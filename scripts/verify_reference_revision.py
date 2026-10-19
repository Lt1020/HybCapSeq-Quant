"""Recalculate an archived batch and verify the scope of reference corrections.

Source detection columns are recovered unchanged from the archived reference
tables. Total input bases are recovered from the saved analysis manifests.
No experimental parameter or alignment observation is inferred or adjusted.
"""

import argparse
import csv
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

from hybcapseq_quant.cli import main as quantify, validate_sample_id
from hybcapseq_quant.core import read_tsv, sha256, write_tsv
from hybcapseq_quant.reference_updates import read_corrections

SOURCE_FIELDS = [
    "sample_ID",
    "Name",
    "description",
    "Length",
    "Segment",
    "Accession",
    "Assembly",
    "Asm_length",
    "kingdom",
    "phylum",
    "tclass",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "RPKMF",
    "read_count",
    "covered_bases",
    "mean_coverage",
    "avg_read_identity",
    "Pi",
    "filtered_reads_in_sample",
]
METHODS = ["sequencing_yield", "depth_calibration"]


def same_number(left, right):
    if left in {"NA", ""} or right in {"NA", ""}:
        return left == right
    return math.isclose(float(left), float(right), rel_tol=1e-10, abs_tol=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-results", type=Path, required=True)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="verify an already recalculated batch without regenerating plots",
    )
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.reuse_existing:
        raise ValueError("Verification output directory must be empty")
    args.outdir.mkdir(parents=True, exist_ok=True)
    inputs = args.outdir / "inputs"
    inputs.mkdir(exist_ok=args.reuse_existing)
    detections = args.outdir / "detection_tables"
    detections.mkdir(exist_ok=args.reuse_existing)
    old = {}
    with ZipFile(args.source_results) as archive:
        manifests = [n for n in archive.namelist() if n.endswith("/run_manifest.json")]
        if len(manifests) != 1:
            raise ValueError("Expected one batch manifest in the source archive")
        prefix = manifests[0].removesuffix("run_manifest.json")

        def read(name):
            return list(
                csv.DictReader(
                    io.StringIO(archive.read(prefix + name).decode("utf-8-sig")),
                    delimiter="\t",
                )
            )

        samples = read("inputs/samples.tsv")
        (inputs / "ercc_standards.tsv").write_bytes(
            archive.read(prefix + "inputs/ercc_standards.tsv")
        )
        for sample in samples:
            sid = sample["sample_id"]
            validate_sample_id(sid)
            source = read(f"depth_calibration/{sid}/{sid}.reference_quantification.tsv")
            write_tsv(
                detections / f"{sid}.detected_virus.info.tsv", source, SOURCE_FIELDS
            )
            original = json.loads(
                archive.read(
                    prefix + f"depth_calibration/{sid}/{sid}.analysis_manifest.json"
                )
            )
            sample["total_effective_bases"] = str(original["total_effective_bases"])
            sample["info_tsv"] = (
                "../detection_tables/" + sid + ".detected_virus.info.tsv"
            )
            sample["fastp_json"] = "NA"
            for method in METHODS:
                old[method, sid] = {
                    "groups": {
                        r["target_group"]: r
                        for r in read(f"{method}/{sid}/{sid}.target_quantification.tsv")
                    },
                    "model": read(f"{method}/{sid}/{sid}.calibration_model.tsv")[0],
                }
    write_tsv(inputs / "samples.tsv", samples)
    results = args.outdir / "results"
    status = (
        0
        if args.reuse_existing
        else quantify(
            [
                "run",
                "--samples",
                str(inputs / "samples.tsv"),
                "--ercc-standards",
                str(inputs / "ercc_standards.tsv"),
                "--db-dir",
                str(args.db_dir),
                "--outdir",
                str(results),
                "--method",
                "both",
                "--reporting-mode",
                "exploratory",
                "--wastewater-unit",
                "copies/L",
            ]
        )
    )
    if status:
        raise RuntimeError("Batch recalculation failed")
    corrected = {c["assembly"] for c in read_corrections()["corrections"]}
    audit = {
        "source_archive_sha256": sha256(args.source_results),
        "sample_count": len(samples),
        "experimental_parameters_changed": False,
        "input_reconstruction": "Unchanged source detection columns; total bases copied from archived analysis manifests; paths made portable.",
        "methods": {},
    }
    for method in METHODS:
        summary = {
            "target_groups": 0,
            "newly_computable_corrected_groups": 0,
            "unchanged_uncorrected_group_estimates": 0,
            "reported_concentrations": 0,
            "remaining_NA": 0,
            "calibration_models_unchanged": True,
        }
        for sample in samples:
            sid = sample["sample_id"]
            p = results / method / sid
            _, new = read_tsv(p / (sid + ".target_quantification.tsv"))
            _, model_rows = read_tsv(p / (sid + ".calibration_model.tsv"))
            model = model_rows[0]
            for key in ["slope", "intercept", "r2", "n_fitted"]:
                if not same_number(
                    model.get(key, "NA"), old[method, sid]["model"].get(key, "NA")
                ):
                    raise AssertionError(
                        "Calibration changed after a reference-only correction"
                    )
            for row in new:
                group = row["target_group"]
                before = old[method, sid]["groups"][group]
                summary["target_groups"] += 1
                if group in corrected:
                    if (
                        before["estimated_extract_copies_ul"] != "NA"
                        or row["estimated_extract_copies_ul"] == "NA"
                    ):
                        raise AssertionError(
                            "Expected a corrected group to become computable"
                        )
                    if (
                        row["reference_length_nt"] != "13859"
                        or "inconsistent_database_assembly_length"
                        in row["quantification_status"]
                    ):
                        raise AssertionError(
                            "Corrected reference length is not used consistently"
                        )
                    summary["newly_computable_corrected_groups"] += 1
                else:
                    for field in [
                        "estimated_reaction_copies_ul",
                        "estimated_extract_copies_ul",
                        "estimated_wastewater_concentration",
                    ]:
                        if not same_number(row[field], before[field]):
                            raise AssertionError("An unrelated target estimate changed")
                    summary["unchanged_uncorrected_group_estimates"] += 1
                summary[
                    (
                        "remaining_NA"
                        if row["concentration_reaction_copies_ul"] == "NA"
                        else "reported_concentrations"
                    )
                ] += 1
            _, generated = read_tsv(p / (sid + ".reference_quantification.tsv"))
            _, source = read_tsv(detections / (sid + ".detected_virus.info.tsv"))
            if len(generated) != len(source) or any(
                a[field] != b[field]
                for a, b in zip(generated, source)
                for field in SOURCE_FIELDS
            ):
                raise AssertionError("Original source columns were changed")
        audit["methods"][method] = summary
    audit["source_columns_preserved"] = True
    (args.outdir / "verification_summary.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
