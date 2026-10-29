#!/usr/bin/env python3
"""Check v3.4 annotation loading and two-method behavior with simulated signals."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from hybcapseq_quant.core import read_tsv, write_tsv  # noqa: E402
from hybcapseq_quant.references import database_catalog  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Simulation output directory must be empty")
    inputs = output / "inputs"
    spec = importlib.util.spec_from_file_location(
        "analytical_demo", ROOT / "examples/quant/generate_demo.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sample, rows = module.create(inputs, sid="SIMULATED_V34")
    catalog, context = database_catalog(args.db_dir)
    _, metadata = read_tsv(args.db_dir / "virus_pathogen_database.all_metadata.tsv")
    by_id = {r["Accession"]: r for r in metadata}
    legacy = next(
        a
        for a, r in catalog.items()
        if not a.startswith(("PCRROI_", "ERCC-"))
        and r["molecule_type"] == "ssRNA"
        and r["assembly_reference_count"] == "1"
        and r["assembly_mass_status"] == "resolved_complete_assembly"
    )
    roi = [a for a in catalog if a.startswith("PCRROI_")]
    rows = [r for r in rows if r["Accession"].startswith("ERCC-")]
    target_fields, _ = read_tsv(ROOT / "src/hybcapseq_quant/data/targets_template.tsv")
    overrides = []
    for acc in [legacy] + roi:
        record = by_id[acc]
        length = int(record["Length"])
        # Analytical RNA-like response only; not a measurement of target chemistry.
        bases = 10000 * length * 340 / module.AVOGADRO * 1e18
        detection = dict(rows[0])
        detection.update(record)
        detection.update(
            {
                "sample_ID": "SIMULATED_V34",
                "description": "SIMULATED SOFTWARE TEST ONLY",
                "mean_coverage": bases / length,
                "covered_bases": length,
                "read_count": max(1, int(bases / 150)),
                "RPKMF": 1,
            }
        )
        rows.append(detection)
        if acc in roi:
            override = dict.fromkeys(target_fields, "NA")
            override.update(
                {
                    "accession": acc,
                    "alignment_qc_verified": "0",
                    "detection_qc_verified": "0",
                    "calibration_transfer_verified": "0",
                }
            )
            overrides.append(override)
    write_tsv(inputs / "SIMULATED_V34.detected_virus.info.tsv", rows)
    write_tsv(inputs / "roi_overrides.tsv", overrides, target_fields)
    result = output / "results"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "hybcapseq_quant",
            "run",
            "--samples",
            str(inputs / "samples.tsv"),
            "--esviritu-dir",
            str(inputs),
            "--db-dir",
            str(args.db_dir),
            "--targets",
            str(inputs / "roi_overrides.tsv"),
            "--outdir",
            str(result),
            "--method",
            "both",
            "--min-r2",
            "0.80",
            "--reporting-mode",
            "exploratory",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    methods = {}
    for method in ["sequencing_yield", "depth_calibration"]:
        matches = list((result / method).rglob("*.target_quantification.tsv"))
        if not matches:
            matches = list(result.rglob("*target_quantification.tsv"))
            matches = [p for p in matches if method in p.parts]
        if len(matches) != 1:
            raise ValueError("Cannot identify simulated method output: " + method)
        _, targets = read_tsv(matches[0])
        by_group = {r["target_group"]: r for r in targets}
        legacy_value = float(
            by_group[by_id[legacy]["Assembly"]]["estimated_extract_copies_ul"]
        )
        if abs(legacy_value / 5000 - 1) > 1e-8:
            raise ValueError("Legacy analytical truth differs")
        roi_values = [by_group[by_id[a]["Assembly"]] for a in roi]
        key = "estimated_extract_copies_ul"
        if method == "sequencing_yield":
            if any(r[key] != "NA" for r in roi_values):
                raise ValueError("ROI inherited a mass model")
        else:
            if any(abs(float(r[key]) / 5000 - 1) > 1e-8 for r in roi_values):
                raise ValueError("Local depth analytical response differs")
        methods[method] = {
            "legacy_extract_copies_ul": legacy_value,
            "roi_rows_checked": len(roi_values),
            "expected_roi_behavior": (
                "mass_unavailable"
                if method == "sequencing_yield"
                else "qualified_depth_response"
            ),
        }
    report = {
        "result": "PASS",
        "data_kind": "explicitly simulated software check; not experimental validation",
        "database_reference_count": context["reference_count"],
        "methods": methods,
    }
    (output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
