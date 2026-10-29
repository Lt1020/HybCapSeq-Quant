#!/usr/bin/env python3
"""Build the reviewed literature-ROI extension of the corrected v3.3 baseline.

No sequence retrieval or oligo design occurs here. Input files are immutable,
checksum-verified reference snapshots. Large alignment assets belong in Releases.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from hybcapseq_quant.reference_updates import (  # noqa: E402
    bundled_reference_files,
    corrected_metadata_bytes,
    read_corrections,
)
from hybcapseq_quant.references import (  # noqa: E402
    CATALOG_NAME,
    MANIFEST_NAME,
    load_catalog,
)
from build_v3_3 import (  # noqa: E402
    COLS,
    FASTA,
    META,
    INDEX,
    append_reference,
    catalog,
    checksum,
    metadata_matches,
    read_meta,
    require,
)


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def read_table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def table_bytes(fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def sequences(path):
    result, current = {}, None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            require(current not in result, "Duplicate reference ID")
            result[current] = ""
        elif line.strip():
            require(current is not None, "Sequence precedes FASTA header")
            result[current] += line.strip()
    return result


def verify_package(package):
    """Audit all frozen package bytes before executing its independent validator."""
    for line in (package / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        path = (package / name).resolve()
        require(path.is_relative_to(package.resolve()), "Unsafe checksum path")
        require(checksum(path) == digest, "Reference package checksum: " + name)
    subprocess.run(
        [
            sys.executable,
            str(package / "scripts/validate_reference_package.py"),
            "--package",
            str(package),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    summary = json.loads((package / "validation_summary.json").read_text())
    require(summary["primary_record_count"] == 21, "Expected 21 primary references")
    require(summary["candidate_record_count"] == 34, "Expected 34 candidate references")
    return summary


def extend_annotations(stage, roi_rows, region_by_id):
    """Preserve legacy assignments; local ROIs do not inherit genome mass models."""
    data, manifest = bundled_reference_files()
    reader = csv.DictReader(io.StringIO(data.decode()), delimiter="\t")
    fields, base_rows = reader.fieldnames, list(reader)
    additions = []
    for metadata in roi_rows:
        acc = metadata["Accession"]
        region = region_by_id[acc]
        row = dict.fromkeys(fields, "NA")
        row.update(
            {
                "accession": acc,
                "assembly": metadata["Assembly"],
                "reference_name": metadata["Name"],
                "length_nt": metadata["Length"],
                "assembly_length_nt": metadata["Asm_length"],
                "family_raw": metadata["family"],
                "family_normalized": metadata["family"].removeprefix("f__"),
                "ncbi_taxid": metadata["TaxID"],
                "ncbi_accession_version": region["source_accession_version"],
                "ncbi_length_nt": region["source_length_nt"],
                "ncbi_retrieval_status": "fixed_source_interval_verified",
                "ncbi_source_url": region["source_record_url"],
                "classification_evidence": "reviewed_literature_reference_interval",
                "annotation_status": "local_reference_requires_target_mass_assignment",
                "notes": (
                    "reference_scope=literature_PCR_ROI;not_a_complete_genome;"
                    "actual_commercial_assay_unconfirmed;"
                    f"source_interval={region['reference_start_1based']}-"
                    f"{region['reference_end_1based']};"
                    "independent_local_response_validation_required"
                ),
                "assembly_reference_count": "1",
                "assembly_mass_status": "unresolved_assembly_member_mass",
            }
        )
        additions.append(row)
    # Keep every prior annotation byte unchanged; append only new local rows.
    extra = table_bytes(fields, additions).split(b"\n", 1)[1]
    (stage / CATALOG_NAME).write_bytes(data + extra)
    manifest.update(
        {
            "database_revision": "v3.4",
            "created_at_utc": json.loads(
                (ROOT / "databases/v3.4/build_config.json").read_text()
            )["created_at_utc"],
            "catalog_sha256": checksum(stage / CATALOG_NAME),
            "database_metadata_sha256": checksum(stage / META),
            "reference_count": len(base_rows) + len(additions),
            "roi_reference_count": len(additions),
            "roi_mass_assignment_policy": "Unavailable by default; no whole-genome model inheritance",
            "reference_scope": "whole_reference_baseline_plus_literature_diagnostic_ROIs",
        }
    )
    manifest["annotation_status_counts"] = dict(
        Counter(r["annotation_status"] for r in base_rows + additions)
    )
    manifest["molecule_type_counts"] = dict(
        Counter(r["molecule_type"] for r in base_rows + additions)
    )
    write_json(stage / MANIFEST_NAME, manifest)
    load_catalog(stage / CATALOG_NAME, stage / META)
    return manifest


def check_index(stage, package, work, binary, threads, expected):
    roi = sequences(package / FASTA)
    ercc = sequences(ROOT / "references/ercc96/data/certified_v1" / FASTA)
    queries = {}
    query_path = work / "index_queries.fna"
    with query_path.open("w") as handle:
        for prefix, records in [("roi", roi), ("ercc", ercc)]:
            for acc, sequence in records.items():
                # Long enough for sr seeds, unambiguous, and present in the input.
                length = min(150, len(sequence))
                query = next(
                    (
                        sequence[i : i + length]
                        for i in range(len(sequence) - length + 1)
                        if set(sequence[i : i + length]) <= set("ACGT")
                    ),
                    None,
                )
                require(query is not None, "No unambiguous index query: " + acc)
                name = prefix + "_check_" + acc
                queries[name] = (acc, prefix)
                handle.write(f">{name}\n{query}\n")
    sam_path = work / "index_queries.sam"
    with sam_path.open("w") as sam, (work / "index_queries.log").open("w") as log:
        subprocess.run(
            [
                str(binary),
                "-a",
                "-x",
                "sr",
                "--secondary=yes",
                "-N",
                "1000",
                "-t",
                str(threads),
                str(stage / INDEX),
                str(query_path),
            ],
            stdout=sam,
            stderr=log,
            check=True,
        )
    indexed, hits, target_hits = {}, defaultdict(set), set()
    with sam_path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if line.startswith("@SQ\t"):
                values = dict(v.split(":", 1) for v in fields[1:])
                require(values["SN"] not in indexed, "Duplicate indexed ID")
                indexed[values["SN"]] = int(values["LN"])
            elif not line.startswith("@") and not int(fields[1]) & 4:
                hits[fields[0]].add(fields[2])
                if queries[fields[0]][0] == fields[2]:
                    target_hits.add(fields[0])
    require(indexed == expected, "Index and FASTA ID/length/order mismatch")
    require(
        target_hits == set(queries),
        "An expected reference was not recovered in index queries",
    )
    audit = [
        {
            "query": name,
            "expected_reference": acc,
            "kind": kind,
            "expected_reference_hit": "true",
            "reported_reference_hits": len(hits[name]),
            "multiple_reference_hits": str(len(hits[name]) > 1).lower(),
        }
        for name, (acc, kind) in queries.items()
    ]
    (stage / "index_query_audit.tsv").write_bytes(table_bytes(list(audit[0]), audit))
    return {
        "index_reference_count": len(indexed),
        "ercc_queries_passed": len(ercc),
        "roi_queries_passed": len(roi),
        "roi_queries_with_multiple_reported_hits": sum(
            a["kind"] == "roi" and a["multiple_reference_hits"] == "true" for a in audit
        ),
        "secondary_alignment_limit": 1000,
        "interpretation": "Reference recovery check; neither uniqueness nor assay specificity validation",
    }


def build(args):
    base, package = args.base_dir.resolve(), args.reference_package.resolve()
    output, work = args.output_dir.resolve(), args.work_dir.resolve()
    require(not output.exists(), "Refuse to overwrite an existing database directory")
    require(args.threads > 0, "Thread count must be positive")
    require(
        subprocess.check_output([str(args.minimap2), "--version"], text=True).strip()
        == "2.31-r1302",
        "Use the pinned minimap2 2.31-r1302 build",
    )
    base_manifest = json.loads(
        (ROOT / "databases/v3.3/esviritu_db_v3.3_manifest.json").read_text()
    )
    for item in base_manifest["files"]:
        # Accept the original v3.3 or its documented v3.3.1 metadata revision.
        if item["name"] == META:
            continue
        require(
            checksum(base / item["name"]) == item["sha256"],
            "Pinned base file differs: " + item["name"],
        )
    specification = read_corrections()
    corrected = corrected_metadata_bytes((base / META).read_bytes(), specification)
    require(
        hashlib.sha256(corrected).hexdigest()
        == specification["corrected_metadata_sha256"],
        "Corrected base metadata differs",
    )
    summary = verify_package(package)
    base_refs, roi_refs = catalog(base / FASTA), catalog(package / FASTA)
    base_rows = list(csv.DictReader(io.StringIO(corrected.decode()), delimiter="\t"))
    roi_rows = read_meta(package / META)
    metadata_matches(base_rows, base_refs)
    metadata_matches(roi_rows, roi_refs)
    require(len(base_refs) == 31003 and len(roi_refs) == 21, "Unexpected input counts")
    require(
        sum(a.startswith("ERCC-") for a in base_refs) == 96, "Base ERCC count differs"
    )
    for field in ["Accession", "Name", "Assembly"]:
        require(
            not ({r[field] for r in base_rows} & {r[field] for r in roi_rows}),
            "Base/ROI collision: " + field,
        )
    require(
        all(
            r["Assembly"] == r["Name"] and r["Length"] == r["Asm_length"]
            for r in roi_rows
        ),
        "ROI grouping/length differs",
    )
    _, regions = read_table(package / "reference_regions.tsv")
    region_by_id = {
        r["reference_id"]: r for r in regions if r["reference_id"] in roi_refs
    }
    require(set(region_by_id) == set(roi_refs), "ROI coordinate annotations incomplete")
    require(
        all(
            r["sequence_check_status"] == "full100_compatible_sites"
            for r in region_by_id.values()
        ),
        "Non-primary reference in merge",
    )
    require(
        sum(r["target"] == "SARS-CoV-2" for r in region_by_id.values()) == 4,
        "Expected exactly four SARS-CoV-2 ROI records",
    )
    work.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="v3.4_", dir=work))
    try:
        append_reference(base / FASTA, package / FASTA, stage / FASTA)
        extra = (package / META).read_bytes().split(b"\n", 1)[1]
        (stage / META).write_bytes(corrected + extra)
        expected = {**base_refs, **roi_refs}
        require(catalog(stage / FASTA) == expected, "Combined FASTA catalog differs")
        metadata_matches(read_meta(stage / META), expected)
        require(
            checksum(stage / FASTA, limit=(base / FASTA).stat().st_size)
            == checksum(base / FASTA),
            "Base FASTA bytes changed",
        )
        require(
            (stage / META).read_bytes().startswith(corrected),
            "Corrected base metadata bytes changed",
        )
        annotations = extend_annotations(stage, roi_rows, region_by_id)
        shutil.copy2(
            ROOT
            / "src/hybcapseq_quant/data/reference/virus_pathogen_database.reference_corrections.json",
            stage,
        )
        with (work / "build_index.log").open("w") as log:
            subprocess.run(
                [
                    str(args.minimap2),
                    "-x",
                    "sr",
                    "-t",
                    str(args.threads),
                    "-d",
                    str(stage / INDEX),
                    str(stage / FASTA),
                ],
                stdout=log,
                stderr=log,
                check=True,
            )
        index_validation = check_index(
            stage, package, work, args.minimap2, args.threads, expected
        )
        roi_bytes = table_bytes(list(regions[0]), [region_by_id[a] for a in roi_refs])
        (stage / "virus_pathogen_database.pcr_roi_regions.tsv").write_bytes(roi_bytes)
        shutil.copytree(package, stage / "pcr_roi_reference_package")
        readme = (ROOT / "databases/v3.4/README.md").read_text()
        readme = readme.replace(
            "../../references/pcr_roi/README_zh.md",
            "pcr_roi_reference_package/README_zh.md",
        )
        readme = readme.replace("../../docs/", "documentation/")
        (stage / "README_v3.4.md").write_text(readme)
        (stage / "documentation").mkdir()
        for name in [
            "WORKFLOW_v3.4_zh.md",
            "METHODS.md",
            "MANUSCRIPT_METHODS.md",
            "TARGET_CORRECTION.md",
        ]:
            shutil.copy2(ROOT / "docs" / name, stage / "documentation" / name)
        report = {
            "database_version": "v3.4",
            "build_date": "2026-10-29",
            "release_kind": "custom literature diagnostic ROI extension; not an upstream EsViritu release",
            "base_database": "v3.3 with verified v3.3.1 metadata and molecular annotations",
            "base_archive_sha256": base_manifest["archive"]["sha256"],
            "base_fasta_sha256": checksum(base / FASTA),
            "base_corrected_metadata_sha256": hashlib.sha256(corrected).hexdigest(),
            "reference_package_manifest_sha256": checksum(package / "SHA256SUMS"),
            "counts": {
                "base_records": len(base_refs),
                "base_ercc_records": 96,
                "added_primary_roi_records": len(roi_refs),
                "combined_records": len(expected),
                "combined_reference_bases": sum(expected.values()),
                "roi_targets_in_main_database": len(
                    {r["target"] for r in region_by_id.values()}
                ),
                "literature_targets_in_audit": summary["source_target_count"],
                "literature_assays_in_audit": summary["assay_record_count"],
                "excluded_roi_candidates_retained_in_audit": 13,
                "sars_cov2_roi_records": 4,
            },
            "minimap2": {"version": "2.31-r1302", "preset": "sr", "k": 21, "w": 11},
            "validation": {
                "base_fasta_bytes_preserved": True,
                "corrected_base_metadata_bytes_preserved": True,
                "only_documented_eight_base_metadata_cells_changed": True,
                "fasta_metadata_index_ids_and_lengths_match": True,
                "molecular_annotation_catalog_validated": True,
                "all_original_molecular_annotation_bytes_preserved": True,
                "reference_package_checksums_and_independent_validator_passed": True,
                "esviritu_experimental_fastq_end_to_end_performed": False,
                "roi_quantitative_transfer_validated": False,
                "actual_commercial_dpcr_assay_confirmed": False,
                **index_validation,
            },
            "roi_measurement_policy": {
                "reference_scope": "published diagnostic region including exactly 100 nt per flank",
                "mass_assignment": "NA by default; explicitly reviewed local-molecule assignments required for yield estimates",
                "depth_estimates": "region response only; exploratory until independently validated",
                "counting": "do not sum ROI and genomic-reference hits or overlapping assay records",
                "calibration_factors": "no inheritance of whole-reference correction factors",
            },
            "files": [
                {
                    "name": p.relative_to(stage).as_posix(),
                    "bytes": p.stat().st_size,
                    "sha256": checksum(p),
                }
                for p in sorted(stage.rglob("*"))
                if p.is_file()
                and not any(part.startswith(".") for part in p.relative_to(stage).parts)
            ],
        }
        if args.expected_manifest:
            locked = json.loads(args.expected_manifest.read_text())
            for key in [
                "counts",
                "files",
                "base_fasta_sha256",
                "base_corrected_metadata_sha256",
                "reference_package_manifest_sha256",
            ]:
                require(
                    report[key] == locked[key],
                    "Rebuild differs from locked manifest: " + key,
                )
        write_json(stage / "esviritu_db_v3.4_manifest.json", report)
        (stage / "SHA256SUMS").write_text(
            "".join(
                checksum(p) + "  " + p.relative_to(stage).as_posix() + "\n"
                for p in sorted(stage.rglob("*"))
                if p.is_file()
                and p != stage / "SHA256SUMS"
                and not any(part.startswith(".") for part in p.relative_to(stage).parts)
            )
        )
        stage.replace(output)
        print(
            json.dumps(
                {
                    "output": str(output),
                    "counts": report["counts"],
                    "index_checks": index_validation,
                },
                indent=2,
            )
        )
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["base-dir", "reference-package", "output-dir", "work-dir", "minimap2"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--expected-manifest", type=Path)
    build(parser.parse_args())


if __name__ == "__main__":
    main()
