#!/usr/bin/env python3
"""Audit supplied, reviewed reference files; never fetch or generate sequences."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

FIELDS = [
    "target",
    "reference_id",
    "source_accession_version",
    "source_url",
    "assay_region",
    "actual_assay_confirmed",
    "reviewed",
    "redistribution_allowed",
    "fasta_path",
    "fasta_sha256",
    "reference_start_1based",
    "reference_end_1based",
    "amplicon_start_1based",
    "amplicon_end_1based",
]


def fasta_lengths(path):
    lengths, current = {}, None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                current = value[1:].split()[0]
                if current in lengths:
                    raise ValueError("Duplicate FASTA identifier")
                lengths[current] = 0
            else:
                if current is None or set(value.upper()) - set("ACGTURYSWKMBDHVN"):
                    raise ValueError("Invalid FASTA structure or nucleotide alphabet")
                lengths[current] += len(value)
    return lengths


def audit(manifest, base_dir=None):
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise ValueError("Reference manifest schema mismatch")
        rows = list(reader)
    problems = []
    if not rows:
        problems.append("no_reviewed_target_records_supplied")
    if base_dir is None:
        problems.append("base_database_directory_not_supplied")
    else:
        for name in [
            "virus_pathogen_database.fna",
            "virus_pathogen_database.all_metadata.tsv",
            "virus_pathogen_database.molecular_annotations.tsv",
            "virus_pathogen_database.molecular_annotations.manifest.json",
        ]:
            if not (base_dir / name).is_file():
                problems.append("missing_base_file:" + name)
    seen = set()
    for row in rows:
        tag = row["target"] + ":" + row["reference_id"]
        if row["reference_id"] in seen:
            problems.append("duplicate_reference_id:" + tag)
        seen.add(row["reference_id"])
        for flag in ["actual_assay_confirmed", "reviewed", "redistribution_allowed"]:
            if row[flag].lower() not in {"true", "1"}:
                problems.append(flag + "_not_confirmed:" + tag)
        for field in [
            "target",
            "reference_id",
            "source_accession_version",
            "source_url",
            "assay_region",
        ]:
            if not row[field].strip():
                problems.append("missing_" + field + ":" + tag)
        path = Path(row["fasta_path"])
        if not path.is_absolute():
            path = manifest.parent / path
        if not path.is_file():
            problems.append("missing_supplied_fasta:" + tag)
            continue
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["fasta_sha256"]:
            problems.append("fasta_checksum_mismatch:" + tag)
            continue
        try:
            lengths = fasta_lengths(path)
            start, end, left, right = [
                int(row[field])
                for field in [
                    "reference_start_1based",
                    "reference_end_1based",
                    "amplicon_start_1based",
                    "amplicon_end_1based",
                ]
            ]
            if not 1 <= start <= left <= right <= end:
                problems.append("invalid_coordinate_order:" + tag)
            if left - start < 100 or end - right < 100:
                problems.append("supplied_reference_flanks_less_than_100nt:" + tag)
            if lengths.get(row["reference_id"]) != end - start + 1:
                problems.append("supplied_reference_length_mismatch:" + tag)
        except (ValueError, KeyError):
            problems.append("invalid_reference_or_coordinates:" + tag)
    return dict(
        candidate_version="v3.4",
        preflight_passed=not problems,
        release_ready=False,
        sequence_asset_created=False,
        n_supplied_records=len(rows),
        blocking_issues=problems,
        further_requirements=[
            "database_metadata_and_molecular_annotation_join",
            "compatible_alignment_index_rebuild",
            "alignment_and_quantification_validation",
        ],
        note="Passing preflight verifies supplied inputs only; it does not build or release a database.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--base-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.manifest, args.base_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return 0 if result["preflight_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
