#!/usr/bin/env python3
"""Reproduce the NIST ERCC96 reference tables using Python 3.10+ only.

Never infer taxonomy, alter bases, trim tails, or include the vector in ERCC96.
Sources are pinned to their downloaded SHA256. Run from any working directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parent
COLS = ["Accession", "description", "Name", "Segment", "kingdom", "phylum",
        "tclass", "order", "family", "genus", "species", "subspecies",
        "Length", "TaxID", "Assembly", "Asm_length"]
BASE = "https://tsapps.nist.gov/srmext/certificates/documents/"
CERT_URL = "https://tsapps.nist.gov/srmext/certificates/2374.pdf"
PINS = {
    "SRM2374_Sequence_v1.FASTA": "337e7121586349129ff9a0004f34d7d3d003c0bd28598e70d76ed48d5ad3778f",
    "SRM2374_ambiguities_resolved_v2.FASTA": "be3b785bf72b5b2b0cb3bcb1e022223cc76fb00ac0d32212a945719799222f8a",
    "SRM2374_ambiguities_resolved_v1.GTF": "4f0e4278811d9dcf5a0ce42bdb8a126816a2546bf16d385f3d541ede1e6b493a",
}
VECTOR = "080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev"
TABLE2_SHA256 = "79624cad5923aab4cf2d30a06b67fa9399650259540d1d840e021e7175ef2dd9"
IUPAC = dict(zip("ACGTRYSWKMBDHVN", ["A", "C", "G", "T", "AG", "CT",
             "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"]))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:]
            require(bool(current) and not any(c.isspace() for c in current),
                    f"Unexpected FASTA header: {current!r}")
            require(current not in records, f"Duplicate FASTA ID: {current}")
            records[current] = ""
        elif line.strip():
            require(current is not None, "Sequence before FASTA header")
            require(not (set(line) - set(IUPAC)), f"Invalid base in {current}")
            records[current] += line
    require(all(records.values()), f"Empty FASTA sequence: {path}")
    return records


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build(out: Path) -> dict:
    src = ROOT / "sources"
    for name, expected in PINS.items():
        require(sha(src / name) == expected, f"Official source SHA256 mismatch: {name}")
    require(sha(src / "NIST_Table2_extracted.tsv") == TABLE2_SHA256, "Extracted Table 2 SHA256 differs")
    if (src / "2374.pdf").exists():
        require(sha(src / "2374.pdf") == "9add6fbc3eb1c801595ce829021f9982d17c762fc64505af1e937172332c9a21", "Certificate PDF SHA256 differs")
    certified = fasta(src / "SRM2374_Sequence_v1.FASTA")
    resolved_all = fasta(src / "SRM2374_ambiguities_resolved_v2.FASTA")
    require(len(certified) == 96, "Certified source must contain exactly 96 records")
    require(all(re.fullmatch(r"ERCC-\d{5}", i) for i in certified), "Non-ERCC certified ID")
    require(set(resolved_all) - set(certified) == {VECTOR}, "Unexpected secondary-source record set")
    require(set(certified) <= set(resolved_all), "Missing resolved ERCC sequence")
    resolved = {i: resolved_all[i] for i in certified}
    require(all(set(s) <= set("ACGT") for s in resolved.values()), "Resolved source has ambiguity")
    require(all(len(certified[i]) == len(resolved[i]) for i in certified), "Version lengths differ")

    gtf = {}
    for line in (src / "SRM2374_ambiguities_resolved_v1.GTF").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        require(len(fields) == 9, "GTF does not have nine columns")
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
        require(fields[0] not in gtf, "Duplicate GTF reference")
        gtf[fields[0]] = (fields, attrs)
    require(set(gtf) == set(resolved_all), "GTF/secondary FASTA record set differs")
    nist_rows = table(src / "NIST_Table2_extracted.tsv")
    nist = {r["ERCC_ID"]: r for r in nist_rows}
    require(len(nist_rows) == len(nist) == 96 and set(nist) == set(certified), "Table 2 IDs differ")
    for i, sequence in certified.items():
        row = nist[i]
        fields, attrs = gtf[i]
        require(int(row["Insert_length_nt"]) == len(sequence), f"Table 2 length differs: {i}")
        require(int(fields[3]) == 1 and int(fields[4]) == len(sequence), f"GTF length differs: {i}")
        require(attrs["gene_id"] == i, f"GTF gene ID differs: {i}")
        require(attrs["transcript_id"] == row["GenBank_accession_reported"], f"GenBank mapping differs: {i}")

    out.mkdir(parents=True, exist_ok=True)
    variants = []
    for label, records, source_name in [
        ("certified_v1", certified, "SRM2374_Sequence_v1.FASTA"),
        ("ambiguities_resolved_v2", resolved, "SRM2374_ambiguities_resolved_v2.FASTA"),
    ]:
        target = out / label
        target.mkdir(exist_ok=True)
        fna = target / "virus_pathogen_database.fna"
        if label == "certified_v1":
            shutil.copyfile(src / source_name, fna)  # Byte-for-byte official original.
        else:
            with fna.open("w", encoding="ascii", newline="\n") as f:
                for i, sequence in records.items():
                    f.write(f">{i}\n")
                    for offset in range(0, len(sequence), 60):
                        f.write(sequence[offset:offset + 60] + "\n")
        meta, provenance = [], []
        for i, sequence in records.items():
            row = dict.fromkeys(COLS, "NA")
            # One ERCC transcript per computational assembly; user-defined grouping.
            row.update(Accession=i, description=i, Name=i, Length=len(sequence),
                       Assembly=i, Asm_length=len(sequence))
            meta.append(row)
            provenance.append({
                "Accession": i, "official_fasta_header": i,
                "GenBank_accession_reported": nist[i]["GenBank_accession_reported"],
                "reference_length_nt": len(sequence),
                "ambiguous_bases": sum(c not in "ACGT" for c in sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "source_filename": source_name, "source_url": BASE + source_name,
                "mapping_source_url": BASE + "SRM2374_ambiguities_resolved_v1.GTF",
            })
        write_tsv(target / "virus_pathogen_database.all_metadata.tsv", COLS, meta)
        write_tsv(target / "sequence_provenance.tsv", list(provenance[0]), provenance)
        check = table(target / "virus_pathogen_database.all_metadata.tsv")
        require(list(check[0]) == COLS, "Metadata schema/order differs")
        require([r["Accession"] for r in check] == list(fasta(fna)), "FASTA/metadata order differs")
        require(all(r["Assembly"] == r["Name"] for r in check), "Assembly/Name differs")
        require(all(r["Asm_length"] == r["Length"] for r in check), "Asm_length/Length differs")
        variants.append({"name": label, "records": len(records),
                         "total_reference_nt": sum(map(len, records.values())),
                         "length_min_nt": min(map(len, records.values())),
                         "length_max_nt": max(map(len, records.values())),
                         "ambiguous_bases": sum(c not in "ACGT" for s in records.values() for c in s)})

    differences = []
    for i, sequence in certified.items():
        for pos, (a, b) in enumerate(zip(sequence, resolved[i]), 1):
            if a != b:
                differences.append({"Accession": i, "position_1based": pos,
                                    "certified_v1_symbol": a, "official_resolved_v2_base": b,
                                    "base_allowed_by_certified_IUPAC": str(b in IUPAC[a]).lower()})
    write_tsv(out / "official_version_differences.tsv", list(differences[0]), differences)
    excluded = [{"source_filename": name, "record_id": VECTOR, "length_nt": len(resolved_all[VECTOR]),
                 "reason": "Not an identifier in the 96-record certified ERCC source; vector reference"}
                for name in ["SRM2374_ambiguities_resolved_v2.FASTA", "SRM2374_ambiguities_resolved_v1.GTF"]]
    write_tsv(out / "excluded_records.tsv", list(excluded[0]), excluded)
    shutil.copyfile(src / "NIST_Table2_extracted.tsv", out / "NIST_Table2_records.tsv")
    warnings = [
        {"code": "OFFICIAL_IUPAC_RESOLUTION_DISAGREEMENT", "accession": "ERCC-00138", "position_1based": 997,
         "detail": "Certified FASTA has R (A/G); official resolved v2 FASTA has T. Both sources retained unchanged."},
        {"code": "OFFICIAL_CERTIFICATE_POSITION_DISAGREEMENT", "accession": "ERCC-00018",
         "detail": "Certified FASTA has its sole R at position 540; certificate Table 2 footnote (e) says 599. No sequence or position was corrected."},
    ]
    incompatibilities = [r for r in differences if r["base_allowed_by_certified_IUPAC"] == "false"]
    require(len(differences) == 11, "Unexpected version difference count")
    require(len(incompatibilities) == 1 and incompatibilities[0]["Accession"] == "ERCC-00138"
            and incompatibilities[0]["position_1based"] == 997, "Unexpected IUPAC discrepancy")
    source_info = [{"filename": n, "url": BASE + n, "sha256": expected,
                    "size_bytes": (src / n).stat().st_size} for n, expected in PINS.items()]
    source_info.append({"filename": "2374.pdf", "url": CERT_URL,
                        "sha256": "9add6fbc3eb1c801595ce829021f9982d17c762fc64505af1e937172332c9a21",
                        "role": "NIST certificate, issued 2017-12-06; Table 2 on pages 7-10"})
    manifest = {
        "dataset": "NIST SRM2374 ERCC96, EsViritu-style truthful annotations",
        "dataset_version": "1.0.1", "retrieved_date_utc": "2026-10-17",
        "sources": source_info, "variants": variants,
        "metadata_columns": COLS, "missing_value": "NA",
        "populated_columns": {"Accession": "Official ERCC identifier (not an NCBI accession)",
                              "description": "Exact official FASTA header without >",
                              "Name": "Official ERCC identifier", "Length": "Counted from the delivered reference",
                              "Assembly": "Copied from Name; user-defined per-transcript computational grouping, not an NCBI assembly accession",
                              "Asm_length": "Copied from Length; one reference per computational assembly"},
        "metadata_schema_source": "https://github.com/cmmr/EsViritu/blob/42aad529d9112451531c81633171bc9a7e1f630b/docs/custom-database.md",
        "validation": {"unique_certified_ids": 96, "metadata_fasta_id_length_match": True,
                       "Assembly_equals_Name_for_all_records": True,
                       "Asm_length_equals_Length_for_all_records": True,
                       "table2_gtf_accession_length_matches": 96,
                       "certified_output_byte_identical_to_source": True,
                       "resolved_ERCC_sequences_identical_to_official_source": True,
                       "extra_vector_records_excluded_from_resolved_reference": 1,
                       "extra_vector_records_excluded_from_ERCC_annotation": 1,
                       "base_differences_between_official_versions": 11,
                       "IUPAC_incompatible_official_resolutions": len(incompatibilities),
                       "esviritu_end_to_end_run_performed": False, "mmi_built": False},
        "warnings": warnings,
        "runtime_note": "Each ERCC record has Assembly=Name and numeric Asm_length=Length, as requested. Assembly is a per-transcript computational grouping label derived from the official ERCC name, not an NCBI assembly accession. The remaining 10 unsupported metadata fields are NA. An EsViritu end-to-end run has not been performed.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    paths = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text("".join(f"{sha(p)}  {p.relative_to(out).as_posix()}\n" for p in paths), encoding="ascii")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "data", help="New or empty output directory")
    args = parser.parse_args()
    if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        parser.error("--out must be a new or empty directory")
    m = build(args.out)
    print(json.dumps({"output": str(args.out.resolve()), "variants": m["variants"],
                      "source_discrepancies_retained": m["warnings"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
