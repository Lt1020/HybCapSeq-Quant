#!/usr/bin/env python3
"""Append unmodified ERCC reference records and rebuild an EsViritu index.

Inputs are checked before writing. The three-file output is staged atomically.
No pathogen sequence or original metadata value is changed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

FASTA = "virus_pathogen_database.fna"
META = "virus_pathogen_database.all_metadata.tsv"
INDEX = "virus_pathogen_database.mmi"
COLS = ["Accession", "description", "Name", "Segment", "kingdom", "phylum",
        "tclass", "order", "family", "genus", "species", "subspecies",
        "Length", "TaxID", "Assembly", "Asm_length"]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def checksum(path, algorithm="sha256", limit=None):
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        while limit is None or limit > 0:
            block = stream.read(1024 * 1024 if limit is None else min(limit, 1024 * 1024))
            if not block:
                break
            h.update(block)
            if limit is not None:
                limit -= len(block)
    return h.hexdigest()


def catalog(path):
    records = {}
    current = None
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith(">"):
                current = line[1:].split()[0]
                require(current not in records, f"Duplicate FASTA ID: {current}")
                records[current] = 0
            elif line.strip():
                require(current is not None, "FASTA sequence before header")
                records[current] += len(line.strip())
    require(records and all(records.values()), "Empty FASTA reference")
    return records


def read_meta(path):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == COLS, f"Unexpected metadata header: {reader.fieldnames}")
        rows = list(reader)
    require(len({r["Accession"] for r in rows}) == len(rows), "Duplicate metadata accession")
    require(all(set(r) == set(COLS) and None not in r.values() for r in rows), "Malformed metadata row")
    return rows


def metadata_matches(rows, references):
    require({r["Accession"] for r in rows} == set(references), "Metadata/FASTA ID sets differ")
    require(all(int(r["Length"]) == references[r["Accession"]] for r in rows), "Metadata/FASTA lengths differ")


def append_reference(base, extra, target):
    with target.open("wb") as out:
        with base.open("rb") as src:
            shutil.copyfileobj(src, out)
        with base.open("rb") as src:
            src.seek(-1, 2)
            if src.read(1) not in (b"\n", b"\r"):
                out.write(b"\n")
        with extra.open("rb") as src:
            shutil.copyfileobj(src, out)


def write_queries(path, source):
    sequences = {}
    current = None
    for line in source.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            sequences[current] = ""
        elif line.strip():
            sequences[current] += line.strip()
    queries = {}
    with path.open("w", encoding="ascii") as out:
        for accession, sequence in sequences.items():
            query = next((sequence[i:i+150] for i in range(len(sequence)-149)
                          if set(sequence[i:i+150]) <= set("ACGT")), None)
            require(query is not None, f"No unambiguous 150 nt verification window: {accession}")
            name = "index_check_" + accession
            queries[name] = accession
            out.write(f">{name}\n{query}\n")
    return queries


def build(args):
    base, ercc, output = args.base_dir.resolve(), args.ercc_dir.resolve(), args.output_dir.resolve()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    require(not output.exists(), "Output directory already exists")
    base_catalog = catalog(base / FASTA)
    ercc_catalog = catalog(ercc / FASTA)
    base_rows, ercc_rows = read_meta(base / META), read_meta(ercc / META)
    metadata_matches(base_rows, base_catalog)
    metadata_matches(ercc_rows, ercc_catalog)
    require(len(base_rows) == 30907, "Unexpected official v3.2.4 record count")
    require(len(ercc_rows) == 96 and all(i.startswith("ERCC-") for i in ercc_catalog), "Expected 96 ERCC records")
    require(not (set(base_catalog) & set(ercc_catalog)), "Base/ERCC accession collision")
    require(all(r["Assembly"] == r["Name"] and r["Asm_length"] == r["Length"] for r in ercc_rows),
            "ERCC Assembly/Name or Asm_length/Length mismatch")
    require(len({r["Assembly"] for r in ercc_rows}) == 96, "ERCC assemblies must be independent")
    expected = {**base_catalog, **ercc_catalog}
    stage = Path(tempfile.mkdtemp(prefix="v3.3_build_", dir=work))
    try:
        append_reference(base / FASTA, ercc / FASTA, stage / FASTA)
        with (stage / META).open("wb") as out:
            with (base / META).open("rb") as src:
                shutil.copyfileobj(src, out)
            with (base / META).open("rb") as src:
                src.seek(-1, 2)
                if src.read(1) not in (b"\n", b"\r"):
                    out.write(b"\n")
            with (ercc / META).open("rb") as src:
                src.readline()  # Preserve the original base header exactly once.
                shutil.copyfileobj(src, out)
        require((stage / FASTA).stat().st_size >= (base / FASTA).stat().st_size + (ercc / FASTA).stat().st_size, "Combined FASTA byte count is incomplete")
        combined = catalog(stage / FASTA)
        combined_rows = read_meta(stage / META)
        require(combined == expected, "Unexpected combined reference catalog")
        metadata_matches(combined_rows, combined)
        require(combined_rows[:len(base_rows)] == base_rows, "Original metadata changed")
        require(combined_rows[len(base_rows):] == ercc_rows, "ERCC metadata changed")
        for name in [FASTA, META]:
            require(checksum(stage / name, limit=(base / name).stat().st_size) == checksum(base / name),
                    f"Original file bytes changed: {name}")
        version = subprocess.check_output([str(args.minimap2), "--version"], text=True).strip()
        index_command = [str(args.minimap2), "-x", "sr", "-t", str(args.threads),
                         "-d", str(stage / INDEX), str(stage / FASTA)]
        with (work / "build_index.log").open("w") as log:
            subprocess.run(index_command, stdout=log, stderr=subprocess.STDOUT, check=True)
        require((stage / INDEX).stat().st_size > 0, "Empty minimap2 index")
        queries = write_queries(work / "ERCC_index_check.fna", ercc / FASTA)
        with (work / "ERCC_index_check.sam").open("w") as sam, (work / "index_check.log").open("w") as log:
            subprocess.run([str(args.minimap2), "-a", "-x", "sr", "-t", str(args.threads),
                            str(stage / INDEX), str(work / "ERCC_index_check.fna")],
                           stdout=sam, stderr=log, check=True)
        indexed_references, matched = {}, set()
        for line in (work / "ERCC_index_check.sam").read_text().splitlines():
            fields = line.split("\t")
            if line.startswith("@SQ\t"):
                values = dict(f.split(":", 1) for f in fields[1:])
                require(values["SN"] not in indexed_references, "Duplicate indexed reference")
                indexed_references[values["SN"]] = int(values["LN"])
            elif not line.startswith("@"):
                flag = int(fields[1])
                if not flag & 4 and fields[2] == queries.get(fields[0]):
                    matched.add(fields[0])
        require(indexed_references == expected, "Index reference IDs/lengths differ from FASTA")
        require(set(queries) == matched, "Not all ERCC verification queries matched their expected references")
        report = {
            "database_version": "v3.3", "release_kind": "custom ERCC extension; not an upstream EsViritu release",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "upstream": {"version": "v3.2.4", "author": "Michael Tisza",
                         "doi": "10.5281/zenodo.17716199", "license": "CC-BY-4.0",
                         "url": "https://zenodo.org/records/17716199/files/esviritu_db_v3.2.4.tar.gz",
                         "archive_md5": "24d85c1ec3cbffff12e921d2f39c91b2"},
            "ercc": {"variant": "NIST SRM2374 certified v1", "records": len(ercc_rows),
                     "fasta_sha256": checksum(ercc / FASTA), "metadata_sha256": checksum(ercc / META),
                     "Assembly_rule": "Assembly = Name", "Asm_length_rule": "Asm_length = Length",
                     "missing_annotation_value": "NA", "reference_bases_added": sum(ercc_catalog.values())},
            "counts": {"upstream_records": len(base_rows), "added_ERCC_records": len(ercc_rows),
                       "combined_records": len(combined), "combined_reference_bases": sum(combined.values()),
                       "metadata_rows": len(combined_rows), "index_reference_records": len(indexed_references)},
            "minimap2": {"version": version, "preset": "sr", "threads": args.threads,
                         "rebuild_command_template": "minimap2 -x sr -d virus_pathogen_database.mmi virus_pathogen_database.fna"},
            "validation": {"FASTA_metadata_ID_length_match": True, "original_FASTA_bytes_preserved": True,
                           "original_metadata_bytes_preserved": True, "ERCC_metadata_preserved": True,
                           "FASTA_index_ID_length_match": True, "ERCC_150nt_checks_expected_reference_hits": len(matched),
                           "EsViritu_end_to_end_test_performed": False},
            "files": [{"name": name, "bytes": (stage / name).stat().st_size,
                       "sha256": checksum(stage / name)} for name in [FASTA, META, INDEX]],
        }
        require(catalog(stage / FASTA) == expected, "Final FASTA catalog changed during indexing")
        require((stage / FASTA).stat().st_size >= (base / FASTA).stat().st_size + (ercc / FASTA).stat().st_size, "Final FASTA is incomplete")
        output.parent.mkdir(parents=True, exist_ok=True)
        stage.rename(output)
        require(catalog(output / FASTA) == expected, "Published FASTA catalog differs")
        for item in report["files"]:
            require((output / item["name"]).stat().st_size == item["bytes"] and checksum(output / item["name"]) == item["sha256"], "Published file changed")
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"output": str(output), "report": str(args.report),
                          "counts": report["counts"], "validation": report["validation"]}, indent=2))
    except BaseException:
        if stage.exists():
            shutil.rmtree(stage)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["base-dir", "ercc-dir", "output-dir", "work-dir", "minimap2", "report"]:
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    require(args.threads > 0, "Threads must be positive")
    build(args)


if __name__ == "__main__":
    main()
