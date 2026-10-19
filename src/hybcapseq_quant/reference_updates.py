"""Verified reference corrections and reversible database revision upgrades."""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile

from . import __version__
from .core import read_tsv, sha256, write_tsv
from .references import (
    CATALOG_NAME,
    MANIFEST_NAME,
    METADATA_NAME,
    PACKAGE_REFERENCE,
    database_catalog,
    load_catalog,
)

CORRECTIONS_NAME = "virus_pathogen_database.reference_corrections.json"
REVISION = "v3.3.1"
FASTA_NAME = "virus_pathogen_database.fna"
INDEX_NAME = "virus_pathogen_database.mmi"


def read_corrections():
    return json.loads((PACKAGE_REFERENCE / CORRECTIONS_NAME).read_text())


def corrected_metadata_bytes(source, specification):
    """Change only documented cells in the checksum-identified metadata file."""
    digest = hashlib.sha256(source).hexdigest()
    if digest == specification.get("corrected_metadata_sha256"):
        return source
    if digest != specification["original_metadata_sha256"]:
        raise ValueError(
            "Reference update does not match this database metadata checksum"
        )
    corrections = specification["corrections"]
    lookup = {
        accession: correction
        for correction in corrections
        for accession in correction["members"]
    }
    lines = source.splitlines(keepends=True)
    header = lines[0].decode("utf-8-sig").rstrip("\r\n").split("\t")
    indexes = {
        name: header.index(name)
        for name in ("Accession", "Assembly", "Length", "Asm_length")
    }
    seen = set()
    result = [lines[0]]
    for line in lines[1:]:
        ending = (
            b"\r\n"
            if line.endswith(b"\r\n")
            else (
                b"\n"
                if line.endswith(b"\n")
                else b"\r" if line.endswith(b"\r") else b""
            )
        )
        body = line[: -len(ending)] if ending else line
        parts = body.split(b"\t")
        accession = parts[indexes["Accession"]].decode()
        if accession in lookup:
            correction = lookup[accession]
            if len(parts) != len(header) or accession in seen:
                raise ValueError("Unexpected or duplicate correction metadata row")
            if parts[indexes["Assembly"]].decode() != correction["assembly"]:
                raise ValueError("Correction assembly identity mismatch")
            if int(parts[indexes["Length"]]) != correction["members"][accession]:
                raise ValueError("Correction reference length mismatch")
            if (
                int(parts[indexes["Asm_length"]])
                != correction["original_assembly_length_nt"]
            ):
                raise ValueError("Correction original assembly length mismatch")
            parts[indexes["Asm_length"]] = str(
                correction["corrected_assembly_length_nt"]
            ).encode()
            line = b"\t".join(parts) + ending
            seen.add(accession)
        result.append(line)
    if seen != set(lookup):
        raise ValueError("Reference correction member set is incomplete")
    return b"".join(result)


def corrected_catalog_bytes(source, specification):
    """Apply the eight verified annotations without rewriting unrelated values."""
    digest = hashlib.sha256(source).hexdigest()
    if digest == specification.get("corrected_catalog_sha256"):
        return source
    if digest != specification["original_catalog_sha256"]:
        raise ValueError("Reference update does not match the bundled base catalog")
    stream = io.StringIO(source.decode())
    reader = csv.DictReader(stream, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)
    lookup = {a: c for c in specification["corrections"] for a in c["members"]}
    seen = set()
    for row in rows:
        accession = row["accession"]
        if accession not in lookup:
            continue
        correction = lookup[accession]
        if (
            row["assembly"] != correction["assembly"]
            or int(row["length_nt"]) != correction["members"][accession]
        ):
            raise ValueError("Catalog correction identity mismatch")
        if int(row["assembly_length_nt"]) != correction["original_assembly_length_nt"]:
            raise ValueError("Catalog correction original length mismatch")
        factor = float(row["molecular_mass_per_nt_g_mol"])
        if row["molecule_type"] != "ssRNA" or factor != 340:
            raise ValueError(
                "Catalog correction chemistry differs from verified evidence"
            )
        row["assembly_length_nt"] = str(correction["corrected_assembly_length_nt"])
        row["assembly_reference_count"] = str(len(correction["members"]))
        row["assembly_molecular_mass_g_mol"] = format(
            factor * correction["corrected_assembly_length_nt"], ".12g"
        )
        row["assembly_mass_status"] = "resolved_complete_assembly"
        note = "reference_correction:" + correction["correction_id"]
        row["notes"] = note if row["notes"] in {"", "NA"} else row["notes"] + ";" + note
        seen.add(accession)
    if seen != set(lookup):
        raise ValueError("Catalog correction member set is incomplete")
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode()


def bundled_reference_files():
    """Return checksum-verified revised catalog bytes and its matching manifest."""
    specification = read_corrections()
    with gzip.open(PACKAGE_REFERENCE / (CATALOG_NAME + ".gz"), "rb") as stream:
        source = stream.read()
    result = corrected_catalog_bytes(source, specification)
    if hashlib.sha256(result).hexdigest() != specification["corrected_catalog_sha256"]:
        raise ValueError("Revised bundled catalog checksum mismatch")
    manifest = specification["corrected_catalog_manifest"]
    if (
        manifest["catalog_sha256"] != specification["corrected_catalog_sha256"]
        or manifest["database_metadata_sha256"]
        != specification["corrected_metadata_sha256"]
    ):
        raise ValueError(
            "Revised reference manifest does not match correction checksums"
        )
    return result, manifest


def applicable_correction(accession, annotation):
    """Authorize only a documented correction present in the selected catalog."""
    if not annotation:
        return None
    for correction in read_corrections()["corrections"]:
        if (
            accession in correction["members"]
            and annotation["assembly"] == correction["assembly"]
            and float(annotation["length_nt"]) == correction["members"][accession]
            and float(annotation["assembly_length_nt"])
            == correction["corrected_assembly_length_nt"]
            and annotation["assembly_mass_status"] == "resolved_complete_assembly"
            and "reference_correction:" + correction["correction_id"]
            in annotation["notes"]
        ):
            return correction
    return None


def verify_fasta(directory, specification):
    fasta = Path(directory) / FASTA_NAME
    if not fasta.is_file():
        raise ValueError(f"Database FASTA is required for a verified update: {fasta}")
    expected = {
        a: n for c in specification["corrections"] for a, n in c["members"].items()
    }
    lengths = {}
    accession, length = None, 0
    digest = hashlib.sha256()
    with fasta.open("rb") as stream:
        for line in stream:
            digest.update(line)
            if line.startswith(b">"):
                if accession in expected:
                    if accession in lengths:
                        raise ValueError("Duplicate correction accession in FASTA")
                    lengths[accession] = length
                accession = line[1:].split()[0].decode()
                length = 0
            else:
                length += len(line.strip())
        if accession in expected:
            if accession in lengths:
                raise ValueError("Duplicate correction accession in FASTA")
            lengths[accession] = length
    if digest.hexdigest() != specification["fasta_sha256"] or lengths != expected:
        raise ValueError(
            "FASTA checksum or correction member lengths differ from verified evidence"
        )
    return {"sha256": digest.hexdigest(), "verified_member_lengths": lengths}


def upgrade_database(directory, check_only=False):
    """Install v3.3.1 with staged validation, backups, rollback and no sequence edits."""
    directory = Path(directory).resolve()
    specification = read_corrections()
    metadata = directory / METADATA_NAME
    if not metadata.is_file():
        raise ValueError(f"Database metadata is required: {metadata}")
    original = metadata.read_bytes()
    revised = corrected_metadata_bytes(original, specification)
    if (
        hashlib.sha256(revised).hexdigest()
        != specification["corrected_metadata_sha256"]
    ):
        raise ValueError(
            "Revised metadata checksum differs from the published correction"
        )
    fasta_evidence = verify_fasta(directory, specification)
    catalog = directory / CATALOG_NAME
    if catalog.exists() and sha256(catalog) not in {
        specification["original_catalog_sha256"],
        specification["corrected_catalog_sha256"],
    }:
        raise ValueError(
            "Existing catalog contains unrecognized changes; no files were replaced"
        )
    report = {
        "database_revision": REVISION,
        "software_version": __version__,
        "status": "update_available",
        "database_directory": str(directory),
        "original_metadata_sha256": hashlib.sha256(original).hexdigest(),
        "corrected_metadata_sha256": specification["corrected_metadata_sha256"],
        "corrected_catalog_sha256": specification["corrected_catalog_sha256"],
        "fasta_evidence": fasta_evidence,
        "alignment_sequences_changed": False,
        "alignment_index_rebuild_required": False,
        "correction_ids": [c["correction_id"] for c in specification["corrections"]],
    }
    if (
        original == revised
        and catalog.is_file()
        and (directory / MANIFEST_NAME).is_file()
    ):
        records, context = database_catalog(directory)
        if context["sha256"] == specification["corrected_catalog_sha256"]:
            return {
                **report,
                "status": "already_current",
                "reference_count": len(records),
            }
    if check_only:
        return {**report, "status": "check_only"}
    lock = directory / ".hybcapseq_quant_upgrade.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise ValueError("Another database upgrade holds the update lock") from exc
    os.close(descriptor)
    try:
        revised_catalog, manifest = bundled_reference_files()
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        backup = directory / "hybcapseq_quant_backups" / (REVISION + "_" + stamp)
        report.update(
            status="updated", completed_at_utc=stamp, backup_directory=str(backup)
        )
        revision_manifest_name = f"esviritu_db_{REVISION}_manifest.json"
        revision_checksums_name = f"esviritu_db_{REVISION}.sha256"
        report_name = "virus_pathogen_database.upgrade_report.json"
        with tempfile.TemporaryDirectory(
            prefix=".reference_update_", dir=directory
        ) as temp:
            staging = Path(temp)
            (staging / METADATA_NAME).write_bytes(revised)
            (staging / CATALOG_NAME).write_bytes(revised_catalog)
            (staging / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
            (staging / CORRECTIONS_NAME).write_text(
                json.dumps(specification, indent=2) + "\n"
            )
            load_catalog(staging / CATALOG_NAME, staging / METADATA_NAME)
            index = directory / INDEX_NAME
            files = [
                {"name": FASTA_NAME, "sha256": fasta_evidence["sha256"]},
                {"name": METADATA_NAME, "sha256": sha256(staging / METADATA_NAME)},
                {"name": CATALOG_NAME, "sha256": sha256(staging / CATALOG_NAME)},
                {"name": MANIFEST_NAME, "sha256": sha256(staging / MANIFEST_NAME)},
                {
                    "name": CORRECTIONS_NAME,
                    "sha256": sha256(staging / CORRECTIONS_NAME),
                },
            ]
            if index.is_file():
                files.append({"name": INDEX_NAME, "sha256": sha256(index)})
            revision_manifest = {
                "database_version": REVISION,
                "source_database_version": "v3.3",
                "reference_count": manifest["reference_count"],
                "reference_corrections": specification["corrections"],
                "alignment_sequences_changed": False,
                "files": files,
            }
            (staging / revision_manifest_name).write_text(
                json.dumps(revision_manifest, indent=2) + "\n"
            )
            (staging / revision_checksums_name).write_text(
                "".join(f"{f['sha256']}  {f['name']}\n" for f in files)
            )
            (staging / report_name).write_text(json.dumps(report, indent=2) + "\n")
            names = [
                METADATA_NAME,
                CATALOG_NAME,
                MANIFEST_NAME,
                CORRECTIONS_NAME,
                revision_manifest_name,
                revision_checksums_name,
                report_name,
            ]
            if metadata.read_bytes() != original:
                raise ValueError("Database metadata changed during update preparation")
            backup.mkdir(parents=True)
            existed = {name: (directory / name).exists() for name in names}
            for name in names:
                target = directory / name
                if existed[name]:
                    shutil.copy2(target, backup / name)
                    shutil.copymode(target, staging / name)
            replaced = []
            try:
                for name in names:
                    os.replace(staging / name, directory / name)
                    replaced.append(name)
                records, context = database_catalog(directory)
            except Exception:
                for name in reversed(replaced):
                    if existed[name]:
                        shutil.copy2(backup / name, directory / name)
                    else:
                        (directory / name).unlink()
                raise
        return {
            **report,
            "reference_count": len(records),
            "catalog_sha256": context["sha256"],
        }
    finally:
        lock.unlink(missing_ok=True)
