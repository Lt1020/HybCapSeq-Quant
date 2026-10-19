"""Database migration and legacy-result checks using synthetic byte fixtures."""

import copy
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from hybcapseq_quant.cli import load_targets, parser, settings_from_args
from hybcapseq_quant.core import read_tsv, reference_groups, sha256
from hybcapseq_quant.reference_updates import (
    CORRECTIONS_NAME,
    corrected_catalog_bytes,
    corrected_metadata_bytes,
    read_corrections,
    upgrade_database,
)
from hybcapseq_quant.references import (
    CATALOG_NAME,
    MANIFEST_NAME,
    METADATA_NAME,
    PACKAGE_REFERENCE,
    load_catalog,
)


class ReferenceUpdateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.published = read_corrections()
        wanted = cls.published["corrections"][0]["members"]
        with gzip.open(PACKAGE_REFERENCE / (CATALOG_NAME + ".gz"), "rt") as source:
            cls.base_rows = [
                r
                for r in csv.DictReader(source, delimiter="\t")
                if r["accession"] in wanted
            ]

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.db = self.root / "database"
        self.db.mkdir()
        self.package = self.root / "package"
        self.package.mkdir()
        spec = copy.deepcopy(self.published)
        correction = spec["corrections"][0]
        # Asm_length is deliberately the last field, as in the real source file.
        header = "Accession\tName\tAssembly\tLength\tAsm_length\r\n"
        body = "".join(
            f"{a}\tTEST_FIXTURE\t{correction['assembly']}\t{n}\t15617\r\n"
            for a, n in correction["members"].items()
        )
        metadata = (header + body).encode()
        (self.db / METADATA_NAME).write_bytes(metadata)
        spec["original_metadata_sha256"] = hashlib.sha256(metadata).hexdigest()
        # Sequence contents are synthetic padding for file-integrity tests only.
        fasta = b"".join(
            b">" + a.encode() + b"\n" + b"A" * n + b"\n"
            for a, n in correction["members"].items()
        )
        (self.db / "virus_pathogen_database.fna").write_bytes(fasta)
        (self.db / "virus_pathogen_database.mmi").write_bytes(b"TEST_INDEX_UNCHANGED")
        spec["fasta_sha256"] = hashlib.sha256(fasta).hexdigest()
        out = io.StringIO(newline="")
        writer = csv.DictWriter(
            out, fieldnames=list(self.base_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(self.base_rows)
        original_catalog = out.getvalue().encode()
        spec["original_catalog_sha256"] = hashlib.sha256(original_catalog).hexdigest()
        (self.db / CATALOG_NAME).write_bytes(original_catalog)
        with gzip.open(self.package / (CATALOG_NAME + ".gz"), "wb") as stream:
            stream.write(original_catalog)
        base_manifest = {
            "schema_version": "1.0",
            "reference_count": 8,
            "catalog_sha256": spec["original_catalog_sha256"],
            "database_metadata_sha256": spec["original_metadata_sha256"],
            "mass_model": "average_reference_position_v1",
        }
        (self.db / MANIFEST_NAME).write_text(json.dumps(base_manifest))
        revised_metadata = corrected_metadata_bytes(metadata, spec)
        revised_catalog = corrected_catalog_bytes(original_catalog, spec)
        spec["corrected_metadata_sha256"] = hashlib.sha256(revised_metadata).hexdigest()
        spec["corrected_catalog_sha256"] = hashlib.sha256(revised_catalog).hexdigest()
        spec["corrected_catalog_manifest"] = {
            **base_manifest,
            "database_revision": "v3.3.1",
            "catalog_sha256": spec["corrected_catalog_sha256"],
            "database_metadata_sha256": spec["corrected_metadata_sha256"],
            "reference_correction_ids": [correction["correction_id"]],
        }
        (self.package / CORRECTIONS_NAME).write_text(json.dumps(spec))
        self.spec = spec
        self.patch = mock.patch(
            "hybcapseq_quant.reference_updates.PACKAGE_REFERENCE", self.package
        )
        self.patch.start()

    def tearDown(self):
        self.patch.stop()
        self.temp.cleanup()

    def fingerprints(self):
        return {
            str(p.relative_to(self.db)): sha256(p)
            for p in self.db.rglob("*")
            if p.is_file()
        }

    def test_upgrade_preserves_sequences_index_and_source_backups(self):
        before = self.fingerprints()
        result = upgrade_database(self.db)
        self.assertEqual(result["status"], "updated")
        records, context = load_catalog(self.db / CATALOG_NAME, self.db / METADATA_NAME)
        self.assertEqual(context["database_revision"], "v3.3.1")
        self.assertTrue(
            all(
                r["assembly_mass_status"] == "resolved_complete_assembly"
                for r in records.values()
            )
        )
        self.assertTrue(
            all(
                float(r["assembly_molecular_mass_g_mol"]) == 4712060
                for r in records.values()
            )
        )
        for name in ["virus_pathogen_database.fna", "virus_pathogen_database.mmi"]:
            self.assertEqual(sha256(self.db / name), before[name])
        for name in [METADATA_NAME, CATALOG_NAME, MANIFEST_NAME]:
            self.assertEqual(
                sha256(Path(result["backup_directory"]) / name), before[name]
            )
        fields, rows = read_tsv(self.db / METADATA_NAME)
        self.assertEqual(len(rows), 8)
        self.assertEqual(fields[-1], "Asm_length")
        self.assertTrue(all(r["Asm_length"] == "13859" for r in rows))
        self.assertEqual((self.db / METADATA_NAME).read_bytes().count(b"\r\n"), 9)
        snapshot = self.fingerprints()
        self.assertEqual(upgrade_database(self.db)["status"], "already_current")
        self.assertEqual(self.fingerprints(), snapshot)

    def test_check_only_performs_no_writes(self):
        before = self.fingerprints()
        self.assertEqual(
            upgrade_database(self.db, check_only=True)["status"], "check_only"
        )
        self.assertEqual(self.fingerprints(), before)

    def test_unknown_metadata_or_fasta_is_not_repaired(self):
        metadata = self.db / METADATA_NAME
        original = metadata.read_bytes()
        metadata.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            upgrade_database(self.db)
        self.assertEqual(metadata.read_bytes(), original + b"\n")
        metadata.write_bytes(original)
        fasta = self.db / "virus_pathogen_database.fna"
        fasta.write_bytes(fasta.read_bytes().replace(b"AAAA", b"CCCC", 1))
        with self.assertRaisesRegex(ValueError, "FASTA"):
            upgrade_database(self.db)
        self.assertEqual(metadata.read_bytes(), original)
        self.assertFalse((self.db / "hybcapseq_quant_backups").exists())

    def test_unrecognized_catalog_is_not_overwritten(self):
        catalog = self.db / CATALOG_NAME
        catalog.write_bytes(catalog.read_bytes() + b"\n")
        before = self.fingerprints()
        with self.assertRaisesRegex(ValueError, "unrecognized changes"):
            upgrade_database(self.db)
        self.assertEqual(before, self.fingerprints())

    def test_failed_replacement_restores_original_files(self):
        before = self.fingerprints()
        replace = os.replace
        calls = 0

        def fail_second(source, destination):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError("TEST simulated write failure")
            return replace(source, destination)

        with mock.patch(
            "hybcapseq_quant.reference_updates.os.replace", side_effect=fail_second
        ):
            with self.assertRaisesRegex(OSError, "simulated"):
                upgrade_database(self.db)
        for name, digest in before.items():
            self.assertEqual(sha256(self.db / name), digest)
        self.assertFalse((self.db / ".hybcapseq_quant_upgrade.lock").exists())

    def corrected_records(self):
        upgrade_database(self.db)
        return load_catalog(self.db / CATALOG_NAME, self.db / METADATA_NAME)[0]

    def test_legacy_annotation_is_preserved_and_effective_reference_is_corrected(self):
        catalog = self.corrected_records()
        rows = [
            {
                "Accession": a,
                "Assembly": r["assembly"],
                "Length": r["length_nt"],
                "Asm_length": "15617",
                "sample_ID": "TEST",
            }
            for a, r in catalog.items()
        ]
        data = {
            "rows": rows,
            "sample": {"sample_id": "TEST"},
            "signals": {
                r["Accession"]: {"aligned_bases": 200, "covered_bases": 100}
                for r in rows
            },
        }
        targets = load_targets(None, [data], catalog)
        group = next(iter(reference_groups(data, targets).values()))
        self.assertEqual(group["reference_length_nt"], 13859)
        self.assertEqual(group["detected_reference_length_nt"], 13859)
        self.assertEqual(group["mean_depth"], 1600 / 13859)
        self.assertTrue(group["reference_length_valid"])
        self.assertEqual(group["molecular_mass_per_nt_g_mol"], 340)
        self.assertTrue(all(r["Asm_length"] == "15617" for r in rows))
        self.assertEqual(len(data["reference_corrections"]), 8)
        self.assertTrue(
            all(r["input_length_was_corrected"] for r in data["reference_corrections"])
        )

    def test_other_length_mismatches_are_not_silently_corrected(self):
        catalog = self.corrected_records()
        a, annotation = next(iter(catalog.items()))
        row = {
            "Accession": a,
            "Assembly": annotation["assembly"],
            "Length": annotation["length_nt"],
            "Asm_length": "15618",
        }
        with self.assertRaisesRegex(ValueError, "Asm_length mismatch"):
            load_targets(None, [{"rows": [row]}], catalog)
        row["Asm_length"] = "15617"
        row["Length"] = str(int(row["Length"]) + 1)
        with self.assertRaisesRegex(ValueError, "Length mismatch"):
            load_targets(None, [{"rows": [row]}], catalog)

    def test_catalog_rejects_semantic_assembly_mass_corruption(self):
        self.corrected_records()
        path = self.db / CATALOG_NAME
        path.write_text(path.read_text().replace("4712060", "4712061"))
        manifest_path = self.db / MANIFEST_NAME
        manifest = json.loads(manifest_path.read_text())
        manifest["catalog_sha256"] = sha256(path)
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "Assembly molecular mass"):
            load_catalog(path, self.db / METADATA_NAME)

    def test_cli_reporting_mode_resolves_extrapolation_defaults(self):
        for mode, expected in [("exploratory", True), ("validated", False)]:
            argv = ["run", "--samples", "TEST.tsv", "--reporting-mode", mode]
            self.assertEqual(
                settings_from_args(parser().parse_args(argv)).allow_extrapolation,
                expected,
            )
            self.assertFalse(
                settings_from_args(
                    parser().parse_args(argv + ["--no-extrapolation"])
                ).allow_extrapolation
            )
            self.assertTrue(
                settings_from_args(
                    parser().parse_args(argv + ["--allow-extrapolation"])
                ).allow_extrapolation
            )


if __name__ == "__main__":
    unittest.main()
