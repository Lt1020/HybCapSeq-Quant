"""Reference preflight uses only nonpathogen synthetic byte fixtures."""

import csv
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "preflight", Path(__file__).parents[1] / "scripts/database/preflight_v3_4.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PreflightTests(unittest.TestCase):
    def test_empty_manifest_is_blocked(self):
        path = (
            Path(__file__).parents[1] / "databases/v3.4/target_references.template.tsv"
        )
        result = MODULE.audit(path)
        self.assertFalse(result["preflight_passed"])
        self.assertFalse(result["release_ready"])
        self.assertIn("no_reviewed_target_records_supplied", result["blocking_issues"])

    def test_input_integrity_does_not_claim_database_release(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            content = b">SIMULATED_REFERENCE\n" + b"ACGT" * 75 + b"\n"
            (root / "supplied.fna").write_bytes(content)
            for name in [
                "virus_pathogen_database.fna",
                "virus_pathogen_database.all_metadata.tsv",
                "virus_pathogen_database.molecular_annotations.tsv",
                "virus_pathogen_database.molecular_annotations.manifest.json",
            ]:
                (root / name).touch()
            row = dict.fromkeys(MODULE.FIELDS, "")
            row.update(
                target="SIMULATED_TARGET",
                reference_id="SIMULATED_REFERENCE",
                source_accession_version="SIMULATED.1",
                source_url="https://example.invalid/synthetic_fixture",
                assay_region="synthetic",
                actual_assay_confirmed="True",
                reviewed="True",
                redistribution_allowed="True",
                fasta_path="supplied.fna",
                fasta_sha256=hashlib.sha256(content).hexdigest(),
                reference_start_1based=1,
                reference_end_1based=300,
                amplicon_start_1based=101,
                amplicon_end_1based=200,
            )
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=MODULE.FIELDS, delimiter="\t"
                )
                writer.writeheader()
                writer.writerow(row)
            result = MODULE.audit(manifest, root)
            self.assertTrue(result["preflight_passed"])
            self.assertFalse(result["release_ready"])
