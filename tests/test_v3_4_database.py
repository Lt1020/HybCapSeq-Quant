"""Contracts for the frozen literature extension and conservative local mass model."""

import csv
import hashlib
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "references/pcr_roi"


class LiteratureExtensionTests(unittest.TestCase):
    def test_primary_scope_and_four_sars_regions(self):
        with (PACKAGE / "virus_pathogen_database.all_metadata.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 21)
        sars = [
            r["Accession"]
            for r in rows
            if r["Accession"].startswith("PCRROI_SARS-CoV-2_")
        ]
        self.assertEqual(len(sars), 4)
        for name in ["CDC2020_N1", "CDC2020_N2", "Corman2020_E", "ChinaCDC2020_ORF1ab"]:
            self.assertTrue(any(name in value for value in sars))
        self.assertTrue(
            all(
                r["Length"] == r["Asm_length"] and r["Assembly"] == r["Name"]
                for r in rows
            )
        )

    def test_frozen_inputs_and_exceptions(self):
        for line in (PACKAGE / "SHA256SUMS").read_text().splitlines():
            expected, name = line.split("  ", 1)
            self.assertEqual(
                hashlib.sha256((PACKAGE / name).read_bytes()).hexdigest(), expected
            )
        with (PACKAGE / "reference_regions.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        primary = {
            r["reference_id"]
            for r in rows
            if r["sequence_check_status"] == "full100_compatible_sites"
        }
        self.assertEqual(len(primary), 21)
        self.assertFalse(any("Rotavirus_A" in value for value in primary))
        for row in rows:
            if row["reference_id"] in primary:
                self.assertEqual(int(row["left_flank_nt"]), 100)
                self.assertEqual(int(row["right_flank_nt"]), 100)


if __name__ == "__main__":
    unittest.main()
