"""Reference integrity checks against the original manufacturer table."""

import csv
import hashlib
import json
import math
from pathlib import Path
import unittest

from hybcapseq_quant.core import NA, STANDARD_PATH, read_tsv


class StandardsTests(unittest.TestCase):
    def test_manufacturer_values_and_checksum(self):
        folder = STANDARD_PATH.parent
        provenance = json.loads((folder / "provenance.json").read_text())
        original = folder / "ERCC_Controls_Analysis.txt"
        self.assertEqual(
            hashlib.sha256(original.read_bytes()).hexdigest(),
            provenance["concentration_source_sha256"],
        )
        with original.open(newline="") as src:
            rows = list(csv.DictReader(src, delimiter="\t"))
        nominal = {
            r["ERCC ID"]: float(r["concentration in Mix 1 (attomoles/ul)"])
            for r in rows
        }
        _, standards = read_tsv(STANDARD_PATH)
        self.assertEqual(set(nominal), {r["ercc_id"] for r in standards})
        for standard in standards:
            expected = nominal[standard["ercc_id"]]
            self.assertEqual(float(standard["stock_attomol_ul"]), expected)
            self.assertTrue(
                math.isclose(
                    float(standard["stock_copies_ul"]),
                    expected * NA * 1e-18,
                    rel_tol=1e-11,
                )
            )

    def test_lengths_match_combined_database_metadata(self):
        path = (
            Path(__file__).parents[1]
            / "references/ercc96/data/certified_v1/virus_pathogen_database.all_metadata.tsv"
        )
        _, metadata = read_tsv(path)
        lengths = {r["Name"]: int(r["Length"]) for r in metadata}
        _, standards = read_tsv(STANDARD_PATH)
        for standard in standards:
            self.assertEqual(int(standard["length_nt"]), lengths[standard["ercc_id"]])

    def test_packaged_and_public_templates_are_identical(self):
        root = Path(__file__).parents[1]
        for name in [
            "samples_template.tsv",
            "targets_template.tsv",
            "input_fields.tsv",
        ]:
            self.assertEqual(
                (STANDARD_PATH.parent / name).read_bytes(),
                (root / "templates" / name).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
