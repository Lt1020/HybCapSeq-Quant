"""Accession identity, evidence precedence and molecular mass verification."""

import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from hybcapseq_quant.cli import load_targets, main
from hybcapseq_quant.core import read_tsv, reference_groups, sha256, write_tsv
from hybcapseq_quant.references import (
    CATALOG_NAME,
    MANIFEST_NAME,
    METADATA_NAME,
    PACKAGE_REFERENCE,
    MASS_FACTORS,
    MASS_MODEL,
    annotate_reference,
    load_catalog,
    install_catalog,
)


class ReferenceAnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        with gzip.open(PACKAGE_REFERENCE / (CATALOG_NAME + ".gz"), "rb") as src:
            (cls.root / CATALOG_NAME).write_bytes(src.read())
        shutil.copyfile(PACKAGE_REFERENCE / MANIFEST_NAME, cls.root / MANIFEST_NAME)
        cls.records, cls.context = load_catalog(cls.root / CATALOG_NAME)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_complete_accession_catalog_and_declared_models(self):
        self.assertEqual(len(self.records), 31003)
        resolved = [
            r
            for r in self.records.values()
            if r["reference_molecular_mass_g_mol"] != "NA"
        ]
        self.assertEqual(len(resolved), 30729)
        self.assertEqual(sum(a.startswith("ERCC-") for a in self.records), 96)
        self.assertEqual(set(r["molecule_type"] for r in resolved), set(MASS_FACTORS))
        for row in resolved:
            self.assertEqual(
                float(row["reference_molecular_mass_g_mol"]),
                float(row["length_nt"]) * MASS_FACTORS[row["molecule_type"]],
            )
        self.assertEqual(self.records["ERCC-00007"]["molecule_type"], "ssRNA")

    def test_prefixed_taxa_and_ictv_genome_priority_preserve_ncbi_chemistry(self):
        row = {
            "Accession": "TEST_1.1",
            "Length": "100",
            "Asm_length": "100",
            "Assembly": "TEST",
            "family": "f__Exampleviridae",
        }
        ncbi = {
            "retrieval_status": "exact_accession_version",
            "slen": 100,
            "moltype": "dna",
            "strand": "ds",
        }
        ictv = {
            ("family", "Exampleviridae"): [
                {"chemistry": "ssRNA", "genome": "ssRNA(+)", "isolate_id": "TEST_ONLY"}
            ]
        }
        result = annotate_reference(row, ncbi, ictv, set())
        self.assertEqual(result["molecule_type"], "ssRNA")
        self.assertEqual(result["ncbi_molecule_type"], "dna")
        self.assertEqual(
            result["annotation_status"], "resolved_ictv_genome_record_chemistry_differs"
        )
        self.assertEqual(result["family_raw"], "f__Exampleviridae")
        self.assertEqual(result["reference_molecular_mass_g_mol"], 34000)
        ncbi["slen"] = 99
        result = annotate_reference(row, ncbi, ictv, set())
        self.assertIsNone(result["reference_molecular_mass_g_mol"])
        self.assertEqual(result["annotation_status"], "ncbi_database_length_mismatch")

    def test_missing_strandedness_never_defaults_to_rna(self):
        row = {"Accession": "TEST_1.1", "Length": "100", "Asm_length": "100"}
        result = annotate_reference(
            row,
            {
                "retrieval_status": "exact_accession_version",
                "slen": 100,
                "moltype": "rna",
            },
            {},
            set(),
        )
        self.assertIsNone(result["reference_molecular_mass_g_mol"])
        self.assertEqual(result["genome_nucleic_acid_type"], "RNA")

    def test_target_loading_uses_exact_accession_and_checks_length(self):
        known = next(
            r
            for r in self.records.values()
            if r["molecule_type"] == "dsRNA"
            and r["assembly_mass_status"] == "resolved_complete_assembly"
        )
        row = {
            "Accession": known["accession"],
            "Length": known["length_nt"],
            "Assembly": known["assembly"],
            "Asm_length": known["assembly_length_nt"],
        }
        target = load_targets(None, [{"rows": [row]}], self.records)[known["accession"]]
        self.assertEqual(float(target["molecular_mass_per_nt_g_mol"]), 680)
        row["Length"] = "1"
        with self.assertRaisesRegex(ValueError, "Length mismatch"):
            load_targets(None, [{"rows": [row]}], self.records)

    def test_reviewed_factor_cannot_conflict_with_inherited_chemistry(self):
        known = next(
            r
            for r in self.records.values()
            if r["molecule_type"] == "ssRNA"
            and not r["accession"].startswith("ERCC-")
            and r["assembly_mass_status"] == "resolved_complete_assembly"
        )
        row = {
            "Accession": known["accession"],
            "Length": known["length_nt"],
            "Assembly": known["assembly"],
            "Asm_length": known["assembly_length_nt"],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "targets.tsv"
            write_tsv(
                path,
                [
                    {
                        "accession": known["accession"],
                        "molecular_mass_per_nt_g_mol": 660,
                        "molecular_mass_source": "SIMULATED conflicting declaration",
                    }
                ],
            )
            with self.assertRaisesRegex(ValueError, "type/factor mismatch"):
                load_targets(path, [{"rows": [row]}], self.records)

    def test_checksum_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shutil.copyfile(self.root / MANIFEST_NAME, root / MANIFEST_NAME)
            (root / CATALOG_NAME).write_text("accession\nFORGED\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                load_catalog(root / CATALOG_NAME)

    def test_install_rejects_unrelated_database_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / METADATA_NAME).write_text("Accession\nTEST\n")
            with self.assertRaisesRegex(ValueError, "checksum"):
                install_catalog(directory)

    def test_unknown_segment_prevents_group_mass_assignment(self):
        rows = [
            {"Accession": a, "Assembly": "GROUP", "Length": 100, "Asm_length": 200}
            for a in ["TEST_A", "TEST_B"]
        ]
        targets = {
            a: {
                "reference_length_nt": 200,
                "molecular_mass_per_nt_g_mol": 680 if a == "TEST_A" else "NA",
            }
            for a in ["TEST_A", "TEST_B"]
        }
        data = {
            "rows": rows,
            "signals": {a: {"aligned_bases": 50, "covered_bases": 50} for a in targets},
        }
        groups = reference_groups(data, targets)
        self.assertIsNone(groups["TEST_A"]["molecular_mass_per_nt_g_mol"])
        self.assertIsNone(groups["TEST_B"]["quantification_molecular_mass_g_mol"])

    def test_unreviewed_taxonomic_labels_do_not_assign_mass(self):
        row = {
            "Accession": "TEST_1",
            "Assembly": "TEST",
            "Length": 100,
            "Asm_length": 100,
            "kingdom": "k__Orthornavirae",
            "family": "f__Flaviviridae",
        }
        result = load_targets(None, [{"rows": [row]}])["TEST_1"]
        self.assertEqual(result["molecular_mass_per_nt_g_mol"], "NA")


if __name__ == "__main__":
    unittest.main()
