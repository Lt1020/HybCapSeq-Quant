"""Analytical checks against independently calculated quantities."""

import contextlib
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from hybcapseq_quant.cli import main, load_targets
from hybcapseq_quant.core import (
    CONC,
    WWCONC,
    METHODS,
    NA,
    STANDARD_PATH,
    Settings,
    calibrate,
    ercc_amounts,
    load_sample,
    quantify_rows,
    read_tsv,
    reference_groups,
    replicate_qc,
    wastewater_factor,
    write_tsv,
)

SPEC = importlib.util.spec_from_file_location(
    "analytical_fixture", Path(__file__).parents[1] / "examples/quant/generate_demo.py"
)
DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEMO)


class QuantificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.sample, self.rows = DEMO.create(self.root)
        _, self.standards = read_tsv(STANDARD_PATH)
        self.path = self.root / "SIMULATED_01.detected_virus.info.tsv"
        self.settings = Settings()

    def tearDown(self):
        self.temp.cleanup()

    def load(self, sample=None):
        return load_sample(
            sample or self.sample, self.path, self.standards, self.root, self.settings
        )

    def quantify(self, data=None, method=METHODS[0], settings=None, targets=None):
        data = data or self.load()
        settings = settings or self.settings
        cal, _, _, _ = calibrate(data, method, settings)
        return quantify_rows(
            data,
            cal,
            targets or load_targets(self.root / "targets.tsv", [data]),
            settings,
        )

    def target(self, rows):
        return next(r for r in rows if r["Accession"] == "SIM_RNA")

    def test_official_transcript_set_and_levels(self):
        self.assertEqual(len(self.standards), 92)
        self.assertEqual(len({r["level"] for r in self.standards}), 22)
        self.assertFalse(
            {r["ercc_id"] for r in self.standards}
            & {"ERCC-00007", "ERCC-00018", "ERCC-00023", "ERCC-00128"}
        )

    def test_nominal_attomole_conversion(self):
        amounts = ercc_amounts(self.standards, self.sample)
        expected = 30000 * 1e-18 * 6.02214076e23 / 1000 * 2
        self.assertAlmostEqual(amounts["ERCC-00130"]["added_copies"], expected)
        self.assertAlmostEqual(
            amounts["ERCC-00130"]["reaction_copies_ul"], expected / 14
        )
        self.assertAlmostEqual(
            amounts["ERCC-00130"]["extract_equivalent_copies_ul"], expected / 2
        )

    def test_both_models_recover_known_molecule_count(self):
        for method in METHODS:
            rows, targets, _ = self.quantify(method=method)
            row = self.target(rows)
            self.assertAlmostEqual(row[CONC], 10000 / 14, places=7)
            self.assertAlmostEqual(
                row["concentration_extract_copies_ul"], 5000, places=7
            )
            self.assertAlmostEqual(
                row[WWCONC], 5000 * 1000 * 80 / 300 * 0.46 / (40 * 0.25), places=6
            )
            self.assertTrue(row["target_validation_verified"])
            self.assertEqual(len(targets), 2)
            for segment in [r for r in rows if r["Accession"].startswith("SIM_SEG_")]:
                self.assertAlmostEqual(segment[CONC], 10000 / 14, places=7)

    def test_yield_uses_arithmetic_mean_not_regression_inversion(self):
        data = self.load()
        data["points"][0]["base_fraction"] *= 1.15
        data["points"][0]["sequencing_yield"] *= 1.15
        calibration, eligible, _, _ = calibrate(data, METHODS[0], self.settings)
        mean_yield = sum(p["sequencing_yield"] for p in eligible) / len(eligible)
        self.assertAlmostEqual(calibration["mean_sequencing_yield"], mean_yield)
        results, _, _ = quantify_rows(
            data,
            calibration,
            load_targets(self.root / "targets.tsv", [data]),
            self.settings,
        )
        signal = data["signals"]["SIM_RNA"]["aligned_bases"]
        expected_mass = data["total_mass"] * signal / data["total_bases"] / mean_yield
        expected_c = expected_mass * 1e-9 * NA / (340 * 10000 * 14)
        self.assertAlmostEqual(self.target(results)[CONC], expected_c, places=7)

    def test_undetected_standard_retains_assigned_input_without_observed_zero(self):
        missing = self.rows[0]["Accession"]
        write_tsv(self.path, self.rows[1:])
        data = self.load()
        point = next(p for p in data["points"] if p["ercc_id"] == missing)
        self.assertFalse(point["detected"])
        self.assertFalse(point["fit_used"])
        self.assertIsNone(point["mean_depth"])
        self.assertGreater(point["added_mass_ng"], 0)
        self.assertEqual(len(data["points"]), 92)

    def test_total_mass_includes_all_spiked_standards(self):
        data = self.load()
        expected = 120 + sum(x["added_mass_ng"] for x in data["amounts"].values())
        self.assertAlmostEqual(data["total_mass"], expected)
        data["points"][0]["fit_used"] = False
        calibration, _, _, _ = calibrate(data, METHODS[0], self.settings)
        self.assertTrue(calibration["passes"])
        self.assertAlmostEqual(data["total_mass"], expected)

    def test_mass_concentration_volume_identity(self):
        with self.assertRaisesRegex(ValueError, "library_input_mass_ng"):
            self.load({**self.sample, "library_input_mass_ng": 100})

    def test_known_standards_have_no_wastewater_concentration(self):
        rows, _, _ = self.quantify()
        standard = next(r for r in rows if r["Accession"] == "ERCC-00130")
        self.assertAlmostEqual(standard[CONC], 2580917.4685714287)
        self.assertEqual(standard["value_origin"], "assigned_standard")
        self.assertIsNone(standard[WWCONC])

    def test_missing_recovery_preserves_reaction_but_not_wastewater(self):
        rows, _, audit = self.quantify(
            self.load({**self.sample, "pmmov_recovery_fraction": "NA"})
        )
        self.assertIsNotNone(self.target(rows)[CONC])
        self.assertIsNone(self.target(rows)[WWCONC])
        self.assertIn("missing", audit["status"])

    def test_liter_conversion(self):
        ml, _ = wastewater_factor(self.sample, "copies/mL")
        liter, _ = wastewater_factor(self.sample, "copies/L")
        self.assertEqual(liter, ml * 1000)

    def test_pmmov_mass_balance_and_dilution(self):
        sample = {
            **self.sample,
            "pmmov_recovery_fraction": "NA",
            "pmmov_raw_extract_copies_ul": 10,
            "pmmov_post_extract_copies_ul": 1000,
            "pmmov_raw_dilution": 1,
            "pmmov_post_dilution": 1,
            "raw_extract_elution_volume_ul": 80,
            "raw_extraction_input_volume_ul": 300,
            "pmmov_above_loq": 1,
        }
        _, audit = wastewater_factor(sample, "copies/mL")
        self.assertAlmostEqual(audit["pmmov_recovery_fraction"], 1.15)
        _, invalid = wastewater_factor(
            {**sample, "pmmov_post_dilution": 10}, "copies/mL"
        )
        self.assertAlmostEqual(invalid["pmmov_recovery_fraction"], 11.5)
        self.assertEqual(invalid["status"], "recovery_outside_acceptance_range")

    def test_missing_total_bases_only_blocks_yield_model(self):
        data = self.load({**self.sample, "total_effective_bases": "NA"})
        self.assertFalse(calibrate(data, METHODS[0], self.settings)[0]["passes"])
        self.assertTrue(calibrate(data, METHODS[1], self.settings)[0]["passes"])

    def test_fastp_stage_and_read_count_are_explicit(self):
        sample = {
            **self.sample,
            "total_effective_bases": "NA",
            "fastp_section": "before_filtering",
        }
        self.path.with_name("SIMULATED_01.fastp.json").write_text(
            json.dumps(
                {
                    "summary": {
                        "before_filtering": {
                            "total_reads": 20000000,
                            "total_bases": 3000000000,
                        },
                        "after_filtering": {"total_reads": 999, "total_bases": 1234},
                    }
                }
            )
        )
        self.assertEqual(self.load(sample)["total_bases"], 3000000000)
        for stage in ["after_filtering", "NA"]:
            with self.assertRaises(ValueError):
                self.load({**sample, "fastp_section": stage})

    def test_unspiked_control_has_no_borrowed_curve(self):
        data = self.load(
            {
                **self.sample,
                "ercc_spike_volume_ul": 0,
                "sample_type": "negative_control",
            }
        )
        for method in METHODS:
            rows, _, audit = self.quantify(data, method=method)
            self.assertTrue(all(r[CONC] is None for r in rows))
            self.assertEqual(audit["status"], "not_applicable")

    def test_invalid_or_duplicate_source_rejected(self):
        write_tsv(self.path, self.rows + [self.rows[0]])
        with self.assertRaises(ValueError):
            self.load()
        write_tsv(self.path, self.rows)
        with self.assertRaises(ValueError):
            self.load({**self.sample, "sample_id": "OTHER"})

    def test_aligned_bases_are_not_coverage_breadth(self):
        signal = self.load()["signals"]["SIM_RNA"]
        self.assertNotAlmostEqual(signal["aligned_bases"], signal["covered_bases"])

    def test_impossible_total_base_accounting_rejected(self):
        with self.assertRaisesRegex(ValueError, "Summed aligned bases"):
            self.load({**self.sample, "total_effective_bases": 100})

    def test_validated_reporting_requires_external_evidence(self):
        data = self.load({**self.sample, "calibration_transfer_verified": 0})
        row = self.target(self.quantify(data)[0])
        self.assertIsNone(row[CONC])
        self.assertIsNone(row["concentration_extract_copies_ul"])
        self.assertIsNotNone(row["estimated_reaction_copies_ul"])
        self.assertIsNotNone(row["estimated_wastewater_concentration"])
        exploratory = self.target(
            self.quantify(data, settings=Settings(reporting_mode="exploratory"))[0]
        )
        self.assertIsNotNone(exploratory[CONC])
        self.assertFalse(exploratory["target_validation_verified"])

    def test_target_specific_loq_overrides_sample_value(self):
        data = self.load()
        targets = load_targets(self.root / "targets.tsv", [data])
        targets["SIM_RNA"]["loq_reaction_copies_ul"] = 1000
        row = self.target(self.quantify(data, targets=targets)[0])
        self.assertIsNone(row[CONC])
        self.assertIn("below_loq", row["quantification_status"])
        targets["SIM_RNA"]["loq_reaction_copies_ul"] = "NA"
        self.assertIsNotNone(self.target(self.quantify(data, targets=targets)[0])[CONC])

    def test_outside_range_requires_explicit_extrapolation(self):
        # Increase target signal while keeping total-base accounting possible.
        for row in self.rows:
            if row["Accession"] == "SIM_RNA":
                row["mean_coverage"] *= 10000
        write_tsv(self.path, self.rows)
        data = self.load()
        for method in METHODS:
            row = self.target(self.quantify(data, method=method)[0])
            self.assertIsNone(row[CONC])
            self.assertIsNotNone(row["estimated_reaction_copies_ul"])
            self.assertIn("outside_calibration_range", row["quantification_status"])
            row = self.target(
                self.quantify(
                    data, method=method, settings=Settings(allow_extrapolation=True)
                )[0]
            )
            self.assertIsNotNone(row[CONC])

    def test_exploratory_low_signal_is_estimated_without_assigning_a_loq(self):
        data = self.load({**self.sample, "loq_reaction_copies_ul": "NA"})
        signal = data["signals"]["SIM_RNA"]
        # Isolate the reporting policy with a positive synthetic sub-range signal.
        signal["aligned_bases"] *= 1e-12
        signal["depth"] *= 1e-12
        for method in METHODS:
            row = self.target(
                self.quantify(
                    data, method=method, settings=Settings(reporting_mode="exploratory")
                )[0]
            )
            self.assertIsNotNone(row[CONC])
            self.assertGreater(row[CONC], 0)
            self.assertFalse(row["within_calibration_range"])
            self.assertIn("outside_calibration_range", row["quantification_status"])
            self.assertIn("loq_unavailable", row["quantification_status"])
            self.assertNotIn("below_loq", row["quantification_status"])
            strict = self.target(
                self.quantify(
                    data,
                    method=method,
                    settings=Settings(
                        reporting_mode="exploratory", allow_extrapolation=False
                    ),
                )[0]
            )
            self.assertIsNone(strict[CONC])
            self.assertEqual(
                strict["estimated_reaction_copies_ul"],
                row["estimated_reaction_copies_ul"],
            )

    def test_duplex_rna_uses_two_strand_mass(self):
        data = self.load()
        targets = load_targets(self.root / "targets.tsv", [data])
        targets["SIM_RNA"]["molecule_type"] = "dsRNA"
        targets["SIM_RNA"]["molecular_mass_per_nt_g_mol"] = 680
        row = self.target(self.quantify(data, targets=targets)[0])
        self.assertAlmostEqual(row["estimated_reaction_copies_ul"], 10000 / 28)

    def test_invalid_assembly_length_suppresses_both_estimators(self):
        data = self.load()
        targets = load_targets(self.root / "targets.tsv", [data])
        targets["SIM_RNA"]["_reference_length_valid"] = False
        for method in METHODS:
            row = self.target(self.quantify(data, method=method, targets=targets)[0])
            self.assertIsNone(row["estimated_reaction_copies_ul"])
            self.assertIn(
                "inconsistent_database_assembly_length", row["quantification_status"]
            )

    def test_molecular_mass_factor_controls_yield_conversion(self):
        data = self.load()
        targets = load_targets(self.root / "targets.tsv", [data])
        targets["SIM_RNA"]["molecular_mass_per_nt_g_mol"] = 660
        row = self.target(self.quantify(data, targets=targets)[0])
        self.assertAlmostEqual(row[CONC], 10000 / 14 * 340 / 660)
        targets["SIM_RNA"]["molecular_mass_per_nt_g_mol"] = "NA"
        row = self.target(self.quantify(data, targets=targets)[0])
        self.assertIsNone(row[CONC])
        self.assertIn("missing_molecular_mass_factor", row["quantification_status"])

    def test_high_r2_does_not_override_slope_criterion(self):
        for row in self.rows:
            if row["Accession"].startswith("ERCC-"):
                row["mean_coverage"] = (
                    float(row["mean_coverage"]) * float(row["Length"])
                ) ** 0.4 / float(row["Length"])
                row["covered_bases"] = min(
                    float(row["covered_bases"]),
                    int(row["mean_coverage"] * row["Length"]),
                )
        write_tsv(self.path, self.rows)
        calibration = calibrate(self.load(), METHODS[0], self.settings)[0]
        self.assertAlmostEqual(calibration["fit"]["r2"], 1)
        self.assertAlmostEqual(calibration["fit"]["slope"], 0.4)
        self.assertFalse(calibration["passes"])

    def test_levels_must_be_adjacent(self):
        data = self.load()
        for point in data["points"]:
            point["fit_used"] = point["level"] in {1, 3, 5, 7, 9}
        calibration = calibrate(data, METHODS[1], self.settings)[0]
        self.assertEqual(calibration["longest_adjacent_levels"], 1)
        self.assertFalse(calibration["passes"])

    def test_replicate_cv_is_diagnostic_unless_requested(self):
        def samples():
            return [
                {
                    "sample": {"replicate_group": "G"},
                    "points": [
                        {
                            "fit_used": True,
                            "sequencing_yield": value,
                            "ercc_id": "ERCC-00002",
                        }
                    ],
                }
                for value in [1, 1, 100]
            ]

        data = samples()
        self.assertFalse(replicate_qc(data)[0]["passes"])
        self.assertTrue(data[0]["points"][0]["fit_used"])
        data = samples()
        replicate_qc(data, Settings(filter_replicate_cv=True))
        self.assertFalse(data[0]["points"][0]["fit_used"])
        data = samples()[:2]
        replicate_qc(data, Settings(filter_replicate_cv=True))
        self.assertTrue(data[0]["points"][0]["fit_used"])

    def test_partial_segment_does_not_shorten_reference(self):
        rows = [r for r in self.rows if r["Accession"] != "SIM_SEG_B"]
        write_tsv(self.path, rows)
        data = self.load()
        groups = reference_groups(data, load_targets(self.root / "targets.tsv", [data]))
        self.assertEqual(groups["SIM_SEG_A"]["reference_length_nt"], 10000)
        self.assertEqual(groups["SIM_SEG_A"]["detected_reference_length_nt"], 4000)
        quantified, _, _ = self.quantify(data)
        row = next(r for r in quantified if r["Accession"] == "SIM_SEG_A")
        self.assertIn("partial_reference_detected", row["quantification_status"])
        self.assertAlmostEqual(row[CONC], 4000 / 14, places=7)

    def test_quality_filtered_bam_blocks(self):
        try:
            import pysam
        except ImportError:
            self.skipTest("Optional BAM dependency is not installed")
        bam_path = self.root / "reads.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "REF", "LN": 100}],
        }
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            for i, (start, flags, mapq, nh) in enumerate(
                [
                    (0, 0, 60, 1),
                    (0, 1024, 60, 1),
                    (0, 256, 60, 1),
                    (0, 0, 0, 1),
                    (0, 0, 60, 2),
                    (1, 2048, 60, 1),
                    (2, 0, 60, 1),
                ]
            ):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "ACGT"
                read.flag = flags
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = mapq
                read.cigarstring = "4M"
                read.set_tag("NH", nh)
                read.query_qualities = pysam.qualitystring_to_array("IIII")
                bam.write(read)
        pysam.index(str(bam_path))
        self.assertEqual(
            main(
                [
                    "collect-bam",
                    "--bam",
                    str(bam_path),
                    "--output",
                    str(self.root / "metrics.tsv"),
                ]
            ),
            0,
        )
        _, rows = read_tsv(self.root / "metrics.tsv")
        self.assertEqual(int(rows[0]["aligned_bases"]), 8)
        self.assertEqual(int(rows[0]["covered_bases"]), 6)
        self.assertEqual(int(rows[0]["read_count"]), 2)
        self.assertTrue((self.root / "metrics.json").is_file())


if __name__ == "__main__":
    unittest.main()
