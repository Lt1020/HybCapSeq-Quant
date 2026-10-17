import contextlib
import csv
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from hybcapseq_quant.cli import main, load_targets
from hybcapseq_quant.core import (CONC, WWCONC, METHODS, NA, STANDARD_PATH, Settings,
    calibrate, ercc_amounts, input_bases, load_sample, quantify_rows, read_tsv,
    reference_groups, replicate_qc, wastewater_factor, write_tsv)

SPEC = importlib.util.spec_from_file_location("demo", Path(__file__).parents[1] / "examples/quant/generate_demo.py")
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

    def load(self, sample=None, settings=None):
        return load_sample(sample or self.sample, self.path, self.standards, self.root, settings or self.settings)

    def test_official_92_transcripts_22_levels_and_excluded_four(self):
        self.assertEqual(len(self.standards), 92)
        self.assertEqual(len({r["level"] for r in self.standards}), 22)
        ids = {r["ERCC_ID"] for r in self.standards}
        self.assertFalse(ids.intersection({"ERCC-00007", "ERCC-00018", "ERCC-00023", "ERCC-00128"}))

    def test_attomoles_dilution_and_reaction_concentration(self):
        amounts = ercc_amounts(self.standards, self.sample)
        # Independent hand calculation for the highest Mix1 standard.
        expected_n = 30000 * 1e-18 * NA / 1000 * 2
        self.assertAlmostEqual(amounts["ERCC-00130"]["added_copies"], expected_n)
        self.assertAlmostEqual(amounts["ERCC-00130"]["reaction_copies_ul"], expected_n / 14)
        self.assertAlmostEqual(amounts["ERCC-00130"]["extract_equivalent_copies_ul"], expected_n / 2)

    def test_both_protocol_methods_recover_known_RNA_and_segmented_genome(self):
        data = self.load()
        targets = load_targets(None, [data])
        for method in METHODS:
            cal, _, _, _ = calibrate(data, method, self.settings)
            self.assertTrue(cal["passes"], cal)
            self.assertAlmostEqual(cal["fit"]["slope"], 1)
            self.assertAlmostEqual(cal["fit"]["r2"], 1)
            results, assemblies, _ = quantify_rows(data, cal, targets, self.settings)
            target = next(r for r in results if r["Accession"] == "SIM_RNA")
            self.assertAlmostEqual(target[CONC], 10000 / 14, places=7)
            self.assertAlmostEqual(target["concentration_extract_copies_ul"], 5000, places=7)
            # Cww = 5000 * (1000*80/300) * 0.46/(40*0.25)
            self.assertAlmostEqual(target[WWCONC], 61333.3333333333, places=6)
            segmented = [r for r in results if r["Accession"].startswith("SIM_SEG_")]
            self.assertEqual(len(segmented), 2)
            for row in segmented:
                self.assertAlmostEqual(row[CONC], 10000 / 14, places=7)
            self.assertEqual(len(assemblies), 2)  # one segmented genome, not two independent viruses

    def test_liter_vs_milliliter_exact_factor(self):
        a, _ = wastewater_factor(self.sample, "copies/mL")
        b, _ = wastewater_factor(self.sample, "copies/L")
        self.assertEqual(b, a * 1000)

    def test_mass_input_includes_DNA_and_RNA_actual_records(self):
        _, records = read_tsv(STANDARD_PATH.with_name("samples_20261017.tsv"))
        self.assertEqual(len(records), 13)
        for row in records:
            mass = float(row["sample_input_mass_ng"])
            mixed = float(row["RNA_conc_ng_ul"]) + float(row["DNA_conc_ng_ul"])
            self.assertAlmostEqual(mass, mixed * float(row["library_extract_ul"]), places=7)
        self.assertGreater(float(records[0]["sample_input_mass_ng"]), 100)
        self.assertEqual(float(records[-1]["ercc_spike_ul"]), 0)

    def test_known_ERCC_never_back_converted_to_wastewater(self):
        data = self.load()
        cal, _, _, _ = calibrate(data, METHODS[0], self.settings)
        results, _, _ = quantify_rows(data, cal, load_targets(None, [data]), self.settings)
        row = next(r for r in results if r["Accession"] == "ERCC-00130")
        self.assertAlmostEqual(row[CONC], 2580917.4685714287)
        self.assertIsNone(row[WWCONC])

    def test_missing_PMMoV_does_not_invent_recovery_one(self):
        sample = {**self.sample, "pmmov_recovery_fraction": "NA"}
        data = self.load(sample)
        cal, _, _, _ = calibrate(data, METHODS[1], self.settings)
        results, _, audit = quantify_rows(data, cal, load_targets(None, [data]), self.settings)
        row = next(r for r in results if r["Accession"] == "SIM_RNA")
        self.assertIsNotNone(row[CONC])
        self.assertIsNone(row[WWCONC])
        self.assertIn("missing", audit["status"])

    def test_PMMoV_formula_uses_mass_recovery_not_concentration_ratio(self):
        sample = {**self.sample, "pmmov_recovery_fraction": "NA", "pmmov_raw_extract_copies_ul": 10,
            "pmmov_post_extract_copies_ul": 1000, "pmmov_raw_dilution": 1, "pmmov_post_dilution": 10,
            "raw_extract_elution_ul": 80, "raw_extract_input_ul": 300, "pmmov_above_loq": 1}
        _, audit = wastewater_factor(sample, "copies/mL")
        self.assertAlmostEqual(audit["pmmov_recovery_fraction"], 11.5)  # excessive recovery is rejected
        self.assertEqual(audit["status"], "recovery_outside_protocol_review_range")
        sample["pmmov_post_dilution"] = 1
        _, audit = wastewater_factor(sample, "copies/mL")
        self.assertAlmostEqual(audit["pmmov_recovery_fraction"], 1.15)

    def test_missing_total_bases_blocks_yield_but_not_depth(self):
        data = self.load({**self.sample, "total_effective_bases": "NA"})
        self.assertFalse(calibrate(data, METHODS[0], self.settings)[0]["passes"])
        self.assertTrue(calibrate(data, METHODS[1], self.settings)[0]["passes"])

    def test_fastp_section_must_match_actual_reads(self):
        sample = {**self.sample, "total_effective_bases": "NA"}
        j = {"summary": {"before_filtering": {"total_reads": 10000000, "total_bases": 123456789},
                         "after_filtering": {"total_reads": 999, "total_bases": 1234}}}
        self.path.with_name("SIMULATED_01.fastp.json").write_text(json.dumps(j))
        data = self.load(sample)
        self.assertEqual(data["total_bases"], 123456789)
        with self.assertRaises(ValueError):
            self.load({**sample, "fastp_section": "after_filtering"})

    def test_no_ERCC_control_does_not_borrow_sample_calibration(self):
        data = self.load({**self.sample, "ercc_spike_ul": 0, "sample_type": "negative_control"})
        for method in METHODS:
            cal, _, _, _ = calibrate(data, method, self.settings)
            self.assertFalse(cal["passes"])
            result, _, _ = quantify_rows(data, cal, load_targets(None, [data]), self.settings)
            self.assertTrue(all(row[CONC] is None for row in result))

    def test_duplicate_reference_and_wrong_sample_rejected(self):
        write_tsv(self.path, self.rows + [self.rows[0]])
        with self.assertRaises(ValueError):
            self.load()
        write_tsv(self.path, self.rows)
        with self.assertRaises(ValueError):
            self.load({**self.sample, "sample_ID": "WRONG"})

    def test_breadth_not_substituted_for_aligned_bases(self):
        data = self.load()
        signal = data["signals"]["SIM_RNA"]
        self.assertGreater(signal["aligned_bases"], signal["covered_bases"] * 20)

    def test_strict_qc_suppresses_unvalidated_estimates(self):
        data = self.load()
        cal, _, _, _ = calibrate(data, METHODS[1], self.settings)
        result, _, _ = quantify_rows(data, cal, load_targets(None, [data]), Settings(strict_qc=True))
        row = next(r for r in result if r["Accession"] == "SIM_RNA")
        self.assertIsNone(row[CONC])
        self.assertIsNotNone(row["estimated_reaction_copies_ul_before_QC"])

    def test_low_fit_R2_retains_known_ERCC_but_suppresses_targets(self):
        data = self.load()
        cal, _, _, _ = calibrate(data, METHODS[1], self.settings)
        cal["passes"] = False
        cal["errors"] = ["R2_below_threshold"]
        result, _, _ = quantify_rows(data, cal, load_targets(None, [data]), self.settings)
        self.assertIsNotNone(result[0][CONC])
        self.assertIsNone(next(r for r in result if r["Accession"] == "SIM_RNA")[CONC])

    def test_replicate_filter_requires_three_observations(self):
        datasets = []
        for i, y in enumerate([1, 1, 100]):
            datasets.append({"sample": {"replicate_group": "G"}, "points": [
                {"fit_used": True, "Y_seq_i": y, "ERCC_ID": "ERCC-00002"}]})
        audit = replicate_qc(datasets)
        self.assertGreater(audit[0]["Y_seq_cv"], 0.35)
        self.assertTrue(all(not d["points"][0]["fit_used"] for d in datasets))

    def test_cli_preserves_original_columns_and_writes_both_methods_and_plots(self):
        with contextlib.redirect_stdout(io.StringIO()):
            code = main(["--samples", str(self.root / "samples.tsv"), "--esviritu-dir", str(self.root),
                         "--sample", "SIMULATED_01", "--outdir", str(self.root / "results")])
        self.assertEqual(code, 0)
        for method in METHODS:
            folder = self.root / "results" / method / "SIMULATED_01"
            fields, rows = read_tsv(next(folder.glob("*.info.quantified.*.tsv")))
            original = list(self.rows[0])
            at = fields.index("RPKMF")
            self.assertEqual(fields[at + 1:at + 3], [CONC, WWCONC])
            for a, b in zip(self.rows, rows):
                for col in original:
                    self.assertEqual(b[col], str(a[col]))
            self.assertGreater((folder / "SIMULATED_01.ERCC_standard_curve.png").stat().st_size, 1000)
            self.assertGreater((folder / "SIMULATED_01.ERCC_standard_curve.pdf").stat().st_size, 1000)
            self.assertEqual(len(json.loads((folder / "SIMULATED_01.quantification_audit.json").read_text())["calibration"]["fit_ERCC_IDs"]), 92)

    def test_actual_FASTQ_base_count_not_nominal_PE150(self):
        p = self.root / "reads.fastq"
        p.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGTAC\n+\nIIIIII\n")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(main(["measure-input-bases", "--sample", "S", "--r1", str(p)]), 0)
        self.assertEqual(json.loads(out.getvalue())["total_effective_bases"], 10)

    def test_out_of_range_target_not_silently_extrapolated(self):
        for row in self.rows:
            if row["Accession"] == "SIM_RNA":
                row["mean_coverage"] *= 1e10
        write_tsv(self.path, self.rows)
        data = self.load()
        for method in METHODS:
            cal, _, _, _ = calibrate(data, method, self.settings)
            result, _, _ = quantify_rows(data, cal, load_targets(None, [data]), self.settings)
            row = next(r for r in result if r["Accession"] == "SIM_RNA")
            self.assertIsNone(row[CONC])
            self.assertIsNotNone(row["estimated_reaction_copies_ul_before_QC"])
            self.assertIn("outside_calibration_range", row["quant_status"])

    def test_molecular_mass_factor_changes_yield_copies_and_unknown_is_NA(self):
        data = self.load()
        targets = load_targets(None, [data])
        cal, _, _, _ = calibrate(data, METHODS[0], self.settings)
        targets["SIM_RNA"]["mass_per_nt_g_mol"] = 660
        result, _, _ = quantify_rows(data, cal, targets, self.settings)
        row = next(r for r in result if r["Accession"] == "SIM_RNA")
        self.assertAlmostEqual(row[CONC], 10000 / 14 * 340 / 660)
        targets["SIM_RNA"]["mass_per_nt_g_mol"] = "NA"
        result, _, _ = quantify_rows(data, cal, targets, self.settings)
        row = next(r for r in result if r["Accession"] == "SIM_RNA")
        self.assertIsNone(row[CONC])
        self.assertIn("missing_molecular_mass_factor", row["quant_status"])

    def test_perfect_R2_does_not_override_wrong_yield_slope(self):
        for row in self.rows:
            if row["Accession"].startswith("ERCC-"):
                n = float(row["mean_coverage"]) * float(row["Length"])
                row["mean_coverage"] = n ** 0.4 / float(row["Length"])
        write_tsv(self.path, self.rows)
        cal, _, _, _ = calibrate(self.load(), METHODS[0], self.settings)
        self.assertAlmostEqual(cal["fit"]["r2"], 1)
        self.assertAlmostEqual(cal["fit"]["slope"], 0.4)
        self.assertFalse(cal["passes"])

    def test_concentration_levels_must_be_adjacent(self):
        data = self.load()
        for point in data["points"]:
            point["fit_used"] = point["level"] in {1, 3, 5, 7, 9}
        cal, _, _, _ = calibrate(data, METHODS[1], self.settings)
        self.assertGreater(cal["n_fitted"], 4)
        self.assertEqual(cal["longest_adjacent_levels"], 1)
        self.assertFalse(cal["passes"])

    def test_BAM_filters_and_exact_covered_blocks(self):
        try:
            import pysam
        except ImportError:
            self.skipTest("optional pysam is not installed")
        bam_path = self.root / "test.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "REF", "LN": 100}]}
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            for i, (start, flag_bits, mapq, nh) in enumerate([
                (0, 0, 60, 1), (0, 1024, 60, 1), (0, 256, 60, 1),
                (0, 0, 0, 1), (0, 0, 60, 2), (1, 2048, 60, 1), (2, 0, 60, 1)]):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "ACGT"
                read.flag = flag_bits
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = mapq
                read.cigarstring = "4M"
                read.set_tag("NH", nh)
                read.query_qualities = pysam.qualitystring_to_array("IIII")
                bam.write(read)
        pysam.index(str(bam_path))
        with contextlib.redirect_stderr(io.StringIO()):
            code = main(["collect-bam", "--bam", str(bam_path), "--output", str(self.root / "metrics.tsv")])
        self.assertEqual(code, 0)
        _, rows = read_tsv(self.root / "metrics.tsv")
        self.assertEqual(int(rows[0]["aligned_bases"]), 8)
        self.assertEqual(int(rows[0]["covered_bases"]), 6)
        self.assertEqual(int(rows[0]["read_count"]), 2)

    def test_batch_mode_processes_all_manifest_rows(self):
        sid = "SIMULATED_02"
        second = {**self.sample, "sample_ID": sid}
        write_tsv(self.root / (sid + ".detected_virus.info.tsv"), [{**r, "sample_ID": sid} for r in self.rows])
        write_tsv(self.root / "batch.tsv", [self.sample, second])
        with contextlib.redirect_stdout(io.StringIO()):
            code = main(["--samples", str(self.root / "batch.tsv"), "--esviritu-dir", str(self.root),
                         "--outdir", str(self.root / "batch_results")])
        self.assertEqual(code, 0)
        _, rows = read_tsv(self.root / "batch_results/run_summary.tsv")
        self.assertEqual(len(rows), 4)
        self.assertEqual({r["sample_ID"] for r in rows}, {"SIMULATED_01", "SIMULATED_02"})


if __name__ == "__main__":
    unittest.main()
