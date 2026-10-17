import contextlib
import importlib.util
import io
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from hybcapseq_quant.cli import choose_table, main
from hybcapseq_quant.core import read_tsv, write_tsv

SPEC = importlib.util.spec_from_file_location("portable_demo", Path(__file__).parents[1] / "examples/quant/generate_demo.py")
DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEMO)


class PortableDirectoriesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def command(self, args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return main([str(v) for v in args])

    def test_default_init_is_generic_not_historical_records(self):
        out = self.root / "inputs"
        self.assertEqual(self.command(["init", "--outdir", out]), 0)
        fields, rows = read_tsv(out / "samples.tsv")
        self.assertIn("reaction_volume_ul", fields)
        self.assertEqual(rows, [])
        self.assertFalse((out / "samples_20261017.tsv").exists())

    def test_init_discovers_current_sample_names_and_leaves_inputs_unknown(self):
        results = self.root / "results for another cohort"
        DEMO.create(results / "nested folder", sid="NEW_SITE_X")
        DEMO.create(results / "another folder", sid="NEW_SITE_Y")
        out = self.root / "inputs"
        self.assertEqual(self.command(["init", "--esviritu-dir", results, "--outdir", out]), 0)
        _, rows = read_tsv(out / "samples.tsv")
        self.assertEqual({r["sample_ID"] for r in rows}, {"NEW_SITE_X", "NEW_SITE_Y"})
        for row in rows:
            self.assertEqual(row["sample_input_mass_ng"], "NA")
            self.assertEqual(row["ercc_dilution"], "NA")
            self.assertFalse(Path(row["info_tsv"]).is_absolute())
            self.assertTrue((out / row["info_tsv"]).is_file())

    def test_historical_13_samples_require_explicit_preset(self):
        out = self.root / "historical inputs"
        self.assertEqual(self.command(["init", "--preset", "20261017", "--outdir", out]), 0)
        _, rows = read_tsv(out / "samples_20261017.tsv")
        self.assertEqual(len(rows), 13)
        self.assertFalse((out / "samples.tsv").exists())

    def test_init_can_accept_entirely_new_manual_sample_ids(self):
        out = self.root / "new inputs"
        self.assertEqual(self.command(["init", "--sample", "FUTURE_01", "FUTURE_02", "--outdir", out]), 0)
        _, rows = read_tsv(out / "samples.tsv")
        self.assertEqual([r["sample_ID"] for r in rows], ["FUTURE_01", "FUTURE_02"])

    def test_moving_batch_and_changing_CWD_preserves_relative_inputs(self):
        original = self.root / "original batch with spaces"
        results = original / "any output name" / "per sample"
        sample, _ = DEMO.create(results, sid="FUTURE_SITE")
        inputs = original / "input files"
        self.assertEqual(self.command(["init", "--esviritu-dir", original / "any output name", "--outdir", inputs]), 0)
        fields, rows = read_tsv(inputs / "samples.tsv")
        rows[0].update(sample)
        write_tsv(inputs / "samples.tsv", rows, fields)
        relocated = self.root / "renamed batch on another path"
        shutil.move(str(original), relocated)
        other_cwd = self.root / "unrelated working directory"
        other_cwd.mkdir()
        previous = Path.cwd()
        try:
            os.chdir(other_cwd)
            code = self.command(["--samples", relocated / "input files/samples.tsv",
                                 "--processed-dir", relocated])
        finally:
            os.chdir(previous)
        self.assertEqual(code, 0)
        root = relocated / "HybCapSeq_Quant_results"
        self.assertFalse((other_cwd / "HybCapSeq_Quant_results").exists())
        for method in ("method1_yield", "method2_depth"):
            _, values = read_tsv(next((root / method / "FUTURE_SITE").glob("*.info.quantified.*.tsv")))
            target = next(v for v in values if v["Accession"] == "SIM_RNA")
            self.assertAlmostEqual(float(target["浓度concentration_ERCC"]), 10000 / 14, places=7)

    def test_current_explicit_result_directory_overrides_saved_old_path(self):
        old = self.root / "previous data"
        new = self.root / "current data"
        DEMO.create(old, sid="SAME_ID")
        DEMO.create(new, sid="SAME_ID")
        old_path = old / "SAME_ID.detected_virus.info.tsv"
        new_path = new / old_path.name
        sample = {"sample_ID": "SAME_ID", "info_tsv": str(old_path)}
        self.assertEqual(choose_table(sample, new, None, self.root), new_path)
        self.assertEqual(choose_table(sample, new, old_path, self.root), old_path)

    def test_explicit_outdir_can_be_outside_processed_directory(self):
        results = self.root / "input data"
        DEMO.create(results)
        processed = self.root / "working location"
        separate = self.root / "independent output disk"
        self.assertEqual(self.command(["--samples", results / "samples.tsv", "--esviritu-dir", results,
                                      "--processed-dir", processed, "--outdir", separate,
                                      "--method", "method2_depth"]), 0)
        self.assertTrue((separate / "run_summary.tsv").is_file())
        self.assertFalse((processed / "HybCapSeq_Quant_results").exists())

    def test_existing_exclusion_file_is_not_overwritten_on_init(self):
        out = self.root / "inputs"
        out.mkdir()
        (out / "ercc_qc.tsv").write_text("important existing exclusions\n")
        self.assertEqual(self.command(["init", "--outdir", out]), 2)
        self.assertEqual((out / "ercc_qc.tsv").read_text(), "important existing exclusions\n")
        self.assertFalse((out / "samples.tsv").exists())

    def test_ambiguous_same_sample_across_batches_is_rejected(self):
        results = self.root / "many batches"
        DEMO.create(results / "batch A", sid="SAME_ID")
        DEMO.create(results / "batch B", sid="SAME_ID")
        self.assertEqual(self.command(["init", "--esviritu-dir", results, "--outdir", self.root / "inputs"]), 2)


if __name__ == "__main__":
    unittest.main()
