"""Command-line contracts, provenance, and relocation checks."""

import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from hybcapseq_quant.cli import choose_table, main
from hybcapseq_quant.core import CONC, WWCONC, METHODS, read_tsv, write_tsv

SPEC = importlib.util.spec_from_file_location(
    "portable_fixture", Path(__file__).parents[1] / "examples/quant/generate_demo.py"
)
DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEMO)


class PortableDirectoriesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def command(self, args):
        with (
            contextlib.redirect_stdout(io.StringIO()),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            return main([str(v) for v in args])

    def test_empty_template_has_no_experimental_defaults(self):
        out = self.root / "inputs"
        self.assertEqual(self.command(["init", "--outdir", out]), 0)
        fields, rows = read_tsv(out / "samples.tsv")
        self.assertIn("library_input_volume_ul", fields)
        self.assertEqual(rows, [])

    def test_discovery_leaves_experimental_fields_unknown(self):
        results = self.root / "results with spaces"
        DEMO.create(results / "first", sid="SAMPLE_A")
        DEMO.create(results / "second", sid="SAMPLE_B")
        out = self.root / "inputs"
        self.assertEqual(
            self.command(["init", "--esviritu-dir", results, "--outdir", out]), 0
        )
        _, rows = read_tsv(out / "samples.tsv")
        self.assertEqual({r["sample_id"] for r in rows}, {"SAMPLE_A", "SAMPLE_B"})
        for row in rows:
            for name in [
                "library_input_mass_ng",
                "ercc_dilution_factor",
                "reaction_volume_ul",
                "fastp_section",
            ]:
                self.assertEqual(row[name], "NA")
            self.assertFalse(Path(row["info_tsv"]).is_absolute())
            self.assertTrue((out / row["info_tsv"]).is_file())

    def test_manual_identifiers(self):
        out = self.root / "inputs"
        self.assertEqual(
            self.command(["init", "--sample", "SAMPLE_A", "SAMPLE_B", "--outdir", out]),
            0,
        )
        _, rows = read_tsv(out / "samples.tsv")
        self.assertEqual([r["sample_id"] for r in rows], ["SAMPLE_A", "SAMPLE_B"])

    def test_relocation_and_cwd_independence(self):
        original = self.root / "project with spaces"
        results = original / "mapping results" / "sample"
        sample, _ = DEMO.create(results, sid="SAMPLE_A")
        inputs = original / "input tables"
        self.assertEqual(
            self.command(
                [
                    "init",
                    "--esviritu-dir",
                    original / "mapping results",
                    "--outdir",
                    inputs,
                ]
            ),
            0,
        )
        fields, rows = read_tsv(inputs / "samples.tsv")
        rows[0].update(sample)
        write_tsv(inputs / "samples.tsv", rows, fields)
        relocated = self.root / "relocated project"
        shutil.move(str(original), relocated)
        unrelated = self.root / "different working directory"
        unrelated.mkdir()
        previous = Path.cwd()
        try:
            os.chdir(unrelated)
            code = self.command(
                [
                    "run",
                    "--samples",
                    relocated / "input tables/samples.tsv",
                    "--targets",
                    relocated / "mapping results/sample/targets.tsv",
                    "--processed-dir",
                    relocated,
                ]
            )
        finally:
            os.chdir(previous)
        self.assertEqual(code, 0)
        root = relocated / "hybcapseq_quant_results"
        self.assertFalse((unrelated / "hybcapseq_quant_results").exists())
        for method in METHODS:
            _, values = read_tsv(
                root / method / "SAMPLE_A/SAMPLE_A.reference_quantification.tsv"
            )
            target = next(v for v in values if v["Accession"] == "SIM_RNA")
            self.assertAlmostEqual(float(target[CONC]), 10000 / 14, places=7)

    def test_yield_requires_catalog_or_explicit_chemistry_manifest(self):
        source = self.root / "mapping"
        DEMO.create(source)
        out = self.root / "output"
        self.assertEqual(
            self.command(
                [
                    "run",
                    "--samples",
                    source / "samples.tsv",
                    "--esviritu-dir",
                    source,
                    "--outdir",
                    out,
                ]
            ),
            2,
        )
        self.assertFalse(out.exists())

    def test_explicit_location_priority(self):
        left, right = self.root / "left", self.root / "right"
        DEMO.create(left, sid="SAMPLE_A")
        DEMO.create(right, sid="SAMPLE_A")
        source = left / "SAMPLE_A.detected_virus.info.tsv"
        sample = {"sample_id": "SAMPLE_A", "info_tsv": str(source)}
        self.assertEqual(
            choose_table(sample, right, None, self.root), right / source.name
        )
        self.assertEqual(choose_table(sample, right, source, self.root), source)

    def test_explicit_outdir_takes_priority(self):
        source = self.root / "mapping"
        DEMO.create(source)
        processed, separate = self.root / "project", self.root / "outputs"
        self.assertEqual(
            self.command(
                [
                    "run",
                    "--samples",
                    source / "samples.tsv",
                    "--targets",
                    source / "targets.tsv",
                    "--esviritu-dir",
                    source,
                    "--processed-dir",
                    processed,
                    "--outdir",
                    separate,
                    "--method",
                    METHODS[1],
                ]
            ),
            0,
        )
        self.assertTrue((separate / "run_summary.tsv").is_file())
        self.assertFalse((processed / "hybcapseq_quant_results").exists())

    def test_existing_inputs_are_not_overwritten(self):
        out = self.root / "inputs"
        out.mkdir()
        (out / "ercc_qc.tsv").write_text("Existing records\n")
        self.assertEqual(self.command(["init", "--outdir", out]), 2)
        self.assertEqual((out / "ercc_qc.tsv").read_text(), "Existing records\n")
        self.assertFalse((out / "samples.tsv").exists())

    def test_ambiguous_same_sample_is_rejected(self):
        mapping = self.root / "mapping"
        DEMO.create(mapping / "batch_a", sid="SAMPLE_A")
        DEMO.create(mapping / "batch_b", sid="SAMPLE_A")
        self.assertEqual(
            self.command(
                ["init", "--esviritu-dir", mapping, "--outdir", self.root / "inputs"]
            ),
            2,
        )

    def test_source_columns_preserved_and_new_fields_follow_rpkmf(self):
        source, out = self.root / "mapping", self.root / "quantification"
        _, original = DEMO.create(source)
        self.assertEqual(
            self.command(
                [
                    "run",
                    "--samples",
                    source / "samples.tsv",
                    "--targets",
                    source / "targets.tsv",
                    "--outdir",
                    out,
                    "--esviritu-dir",
                    source,
                ]
            ),
            0,
        )
        for method in METHODS:
            folder = out / method / "SIMULATED_01"
            fields, values = read_tsv(
                folder / "SIMULATED_01.reference_quantification.tsv"
            )
            index = fields.index("RPKMF")
            self.assertEqual(fields[index + 1 : index + 3], [CONC, WWCONC])
            self.assertTrue(all(c.isascii() for c in fields))
            for row, annotated in zip(original, values):
                for key, value in row.items():
                    self.assertEqual(annotated[key], str(value))
            for figure in ["calibration_curve", "reaction_concentration_curve"]:
                for suffix in ["png", "pdf"]:
                    self.assertGreater(
                        (folder / f"SIMULATED_01.{figure}.{suffix}").stat().st_size,
                        1000,
                    )
            _, points = read_tsv(folder / "SIMULATED_01.calibration_points.tsv")
            self.assertEqual(len(points), 92)
            self.assertTrue(all(p["included_in_model"] == "true" for p in points))
            self.assertLess(max(abs(float(p["residual_log10"])) for p in points), 1e-8)
        manifest = json.loads((out / "run_manifest.json").read_text())
        self.assertEqual(len(manifest["input_files"]["samples"]["sha256"]), 64)
        self.assertIn("settings", manifest)

    def test_batch_mode_and_single_sample_selection(self):
        source = self.root / "mapping"
        a, _ = DEMO.create(source / "first", sid="SAMPLE_A")
        b, _ = DEMO.create(source / "second", sid="SAMPLE_B")
        write_tsv(self.root / "samples.tsv", [a, b])
        common = [
            "run",
            "--samples",
            self.root / "samples.tsv",
            "--targets",
            source / "first/targets.tsv",
            "--esviritu-dir",
            source,
        ]
        self.assertEqual(self.command(common + ["--outdir", self.root / "all"]), 0)
        _, rows = read_tsv(self.root / "all/run_summary.tsv")
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            self.command(
                common + ["--sample", "SAMPLE_B", "--outdir", self.root / "one"]
            ),
            0,
        )
        _, rows = read_tsv(self.root / "one/run_summary.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["sample_id"] for r in rows}, {"SAMPLE_B"})

    def test_existing_results_are_not_mixed(self):
        source = self.root / "mapping"
        DEMO.create(source)
        out = self.root / "output"
        out.mkdir()
        (out / "existing.txt").write_text("Existing analysis\n")
        self.assertEqual(
            self.command(
                [
                    "run",
                    "--samples",
                    source / "samples.tsv",
                    "--targets",
                    source / "targets.tsv",
                    "--esviritu-dir",
                    source,
                    "--outdir",
                    out,
                ]
            ),
            2,
        )
        self.assertEqual(list(out.iterdir()), [out / "existing.txt"])

    def test_fastq_uses_actual_base_lengths(self):
        fastq = self.root / "reads.fastq"
        fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGTAC\n+\nIIIIII\n")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(
                main(
                    ["measure-input-bases", "--sample", "SAMPLE_A", "--r1", str(fastq)]
                ),
                0,
            )
        self.assertEqual(json.loads(out.getvalue())["total_effective_bases"], 10)

    def test_unknown_fields_and_invalid_identifiers_are_rejected(self):
        source = self.root / "mapping"
        sample, _ = DEMO.create(source)
        sample["unexpected_field"] = 1
        write_tsv(source / "samples.tsv", [sample])
        self.assertEqual(
            self.command(
                ["run", "--samples", source / "samples.tsv", "--esviritu-dir", source]
            ),
            2,
        )
        self.assertEqual(
            self.command(
                ["init", "--sample", "../SAMPLE", "--outdir", self.root / "inputs"]
            ),
            2,
        )

    def test_invalid_settings_are_rejected(self):
        source = self.root / "mapping"
        DEMO.create(source)
        for args in [
            ["--min-r2", "1.1"],
            ["--min-levels", "3"],
            ["--slope-min", "2", "--slope-max", "1"],
        ]:
            self.assertEqual(
                self.command(
                    [
                        "run",
                        "--samples",
                        source / "samples.tsv",
                        "--targets",
                        source / "targets.tsv",
                        "--esviritu-dir",
                        source,
                    ]
                    + args
                ),
                2,
            )

    def test_target_template_contains_only_virus_references(self):
        source = self.root / "mapping"
        DEMO.create(source)
        out = self.root / "targets.tsv"
        self.assertEqual(
            self.command(
                ["prepare-targets", "--esviritu-dir", source, "--output", out]
            ),
            0,
        )
        _, rows = read_tsv(out)
        self.assertEqual(len(rows), 3)
        self.assertFalse(any(r["accession"].startswith("ERCC-") for r in rows))


if __name__ == "__main__":
    unittest.main()
