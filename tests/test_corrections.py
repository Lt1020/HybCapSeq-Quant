"""Bias correction checked against independent numerical and grouping truths."""

import copy
import csv
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from hybcapseq_quant.cli import main, parser, settings_from_args
from hybcapseq_quant.core import Settings, calibrate
from hybcapseq_quant.corrections import (
    analyze,
    eligible_rows,
    factor,
    feasibility,
    write_rows,
)


def fixture():
    rows = []
    for parent in range(4):
        for extraction in range(2):
            for dilution in range(3):
                for method, coefficient in [
                    ("depth_calibration", 10),
                    ("sequencing_yield", 5),
                ]:
                    truth = (parent + 1) * 100 * (extraction + 1)
                    rows.append(
                        dict(
                            batch="synthetic",
                            sample_id=f"P{parent}_E{extraction}_D{dilution}",
                            parent=f"P{parent}",
                            extraction_id=f"P{parent}_E{extraction}",
                            method=method,
                            target="SIMULATED_TARGET",
                            seq_extract_copies_ul=str(truth / coefficient),
                            dpcr_extract_copies_ul=str(truth),
                            complete_in_range="True",
                            dpcr_comparator_eligible="True",
                            model_accepted="True",
                        )
                    )
    return rows


class CorrectionTests(unittest.TestCase):
    def test_known_factor_and_heldout_prediction(self):
        result = analyze(fixture(), bootstrap=30)
        parameter = next(
            row
            for row in result["parameters"]
            if row["fit_scope"] == "pooled"
            and row["correction_profile"] == "depth_calibration"
        )
        self.assertAlmostEqual(parameter["factor"], 10)
        self.assertEqual(
            parameter["validation_status"], "all_heldout_pairs_within_30pct"
        )
        for row in result["predictions"]:
            self.assertNotIn(
                row["batch"] + "|" + row["parent"],
                row["training_parent_ids"].split(";"),
            )
            if row["correction_profile"] != "shared_depth_yield":
                self.assertAlmostEqual(
                    row["corrected_concentration"], float(row["dpcr_extract_copies_ul"])
                )

    def test_method_specific_factors_cannot_always_be_shared(self):
        result = analyze(fixture(), bootstrap=10)
        shared = next(
            row
            for row in result["feasibility"]
            if row["fit_scope"] == "pooled"
            and row["correction_profile"] == "shared_depth_yield"
        )
        self.assertFalse(shared["all_pairs_30pct_feasible"])
        self.assertGreater(
            shared["feasible_factor_lower"], shared["feasible_factor_upper"]
        )

    def test_repeated_dilutions_do_not_change_parent_weights(self):
        rows, _ = eligible_rows(fixture())
        depth = [row for row in rows if row["method"] == "depth_calibration"]
        baseline = factor(depth)
        replicated = depth + [dict(depth[0]) for _ in range(200)]
        self.assertAlmostEqual(factor(replicated), baseline)

    def test_zero_and_ineligible_comparators_are_audited(self):
        rows = fixture()
        rows[0]["dpcr_extract_copies_ul"] = "0"
        rows[1]["dpcr_comparator_eligible"] = "False"
        valid, audit = eligible_rows(rows)
        self.assertEqual(len(audit), len(rows))
        self.assertEqual(len(valid), len(rows) - 2)
        self.assertIn("nonpositive_or_missing_dpcr", audit[0]["correction_exclusion"])
        self.assertEqual(rows[0]["dpcr_extract_copies_ul"], "0")

    def test_duplicate_keys_fail(self):
        rows = fixture()
        with self.assertRaises(ValueError):
            eligible_rows(rows + [copy.deepcopy(rows[0])])

    def test_three_parents_do_not_produce_a_four_parent_validation(self):
        result = analyze(
            [row for row in fixture() if row["parent"] != "P3"], bootstrap=10
        )
        self.assertFalse(result["predictions"])
        self.assertTrue(
            all(
                row["validation_status"] == "insufficient_validation"
                for row in result["parameters"]
            )
        )

    def test_cli_round_trip_preserves_measurements(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pairs = root / "pairs.tsv"
            write_rows(pairs, fixture())
            self.assertEqual(
                main(
                    [
                        "fit-correction",
                        "--pairs",
                        str(pairs),
                        "--outdir",
                        str(root / "fit"),
                        "--bootstrap",
                        "20",
                    ]
                ),
                0,
            )
            self.assertEqual(
                main(
                    [
                        "apply-correction",
                        "--pairs",
                        str(pairs),
                        "--parameters",
                        str(root / "fit/correction_parameters.tsv"),
                        "--fit-scope",
                        "pooled",
                        "--output",
                        str(root / "corrected.tsv"),
                    ]
                ),
                0,
            )
            with (root / "corrected.tsv").open() as handle:
                corrected = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["seq_extract_copies_ul"] for row in corrected],
                [row["seq_extract_copies_ul"] for row in fixture()],
            )
            self.assertTrue(
                all(
                    row["correction_status"] == "evaluated_coefficient"
                    for row in corrected
                )
            )

    def test_default_and_optional_historical_slope_gate(self):
        args = parser().parse_args(["run", "--samples", "samples.tsv"])
        settings = settings_from_args(args)
        self.assertEqual(settings.min_r2, 0.80)
        self.assertFalse(settings.enforce_yield_slope_range)
        args = parser().parse_args(
            [
                "run",
                "--samples",
                "samples.tsv",
                "--min-r2",
                ".95",
                "--enforce-yield-slope-range",
            ]
        )
        settings = settings_from_args(args)
        self.assertEqual(settings.min_r2, 0.95)
        self.assertTrue(settings.enforce_yield_slope_range)

    def test_model_acceptance_boundary_and_nonpositive_slope(self):
        points = [
            dict(
                fit_used=True,
                base_fraction=0.01 * n,
                mass_fraction=0.001 * n,
                normalized_depth_per_ul=float(n),
                extract_equivalent_copies_ul=float(n),
                level=n,
                sequencing_yield=10.0,
                ercc_id=f"SIM{n}",
            )
            for n in range(1, 5)
        ]
        data = dict(points=points, total_bases=1000)
        for method in ["sequencing_yield", "depth_calibration"]:
            for r2, slope, accepted in [
                (0.80, 1.0, True),
                (0.799, 1.0, False),
                (0.99, 0.0, False),
                (0.99, -1.0, False),
            ]:
                with mock.patch(
                    "hybcapseq_quant.core.regress",
                    return_value=dict(r2=r2, slope=slope, intercept=0),
                ):
                    result = calibrate(data, method, Settings())[0]
                self.assertEqual(result["passes"], accepted)
