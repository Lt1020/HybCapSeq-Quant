"""Grouped, out-of-sample assessment of target-specific multiplicative bias.

This module does not fit detection thresholds or infer missing concentrations.
Concentrations must refer to the same specimen and measurement layer.
"""

import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics as st
from collections import defaultdict

from . import __version__

REQUIRED = {
    "batch",
    "sample_id",
    "parent",
    "extraction_id",
    "method",
    "target",
    "seq_extract_copies_ul",
    "dpcr_extract_copies_ul",
    "complete_in_range",
    "dpcr_comparator_eligible",
    "model_accepted",
}


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError(
                "Missing correction fields: "
                + ", ".join(sorted(REQUIRED - set(reader.fieldnames or [])))
            )
        return list(reader)


def write_rows(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "NA" if row.get(key) is None else row.get(key) for key in fields}
            )


def numeric(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def yes(value):
    return str(value).strip().lower() in {"true", "1"}


def parent_key(row):
    return row["batch"] + "|" + row["parent"]


def eligible_rows(rows):
    accepted, audit = [], []
    seen = set()
    for row in rows:
        key = (row["batch"], row["sample_id"], row["method"], row["target"])
        if key in seen:
            raise ValueError("Duplicate sample/method/target: " + str(key))
        seen.add(key)
        x, y = numeric(row["seq_extract_copies_ul"]), numeric(
            row["dpcr_extract_copies_ul"]
        )
        reasons = []
        if not yes(row["dpcr_comparator_eligible"]):
            reasons.append("ineligible_dpcr_comparator")
        if not yes(row["model_accepted"]):
            reasons.append("calibration_not_accepted")
        if not yes(row["complete_in_range"]):
            reasons.append("incomplete_or_outside_calibration_range")
        if x is None or x <= 0:
            reasons.append("nonpositive_or_missing_sequencing")
        if y is None or y <= 0:
            reasons.append("nonpositive_or_missing_dpcr")
        item = dict(
            row,
            correction_input_eligible=not reasons,
            correction_exclusion=";".join(reasons) or "eligible",
        )
        audit.append(item)
        if not reasons:
            accepted.append(dict(row, x=x, y=y))
    return accepted, audit


def parent_log_ratios(rows):
    # Collapse repeated ERCC dilutions within extraction and method, then give
    # methods/extractions equal weight within a parent and parents equal weight.
    extraction = defaultdict(list)
    for row in rows:
        extraction[(parent_key(row), row["extraction_id"], row["method"])].append(
            math.log10(row["y"] / row["x"])
        )
    parent = defaultdict(list)
    for (key, _, _), values in extraction.items():
        parent[key].append(st.median(values))
    return {key: st.median(values) for key, values in parent.items()}


def factor(rows):
    logs = parent_log_ratios(rows)
    return 10 ** st.median(logs.values()) if logs else None


def factor_interval(rows, repeats=1000):
    logs = list(parent_log_ratios(rows).values())
    if len(logs) < 2 or repeats < 2:
        return None, None
    seed = int(hashlib.sha256(repr(sorted(logs)).encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)
    draws = sorted(
        10 ** st.median(rng.choices(logs, k=len(logs))) for _ in range(repeats)
    )
    return draws[int((repeats - 1) * 0.025)], draws[int((repeats - 1) * 0.975)]


def metrics(rows, concentration="corrected_concentration"):
    if not rows:
        return {"n_pairs": 0, "n_parents": 0}
    x = [math.log10(row["y"]) for row in rows]
    y = [math.log10(row[concentration]) for row in rows]
    ape = [100 * abs(row[concentration] / row["y"] - 1) for row in rows]
    error = [b - a for a, b in zip(x, y)]
    ccc = r2 = None
    if len(rows) > 1:
        vx, vy = st.variance(x), st.variance(y)
        cov = sum((a - st.mean(x)) * (b - st.mean(y)) for a, b in zip(x, y)) / (
            len(x) - 1
        )
        denominator = vx + vy + (st.mean(x) - st.mean(y)) ** 2
        ccc = 2 * cov / denominator if denominator else None
        r2 = cov**2 / (vx * vy) if vx > 0 and vy > 0 else None
    return dict(
        n_pairs=len(rows),
        n_parents=len({parent_key(row) for row in rows}),
        median_APE_pct=st.median(ape),
        max_APE_pct=max(ape),
        fraction_APE_within_30=sum(value <= 30 + 1e-10 for value in ape) / len(ape),
        mean_log10_bias=st.mean(error),
        rmse_log10=math.sqrt(st.mean(value**2 for value in error)),
        Lin_CCC_log10=ccc,
        r2_log10=r2,
    )


def feasibility(rows):
    # Exact necessary and sufficient intersection for a fixed positive k.
    lower = max((0.7 * row["y"] / row["x"] for row in rows), default=None)
    upper = min((1.3 * row["y"] / row["x"] for row in rows), default=None)
    return dict(
        all_pairs_30pct_feasible=lower <= upper if rows else None,
        feasible_factor_lower=lower,
        feasible_factor_upper=upper,
    )


def prediction(row, coefficient, scheme, training, held_out, scope, profile):
    train_ids = sorted({parent_key(value) for value in training})
    if parent_key(row) in train_ids:
        raise ValueError("Training/validation parent leakage")
    return dict(
        row,
        fit_scope=scope,
        correction_profile=profile,
        validation_scheme=scheme,
        held_out=held_out,
        training_parent_ids=";".join(train_ids),
        training_parent_count=len(train_ids),
        correction_factor=coefficient,
        original_concentration=row["x"],
        corrected_concentration=row["x"] * coefficient,
        before_APE_pct=100 * abs(row["x"] / row["y"] - 1),
        after_APE_pct=100 * abs(row["x"] * coefficient / row["y"] - 1),
    )


def analyze(rows, min_training_parents=3, bootstrap=1000):
    usable, audit = eligible_rows(rows)
    parameters, feasibility_rows, predictions, gaps = [], [], [], []
    batches = sorted({row["batch"] for row in rows})
    targets = sorted({row["target"] for row in rows})
    methods = sorted({row["method"] for row in rows})
    profiles = [(method, {method}) for method in methods]
    profiles.append(("shared_depth_yield", {"depth_calibration", "sequencing_yield"}))
    for target in targets:
        for profile, members in profiles:
            all_rows = [
                row
                for row in usable
                if row["target"] == target and row["method"] in members
            ]
            for scope in ["pooled", *batches]:
                group = [
                    row
                    for row in all_rows
                    if scope == "pooled" or row["batch"] == scope
                ]
                coefficient = factor(group)
                parents = sorted({parent_key(row) for row in group})
                lo, hi = factor_interval(group, bootstrap)
                item = dict(
                    target=target,
                    correction_profile=profile,
                    fit_scope=scope,
                    factor=coefficient,
                    factor_cluster_CI_low=lo,
                    factor_cluster_CI_high=hi,
                    n_training_pairs=len(group),
                    n_training_parents=len(parents),
                    training_parent_ids=";".join(parents),
                    concentration_layer="extract_copies_per_ul",
                    unit="copies/uL",
                    training_min_concentration=min(
                        (row["x"] for row in group), default=None
                    ),
                    training_max_concentration=max(
                        (row["x"] for row in group), default=None
                    ),
                    coefficient_status="descriptive_fit" if group else "not_estimable",
                )
                feasibility_rows.append(
                    dict(
                        target=target,
                        correction_profile=profile,
                        fit_scope=scope,
                        n_pairs=len(group),
                        **feasibility(group),
                    )
                )
                parameters.append(item)
                for held_out in parents:
                    train = [row for row in group if parent_key(row) != held_out]
                    test = [row for row in group if parent_key(row) == held_out]
                    count = len({parent_key(row) for row in train})
                    if count < min_training_parents:
                        gaps.append(
                            dict(
                                target=target,
                                correction_profile=profile,
                                fit_scope=scope,
                                held_out=held_out,
                                reason="insufficient_training_parents",
                                n_training_parents=count,
                                n_test_pairs=len(test),
                            )
                        )
                        continue
                    k = factor(train)
                    predictions.extend(
                        prediction(
                            row,
                            k,
                            "leave_one_parent_out",
                            train,
                            held_out,
                            scope,
                            profile,
                        )
                        for row in test
                    )
            for train_batch in batches:
                train = [row for row in all_rows if row["batch"] == train_batch]
                if len({parent_key(row) for row in train}) < min_training_parents:
                    continue
                k = factor(train)
                test = [row for row in all_rows if row["batch"] != train_batch]
                predictions.extend(
                    prediction(
                        row, k, "cross_batch", train, row["batch"], train_batch, profile
                    )
                    for row in test
                )
    grouped = defaultdict(list)
    for row in predictions:
        grouped[
            (
                row["target"],
                row["correction_profile"],
                row["fit_scope"],
                row["validation_scheme"],
                row["batch"],
                row["method"],
            )
        ].append(row)
    summaries = []
    for key, group in sorted(grouped.items()):
        meta = dict(
            zip(
                [
                    "target",
                    "correction_profile",
                    "fit_scope",
                    "validation_scheme",
                    "test_batch",
                    "method",
                ],
                key,
            )
        )
        for stage, field in [
            ("before", "original_concentration"),
            ("after", "corrected_concentration"),
        ]:
            summaries.append(dict(**meta, stage=stage, **metrics(group, field)))
    for item in parameters:
        holdout = [
            row
            for row in predictions
            if row["target"] == item["target"]
            and row["correction_profile"] == item["correction_profile"]
            and row["fit_scope"] == item["fit_scope"]
            and row["validation_scheme"] == "leave_one_parent_out"
        ]
        m = metrics(holdout)
        sufficient = (
            m["n_pairs"] == item["n_training_pairs"]
            and m["n_parents"] >= min_training_parents + 1
        )
        passed = sufficient and m.get("max_APE_pct", math.inf) <= 30 + 1e-10
        item.update(
            heldout_n_pairs=m["n_pairs"],
            heldout_n_parents=m["n_parents"],
            heldout_median_APE_pct=m.get("median_APE_pct"),
            heldout_max_APE_pct=m.get("max_APE_pct"),
            heldout_fraction_APE_within_30=m.get("fraction_APE_within_30"),
            validation_status=(
                "all_heldout_pairs_within_30pct"
                if passed
                else (
                    "heldout_target_not_met"
                    if sufficient
                    else "insufficient_validation"
                )
            ),
            formal_assay_validation=False,
        )
    return dict(
        parameters=parameters,
        feasibility=feasibility_rows,
        predictions=predictions,
        summary=summaries,
        input_audit=audit,
        validation_gaps=gaps,
    )


def fit_command(args):
    if args.min_training_parents < 3 or args.bootstrap < 2:
        raise ValueError(
            "Require at least three training parents and two bootstrap draws"
        )
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError("Correction output directory must be empty")
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = analyze(read_rows(args.pairs), args.min_training_parents, args.bootstrap)
    for name, rows in data.items():
        # Input x/y are helper copies; retain original source columns unchanged.
        clean = [
            {key: value for key, value in row.items() if key not in {"x", "y"}}
            for row in rows
        ]
        write_rows(args.outdir / ("correction_" + name + ".tsv"), clean)
    manifest = dict(
        software_version=__version__,
        input_sha256=hashlib.sha256(args.pairs.read_bytes()).hexdigest(),
        command_line=args.command_line,
        min_training_parents=args.min_training_parents,
        bootstrap_repeats=args.bootstrap,
        error_target_pct=30,
        estimator="median_parent_log_ratio_after_extraction_method_collapse",
        scope="exploratory_bias_assessment_not_formal_assay_validation",
        n_input_rows=len(data["input_audit"]),
        n_validation_predictions=len(data["predictions"]),
    )
    (args.outdir / "correction_run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return 0


def apply_command(args):
    with args.parameters.open(encoding="utf-8-sig", newline="") as handle:
        coefficients = list(csv.DictReader(handle, delimiter="\t"))
    selected = {}
    for row in coefficients:
        if (
            row["fit_scope"] != args.fit_scope
            or row["correction_profile"] == "shared_depth_yield"
        ):
            continue
        key = row["target"], row["correction_profile"]
        if key in selected:
            raise ValueError("Duplicate correction coefficient: " + str(key))
        selected[key] = row
    source = read_rows(args.pairs)
    output = []
    for row in source:
        parameter = selected.get((row["target"], row["method"]))
        value = numeric(row["seq_extract_copies_ul"])
        coefficient = numeric(parameter.get("factor")) if parameter else None
        status = "missing_coefficient"
        corrected = None
        if parameter and coefficient is not None:
            lo, hi = numeric(parameter.get("training_min_concentration")), numeric(
                parameter.get("training_max_concentration")
            )
            accepted = (
                parameter.get("validation_status") == "all_heldout_pairs_within_30pct"
            )
            status = "not_validated_for_30pct"
            if value is None or value <= 0:
                status = "nonpositive_or_missing_input"
            elif not yes(row["model_accepted"]) or not yes(row["complete_in_range"]):
                status = "invalid_or_out_of_range_source_estimate"
            elif lo is None or hi is None or not lo <= value <= hi:
                status = "outside_correction_training_range"
            elif accepted or args.allow_exploratory:
                corrected = value * coefficient
                status = (
                    "evaluated_coefficient" if accepted else "exploratory_coefficient"
                )
        output.append(
            dict(
                row,
                correction_factor=coefficient,
                corrected_extract_copies_ul=corrected,
                correction_status=status,
                correction_fit_scope=args.fit_scope,
            )
        )
    if args.output.exists():
        raise ValueError("Correction output file already exists")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_rows(args.output, output)
    provenance = dict(
        software_version=__version__,
        pairs_sha256=hashlib.sha256(args.pairs.read_bytes()).hexdigest(),
        parameters_sha256=hashlib.sha256(args.parameters.read_bytes()).hexdigest(),
        fit_scope=args.fit_scope,
        allow_exploratory=args.allow_exploratory,
        command_line=args.command_line,
        input_values_preserved=True,
    )
    args.output.with_suffix(args.output.suffix + ".manifest.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return 0


def add_commands(sub):
    fit = sub.add_parser(
        "fit-correction",
        help="assess per-target multiplicative bias with grouped validation",
    )
    fit.add_argument("--pairs", type=Path, required=True)
    fit.add_argument("--outdir", type=Path, required=True)
    fit.add_argument("--min-training-parents", type=int, default=3)
    fit.add_argument("--bootstrap", type=int, default=1000)
    fit.set_defaults(func=fit_command)
    apply = sub.add_parser(
        "apply-correction",
        help="append correction estimates without changing source values",
    )
    apply.add_argument("--pairs", type=Path, required=True)
    apply.add_argument("--parameters", type=Path, required=True)
    apply.add_argument(
        "--fit-scope", required=True, help="pooled or an exact training batch label"
    )
    apply.add_argument("--allow-exploratory", action="store_true")
    apply.add_argument("--output", type=Path, required=True)
    apply.set_defaults(func=apply_command)
