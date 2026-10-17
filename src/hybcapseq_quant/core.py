"""Equations 5–22 of the supplied protocol, with explicit concentration layers.

No RPKMF, covered-bases breadth, or nominal read length is substituted for
aligned bases. All input detection rows survive, in their original order.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics as stats
from dataclasses import dataclass
from pathlib import Path

NA = 6.02214076e23
MISSING = {"", "NA", "N/A", "NAN", "NONE", "NULL", "/"}
ERCC_PATTERN = re.compile(r"(?<![A-Za-z0-9])ERCC[-_](\d{5})(?!\d)", re.I)
METHODS = ("method1_yield", "method2_depth")
CONC = "浓度concentration_ERCC"
WWCONC = "污水浓度concentration_wastewater"
STANDARD_PATH = Path(__file__).parent / "data" / "ercc92_mix1.tsv"


def number(value, *, name="value", positive=False, nonnegative=False):
    if value is None or str(value).strip().upper() in MISSING:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}: invalid number {value!r}") from exc
    if not math.isfinite(v) or (positive and v <= 0) or (nonnegative and v < 0):
        raise ValueError(f"{name}: invalid value {value!r}")
    return v


def flag(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate TSV column: {path}")
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"Too many fields; input must be tab-delimited: {path}")
    return fields, rows


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".12g") if math.isfinite(value) else "NA"
    return str(value)


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k)) for k in fields})


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as src:
        while block := src.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def ercc_id(row):
    ids = set()
    for key in ("Accession", "Name", "Assembly"):
        ids.update("ERCC-" + m for m in ERCC_PATTERN.findall(str(row.get(key, ""))))
    if len(ids) > 1:
        raise ValueError(f"Conflicting ERCC identifiers: {row}")
    return next(iter(ids), None)


def ercc_amounts(standards, sample):
    dilution = number(sample.get("ercc_dilution"), name="ercc_dilution", positive=True)
    spike = number(sample.get("ercc_spike_ul"), name="ercc_spike_ul", nonnegative=True)
    reaction = number(sample.get("reaction_volume_ul"), name="reaction_volume_ul", positive=True)
    vlib = number(sample.get("library_extract_ul"), name="library_extract_ul", positive=True)
    if None in (dilution, spike, reaction, vlib):
        raise ValueError("ercc_dilution, ercc_spike_ul, reaction_volume_ul and library_extract_ul are required")
    if vlib + spike > reaction + 1e-8:
        raise ValueError("Extract plus ERCC volumes exceed starting reaction volume")
    if sample.get("ercc_mix", "Mix1") not in {"Mix1", "Mix 1"}:
        raise ValueError("This concentration file is for commercial ERCC Mix1")
    result = {}
    for item in standards:
        eid = item["ERCC_ID"]
        stock = number(item["mix1_attomol_ul"], name=eid, positive=True) * NA * 1e-18
        length = number(item["Length"], name=eid + " Length", positive=True)
        copies = stock / dilution * spike
        result[eid] = {
            "ERCC_ID": eid, "Length": length,
            "level": int(item["level"]), "stock_attomol_ul": float(item["mix1_attomol_ul"]),
            "stock_copies_ul": stock, "working_copies_ul": stock / dilution,
            "added_copies": copies, "reaction_copies_ul": copies / reaction,
            "extract_equivalent_copies_ul": copies / vlib,
            "added_mass_ng": copies * length * 340 / NA * 1e9,
            "fit_include": flag(item.get("fit_include", "1")),
        }
    if len(result) != len(standards):
        raise ValueError("Duplicate ERCC standard ID")
    return result


def regress(x, y):
    if len(x) < 2 or len(set(x)) < 2:
        return None
    xm, ym = stats.fmean(x), stats.fmean(y)
    xx = sum((a - xm) ** 2 for a in x)
    slope = sum((a - xm) * (b - ym) for a, b in zip(x, y)) / xx
    intercept = ym - slope * xm
    residuals = [b - (intercept + slope * a) for a, b in zip(x, y)]
    yy = sum((b - ym) ** 2 for b in y)
    r2 = 1 - sum(v * v for v in residuals) / yy if yy else None
    return {"slope": slope, "intercept": intercept, "r2": r2,
            "n": len(x), "x_min": min(x), "x_max": max(x),
            "equation": f"y = {intercept:.10g} + {slope:.10g} * x"}


def adjacent_run(levels):
    best = current = 0
    previous = None
    for level in sorted(set(levels)):
        current = current + 1 if previous is not None and level == previous + 1 else 1
        best = max(best, current)
        previous = level
    return best


@dataclass
class Settings:
    min_ercc_reads: int = 1
    min_ercc_breadth: float = 0.0
    min_r2: float = 0.95
    min_levels: int = 4
    slope_min: float = 0.8
    slope_max: float = 1.2
    allow_extrapolation: bool = False
    strict_qc: bool = False
    ascii_columns: bool = False
    wastewater_unit: str = "copies/mL"


def resolve_path(value, base):
    if value is None or str(value).strip().upper() in MISSING:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(base) / path


def input_bases(sample, table_path, rows, base):
    supplied = number(sample.get("total_effective_bases"), name="total_effective_bases", positive=True)
    if supplied is not None:
        return supplied, "sample_manifest"
    path = resolve_path(sample.get("fastp_json"), base)
    if path is None:
        path = table_path.with_name(sample["sample_ID"] + ".fastp.json")
    if not path.exists():
        return None, "not_available"
    section = sample.get("fastp_section", "before_filtering")
    if section not in {"before_filtering", "after_filtering"}:
        raise ValueError("fastp_section must be before_filtering or after_filtering")
    data = json.loads(path.read_text())["summary"][section]
    count = number(data.get("total_reads"), name="fastp total_reads", positive=True)
    recorded = {number(r.get("filtered_reads_in_sample"), positive=True) for r in rows}
    recorded.discard(None)
    if recorded and (len(recorded) != 1 or count not in recorded):
        raise ValueError(f"fastp {section} read count does not match EsViritu filtered_reads_in_sample: {path}")
    return number(data.get("total_bases"), name="fastp total_bases", positive=True), f"{path}:{section}"


def load_sample(sample, table_path, standards, base, settings, ercc_overrides=None):
    fields, rows = read_tsv(table_path)
    required = {"sample_ID", "Name", "Length", "Accession", "RPKMF", "read_count", "covered_bases", "mean_coverage"}
    if not required.issubset(fields):
        raise ValueError(f"Missing detection columns: {sorted(required - set(fields))}: {table_path}")
    seen = set()
    for row in rows:
        if row["sample_ID"] != sample["sample_ID"]:
            raise ValueError(f"sample_ID does not match manifest: {table_path}")
        acc = row["Accession"]
        if acc in seen:
            raise ValueError(f"Duplicate reference would be double-counted: {acc}")
        seen.add(acc)
        length = number(row["Length"], positive=True)
        depth = number(row["mean_coverage"], nonnegative=True)
        breadth = number(row["covered_bases"], nonnegative=True)
        reads = number(row["read_count"], nonnegative=True)
        if None in (length, depth, breadth, reads) or breadth > length:
            raise ValueError(f"Invalid alignment statistics: {acc}")
    amounts = ercc_amounts(standards, sample)
    mass = number(sample.get("sample_input_mass_ng"), name="sample_input_mass_ng", positive=True)
    cna = number(sample.get("nucleic_acid_conc_ng_ul"), name="nucleic_acid_conc_ng_ul", positive=True)
    vlib = number(sample["library_extract_ul"], positive=True)
    if mass is None or cna is None:
        raise ValueError("sample_input_mass_ng and nucleic_acid_conc_ng_ul are required (RNA + DNA)")
    if not math.isclose(mass, cna * vlib, rel_tol=1e-5, abs_tol=1e-7):
        raise ValueError("sample_input_mass_ng != nucleic_acid_conc_ng_ul * library_extract_ul")
    total_mass = mass + sum(v["added_mass_ng"] for v in amounts.values())
    total_bases, total_source = input_bases(sample, table_path, rows, base)
    evidence_path = resolve_path(sample.get("alignment_metrics_tsv"), base)
    metrics = {}
    if evidence_path:
        _, observations = read_tsv(evidence_path)
        for item in observations:
            if item["Accession"] in metrics:
                raise ValueError("Duplicate alignment metrics accession")
            metrics[item["Accession"]] = item
        if set(seen) - set(metrics):
            raise ValueError("Alignment metrics missing detected references")
    signals = {}
    for row in rows:
        acc, length = row["Accession"], float(row["Length"])
        if metrics:
            item = metrics[acc]
            if not math.isclose(float(item["Length"]), length):
                raise ValueError(f"Alignment/reference length mismatch: {acc}")
            n = number(item["aligned_bases"], nonnegative=True)
            if n is None:
                raise ValueError("Missing aligned_bases in metrics")
            depth = n / length
            covered = number(item["covered_bases"], nonnegative=True)
            reads = number(item["read_count"], nonnegative=True)
        else:
            depth = float(row["mean_coverage"])
            n = depth * length
            covered, reads = float(row["covered_bases"]), float(row["read_count"])
        if covered is None or reads is None or covered > length:
            raise ValueError(f"Invalid alignment metrics: {acc}")
        signals[acc] = {"aligned_bases": n, "depth": depth, "breadth": covered / length,
                        "covered_bases": covered, "reads": reads}
    points = []
    observed_ercc = {}
    for row in rows:
        eid = ercc_id(row)
        if eid:
            if eid in observed_ercc:
                raise ValueError(f"Multiple references for the same ERCC: {eid}")
            observed_ercc[eid] = row
    for eid, amount in amounts.items():
        row = observed_ercc.get(eid)
        point = dict(amount)
        point.update({"sample_ID": sample["sample_ID"], "detected": row is not None,
                      "aligned_bases": 0.0, "mean_depth": 0.0, "breadth": 0.0,
                      "read_count": 0.0, "fit_exclusion": "not_detected", "fit_used": False,
                      "mass_fraction": amount["added_mass_ng"] / total_mass,
                      "base_fraction": None, "Y_seq_i": None,
                      "O_normalized": None, "replicate_Y_cv": None})
        if row:
            if float(row["Length"]) != amount["Length"]:
                raise ValueError(f"ERCC length differs from v3.3 reference: {eid}")
            signal = signals[row["Accession"]]
            reasons = []
            if amount["added_copies"] <= 0:
                reasons.append("no_ERCC_spiked")
            if signal["aligned_bases"] <= 0:
                reasons.append("zero_aligned_bases")
            if signal["reads"] < settings.min_ercc_reads:
                reasons.append("below_min_ercc_reads")
            if signal["breadth"] < settings.min_ercc_breadth:
                reasons.append("below_min_ercc_breadth")
            if not amount["fit_include"]:
                reasons.append("excluded_in_standards")
            override = (ercc_overrides or {}).get((sample["sample_ID"], eid))
            if override and not flag(override.get("use_for_fit")):
                if not override.get("reason"):
                    raise ValueError("An ERCC exclusion requires a reason")
                reasons.append("manual:" + override["reason"])
            loq = number(sample.get("loq_reaction_copies_ul"), nonnegative=True)
            if loq is not None and amount["reaction_copies_ul"] < loq:
                reasons.append("below_supplied_LoQ")
            bf = signal["aligned_bases"] / total_bases if total_bases else None
            mf = amount["added_mass_ng"] / total_mass
            point.update({"aligned_bases": signal["aligned_bases"], "mean_depth": signal["depth"],
                          "breadth": signal["breadth"], "read_count": signal["reads"],
                          "fit_exclusion": ";".join(reasons), "fit_used": not reasons,
                          "base_fraction": bf, "Y_seq_i": bf / mf if bf is not None and mf > 0 else None,
                          "O_normalized": signal["depth"] * cna / mass})
        points.append(point)
    return {"sample": sample, "path": table_path, "fields": fields, "rows": rows,
            "amounts": amounts, "points": points, "signals": signals,
            "total_mass": total_mass, "total_bases": total_bases, "total_bases_source": total_source,
            "input_sha256": sha256(table_path), "alignment_metrics_sha256": sha256(evidence_path) if evidence_path else None,
            "evidence_basis": "alignment_metrics" if metrics else "EsViritu_full_reference_mean_depth_times_Length"}


def replicate_qc(samples):
    groups = {}
    for data in samples:
        group = data["sample"].get("replicate_group", "NA")
        if group.upper() in MISSING:
            continue
        for point in data["points"]:
            if point["fit_used"] and point["Y_seq_i"] is not None and point["Y_seq_i"] > 0:
                groups.setdefault((data["sample"].get("batch_ID", "NA"), group, point["ERCC_ID"]), []).append(point)
    results = []
    for (batch, group, eid), points in groups.items():
        ys = [p["Y_seq_i"] for p in points]
        cv = stats.stdev(ys) / stats.fmean(ys) if len(ys) >= 2 else None
        assessed = len(ys) >= 3
        for p in points:
            p["replicate_Y_cv"] = cv
            if assessed and cv > 0.35:
                p["fit_used"] = False
                p["fit_exclusion"] = "replicate_Y_CV_above_35pct"
        results.append({"batch_ID": batch, "replicate_group": group, "ERCC_ID": eid, "n": len(ys),
                        "Y_seq_cv": cv, "three_replicates_assessed": assessed,
                        "passes": assessed and cv <= 0.35})
    return results


def calibrate(data, method, settings):
    eligible = [p for p in data["points"] if p["fit_used"]]
    errors = []
    if method == "method1_yield":
        eligible = [p for p in eligible if p["base_fraction"] is not None and p["base_fraction"] > 0]
        x = [math.log10(p["mass_fraction"]) for p in eligible]
        y = [math.log10(p["base_fraction"]) for p in eligible]
        if data["total_bases"] is None:
            errors.append("missing_total_effective_bases")
        axes = ("log10(ERCC input mass / total input mass)", "log10(ERCC aligned bases / total effective bases)")
    else:
        eligible = [p for p in eligible if p["O_normalized"] is not None and p["O_normalized"] > 0]
        x = [math.log10(p["O_normalized"]) for p in eligible]
        y = [math.log10(p["extract_equivalent_copies_ul"]) for p in eligible]
        axes = ("log10(depth * mixed-NA concentration / sample input mass)", "log10(known ERCC extract-equivalent copies/uL)")
    fit = regress(x, y)
    run = adjacent_run(p["level"] for p in eligible)
    if run < settings.min_levels:
        errors.append("fewer_than_4_adjacent_levels" if settings.min_levels == 4 else "insufficient_adjacent_levels")
    if fit is None:
        errors.append("insufficient_regression_points")
    else:
        if fit["r2"] is None or fit["r2"] < settings.min_r2:
            errors.append("R2_below_threshold")
        if fit["slope"] <= 0:
            errors.append("nonpositive_slope")
        if method == "method1_yield" and not settings.slope_min <= fit["slope"] <= settings.slope_max:
            errors.append("yield_slope_outside_protocol_range")
    yields = [p["Y_seq_i"] for p in eligible] if method == "method1_yield" else []
    result = {"method": method, "fit": fit, "passes": not errors, "errors": errors,
              "n_fitted": len(eligible), "distinct_levels": len(set(p["level"] for p in eligible)),
              "longest_adjacent_levels": run, "x_label": axes[0], "y_label": axes[1],
              "Y_seq_mean": stats.fmean(yields) if yields else None,
              "Y_seq_median": stats.median(yields) if yields else None,
              "fit_ERCC_IDs": [p["ERCC_ID"] for p in eligible]}
    return result, eligible, x, y


def wastewater_factor(sample, unit):
    if sample.get("sample_type", "wastewater") != "wastewater":
        return None, {"status": "not_a_wastewater_sample", "unit": unit}
    names = ("extract_elution_ul", "extract_input_concentrate_ul", "concentrate_volume_ml", "wastewater_volume_ml")
    v = {k: number(sample.get(k), name=k, positive=True) for k in names}
    recovery = number(sample.get("pmmov_recovery_fraction"), name="pmmov_recovery_fraction", positive=True)
    audit = {"unit": unit, **v}
    if recovery is None:
        q0 = number(sample.get("pmmov_raw_extract_copies_ul"), nonnegative=True)
        q1 = number(sample.get("pmmov_post_extract_copies_ul"), nonnegative=True)
        d0 = number(sample.get("pmmov_raw_dilution"), positive=True)
        d1 = number(sample.get("pmmov_post_dilution"), positive=True)
        ve0 = number(sample.get("raw_extract_elution_ul"), positive=True)
        vx0 = number(sample.get("raw_extract_input_ul"), positive=True)
        loq_pass = flag(sample.get("pmmov_above_loq"))
        if None not in (q0, q1, d0, d1, ve0, vx0, *v.values()) and q0 > 0 and q1 > 0 and loq_pass:
            c0 = 1000 * q0 * d0 * ve0 / vx0
            c1 = 1000 * q1 * d1 * v["extract_elution_ul"] / v["extract_input_concentrate_ul"]
            recovery = c1 * v["concentrate_volume_ml"] / (c0 * v["wastewater_volume_ml"])
            audit.update({"pmmov_raw_copies_ml": c0, "pmmov_concentrate_copies_ml": c1,
                          "recovery_source": "protocol_equations_5_to_8"})
        else:
            audit["recovery_source"] = "not_available_or_PMMoV_below_LoQ"
    else:
        audit["recovery_source"] = "sample_manifest"
    audit["pmmov_recovery_fraction"] = recovery
    missing = [k for k, value in v.items() if value is None]
    if recovery is None:
        missing.append("pmmov_recovery_fraction_or_complete_dPCR_inputs")
    if missing:
        audit["status"] = "missing:" + ",".join(missing)
        return None, audit
    if not 0.01 <= recovery <= 1.50:
        audit["status"] = "recovery_outside_protocol_review_range"
        return None, audit
    multiplier = 1000 * v["extract_elution_ul"] / v["extract_input_concentrate_ul"]
    multiplier *= v["concentrate_volume_ml"] / (v["wastewater_volume_ml"] * recovery)
    if unit == "copies/L":
        multiplier *= 1000
    audit.update({"status": "available_PMMoV_proxy_correction", "extract_to_wastewater_factor": multiplier})
    return multiplier, audit


def reference_groups(data, targets):
    groups = {}
    for row in data["rows"]:
        if ercc_id(row):
            continue
        acc = row["Accession"]
        target = targets.get(acc, {})
        assembly = row.get("Assembly", "NA")
        default_id = assembly if str(assembly).strip().upper() not in MISSING else acc
        group_id = target.get("quant_group") or default_id
        if str(group_id).upper() in MISSING:
            group_id = default_id
        groups.setdefault(group_id, []).append((row, target))
    result = {}
    for gid, items in groups.items():
        lengths = set()
        masses = set()
        for row, target in items:
            length = number(target.get("quant_reference_length"), positive=True)
            if length is None:
                length = number(row.get("Asm_length"), positive=True) or float(row["Length"])
            lengths.add(length)
            mass = number(target.get("mass_per_nt_g_mol"), positive=True)
            if mass is not None:
                masses.add(mass)
        if len(lengths) != 1 or len(masses) > 1:
            raise ValueError(f"Inconsistent quantification reference metadata for {gid}")
        length = next(iter(lengths))
        if sum(float(r["Length"]) for r, _ in items) > length + 1e-6:
            raise ValueError(f"Observed reference lengths exceed locked genome length: {gid}")
        n = sum(data["signals"][r["Accession"]]["aligned_bases"] for r, _ in items)
        covered = sum(data["signals"][r["Accession"]]["covered_bases"] for r, _ in items)
        item = {"quant_group": gid, "quant_reference_length": length, "aligned_bases": n,
                "mean_depth": n / length, "breadth": covered / length,
                "mass_per_nt_g_mol": next(iter(masses), None),
                "detected_accessions": ",".join(r["Accession"] for r, _ in items),
                "detected_reference_length": sum(float(r["Length"]) for r, _ in items)}
        for row, _ in items:
            result[row["Accession"]] = item
    return result


def quantify_rows(data, calibration, targets, settings):
    sample, method = data["sample"], calibration["method"]
    vlib = float(sample["library_extract_ul"])
    reaction = float(sample["reaction_volume_ul"])
    cna = float(sample["nucleic_acid_conc_ng_ul"])
    mass = float(sample["sample_input_mass_ng"])
    factor, recovery_audit = wastewater_factor(sample, settings.wastewater_unit)
    groups = reference_groups(data, targets)
    results, group_results = [], {}
    for row in data["rows"]:
        eid = ercc_id(row)
        copy = dict(row)
        status = []
        ext = reaction_c = ww = estimated = None
        group = None
        if eid:
            known = data["amounts"].get(eid)
            if known and known["added_copies"] > 0:
                reaction_c = known["reaction_copies_ul"]
                status = ["ERCC_known_input"]
            elif known:
                status = ["ERCC_detected_but_not_spiked"]
            else:
                status = ["ERCC_not_in_commercial_Mix1"]
        else:
            group = groups[row["Accession"]]
            if group["detected_reference_length"] < group["quant_reference_length"]:
                status.append("partial_quantification_reference_detected")
            n, length = group["aligned_bases"], group["quant_reference_length"]
            if n <= 0:
                status.append("no_quantifiable_aligned_bases")
            if not calibration["passes"]:
                status.append("calibration_failed:" + ",".join(calibration["errors"]))
            if calibration["passes"] and n > 0:
                if method == "method1_yield":
                    mw = group["mass_per_nt_g_mol"]
                    if mw is None:
                        status.append("missing_molecular_mass_factor")
                    else:
                        inferred_mass = data["total_mass"] * (n / data["total_bases"]) / calibration["Y_seq_mean"]
                        molecules = inferred_mass * 1e-9 * NA / (mw * length)
                        ext = molecules / vlib
                        log_input = math.log10(inferred_mass / data["total_mass"])
                else:
                    observed = group["mean_depth"] * cna / mass
                    log_input = math.log10(observed)
                    fit = calibration["fit"]
                    ext = 10 ** (fit["intercept"] + fit["slope"] * log_input)
                if ext is not None:
                    estimated = ext * vlib / reaction
                    fit = calibration["fit"]
                    outside = log_input < fit["x_min"] - 1e-9 or log_input > fit["x_max"] + 1e-9
                    if outside:
                        status.append("outside_calibration_range")
                    if not outside or settings.allow_extrapolation:
                        reaction_c = estimated
                    loq = number(sample.get("loq_reaction_copies_ul"), nonnegative=True)
                    if loq is not None and estimated < loq:
                        status.append("below_supplied_LoQ")
                        reaction_c = None
                    verified = all(flag(sample.get(k)) for k in ("alignment_qc_verified", "detection_qc_verified"))
                    if not verified:
                        status.append("alignment_or_detection_QC_unverified")
                    if loq is None:
                        status.append("LoQ_not_validated")
                    if settings.strict_qc and (not verified or loq is None):
                        reaction_c = None
                    if reaction_c is not None and factor is not None:
                        ww = ext * factor
                    if not status:
                        status.append("quantified_with_supplied_QC_and_LoQ")
            if factor is None:
                status.append("wastewater_conversion_unavailable")
        copy.update({CONC: reaction_c, WWCONC: ww,
                     "quant_method": method, "concentration_ERCC_unit": "copies/ul_reaction",
                     "concentration_wastewater_unit": settings.wastewater_unit,
                     "estimated_reaction_copies_ul_before_QC": estimated,
                     "concentration_extract_copies_ul": ext,
                     "quant_group": group["quant_group"] if group else eid,
                     "quant_reference_length": group["quant_reference_length"] if group else None,
                     "aligned_bases_for_quant": group["aligned_bases"] if group else data["signals"][row["Accession"]]["aligned_bases"],
                     "quant_status": ";".join(status),
                     "molecular_mass_factor_g_mol_nt": group["mass_per_nt_g_mol"] if group else 340,
                     "pmmov_recovery_fraction": recovery_audit.get("pmmov_recovery_fraction")})
        results.append(copy)
        if group:
            group_results[group["quant_group"]] = {"sample_ID": sample["sample_ID"], **group,
                CONC: reaction_c, WWCONC: ww, "concentration_extract_copies_ul": ext,
                "quant_status": copy["quant_status"], "wastewater_unit": settings.wastewater_unit}
    return results, list(group_results.values()), recovery_audit


def plot_curve(path, x, y, fit, xlabel, ylabel, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.7, 5.2), constrained_layout=True)
    ax.scatter(x, y, s=24, color="#2563A6", alpha=0.75, edgecolors="none")
    if fit:
        ends = [min(x), max(x)]
        ax.plot(ends, [fit["intercept"] + fit["slope"] * v for v in ends], color="#B44B35", lw=1.8)
        r2 = "NA" if fit["r2"] is None else f'{fit["r2"]:.4f}'
        ax.text(0.03, 0.97, f'y = {fit["intercept"]:.4f} + {fit["slope"]:.4f} x\n$R^2$ = {r2}; n = {len(x)}',
                transform=ax.transAxes, va="top", fontsize=10)
    else:
        ax.text(0.5, 0.5, "Insufficient calibration evidence", transform=ax.transAxes, ha="center")
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.spines[["top", "right"]].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(str(path) + "." + ext, dpi=300)
    plt.close(fig)
