"""Mass-yield and normalized-depth calibration with explicit measurement units."""

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
METHODS = ("sequencing_yield", "depth_calibration")
CONC = "concentration_reaction_copies_ul"
WWCONC = "concentration_wastewater"
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
    if value is None or str(value).strip().upper() in MISSING:
        return False
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate TSV column: {path}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(
            f"Inconsistent field count; input must be tab-delimited: {path}"
        )
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
        writer = csv.DictWriter(
            out, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
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
    dilution = number(
        sample.get("ercc_dilution_factor"), name="ercc_dilution_factor", positive=True
    )
    spike = number(
        sample.get("ercc_spike_volume_ul"),
        name="ercc_spike_volume_ul",
        nonnegative=True,
    )
    reaction = number(
        sample.get("reaction_volume_ul"), name="reaction_volume_ul", positive=True
    )
    vlib = number(
        sample.get("library_input_volume_ul"),
        name="library_input_volume_ul",
        positive=True,
    )
    if None in (dilution, spike, reaction, vlib):
        raise ValueError(
            "ercc_dilution_factor, ercc_spike_volume_ul, reaction_volume_ul and library_input_volume_ul are required"
        )
    if vlib + spike > reaction + 1e-8:
        raise ValueError("Extract plus ERCC volumes exceed starting reaction volume")
    if sample.get("ercc_mix") != "Mix1":
        raise ValueError("This concentration file is for commercial ERCC Mix1")
    result = {}
    for item in standards:
        eid = item["ercc_id"]
        stock = number(item["stock_attomol_ul"], name=eid, positive=True) * NA * 1e-18
        length = number(item["length_nt"], name=eid + " length_nt", positive=True)
        copies = stock / dilution * spike
        result[eid] = {
            "ercc_id": eid,
            "length_nt": length,
            "level": int(item["level"]),
            "stock_attomol_ul": float(item["stock_attomol_ul"]),
            "stock_copies_ul": stock,
            "working_copies_ul": stock / dilution,
            "added_copies": copies,
            "reaction_copies_ul": copies / reaction,
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
    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "rmse_log10": math.sqrt(sum(v * v for v in residuals) / len(x)),
        "residual_degrees_of_freedom": len(x) - 2,
        "n": len(x),
        "x_min": min(x),
        "x_max": max(x),
        "equation": f"y = {intercept:.10g} + {slope:.10g} * x",
    }


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
    reporting_mode: str = "validated"
    filter_replicate_cv: bool = False
    max_replicate_cv: float = 0.35
    wastewater_unit: str = "copies/mL"


def resolve_path(value, base):
    if value is None or str(value).strip().upper() in MISSING:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(base) / path


def input_bases(sample, table_path, rows, base):
    supplied = number(
        sample.get("total_effective_bases"), name="total_effective_bases", positive=True
    )
    if supplied is not None:
        return supplied, "sample_manifest"
    path = resolve_path(sample.get("fastp_json"), base)
    if path is None:
        path = table_path.with_name(sample["sample_id"] + ".fastp.json")
    if not path.exists():
        return None, "not_available"
    section = sample.get("fastp_section")
    if section not in {"before_filtering", "after_filtering"}:
        raise ValueError("fastp_section must be before_filtering or after_filtering")
    data = json.loads(path.read_text())["summary"][section]
    count = number(data.get("total_reads"), name="fastp total_reads", positive=True)
    recorded = {number(r.get("filtered_reads_in_sample"), positive=True) for r in rows}
    recorded.discard(None)
    if recorded and (len(recorded) != 1 or count not in recorded):
        raise ValueError(
            f"fastp {section} read count does not match EsViritu filtered_reads_in_sample: {path}"
        )
    return (
        number(data.get("total_bases"), name="fastp total_bases", positive=True),
        f"{path}:{section}",
    )


def load_sample(sample, table_path, standards, base, settings, ercc_overrides=None):
    fields, rows = read_tsv(table_path)
    if sample.get("sample_type") not in {
        "wastewater",
        "negative_control",
        "reference_material",
    }:
        raise ValueError(
            "sample_type must be wastewater, negative_control or reference_material"
        )
    required = {
        "sample_ID",
        "Name",
        "Length",
        "Accession",
        "RPKMF",
        "read_count",
        "covered_bases",
        "mean_coverage",
    }
    if not required.issubset(fields):
        raise ValueError(
            f"Missing detection columns: {sorted(required - set(fields))}: {table_path}"
        )
    if CONC in fields or WWCONC in fields or "quantification_method" in fields:
        raise ValueError("Provide the original EsViritu detection table")
    seen = set()
    for row in rows:
        if row["sample_ID"] != sample["sample_id"]:
            raise ValueError(
                f"Source sample_ID does not match manifest sample_id: {table_path}"
            )
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
    mass = number(
        sample.get("library_input_mass_ng"), name="library_input_mass_ng", positive=True
    )
    cna = number(
        sample.get("extract_nucleic_acid_concentration_ng_ul"),
        name="extract_nucleic_acid_concentration_ng_ul",
        positive=True,
    )
    vlib = number(sample["library_input_volume_ul"], positive=True)
    if mass is None or cna is None:
        raise ValueError(
            "library_input_mass_ng and extract_nucleic_acid_concentration_ng_ul are required (RNA + DNA)"
        )
    if not math.isclose(mass, cna * vlib, rel_tol=1e-5, abs_tol=1e-7):
        raise ValueError(
            "library_input_mass_ng != extract_nucleic_acid_concentration_ng_ul * library_input_volume_ul"
        )
    total_mass = mass + sum(v["added_mass_ng"] for v in amounts.values())
    total_bases, total_source = input_bases(sample, table_path, rows, base)
    evidence_path = resolve_path(sample.get("alignment_metrics_tsv"), base)
    metrics = {}
    if evidence_path:
        _, observations = read_tsv(evidence_path)
        for item in observations:
            if item["accession"] in metrics:
                raise ValueError("Duplicate alignment metrics accession")
            metrics[item["accession"]] = item
        if set(seen) - set(metrics):
            raise ValueError("Alignment metrics missing detected references")
    signals = {}
    for row in rows:
        acc, length = row["Accession"], float(row["Length"])
        if metrics:
            item = metrics[acc]
            if not math.isclose(float(item["length_nt"]), length):
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
        if (
            covered is None
            or reads is None
            or covered > length
            or covered > n + max(1e-6, n * 1e-6)
        ):
            raise ValueError(f"Invalid alignment metrics: {acc}")
        signals[acc] = {
            "aligned_bases": n,
            "depth": depth,
            "breadth": covered / length,
            "covered_bases": covered,
            "reads": reads,
        }
    if total_bases and sum(
        v["aligned_bases"] for v in signals.values()
    ) > total_bases * (1 + 1e-6):
        raise ValueError(
            "Summed aligned bases exceed total_effective_bases; review input stage and alignment accounting"
        )
    fastp_path = resolve_path(sample.get("fastp_json"), base)
    if fastp_path is None:
        fastp_path = table_path.with_name(sample["sample_id"] + ".fastp.json")
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
        point.update(
            {
                "sample_id": sample["sample_id"],
                "detected": row is not None,
                "aligned_bases": None,
                "mean_depth": None,
                "breadth": None,
                "read_count": None,
                "fit_exclusion": "not_detected",
                "fit_used": False,
                "mass_fraction": amount["added_mass_ng"] / total_mass,
                "base_fraction": None,
                "sequencing_yield": None,
                "normalized_depth_per_ul": None,
                "replicate_yield_cv": None,
            }
        )
        if row:
            if float(row["Length"]) != amount["length_nt"]:
                raise ValueError(f"ERCC length differs from v3.3 reference: {eid}")
            signal = signals[row["Accession"]]
            reasons = []
            if amount["added_copies"] <= 0:
                reasons.append("no_ercc_spiked")
            if signal["aligned_bases"] <= 0:
                reasons.append("zero_aligned_bases")
            if signal["reads"] < settings.min_ercc_reads:
                reasons.append("below_minimum_ercc_reads")
            if signal["breadth"] < settings.min_ercc_breadth:
                reasons.append("below_min_ercc_breadth")
            if not amount["fit_include"]:
                reasons.append("excluded_in_standards")
            override = (ercc_overrides or {}).get((sample["sample_id"], eid))
            if override and not flag(override.get("use_for_fit")):
                if not override.get("reason"):
                    raise ValueError("An ERCC exclusion requires a reason")
                reasons.append("manual:" + override["reason"])
            bf = signal["aligned_bases"] / total_bases if total_bases else None
            mf = amount["added_mass_ng"] / total_mass
            point.update(
                {
                    "aligned_bases": signal["aligned_bases"],
                    "mean_depth": signal["depth"],
                    "breadth": signal["breadth"],
                    "read_count": signal["reads"],
                    "fit_exclusion": ";".join(reasons),
                    "fit_used": not reasons,
                    "base_fraction": bf,
                    "sequencing_yield": bf / mf if bf is not None and mf > 0 else None,
                    "normalized_depth_per_ul": signal["depth"] * cna / mass,
                }
            )
        points.append(point)
    return {
        "sample": sample,
        "path": table_path,
        "fields": fields,
        "rows": rows,
        "amounts": amounts,
        "points": points,
        "signals": signals,
        "total_mass": total_mass,
        "total_bases": total_bases,
        "total_bases_source": total_source,
        "input_sha256": sha256(table_path),
        "fastp_json_sha256": (
            sha256(fastp_path)
            if fastp_path.is_file() and total_source != "sample_manifest"
            else None
        ),
        "alignment_metrics_sha256": sha256(evidence_path) if evidence_path else None,
        "evidence_basis": (
            "alignment_metrics"
            if metrics
            else "esviritu_full_reference_mean_depth_times_length"
        ),
    }


def replicate_qc(samples, settings=None):
    settings = settings or Settings()
    groups = {}
    for data in samples:
        group = data["sample"].get("replicate_group", "NA")
        if group.upper() in MISSING:
            continue
        for point in data["points"]:
            if (
                point["fit_used"]
                and point["sequencing_yield"] is not None
                and point["sequencing_yield"] > 0
            ):
                groups.setdefault(
                    (data["sample"].get("batch_id", "NA"), group, point["ercc_id"]), []
                ).append(point)
    results = []
    for (batch, group, eid), points in groups.items():
        ys = [p["sequencing_yield"] for p in points]
        cv = stats.stdev(ys) / stats.fmean(ys) if len(ys) >= 2 else None
        assessed = len(ys) >= 3
        for p in points:
            p["replicate_yield_cv"] = cv
            if (
                assessed
                and settings.filter_replicate_cv
                and cv > settings.max_replicate_cv
            ):
                p["fit_used"] = False
                p["fit_exclusion"] = "replicate_yield_cv_above_threshold"
        results.append(
            {
                "batch_id": batch,
                "replicate_group": group,
                "ercc_id": eid,
                "n": len(ys),
                "sequencing_yield_cv": cv,
                "three_replicates_assessed": assessed,
                "filter_applied": assessed
                and settings.filter_replicate_cv
                and cv > settings.max_replicate_cv,
                "passes": assessed and cv <= settings.max_replicate_cv,
            }
        )
    return results


def calibrate(data, method, settings):
    eligible = [p for p in data["points"] if p["fit_used"]]
    errors = []
    if method == "sequencing_yield":
        eligible = [
            p
            for p in eligible
            if p["base_fraction"] is not None and p["base_fraction"] > 0
        ]
        x = [math.log10(p["mass_fraction"]) for p in eligible]
        y = [math.log10(p["base_fraction"]) for p in eligible]
        if data["total_bases"] is None:
            errors.append("missing_total_effective_bases")
        axes = (r"$\log_{10}(m_i / M_{total})$", r"$\log_{10}(n_i / n_{total})$")
    else:
        eligible = [
            p
            for p in eligible
            if p["normalized_depth_per_ul"] is not None
            and p["normalized_depth_per_ul"] > 0
        ]
        x = [math.log10(p["normalized_depth_per_ul"]) for p in eligible]
        y = [math.log10(p["extract_equivalent_copies_ul"]) for p in eligible]
        axes = (
            r"$\log_{10}(O_i / 1\;\mathrm{\mu L}^{-1})$",
            "log10(extract-equivalent concentration / 1 copy per uL)",
        )
    fit = regress(x, y)
    run = adjacent_run(p["level"] for p in eligible)
    if run < settings.min_levels:
        errors.append(
            "fewer_than_4_adjacent_levels"
            if settings.min_levels == 4
            else "insufficient_adjacent_levels"
        )
    if fit is None:
        errors.append("insufficient_regression_points")
    else:
        if fit["r2"] is None or fit["r2"] < settings.min_r2:
            errors.append("r2_below_threshold")
        if fit["slope"] <= 0:
            errors.append("nonpositive_slope")
        if (
            method == "sequencing_yield"
            and not settings.slope_min <= fit["slope"] <= settings.slope_max
        ):
            errors.append("yield_slope_outside_acceptance_range")
    yields = (
        [p["sequencing_yield"] for p in eligible]
        if method == "sequencing_yield"
        else []
    )
    result = {
        "method": method,
        "fit": fit,
        "passes": not errors,
        "errors": errors,
        "n_fitted": len(eligible),
        "distinct_levels": len(set(p["level"] for p in eligible)),
        "longest_adjacent_levels": run,
        "x_label": axes[0],
        "y_label": axes[1],
        "mean_sequencing_yield": stats.fmean(yields) if yields else None,
        "median_sequencing_yield": stats.median(yields) if yields else None,
        "fitted_ercc_ids": [p["ercc_id"] for p in eligible],
    }
    return result, eligible, x, y


def wastewater_factor(sample, unit):
    if sample.get("sample_type") != "wastewater":
        return None, {"status": "not_applicable", "unit": unit}
    names = (
        "extract_elution_volume_ul",
        "extraction_input_volume_ul",
        "concentrate_volume_ml",
        "wastewater_volume_ml",
    )
    v = {k: number(sample.get(k), name=k, positive=True) for k in names}
    recovery = number(
        sample.get("pmmov_recovery_fraction"),
        name="pmmov_recovery_fraction",
        positive=True,
    )
    audit = {"unit": unit, **v}
    if recovery is None:
        q0 = number(sample.get("pmmov_raw_extract_copies_ul"), nonnegative=True)
        q1 = number(sample.get("pmmov_post_extract_copies_ul"), nonnegative=True)
        d0 = number(sample.get("pmmov_raw_dilution"), positive=True)
        d1 = number(sample.get("pmmov_post_dilution"), positive=True)
        ve0 = number(sample.get("raw_extract_elution_volume_ul"), positive=True)
        vx0 = number(sample.get("raw_extraction_input_volume_ul"), positive=True)
        loq_pass = flag(sample.get("pmmov_above_loq"))
        if (
            None not in (q0, q1, d0, d1, ve0, vx0, *v.values())
            and q0 > 0
            and q1 > 0
            and loq_pass
        ):
            c0 = 1000 * q0 * d0 * ve0 / vx0
            c1 = (
                1000
                * q1
                * d1
                * v["extract_elution_volume_ul"]
                / v["extraction_input_volume_ul"]
            )
            recovery = (
                c1 * v["concentrate_volume_ml"] / (c0 * v["wastewater_volume_ml"])
            )
            audit.update(
                {
                    "pmmov_raw_copies_ml": c0,
                    "pmmov_concentrate_copies_ml": c1,
                    "recovery_source": "pmmov_mass_balance",
                }
            )
        else:
            audit["recovery_source"] = "missing_or_below_pmmov_loq"
    else:
        audit["recovery_source"] = "sample_manifest"
    audit["pmmov_recovery_fraction"] = recovery
    missing = [k for k, value in v.items() if value is None]
    if recovery is None:
        missing.append("pmmov_recovery_fraction_or_complete_dpcr_inputs")
    if missing:
        audit["status"] = "missing:" + ",".join(missing)
        return None, audit
    if not 0.01 <= recovery <= 1.50:
        audit["status"] = "recovery_outside_acceptance_range"
        return None, audit
    multiplier = 1000 * v["extract_elution_volume_ul"] / v["extraction_input_volume_ul"]
    multiplier *= v["concentrate_volume_ml"] / (v["wastewater_volume_ml"] * recovery)
    if unit == "copies/L":
        multiplier *= 1000
    audit.update(
        {
            "status": "available_pmmov_proxy_correction",
            "extract_to_wastewater_factor": multiplier,
        }
    )
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
        group_id = target.get("target_group") or default_id
        if str(group_id).upper() in MISSING:
            group_id = default_id
        groups.setdefault(group_id, []).append((row, target))
    result = {}
    for gid, items in groups.items():
        lengths = set()
        masses = set()
        for row, target in items:
            length = number(target.get("reference_length_nt"), positive=True)
            if length is None:
                length = number(row.get("Asm_length"), positive=True) or float(
                    row["Length"]
                )
            lengths.add(length)
            mass = number(target.get("molecular_mass_per_nt_g_mol"), positive=True)
            if mass is not None:
                masses.add(mass)
        if len(lengths) != 1 or len(masses) > 1:
            raise ValueError(
                f"Inconsistent quantification reference metadata for {gid}"
            )
        length = next(iter(lengths))
        if sum(float(r["Length"]) for r, _ in items) > length + 1e-6:
            raise ValueError(
                f"Observed reference lengths exceed locked genome length: {gid}"
            )
        n = sum(data["signals"][r["Accession"]]["aligned_bases"] for r, _ in items)
        covered = sum(
            data["signals"][r["Accession"]]["covered_bases"] for r, _ in items
        )
        item = {
            "target_group": gid,
            "reference_length_nt": length,
            "aligned_bases": n,
            "mean_depth": n / length,
            "breadth": covered / length,
            "molecular_mass_per_nt_g_mol": next(iter(masses), None),
            "detected_accessions": ",".join(r["Accession"] for r, _ in items),
            "detected_reference_length_nt": sum(float(r["Length"]) for r, _ in items),
        }
        for row, _ in items:
            result[row["Accession"]] = item
    return result


def quantify_rows(data, calibration, targets, settings):
    """Return source annotations, one result per target group, and recovery evidence.

    Estimated values describe the calibration model. Reportable concentrations
    additionally obey the selected reporting policy and supplied validation.
    """
    sample, method = data["sample"], calibration["method"]
    vlib = float(sample["library_input_volume_ul"])
    reaction = float(sample["reaction_volume_ul"])
    cna = float(sample["extract_nucleic_acid_concentration_ng_ul"])
    mass = float(sample["library_input_mass_ng"])
    factor, recovery_audit = wastewater_factor(sample, settings.wastewater_unit)
    groups = reference_groups(data, targets)
    results, group_results = [], {}

    def resolved(target, key):
        value = target.get(key)
        return (
            sample.get(key)
            if value is None or str(value).strip().upper() in MISSING
            else value
        )

    for row in data["rows"]:
        eid = ercc_id(row)
        result = dict(row)
        flags = []
        ext = estimate = estimated_ww = reported = reported_ext = reported_ww = None
        group = None
        reportable = False
        within_range = None
        loq = None
        validated = False
        if eid:
            known = data["amounts"].get(eid)
            if known and known["added_copies"] > 0:
                reported = known["reaction_copies_ul"]
                flags = ["known_ercc_input"]
            else:
                flags = [
                    "ercc_detected_without_spike" if known else "ercc_absent_from_mix1"
                ]
            value_origin = "assigned_standard"
        else:
            value_origin = "model_estimate"
            group = groups[row["Accession"]]
            if group["detected_reference_length_nt"] < group["reference_length_nt"]:
                flags.append("partial_reference_detected")
            n, length = group["aligned_bases"], group["reference_length_nt"]
            if n <= 0:
                flags.append("no_quantifiable_aligned_bases")
            if not calibration["passes"]:
                flags.append("calibration_failed")
            if calibration["passes"] and n > 0:
                if method == "sequencing_yield":
                    mw = group["molecular_mass_per_nt_g_mol"]
                    if mw is None:
                        flags.append("missing_molecular_mass_factor")
                    else:
                        inferred_mass = (
                            data["total_mass"]
                            * (n / data["total_bases"])
                            / calibration["mean_sequencing_yield"]
                        )
                        ext = inferred_mass * 1e-9 * NA / (mw * length * vlib)
                        log_input = math.log10(inferred_mass / data["total_mass"])
                else:
                    log_input = math.log10(group["mean_depth"] * cna / mass)
                    fit = calibration["fit"]
                    ext = 10 ** (fit["intercept"] + fit["slope"] * log_input)
                if ext is not None:
                    estimate = ext * vlib / reaction
                    estimated_ww = ext * factor if factor is not None else None
                    fit = calibration["fit"]
                    within_range = (
                        fit["x_min"] - 1e-9 <= log_input <= fit["x_max"] + 1e-9
                    )
                    if not within_range:
                        flags.append("outside_calibration_range")
                    member_targets = [
                        targets.get(a, {})
                        for a in group["detected_accessions"].split(",")
                    ]
                    loqs = [
                        number(
                            resolved(t, "loq_reaction_copies_ul"),
                            name="loq_reaction_copies_ul",
                            positive=True,
                        )
                        for t in member_targets
                    ]
                    loq = (
                        max(loqs) if loqs and all(v is not None for v in loqs) else None
                    )
                    alignment_verified = all(
                        flag(resolved(t, "alignment_qc_verified"))
                        for t in member_targets
                    )
                    detection_verified = all(
                        flag(resolved(t, "detection_qc_verified"))
                        for t in member_targets
                    )
                    transfer_verified = all(
                        flag(resolved(t, "calibration_transfer_verified"))
                        for t in member_targets
                    )
                    if not alignment_verified:
                        flags.append("alignment_qc_unverified")
                    if not detection_verified:
                        flags.append("detection_qc_unverified")
                    if not transfer_verified:
                        flags.append("calibration_transfer_unverified")
                    if loq is None:
                        flags.append("loq_unavailable")
                    elif estimate < loq:
                        flags.append("below_loq")
                    validated = (
                        alignment_verified
                        and detection_verified
                        and transfer_verified
                        and loq is not None
                    )
                    policy_pass = validated or settings.reporting_mode == "exploratory"
                    range_pass = within_range or settings.allow_extrapolation
                    above_loq = loq is None or estimate >= loq
                    reportable = policy_pass and range_pass and above_loq
                    if reportable:
                        reported, reported_ext = estimate, ext
                        reported_ww = estimated_ww
            if factor is None and sample["sample_type"] == "wastewater":
                flags.append("wastewater_conversion_unavailable")
            if not flags:
                flags.append("accepted")
        values = {
            CONC: reported,
            WWCONC: reported_ww,
            "quantification_method": method,
            "reaction_concentration_unit": "copies/uL",
            "wastewater_concentration_unit": settings.wastewater_unit,
            "estimated_reaction_copies_ul": estimate,
            "estimated_extract_copies_ul": ext,
            "estimated_wastewater_concentration": estimated_ww,
            "concentration_extract_copies_ul": reported_ext,
            "target_group": group["target_group"] if group else eid,
            "reference_length_nt": (
                group["reference_length_nt"] if group else float(row["Length"])
            ),
            "quantification_aligned_bases": (
                group["aligned_bases"]
                if group
                else data["signals"][row["Accession"]]["aligned_bases"]
            ),
            "quantification_status": ";".join(flags),
            "value_origin": value_origin,
            "reporting_mode": settings.reporting_mode,
            "target_validation_verified": validated,
            "reportable_target_concentration": reportable,
            "within_calibration_range": within_range,
            "loq_reaction_copies_ul": loq,
            "molecular_mass_per_nt_g_mol": (
                group["molecular_mass_per_nt_g_mol"] if group else 340
            ),
            "pmmov_recovery_fraction": recovery_audit.get("pmmov_recovery_fraction"),
        }
        result.update(values)
        results.append(result)
        if group:
            group_results[group["target_group"]] = {
                "sample_id": sample["sample_id"],
                **group,
                **values,
            }
    return results, list(group_results.values()), recovery_audit


def plot_curve(path, x, y, fit, xlabel, ylabel, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.7, 5.2), constrained_layout=True)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 10})
    ax.scatter(x, y, s=24, color="#2563A6", alpha=0.75, edgecolors="none")
    if fit:
        ends = [min(x), max(x)]
        ax.plot(
            ends,
            [fit["intercept"] + fit["slope"] * v for v in ends],
            color="#B44B35",
            lw=1.8,
        )
        r2 = "NA" if fit["r2"] is None else f'{fit["r2"]:.4f}'
        ax.text(
            0.03,
            0.97,
            f'y = {fit["intercept"]:.4f} + {fit["slope"]:.4f} x\n$R^2$ = {r2}; n = {len(x)}',
            transform=ax.transAxes,
            va="top",
            fontsize=10,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "Insufficient calibration evidence",
            transform=ax.transAxes,
            ha="center",
        )
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.spines[["top", "right"]].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(str(path) + "." + ext, dpi=300, bbox_inches="tight")
    plt.close(fig)
