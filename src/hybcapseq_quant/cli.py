"""Portable command-line interface for ERCC-calibrated quantification."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import math
import os
import platform
import re
import shutil
import sys
from pathlib import Path

from . import __version__
from .references import (
    CATALOG_NAME,
    MANIFEST_NAME,
    METADATA_NAME,
    MASS_FACTORS,
    MASS_MODEL,
    annotation_summary,
    build_catalog,
    database_catalog,
    install_catalog,
    load_catalog,
)
from .core import (
    CONC,
    WWCONC,
    METHODS,
    STANDARD_PATH,
    Settings,
    calibrate,
    ercc_amounts,
    flag,
    load_sample,
    number,
    plot_curve,
    quantify_rows,
    read_tsv,
    replicate_qc,
    resolve_path,
    sha256,
    write_tsv,
)

DATA = Path(__file__).parent / "data"
SAMPLE_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def validate_sample_id(value):
    if not SAMPLE_PATTERN.fullmatch(value) or value in {".", ".."}:
        raise ValueError(
            "sample_id must start with an alphanumeric character and contain only letters, digits, periods, underscores or hyphens"
        )


def choose_table(sample, directory, explicit, base):
    if explicit:
        path = Path(explicit)
    elif directory is not None:
        directory = Path(directory)
        if not directory.is_dir():
            raise ValueError(f"EsViritu results directory not found: {directory}")
        name = sample["sample_id"] + ".detected_virus.info.tsv"
        matches = sorted(directory.rglob(name))
        if len(matches) != 1:
            raise ValueError(
                f"Expected one {name}, found {len(matches)} under {directory}"
            )
        return matches[0].resolve()
    else:
        path = resolve_path(sample.get("info_tsv"), base)
        if path is None:
            raise ValueError("Supply --esviritu-dir, --input-tsv or manifest info_tsv")
    if not path.is_file():
        raise ValueError(f"Detection table not found: {path}")
    return path.resolve()


def output_directory(args):
    if args.outdir is not None:
        return args.outdir.resolve()
    if args.processed_dir is not None:
        return (args.processed_dir / "hybcapseq_quant_results").resolve()
    return Path("hybcapseq_quant_results").resolve()


def load_targets(path, datasets, catalog=None):
    """Use accession-resolved genome chemistry, with explicit reviewed overrides."""
    targets = {}
    catalog = catalog or {}
    for data in datasets:
        for row in data["rows"]:
            from .core import ercc_id

            if ercc_id(row):
                continue
            acc = row["Accession"]
            annotation = catalog.get(acc, {})
            if catalog and not annotation:
                raise ValueError(
                    f"Accession absent from molecular annotation catalog: {acc}"
                )
            if annotation:
                for source_key, annotation_key in (
                    ("Length", "length_nt"),
                    ("Asm_length", "assembly_length_nt"),
                ):
                    if float(row[source_key]) != float(annotation[annotation_key]):
                        raise ValueError(
                            f"Detection table/catalog {source_key} mismatch: {acc}"
                        )
                if row.get("Assembly") != annotation["assembly"]:
                    raise ValueError(
                        f"Detection table/catalog assembly mismatch: {acc}"
                    )
            complete = (
                annotation.get("assembly_mass_status") == "resolved_complete_assembly"
            )
            target = {
                "accession": acc,
                "target_group": row.get("Assembly", "NA"),
                "reference_length_nt": row.get("Asm_length") or row["Length"],
                "molecular_mass_per_nt_g_mol": (
                    annotation.get("molecular_mass_per_nt_g_mol", "NA")
                    if complete
                    else "NA"
                ),
                "molecule_type": annotation.get("molecule_type", "NA"),
                "molecular_mass_source": annotation.get(
                    "classification_evidence", "NA"
                ),
                "_annotation_status": (
                    annotation.get("annotation_status")
                    if complete
                    else ";".join(
                        str(annotation.get(k, "molecular_annotation_unavailable"))
                        for k in ("annotation_status", "assembly_mass_status")
                    )
                ),
                "_mass_model": annotation.get("mass_model", "NA"),
                "_reference_length_valid": annotation.get("assembly_mass_status")
                != "inconsistent_assembly_reference_length",
            }
            if acc in targets and targets[acc] != target:
                raise ValueError(f"Conflicting target metadata across samples: {acc}")
            targets[acc] = target
    if path:
        fields, custom = read_tsv(path)
        target_fields, _ = read_tsv(DATA / "targets_template.tsv")
        if "accession" not in fields or set(fields) - set(target_fields):
            raise ValueError(
                "Target manifest requires accession and documented target fields only"
            )
        seen = set()
        for row in custom:
            acc = row["accession"]
            if acc in seen:
                raise ValueError("Duplicate accession in target manifest")
            seen.add(acc)
            target = targets.setdefault(acc, {"accession": acc})
            for key, value in row.items():
                if str(value).strip().upper() not in {"", "NA", "N/A", "NONE", "NULL"}:
                    target[key] = value
            explicit_factor = number(
                row.get("molecular_mass_per_nt_g_mol"), positive=True
            )
            explicit_type = row.get("molecule_type")
            if explicit_factor is not None or explicit_type in MASS_FACTORS:
                factor = explicit_factor or MASS_FACTORS[explicit_type]
                effective_type = target.get("molecule_type")
                if (
                    effective_type in MASS_FACTORS
                    and factor != MASS_FACTORS[effective_type]
                ):
                    raise ValueError(f"Reviewed molecular type/factor mismatch: {acc}")
                source = row.get("molecular_mass_source", "NA")
                if str(source).strip().upper() in {"", "NA"}:
                    raise ValueError(
                        f"Reviewed molecular mass requires a source: {acc}"
                    )
                target["molecular_mass_per_nt_g_mol"] = factor
                target["_annotation_status"] = "reviewed_target_manifest_assignment"
                target["_mass_model"] = (
                    MASS_MODEL
                    if effective_type in MASS_FACTORS
                    else "reviewed_custom_factor"
                )
    return targets


def calibration_points(data, calibration, x, y):
    included = dict(zip(calibration["fitted_ercc_ids"], zip(x, y)))
    rows = []
    for point in data["points"]:
        row = dict(point)
        pair = included.get(point["ercc_id"])
        row["included_in_model"] = pair is not None
        row["predictor_log10"] = pair[0] if pair else None
        row["response_log10"] = pair[1] if pair else None
        fit = calibration["fit"]
        prediction = fit["intercept"] + fit["slope"] * pair[0] if fit and pair else None
        row["fitted_response_log10"] = prediction
        row["residual_log10"] = pair[1] - prediction if prediction is not None else None
        rows.append(row)
    return rows


def reaction_display_model(data, calibration, eligible):
    """Express the primary model on reaction concentration versus mean depth."""
    if not calibration["fit"] or not eligible:
        return [], [], None
    sample = data["sample"]
    x = [math.log10(p["mean_depth"]) for p in eligible]
    y = [math.log10(p["reaction_copies_ul"]) for p in eligible]
    if calibration["method"] == "depth_calibration":
        slope = calibration["fit"]["slope"]
        intercept = calibration["fit"]["intercept"]
        intercept += slope * math.log10(
            float(sample["extract_nucleic_acid_concentration_ng_ul"])
            / float(sample["library_input_mass_ng"])
        )
        intercept += math.log10(
            float(sample["library_input_volume_ul"])
            / float(sample["reaction_volume_ul"])
        )
    elif calibration["mean_sequencing_yield"]:
        slope = 1.0
        intercept = math.log10(
            data["total_mass"]
            * 1e-9
            * 6.02214076e23
            / (
                data["total_bases"]
                * calibration["mean_sequencing_yield"]
                * 340
                * float(sample["reaction_volume_ul"])
            )
        )
    else:
        return x, y, None
    variation = sum((v - sum(y) / len(y)) ** 2 for v in y)
    r2 = (
        1 - sum((v - intercept - slope * u) ** 2 for u, v in zip(x, y)) / variation
        if variation
        else None
    )
    return (
        x,
        y,
        {
            "intercept": intercept,
            "slope": slope,
            "r2": r2,
            "equation": f"log10(C_reaction / 1 copy per uL) = {intercept:.12g} + {slope:.12g} * log10(mean_depth)",
            "derivation": "algebraic transformation of the primary model; no additional fit",
            "molecule_type": "RNA",
            "molecular_mass_per_nt_g_mol": 340,
        },
    )


def settings_from_args(args):
    settings = Settings(
        min_ercc_reads=args.min_ercc_reads,
        min_ercc_breadth=args.min_ercc_breadth,
        min_r2=args.min_r2,
        min_levels=args.min_levels,
        slope_min=args.slope_min,
        slope_max=args.slope_max,
        allow_extrapolation=args.allow_extrapolation,
        reporting_mode=args.reporting_mode,
        filter_replicate_cv=args.filter_replicate_cv,
        max_replicate_cv=args.max_replicate_cv,
        wastewater_unit=args.wastewater_unit,
    )
    if settings.min_levels < 4 or settings.min_ercc_reads < 1:
        raise ValueError(
            "Require at least four adjacent concentration levels and one ERCC read"
        )
    if not 0 <= settings.min_ercc_breadth <= 1 or not 0 <= settings.min_r2 <= 1:
        raise ValueError("Breadth and minimum R-squared must be within [0, 1]")
    if (
        not 0 < settings.slope_min <= settings.slope_max
        or not settings.max_replicate_cv > 0
    ):
        raise ValueError(
            "Slope bounds and maximum replicate CV must be positive and ordered"
        )
    return settings


def run(args):
    base = args.samples.resolve().parent
    fields, samples = read_tsv(args.samples)
    schema_fields, _ = read_tsv(DATA / "samples_template.tsv")
    if "sample_id" not in fields or not samples:
        raise ValueError("Sample manifest requires sample_id and at least one sample")
    if set(fields) - set(schema_fields):
        raise ValueError(
            f"Unknown sample fields: {sorted(set(fields) - set(schema_fields))}"
        )
    if len({s["sample_id"] for s in samples}) != len(samples):
        raise ValueError("Duplicate sample_id in sample manifest")
    for sample in samples:
        validate_sample_id(sample["sample_id"])
    if args.sample:
        missing = set(args.sample) - {s["sample_id"] for s in samples}
        if missing:
            raise ValueError(f"Unknown sample IDs: {sorted(missing)}")
        samples = [s for s in samples if s["sample_id"] in set(args.sample)]
    if args.input_tsv and len(samples) != 1:
        raise ValueError("--input-tsv requires exactly one selected sample")
    standard_fields, standards = read_tsv(args.ercc_standards)
    required_standard_fields = {
        "ercc_id",
        "stock_attomol_ul",
        "length_nt",
        "level",
        "fit_include",
    }
    if not required_standard_fields.issubset(standard_fields):
        raise ValueError(
            "ERCC standards require ercc_id, stock_attomol_ul, length_nt, level and fit_include"
        )
    for standard in standards:
        for key in ("stock_attomol_ul", "length_nt", "level"):
            if number(standard[key], name=key, positive=True) is None:
                raise ValueError(f"Missing required standard value: {key}")
        flag(standard["fit_include"])
    _, expected_standards = read_tsv(STANDARD_PATH)
    expected = {r["ercc_id"] for r in expected_standards}
    if len(standards) != 92 or {r["ercc_id"] for r in standards} != expected:
        raise ValueError(
            "Commercial Mix 1 standards require the 92 documented transcript identifiers"
        )
    # Nominal level ordering is determined by concentration, not editable ordinal labels.
    levels = {
        v: i + 1
        for i, v in enumerate(sorted({float(r["stock_attomol_ul"]) for r in standards}))
    }
    for standard in standards:
        if int(standard["level"]) != levels[float(standard["stock_attomol_ul"])]:
            raise ValueError("ERCC concentration level ordering is inconsistent")
    settings = settings_from_args(args)
    overrides = {}
    if args.ercc_qc:
        _, exclusions = read_tsv(args.ercc_qc)
        for row in exclusions:
            key = (row["sample_id"], row["ercc_id"])
            if key in overrides or row["ercc_id"] not in expected:
                raise ValueError("Duplicate or unknown ERCC QC record")
            overrides[key] = row
    datasets = [
        load_sample(
            s,
            choose_table(s, args.esviritu_dir, args.input_tsv, base),
            standards,
            base,
            settings,
            overrides,
        )
        for s in samples
    ]
    repeat_audit = replicate_qc(datasets, settings)
    catalog, catalog_context = (
        database_catalog(args.db_dir) if args.db_dir else ({}, {})
    )
    targets = load_targets(args.targets, datasets, catalog)
    if (
        args.method in {"both", "sequencing_yield"}
        and not catalog
        and args.targets is None
    ):
        raise ValueError(
            "Sequencing-yield quantification requires --db-dir with molecular annotations "
            "or an explicitly reviewed --targets manifest. Run install-references first."
        )
    root = output_directory(args)
    if root.exists() and any(root.iterdir()):
        raise ValueError(
            f"Output directory must be empty to prevent mixing analyses: {root}"
        )
    # Validate genome grouping before producing any result files.
    from .core import reference_groups

    for data in datasets:
        reference_groups(data, targets)
    root.mkdir(parents=True, exist_ok=True)
    write_tsv(root / "inputs" / "samples.tsv", samples)
    write_tsv(
        root / "inputs" / "targets.tsv",
        list(targets.values()),
        read_tsv(DATA / "targets_template.tsv")[0],
    )
    if catalog:
        used = {r["Accession"] for d in datasets for r in d["rows"]}
        write_tsv(
            root / "inputs" / "molecular_annotations.tsv",
            [catalog[a] for a in sorted(used)],
        )
    shutil.copyfile(args.ercc_standards, root / "inputs" / "ercc_standards.tsv")
    write_tsv(
        root / "replicate_diagnostics.tsv",
        repeat_audit,
        [
            "batch_id",
            "replicate_group",
            "ercc_id",
            "n",
            "sequencing_yield_cv",
            "three_replicates_assessed",
            "filter_applied",
            "passes",
        ],
    )
    summary = []
    for data in datasets:
        sid = data["sample"]["sample_id"]
        write_tsv(
            root / "inputs" / f"{sid}.ercc_inputs.tsv", list(data["amounts"].values())
        )
        for method in METHODS if args.method == "both" else (args.method,):
            directory = root / method / sid
            directory.mkdir(parents=True)
            calibration, eligible, x, y = calibrate(data, method, settings)
            results, assemblies, recovery = quantify_rows(
                data, calibration, targets, settings
            )
            primary = list(data["fields"])
            if any(col in primary for col in (CONC, WWCONC, "quantification_method")):
                raise ValueError("Provide the original EsViritu table")
            at = primary.index("RPKMF") + 1
            primary[at:at] = [CONC, WWCONC]
            extras = list(
                dict.fromkeys(k for r in results for k in r if k not in primary)
            )
            write_tsv(
                directory / f"{sid}.reference_quantification.tsv",
                results,
                primary + extras,
            )
            write_tsv(
                directory / f"{sid}.target_quantification.tsv",
                assemblies,
                (
                    list(assemblies[0])
                    if assemblies
                    else ["sample_id", "target_group", CONC, WWCONC]
                ),
            )
            write_tsv(
                directory / f"{sid}.calibration_points.tsv",
                calibration_points(data, calibration, x, y),
            )
            write_tsv(
                directory / f"{sid}.wastewater_conversion.tsv",
                [{"sample_id": sid, **recovery}],
            )
            fit = calibration["fit"] or {}
            write_tsv(
                directory / f"{sid}.calibration_model.tsv",
                [
                    {
                        "sample_id": sid,
                        "method": method,
                        **fit,
                        "accepted": calibration["passes"],
                        "failure_codes": ";".join(calibration["errors"]),
                        "n_fitted": calibration["n_fitted"],
                        "distinct_levels": calibration["distinct_levels"],
                        "adjacent_levels": calibration["longest_adjacent_levels"],
                        "mean_sequencing_yield": calibration["mean_sequencing_yield"],
                        "median_sequencing_yield": calibration[
                            "median_sequencing_yield"
                        ],
                        "predictor_definition": calibration["x_label"],
                        "response_definition": calibration["y_label"],
                        "fitted_ercc_ids": ",".join(calibration["fitted_ercc_ids"]),
                    }
                ],
            )
            plot_curve(
                directory / f"{sid}.calibration_curve",
                x,
                y,
                calibration["fit"],
                calibration["x_label"],
                calibration["y_label"],
                sid + " | " + method,
            )
            rx, ry, display_model = reaction_display_model(data, calibration, eligible)
            plot_curve(
                directory / f"{sid}.reaction_concentration_curve",
                rx,
                ry,
                display_model,
                "log10(mean depth)",
                "log10(starting reaction concentration / 1 copy per uL)",
                sid + " | reaction concentration",
            )
            write_tsv(
                directory / f"{sid}.reaction_concentration_model.tsv",
                [display_model] if display_model else [],
                [
                    "intercept",
                    "slope",
                    "r2",
                    "equation",
                    "derivation",
                    "molecule_type",
                    "molecular_mass_per_nt_g_mol",
                ],
            )
            metadata = {
                "sample": data["sample"],
                "method": method,
                "software_version": __version__,
                "molecular_annotation_catalog": catalog_context,
                "python_version": platform.python_version(),
                "source_info_tsv": str(data["path"]),
                "source_info_sha256": data["input_sha256"],
                "alignment_basis": data["evidence_basis"],
                "alignment_metrics_sha256": data["alignment_metrics_sha256"],
                "fastp_json_sha256": data["fastp_json_sha256"],
                "total_effective_bases": data["total_bases"],
                "total_bases_source": data["total_bases_source"],
                "total_input_mass_ng": data["total_mass"],
                "calibration": calibration,
                "wastewater_conversion": recovery,
                "settings": vars(settings),
            }
            (directory / f"{sid}.analysis_manifest.json").write_text(
                json.dumps(metadata, indent=2, allow_nan=False) + "\n"
            )
            summary.append(
                {
                    "sample_id": sid,
                    "method": method,
                    "calibration_accepted": calibration["passes"],
                    "calibration_failure_codes": ";".join(calibration["errors"]),
                    "slope": fit.get("slope"),
                    "intercept": fit.get("intercept"),
                    "r2": fit.get("r2"),
                    "n_fitted": len(eligible),
                    "n_reportable_targets": sum(
                        r["reportable_target_concentration"] for r in assemblies
                    ),
                    "n_estimated_targets": sum(
                        r["estimated_reaction_copies_ul"] is not None
                        for r in assemblies
                    ),
                    "n_outside_calibration_range": sum(
                        r["within_calibration_range"] is False for r in assemblies
                    ),
                    "n_missing_molecular_mass": sum(
                        r["molecular_mass_per_nt_g_mol"] is None for r in assemblies
                    ),
                    "n_validated_targets": sum(
                        r["target_validation_verified"] for r in assemblies
                    ),
                    "mean_sequencing_yield": calibration["mean_sequencing_yield"],
                    "wastewater_status": recovery["status"],
                    "output_directory": str(directory),
                }
            )
            print(
                f"{sid}\t{method}\tcalibration_accepted={calibration['passes']}\t{directory}"
            )
    write_tsv(root / "run_summary.tsv", summary)
    from importlib.metadata import version

    input_files = {
        "samples": args.samples,
        "ercc_standards": args.ercc_standards,
        "targets": args.targets,
        "ercc_qc": args.ercc_qc,
    }
    provenance = {
        "schema_version": "1.0",
        "software_version": __version__,
        "started_at_utc": args.started_at_utc,
        "completed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command_line": args.command_line,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": {"matplotlib": version("matplotlib")},
        "input_files": {
            k: {"path": str(v.resolve()), "sha256": sha256(v)}
            for k, v in input_files.items()
            if v
        },
        "molecular_annotation_catalog": catalog_context,
        "reference_provenance": json.loads((DATA / "provenance.json").read_text()),
        "settings": vars(settings),
        "selected_sample_ids": [s["sample_id"] for s in samples],
        "replicate_cv_filter_enabled": settings.filter_replicate_cv,
    }
    (root / "run_manifest.json").write_text(
        json.dumps(provenance, indent=2, allow_nan=False) + "\n"
    )
    return 0


def initialize(args):
    names = [
        "ercc_standards.tsv",
        "samples.tsv",
        "input_fields.tsv",
        "ercc_qc.tsv",
        "targets_template.tsv",
    ]
    for name in names:
        if (args.outdir / name).exists():
            raise ValueError(f"Refusing to overwrite input file: {args.outdir / name}")
    fields, _ = read_tsv(DATA / "samples_template.tsv")
    discovered = {}
    if args.esviritu_dir is not None:
        if not args.esviritu_dir.is_dir():
            raise ValueError(
                f"EsViritu results directory not found: {args.esviritu_dir}"
            )
        for path in sorted(args.esviritu_dir.rglob("*.detected_virus.info.tsv")):
            sid = path.name.removesuffix(".detected_virus.info.tsv")
            if sid in discovered:
                raise ValueError(
                    f"Multiple result files for sample {sid}; select one batch directory"
                )
            header, rows = read_tsv(path)
            if "sample_ID" not in header or any(r["sample_ID"] != sid for r in rows):
                raise ValueError(f"Filename/sample_ID mismatch: {path}")
            discovered[sid] = path
        if not discovered:
            raise ValueError("No original EsViritu detection tables found")
    ids = args.sample if args.sample else list(discovered)
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs requested")
    samples = []
    for sid in ids:
        validate_sample_id(sid)
        if args.esviritu_dir is not None and sid not in discovered:
            raise ValueError(f"Requested sample not found: {sid}")
        row = dict.fromkeys(fields, "NA")
        row.update(
            {
                "sample_id": sid,
                "ercc_mix": "Mix1",
                "alignment_qc_verified": "0",
                "detection_qc_verified": "0",
                "calibration_transfer_verified": "0",
                "pmmov_above_loq": "0",
            }
        )
        if sid in discovered:
            row["info_tsv"] = os.path.relpath(
                discovered[sid].resolve(), args.outdir.resolve()
            )
        samples.append(row)
    args.outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(STANDARD_PATH, args.outdir / "ercc_standards.tsv")
    shutil.copyfile(DATA / "input_fields.tsv", args.outdir / "input_fields.tsv")
    shutil.copyfile(DATA / "targets_template.tsv", args.outdir / "targets_template.tsv")
    write_tsv(args.outdir / "samples.tsv", samples, fields)
    write_tsv(
        args.outdir / "ercc_qc.tsv",
        [],
        ["sample_id", "ercc_id", "use_for_fit", "reason"],
    )
    print(args.outdir)
    return 0


def prepare_targets(args):
    datasets = []
    for path in sorted(args.esviritu_dir.rglob("*.detected_virus.info.tsv")):
        _, rows = read_tsv(path)
        datasets.append({"rows": rows})
    if not datasets:
        raise ValueError("No original EsViritu detection tables found")
    catalog, _ = database_catalog(args.db_dir) if args.db_dir else ({}, {})
    targets = load_targets(None, datasets, catalog)
    if args.db_metadata:
        _, reference = read_tsv(args.db_metadata)
        reference = {r["Accession"]: r for r in reference}
        for acc, target in targets.items():
            if acc not in reference:
                raise ValueError(f"Accession absent from database metadata: {acc}")
            target["reference_length_nt"] = reference[acc]["Asm_length"]
    if args.output.exists():
        raise ValueError(f"Refusing to overwrite target manifest: {args.output}")
    fields, _ = read_tsv(DATA / "targets_template.tsv")
    write_tsv(args.output, list(targets.values()), fields)
    print(args.output)
    return 0


def measure_bases(args):
    total = reads = 0
    for path in [args.r1, args.r2]:
        if not path:
            continue
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as src:
            while True:
                header = src.readline()
                if not header:
                    break
                sequence, plus, quality = (
                    src.readline().rstrip(),
                    src.readline(),
                    src.readline().rstrip(),
                )
                if (
                    not header.startswith("@")
                    or not plus.startswith("+")
                    or not sequence
                    or len(sequence) != len(quality)
                ):
                    raise ValueError(f"Invalid four-line FASTQ record in {path}")
                total += len(sequence)
                reads += 1
    print(
        json.dumps(
            {
                "sample_id": args.sample,
                "total_effective_bases": total,
                "individual_reads": reads,
            }
        )
    )
    return 0


def collect_bam(args):
    try:
        import pysam
    except ImportError as exc:
        raise ValueError(
            "Install optional BAM support: python -m pip install '.[bam]'"
        ) from exc
    if args.min_mapq < 0:
        raise ValueError("Minimum mapping quality must be nonnegative")
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        rows = []
        for name, length in zip(bam.references, bam.lengths):
            events, count, aligned = {}, 0, 0
            for read in bam.fetch(name):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_qcfail
                    or read.is_duplicate
                ):
                    continue
                if read.mapping_quality < args.min_mapq or (
                    read.has_tag("NH") and read.get_tag("NH") != 1
                ):
                    continue
                count += 1
                for start, end in read.get_blocks():
                    aligned += end - start
                    events[start] = events.get(start, 0) + 1
                    events[end] = events.get(end, 0) - 1
            active = covered = previous = 0
            for position, delta in sorted(events.items()):
                if active:
                    covered += position - previous
                active += delta
                previous = position
            rows.append(
                {
                    "accession": name,
                    "length_nt": length,
                    "aligned_bases": aligned,
                    "mean_depth": aligned / length,
                    "covered_bases": covered,
                    "read_count": count,
                }
            )
    write_tsv(args.output, rows)
    evidence = {
        "bam_sha256": sha256(args.bam),
        "minimum_mapping_quality": args.min_mapq,
        "excluded_flags": [
            "unmapped",
            "secondary",
            "supplementary",
            "qc_fail",
            "duplicate",
        ],
        "nh_rule": "NH must equal 1 when present; absent NH is accepted subject to MAPQ",
        "duplicate_rule": "Only pre-existing duplicate flags are used; duplicates are not inferred",
        "pysam_version": pysam.__version__,
    }
    args.output.with_suffix(".json").write_text(json.dumps(evidence, indent=2) + "\n")
    return 0


def manage_references(args):
    if args.command == "install-references":
        records, context = install_catalog(args.db_dir, replace=args.replace)
    elif args.command == "validate-db":
        records, context = database_catalog(args.db_dir)
    else:
        build_catalog(
            args.db_dir / METADATA_NAME,
            args.ictv_vmr,
            args.cache_dir,
            args.output_dir or args.db_dir,
            workers=args.workers,
            ncbi_snapshot=args.ncbi_snapshot,
        )
        output = args.output_dir or args.db_dir
        records, context = load_catalog(
            output / CATALOG_NAME, args.db_dir / METADATA_NAME
        )
    print(json.dumps(annotation_summary(records, context), indent=2))
    return 0


def parser():
    p = argparse.ArgumentParser(
        prog="hybcapseq_quant",
        description="ERCC-calibrated quantification of EsViritu results",
    )
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="quantify selected samples or all manifest samples")
    r.add_argument("--samples", type=Path, required=True)
    r.add_argument("--sample", nargs="+", help="exact sample IDs; omit for batch mode")
    r.add_argument("--esviritu-dir", type=Path)
    r.add_argument("--input-tsv", type=Path)
    r.add_argument(
        "--processed-dir",
        type=Path,
        help="working directory; default output is its hybcapseq_quant_results subdirectory",
    )
    r.add_argument(
        "--outdir",
        type=Path,
        help="explicit output directory; takes priority over --processed-dir",
    )
    r.add_argument("--ercc-standards", type=Path, default=STANDARD_PATH)
    r.add_argument("--targets", type=Path)
    r.add_argument(
        "--db-dir",
        type=Path,
        help="alignment database directory containing accession molecular annotations",
    )
    r.add_argument("--ercc-qc", type=Path)
    r.add_argument("--method", choices=["both", *METHODS], default="both")
    r.add_argument(
        "--reporting-mode", choices=["validated", "exploratory"], default="validated"
    )
    r.add_argument(
        "--wastewater-unit", choices=["copies/mL", "copies/L"], default="copies/mL"
    )
    r.add_argument("--min-ercc-reads", type=int, default=1)
    r.add_argument("--min-ercc-breadth", type=float, default=0)
    r.add_argument("--min-r2", type=float, default=0.95)
    r.add_argument("--min-levels", type=int, default=4)
    r.add_argument("--slope-min", type=float, default=0.8)
    r.add_argument("--slope-max", type=float, default=1.2)
    r.add_argument("--allow-extrapolation", action="store_true")
    r.add_argument(
        "--filter-replicate-cv",
        action="store_true",
        help="exclude ERCCs exceeding the predeclared replicate CV threshold",
    )
    r.add_argument("--max-replicate-cv", type=float, default=0.35)
    r.set_defaults(func=run)
    i = sub.add_parser("init", help="create input templates for the current sample set")
    i.add_argument("--outdir", type=Path, required=True)
    i.add_argument("--esviritu-dir", type=Path)
    i.add_argument("--sample", nargs="+")
    i.set_defaults(func=initialize)
    t = sub.add_parser(
        "prepare-targets", help="prepare target reference and validation metadata"
    )
    t.add_argument("--esviritu-dir", type=Path, required=True)
    t.add_argument("--db-metadata", type=Path)
    t.add_argument("--db-dir", type=Path)
    t.add_argument("--output", type=Path, required=True)
    t.set_defaults(func=prepare_targets)
    b = sub.add_parser(
        "measure-input-bases", help="count actual bases in analysis-input FASTQs"
    )
    b.add_argument("--sample", required=True)
    b.add_argument("--r1", type=Path, required=True)
    b.add_argument("--r2", type=Path)
    b.set_defaults(func=measure_bases)
    c = sub.add_parser(
        "collect-bam", help="collect alignment evidence from an indexed BAM"
    )
    c.add_argument("--bam", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True)
    c.add_argument("--min-mapq", type=int, default=20)
    c.set_defaults(func=collect_bam)
    for command in ("install-references", "validate-db", "build-molecular-annotations"):
        ref = sub.add_parser(
            command, help="install, validate or rebuild accession molecular annotations"
        )
        ref.add_argument("--db-dir", type=Path, required=True)
        if command == "install-references":
            ref.add_argument("--replace", action="store_true")
        if command == "build-molecular-annotations":
            ref.add_argument("--ictv-vmr", type=Path, required=True)
            ref.add_argument("--cache-dir", type=Path, required=True)
            ref.add_argument(
                "--ncbi-snapshot",
                type=Path,
                help="rebuild offline from a saved NCBI summary snapshot",
            )
            ref.add_argument("--output-dir", type=Path)
            ref.add_argument("--workers", type=int, choices=[1, 2, 3], default=3)
        ref.set_defaults(func=manage_references)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    command_line = ["hybcapseq_quant", *argv]
    commands = {
        "run",
        "init",
        "prepare-targets",
        "measure-input-bases",
        "collect-bam",
        "install-references",
        "validate-db",
        "build-molecular-annotations",
    }
    if (
        argv
        and argv[0] not in commands
        and argv[0] not in {"-h", "--help", "--version"}
    ):
        argv.insert(0, "run")
    args = parser().parse_args(argv)
    args.started_at_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    args.command_line = command_line
    try:
        return args.func(args)
    except (ValueError, OSError, KeyError, OverflowError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
