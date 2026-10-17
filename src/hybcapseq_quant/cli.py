"""Command-line interface for protocol-based ERCC quantification."""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import platform
import shutil
import statistics
import sys
from pathlib import Path

from . import __version__
from .core import (CONC, WWCONC, METHODS, STANDARD_PATH, Settings, calibrate,
                   ercc_amounts, flag, fmt, load_sample, number, plot_curve,
                   quantify_rows, read_tsv, replicate_qc, resolve_path, sha256,
                   write_tsv)

DATA = Path(__file__).parent / "data"


def choose_table(sample, directory, explicit, base):
    # Explicit current-run locations take priority over a saved manifest path.
    if explicit:
        path = Path(explicit)
    elif directory is not None:
        directory = Path(directory)
        if not directory.is_dir():
            raise ValueError(f"EsViritu results directory not found: {directory}")
        name = sample["sample_ID"] + ".detected_virus.info.tsv"
        matches = list(directory.rglob(name))
        if len(matches) != 1:
            raise ValueError(f"Expected one {name}, found {len(matches)} under {directory}")
        return matches[0]
    else:
        path = resolve_path(sample.get("info_tsv"), base)
        if path is None:
            raise ValueError("Supply --esviritu-dir, --input-tsv or manifest info_tsv")
    if not path.is_file():
        raise ValueError(f"Detection table not found: {path}")
    return path


def output_directory(args):
    if args.outdir is not None:
        return args.outdir
    if args.processed_dir is not None:
        return args.processed_dir / "HybCapSeq_Quant_results"
    return Path("HybCapSeq_Quant_results")


def load_targets(path, datasets):
    _, rules = read_tsv(DATA / "molecule_families.tsv")
    rules = {r["family"]: r for r in rules}
    targets = {}
    for data in datasets:
        for row in data["rows"]:
            acc = row["Accession"]
            rule = rules.get(row.get("family"), {})
            if not rule and row.get("kingdom") == "Orthornavirae":
                rule = {"molecule_type": "RNA", "mass_per_nt_g_mol": "340",
                        "source": "Orthornavirae RNA reference; protocol equation 14 uses 340 per reference nt"}
            targets[acc] = {"Accession": acc, "quant_group": row.get("Assembly", "NA"),
                "quant_reference_length": row.get("Asm_length") or row["Length"],
                "mass_per_nt_g_mol": rule.get("mass_per_nt_g_mol", "NA"),
                "molecule_type": rule.get("molecule_type", "NA"),
                "mass_factor_source": rule.get("source", "not_available")}
    if path:
        _, custom = read_tsv(path)
        seen = set()
        for row in custom:
            if row["Accession"] in seen:
                raise ValueError("Duplicate accession in targets manifest")
            seen.add(row["Accession"])
            targets.setdefault(row["Accession"], {}).update(row)
            if str(row.get("mass_per_nt_g_mol", "NA")).upper() in {"", "NA"}:
                known_type = {"RNA": "340", "dsDNA": "660"}.get(row.get("molecule_type"))
                if known_type:
                    targets[row["Accession"]]["mass_per_nt_g_mol"] = known_type
    return targets


def run(args):
    base = args.samples.resolve().parent
    _, samples = read_tsv(args.samples)
    if not samples or len({s["sample_ID"] for s in samples}) != len(samples):
        raise ValueError("Sample manifest is empty or contains duplicate sample_ID")
    for s in samples:
        sid = s["sample_ID"]
        if Path(sid).name != sid or sid in {".", ".."} or "\\" in sid:
            raise ValueError("sample_ID cannot contain path separators")
    if args.sample:
        wanted = set(args.sample)
        missing = wanted - {s["sample_ID"] for s in samples}
        if missing:
            raise ValueError(f"Unknown sample IDs: {sorted(missing)}")
        samples = [s for s in samples if s["sample_ID"] in wanted]
    if args.input_tsv and len(samples) != 1:
        raise ValueError("--input-tsv requires exactly one selected sample")
    _, standards = read_tsv(args.ercc_standards)
    if len(standards) != 92:
        raise ValueError("Commercial Mix 1 standards must contain exactly 92 transcripts")
    settings = Settings(min_ercc_reads=args.min_ercc_reads, min_ercc_breadth=args.min_ercc_breadth,
        min_r2=args.min_r2, min_levels=args.min_levels, slope_min=args.slope_min, slope_max=args.slope_max,
        allow_extrapolation=args.allow_extrapolation, strict_qc=args.strict_qc,
        ascii_columns=args.ascii_columns, wastewater_unit=args.wastewater_unit)
    if settings.min_levels < 4 or settings.min_ercc_reads < 1 or not 0 <= settings.min_ercc_breadth <= 1:
        raise ValueError("Require >=4 levels, >=1 ERCC reads and breadth in [0,1]")
    overrides = {}
    if args.ercc_qc:
        _, exclusions = read_tsv(args.ercc_qc)
        for row in exclusions:
            key = (row["sample_ID"], row["ERCC_ID"])
            if key in overrides:
                raise ValueError("Duplicate sample/ERCC QC override")
            overrides[key] = row
    datasets = [load_sample(s, choose_table(s, args.esviritu_dir, args.input_tsv, base),
                           standards, base, settings, overrides) for s in samples]
    repeat_audit = replicate_qc(datasets)
    targets = load_targets(args.targets, datasets)
    root = output_directory(args)
    root.mkdir(parents=True, exist_ok=True)
    write_tsv(root / "inputs" / "samples_used.tsv", samples)
    write_tsv(root / "inputs" / "targets_resolved.tsv", list(targets.values()))
    shutil.copyfile(args.ercc_standards, root / "inputs" / "ercc92_mix1.tsv")
    write_tsv(root / "replicate_ERCC_Yseq_CV.tsv", repeat_audit,
        ["batch_ID", "replicate_group", "ERCC_ID", "n", "Y_seq_cv", "three_replicates_assessed", "passes"])
    summary = []
    for data in datasets:
        sid = data["sample"]["sample_ID"]
        if Path(sid).name != sid or sid in {".", ".."} or "\\" in sid:
            raise ValueError("sample_ID cannot contain path separators")
        write_tsv(root / "inputs" / (sid + ".ERCC92_inputs.tsv"), list(data["amounts"].values()))
        for method in (METHODS if args.method == "both" else (args.method,)):
            directory = root / method / sid
            directory.mkdir(parents=True, exist_ok=True)
            calibration, eligible, x, y = calibrate(data, method, settings)
            results, assemblies, recovery = quantify_rows(data, calibration, targets, settings)
            primary = list(data["fields"])
            if CONC in primary or WWCONC in primary or "concentration_ERCC" in primary:
                raise ValueError("Input is already quantified; provide the original EsViritu table")
            at = primary.index("RPKMF") + 1
            primary[at:at] = [CONC, WWCONC]
            extras = list(dict.fromkeys(k for r in results for k in r if k not in primary))
            cols = primary + extras
            if settings.ascii_columns:
                mapping = {CONC: "concentration_ERCC", WWCONC: "concentration_wastewater"}
                results = [{mapping.get(k, k): v for k, v in r.items()} for r in results]
                assemblies = [{mapping.get(k, k): v for k, v in r.items()} for r in assemblies]
                cols = [mapping.get(k, k) for k in cols]
            write_tsv(directory / f"{sid}.detected_virus.info.quantified.{method}.tsv", results, cols)
            write_tsv(directory / f"{sid}.assembly.quantified.{method}.tsv", assemblies,
                      list(assemblies[0]) if assemblies else ["sample_ID", "quant_group"])
            write_tsv(directory / f"{sid}.ERCC_calibration_points.tsv", data["points"])
            write_tsv(directory / f"{sid}.wastewater_conversion.tsv", [{"sample_ID": sid, **recovery}])
            fit = calibration["fit"] or {}
            write_tsv(directory / f"{sid}.ERCC_regression.tsv", [{"sample_ID": sid, "method": method,
                **fit, "passes": calibration["passes"], "errors": ";".join(calibration["errors"]),
                "n_fitted": calibration["n_fitted"], "distinct_levels": calibration["distinct_levels"],
                "adjacent_levels": calibration["longest_adjacent_levels"],
                "Y_seq_mean": calibration["Y_seq_mean"], "Y_seq_median": calibration["Y_seq_median"],
                "x_definition": calibration["x_label"], "y_definition": calibration["y_label"],
                "fit_ERCC_IDs": ",".join(calibration["fit_ERCC_IDs"])}])
            plot_curve(directory / f"{sid}.ERCC_standard_curve", x, y, calibration["fit"],
                       calibration["x_label"], calibration["y_label"], sid + " | " + method)
            # Display the same model on this sample's starting-reaction scale.
            rx = [math.log10(p["mean_depth"]) for p in eligible]
            ry = [math.log10(p["reaction_copies_ul"]) for p in eligible]
            display_fit = None
            if calibration["fit"] and eligible:
                s = data["sample"]
                if method == "method2_depth":
                    b = fit["slope"]
                    a = fit["intercept"] + b * math.log10(float(s["nucleic_acid_conc_ng_ul"]) / float(s["sample_input_mass_ng"]))
                    a += math.log10(float(s["library_extract_ul"]) / float(s["reaction_volume_ul"]))
                else:
                    b = 1.0
                    a = math.log10(data["total_mass"] * 1e-9 * 6.02214076e23 /
                        (data["total_bases"] * calibration["Y_seq_mean"] * 340 * float(s["reaction_volume_ul"])))
                yy = sum((v - sum(ry) / len(ry)) ** 2 for v in ry)
                dr2 = 1 - sum((v - (a + b * u)) ** 2 for u, v in zip(rx, ry)) / yy if yy else None
                display_fit = {"intercept": a, "slope": b, "r2": dr2}
                write_tsv(directory / f"{sid}.reaction_display_equation.tsv", [{**display_fit,
                    "equation": f"log10(C_reaction) = {a:.12g} + {b:.12g} * log10(mean_depth)",
                    "note": "algebraic display of the primary model; not a replacement fit"}])
            plot_curve(directory / f"{sid}.ERCC_reaction_concentration_curve", rx, ry, display_fit,
                       "log10(ERCC mean depth)", "log10(known ERCC copies/uL in starting reaction)", sid + " | reaction scale")
            metadata = {"sample": data["sample"], "method": method, "software_version": __version__,
                "python": platform.python_version(), "source_info_tsv": str(data["path"]),
                "source_info_sha256": data["input_sha256"], "alignment_basis": data["evidence_basis"],
                "alignment_metrics_sha256": data["alignment_metrics_sha256"],
                "total_effective_bases": data["total_bases"], "total_bases_source": data["total_bases_source"],
                "total_input_mass_including_ERCC_ng": data["total_mass"], "calibration": calibration,
                "wastewater_conversion": recovery, "settings": vars(settings)}
            (directory / f"{sid}.quantification_audit.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
            summary.append({"sample_ID": sid, "method": method, "calibration_pass": calibration["passes"],
                "calibration_errors": ";".join(calibration["errors"]), "slope": fit.get("slope"),
                "intercept": fit.get("intercept"), "r2": fit.get("r2"), "n_fitted": len(eligible),
                "Y_seq_mean": calibration["Y_seq_mean"], "Y_seq_median": calibration["Y_seq_median"],
                "wastewater_status": recovery["status"], "output_directory": str(directory)})
            print(f"{sid}\t{method}\tcurve_pass={calibration['passes']}\t{directory}")
    write_tsv(root / "run_summary.tsv", summary)
    library_groups = {}
    for row in summary:
        record = next(s for s in samples if s["sample_ID"] == row["sample_ID"])
        group = record.get("replicate_group", "NA")
        if row["method"] == "method1_yield" and row["Y_seq_mean"] and group != "NA":
            library_groups.setdefault((record.get("batch_ID", "NA"), group), []).append(row["Y_seq_mean"])
    library_cv = []
    for (batch, group), values in library_groups.items():
        cv = statistics.stdev(values) / statistics.fmean(values) if len(values) > 1 else None
        library_cv.append({"batch_ID": batch, "replicate_group": group, "n": len(values),
            "library_Y_seq_mean": statistics.fmean(values), "library_Y_seq_CV": cv,
            "passes": len(values) >= 3 and cv <= 0.35})
    write_tsv(root / "replicate_library_Yseq_CV.tsv", library_cv,
        ["batch_ID", "replicate_group", "n", "library_Y_seq_mean", "library_Y_seq_CV", "passes"])
    provenance = {"version": __version__, "samples_sha256": sha256(args.samples),
        "ercc_standards_sha256": sha256(args.ercc_standards),
        "targets_sha256": sha256(args.targets) if args.targets else None,
        "ercc_qc_sha256": sha256(args.ercc_qc) if args.ercc_qc else None,
        "formula_provenance": json.loads((DATA / "provenance.json").read_text()),
        "note": "Input info TSVs alone do not certify unique alignment, deduplication, detection validity or LoQ."}
    (root / "run_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n")
    return 0


def initialize(args):
    if args.preset == "20261017" and (args.esviritu_dir is not None or args.sample):
        raise ValueError("Historical preset cannot be combined with new-sample discovery")
    sample_name = "samples_20261017.tsv" if args.preset == "20261017" else "samples.tsv"
    names = ["ercc92_mix1.tsv", sample_name, "input_fields.tsv", "ercc_qc.tsv"]
    # Check all destinations before writing, including an existing exclusions file.
    for name in names:
        if (args.outdir / name).exists():
            raise ValueError(f"Refusing to overwrite input file: {args.outdir / name}")
    fields, _ = read_tsv(DATA / "samples_template.tsv")
    samples = []
    discovered = {}
    if args.esviritu_dir is not None:
        if not args.esviritu_dir.is_dir():
            raise ValueError(f"EsViritu results directory not found: {args.esviritu_dir}")
        for path in sorted(args.esviritu_dir.rglob("*.detected_virus.info.tsv")):
            sid = path.name.removesuffix(".detected_virus.info.tsv")
            if sid in discovered:
                raise ValueError(f"Multiple result files for sample {sid}; choose a single batch directory")
            header, rows = read_tsv(path)
            if "sample_ID" not in header or any(r["sample_ID"] != sid for r in rows):
                raise ValueError(f"Filename/sample_ID mismatch: {path}")
            discovered[sid] = path
        if not discovered:
            raise ValueError("No original EsViritu detection tables found")
    ids = args.sample if args.sample else list(discovered)
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs requested")
    for sid in ids:
        if Path(sid).name != sid or sid in {".", ".."} or "\\" in sid:
            raise ValueError("sample_ID cannot contain path separators")
        if args.esviritu_dir is not None and sid not in discovered:
            raise ValueError(f"Requested sample not found in results directory: {sid}")
        row = {key: "NA" for key in fields}
        row.update({"sample_ID": sid, "ercc_mix": "Mix1", "fastp_section": "before_filtering",
                    "alignment_qc_verified": "0", "detection_qc_verified": "0", "pmmov_above_loq": "0",
                    "notes": "New-batch template; fill actual sample type, inputs, ERCC addition and recovery parameters"})
        if sid in discovered:
            row["info_tsv"] = os.path.relpath(discovered[sid].resolve(), args.outdir.resolve())
        samples.append(row)
    args.outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(STANDARD_PATH, args.outdir / "ercc92_mix1.tsv")
    shutil.copyfile(DATA / "input_fields.tsv", args.outdir / "input_fields.tsv")
    if args.preset == "20261017":
        shutil.copyfile(DATA / "samples_20261017.tsv", args.outdir / sample_name)
    else:
        write_tsv(args.outdir / sample_name, samples, fields)
    write_tsv(args.outdir / "ercc_qc.tsv", [], ["sample_ID", "ERCC_ID", "use_for_fit", "reason"])
    print(args.outdir)
    return 0


def prepare_targets(args):
    datasets = []
    for path in args.esviritu_dir.rglob("*.detected_virus.info.tsv"):
        _, rows = read_tsv(path)
        datasets.append({"rows": rows})
    if not datasets:
        raise ValueError("No original EsViritu detection tables found")
    targets = load_targets(None, datasets)
    if args.db_metadata:
        _, reference = read_tsv(args.db_metadata)
        reference = {r["Accession"]: r for r in reference}
        for acc, target in targets.items():
            row = reference.get(acc)
            if row is None:
                raise ValueError(f"Accession absent from database metadata: {acc}")
            target["quant_reference_length"] = row["Asm_length"]
    write_tsv(args.output, list(targets.values()))
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
                sequence, plus, quality = src.readline().rstrip(), src.readline(), src.readline().rstrip()
                if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                    raise ValueError(f"Invalid four-line FASTQ record in {path}")
                total += len(sequence)
                reads += 1
    print(json.dumps({"sample_ID": args.sample, "total_effective_bases": total, "individual_reads": reads}))
    return 0


def collect_bam(args):
    try:
        import pysam
    except ImportError as exc:
        raise ValueError("Install optional BAM support: python -m pip install '.[bam]'") from exc
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        rows = []
        for name, length in zip(bam.references, bam.lengths):
            events, count, aligned = {}, 0, 0
            for read in bam.fetch(name):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate:
                    continue
                if read.mapping_quality < args.min_mapq or (read.has_tag("NH") and read.get_tag("NH") != 1):
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
            rows.append({"Accession": name, "Length": length, "aligned_bases": aligned,
                         "mean_depth": aligned / length, "covered_bases": covered, "read_count": count})
    write_tsv(args.output, rows)
    print("Excluded secondary/supplementary/QC-fail/duplicate-flagged alignments; "
          "this does not identify unmarked PCR duplicates or validate coverage uniformity.", file=sys.stderr)
    return 0


def parser():
    p = argparse.ArgumentParser(prog="hybcapseq_quant", description="ERCC92 two-method calibration of EsViritu TSV results")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="quantify one selected sample or all manifest samples")
    r.add_argument("--samples", type=Path, required=True)
    r.add_argument("--sample", nargs="+", help="exact sample IDs; omit for batch mode")
    r.add_argument("--esviritu-dir", type=Path)
    r.add_argument("--input-tsv", type=Path)
    r.add_argument("--processed-dir", type=Path, help="current batch working directory; default output is its HybCapSeq_Quant_results subdirectory")
    r.add_argument("--outdir", type=Path, help="explicit output directory; overrides --processed-dir output default")
    r.add_argument("--ercc-standards", type=Path, default=STANDARD_PATH)
    r.add_argument("--targets", type=Path)
    r.add_argument("--ercc-qc", type=Path)
    r.add_argument("--method", choices=["both", *METHODS], default="both")
    r.add_argument("--wastewater-unit", choices=["copies/mL", "copies/L"], default="copies/mL")
    r.add_argument("--min-ercc-reads", type=int, default=1)
    r.add_argument("--min-ercc-breadth", type=float, default=0)
    r.add_argument("--min-r2", type=float, default=0.95)
    r.add_argument("--min-levels", type=int, default=4)
    r.add_argument("--slope-min", type=float, default=0.8)
    r.add_argument("--slope-max", type=float, default=1.2)
    r.add_argument("--allow-extrapolation", action="store_true")
    r.add_argument("--strict-qc", action="store_true", help="suppress target values unless supplied alignment/detection QC and LoQ are verified")
    r.add_argument("--ascii-columns", action="store_true")
    r.set_defaults(func=run)
    i = sub.add_parser("init", help="create a generic new-batch template, or explicitly select historical records")
    i.add_argument("--outdir", type=Path, required=True)
    i.add_argument("--esviritu-dir", type=Path, help="discover sample IDs from this actual results directory")
    i.add_argument("--sample", nargs="+", help="sample IDs to initialize; may filter discovered results")
    i.add_argument("--preset", choices=["generic", "20261017"], default="generic",
                   help="generic by default; 20261017 explicitly copies the actual 13 historical records")
    i.set_defaults(func=initialize)
    t = sub.add_parser("prepare-targets", help="generate reviewable genome-length/molecular-mass input")
    t.add_argument("--esviritu-dir", type=Path, required=True)
    t.add_argument("--db-metadata", type=Path)
    t.add_argument("--output", type=Path, required=True)
    t.set_defaults(func=prepare_targets)
    b = sub.add_parser("measure-input-bases", help="count actual bases in the FASTQs used by EsViritu")
    b.add_argument("--sample", required=True)
    b.add_argument("--r1", type=Path, required=True)
    b.add_argument("--r2", type=Path)
    b.set_defaults(func=measure_bases)
    c = sub.add_parser("collect-bam", help="optional explicit alignment-base evidence from an indexed BAM")
    c.add_argument("--bam", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True)
    c.add_argument("--min-mapq", type=int, default=20)
    c.set_defaults(func=collect_bam)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    commands = {"run", "init", "prepare-targets", "measure-input-bases", "collect-bam"}
    if argv and argv[0] not in commands and argv[0] not in {"-h", "--help", "--version"}:
        argv.insert(0, "run")
    p = parser()
    args = p.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError, KeyError, OverflowError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
