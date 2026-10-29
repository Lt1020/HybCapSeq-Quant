# HybCapSeq-Quant

HybCapSeq-Quant estimates nucleic-acid abundance from ERCC-calibrated hybrid-capture sequencing analyzed with EsViritu. The sequencing-yield and normalized-depth estimators separately produce extract and starting-library-mixture concentrations. Wastewater conversion additionally requires actual experimental volumes and a measured recovery proxy.

Software **2.2.0**; database **[v3.4](databases/v3.4/README.md)** adds 21 reviewed literature PCR reference regions to the corrected v3.3 baseline, with 31,024 total records and a rebuilt index. The four added SARS-CoV-2 regions are N1, N2, E and ORF1ab. Database and software versions are independent.

[中文说明](README_zh.md) · [Methods](docs/METHODS.md) · [Fields](docs/DATA_DICTIONARY.md) · [Verification](VALIDATION.md)

## Install

Use Linux and Python 3.10–3.12. Record the exact commit, dependencies and reference checksums for publication.

```bash
git clone https://github.com/Lt1020/HybCapSeq-Quant.git
cd HybCapSeq-Quant
conda env create -f environment.yml
conda activate hybcapseq_quant
python -m pip install .
hybcapseq_quant --version
```

BAM-derived signal extraction is optional: `python -m pip install '.[bam]'`.

## Database and inputs

Download and verify the [complete v3.4 database](https://github.com/Lt1020/HybCapSeq-Quant/releases/tag/db-v3.4) using the commands in its [database guide](databases/v3.4/README.md). Set `DB_DIR` to the extracted `v3.4` directory:

```bash
DB_DIR="$PWD/v3.4"
hybcapseq_quant validate-db --db-dir "$DB_DIR"
hybcapseq_quant init --esviritu-dir esviritu_results --outdir inputs
hybcapseq_quant prepare-targets \
  --esviritu-dir esviritu_results --db-dir "$DB_DIR" \
  --output inputs/targets.tsv
```

v3.4 includes matching molecular annotations. Re-run EsViritu against it; new sequences require new alignments. The legacy `upgrade-db` command applies only to existing v3.3 installations and is unnecessary for v3.4. Existing v3.3 users can retain that database or follow the [v3.3.1 annotation update](docs/UPGRADE_2.1.0.md).

PCRROI records are published diagnostic windows, not confirmed commercial dPCR kit regions. Genomic and overlapping ROI counts are not additive. Local molecular-mass assignments remain unavailable by default; the release does not assert validated local absolute concentrations. See [reference scope and counting](databases/v3.4/README.md#counting-and-quantitative-interpretation) and the [Chinese workflow](docs/WORKFLOW_v3.4_zh.md).

Review every sample input: actual extract volume, combined DNA+RNA input mass, starting mixture volume, ERCC stock dilution and added volume. Unknown quantities remain `NA`; templates do not assume experiment volumes or recovery. The catalog has 96 ERCC references; only the 92 manufacturer-assigned Mix 1 transcripts receive concentrations. Yield estimation requires analysis-input bases or the correct matching fastp section. See the [field dictionary](docs/DATA_DICTIONARY.md).

## Quantify

```bash
hybcapseq_quant run \
  --samples inputs/samples.tsv --esviritu-dir esviritu_results \
  --db-dir "$DB_DIR" --targets inputs/targets.tsv \
  --outdir quantification_results --method both \
  --min-r2 0.80 --wastewater-unit copies/L
```

Both models require positive slopes, R² ≥0.80 and at least four adjacent assigned levels. The yield slope range 0.8–1.2 is diagnostic only. Explicitly use `--min-r2 0.95 --enforce-yield-slope-range` to reproduce the former acceptance policy. Actual settings and diagnostic codes are exported.

Default `validated` reporting requires external validation declarations and LoQ. Mathematical estimates and reportable concentrations are distinct. `--reporting-mode exploratory` permits exploratory estimates; it does not establish accuracy or analytical limits. `--no-extrapolation` restricts estimates to calibration range.

Method-specific outputs include reference/target concentration tables, calibration points/models and plots, wastewater conversion, and manifests. Source fields are preserved. Use the target-level table; do not sum repeated segment rows. Nonempty output directories are rejected. RPKMF is preserved as an upstream metric; an RPKMF estimator is not implemented in the main `run` command.

## Assess target-specific bias

```bash
hybcapseq_quant fit-correction \
  --pairs reviewed_dpcr_pairs.tsv --outdir correction_assessment
hybcapseq_quant apply-correction \
  --pairs reviewed_dpcr_pairs.tsv \
  --parameters correction_assessment/correction_parameters.tsv \
  --fit-scope pooled --output corrected_pairs.tsv
```

The reviewed pair schema and equations are in [TARGET_CORRECTION.md](docs/TARGET_CORRECTION.md). The estimator balances parents after collapsing repeated dilutions within extraction and method. It evaluates method-specific and shared depth/yield factors, exact feasibility of all observed APE values ≤30%, grouped held-out predictions, and cross-batch transfer. Each validation fit requires at least three training parents. Ineligible, missing and zero comparisons remain auditable.

Application appends estimates and preserves source values. By default a coefficient needs complete leave-parent-out coverage with every held-out APE ≤30%; other coefficients require `--allow-exploratory`. Invalid source calibrations and values outside the factor training range remain unfilled. Passing these data checks is not formal assay validation or a guarantee for future specimens.

## Verify and reproduce

```bash
python -m unittest discover -s tests -v
python examples/quant/generate_demo.py
hybcapseq_quant run --samples examples/quant/samples.tsv \
  --esviritu-dir examples/quant --targets examples/quant/targets.tsv \
  --outdir analytical_example_results
python examples/correction/generate_demo.py
hybcapseq_quant fit-correction --pairs examples/correction/simulated_pairs.tsv \
  --outdir simulated_correction_results
```

Examples contain explicitly simulated measurements with analytical truth. [Software verification](VALIDATION.md) checks computations and data contracts; experimental accuracy, recovery, inhibition, calibration transfer, LoD and LoQ require independent evidence.

## Documentation and attribution

[Manuscript methods](docs/MANUSCRIPT_METHODS.md) describe implemented computations. [CHANGELOG.md](CHANGELOG.md) records migration and version changes. Cite the exact software commit and actual database revision; distinguish the released reference assets from experimental validation of a particular measurement model.

The existing MIT software license is retained. Third-party references retain their attribution and distribution conditions. This public repository contains code, documentation and simulated examples; individual experimental measurements and specimen identifiers are not distributed.
