# HybCapSeq Quant 2.2.0

This release adds grouped assessment of target-specific multiplicative bias and updates the research calibration defaults. The active reference database remains v3.3.1.

## Calibration policy

Both estimators default to R² ≥0.80 with a positive slope and at least four adjacent assigned levels. Yield slopes outside 0.8–1.2 are reported as diagnostics. Explicit `--min-r2 0.95 --enforce-yield-slope-range` reproduces the previous acceptance policy. Reporting-mode, interpolation, LoQ and external-validation requirements are unchanged.

## Target correction

`fit-correction` estimates method-specific and shared depth/yield factors, checks exact feasibility of all observed APE values within 30%, and exports leave-parent-out and cross-batch predictions. It collapses repeated ERCC conditions within extraction and method before balancing parent specimens. Each independent training split requires at least three parents. Missing or ineligible comparisons remain in the audit.

`apply-correction` appends estimates while preserving original concentrations. Factors without complete held-out performance within 30% require explicit exploratory application. Invalid source calibration and values outside a factor's training interval remain missing. These data checks do not certify analytical accuracy or guarantee future results.

## Database and reproducibility

The v3.3.1 annotation update is retained. The v3.4 directory contains a preparation specification and input preflight only: reviewed assay-region references and the complete alignment database have not been supplied. No v3.4 FASTA, index or database release is provided, and no local-region analysis has been performed.

English/Chinese instructions, computational methods, manuscript methods, correction field definitions and synthetic examples accompany the code. Actual experimental concentrations and specimen identifiers are not published in the repository. See `VALIDATION.md` for numerical checks and their limits. The included v3.3.1 update asset is annotation-only; it is not a v3.4 alignment database.
