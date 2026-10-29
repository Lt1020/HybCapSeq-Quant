# Change log

## Database v3.4 — 2026-10-29

- Append 21 reviewed diagnostic ROI sequences to v3.3, retaining the verified v3.3.1 metadata correction: 31,024 total references, including 96 ERCCs.
- Preserve four SARS-CoV-2 ROI records (N1, N2, E, ORF1ab); retain the full 35-assay audit and 13 excluded candidate windows.
- Rebuild and verify the minimap2 short-read index, add compatible local-reference annotations, and publish the complete database with checksums.
- Document overlapping-reference counting, local versus genomic measurement definitions, installation and manuscript methods. Software remains 2.2.0; calibration defaults remain R²≥0.80 and diagnostic-only yield-slope range.


## 2.2.0

- Set both calibration defaults to R² ≥0.80; retain positive slopes and four adjacent levels.
- Make the yield slope range 0.8–1.2 diagnostic. `--enforce-yield-slope-range` restores the gate; add `--min-r2 0.95` for the former default policy.
- Add grouped target-specific correction, exact ±30% feasibility, held-out/cross-batch assessment, and application preserving original values.
- Add simulated correction examples and reference-input preflight.
- Update English/Chinese documentation and manuscript methods.
- Retain active database v3.3.1. Database v3.4 preparation does not include sequence assets, a new index, or a completed release.

## 2.1.0

See [release notes](docs/RELEASE_NOTES_2.1.0.md) and [migration](docs/UPGRADE_2.1.0.md) for the accession annotation, v3.3.1 correction and exploratory reporting policy.
