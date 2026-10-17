# Quantification methods

## Scope and measurement definition

HybCapSeq-Quant estimates target nucleic-acid copies from ERCC-calibrated hybrid-capture sequencing. Its two estimators operate on the same observed reference signal but use different calibration equations. They are reported separately and are not averaged. Copies refer to reference-equivalent nucleic-acid molecules; sequencing does not establish infectious-particle counts.

ERCC is added to the library reaction after nucleic-acid extraction. Its assigned concentration therefore calibrates library input and sequencing response. Back-conversion to wastewater additionally requires extraction/concentration volumes and an independently determined recovery proxy.

## Notation and units

| Symbol | Definition | Unit |
|---|---|---|
| $N_A$ | Avogadro constant, exactly $6.02214076\times10^{23}$ | mol⁻¹ |
| $s_i$ | Manufacturer stock assignment for ERCC transcript $i$ | amol/µL |
| $D$ | Total stock dilution factor | dimensionless |
| $V_s$ | Volume of diluted ERCC added | µL |
| $V_l$ | Extract aliquot entering the library reaction | µL |
| $V_r$ | Total starting reaction volume | µL |
| $M_l$ | Mixed DNA plus RNA input mass before ERCC addition | ng |
| $C_{NA}$ | Mixed nucleic-acid concentration of the extract | ng/µL |
| $L_i,L_v$ | ERCC transcript length or locked complete target reference length | nt or bp |
| $n_i,n_v$ | Cumulative aligned read bases | bases |
| $n_T$ | Total measured bases in the analysis-input read set | bases |
| $k_v$ | Approximate molecular mass per reference position | g mol⁻¹ nt⁻¹ or g mol⁻¹ bp⁻¹ |

The input identity $M_l=C_{NA}V_l$ is checked with relative tolerance $10^{-5}$ and absolute tolerance $10^{-7}$ ng. For RNA, $k=340$; for double-stranded DNA, $k=660$. These are average molecular-mass approximations. The RNA factor is not assigned to unknown DNA chemistries.

## ERCC assignments

The bundled table contains the manufacturer's 92 Mix 1 transcript assignments and 22 distinct nominal concentration levels. NIST reference identifiers ERCC-00007, ERCC-00018, ERCC-00023 and ERCC-00128 do not receive Mix 1 concentration assignments. Transcript lengths match the certified-v1 sequence representation in the combined reference database.

For each spiked transcript:

$$
N_i=s_i10^{-18}N_A V_s/D,
\qquad
m_i=N_i L_i(340)10^9/N_A.
$$

The known starting-reaction and extract-equivalent concentrations are:

$$
C_{r,i}=N_i/V_r,\qquad A_i=N_i/V_l.
$$

The latter is a mathematical equivalent for calibration, not a claim that ERCC was present in the extracted specimen. Total library mass is

$$
M_T=M_l+\sum_{i=1}^{92}m_i.
$$

All spiked ERCC mass is included, including transcripts that are undetected or excluded from fitting. Nominal stock assignments do not include lot-specific measurement, dilution or pipetting uncertainty.

## Alignment signal

For the verified EsViritu v1.1.2 coverage implementation, `mean_coverage` is the mean over the full reference, including zero-depth positions. Consequently, the summary-derived signal is

$$
n_j=\mathrm{mean\_coverage}_j L_j.
$$

`covered_bases` counts distinct covered reference positions and is used for breadth, not cumulative aligned bases. `read_count` and `RPKMF` are retained as upstream observations; neither replaces $n_j$ in these equations. Summary-derived alignment bases may be approximate if upstream depth values have been rounded.

An optional alignment-metrics table can replace the summary signal. `collect-bam` counts `get_blocks()` aligned bases after excluding unmapped, secondary, supplementary, QC-failed and duplicate-flagged alignments, applying a MAPQ threshold, and requiring NH=1 when NH exists. It does not identify unmarked PCR duplicates, validate UMI processing, or prove unique mapping when NH is absent. Overlapping mates remain separate read-base observations; no fragment-level overlap correction is applied. Total-input bases must use a compatible read-counting convention.

Measured total-input bases can be supplied directly or read from an explicitly selected fastp report section. When fastp is used, its total read count must match the EsViritu sample read count. The sum of reference-aligned bases must not exceed total input bases beyond numerical tolerance. A matching read count alone does not establish that the correct filtering stage was selected.

## Sequencing-yield estimator

For every eligible ERCC transcript, let

$$
f_i=m_i/M_T,\qquad q_i=n_i/n_T,\qquad Y_i=q_i/f_i.
$$

Unweighted ordinary least squares is fitted in the specified direction:

$$
\log_{10}q_i=a+b\log_{10}f_i.
$$

The regression evaluates response proportionality. The concentration estimator uses the arithmetic mean

$$
\bar Y=\frac{1}{K}\sum_{i\in\mathcal I}Y_i,
$$

where $\mathcal I$ is the eligible transcript set and $K$ is its size. The median is also reported as a diagnostic. Regression inversion is not substituted for this estimator.

For target $v$:

$$
\widehat m_v=M_T\frac{n_v/n_T}{\bar Y},
\quad
\widehat N_v=\frac{\widehat m_v10^{-9}N_A}{k_v L_v},
\quad
\widehat C_{e,v}=\widehat N_v/V_l,
\quad
\widehat C_{r,v}=\widehat N_v/V_r.
$$

The calibration range is the fitted ERCC log mass-fraction interval. A target's inferred mass fraction is checked against that interval. If the molecular-mass factor is unavailable, the yield concentration is unavailable.

## Normalized-depth estimator

For each eligible ERCC transcript:

$$
d_i=n_i/L_i,\qquad O_i=d_i C_{NA}/M_l.
$$

$O_i$ has units µL⁻¹. Numerical logarithms use concentrations expressed in copies/µL and $O_i$ expressed in µL⁻¹. The model is

$$
\log_{10}\!\left(\frac{A_i}{1\ \mathrm{copy}/\mu\mathrm L}\right)
=\alpha+\beta\log_{10}\!\left(\frac{O_i}{1\ \mu\mathrm L^{-1}}\right).
$$

For a target reference group:

$$
d_v=n_v/L_v,\quad O_v=d_v C_{NA}/M_l,
\quad
\widehat C_{e,v}=10^{\alpha+\beta\log_{10}O_v},
\quad
\widehat C_{r,v}=\widehat C_{e,v}V_l/V_r.
$$

The fitted ERCC log-normalized-depth interval defines the interpolation range. ERCC RNA calibration is empirical; transfer to targets, including DNA targets, must be supported independently.

## Reference groups and segmented targets

The default quantification group is the upstream `Assembly`, falling back to `Accession` when unavailable. Users can supply reviewed group identifiers and complete reference lengths. Aligned bases are summed across detected segments, while depth and molecular mass use the locked complete reference length. Missing segments do not shorten the denominator. Partial reference detection is flagged.

Every source segment row receives its group's concentration. `target_quantification.tsv` contains each group once and is the appropriate table for target-level summaries. Repeated segment rows must not be added. Alternative assemblies, strains and shared reference regions require explicit review; the software does not automatically merge them into a species abundance or redistribute ambiguous reads.

## Wastewater conversion and recovery

Let $V_{e1}$ be concentrate-extract elution volume (µL), $V_{x1}$ the actual concentrate aliquot extracted (µL), $V_c$ total concentrate volume (mL), and $V_0$ original wastewater volume (mL). With recovery fraction $r$:

$$
\widehat C_{c,v}=1000\widehat C_{e,v}V_{e1}/V_{x1},
\qquad
\widehat C_{w,v}=\widehat C_{c,v}\frac{V_c}{V_0r}
\quad[\mathrm{copies/mL}].
$$

Multiply by 1000 for copies/L. The concentrate aliquot excludes lysis reagents; the elution volume is the actual recovered volume.

The supplied recovery fraction may be used directly. Alternatively, for original-water and concentrate dPCR extract measurements $q_0,q_1$, dilution factors $D_0,D_1$, and original-water extraction volumes $V_{e0},V_{x0}$:

$$
C_0=1000q_0D_0V_{e0}/V_{x0},\qquad
C_1=1000q_1D_1V_{e1}/V_{x1},
\qquad
r=\frac{C_1V_c}{C_0V_0}.
$$

This is mass recovery, not the concentrate-to-original concentration ratio. Both dPCR measurements must satisfy external quantification/inhibition criteria before this route is enabled. The software accepts recovery fractions in the predeclared review interval 0.01–1.50 and flags values outside it. Values above one may reflect measurement variability and require scientific review; acceptance does not prove recovery accuracy. Missing recovery is not replaced by one. PMMoV is a proxy whose relationship to each target's recovery remains an experimental assumption.

## Calibration eligibility and reporting

A transcript must be detected, spiked, have positive aligned signal, satisfy the configured read and breadth thresholds, and not be excluded in the standards or ERCC QC table. All 92 standards, including undetected transcripts, remain in the point audit. ERCC exclusions require a reason. An undetected transcript is not assigned a fabricated observed depth and does not enter a logarithmic fit. Target LoQ values are not used to infer ERCC eligibility.

Default model acceptance criteria are at least four adjacent nominal concentration levels, R²≥0.95, positive slope, and a yield-model slope of 0.8–1.2. These are configurable analysis criteria, not universal analytical-performance guarantees. Eligible transcripts have equal regression weight; concentration levels are not collapsed, and no maximum-R² subset search is performed. The reported RMSE is $\sqrt{\mathrm{SSE}/K}$ in log10 units; it is not a prediction interval. Regression treats predictors as fixed and does not model measurement error, censoring bias or heteroscedasticity.

Replicate diagnostics use the sample standard deviation divided by the arithmetic mean of each transcript's $Y_i$, grouped by batch and technical-replicate group. Three observations are required for threshold assessment. Filtering is disabled by default. If explicitly enabled, transcripts above the configured CV threshold are excluded from both models. This rule requires total-input bases and a declared replicate set; selection of different samples can change the assessed set and must be reported. Diagnostics are not a standalone validation of concentration precision.

Two concentration layers are exported:

1. `estimated_*`: the mathematical model estimates, retained with qualification flags when an accepted model exists.
2. `concentration_*`: values satisfying the selected reporting policy. Default `validated` reporting requires externally verified alignment QC, detection/background QC, calibration transfer and a positive independently established target LoQ. Target-specific non-missing inputs override sample-level inputs. Values below LoQ are suppressed. Missing validation declarations or LoQ do not prevent export of clearly identified estimates.

`exploratory` reporting relaxes the external-validation requirement but retains status flags. Calibration acceptance, positive signal, molecular-mass requirements, supplied LoQ and range checks still apply. Extrapolation requires explicit opt-in and remains flagged. A declaration set to true records user-provided evidence; the software does not independently certify the experiment.

Known ERCC rows contain assigned reaction concentrations, not predicted target concentrations, and are identified as `assigned_standard`. Their wastewater concentration is always unavailable. An unspiked sample receives no borrowed calibration model.

## Reproducibility and analytical validation

Input tables use UTF-8 TSV, explicit units and `NA` for missing values. Source rows and fields are retained. Inputs and model settings are snapshotted; run and sample manifests record checksums, software versions, alignment basis and equation diagnostics. Nonempty output directories are rejected to prevent mixed runs. The combined database and upstream software version should be fixed and reported separately.

The implementation provides point estimates. It does not propagate uncertainty in stock concentration, pipetting, extraction volume, recovery or calibration transfer. Validation of absolute accuracy, target-specific LoD/LoQ, capture bias, inhibition, background, recovery and reproducibility requires independent experimental evidence, preferably including orthogonal quantification. High calibration R² and passing software tests do not establish these properties.

## Primary reference sources

- [EsViritu](https://github.com/cmmr/EsViritu) and [coverage implementation, v1.1.2](https://github.com/cmmr/EsViritu/blob/v1.1.2/src/EsViritu/esv_funcs.py).
- [Manufacturer ERCC concentration assignments](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt), distributed verbatim with a recorded SHA256.
- [NIST SRM 2374 reference record](https://www-s.nist.gov/srmors/view_detail.cfm?srm=2374), with sequence provenance retained in `references/ercc96`.
