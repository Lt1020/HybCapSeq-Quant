# Target-specific bias assessment

## Input and estimand

The correction is an empirical multiplicative adjustment of extract-equivalent concentration. It does not change aligned reads, establish LoD/LoQ or identify a physical cause of bias. Comparisons require the same specimen, concentration layer and units.

The UTF-8 TSV requires `batch`, `sample_id`, `parent`, `extraction_id`, `method`, `target`, `seq_extract_copies_ul`, `dpcr_extract_copies_ul`, `complete_in_range`, `dpcr_comparator_eligible`, and `model_accepted`. True flags are `True` or `1`. All three eligibility flags must be true and both concentrations positive for fitting. Eligibility must be established independently of the observed error. Batch/sample/method/target keys must be unique. Additional columns are preserved. An external method label such as RPKMF can be assessed separately; this does not add that estimator to `run`.

## Coefficient

For sequencing estimate x and dPCR concentration y, calculate log10(y/x). Take its median across dilution observations within batch/parent/extraction/method; take the median of extraction-method summaries within parent; take the median across parents. The factor is 10 raised to this final median. Shared depth/yield profiles summarize those two methods together, without RPKMF. Additional dilution observations do not automatically increase a parent's weight.

Fit scopes include each batch and all batches pooled. Parent identifiers are qualified by batch. Users must choose identifiers reflecting actual specimen independence; repeated collection from a site is not assumed to represent the same physical specimen. A deterministic parent bootstrap (1,000 draws by default) provides a conditional factor interval; fewer than two parents produces no interval. This is not a prediction interval and does not propagate all reference or laboratory uncertainties.

## Exact 30% feasibility

A fixed positive K gives every positive pair APE ≤30% if and only if the intervals [0.7 y_i/x_i, 1.3 y_i/x_i] intersect. The tool reports the maximum lower bound, minimum upper bound and feasibility flag. It does not adjust factors per sample. Feasibility on existing data does not guarantee future performance.

## Independent prediction

Leave-one-parent-out excludes every observation of the test parent, including repeated extractions and ERCC conditions. Each training split requires at least three parents. Insufficient splits are recorded; in-sample fallback is not used. Cross-batch assessment trains on one batch and predicts other batches when training is sufficient. Predictions retain training IDs, coefficients and errors.

Summaries report n, parent count, median/maximum APE, the fraction within 30%, mean log10 bias, log10 RMSE, R² and Lin CCC on log10 concentrations. Library-target observations are not claimed to be independent biological samples. CCC uses sample covariance and variances.

`all_heldout_pairs_within_30pct` requires complete eligible leave-parent-out coverage, at least four evaluated parents, and all held-out APE values ≤30%. `heldout_target_not_met` means adequate evaluation failed this criterion; `insufficient_validation` indicates missing or insufficient independent evaluation. `formal_assay_validation` remains false for all empirical coefficients.

## Application

`apply-correction` selects the exact fit scope and method-specific factor. Shared profiles are assessed but not implicitly substituted. A positive source estimate, accepted calibration, complete in-range composition, and a raw estimate inside the factor's training interval are required. The future observation does not need a positive dPCR result; its dPCR field may be `NA`.

Unless `--allow-exploratory` is explicit, only `all_heldout_pairs_within_30pct` factors are applied. Other values remain missing with a reason. Source values are preserved. A manifest fingerprints input, coefficients, software and options. Application to training observations is descriptive; accuracy claims must use held-out predictions.

## Outputs

| File | Content |
|---|---|
| `correction_parameters.tsv` | Factors, conditional intervals, ranges, counts and status |
| `correction_feasibility.tsv` | Exact observed-data ±30% feasibility |
| `correction_predictions.tsv` | Held-out observations and training provenance |
| `correction_summary.tsv` | Before/after held-out metrics |
| `correction_input_audit.tsv` | Every original pair and exclusions |
| `correction_validation_gaps.tsv` | Insufficient validation splits |
| `correction_run_manifest.json` | Settings and source fingerprint |

Experimental measurements are not distributed in the public repository. Examples are explicitly simulated.
