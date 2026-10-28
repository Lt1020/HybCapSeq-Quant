# Software verification

## Analytical fixture

`examples/quant/generate_demo.py` generates explicitly simulated reference observations under a continuous, noiseless response model:

$$
n_j=10^9m_j,
$$

where mass is expressed in ng. The fixture includes the 92 commercial Mix 1 standards, a single RNA reference with 10,000 input molecules, and a two-segment RNA reference with the same molecule count. Segment lengths sum to the locked complete reference length. Read-count and breadth columns support input-contract tests; the continuous signal is an analytical fixture, not simulated raw sequencing reads or experimental data.

With a 14 µL reaction, a 2 µL extract aliquot and the supplied synthetic conversion parameters, the expected target concentrations are:

| Quantity | Analytical value |
|---|---:|
| Reaction concentration | 10,000 / 14 copies/µL |
| Extract concentration | 5,000 copies/µL |
| Wastewater concentration | 61,333.3333333 copies/mL |

These values are recovered independently by both estimators. Validation flags and LoQ in the fixture are synthetic assertions used to exercise reporting; they do not represent laboratory validation.

## Verification domains

The automated suite checks official stock assignments and original-source checksums, reference lengths, dilution and Avogadro conversion, inclusion of all spiked standard mass, arithmetic-mean yield estimation, normalized-depth calibration, target molecular mass, segment aggregation, mass-balance recovery and unit conversion. It checks positive yield slopes outside the diagnostic interval, rejection when the historical gate is explicitly enabled, adjacent-level requirements, non-extrapolation, missing inputs, reporting policy, target-specific LoQ, unspiked controls and optional replicate filtering.

Command-line checks cover single-sample and batch execution, source-column preservation, English output fields, input/output conflict detection, relative paths, relocation, directory names with spaces, independent working directories, installed schemas, provenance and figure generation. Optional BAM checks use a constructed alignment file containing known duplicate, secondary, supplementary and multi-mapping records.

GitHub Actions runs the suite and analytical example under Python 3.10, 3.11 and 3.12. Each job retains its generated result files and plots as an artifact. The status for a particular analysis must be taken from its exact commit, not inferred from the project badge alone.

## Molecular annotation verification

The suite validates the packaged full-catalog checksum, 31,003 unique accession keys, mass arithmetic and declared chemistry factors. It checks exact source length, preservation of NCBI/ICTV disagreement, prefix-normalized taxon lookup, unassigned unknown strandedness, accession-specific dsRNA assignment, rejection of unrelated database metadata and checksum tampering, and refusal to derive a group mass from only its known segment. The analytical fixture supplies an explicitly simulated chemistry manifest rather than relying on viral taxonomy labels.

These checks establish the integrity of the distributed source-to-model mapping; they do not independently validate every deposited taxonomic assignment or establish chemical composition in a study sample.

## Reference revision and reporting-mode verification

Version 2.1.0 tests staged database upgrades, byte-level preservation of unrelated metadata and line endings, unchanged FASTA/index files, original-file backups, no-write inspection, idempotence, rejection of unknown changes, recovery from a simulated write failure, exact legacy-reference compatibility and continued rejection of other mismatches. Catalog validation independently recomputes full-group length, member count, status and mass arithmetic. Exploratory-mode tests verify that a positive out-of-range model estimate is retained without assigning an ERCC endpoint as a LoQ; explicit non-extrapolation remains available.

The aggregate archived-batch regression in `validation/reference_revision_2.1.0.json` covers 13 libraries. Each estimator has 1,603 sample/reference-group observations: seven documented corrected groups become computable, while estimates for the other 1,596 groups and all calibration coefficients remain unchanged. Original detection fields remain byte-equivalent as TSV cell values. These are software regression results on existing observations, not a new experimental accuracy validation. `scripts/verify_reference_revision.py` reproduces this audit from the original result archive and a revised reference directory; no sample measurements are distributed in this repository.

## Version 2.2.0 additions

Correction fixtures have independent analytical factors of 10 and 5 for the depth and yield observations. Tests check exact factor recovery, incompatible shared factors, unchanged parent weights under extra dilution observations, grouped training/test separation, explicit exclusion audits, duplicate-key rejection, insufficient training parents, source-value preservation and end-to-end command execution. Policy tests check the 0.80 default and explicit legacy slope gate. Reference-preflight fixtures are nonpathogen synthetic padding and verify that input checks do not claim a database release. No experimental target measurements enter public test fixtures.

Local checks and GitHub Actions results refer to a specific commit. Optional BAM tests are skipped when pysam is absent; that skip must remain visible. Passing source checks does not establish empirical accuracy. Database v3.4 remains blocked until actual assay-region information and reviewed reference files are supplied; no v3.4 sequence/index validation has been performed.

## Reproduction commands

```bash
python -m pip install '.[bam]'
python -m unittest discover -s tests -v
python examples/quant/generate_demo.py
hybcapseq_quant run \
  --samples examples/quant/samples.tsv \
  --esviritu-dir examples/quant \
  --targets examples/quant/targets.tsv \
  --outdir analytical_example_results
```

The output directory must not already contain another analysis. Tests operate in temporary directories.

## Interpretation

These checks establish mathematical and software behavior for the tested inputs. They do not establish target-specific trueness, precision, LoD, LoQ, unbiased RNA-to-DNA calibration transfer, or representativeness of PMMoV recovery. No experimental performance metric is inferred from the fixture.
