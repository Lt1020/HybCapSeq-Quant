# Software verification

## Analytical fixture

`examples/quant/generate_demo.py` generates explicitly simulated reference observations under a continuous, noiseless response model:

\[
n_j=10^9m_j,
\]

where mass is expressed in ng. The fixture includes the 92 commercial Mix 1 standards, a single RNA reference with 10,000 input molecules, and a two-segment RNA reference with the same molecule count. Segment lengths sum to the locked complete reference length. Read-count and breadth columns support input-contract tests; the continuous signal is an analytical fixture, not simulated raw sequencing reads or experimental data.

With a 14 µL reaction, a 2 µL extract aliquot and the supplied synthetic conversion parameters, the expected target concentrations are:

| Quantity | Analytical value |
|---|---:|
| Reaction concentration | 10,000 / 14 copies/µL |
| Extract concentration | 5,000 copies/µL |
| Wastewater concentration | 61,333.3333333 copies/mL |

These values are recovered independently by both estimators. Validation flags and LoQ in the fixture are synthetic assertions used to exercise reporting; they do not represent laboratory validation.

## Verification domains

The automated suite checks official stock assignments and original-source checksums, reference lengths, dilution and Avogadro conversion, inclusion of all spiked standard mass, arithmetic-mean yield estimation, normalized-depth calibration, target molecular mass, segment aggregation, mass-balance recovery and unit conversion. It also checks failed slope criteria despite perfect R², adjacent-level requirements, non-extrapolation, missing inputs, reporting policy, target-specific LoQ, unspiked controls and optional replicate filtering.

Command-line checks cover single-sample and batch execution, source-column preservation, English output fields, input/output conflict detection, relative paths, relocation, directory names with spaces, independent working directories, installed schemas, provenance and figure generation. Optional BAM checks use a constructed alignment file containing known duplicate, secondary, supplementary and multi-mapping records.

GitHub Actions runs the suite and analytical example under Python 3.10, 3.11 and 3.12. Each job retains its generated result files and plots as an artifact. The status for a particular analysis must be taken from its exact commit, not inferred from the project badge alone.

## Reproduction

```bash
python -m pip install '.[bam]'
python -m unittest discover -s tests -v
python examples/quant/generate_demo.py
hybcapseq_quant run \
  --samples examples/quant/samples.tsv \
  --esviritu-dir examples/quant \
  --outdir analytical_example_results
```

The output directory must not already contain another analysis. Tests operate in temporary directories.

## Interpretation

These checks establish mathematical and software behavior for the tested inputs. They do not establish target-specific trueness, precision, LoD, LoQ, unbiased RNA-to-DNA calibration transfer, or representativeness of PMMoV recovery. No experimental performance metric is inferred from the fixture.
