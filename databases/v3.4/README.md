# EsViritu custom database v3.4

Version v3.4 appends **21 reviewed literature diagnostic reference regions** to the 31,003-record v3.3 database, retains the verified v3.3.1 metadata correction, and rebuilds the short-read minimap2 index. The combined database contains **31,024 records**, including 96 ERCC references. This is a project database revision, independent of HybCapSeq-Quant software 2.2.0 and upstream EsViritu versions.

The four added SARS-CoV-2 references represent CDC N1, CDC N2, Corman E and China CDC ORF1ab. These counts refer to added PCRROI records; the original genomic-reference records remain present.

## Download and installation

The complete alignment database is distributed as a [Release asset](https://github.com/Lt1020/HybCapSeq-Quant/releases/tag/db-v3.4), rather than committed to Git history.

```bash
curl --fail --location --output esviritu_db_v3.4.tar.gz \
  https://github.com/Lt1020/HybCapSeq-Quant/releases/download/db-v3.4/esviritu_db_v3.4.tar.gz
curl --fail --location --output esviritu_db_v3.4.sha256 \
  https://github.com/Lt1020/HybCapSeq-Quant/releases/download/db-v3.4/esviritu_db_v3.4.sha256
sha256sum -c esviritu_db_v3.4.sha256
tar -xzf esviritu_db_v3.4.tar.gz
cd v3.4
sha256sum -c SHA256SUMS
cd ..
hybcapseq_quant validate-db --db-dir "$PWD/v3.4"
```

Point the installed EsViritu `--db` option to the extracted `v3.4` directory. Re-run alignment when moving from v3.3/v3.3.1: an annotation update cannot create observations against new sequences. Do not run the legacy `upgrade-db` or `install-references` commands on v3.4; the matching annotations are included.

## Contents and provenance

The original filenames remain `virus_pathogen_database.fna`, `virus_pathogen_database.all_metadata.tsv` and `virus_pathogen_database.mmi`. Metadata retains its 16-column schema and a single header. Molecular annotations, the database manifest, coordinate mapping, query audit, checksums and the complete frozen PCR reference package accompany these files.

The original v3.3 FASTA is retained byte-for-byte as a prefix. Only the eight already documented Rochester influenza B `Asm_length` cells receive the v3.3.1 correction; the corrected metadata is then retained as a prefix. All 96 ERCC sequences and their grouping are preserved. The index is rebuilt with minimap2 **2.31-r1302**, preset `sr`, k=21 and w=11. The original index is not reused after appending sequences.

The frozen literature package covers 12 targets and 35 published assays, yielding 34 unique sequence windows. Its main set contains 21 records across 11 targets with exactly 100 nucleotides outside each amplicon boundary and compatible selected published oligo sites. Seven boundary-limited and six mismatch candidates remain in `pcr_roi_reference_package/exceptions/`; they are not appended to the alignment FASTA. All three reviewed rotavirus A windows are boundary-limited, so no rotavirus PCRROI is added. The original rotavirus references remain available.

FASTA identifiers follow `PCRROI_<target>_<assay>_<accession.version>_<start>_<end>`. PCRROI means *PCR region of interest*, a project prefix. Coordinates are plus-strand, one-based inclusive source coordinates. The amplicon includes both primer-binding sites; 100 nt are added outside those boundaries. Every emitted base is retained from a fixed NCBI accession version. No padding, oligo substitution or cross-strain sequence joining is used. `Length=Asm_length` equals the actual local fragment length; each unique ROI has its own `Name=Assembly`.

The [source package](../../references/pcr_roi/README_zh.md) retains the literature, accession roles, coordinates, oligo audit, exceptions and reconstruction scripts. A published assay is not evidence that a commercial dPCR kit uses the same region. Actual commercial assay identity and experimental calibration transfer remain unconfirmed.

## Counting and quantitative interpretation

Genomic references and their local regions overlap. Index verification explicitly permits secondary alignments and records how many reference hits are reported. An expected-reference hit verifies recoverability, not unique mapping, sensitivity or assay specificity. Competitive primary alignments may be assigned to a genomic reference instead of a PCRROI. A zero PCRROI count therefore does not establish absence of the region.

Do not add genomic and ROI counts, sum overlapping assays, or count PCRROI labels as additional viral species. Preserve reference scope in all comparisons. For region-restricted abundance, define a read/fragment assignment policy from compatible alignments before analysis and audit ambiguous mappings; the database alone does not implement a new deduplicated ROI counting method.

The molecular sidecar retains all original genome-reference assignments. Added ROI molecular-mass factors are `NA` by default, because a local response must not inherit a whole-genome model. Software 2.2.0 can retain their observed signals and compute explicitly qualified depth estimates. Yield concentrations for these ROIs remain unavailable unless a reviewed `--targets` table supplies an explicit local-molecule mass assignment. Any such estimate is a local-reference equivalent, not an independently validated genome concentration. Calibration-transfer declarations and LoQ must refer to that same measurement definition. Existing whole-reference correction factors are not transferred to a new region or new alignment run.

Both implemented calibration models continue to require positive slope, R² ≥0.80 and at least four adjacent assigned ERCC levels. The yield slope interval 0.8–1.2 remains diagnostic by default. RPKMF remains an observed metric, not a third estimator in `run`. See the [computational workflow](../../docs/WORKFLOW_v3.4_zh.md) and [methods](../../docs/METHODS.md).

## Rebuild and verification

With the pinned v3.3 database extracted locally and minimap2 2.31 available:

```bash
python scripts/database/build_v3_4.py \
  --base-dir /path/to/v3.3 \
  --reference-package references/pcr_roi \
  --output-dir /path/to/new/v3.4 \
  --work-dir /path/to/build-work \
  --minimap2 /path/to/minimap2 \
  --expected-manifest databases/v3.4/esviritu_db_v3.4_manifest.json
hybcapseq_quant validate-db --db-dir /path/to/new/v3.4
```

The build rejects unrecognized base checksums, altered source-package files, identifier/group collisions, inconsistent metadata and annotation lengths, wrong primary-set membership, and an index whose references differ from the FASTA. All 96 ERCC verification queries and 21 ROI verification queries must recover their expected references. Release construction additionally verifies every archived file after decompression. The manifest distinguishes these computational checks from the experimental EsViritu/dPCR validation that has not been performed.

## Attribution

The base viral database is Michael Tisza's *Virus Pathogen Database for EsViritu*, v3.2.4, [DOI 10.5281/zenodo.17716199](https://doi.org/10.5281/zenodo.17716199), CC BY 4.0. Changes comprise ERCC and literature-region additions, the documented v3.3.1 metadata correction and a rebuilt index. ERCC sequences are from NIST SRM 2374 certified v1. Public diagnostic method citations are recorded individually in `literature_assays.tsv`; article PDFs are not redistributed. NCBI molecular sequence data are used with accession/version attribution under its [molecular data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/). NCBI does not itself restrict these data, but this is not a claim that third-party rights do not exist. The software MIT license does not relicense third-party reference material.
