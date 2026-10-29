# EsViritu custom database v3.4

This release combines the original v3.3 reference sequences, the verified v3.3.1 metadata correction and 21 reviewed literature diagnostic-region references.

- 31,003 baseline records + 21 PCRROI records = 31,024 references; all 96 ERCC sequences are retained.
- Four added SARS-CoV-2 regions: N1, N2, E and ORF1ab.
- Original FASTA bytes are preserved. Only the eight previously documented influenza B assembly-length cells are revised before the local metadata is appended.
- The complete literature audit retains 35 assays and 34 unique windows; 13 boundary/mismatch candidates are excluded from the alignment FASTA. No new rotavirus A ROI is added because its reviewed windows lack a complete flank.
- Minimap2 2.31-r1302 short-read index, k=21 and w=11, is rebuilt. FASTA, metadata and index reference IDs/lengths match. All 96 ERCC and 21 ROI verification queries recover expected references.
- All 21 ROI verification queries have multiple reference hits. Genomic and overlapping ROI observations must not be summed or treated as independent measurements.
- Molecular annotations retain all baseline assignments; local ROI mass factors remain NA unless explicitly assigned in a reviewed target manifest. Local depth predictions are qualified response estimates, not validated genome concentrations.
- Complete source provenance, coordinates, literature, oligo audit, independent validators and English/Chinese computational methods accompany the archive.

Download `esviritu_db_v3.4.tar.gz` and `esviritu_db_v3.4.sha256`, run `sha256sum -c esviritu_db_v3.4.sha256`, then extract the archive. Point EsViritu `--db` to the extracted `v3.4` directory. New reference sequences require new alignments. No legacy annotation-only upgrade is required for this archive.

Software remains 2.2.0. Both primary calibration models retain R²≥0.80, positive slopes and at least four adjacent ERCC levels; the yield slope range 0.8–1.2 is diagnostic only by default. See the [workflow](https://github.com/Lt1020/HybCapSeq-Quant/blob/main/docs/WORKFLOW_v3.4_zh.md) and [manuscript methods](https://github.com/Lt1020/HybCapSeq-Quant/blob/main/docs/MANUSCRIPT_METHODS.md).

Database consistency and simulated numerical checks do not establish experimental sensitivity, absolute accuracy, LoD/LoQ, commercial assay identity or regional calibration transfer. The original viral database retains Michael Tisza's CC BY 4.0 attribution (DOI 10.5281/zenodo.17716199); ERCC and public molecular reference sources retain their own attribution and applicable terms.
