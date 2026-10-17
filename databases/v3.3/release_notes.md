Official EsViritu DB v3.2.4 + 96 NIST SRM2374 certified-v1 ERCC references.

- 30,907 original reference records + 96 ERCC records = 31,003 records.
- Original viral FASTA and metadata bytes are preserved; ERCC records are appended.
- ERCC Assembly = Name; Asm_length = Length. Unsupported ERCC metadata remains NA.
- The three original database filenames and 16-column metadata schema are preserved.
- Short-read minimap2 index rebuilt using minimap2 2.31-r1302, k=21, w=11.
- All FASTA/metadata/index IDs and lengths match; all 96 ERCC verification queries match their expected references.
- Archive SHA256: `86e3f472b5376bb75c0baac29f460f3e13f4766f091659107c17e04451198e2c`.
- Archive size: 447,693,545 bytes. The release archive is byte-identical to the delivered v3.3 package.

This is a custom ERCC extension, not an upstream EsViritu release. No complete EsViritu run on the experimental FASTQ files is claimed.

Upstream: Michael Tisza, *Virus Pathogen Database for EsViritu*, v3.2.4, DOI: [10.5281/zenodo.17716199](https://doi.org/10.5281/zenodo.17716199), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The changes are the addition of the ERCC records and regeneration of the index.

ERCC source: [NIST SRM 2374 certified sequence v1](https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA).

Download `esviritu_db_v3.3.tar.gz`, verify it with the SHA256 file, and extract it:

```bash
sha256sum -c esviritu_db_v3.3.sha256
tar -xzf esviritu_db_v3.3.tar.gz
```

Point EsViritu `--db` to the extracted `v3.3` directory. See `README_v3.3.md` and `esviritu_db_v3.3_manifest.json` for provenance and validation details.
