# Database v3.4 preparation

**Status: incomplete; no v3.4 sequence database or alignment index is released.**

The intended extension adds reviewed references associated with the actual digital-PCR assay regions. The current input material identifies assay targets but does not supply their confirmed regions, reviewed local-reference FASTA files, or redistribution permissions. A commonly published assay must not be substituted for an undisclosed commercial assay.

The validated v3.3.1 metadata correction and its accession-resolved molecular annotations remain the active reference baseline. This directory is a preparation record, not a replacement database. Software 2.2.0 does not imply that database v3.4 exists.

## Required inputs

Complete `target_references.template.tsv` using existing reviewed files. Coordinates are one-based inclusive on the source accession. `reference_id` must match the supplied FASTA record exactly. Record the actual assay identity, source accession version, source URL, checksum, review status, and redistribution status. The flanking-region check only audits the supplied annotation; it does not obtain, infer, or generate sequence content.

```bash
python scripts/database/preflight_v3_4.py \
  --manifest reviewed_target_references.tsv \
  --base-dir /path/to/current_database \
  --output reference_preflight.json
```

An empty template or missing inputs produces exit code 2 and an explicit list of blocking issues. Even a passing input audit does not establish release readiness. A completed extension additionally needs compatible metadata, reviewed molecular annotations, unique identifiers, a matching rebuilt index, and end-to-end validation. Local-region signals require a separately defined measurement model and must not inherit whole-genome lengths or correction factors without validation.

## Release conditions

Only after all sequence, metadata, index, licensing, and functional checks pass should a `db-v3.4` release be created. The release must include file sizes, SHA256 checksums, build versions, changes relative to v3.3.1, and a tested installation command. Until then, continue to use the existing database and do not label old alignments as v3.4 results.
