# Molecular annotation sidecar

## Identity and interpretation

`virus_pathogen_database.molecular_annotations.tsv` contains one row per exact upstream `Accession`. It is independent of the original metadata schema. The companion manifest binds its checksum, source workbook, saved NCBI summaries and matching database metadata. The packaged compressed table is installed beside the alignment database using `install-references --db-dir`.

The selected chemistry describes a viral genome-equivalent or explicitly identified ERCC RNA transcript. The deposited NCBI record and the ICTV genome category remain separate observations. Chemistry assignments do not establish the physical composition of an individual sample, calibration transfer or experimental quantification accuracy.

## Fields

| Field | Definition |
|---|---|
| `accession` | Exact database identifier, including version where present |
| `assembly` | Unchanged database `Assembly` |
| `reference_name` | Unchanged database `Name` |
| `length_nt` | Represented-strand reference length; numerical bp length for duplex genomes |
| `assembly_length_nt` | Declared reference-group length, including published checksum-bound corrections; historical input values remain in source detection fields |
| `family_raw` | Original family annotation, including rank prefix |
| `family_normalized` | Prefix-normalized lookup value; missing remains NA |
| `ncbi_taxid` | Taxonomy identifier in the retrieved accession summary |
| `ncbi_accession_version` | Exact returned accession version; a newer version is not substituted |
| `ncbi_length_nt` | Length in the retrieved NCBI summary |
| `ncbi_molecule_type` | Original NCBI record molecule type, distinct from selected genome chemistry |
| `ncbi_strandedness` | Original NCBI record strand annotation; absent values remain NA |
| `ncbi_biomol` | Original NCBI biomolecule annotation |
| `ncbi_record_update` | Update date reported by NCBI |
| `ncbi_retrieval_status` | Exact-version retrieval status |
| `ncbi_source_url` | Nuccore accession URL |
| `ncbi_response_sha256` | SHA256 of the cached complete ESummary response batch |
| `ictv_genome_type` | VMR genome categories supporting the assignment, including polarity/RT notation |
| `ictv_match_rank` | Closest successful matching level: accession, species, genus, family or order |
| `ictv_match_taxon` | Matched accession or taxon value |
| `ictv_isolate_ids` | Supporting VMR isolate IDs for small evidence sets; NA for sets larger than eight records |
| `ictv_supporting_record_count` | Number of supporting VMR records |
| `ictv_evidence_sha256` | SHA256 of the ordered genome/chemistry/isolate evidence, serialized with sorted JSON keys |
| `ictv_version` | Fixed official VMR version |
| `ictv_source_url` | Official version-specific workbook URL |
| `genome_nucleic_acid_type` | RNA or DNA when established; can remain known when strandedness is unresolved |
| `molecule_type` | Selected ssRNA, dsRNA, ssDNA, dsDNA, partially_dsDNA or NA |
| `strandedness` | ss, ds, partially_ds or NA |
| `molecular_mass_per_nt_g_mol` | Declared average mass factor per reference position, g/mol; NA if unresolved |
| `reference_molecular_mass_g_mol` | Approximate theoretical molar mass of this individual reference, factor × length |
| `reference_mass_g_per_molecule` | Individual-reference theoretical mass divided by Avogadro constant, g |
| `reference_mass_ng_per_molecule` | Individual-reference mass in ng |
| `mass_model` | `average_reference_position_v1` when a mass is assigned |
| `mass_model_source` | Source documenting the average conversion convention; dsRNA is an explicit two-strand derivation |
| `classification_evidence` | Evidence route used for selected genome/transcript chemistry |
| `annotation_status` | Resolution, evidence-disagreement or missing-mass reason |
| `notes` | Interpretation qualifications, including noncommercial ERCC quantity status |
| `assembly_reference_count` | Number of references in the complete database assembly |
| `assembly_molecular_mass_g_mol` | Sum of reference masses if complete length and chemistry checks pass; otherwise NA |
| `assembly_mass_status` | Group-level completeness/coherence status |

A reference row and an assembly row are different measurement objects. A detected segment may have a known individual-reference mass while its complete group mass remains unavailable. Quantification uses complete group length and mass. When an assembly appears in multiple source rows, its concentration is repeated; these concentrations must not be summed.

## Status interpretation

| Status | Interpretation |
|---|---|
| `resolved_ictv_genome` | ICTV genome chemistry assigned after exact NCBI accession/length verification |
| `resolved_ictv_genome_record_chemistry_differs` | Genome chemistry assigned; deposited NCBI record has a different RNA/DNA representation |
| `resolved_ictv_genome_record_strandedness_differs` | Genome chemistry assigned; deposited record has different strandedness |
| `resolved_ncbi_accession` | Explicit NCBI molecule and strand annotations used when an ICTV type is unavailable |
| `assigned_ERCC_transcript_reference` | RNA transcript mass model; concentration assignment is a separate Mix 1 decision |
| `partially_double_stranded_DNA_requires_material_assignment` | Partial duplex completion prevents a unique genome mass |
| `unresolved_molecular_type` | Strand count/type is not established; no default mass is assigned |
| `exact_accession_version_unavailable` | Requested version was not returned; no newer record is silently used |
| `ncbi_database_length_mismatch` | NCBI and database lengths differ; mass assignment is blocked |
| `resolved_complete_assembly` | Sum of member lengths equals declared complete length, and member mass factors agree |
| `unresolved_assembly_member_mass` | At least one full-database assembly member lacks mass |
| `inconsistent_assembly_reference_length` | Sum of full-database reference lengths differs from declared complete length |
| `inconsistent_assembly_molecular_chemistry` | Complete-group members do not share a coherent mass factor |

An explicit sourced target manifest can supply a material-specific mass model. This records a user declaration, not automatic experimental confirmation. A catalog length inconsistency remains blocked in both estimators until the database/reference definition has been reviewed and a consistent catalog generated.

## Reproducibility

The complete saved NCBI summary snapshot contains one JSON object per viral accession and the checksum of the retrieval batch. It is a selected-field provenance snapshot, not a verbatim archival copy of every GenBank feature record. Offline rebuilding accepts this snapshot and the checksum-verified official VMR workbook. The source URLs and hashes permit independent verification. Creation timestamps are provenance metadata and are not part of the quantitative model.
