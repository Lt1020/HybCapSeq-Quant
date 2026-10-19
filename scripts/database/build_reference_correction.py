"""Reproduce the v3.3.1 correction from the unchanged v3.3 reference files."""

import argparse
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path

from hybcapseq_quant.reference_updates import (
    CORRECTIONS_NAME,
    corrected_catalog_bytes,
    corrected_metadata_bytes,
    verify_fasta,
)
from hybcapseq_quant.references import (
    CATALOG_NAME,
    MANIFEST_NAME,
    METADATA_NAME,
    PACKAGE_REFERENCE,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--created-at", help="UTC timestamp for byte-identical provenance reproduction"
    )
    args = parser.parse_args()
    correction = {
        "correction_id": "rochester_0093_2013_reference_length_v3.3.1",
        "reference_name": "Influenza B virus (B/Rochester/0093/2013)",
        "assembly": "set:c1e60ab7f8b946bb1afa307f125370c361784dee",
        "original_assembly_length_nt": 15617,
        "corrected_assembly_length_nt": 13859,
        "members": {
            "KX351451.1": 1942,
            "KX351452.1": 2384,
            "KX351450.1": 2292,
            "KX254829.1": 1758,
            "KX351455.1": 1832,
            "KX351454.1": 1549,
            "KX254832.1": 1012,
            "KX351453.1": 1090,
        },
        "basis": "Sum of the eight distinct reference records verified against the existing FASTA and exact-accession NCBI length annotations.",
        "reference_length_interpretation": "Length of the represented reference collection; not a new claim about complete natural genome length.",
        "source_urls": [
            f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
            for acc in [
                "KX351451.1",
                "KX351452.1",
                "KX351450.1",
                "KX254829.1",
                "KX351455.1",
                "KX351454.1",
                "KX254832.1",
                "KX351453.1",
            ]
        ],
    }
    specification = {
        "schema_version": "1.0",
        "database_revision": "v3.3.1",
        "original_metadata_sha256": "bb406d8b177deede56ba3642fa2c5e7abfda56911f3fb3a91893a52da0a39ea3",
        "original_catalog_sha256": "efa18296a91ee008818c97e046e60cb6d6c8648697e27a9891fc33ff812e8c59",
        "fasta_sha256": "16a4a1b616d3cc9c45bd43f446aa561a4e6868b8f0b9861819b99cf402a2c8ed",
        "corrections": [correction],
    }
    assert (
        sum(correction["members"].values())
        == correction["corrected_assembly_length_nt"]
    )
    verify_fasta(args.source_db, specification)
    metadata = corrected_metadata_bytes(
        (args.source_db / METADATA_NAME).read_bytes(), specification
    )
    with gzip.open(PACKAGE_REFERENCE / (CATALOG_NAME + ".gz"), "rb") as stream:
        catalog = corrected_catalog_bytes(stream.read(), specification)
    specification["corrected_metadata_sha256"] = hashlib.sha256(metadata).hexdigest()
    specification["corrected_catalog_sha256"] = hashlib.sha256(catalog).hexdigest()
    manifest = json.loads((PACKAGE_REFERENCE / MANIFEST_NAME).read_text())
    manifest.update(
        created_at_utc=args.created_at or dt.datetime.now(dt.timezone.utc).isoformat(),
        catalog_sha256=specification["corrected_catalog_sha256"],
        database_metadata_sha256=specification["corrected_metadata_sha256"],
        database_revision="v3.3.1",
        source_database_metadata_sha256=specification["original_metadata_sha256"],
        source_catalog_sha256=specification["original_catalog_sha256"],
        reference_correction_ids=[correction["correction_id"]],
    )
    specification["corrected_catalog_manifest"] = manifest
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / CORRECTIONS_NAME).write_text(
        json.dumps(specification, indent=2) + "\n"
    )
    print(
        json.dumps(
            {k: v for k, v in specification.items() if k.endswith("sha256")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
