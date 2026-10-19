"""Build deterministic, sequence-free v3.3.1 reference update assets."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

from hybcapseq_quant.reference_updates import (
    CORRECTIONS_NAME,
    bundled_reference_files,
    read_corrections,
)
from hybcapseq_quant.references import CATALOG_NAME, MANIFEST_NAME


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalog, manifest = bundled_reference_files()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", filename="", mtime=0) as dest:
        dest.write(catalog)
    payloads = {
        CATALOG_NAME + ".gz": buffer.getvalue(),
        MANIFEST_NAME: (json.dumps(manifest, indent=2) + "\n").encode(),
        CORRECTIONS_NAME: (json.dumps(read_corrections(), indent=2) + "\n").encode(),
        "README.md": (
            Path(__file__).resolve().parents[2] / "docs/UPGRADE_2.1.0.md"
        ).read_bytes(),
    }
    for name, body in payloads.items():
        if name != "README.md":
            (args.output_dir / name).write_bytes(body)
    archive_path = args.output_dir / "esviritu_db_v3.3.1_update.tar.gz"
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(
                fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT
            ) as archive:
                for name in sorted(payloads):
                    body = payloads[name]
                    item = tarfile.TarInfo("reference_update/" + name)
                    item.size = len(body)
                    item.mode = 0o644
                    item.mtime = 0
                    archive.addfile(item, io.BytesIO(body))
    with tarfile.open(archive_path, "r:gz") as archive:
        for name, body in payloads.items():
            if archive.extractfile("reference_update/" + name).read() != body:
                raise ValueError("Reference update archive round-trip differs")
    print(
        json.dumps(
            {
                "archive": str(archive_path),
                "sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest(),
                "contains_alignment_sequences": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
