#!/usr/bin/env python3
"""Recreate the delivered v3.3 archive byte-for-byte, with pinned gzip metadata."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
import zlib


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--database-dir", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--recipe", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    recipe = json.loads(a.recipe.read_text())
    if sys.version_info[:2] != (3, 12):
        raise ValueError("Use Python 3.12 for identical tar/gzip serialization")
    if zlib.ZLIB_RUNTIME_VERSION != recipe["zlib_runtime_version"]:
        raise ValueError("zlib runtime version differs from original archive")
    expected = {f["name"]: f for f in manifest["files"]}
    for name, info in expected.items():
        file = a.database_dir / name
        if file.stat().st_size != info["bytes"] or sha(file) != info["sha256"]:
            raise ValueError(f"Database file differs from delivered v3.3: {name}")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("wb") as raw:
        with gzip.GzipFile(filename=recipe["gzip_filename"], fileobj=raw, mode="wb",
                           compresslevel=recipe["compresslevel"], mtime=recipe["gzip_mtime"]) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for record in recipe["members"]:
                    member = tarfile.TarInfo(record["name"])
                    for key, value in record.items():
                        if key not in {"name", "type"}:
                            setattr(member, key, value)
                    member.type = record["type"].encode("ascii")
                    path = PurePosixPath(member.name)
                    if path.parts[0] != "v3.3" or ".." in path.parts:
                        raise ValueError("Unexpected archive member")
                    if member.isdir():
                        archive.addfile(member)
                    elif member.isfile() and len(path.parts) == 2 and path.name in expected:
                        with (a.database_dir / path.name).open("rb") as src:
                            archive.addfile(member, src)
                    else:
                        raise ValueError("Unexpected archive member type or name")
    digest = sha(a.output)
    if a.output.stat().st_size != manifest["archive"]["bytes"] or digest != manifest["archive"]["sha256"]:
        raise ValueError(f"Archive differs from original delivered package: {digest}")
    print(json.dumps({"archive": str(a.output), "bytes": a.output.stat().st_size,
                      "sha256": digest, "identical_to_delivered_archive": True}, indent=2))


if __name__ == "__main__":
    main()
