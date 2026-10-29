#!/usr/bin/env python3
"""Package a verified database and check every archived byte after decompression."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.database_dir.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    expected = {}
    for line in (source / "SHA256SUMS").read_text().splitlines():
        sha, name = line.split("  ", 1)
        path = (source / name).resolve()
        if not path.is_relative_to(source) or digest(path) != sha:
            raise ValueError("Database checksum mismatch: " + name)
        expected[name] = sha
    expected["SHA256SUMS"] = digest(source / "SHA256SUMS")
    paths = {
        p.relative_to(source).as_posix(): p
        for p in source.rglob("*")
        if p.is_file()
        and not any(part.startswith(".") for part in p.relative_to(source).parts)
    }
    if set(paths) != set(expected):
        raise ValueError("Unexpected or missing database package files")
    archive_path = output / "esviritu_db_v3.4.tar.gz"
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(
            fileobj=raw, filename="", mode="wb", mtime=0, compresslevel=6
        ) as gz:
            with tarfile.open(
                fileobj=gz, mode="w", format=tarfile.PAX_FORMAT
            ) as archive:
                for name, path in sorted(paths.items()):
                    member = tarfile.TarInfo("v3.4/" + name)
                    member.size = path.stat().st_size
                    member.mode = 0o644
                    member.mtime = 0
                    with path.open("rb") as handle:
                        archive.addfile(member, handle)
    checked = set()
    with tarfile.open(archive_path, "r|gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.startswith("v3.4/"):
                raise ValueError("Unexpected archive member")
            name = member.name.removeprefix("v3.4/")
            if name not in expected or name in checked:
                raise ValueError("Unlisted or duplicate archive member")
            h = hashlib.sha256()
            with archive.extractfile(member) as handle:
                while block := handle.read(1024 * 1024):
                    h.update(block)
            if h.hexdigest() != expected[name]:
                raise ValueError("Archived content changed: " + name)
            checked.add(name)
    if checked != set(expected):
        raise ValueError("Incomplete database archive")
    manifest = json.loads((source / "esviritu_db_v3.4_manifest.json").read_text())
    manifest["archive"] = {
        "filename": archive_path.name,
        "bytes": archive_path.stat().st_size,
        "sha256": digest(archive_path),
        "directory": "v3.4",
        "roundtrip_all_member_sha256_verified": True,
        "checked_file_count": len(checked),
    }
    (output / "esviritu_db_v3.4_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    (output / "esviritu_db_v3.4.sha256").write_text(
        manifest["archive"]["sha256"] + "  " + archive_path.name + "\n"
    )
    shutil.copy2(source / "README_v3.4.md", output)
    print(json.dumps(manifest["archive"], indent=2))


if __name__ == "__main__":
    main()
