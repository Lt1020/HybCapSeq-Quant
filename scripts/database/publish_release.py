#!/usr/bin/env python3
"""Publish only verified v3.3 assets; never overwrite an existing different asset."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

TAG = "db-v3.3"
TITLE = "EsViritu DB v3.3 + ERCC96"


def gh(*args):
    return subprocess.check_output(["gh", *args], text=True)


def release(repo):
    result = subprocess.run(["gh", "api", f"repos/{repo}/releases/tags/{TAG}"],
                            text=True, capture_output=True)
    if result.returncode == 0:
        return json.loads(result.stdout)
    if "404" in result.stderr:
        return None
    raise RuntimeError("Could not query release; check GitHub Actions permissions")


def digest(file):
    h = hashlib.sha256()
    with file.open("rb") as src:
        while block := src.read(1024 * 1024):
            h.update(block)
    return "sha256:" + h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assets", type=Path, required=True)
    p.add_argument("--notes", type=Path, required=True)
    a = p.parse_args()
    repo = os.environ["GITHUB_REPOSITORY"]
    if repo != "Lt1020/HybCapSeq-Quant":
        raise ValueError("Publishing is restricted to the authorized repository")
    names = ["esviritu_db_v3.3.tar.gz", "esviritu_db_v3.3_manifest.json",
             "README_v3.3.md", "esviritu_db_v3.3.sha256"]
    expected = {name: {"size": (a.assets / name).stat().st_size,
                       "digest": digest(a.assets / name)} for name in names}
    manifest = json.loads((a.assets / names[1]).read_text())
    if expected[names[0]]["digest"] != "sha256:" + manifest["archive"]["sha256"]:
        raise ValueError("Unverified database archive")
    current = release(repo)
    if current is None:
        gh("release", "create", TAG, "--repo", repo, "--target", os.environ["GITHUB_SHA"],
           "--title", TITLE, "--notes-file", str(a.notes), "--draft")
        current = release(repo)
    if current["name"] != TITLE:
        raise ValueError("Existing tag is used by a different release")
    for name in names:
        assets = {x["name"]: x for x in current["assets"]}
        existing = assets.get(name)
        if existing:
            if any(existing.get(k) != v for k, v in expected[name].items()):
                raise ValueError(f"Existing asset differs; refusing to overwrite: {name}")
        else:
            if not current["draft"]:
                raise ValueError("Published release is missing an expected asset")
            gh("release", "upload", TAG, str(a.assets / name), "--repo", repo)
            current = release(repo)
    assets = {x["name"]: x for x in current["assets"]}
    for name, fields in expected.items():
        if any(assets[name].get(k) != v for k, v in fields.items()):
            raise ValueError(f"Uploaded asset verification failed: {name}")
    if current["draft"]:
        gh("release", "edit", TAG, "--repo", repo, "--draft=false")
        current = release(repo)
    if current["draft"]:
        raise ValueError("Release was not published")
    print(json.dumps({"release_url": current["html_url"], "assets": [
        {k: x[k] for k in ["name", "size", "digest", "browser_download_url"]}
        for x in current["assets"]]}, indent=2))


if __name__ == "__main__":
    main()
