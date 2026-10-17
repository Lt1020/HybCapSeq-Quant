#!/usr/bin/env python3
"""Publish verified assets through release IDs, including unpublished drafts."""
import argparse
import hashlib
import http.client
import json
import os
from pathlib import Path
import time
from urllib.parse import urlencode, urlsplit

TAG = "db-v3.3"
TITLE = "EsViritu DB v3.3 + ERCC96"


def request(method, host, path, body=None, content_type="application/json", size=None):
    headers = {"Authorization": "Bearer " + os.environ["GH_TOKEN"],
               "Accept": "application/vnd.github+json", "X-GitHub-Api-Version": "2022-11-28",
               "User-Agent": "HybCapSeq-Quant-db-publisher"}
    if body is not None:
        headers["Content-Type"] = content_type
        headers["Content-Length"] = str(size if size is not None else len(body))
    connection = http.client.HTTPSConnection(host, timeout=600)
    try:
        connection.request(method, path, body=body, headers=headers)
        response = connection.getresponse()
        raw = response.read()
        data = json.loads(raw) if raw else None
        if response.status not in (200, 201):
            message = data.get("message", "request failed") if isinstance(data, dict) else "request failed"
            raise RuntimeError(f"GitHub HTTP {response.status}: {message}")
        return data
    finally:
        connection.close()


def api(method, path, data=None):
    return request(method, "api.github.com", path,
                   json.dumps(data).encode("utf-8") if data is not None else None)


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
    base = f"/repos/{repo}/releases"
    matches = [r for r in api("GET", base + "?per_page=100") if r["tag_name"] == TAG]
    if len(matches) > 1:
        raise ValueError("Multiple releases use the expected tag")
    current = matches[0] if matches else api("POST", base, {
        "tag_name": TAG, "target_commitish": os.environ["GITHUB_SHA"],
        "name": TITLE, "body": a.notes.read_text(), "draft": True})
    if current["name"] != TITLE:
        raise ValueError("Existing tag is used by a different release")
    detail = base + "/" + str(current["id"])
    if current["draft"] and current["target_commitish"] != os.environ["GITHUB_SHA"]:
        current = api("PATCH", detail, {"target_commitish": os.environ["GITHUB_SHA"]})
    for name in names:
        assets = {x["name"]: x for x in current["assets"]}
        existing = assets.get(name)
        if existing:
            if existing["size"] != expected[name]["size"]:
                raise ValueError(f"Existing asset differs; refusing to overwrite: {name}")
        else:
            if not current["draft"]:
                raise ValueError("Published release is missing an expected asset")
            upload = urlsplit(current["upload_url"].split("{", 1)[0])
            if upload.scheme != "https" or upload.hostname != "uploads.github.com":
                raise ValueError("Unexpected GitHub upload origin")
            if upload.path != detail + "/assets":
                raise ValueError("Unexpected GitHub upload path")
            with (a.assets / name).open("rb") as src:
                request("POST", upload.hostname, upload.path + "?" + urlencode({"name": name}),
                        src, "application/octet-stream", expected[name]["size"])
        for attempt in range(12):
            current = api("GET", detail)
            item = next(x for x in current["assets"] if x["name"] == name)
            if item.get("digest") is not None:
                break
            time.sleep(5)
        if any(item.get(k) != v for k, v in expected[name].items()):
            raise ValueError(f"Uploaded asset verification failed: {name}")
    if current["draft"]:
        current = api("PATCH", detail, {"draft": False})
    if current["draft"]:
        raise ValueError("Release was not published")
    print(json.dumps({"release_url": current["html_url"], "assets": [
        {k: x[k] for k in ["name", "size", "digest", "browser_download_url"]}
        for x in current["assets"]]}, indent=2))


if __name__ == "__main__":
    main()
