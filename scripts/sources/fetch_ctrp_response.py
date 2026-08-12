"""Fetch the CTRPv2 drug-response data the pipeline trains on, pinned to one Zenodo version.

**What this downloads.** DrEval's reprocessed CTRPv2: the raw CTRPv2 dose-response measurements,
normalised per replicate against the no-drug control and re-fitted with CurveCurator, published as part
of the drevalpy benchmark suite.

    Ammar, J. et al. (daisybio). *DrEval* / ``drevalpy`` -- see
    ``papers/DrEval_s41467-026-72903-w.pdf``, Methods, "Benchmark data". Datasets on Zenodo,
    concept DOI 10.5281/zenodo.12633909.

**Why not ``drevalpy.datasets.loader.load_ctrpv2()``.** That function resolves the *concept* DOI to
whatever the latest release happens to be, and passes ``redownload=True``, so it re-downloads ~518 MB
on every call and silently returns different data as the upstream record is updated. A target whose
version cannot be named is not a citable source, so this script pins a concrete record instead and
verifies the bytes against the checksums that record publishes.

**Why the record also carries Cellosaurus.** ``meta.zip`` ships a Cellosaurus release, which
:mod:`scripts.sources.cellosaurus` uses to give every cell line a persistent accession.
Cellosaurus publishes no stable per-release download URL of its own, so pinning it by this record is
what makes that resolution reproducible.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from datetime import date
from pathlib import Path

import requests

from scripts.layout import FETCH_CACHE_DIRNAME, ZENODO_RESPONSE_RECORD

#: Which record, and where it lands, are configuration and live in ``layout``; the checksums are a
#: property of that record and live here, next to the code that verifies them.
ZENODO_RECORD = ZENODO_RESPONSE_RECORD
CACHE_DIRNAME = FETCH_CACHE_DIRNAME

#: Archives we need, with the MD5s the record publishes. ``CTRPv2.zip`` holds the response table;
#: ``meta.zip`` holds the Cellosaurus release and the tissue mapping. Bumping
#: ``ZENODO_RESPONSE_RECORD`` without updating these will (correctly) fail verification.
FILES: dict[str, str] = {
    "CTRPv2.zip": "a3a6cc49ed57ba1c9b163e56b88cce2f",
    "meta.zip": "658e75f9a0afa95b3c3c88f1f41e305f",
}

_API = "https://zenodo.org/api/records"


def _md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def fetch_ctrp_response(metadata_dir: str | Path, *, force: bool = False) -> Path:
    """Download, verify and extract the pinned record into ``metadata_dir/CACHE_DIRNAME``.

    Idempotent: an archive already present with the right MD5 is neither re-downloaded nor
    re-extracted, so this is safe to call at the head of the pipeline on every run. Returns the
    cache directory, and writes ``provenance.json`` into it recording the record, its DOI, its
    publication date, the retrieval date and the verified checksums -- so the version is readable
    from the data rather than only from this file.

    :raises RuntimeError: if a downloaded archive does not match its published MD5.
    """
    dest = Path(metadata_dir) / CACHE_DIRNAME
    dest.mkdir(parents=True, exist_ok=True)

    meta = requests.get(f"{_API}/{ZENODO_RECORD}", timeout=120)
    meta.raise_for_status()
    record = meta.json()
    published = record["metadata"]["publication_date"]
    print(f"Zenodo record {ZENODO_RECORD} -- {record['title']!r}, "
          f"DOI {record['doi']}, published {published}")

    for name, expected in FILES.items():
        archive = dest / name
        if archive.exists() and not force and _md5(archive) == expected:
            print(f"  {name}: cached, MD5 verified -- skipping download")
        else:
            url = f"{_API}/{ZENODO_RECORD}/files/{name}/content"
            print(f"  {name}: downloading from {url}")
            with requests.get(url, stream=True, timeout=600) as resp:
                resp.raise_for_status()
                with archive.open("wb") as fh:
                    for block in resp.iter_content(chunk_size=1 << 20):
                        fh.write(block)
            got = _md5(archive)
            if got != expected:
                archive.unlink()
                raise RuntimeError(
                    f"{name} failed verification: expected MD5 {expected}, got {got}. "
                    f"The record may have been revised; do not proceed on unverified data."
                )
            print(f"  {name}: downloaded, MD5 verified")

        with zipfile.ZipFile(archive) as z:
            members = [m for m in z.namelist() if not m.startswith("__MACOSX/")]
            if force or not all((dest / m).exists() for m in members):
                print(f"  {name}: extracting {len(members)} members")
                z.extractall(dest, members=members)

    (dest / "provenance.json").write_text(json.dumps({
        "zenodo_record": ZENODO_RECORD,
        "doi": record["doi"],
        "title": record["title"],
        "publication_date": published,
        "retrieved": date.today().isoformat(),
        "files": {name: {"md5": expected} for name, expected in FILES.items()},
        "fetched_by": "scripts/sources/fetch_ctrp_response.py",
    }, indent=2) + "\n")
    print(f"Cached at {dest}")
    return dest
