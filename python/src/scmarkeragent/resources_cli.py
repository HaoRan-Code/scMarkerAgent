"""Fetch and verify the curated marker resource.

The resource is hundreds of megabytes and lives in a versioned archive rather than in the
package, so getting it is a step a user has to take. It is also the thing an annotation is
measured against, which is why every file is checked against a shipped SHA-256 index
before the directory is declared usable: a truncated download is otherwise indistinguishable
from a smaller database, and the run that follows would quietly annotate against half a
resource.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

from .config import PACKAGE_ROOT

RESOURCE_INDEX = PACKAGE_ROOT / "resources" / "resource_index.json"
CHUNK = 8 << 20


def load_index(path: Path = RESOURCE_INDEX) -> dict:
    if not path.is_file():
        raise SystemExit(
            f"resource index not found: {path}\n"
            "  This package was built without required_files/resource_index.json."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(dest: Path, index: dict, quiet: bool = False) -> list[str]:
    """Every indexed file present, the right size and the right hash. Returns problems."""
    problems: list[str] = []
    for name, expected in sorted(index["files"].items()):
        path = dest / name
        if not path.is_file():
            problems.append(f"missing: {name}")
            continue
        size = path.stat().st_size
        if size != expected["bytes"]:
            problems.append(
                f"wrong size: {name} ({size} bytes, expected {expected['bytes']})"
            )
            continue
        if not quiet:
            print(f"  checking {name} ...", flush=True)
        if sha256(path) != expected["sha256"]:
            problems.append(f"checksum mismatch: {name}")
    return problems


def _download(url: str, target: Path) -> None:
    print(f"downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response:  # noqa: S310 - user-supplied archive URL
        total = int(response.headers.get("Content-Length") or 0)
        done = 0
        with target.open("wb") as handle:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                handle.write(chunk)
                done += len(chunk)
                if total:
                    print(f"\r  {done / total:6.1%}", end="", flush=True)
    print()


def _unpack(archive: Path, dest: Path, index: dict) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as staging_name:
        staging = Path(staging_name)
        with tarfile.open(archive) as tar:
            # Refuse absolute paths and parent traversal rather than trusting the archive.
            for member in tar.getmembers():
                target = (staging / member.name).resolve()
                if not str(target).startswith(str(staging.resolve())):
                    raise SystemExit(f"archive member escapes the destination: {member.name}")
            tar.extractall(staging)  # noqa: S202 - members validated above
        # The bundle may be archived with or without a top-level directory.
        roots = [staging] + [p for p in staging.iterdir() if p.is_dir()]
        first = next(iter(index["files"]))
        source = next((root for root in roots if (root / first).is_file()), None)
        if source is None:
            raise SystemExit(
                f"archive does not contain the expected bundle (looked for {first})"
            )
        for name in index["files"]:
            target = dest / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source / name), str(target))


def download_resources(dest: Path, url: str | None, index_path: Path = RESOURCE_INDEX) -> int:
    index = load_index(index_path)
    dest = Path(dest).expanduser()
    human = index["total_bytes"] / 1e9

    if not verify(dest, index, quiet=True):
        print(f"resource already complete at {dest} ({index['bundle']})")
        return 0

    if not url:
        url = index.get("archive_url") or ""
    if not url:
        raise SystemExit(
            "no archive URL is configured yet.\n"
            f"  Bundle {index['bundle']} is ~{human:.1f} GB across "
            f"{index['file_count']} files.\n"
            "  Pass --url <archive>, or download the bundle by hand and unpack it into\n"
            f"  {dest}\n"
            "  See required_files/README.md for the archive DOI."
        )

    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "scmarkeragent-resources.tar.gz"
        _download(url, archive)
        print("unpacking ...", flush=True)
        _unpack(archive, dest, index)

    print("verifying ...", flush=True)
    problems = verify(dest, index)
    if problems:
        for problem in problems:
            print(f"  {problem}")
        raise SystemExit("resource verification failed; the download is not usable")
    print(f"resource ready at {dest} ({index['bundle']}, {human:.1f} GB)")
    print(f"  use it with: --resource-dir {dest}")
    return 0


def verify_resources(dest: Path, index_path: Path = RESOURCE_INDEX) -> int:
    index = load_index(index_path)
    dest = Path(dest).expanduser()
    problems = verify(dest, index)
    if problems:
        for problem in problems:
            print(f"  {problem}")
        print(f"resource at {dest} is NOT usable ({len(problems)} problem(s))")
        return 2
    print(f"resource at {dest} verified ({index['bundle']}, {index['file_count']} files)")
    return 0
