#!/usr/bin/env python3
"""Verify the repository canonical-file SHA-256 manifest."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SHA256SUMS.txt"
LINE_RE = re.compile(r"^([0-9a-f]{64})  (.+)$")


def fail(message: str) -> None:
    print(f"VERDICT: SHA-256 MANIFEST FAILED: {message}", file=sys.stderr)
    raise SystemExit(1)


def tracked_inventory() -> set[str]:
    raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT)
    return {
        item.decode("utf-8")
        for item in raw.split(b"\0")
        if item and item.decode("utf-8") != "SHA256SUMS.txt"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=str(MANIFEST))
    args = parser.parse_args()
    manifest = Path(args.manifest)
    if not manifest.is_file():
        fail(f"manifest is missing: {manifest}")

    seen: set[str] = set()
    checked = 0
    for line_number, raw in enumerate(
        manifest.read_text(encoding="utf-8").splitlines(), 1
    ):
        if not raw or raw.startswith("#"):
            continue
        match = LINE_RE.fullmatch(raw)
        if match is None:
            fail(f"malformed line {line_number}")
        expected, rel_text = match.groups()
        if "\\" in rel_text:
            fail(f"nonportable path on line {line_number}: {rel_text}")
        rel = Path(rel_text)
        if rel.is_absolute() or ".." in rel.parts or rel_text.startswith("/"):
            fail(f"unsafe path on line {line_number}: {rel_text}")
        normalized = rel.as_posix()
        if normalized in seen:
            fail(f"duplicate path: {normalized}")
        seen.add(normalized)
        path = ROOT / rel
        if not path.is_file():
            fail(f"missing file: {normalized}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            fail(f"hash mismatch: {normalized}")
        checked += 1

    if checked == 0:
        fail("manifest contains no files")
    expected_paths = tracked_inventory()
    if seen != expected_paths:
        missing = sorted(expected_paths - seen)
        extra = sorted(seen - expected_paths)
        detail = []
        if missing:
            detail.append("missing tracked paths: " + ", ".join(missing))
        if extra:
            detail.append("untracked paths: " + ", ".join(extra))
        fail("; ".join(detail))
    print("VERDICT: SHA-256 MANIFEST VERIFIED")
    print(f"files = {checked}")


if __name__ == "__main__":
    main()
