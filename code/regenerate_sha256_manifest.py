#!/usr/bin/env python3
"""Regenerate SHA256SUMS.txt from the tracked public repository tree."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SHA256SUMS.txt"


def main() -> None:
    raw = subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=ROOT
    ).decode("utf-8")
    paths = sorted(path for path in raw.split("\0") if path and path != "SHA256SUMS.txt")
    lines = []
    for rel in paths:
        if "\\" in rel or rel.startswith("/") or ".." in Path(rel).parts:
            raise SystemExit(f"unsafe tracked path: {rel}")
        path = ROOT / rel
        if not path.is_file():
            raise SystemExit(f"tracked path is not a file: {rel}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {rel}")
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {MANIFEST.relative_to(ROOT)} with {len(lines)} files")


if __name__ == "__main__":
    main()
