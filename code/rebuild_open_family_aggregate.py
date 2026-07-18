#!/usr/bin/env python3
"""Rebuild a semantically validated, content-addressed open-family aggregate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from report_validation import ReportValidationError, build_combined


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("report_dir")
    ap.add_argument("--config", default="")
    ns = ap.parse_args()
    report_dir = Path(ns.report_dir)
    existing = report_dir / "combined.json"
    if not ns.config:
        raise SystemExit("supply --config so the aggregate is bound to the config bytes")
    config_path = Path(ns.config).resolve()
    try:
        combined = build_combined(report_dir.resolve(), config_path)
    except ReportValidationError as exc:
        raise SystemExit(f"semantic report validation failed: {exc}") from exc
    existing.write_text(json.dumps(combined, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(combined, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
