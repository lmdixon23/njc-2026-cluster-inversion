#!/usr/bin/env python3
"""Semantic validation and binding for public certificate reports.

The checksum manifest protects bytes.  This module separately checks that each
proof report says what its top-level verdict claims, matches the selected
configuration and witness, covers the expected domain, and is bound to every
data export on which it depends.
"""
from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path

EXPECTED_POSITIVE = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3]]
EXPECTED_NEGATIVE = [[2, 3]]
EXPECTED_REPORTS = {
    "compact",
    "compact_arb",
    "compact_independent",
    "compact_leaf_paths",
    "geometry",
    "tail",
    "tail_independent",
}
DATA_EXPORTS = {"compact_leaf_paths"}
SPLIT_RULE = "longer-side; x on ties; child0=lower/left, child1=upper/right"
NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


class ReportValidationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReportValidationError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportValidationError(f"cannot read {path}: {exc}") from exc
    require(isinstance(value, dict), f"{path}: top level must be an object")
    return value


def as_fraction(value, label: str) -> F:
    try:
        return F(Decimal(str(value))) if "/" not in str(value) else F(str(value))
    except Exception as exc:
        raise ReportValidationError(f"{label}: invalid rational value {value}") from exc


def first_decimal(value, label: str) -> Decimal:
    match = NUMBER_RE.search(str(value))
    require(match is not None, f"{label}: no numeric endpoint")
    return Decimal(match.group(0))


def config_deltas(configuration: dict) -> tuple[F, F, F]:
    common = configuration.get("delta", "0")
    return (
        as_fraction(configuration.get("delta_A", common), "config delta_A"),
        as_fraction(configuration.get("delta_B", common), "config delta_B"),
        as_fraction(configuration.get("delta_c", common), "config delta_c"),
    )


def abc_payload_sha256(witness_path: Path) -> str:
    data = load_json(witness_path)
    payload = json.dumps(
        {
            "A": [[str(x) for x in row] for row in data["A"]],
            "B": [[str(x) for x in row] for row in data["B"]],
            "c": [str(x) for x in data["c"]],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_signs(report: dict, label: str) -> None:
    signs = report.get("minor_signs")
    require(isinstance(signs, dict), f"{label}: missing minor_signs")
    require(signs.get("positive") == EXPECTED_POSITIVE, f"{label}: wrong positive pairs")
    require(signs.get("negative") == EXPECTED_NEGATIVE, f"{label}: wrong negative pair")
    require(signs.get("ambiguous") == [], f"{label}: ambiguous minor sign")


def validate_leaf_paths(report: dict, expected_R: int) -> dict:
    require(report.get("R") == str(expected_R), "compact_leaf_paths: wrong R")
    require(report.get("split_rule") == SPLIT_RULE, "compact_leaf_paths: wrong split rule")
    paths = report.get("paths")
    require(isinstance(paths, list) and paths, "compact_leaf_paths: empty path list")
    require(all(isinstance(p, str) and set(p) <= {"0", "1"} for p in paths),
            "compact_leaf_paths: malformed path")
    ordered = sorted(paths)
    require(len(set(ordered)) == len(ordered), "compact_leaf_paths: duplicate path")
    require(all(not ordered[i + 1].startswith(ordered[i]) for i in range(len(ordered) - 1)),
            "compact_leaf_paths: paths are not prefix free")
    max_depth = max(map(len, ordered))
    numerator = sum(1 << (max_depth - len(path)) for path in ordered)
    denominator = 1 << max_depth
    require(numerator == denominator, "compact_leaf_paths: exact dyadic volume does not close")
    return {
        "count": len(ordered),
        "max_depth": max_depth,
        "normalized_numerator": numerator,
        "normalized_denominator": denominator,
    }


def validate_compact(report: dict, configuration: dict, witness_sha: str, path_info: dict) -> None:
    label = "compact"
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, f"{label}: witness hash mismatch")
    require(report.get("R") == int(configuration["R"]), f"{label}: wrong R")
    require(report.get("grid") == int(configuration["grid"]), f"{label}: wrong grid")
    require(as_fraction(report.get("hmin"), f"{label} hmin") == as_fraction(configuration["hmin"], "config hmin"),
            f"{label}: wrong hmin")
    expected_deltas = config_deltas(configuration)
    actual_deltas = tuple(as_fraction(report.get(key), f"{label} {key}") for key in ("delta_A", "delta_B", "delta_c"))
    require(actual_deltas == expected_deltas, f"{label}: parameter deltas do not match config")
    validate_signs(report, label)
    interior = report.get("interior")
    require(isinstance(interior, dict), f"{label}: missing interior result")
    require(interior.get("passed") is True, f"{label}: interior passed is false")
    require(interior.get("paused") is False, f"{label}: interior is paused")
    require(interior.get("unresolved") == 0, f"{label}: unresolved leaves remain")
    require(interior.get("budget_exhausted") is False, f"{label}: budget exhausted")
    require(interior.get("first_failures") == [], f"{label}: recorded failures are nonempty")
    processed = interior.get("processed")
    leaves = interior.get("certified_leaves")
    require(isinstance(processed, int) and isinstance(leaves, int) and leaves > 0,
            f"{label}: invalid processed or leaf count")
    require(processed == 2 * leaves - 1, f"{label}: full binary tree count is inconsistent")
    require(leaves == path_info["count"], f"{label}: leaf export count mismatch")


def validate_compact_independent(report: dict, configuration: dict, path_info: dict) -> None:
    label = "compact_independent"
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    expected_deltas = config_deltas(configuration)
    actual_deltas = tuple(as_fraction(report.get(key), f"{label} {key}") for key in ("delta_A", "delta_B", "delta_c"))
    require(actual_deltas == expected_deltas, f"{label}: parameter deltas do not match config")
    require(report.get("leaf_count") == path_info["count"], f"{label}: leaf count mismatch")
    require(report.get("checked") == path_info["count"], f"{label}: checked count mismatch")
    require(report.get("prefix_free") is True, f"{label}: prefix-free check failed")
    require(report.get("exact_volume_closure") is True, f"{label}: volume closure failed")
    require(report.get("failures") == [], f"{label}: failures are nonempty")
    volume = report.get("volume", {})
    require(str(volume.get("normalized_numerator")) == str(volume.get("normalized_denominator")),
            f"{label}: reported volume does not close")


def validate_geometry(report: dict, witness_sha: str) -> None:
    label = "geometry"
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, f"{label}: witness hash mismatch")
    checks = report.get("checks")
    require(isinstance(checks, dict) and checks and all(value is True for value in checks.values()),
            f"{label}: a required predicate failed")
    require(report.get("segment_failures") == [], f"{label}: segment failures are nonempty")


def validate_tail(report: dict, configuration: dict, witness_sha: str, payload_sha: str) -> None:
    label = "tail"
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, f"{label}: witness hash mismatch")
    require(report.get("ABC_payload_sha256") == payload_sha, f"{label}: payload hash mismatch")
    require(report.get("R") == int(configuration["R"]), f"{label}: wrong R")
    require(report.get("grid") == int(configuration["grid"]), f"{label}: wrong grid")
    expected_deltas = config_deltas(configuration)
    actual_deltas = tuple(as_fraction(report.get(key), f"{label} {key}") for key in ("delta_A", "delta_B", "delta_c"))
    require(actual_deltas == expected_deltas, f"{label}: parameter deltas do not match config")
    require(report.get("minor_sign_preservation") is True, f"{label}: minor signs are not preserved")
    require(first_decimal(report.get("tail_margin_residual_lower"), f"{label} residual") > 0,
            f"{label}: residual margin is not positive")


def validate_tail_independent(report: dict, configuration: dict) -> None:
    label = "tail_independent"
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    expected_deltas = config_deltas(configuration)
    actual_deltas = tuple(as_fraction(report.get(key), f"{label} {key}") for key in ("delta_A", "delta_B", "delta_c"))
    require(actual_deltas == expected_deltas, f"{label}: parameter deltas do not match config")
    require(report.get("minor_sign_preservation") is True, f"{label}: minor signs are not preserved")
    require(first_decimal(report.get("residual_lower"), f"{label} residual") > 0,
            f"{label}: residual margin is not positive")


def validate_compact_arb(
    report: dict,
    configuration: dict,
    config_path: Path,
    witness_sha: str,
    compact_path: Path,
    leaf_path: Path,
    path_info: dict,
) -> None:
    label = "compact_arb"
    require(report.get("schema") == "arb-compact-replay-v1", f"{label}: wrong schema")
    require(report.get("verdict") == "PASS", f"{label}: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, f"{label}: witness hash mismatch")
    require(report.get("configuration_sha256") == file_sha256(config_path), f"{label}: config hash mismatch")
    require(report.get("compact_report_sha256") == file_sha256(compact_path), f"{label}: compact hash mismatch")
    require(report.get("leaf_paths_sha256") == file_sha256(leaf_path), f"{label}: leaf export hash mismatch")
    require(report.get("R") == int(configuration["R"]), f"{label}: wrong R")
    expected_deltas = config_deltas(configuration)
    actual_deltas = tuple(as_fraction(report.get(key), f"{label} {key}") for key in ("delta_A", "delta_B", "delta_c"))
    require(actual_deltas == expected_deltas, f"{label}: parameter deltas do not match config")
    require(isinstance(report.get("precision_bits"), int) and report["precision_bits"] >= 128,
            f"{label}: precision is too low")
    require(report.get("leaf_count") == path_info["count"], f"{label}: leaf count mismatch")
    require(report.get("checked") == path_info["count"], f"{label}: checked count mismatch")
    require(report.get("prefix_free") is True, f"{label}: prefix-free check failed")
    require(report.get("exact_volume_closure") is True, f"{label}: volume closure failed")
    require(report.get("failures") == [], f"{label}: failures are nonempty")
    require(first_decimal(report.get("minimum_leaf_lower"), f"{label} minimum") > 0,
            f"{label}: minimum leaf lower bound is not positive")


def resolve_witness(config_path: Path, configuration: dict) -> Path:
    witness = Path(configuration["witness"])
    return witness if witness.is_absolute() else (config_path.parent / witness).resolve()


def validate_open_family_reports(report_dir: Path, config_path: Path) -> tuple[dict, dict]:
    configuration = load_json(config_path)
    actual_names = {path.stem for path in report_dir.glob("*.json") if path.name != "combined.json"}
    require(actual_names == EXPECTED_REPORTS,
            f"{report_dir}: report set mismatch; missing={sorted(EXPECTED_REPORTS-actual_names)} extra={sorted(actual_names-EXPECTED_REPORTS)}")
    paths = {name: report_dir / f"{name}.json" for name in EXPECTED_REPORTS}
    reports = {name: load_json(path) for name, path in paths.items()}
    witness_path = resolve_witness(config_path, configuration)
    witness_sha = file_sha256(witness_path)
    payload_sha = abc_payload_sha256(witness_path)
    path_info = validate_leaf_paths(reports["compact_leaf_paths"], int(configuration["R"]))
    validate_compact(reports["compact"], configuration, witness_sha, path_info)
    validate_compact_independent(reports["compact_independent"], configuration, path_info)
    validate_geometry(reports["geometry"], witness_sha)
    validate_tail(reports["tail"], configuration, witness_sha, payload_sha)
    validate_tail_independent(reports["tail_independent"], configuration)
    validate_compact_arb(
        reports["compact_arb"], configuration, config_path, witness_sha,
        paths["compact"], paths["compact_leaf_paths"], path_info,
    )
    return configuration, path_info


def build_combined(report_dir: Path, config_path: Path) -> dict:
    configuration, path_info = validate_open_family_reports(report_dir, config_path)
    reports = {}
    for name in sorted(EXPECTED_REPORTS):
        path = report_dir / f"{name}.json"
        obj = load_json(path)
        reports[name] = {
            "role": "data_export" if name in DATA_EXPORTS else "check_report",
            "verdict": None if name in DATA_EXPORTS else obj.get("verdict"),
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
        }
    return {
        "schema": "open-family-combined-v3",
        "verdict": "PASS",
        "configuration": configuration,
        "configuration_sha256": file_sha256(config_path),
        "leaf_count": path_info["count"],
        "reports": reports,
    }


def validate_combined(report_dir: Path, config_path: Path) -> None:
    actual = load_json(report_dir / "combined.json")
    expected = build_combined(report_dir, config_path)
    require(actual == expected, f"{report_dir / 'combined.json'}: aggregate content is stale or inconsistent")


def validate_point_tail(report: dict, witness_sha: str, R: int, grid: int) -> None:
    require(report.get("verdict") == "PASS", "point_tail_recomputed: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, "point_tail_recomputed: witness hash mismatch")
    require(report.get("R") == R and report.get("grid") == grid, "point_tail_recomputed: wrong grid context")
    validate_signs(report, "point_tail_recomputed")
    tail = report.get("tail", {})
    require(tail.get("passed") is True, "point_tail_recomputed: tail failed")
    require(first_decimal(tail.get("worst_margin"), "point tail margin") > 0,
            "point_tail_recomputed: margin is not positive")


def validate_point_tail_arb(report_path: Path, point_path: Path, witness_path: Path, R: int, grid: int) -> None:
    report = load_json(report_path)
    witness_sha = file_sha256(witness_path)
    require(report.get("schema") == "arb-point-tail-replay-v1", "point_tail_arb: wrong schema")
    require(report.get("verdict") == "PASS", "point_tail_arb: verdict is not PASS")
    require(report.get("witness_sha256") == witness_sha, "point_tail_arb: witness hash mismatch")
    require(report.get("point_report_sha256") == file_sha256(point_path), "point_tail_arb: point report hash mismatch")
    require(report.get("R") == R and report.get("grid") == grid, "point_tail_arb: wrong grid context")
    require(isinstance(report.get("precision_bits"), int) and report["precision_bits"] >= 128,
            "point_tail_arb: precision is too low")
    require(report.get("checked") == grid, "point_tail_arb: incomplete grid coverage")
    require(report.get("failures") == [], "point_tail_arb: failures are nonempty")
    require(first_decimal(report.get("worst_margin_lower"), "point_tail_arb margin") > 0,
            "point_tail_arb: margin is not positive")
