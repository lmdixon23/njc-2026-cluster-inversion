#!/usr/bin/env python3
"""Independent Arb replay of every exported compact parameter leaf.

This implementation shares no interval helper with the mpmath generator.  It
reconstructs each dyadic box from its binary path, rebuilds parameter intervals
from the exact JSON witness, and requires a strictly positive Arb lower ball for
the determinant on every leaf.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import time
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx

_B = None
_C = None
_M = None
_R = None


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ball_interval(lo: F, hi: F):
    if lo > hi:
        raise ValueError("reversed interval")
    midpoint = (lo + hi) / 2
    radius = (hi - lo) / 2
    return arb(fraction_text(midpoint), fraction_text(radius))


def load_exact(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=F, parse_int=F)


def determinant(u, v):
    return u[0] * v[1] - u[1] * v[0]


def build_parameter_balls(data: dict, delta_a: F, delta_b: F, delta_c: F):
    A = [[ball_interval(F(value) - delta_a, F(value) + delta_a) for value in row] for row in data["A"]]
    B = [[ball_interval(F(value) - delta_b, F(value) + delta_b) for value in row] for row in data["B"]]
    C = [ball_interval(F(value) - delta_c, F(value) + delta_c) for value in data["c"]]
    columns = [(A[0][index], A[1][index]) for index in range(4)]
    minors = {
        (i, j): determinant(columns[i], columns[j]) * determinant(B[i], B[j])
        for i in range(4) for j in range(i + 1, 4)
    }
    return B, C, minors


def initialize_worker(witness: str, delta_a: str, delta_b: str, delta_c: str, R: int, precision: int):
    global _B, _C, _M, _R
    ctx.prec = precision
    data = load_exact(Path(witness))
    _B, _C, _M = build_parameter_balls(data, F(delta_a), F(delta_b), F(delta_c))
    _R = F(R)


def path_box(path: str):
    x1 = y1 = -_R
    x2 = y2 = _R
    for bit in path:
        width = x2 - x1
        height = y2 - y1
        if width >= height:
            midpoint = (x1 + x2) / 2
            if bit == "0":
                x2 = midpoint
            else:
                x1 = midpoint
        else:
            midpoint = (y1 + y2) / 2
            if bit == "0":
                y2 = midpoint
            else:
                y1 = midpoint
    return x1, y1, x2, y2


def sigma_prime_bounds(value):
    # Evaluating exp(value)+exp(-value) over a wide interval can lose the
    # positive lower endpoint through ball wrapping.  The stable scalar form
    # e/(1+e)^2 with e=exp(-|t|), combined with endpoint monotonicity, avoids
    # overflow and encloses the complete interval range.
    def scalar(endpoint):
        exponential = (-abs(endpoint)).exp()
        return exponential / (1 + exponential) ** 2
    left = scalar(value.lower())
    right = scalar(value.upper())
    lower = left.min(right).lower()
    upper = left.max(right)
    if value.contains(0):
        upper = upper.max(arb(1) / 4)
    return lower, upper.upper()


def check_path(path: str):
    x1, y1, x2, y2 = path_box(path)
    X = ball_interval(x1, x2)
    Y = ball_interval(y1, y2)
    weights = [sigma_prime_bounds(_B[i][0] * X + _B[i][1] * Y + _C[i]) for i in range(4)]
    value = arb(0)
    for (i, j), minor in _M.items():
        if minor > 0:
            value += minor.lower() * weights[i][0] * weights[j][0]
        elif minor < 0:
            value += minor.lower() * weights[i][1] * weights[j][1]
        else:
            return path, "minor sign ambiguity", False, float("-inf")
    return path, value.str(40, radius=True), bool(value > 0), float(value.lower())


def validate_partition(paths: list[str]):
    if not paths or any(set(path) - {"0", "1"} for path in paths):
        return False, False, 0
    ordered = sorted(paths)
    prefix_free = len(set(ordered)) == len(ordered) and all(
        not ordered[i + 1].startswith(ordered[i]) for i in range(len(ordered) - 1)
    )
    max_depth = max(map(len, ordered))
    numerator = sum(1 << (max_depth - len(path)) for path in ordered)
    volume_closed = numerator == 1 << max_depth
    return prefix_free, volume_closed, max_depth


def config_deltas(configuration: dict):
    common = F(str(configuration.get("delta", "0")))
    return (
        F(str(configuration.get("delta_A", common))),
        F(str(configuration.get("delta_B", common))),
        F(str(configuration.get("delta_c", common))),
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--chunksize", type=int, default=256)
    parser.add_argument("--json-out", default="")
    args = parser.parse_args()
    require_precision = args.precision >= 128
    if not require_precision:
        raise SystemExit("precision must be at least 128 bits")

    config_path = Path(args.config).resolve()
    configuration = json.loads(config_path.read_text(encoding="utf-8"))
    witness_path = (config_path.parent / configuration["witness"]).resolve()
    report_dir = root / "results" / configuration["name"]
    leaf_path = report_dir / "compact_leaf_paths.json"
    compact_path = report_dir / "compact.json"
    output_path = Path(args.json_out) if args.json_out else report_dir / "compact_arb.json"
    leaf_report = json.loads(leaf_path.read_text(encoding="utf-8"))
    paths = leaf_report["paths"]
    R = int(configuration["R"])
    delta_a, delta_b, delta_c = config_deltas(configuration)
    prefix_free, volume_closed, max_depth = validate_partition(paths)
    workers = args.workers or int(configuration.get("workers", max(1, min(8, mp.cpu_count()))))

    failures = []
    checked = 0
    minimum = None
    minimum_float = None
    started = time.time()
    context = mp.get_context("spawn")
    initargs = (
        str(witness_path), fraction_text(delta_a), fraction_text(delta_b),
        fraction_text(delta_c), R, args.precision,
    )
    with context.Pool(workers, initializer=initialize_worker, initargs=initargs) as pool:
        for path, interval, positive, lower_float in pool.imap_unordered(check_path, paths, chunksize=args.chunksize):
            checked += 1
            if minimum_float is None or lower_float < minimum_float:
                minimum_float = lower_float
                minimum = interval
            if not positive and len(failures) < 10:
                failures.append({"path": path, "determinant_ball": interval})

    passed = prefix_free and volume_closed and checked == len(paths) and not failures
    report = {
        "schema": "arb-compact-replay-v1",
        "verdict": "PASS" if passed else "FAIL",
        "arithmetic": "python-flint Arb independent ball replay",
        "precision_bits": args.precision,
        "witness_sha256": sha256(witness_path),
        "configuration_sha256": sha256(config_path),
        "compact_report_sha256": sha256(compact_path),
        "leaf_paths_sha256": sha256(leaf_path),
        "delta_A": fraction_text(delta_a),
        "delta_B": fraction_text(delta_b),
        "delta_c": fraction_text(delta_c),
        "R": R,
        "leaf_count": len(paths),
        "checked": checked,
        "prefix_free": prefix_free,
        "exact_volume_closure": volume_closed,
        "maximum_path_depth": max_depth,
        "minimum_leaf_lower": minimum,
        "failures": failures,
        "workers": workers,
        "elapsed_seconds": time.time() - started,
    }
    text = json.dumps(report, indent=2) + "\n"
    output_path.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
