#!/usr/bin/env python3
"""Independent Arb replay of every fixed-radius angular tail direction."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import multiprocessing as mp
import time
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx
from report_validation import validate_planar_four, validate_fixed_context, require




def _strip_timing(payload):
    """Return a copy with wall-clock fields removed, recursively.

    Persisted artifacts are content-addressed in SHA256SUMS.txt, so their bytes
    must depend only on the mathematics.  Timings stay on stdout.
    """
    if isinstance(payload, dict):
        return {k: _strip_timing(v) for k, v in payload.items()
                if k != "elapsed_seconds"}
    if isinstance(payload, list):
        return [_strip_timing(v) for v in payload]
    return payload

def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_ball(value: F):
    return arb(fraction_text(value))


def determinant(u, v):
    return u[0] * v[1] - u[1] * v[0]


def load_witness(path: Path):
    raw = path.read_bytes()
    data = json.loads(raw.decode("utf-8"), parse_float=F, parse_int=F)
    A = [[F(value) for value in row] for row in data["A"]]
    B = [[F(value) for value in row] for row in data["B"]]
    c = [F(value) for value in data["c"]]
    validate_planar_four(A, B, c)
    columns = [(A[0][index], A[1][index]) for index in range(4)]
    minors = {
        pair: determinant(columns[pair[0]], columns[pair[1]]) * determinant(B[pair[0]], B[pair[1]])
        for pair in itertools.combinations(range(4), 2)
    }
    return B, c, minors, hashlib.sha256(raw).hexdigest()


def worker(arguments):
    witness, R, grid, start, end, precision = arguments
    ctx.prec = precision
    B, c, minors, _ = load_witness(Path(witness))
    positive = [pair for pair, value in minors.items() if value > 0]
    negative = [pair for pair, value in minors.items() if value < 0]
    if negative != [(2, 3)] or len(positive) != 5:
        return {"passed": False, "checked": 0, "failures": [{"reason": "minor sign pattern"}]}
    bad_pair = negative[0]
    bad_magnitude = -minors[bad_pair]
    K = {}
    L = {}
    for pair in positive:
        K[pair] = (16 * exact_ball(bad_magnitude) / exact_ball(minors[pair])).log()
        K[pair] += sum(exact_ball(abs(c[index])) for index in (*bad_pair, *pair))
        L[pair] = sum(
            (exact_ball(B[index][0]) ** 2 + exact_ball(B[index][1]) ** 2).sqrt()
            for index in (*bad_pair, *pair)
        )
    half_arc = arb.pi() / grid * (1 + arb(1) / (1 << 30))
    step = 2 * arb.pi() / grid
    worst_float = None
    worst_ball = None
    worst_index = None
    worst_pair = None
    certified_lower = None
    checked = 0
    failures = []

    def rho(index, cosine, sine):
        return abs(exact_ball(B[index][0]) * cosine + exact_ball(B[index][1]) * sine)

    for index in range(start, end):
        theta = step * index
        cosine = theta.cos()
        sine = theta.sin()
        negative_rates = rho(bad_pair[0], cosine, sine) + rho(bad_pair[1], cosine, sine)
        best = None
        best_float = None
        best_pair = None
        for pair in positive:
            gap = negative_rates - rho(pair[0], cosine, sine) - rho(pair[1], cosine, sine)
            gap -= L[pair] * half_arc
            # Extrapolation to r>=R requires the corrected arc slope, not
            # merely a positive radius-R margin when K can be negative.
            if not gap > 0:
                continue
            margin = R * gap - K[pair]
            lower_float = float(margin.lower())
            if best_float is None or lower_float > best_float:
                best = margin
                best_float = lower_float
                best_pair = pair
        checked += 1
        if best is None or not best > 0:
            failures.append({
                "grid_index": index,
                "best_margin_ball": None if best is None else best.str(40, radius=True),
            })
            break
        certified_lower = best.lower() if certified_lower is None else certified_lower.min(best.lower())
        if worst_float is None or best_float < worst_float:
            worst_float = best_float
            worst_ball = best.str(40, radius=True)
            worst_index = index
            worst_pair = list(best_pair)
    return {
        "passed": not failures and checked == end - start,
        "checked": checked,
        "worst_margin_lower": None if certified_lower is None else certified_lower.str(40, radius=True),
        "worst_grid_index": worst_index,
        "worst_pair": worst_pair,
        "failures": failures,
        "start": start,
        "end": end,
        "worst_float": worst_float,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--witness", default=str(root / "data" / "witness.json"))
    parser.add_argument("--point-report", default=str(root / "results" / "point_tail_recomputed.json"))
    parser.add_argument("--R", type=int, default=582)
    parser.add_argument("--grid", type=int, default=200000)
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--workers", type=int, default=max(1, min(8, mp.cpu_count())))
    parser.add_argument("--json-out", default=str(root / "results" / "point_tail_arb.json"))
    args = parser.parse_args()
    if args.precision < 128:
        raise SystemExit("precision must be at least 128 bits")
    witness_path = Path(args.witness).resolve()
    point_path = Path(args.point_report).resolve()
    _, _, _, witness_sha = load_witness(witness_path)
    validate_fixed_context(witness_sha, args.R, args.grid)
    require(args.workers > 0, 'workers must be positive')
    chunks = []
    for worker_index in range(args.workers):
        start = args.grid * worker_index // args.workers
        end = args.grid * (worker_index + 1) // args.workers
        if start < end:
            chunks.append((str(witness_path), args.R, args.grid, start, end, args.precision))
    started = time.time()
    context = mp.get_context("spawn")
    with context.Pool(len(chunks)) as pool:
        parts = pool.map(worker, chunks)
    checked = sum(part["checked"] for part in parts)
    failures = [failure for part in parts for failure in part["failures"]]
    passed = checked == args.grid and not failures and all(part["passed"] for part in parts)
    eligible = [part for part in parts if part["worst_float"] is not None]
    worst = min(eligible, key=lambda part: part["worst_float"]) if eligible else {}
    certified_lower = None
    ctx.prec = args.precision
    for part in eligible:
        bound = arb(part['worst_margin_lower']).lower()
        certified_lower = bound if certified_lower is None else certified_lower.min(bound)
    report = {
        "schema": "arb-point-tail-replay-v1",
        "verdict": "PASS" if passed else "FAIL",
        "arithmetic": "python-flint Arb independent angular-grid replay",
        "precision_bits": args.precision,
        "witness_sha256": witness_sha,
        "point_report_sha256": hashlib.sha256(point_path.read_bytes()).hexdigest(),
        "R": args.R,
        "grid": args.grid,
        "checked": checked,
        "worst_margin_lower": None if certified_lower is None else certified_lower.str(40, radius=True),
        "worst_grid_index": worst.get("worst_grid_index"),
        "worst_pair": worst.get("worst_pair"),
        "failures": failures[:10],
        "workers": len(chunks),
        "elapsed_seconds": time.time() - started,
    }
    text = json.dumps(report, indent=2) + "\n"
    Path(args.json_out).write_text(json.dumps(_strip_timing(report), indent=2) + "\n", encoding="utf-8", newline="\n")
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
