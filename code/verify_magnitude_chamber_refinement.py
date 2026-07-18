#!/usr/bin/env python3
"""Exact regression test for the refined magnitude-chamber theorem.

The two rational magnitude vectors below have identical signs on every endpoint
orientation predicate used by the original proof and the same crossing set, but
different winding spectra.  Crossing-order comparisons separate them, confirming
why the refined linear walls are necessary.
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "code"))
from asymptotic_polygon import det, jsonable, on_segment, sub  # noqa: E402
from polygon_winding import line_intersection, spectrum  # noqa: E402

DIRECTIONS = [(Q(6), Q(1)), (Q(-4), Q(4)), (Q(-3), Q(-5)), (Q(2), Q(4)), (Q(5), Q(-2))]
LAMBDA_A = [Q(49, 16), Q(35, 16), Q(15, 4), Q(4), Q(43, 16)]
LAMBDA_B = [Q(9, 4), Q(7, 8), Q(2), Q(39, 16), Q(11, 8)]


def polygon(magnitudes):
    half = [(magnitudes[i] * DIRECTIONS[i][0], magnitudes[i] * DIRECTIONS[i][1]) for i in range(5)]
    edges = half + [(-x, -y) for x, y in half]
    vertices = [(Q(0), Q(0))]
    for e in edges:
        vertices.append((vertices[-1][0] + e[0], vertices[-1][1] + e[1]))
    if vertices[-1] != vertices[0]:
        raise RuntimeError("polygon does not close")
    return vertices[:-1]


def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0


def endpoint_occurrence_signs(vertices):
    """Directed edge versus every nonincident vertex orientation sign."""
    n = len(vertices)
    out = []
    for r in range(n):
        er = sub(vertices[(r + 1) % n], vertices[r])
        for s in range(n):
            if s in (r, (r + 1) % n):
                continue
            value = det(er, sub(vertices[s], vertices[r]))
            if value == 0:
                raise RuntimeError(f"nongeneric endpoint predicate at edge {r}, vertex {s}")
            out.append(sign(value))
    return tuple(out)


def crossing_pairs(vertices):
    n = len(vertices)
    out = []
    for r, s in combinations(range(n), 2):
        if s == r + 1 or (r == 0 and s == n - 1):
            continue
        p = line_intersection(vertices[r], vertices[(r + 1) % n], vertices[s], vertices[(s + 1) % n])
        if p is not None:
            out.append((r, s))
    return tuple(out)


def line_parameter(vertices, r, s):
    """Parameter t on the oriented line of edge r at its intersection with edge s."""
    er = sub(vertices[(r + 1) % len(vertices)], vertices[r])
    es = sub(vertices[(s + 1) % len(vertices)], vertices[s])
    denominator = det(er, es)
    if denominator == 0:
        return None
    return det(sub(vertices[s], vertices[r]), es) / denominator


def crossing_order_signs(vertices):
    n = len(vertices)
    values = []
    labels = []
    for r in range(n):
        partners = []
        for s in range(n):
            if s == r or s == (r + 1) % n or r == (s + 1) % n:
                continue
            t = line_parameter(vertices, r, s)
            if t is not None:
                partners.append((s, t))
        for (s, t_rs), (u, t_ru) in combinations(partners, 2):
            difference = t_rs - t_ru
            if difference == 0:
                # Some equalities are structural and produce the zero form; they do
                # not define a chamber wall and are omitted.
                continue
            values.append(sign(difference))
            labels.append((r, s, u))
    return tuple(values), tuple(labels)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-out", default="results/structural/magnitude_chamber_refinement_witness.json")
    ns = ap.parse_args()

    va, vb = polygon(LAMBDA_A), polygon(LAMBDA_B)
    old_a, old_b = endpoint_occurrence_signs(va), endpoint_occurrence_signs(vb)
    order_a, labels_a = crossing_order_signs(va)
    order_b, labels_b = crossing_order_signs(vb)
    if labels_a != labels_b:
        raise RuntimeError("comparison indexing changed")

    spec_a = spectrum(va)["winding_values"]
    spec_b = spectrum(vb)["winding_values"]
    pairs_a, pairs_b = crossing_pairs(va), crossing_pairs(vb)
    changed_order = [list(labels_a[i]) for i, (x, y) in enumerate(zip(order_a, order_b)) if x != y]

    checks = {
        "endpoint_sign_vectors_identical": old_a == old_b,
        "endpoint_predicates_nonzero": 0 not in old_a and 0 not in old_b,
        "endpoint_predicate_occurrences": len(old_a) == 80,
        "crossing_sets_identical": pairs_a == pairs_b,
        "crossing_order_refinement_separates": bool(changed_order),
        "spectrum_A": spec_a == [-2, -1, 0, 1],
        "spectrum_B": spec_b == [-1, 0, 1],
    }
    verdict = "PASS" if all(checks.values()) else "FAIL"
    out = {
        "description": "Counterexample to endpoint-orientation signs as a complete winding-spectrum chamber invariant",
        "directions": DIRECTIONS,
        "lambda_A": LAMBDA_A,
        "lambda_B": LAMBDA_B,
        "endpoint_predicate_occurrences": len(old_a),
        "crossing_pairs": pairs_a,
        "spectrum_A": spec_a,
        "spectrum_B": spec_b,
        "changed_crossing_order_comparisons": changed_order,
        "checks": checks,
        "verdict": verdict,
    }
    output = ROOT / ns.json_out
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(jsonable(out), indent=2) + "\n", encoding="utf-8")
    print(json.dumps(jsonable(out), indent=2))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
