#!/usr/bin/env python3
"""Exact certificate for the witness asymptotic-polygon theorem package."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F

from asymptotic_polygon import (
    as_vec,
    best_integer_functional,
    build_polygon,
    columns_from_matrix,
    fraction_text,
    graph_gap,
    jsonable,
    mixed_minors,
    monotone_one_sided,
    polygon_signed_double_area,
    polygon_simplicity,
    strict_convexity_by_turns,
    uniform_geometry_radius,
    uniform_row_radius,
    load_exact_json,
)

EXPECTED_SEQUENCE = [(2, 1), (3, 1), (0, 1), (1, -1), (2, -1), (3, -1), (0, -1), (1, 1)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--witness", default=str(Path(__file__).resolve().parents[1] / "data" / "witness.json"))
    parser.add_argument("--json-out", default="")
    args = parser.parse_args()

    path = Path(args.witness)
    raw = path.read_bytes()
    data = load_exact_json(path)
    A, B = data["A"], data["B"]
    poly = build_polygon(A, B)

    sequence = [(cross.row_index, cross.epsilon) for cross in poly.crossings]
    sequence_ok = sequence == EXPECTED_SEQUENCE
    simple, failures = polygon_simplicity(poly.vertices)
    double_area = polygon_signed_double_area(poly.vertices)
    convex, orientation = strict_convexity_by_turns(poly.edges)
    minors = mixed_minors(A, B)
    minor_signs = {key: 1 if value > 0 else -1 if value < 0 else 0 for key, value in minors.items()}

    q = best_integer_functional(poly.half_edges, limit=80)
    if q is None:
        raise RuntimeError("No small integer monotonicity functional found")
    one_sided = monotone_one_sided(poly.prefixes, q)
    gap = graph_gap(poly.prefixes, q)
    geom_radius = uniform_geometry_radius(poly.prefixes, q) if one_sided["passed"] else None
    row_radius = uniform_row_radius([as_vec(row) for row in B])

    checks = {
        "crossing_sequence": sequence_ok,
        "simple_exact_segment_formula": simple,
        "monotone_one_sided": one_sided["passed"],
        "monotone_graph_gap": gap["passed"],
        "nonconvex": not convex,
        "one_negative_mixed_minor": sum(value < 0 for value in minors.values()) == 1,
    }
    passed = all(checks.values())

    report = {
        "verdict": "PASS" if passed else "FAIL",
        "witness_sha256": hashlib.sha256(raw).hexdigest(),
        "checks": checks,
        "crossing_sequence": sequence,
        "vertices": poly.vertices,
        "signed_double_area": double_area,
        "signed_area": double_area / 2,
        "mixed_minors": minors,
        "mixed_minor_signs": minor_signs,
        "strictly_convex": convex,
        "turn_orientation": orientation,
        "monotone_one_sided": one_sided,
        "graph_gap": gap,
        "geometry_open_radius": geom_radius,
        "row_order_open_radius_inf": row_radius,
        "segment_failures": failures,
    }

    print(f"witness_sha256 = {report['witness_sha256']}")
    print("crossing_sequence = " + ", ".join(("+" if eps > 0 else "-") + f"a{idx}" for idx, eps in sequence))
    for name, ok in checks.items():
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    print(f"q = ({fraction_text(q[0])}, {fraction_text(q[1])})")
    print("q(g_k) = " + ", ".join(fraction_text(v) for v in one_sided["q_values"]))
    print("det(S,p_k) = " + ", ".join(fraction_text(v) for v in one_sided["side_values"]))
    print("graph_gap = " + ", ".join(f"x={float(x):.9g}:H={float(h):.9g}" for x, h in gap["breakpoint_values"]))
    if geom_radius:
        print(f"delta_A_edge_inf >= {float(geom_radius['delta_edge_inf_lower']):.12g}")
    print(f"delta_B_row_inf >= {float(row_radius):.12g}")
    print(f"VERDICT: {report['verdict']}")

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(jsonable(report), indent=2) + "\n", encoding="utf-8")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
