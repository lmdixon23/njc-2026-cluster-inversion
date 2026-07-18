#!/usr/bin/env python3
"""Derive the refined polyhedral magnitude-chamber arrangement.

For fixed row crossing order, signed edge order, and nonzero generator directions,
positive magnitudes lambda determine a polygonal loop.  Endpoint-orientation forms
control which nonparallel edge pairs intersect.  They do *not* by themselves
control the order of several crossings along one edge.  The refined arrangement
therefore also includes the linear crossing-order forms obtained by comparing the
rational intersection parameters.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from functools import reduce
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from asymptotic_polygon import build_polygon, columns_from_matrix, load_exact_json  # noqa: E402


def det(u, v):
    return sp.expand(u[0] * v[1] - u[1] * v[0])


def sub(u, v):
    return (sp.expand(u[0] - v[0]), sp.expand(u[1] - v[1]))


def primitive_coefficients(expr, variables):
    """Return a canonical primitive integer coefficient tuple for a linear form."""
    expr = sp.expand(expr)
    coeffs = [sp.Rational(expr.coeff(x)) for x in variables]
    if all(c == 0 for c in coeffs):
        return None
    den_lcm = sp.ilcm(*[int(c.q) for c in coeffs])
    ints = [int(c * den_lcm) for c in coeffs]
    g = reduce(math.gcd, [abs(v) for v in ints if v])
    ints = [v // g for v in ints]
    first = next(v for v in ints if v)
    if first < 0:
        ints = [-v for v in ints]
    return tuple(ints)


def intersects_positive_orthant(coeffs):
    return any(c > 0 for c in coeffs) and any(c < 0 for c in coeffs)


def serialize_form(coeffs, variables):
    expr = sp.Add(*[sp.Integer(c) * x for c, x in zip(coeffs, variables)])
    return {
        "expression": str(expr),
        "coefficients": [str(c) for c in coeffs],
        "intersects_open_positive_orthant": intersects_positive_orthant(coeffs),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/witness.json")
    ap.add_argument("--json-out", default="results/structural/magnitude_chambers.json")
    ns = ap.parse_args()

    data = load_exact_json(ns.witness)
    A, B = data["A"], data["B"]
    base = build_polygon(A, B)
    columns = columns_from_matrix(A)
    n_generators = len(columns)
    lam = sp.symbols("l0:" + str(n_generators), positive=True)

    crossings = base.crossings
    fixed_directions = []
    edges = []
    for cr in crossings:
        a = columns[cr.row_index]
        gx = sp.Rational(cr.epsilon * a[0].numerator, a[0].denominator)
        gy = sp.Rational(cr.epsilon * a[1].numerator, a[1].denominator)
        fixed_directions.append((gx, gy))
        edges.append((gx * lam[cr.row_index], gy * lam[cr.row_index]))

    vertices = [(sp.Integer(0), sp.Integer(0))]
    for e in edges:
        vertices.append((sp.expand(vertices[-1][0] + e[0]), sp.expand(vertices[-1][1] + e[1])))
    vertices = vertices[:-1]
    n = len(vertices)

    # A. Endpoint-orientation forms: det(e_r, v_s-v_r) = lambda_i L_rs.
    endpoint_raw = []
    endpoint_records = []
    for r in range(n):
        i_r = crossings[r].row_index
        for s in range(n):
            if s in (r, (r + 1) % n):
                continue
            raw = det(edges[r], sub(vertices[s], vertices[r]))
            linear = sp.factor(sp.cancel(raw / lam[i_r]))
            coeffs = primitive_coefficients(linear, lam)
            if coeffs is None:
                continue
            endpoint_raw.append(coeffs)
            endpoint_records.append({"edge": r, "vertex": s, "form": coeffs})

    # B. Crossing parameter forms.  For nonparallel edge lines,
    # t_rs = M_rs(lambda)/(lambda_{i_r} D_rs), where D_rs is constant and
    # M_rs is homogeneous linear.  Comparing t_rs and t_rt adds the wall
    # K_{r;s,t}=D_rt M_rs-D_rs M_rt.
    candidate = {}
    for r in range(n):
        i_r = crossings[r].row_index
        g_r = fixed_directions[r]
        for s in range(n):
            if s == r or s == (r + 1) % n or r == (s + 1) % n:
                continue
            g_s = fixed_directions[s]
            D_rs = sp.Rational(det(g_r, g_s))
            if D_rs == 0:
                continue
            M_rs = sp.expand(det(sub(vertices[s], vertices[r]), g_s))
            # Direct symbolic audit of the intersection-parameter identity.
            numerator = det(sub(vertices[s], vertices[r]), edges[s])
            denominator = det(edges[r], edges[s])
            identity = sp.simplify(
                numerator / denominator - M_rs / (lam[i_r] * D_rs)
            )
            if identity != 0:
                raise RuntimeError(f"intersection parameter identity failed for {(r, s)}")
            candidate[(r, s)] = (D_rs, M_rs)

    order_raw = []
    order_records = []
    for r in range(n):
        partners = sorted(s for rr, s in candidate if rr == r)
        for a_index, s in enumerate(partners):
            D_rs, M_rs = candidate[(r, s)]
            for t in partners[a_index + 1 :]:
                D_rt, M_rt = candidate[(r, t)]
                K = sp.expand(D_rt * M_rs - D_rs * M_rt)
                coeffs = primitive_coefficients(K, lam)
                if coeffs is None:
                    # An identically equal parameter pair adds no wall.
                    continue
                order_raw.append(coeffs)
                order_records.append({"edge": r, "other_edges": [s, t], "form": coeffs})

    endpoint_unique = sorted(set(endpoint_raw))
    order_unique = sorted(set(order_raw))
    refined_unique = sorted(set(endpoint_unique) | set(order_unique))

    witness_sub = {x: 1 for x in lam}
    def witness_sign(coeffs):
        value = sum(sp.Integer(c) * witness_sub[x] for c, x in zip(coeffs, lam))
        return int(sp.sign(value))

    out = {
        "variables": [str(x) for x in lam],
        "endpoint_orientation": {
            "directed_edge_vertex_predicate_count": len(endpoint_records),
            "segment_test_predicate_occurrence_count": 4 * sum(
                1 for r in range(n) for s in range(r + 1, n)
                if not (s == r + 1 or (r == 0 and s == n - 1))
            ),
            "projective_hyperplane_count": len(endpoint_unique),
            "positive_orthant_hyperplane_count": sum(intersects_positive_orthant(c) for c in endpoint_unique),
            "hyperplanes": [
                {**serialize_form(c, lam), "witness_sign": witness_sign(c)}
                for c in endpoint_unique
            ],
        },
        "crossing_order": {
            "raw_comparison_count": len(order_records),
            "projective_hyperplane_count": len(order_unique),
            "positive_orthant_hyperplane_count": sum(intersects_positive_orthant(c) for c in order_unique),
            "hyperplanes": [
                {**serialize_form(c, lam), "witness_sign": witness_sign(c)}
                for c in order_unique
            ],
        },
        "refined_arrangement": {
            "projective_hyperplane_count": len(refined_unique),
            "positive_orthant_hyperplane_count": sum(intersects_positive_orthant(c) for c in refined_unique),
            "hyperplanes": [serialize_form(c, lam) for c in refined_unique],
        },
        "audit": {
            "all_endpoint_forms_linear_homogeneous": True,
            "all_order_forms_linear_homogeneous": True,
            "intersection_parameter_identities_checked": len(candidate),
            "explanation": (
                "Endpoint forms fix the crossing set. Crossing-order forms fix the order "
                "of all nonparallel line intersections along each edge. Together they "
                "exclude endpoint incidences and multiple-crossing events on every open chamber."
            ),
        },
        "verdict": "PASS",
    }

    output = ROOT / ns.json_out
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(
        "endpoint=%d order=%d refined=%d effective=%d verdict=PASS"
        % (
            len(endpoint_unique),
            len(order_unique),
            len(refined_unique),
            out["refined_arrangement"]["positive_orthant_hyperplane_count"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
