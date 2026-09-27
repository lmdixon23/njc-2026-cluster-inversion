#!/usr/bin/env python3
"""Exact regression for Proposition `prop:coorddep` in paper/main.tex.

This checks the mixed minors, the six attainable sector vertices, simplicity,
signed area, and the rational inequalities used in the analytic proof that the
chosen segment has negative averaged-Jacobian determinant.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from asymptotic_polygon import polygon_simplicity, polygon_signed_double_area

A = [(F(0), F(1)), (F(1), F(0)), (F(-10), F(-10))]  # columns
B = [(F(1), F(0)), (F(0), F(1)), (F(3,10), F(3,10))]  # rows


def det(u, v):
    return u[0]*v[1] - u[1]*v[0]


def add(*vs):
    return (sum(v[0] for v in vs), sum(v[1] for v in vs))


def main() -> int:
    mixed = {}
    for i, j in [(0,1),(0,2),(1,2)]:
        mixed[(i,j)] = det(A[i], A[j]) * det(B[i], B[j])
    assert mixed == {(0,1): F(-1), (0,2): F(3), (1,2): F(3)}

    # One exact integer direction from each open sector of the three row walls.
    directions = [(1,1),(-1,2),(-2,1),(-1,-1),(1,-2),(2,-1)]
    vertices = []
    for u in directions:
        active = [A[i] for i,b in enumerate(B) if b[0]*u[0] + b[1]*u[1] > 0]
        vertices.append(add(*active) if active else (F(0),F(0)))
    expected = [
        (F(-9),F(-9)), (F(-9),F(-10)), (F(1),F(0)),
        (F(0),F(0)), (F(0),F(1)), (F(-10),F(-9)),
    ]
    assert vertices == expected
    simple, failures = polygon_simplicity(vertices)
    assert simple and not failures
    assert polygon_signed_double_area(vertices) == F(38)

    # The manuscript uses e>2 to obtain these rational bounds:
    # d0>999/80000, d1>999/140000, d2<1/1152,
    # while d0<1/80 and d1<1/140.
    det_upper = -F(999,80000)*F(999,140000) + 3*F(1,1152)*(F(1,80)+F(1,140))
    assert det_upper == -F(182179,4800000000)
    assert det_upper < 0

    print('VERDICT: CLOSED-WIDTH COORDINATE-DEPENDENCE EXAMPLE VERIFIED')
    print('mixed minors =', mixed)
    print('vertices =', vertices)
    print('signed area = 19')
    print('averaged determinant rational upper bound =', det_upper)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
