#!/usr/bin/env python3
"""Exact counterexample: determinant signs do not determine polygon simplicity."""
from fractions import Fraction as F

from asymptotic_polygon import add, det, polygon_simplicity, scale


def vertices_from_half(edges):
    full = list(edges) + [scale(F(-1), edge) for edge in edges]
    vertices = [(F(0), F(0))]
    for edge in full:
        vertices.append(add(vertices[-1], edge))
    assert vertices[-1] == (0, 0)
    return vertices[:-1]


def sign_pattern(edges):
    return tuple(1 if det(edges[i], edges[j]) > 0 else -1 if det(edges[i], edges[j]) < 0 else 0
                 for i in range(len(edges)) for j in range(i + 1, len(edges)))


def main() -> int:
    directions = [(F(1), F(-3)), (F(1), F(1)), (F(1), F(0)), (F(1), F(3))]
    simple_edges = directions
    crossed_edges = [directions[0], scale(F(7), directions[1]), scale(F(7), directions[2]), directions[3]]
    simple, fail_simple = polygon_simplicity(vertices_from_half(simple_edges))
    crossed, fail_crossed = polygon_simplicity(vertices_from_half(crossed_edges))
    same_signs = sign_pattern(simple_edges) == sign_pattern(crossed_edges)
    passed = simple and not crossed and same_signs
    print(f"[{'PASS' if simple else 'FAIL'}] unit-length realization is simple")
    print(f"[{'PASS' if not crossed else 'FAIL'}] rescaled realization self-intersects")
    print(f"[{'PASS' if same_signs else 'FAIL'}] pairwise determinant sign pattern is unchanged")
    if fail_crossed:
        for item in fail_crossed:
            if item.get("reason") == "nonadjacent-edge intersection":
                print(f"crossing edges {item['edges']} kind={item['kind']}")
    print(f"VERDICT: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
