#!/usr/bin/env python3
"""Exact geometry for planar saturating-ridge asymptotic polygons.

All decision-path arithmetic uses fractions.Fraction. Floating point is used only
for human-readable diagnostics and ranking candidate linear functionals.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cmp_to_key
from pathlib import Path
from typing import Iterable, Sequence

Vec = tuple[F, F]


def load_exact_json(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle, parse_float=F, parse_int=F)


def as_vec(v: Sequence[F]) -> Vec:
    if len(v) != 2:
        raise ValueError("Expected a planar vector")
    return F(v[0]), F(v[1])


def add(a: Vec, b: Vec) -> Vec:
    return a[0] + b[0], a[1] + b[1]


def sub(a: Vec, b: Vec) -> Vec:
    return a[0] - b[0], a[1] - b[1]


def scale(c: F, a: Vec) -> Vec:
    return c * a[0], c * a[1]


def dot(a: Vec, b: Vec) -> F:
    return a[0] * b[0] + a[1] * b[1]


def det(a: Vec, b: Vec) -> F:
    return a[0] * b[1] - a[1] * b[0]


def norm_inf(a: Vec) -> F:
    return max(abs(a[0]), abs(a[1]))


def sign(x: F) -> int:
    return 1 if x > 0 else -1 if x < 0 else 0


def r90(a: Vec) -> Vec:
    return -a[1], a[0]


def half_plane(v: Vec) -> int:
    x, y = v
    return 0 if y > 0 or (y == 0 and x >= 0) else 1


def angle_compare(left: tuple[Vec, int, int], right: tuple[Vec, int, int]) -> int:
    a, b = left[0], right[0]
    ha, hb = half_plane(a), half_plane(b)
    if ha != hb:
        return -1 if ha < hb else 1
    cross = det(a, b)
    if cross > 0:
        return -1
    if cross < 0:
        return 1
    # Collinear roots should only be antipodal copies of one row or a degeneracy.
    # Stable tie-breaking keeps output deterministic; callers separately reject
    # parallel distinct rows.
    return (left[1] > right[1]) - (left[1] < right[1])


@dataclass(frozen=True)
class Crossing:
    root: Vec
    row_index: int
    epsilon: int


@dataclass
class PolygonData:
    crossings: list[Crossing]
    edges: list[Vec]
    vertices: list[Vec]
    half_edges: list[Vec]
    prefixes: list[Vec]
    endpoint: Vec


def pairwise_nonparallel(rows: Sequence[Vec]) -> bool:
    return all(det(rows[i], rows[j]) != 0 for i in range(len(rows)) for j in range(i + 1, len(rows)))


def crossing_order(rows: Sequence[Vec]) -> list[Crossing]:
    if any(row == (0, 0) for row in rows):
        raise ValueError("All rows of B must be nonzero")
    if not pairwise_nonparallel(rows):
        raise ValueError("Rows of B must be pairwise nonparallel")

    raw: list[tuple[Vec, int, int]] = []
    for i, b in enumerate(rows):
        # Two roots of b dot u = 0.  The derivative along increasing angle is
        # b dot R90(u); its sign is the activation switch direction.
        base = (-b[1], b[0])
        for root in (base, scale(F(-1), base)):
            derivative = dot(b, r90(root))
            if derivative == 0:
                raise AssertionError("Impossible zero derivative at a nonzero row")
            raw.append((root, i, sign(derivative)))
    raw.sort(key=cmp_to_key(angle_compare))
    return [Crossing(root=r, row_index=i, epsilon=e) for r, i, e in raw]


def columns_from_matrix(A: Sequence[Sequence[F]]) -> list[Vec]:
    if len(A) != 2 or len(A[0]) != len(A[1]):
        raise ValueError("A must be a 2 by N matrix")
    return [(F(A[0][j]), F(A[1][j])) for j in range(len(A[0]))]


def build_polygon(A: Sequence[Sequence[F]], B: Sequence[Sequence[F]], saturation_span: F = F(1)) -> PolygonData:
    columns = columns_from_matrix(A)
    rows = [as_vec(row) for row in B]
    if len(columns) != len(rows):
        raise ValueError("A and B have incompatible hidden width")
    crossings = crossing_order(rows)
    edges = [scale(saturation_span * crossing.epsilon, columns[crossing.row_index]) for crossing in crossings]
    vertices: list[Vec] = [(F(0), F(0))]
    for edge in edges:
        vertices.append(add(vertices[-1], edge))
    if vertices[-1] != (0, 0):
        raise AssertionError("The full signed edge sequence must close")
    vertices = vertices[:-1]
    n = len(columns)
    half_edges = edges[:n]
    prefixes = [(F(0), F(0))]
    for edge in half_edges:
        prefixes.append(add(prefixes[-1], edge))
    endpoint = prefixes[-1]
    return PolygonData(crossings, edges, vertices, half_edges, prefixes, endpoint)


def polygon_signed_double_area(vertices: Sequence[Vec]) -> F:
    total = F(0)
    for i, p in enumerate(vertices):
        total += det(p, vertices[(i + 1) % len(vertices)])
    return total


def on_segment(point: Vec, a: Vec, b: Vec) -> bool:
    if det(sub(b, a), sub(point, a)) != 0:
        return False
    return (
        min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
        and min(a[1], b[1]) <= point[1] <= max(a[1], b[1])
    )


def segment_intersection_kind(a: Vec, b: Vec, c: Vec, d: Vec) -> str:
    """Return none, endpoint, proper, or overlap for two closed segments."""
    u, v = sub(b, a), sub(d, c)
    o1, o2 = det(u, sub(c, a)), det(u, sub(d, a))
    o3, o4 = det(v, sub(a, c)), det(v, sub(b, c))

    if o1 == o2 == o3 == o4 == 0:
        axis = 0 if max(abs(u[0]), abs(v[0])) >= max(abs(u[1]), abs(v[1])) else 1
        lo = max(min(a[axis], b[axis]), min(c[axis], d[axis]))
        hi = min(max(a[axis], b[axis]), max(c[axis], d[axis]))
        if hi < lo:
            return "none"
        if hi == lo:
            return "endpoint"
        return "overlap"

    if o1 * o2 < 0 and o3 * o4 < 0:
        return "proper"
    if (
        (o1 == 0 and on_segment(c, a, b))
        or (o2 == 0 and on_segment(d, a, b))
        or (o3 == 0 and on_segment(a, c, d))
        or (o4 == 0 and on_segment(b, c, d))
    ):
        return "endpoint"
    return "none"


def polygon_simplicity(vertices: Sequence[Vec]) -> tuple[bool, list[dict]]:
    m = len(vertices)
    failures: list[dict] = []
    if m < 3:
        return False, [{"reason": "fewer than three vertices"}]
    if len(set(vertices)) != m:
        failures.append({"reason": "repeated vertex"})
    for i in range(m):
        if vertices[i] == vertices[(i + 1) % m]:
            failures.append({"reason": "zero edge", "edge": i})

    for i in range(m):
        a, b = vertices[i], vertices[(i + 1) % m]
        for j in range(i + 1, m):
            adjacent = j == i + 1 or (i == 0 and j == m - 1)
            c, d = vertices[j], vertices[(j + 1) % m]
            kind = segment_intersection_kind(a, b, c, d)
            if adjacent:
                if kind in {"proper", "overlap", "none"}:
                    failures.append({"reason": "bad adjacent-edge incidence", "edges": [i, j], "kind": kind})
            elif kind != "none":
                failures.append({"reason": "nonadjacent-edge intersection", "edges": [i, j], "kind": kind})
    return not failures, failures


def mixed_minors(A: Sequence[Sequence[F]], B: Sequence[Sequence[F]]) -> dict[str, F]:
    columns = columns_from_matrix(A)
    rows = [as_vec(row) for row in B]
    return {
        f"{i}{j}": det(columns[i], columns[j]) * det(rows[i], rows[j])
        for i in range(len(columns)) for j in range(i + 1, len(columns))
    }


def strict_convexity_by_turns(edges: Sequence[Vec]) -> tuple[bool, int]:
    turns = [sign(det(edges[i], edges[(i + 1) % len(edges)])) for i in range(len(edges))]
    nonzero = [s for s in turns if s]
    if len(nonzero) != len(turns):
        return False, 0
    orientation = nonzero[0]
    return all(s == orientation for s in nonzero), orientation


def best_integer_functional(edges: Sequence[Vec], limit: int = 60) -> Vec | None:
    best: tuple[float, int, Vec] | None = None
    for qx in range(-limit, limit + 1):
        for qy in range(-limit, limit + 1):
            if qx == 0 and qy == 0:
                continue
            common = math.gcd(abs(qx), abs(qy))
            q = (F(qx // common), F(qy // common))
            vals = [dot(q, edge) for edge in edges]
            if min(vals) <= 0:
                continue
            score = float(min(vals)) / math.hypot(qx, qy)
            complexity = abs(qx) + abs(qy)
            candidate = (score, -complexity, q)
            if best is None or candidate > best:
                best = candidate
    return None if best is None else best[2]


def monotone_one_sided(prefixes: Sequence[Vec], q: Vec) -> dict:
    edges = [sub(prefixes[k], prefixes[k - 1]) for k in range(1, len(prefixes))]
    q_values = [dot(q, edge) for edge in edges]
    S = prefixes[-1]
    side_values = [det(S, prefixes[k]) for k in range(1, len(prefixes) - 1)]
    side_signs = {sign(v) for v in side_values}
    passed = min(q_values) > 0 and len(side_signs) == 1 and 0 not in side_signs and S != (0, 0)
    return {
        "passed": passed,
        "q": q,
        "q_values": q_values,
        "side_values": side_values,
        "side_sign": next(iter(side_signs)) if len(side_signs) == 1 else 0,
    }


def piecewise_linear_value(xs: Sequence[F], ys: Sequence[F], x: F) -> F:
    if not (xs[0] <= x <= xs[-1]):
        raise ValueError("Point lies outside graph interval")
    if x == xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            if x == xs[i]:
                return ys[i]
            width = xs[i + 1] - xs[i]
            return ys[i] + (x - xs[i]) * (ys[i + 1] - ys[i]) / width
    raise AssertionError("Unreachable interpolation state")


def graph_gap(prefixes: Sequence[Vec], q: Vec, r: Vec | None = None) -> dict:
    if r is None:
        r = r90(q)
    if det(q, r) == 0:
        raise ValueError("q and r must be independent")
    xs = [dot(q, point) for point in prefixes]
    ys = [dot(r, point) for point in prefixes]
    if any(xs[i + 1] <= xs[i] for i in range(len(xs) - 1)):
        return {"passed": False, "reason": "chain is not strictly q-monotone", "q": q, "r": r}
    T, Y = xs[-1], ys[-1]
    breakpoints = sorted({x for x in xs[1:-1]} | {T - x for x in xs[1:-1]})
    values: list[tuple[F, F]] = []
    for beta in breakpoints:
        if 0 < beta < T:
            h = piecewise_linear_value(xs, ys, beta) + piecewise_linear_value(xs, ys, T - beta) - Y
            values.append((beta, h))
    signs = {sign(value) for _, value in values}
    passed = bool(values) and len(signs) == 1 and 0 not in signs
    return {"passed": passed, "q": q, "r": r, "breakpoint_values": values, "sign": next(iter(signs)) if len(signs) == 1 else 0}


def uniform_geometry_radius(prefixes: Sequence[Vec], q: Vec) -> dict:
    """Conservative L-infinity radius for perturbing every half-edge generator.

    If ||e_k||_inf <= delta for every k and delta is below the returned lower
    bound, the same q-monotonicity and one-sided determinant signs survive.
    """
    n = len(prefixes) - 1
    edges = [sub(prefixes[k], prefixes[k - 1]) for k in range(1, len(prefixes))]
    S = prefixes[-1]
    qnorm1 = abs(q[0]) + abs(q[1])
    q_cap = min(dot(q, edge) / qnorm1 for edge in edges)
    side = [det(S, prefixes[k]) for k in range(1, n)]
    if not side or len({sign(v) for v in side}) != 1 or 0 in {sign(v) for v in side}:
        raise ValueError("One-sided criterion does not hold")

    def side_ok(delta: F) -> bool:
        for k, value in enumerate(side, start=1):
            bound = 2 * n * delta * norm_inf(prefixes[k]) + 2 * k * delta * norm_inf(S) + 2 * n * k * delta * delta
            if bound >= abs(value):
                return False
        return True

    lo, hi = F(0), q_cap
    for _ in range(100):
        mid = (lo + hi) / 2
        if side_ok(mid):
            lo = mid
        else:
            hi = mid
    return {"delta_edge_inf_lower": lo, "q_cap": q_cap, "side_cap_lower": lo}


def uniform_row_radius(rows: Sequence[Vec]) -> F:
    """Conservative L-infinity row perturbation radius preserving line order."""
    if not pairwise_nonparallel(rows):
        raise ValueError("Rows are already degenerate")
    upper = min(norm_inf(row) for row in rows)

    def ok(delta: F) -> bool:
        if delta >= upper:
            return False
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                bound = 2 * delta * norm_inf(rows[j]) + 2 * delta * norm_inf(rows[i]) + 2 * delta * delta
                if bound >= abs(det(rows[i], rows[j])):
                    return False
        return True

    lo, hi = F(0), upper
    for _ in range(100):
        mid = (lo + hi) / 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def fraction_text(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def jsonable(value):
    if isinstance(value, F):
        return fraction_text(value)
    if isinstance(value, tuple):
        return [jsonable(item) for item in value]
    if isinstance(value, list):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    return value
