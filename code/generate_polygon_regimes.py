#!/usr/bin/env python3
"""Generate the publication polygon-regime figure from canonical exact data."""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "pdf.compression": 9,
})

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]


def exact_witness_vertices(witness_path: Path):
    data = json.loads(witness_path.read_text(encoding="utf-8"), parse_float=F, parse_int=F)
    columns = [(F(data["A"][0][i]), F(data["A"][1][i])) for i in range(4)]
    active_sets = [(), (1,), (1, 2), (1, 2, 3), (0, 1, 2, 3), (0, 2, 3), (0, 3), (0,)]
    return [
        (sum((columns[i][0] for i in active), F(0)), sum((columns[i][1] for i in active), F(0)))
        for active in active_sets
    ]


def obstruction_vertices():
    half_edges = [(F(2), F(3)), (F(-4), F(-4)), (F(5), F(-3)), (F(-1), F(2))]
    edges = half_edges + [(-x, -y) for x, y in half_edges]
    vertices = [(F(0), F(0))]
    for edge in edges:
        vertices.append((vertices[-1][0] + edge[0], vertices[-1][1] + edge[1]))
    return vertices[:-1]


def exact_crossing(a, b, c, d):
    ux, uy = b[0] - a[0], b[1] - a[1]
    vx, vy = d[0] - c[0], d[1] - c[1]
    denominator = ux * vy - uy * vx
    if denominator == 0:
        return None
    wx, wy = c[0] - a[0], c[1] - a[1]
    t = (wx * vy - wy * vx) / denominator
    s = (wx * uy - wy * ux) / denominator
    if 0 < t < 1 and 0 < s < 1:
        return a[0] + t * ux, a[1] + t * uy
    return None


def transverse_crossings(vertices):
    points = []
    count = len(vertices)
    for i in range(count):
        for j in range(i + 1, count):
            if j == i + 1 or (i == 0 and j == count - 1):
                continue
            point = exact_crossing(vertices[i], vertices[(i + 1) % count], vertices[j], vertices[(j + 1) % count])
            if point is not None and point not in points:
                points.append(point)
    return points


def to_array(vertices) -> np.ndarray:
    return np.array([[float(x), float(y)] for x, y in vertices], dtype=float)


def close_polygon(poly: np.ndarray) -> np.ndarray:
    return np.vstack([poly, poly[0]])


def center_and_align(poly: np.ndarray) -> np.ndarray:
    center = (poly[0] + poly[len(poly) // 2]) / 2
    centered = poly - center
    axis = poly[len(poly) // 2] - poly[0]
    angle = np.arctan2(axis[1], axis[0])
    cosine = np.cos(-angle)
    sine = np.sin(-angle)
    rotation = np.array([[cosine, -sine], [sine, cosine]])
    return centered @ rotation.T


def winding_number(poly: np.ndarray, point: np.ndarray) -> int:
    total_angle = 0.0
    closed = close_polygon(poly)
    for start, end in zip(closed[:-1], closed[1:]):
        u = start - point
        v = end - point
        total_angle += np.arctan2(u[0] * v[1] - u[1] * v[0], u @ v)
    return int(np.rint(total_angle / (2 * np.pi)))


def shade_positive_winding(ax, poly, padding, resolution=180):
    lower = poly.min(axis=0) - padding
    upper = poly.max(axis=0) + padding
    xs = np.linspace(lower[0], upper[0], resolution)
    ys = np.linspace(lower[1], upper[1], resolution)
    values = np.array([[winding_number(poly, np.array([x, y])) for x in xs] for y in ys])
    positive = sorted(value for value in set(values.ravel()) if value > 0)
    if positive:
        levels = np.arange(0.5, max(positive) + 1.5, 1.0)
        ax.contourf(xs, ys, values, levels=levels, alpha=0.15)
    return lower, upper


def add_orientation_arrows(ax, poly):
    closed = close_polygon(poly)
    for edge_index in (1, len(poly) // 2 + 1):
        start = closed[edge_index]
        end = closed[edge_index + 1]
        ax.add_patch(FancyArrowPatch(
            start + 0.34 * (end - start), start + 0.68 * (end - start),
            arrowstyle="-|>", mutation_scale=10, linewidth=1.1,
        ))


def plot_panel(ax, poly, title, labels, padding, crossings=None):
    lower, upper = shade_positive_winding(ax, poly, padding)
    closed = close_polygon(poly)
    ax.plot(closed[:, 0], closed[:, 1], marker="o", linewidth=2)
    add_orientation_arrows(ax, poly)
    if crossings is not None and len(crossings):
        points = np.asarray(crossings, dtype=float)
        ax.scatter(points[:, 0], points[:, 1], marker="x", s=38, linewidths=1.6,
                   color="#9b1c31", zorder=6)
    for position, label in labels:
        ax.text(position[0], position[1], label, ha="center", va="center", fontsize=9.5,
                bbox={"boxstyle": "circle,pad=0.16", "alpha": 0.62})
    ax.set_xlim(lower[0], upper[0])
    ax.set_ylim(lower[1], upper[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(1)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--witness", default=str(ROOT / "data" / "witness.json"))
    parser.add_argument("--pdf-out", default=str(ROOT / "paper" / "figures" / "polygon_regimes.pdf"))
    parser.add_argument("--png-out", default="")
    args = parser.parse_args()
    witness_path = Path(args.witness).resolve()
    pdf_path = Path(args.pdf_out).resolve()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    convex = to_array([(0, 0), (4, 3), (8, 4), (12, 3), (16, 0), (12, -3), (8, -4), (4, -3)])
    witness = center_and_align(to_array(exact_witness_vertices(witness_path)))
    obstruction_exact = obstruction_vertices()
    obstruction = to_array(obstruction_exact)
    crossings = to_array(transverse_crossings(obstruction_exact))

    figure, axes = plt.subplots(1, 3, figsize=(12.2, 4.05), constrained_layout=True)
    plot_panel(axes[0], convex, "(a) convex, index one",
               [(np.array([8.0, 0.0]), "$1$"), (np.array([8.0, 4.8]), "$0$")], 0.75)
    plot_panel(axes[1], witness, "(b) nonconvex, index one",
               [(np.array([0.0, 0.0]), "$1$"), (np.array([0.0, 1.1]), "$0$")], 0.30)
    plot_panel(axes[2], obstruction, "(c) self-intersecting, index two",
               [(np.array([0.12, 0.72]), "$2$"), (np.array([1.0, -1.15]), "$1$"),
                (np.array([3.65, 2.55]), "$0$")], 0.55, crossings=crossings)

    pdf_metadata = {
        "Title": "Three asymptotic polygon winding regimes",
        "Author": "Logan M. Dixon",
        "Creator": "code/generate_polygon_regimes.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(pdf_path, bbox_inches="tight", metadata=pdf_metadata)
    if args.png_out:
        png_path = Path(args.png_out).resolve()
        png_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(png_path, dpi=220, bbox_inches="tight", metadata={"Software": "code/generate_polygon_regimes.py"})
        print(f"Wrote: {png_path}")
    plt.close(figure)
    print(f"witness_sha256 = {hashlib.sha256(witness_path.read_bytes()).hexdigest()}")
    print(f"Wrote: {pdf_path}")


if __name__ == "__main__":
    main()
