"""Beginner figure: geometry changes (shape kept) vs topology changes (shape changed)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(1)
n = 220
base = rng.normal(size=(n, 2))
th = rng.uniform(0, 2 * np.pi, n)
ring = np.c_[np.cos(th), np.sin(th)] * 1.6 + rng.normal(scale=0.12, size=(n, 2))
split = rng.normal(scale=0.45, size=(n, 2)) + np.c_[rng.choice([-1.4, 1.4], n), np.zeros(n)]
panels = [("Control", base, "the starting cloud"),
          ("Moved", base + [1.6, 0.8], "translation"),
          ("Grown", base * 1.45, "scaling"),
          ("Stretched", base * [1.9, 0.45], "anisotropy"),
          ("Ring", ring, "a hole appears (H1)"),
          ("Split", split, "two groups (H0)")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
fig, axes = plt.subplots(1, 6, figsize=(10, 2.3))
for ax, (t, P, sub) in zip(axes, panels):
    ax.scatter(P[:, 0], P[:, 1], s=4, c="black", linewidths=0)
    ax.set_xlim(-4.2, 4.2); ax.set_ylim(-4.2, 4.2); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#999999")
    ax.set_title(t, fontsize=9, fontweight="bold")
    ax.set_xlabel(sub, fontsize=7.5)
fig.text(0.215 + 0.09, 0.99, "Geometry changes: the shape is the same", ha="center", va="top", fontsize=8.5, style="italic")
fig.text(0.79, 0.99, "Topology changes: the shape is different", ha="center", va="top", fontsize=8.5, style="italic")
fig.add_artist(plt.Line2D([0.665, 0.665], [0.08, 0.9], color="black", linewidth=1))
fig.tight_layout(rect=(0, 0, 1, 0.9))
out = Path(__file__).resolve().parent / "outputs/concept_geometry_vs_topology.png"
fig.savefig(out, dpi=220, bbox_inches="tight")
print(out)
