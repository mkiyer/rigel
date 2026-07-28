"""Figure: the RNA strand overdispersion, before and after the per-junction SJ strand table.

Two panels, four real cfRNA libraries:

  A — the FITTED od_r that reaches the strand likelihood.  OLD sits on the Beta(2,2) ceiling on
      4/4 libraries (it carries no information about any of them); NEW lands inside the
      0.001–0.016 band measured independently from deep junctions before this design existed.
  B — the RAW pooled MoM, before the ceiling clips it.  OLD is 10.7–79.9 — mathematically
      impossible for an intraclass correlation, which is bounded by 1 — so the ceiling was the
      only thing making the output finite.  NEW is admissible on every library, unclipped.

    OMP_NUM_THREADS=1 python scratchpad/sj_table_figure.py
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("docs/calibration/figures/sj_strand_table_od_r.png")

# Validated categorical slots 1-3 (all-pairs, light mode) from the data-viz reference palette.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#b8b7b0"
SURFACE = "#fcfcfb"

# Measured 2026-07-28 on ~/Downloads/rigel_runs/cfrna/_calib_cache (scratchpad/sj_table_od_r.py).
LIB = ["LBX0190", "LBX0588", "MO_3021", "vcap"]
KAPPA = [0.002314, 0.012032, 0.002034, 0.000060]
N_JUNC = [9_526, 4_506, 67_935, 155_025]
OLD_FIT = [0.2000, 0.2000, 0.2000, 0.2000]  # all four clipped to the Beta(2,2) ceiling
NEW_FIT = [0.00858, 0.01368, 0.00736, 0.01341]
OLD_RAW = [10.728, 13.547, 79.882, 12.614]
NEW_RAW = [0.00857, 0.01364, 0.00736, 0.01341]
# The INDEPENDENT route: deep junctions only, re-centred on their own mean
# (strand_overdispersion_design.md §3, measured before this design existed).
DEEP = [0.0035, 0.0020, 0.0023, 0.0011]

BAND_LO, BAND_HI = 0.001, 0.016
CEILING = 0.2

fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(13.0, 5.6), facecolor=SURFACE, gridspec_kw={"wspace": 0.30}
)
fig.subplots_adjust(top=0.80, bottom=0.30, left=0.135, right=0.965)
y = np.arange(len(LIB))[::-1]


def style(ax, xlabel):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{n}\n" + r"$\kappa$=" + f"{k:.5f}  ·  {j:,} junc" for n, k, j in zip(LIB, KAPPA, N_JUNC)],
        fontsize=8.5,
        color=INK2,
    )
    ax.tick_params(axis="x", labelsize=8.5, colors=INK2, length=3)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(xlabel, fontsize=9, color=INK2, labelpad=8)
    ax.grid(axis="x", color=MUTED, lw=0.5, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.set_ylim(-0.85, len(LIB) - 0.25)


# ---------------------------------------------------------------- A: the fitted od_r
axA.axvspan(BAND_LO, BAND_HI, color=AQUA, alpha=0.13, zorder=0, lw=0)
axA.axvline(CEILING, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=1)
for yi, o, n, d in zip(y, OLD_FIT, NEW_FIT, DEEP):
    axA.plot([n, o], [yi, yi], color=MUTED, lw=2, zorder=2, solid_capstyle="round")
    axA.plot(o, yi, "o", ms=9, color=ORANGE, mec=SURFACE, mew=2, zorder=4)
    axA.plot(n, yi, "o", ms=9, color=BLUE, mec=SURFACE, mew=2, zorder=4)
    axA.plot(d, yi, "D", ms=7, color=AQUA, mec=SURFACE, mew=1.6, zorder=3)
    axA.annotate(
        f"{n:.4f}", (n, yi), textcoords="offset points", xytext=(0, 11),
        ha="center", fontsize=8, color=INK,
    )
style(axA, "fitted RNA strand overdispersion  $od_r$  (log scale)")
axA.set_xlim(4e-4, 0.6)
axA.text(
    CEILING, len(LIB) - 0.42, " Beta(2,2)\n ceiling",
    fontsize=8, color=INK2, va="top", ha="left",
)
axA.text(
    np.sqrt(BAND_LO * BAND_HI), -0.72, "independent deep-junction range 0.001–0.016",
    fontsize=8, color=INK2, ha="center",
)
axA.set_title(
    "A   $od_r$ was pinned to its guard on 4/4 libraries; now it is a measurement",
    fontsize=10.5, color=INK, loc="left", pad=12,
)

# ---------------------------------------------------------------- B: the raw MoM
axB.axvspan(1.0, 300, color=ORANGE, alpha=0.10, zorder=0, lw=0)
axB.axvline(1.0, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=1)
for yi, o, n in zip(y, OLD_RAW, NEW_RAW):
    axB.plot([n, o], [yi, yi], color=MUTED, lw=2, zorder=2, solid_capstyle="round")
    axB.plot(o, yi, "o", ms=9, color=ORANGE, mec=SURFACE, mew=2, zorder=4)
    axB.plot(n, yi, "o", ms=9, color=BLUE, mec=SURFACE, mew=2, zorder=4)
    axB.annotate(
        f"{o:.1f}", (o, yi), textcoords="offset points", xytext=(0, 11),
        ha="center", fontsize=8, color=INK,
    )
style(axB, "raw pooled method-of-moments, before clipping  (log scale)")
axB.set_xlim(4e-3, 300)
axB.text(
    1.15, -0.72, "an intraclass correlation cannot exceed 1",
    fontsize=8, color=INK2, ha="left",
)
axB.set_title(
    "B   the old estimate was not a dispersion at all",
    fontsize=10.5, color=INK, loc="left", pad=12,
)

handles = [
    plt.Line2D([], [], marker="o", ls="", ms=9, color=ORANGE, mec=SURFACE, mew=2,
               label="OLD — accumulator boundary spliced sides (ANNOT + UNANNOT + IMPLICIT)"),
    plt.Line2D([], [], marker="o", ls="", ms=9, color=BLUE, mec=SURFACE, mew=2,
               label="NEW — per-junction SJ strand table (the population $\\kappa$ is the marginal of)"),
    plt.Line2D([], [], marker="D", ls="", ms=7, color=AQUA, mec=SURFACE, mew=1.6,
               label="independent route — deep junctions (depth $\\geq$ 1000), re-centred on their own mean"),
]
fig.legend(
    handles=handles, loc="lower left", ncol=1, frameon=False, fontsize=8.5,
    bbox_to_anchor=(0.135, 0.015), labelcolor=INK2, handletextpad=0.6,
    borderaxespad=0.0,
)
fig.suptitle(
    "The RNA strand Beta-Binomial: both halves now come from one population",
    fontsize=13, color=INK, x=0.012, ha="left", y=0.955, fontweight="medium",
)
fig.text(
    0.012, 0.895,
    "four real cfRNA libraries · fitted 2026-07-28 · docs/calibration/sj_strand_table_design.md",
    fontsize=8.5, color=INK2, ha="left",
)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=170, facecolor=SURFACE)
print(f"wrote {OUT}")
