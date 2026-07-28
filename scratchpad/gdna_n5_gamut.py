"""N5 — THE FULL GAMUT: every condition's landscape, fit vs reference, on a log-y axis.

Owner request: see all conditions, and specifically that the zero-gDNA and LOW-gDNA cases behave. Those are
where a landscape can invent enrichment that is not there, and they are the hardest regime because the
depleted mode is at the resolution wall and there is barely any gDNA to place.

Each panel: the FIT (regions-only substrate, kNN 0.5) against the CORRECTED reference (population
resolution) and, in thin orange, the PRODUCTION reference that N5 replaced — so the comb and its repair are
both visible everywhere, not just on the three conditions the summary figure shows.

Panel titles are colour-coded by gDNA level: ZERO in red, LOW (gdna1/gdna5) in orange, the rest in black.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n5_gamut.py [--suite ambig|quick]
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n5_roughness import region_substrate, roughness  # noqa: E402

KNN = 0.5
LEVEL_COLOR = {"none": "#C00000", "gdna1": "#D55E00", "gdna5": "#D55E00"}


def panel_grid(n):
    cols = 8 if n > 16 else 4
    return int(np.ceil(n / cols)), cols


def main(suite="ambig"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scen = sorted(L.load_scenarios(suite), key=lambda s: (s["group"][1], s["group"][0], s["group"][2]))
    rows, cols = panel_grid(len(scen))
    fig, axes = plt.subplots(rows, cols, figsize=(3.05 * cols, 2.75 * rows), squeeze=False)
    for ax, s in zip(axes.ravel(), scen):
        mk = L.masks(s)
        sel = region_substrate(s, mk)
        fit = L.recipe(s, sel=sel, w=L.recipe_weights(s, sel, mk), knn_scale=KNN)
        corr = L.oracle_landscape(s, knn_scale=KNN)
        comb = L.oracle_landscape(s)
        if fit is not None:
            ax.fill_between(L.GRID, 1e-9, fit, color="#56B4E9", alpha=0.55, lw=0, label="FIT")
        if comb is not None:
            ax.plot(L.GRID, comb, color="#D55E00", lw=0.8, alpha=0.85, label="ref: PRODUCTION")
        if corr is not None:
            ax.plot(L.GRID, corr, color="0.12", lw=1.7, label="ref: CORRECTED")
        cap, lvl, ss, nr = s["group"]
        col = LEVEL_COLOR.get(lvl, "0.1")
        ax.set_title(f"{lvl} · cap {cap} · ss {ss} · {nr}", fontsize=7.5, color=col,
                     fontweight="bold" if lvl in LEVEL_COLOR else "normal")
        ax.set_yscale("log")
        ax.set_ylim(2e-6, 0.5)
        ax.set_xlim(-5, 2.5)
        ax.axvspan(0.0, 2.5, color="0.86", alpha=0.35, zorder=0)
        ax.tick_params(labelsize=6.5)
        ax.grid(alpha=0.22, lw=0.5)
        ax.set_axisbelow(True)
        if fit is not None and corr is not None:
            hi = L.GRID > L.two_component(corr)["split"]
            ax.text(0.03, 0.05, f"enr {float(fit[hi].sum()):.3f} / {float(corr[hi].sum()):.3f}\n"
                                f"rough {roughness(fit, 1.0, 2.5):.1f} vs {roughness(corr, 1.0, 2.5):.1f}",
                    transform=ax.transAxes, fontsize=6, color="0.25", va="bottom")
    for ax in axes.ravel()[len(scen):]:
        ax.axis("off")
    axes[0][0].legend(fontsize=6.5, loc="upper left", framealpha=0.9)
    fig.suptitle(f"N5 — the full gamut, {suite} suite ({len(scen)} conditions). FIT (regions only, kNN 0.5) "
                 f"vs the CORRECTED reference; thin orange is the PRODUCTION reference N5 replaced.\n"
                 f"Titles in RED = zero gDNA (nothing to find), ORANGE = low gDNA (the hardest regime). "
                 f"Shaded band is log10 rho > 0, where the comb lived.", fontsize=12)
    fig.supxlabel("log$_{10}$ $\\rho_g$", fontsize=10)
    fig.supylabel("landscape density (log scale)", fontsize=10)
    fig.tight_layout(rect=(0.012, 0.015, 1, 0.945 if rows > 4 else 0.90))
    out = f"/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_n5_gamut_{suite}.png"
    fig.savefig(out, dpi=100)
    print(f"wrote {out}")


if __name__ == "__main__":
    su = sys.argv[sys.argv.index("--suite") + 1] if "--suite" in sys.argv else "ambig"
    main(su)
