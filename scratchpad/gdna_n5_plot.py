"""N5 — the figure: why the landscape combs above log10 rho = 0, and what fixes it.

Top row    — three conditions on a log-y axis: the PRODUCTION reference (a comb), the CORRECTED reference
             (population resolution), and the production fit. The fit was already protected by W2's kNN
             term; the reference was not — and the reference is the scoring instrument.
Bottom row — the mechanism (kernel width collapsing as rho^(-1/2) and going sub-grid), who that turns into
             a delta, and the shape-based bandwidth selection against ground truth.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n5_plot.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds  # noqa: E402
from gdna_n5_bandwidth_diag import BANDS  # noqa: E402
from gdna_n5_reference import SCALES, truth_stats  # noqa: E402
from gdna_n5_roughness import region_substrate, roughness  # noqa: E402

_EPS = 1e-12
OUT = "/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_n5_bandwidth.png"
KNN = 0.5
PICKS = ["gdna_gdna300_ss_0.50_nrna_none_capture_on",
         "gdna_gdna100_ss_0.99_nrna_none_capture_on",
         "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong"]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ss = conds()
    scen = {s["cond"]: s for s in ss}
    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.4))

    # ── top row: the reference, before and after ───────────────────────────────────────────────────
    for ax, cond in zip(axes[0], PICKS):
        s = scen[cond]
        mk = L.masks(s)
        sel = region_substrate(s, mk)
        fit = L.recipe(s, sel=sel, w=L.recipe_weights(s, sel, mk), knn_scale=KNN)
        ax.fill_between(L.GRID, 1e-9, fit, color="#56B4E9", alpha=0.55, lw=0,
                        label=f"FIT, kNN {KNN} (regions only)")
        ax.plot(L.GRID, L.oracle_landscape(s), color="#D55E00", lw=1.1,
                label="reference, PRODUCTION (measurement kernel)")
        ax.plot(L.GRID, L.oracle_landscape(s, knn_scale=KNN), color="0.12", lw=2.2,
                label="reference, CORRECTED (population resolution)")
        ax.set_yscale("log")
        ax.set_ylim(2e-6, 0.4)
        ax.set_xlim(-5, 2.5)
        ax.axvspan(0.0, 2.5, color="0.85", alpha=0.35, zorder=0)
        ax.set_title(" · ".join(s["group"]), fontsize=9.5)
        ax.set_xlabel("log$_{10}$ $\\rho_g$", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0][0].legend(fontsize=7.5, loc="lower left", framealpha=0.93)
    axes[0][0].set_ylabel("landscape density (log)", fontsize=9)

    # ── bottom-left: the mechanism ─────────────────────────────────────────────────────────────────
    ax = axes[1][0]
    mid = [0.5 * (lo + hi) for lo, hi in BANDS]
    sd_f, sd_o, kh = [], [], []
    for lo, hi in BANDS:
        a_, b_, c_ = [], [], []
        for s in ss:
            mk = L.masks(s)
            sel = region_substrate(s, mk)
            g, E = s["g_hat"][sel], s["eff"][sel]
            a = np.log10(np.maximum(g, 1.0)) - np.log10(np.maximum(E, _EPS))
            m = (a > lo) & (a <= hi)
            if m.sum() >= 3:
                a_.append(np.median(1.0 / (np.sqrt(np.maximum(g[m], 1.0)) * L.LN10)))
                c_.append(np.median(L.knn_widths(a, KNN)[m]))
            osel = L.live_region(s)
            G, Eo = s["G"][osel], s["eff"][osel]
            ao = np.log10(np.maximum(G, 1.0)) - np.log10(np.maximum(Eo, _EPS))
            mo = (ao > lo) & (ao <= hi)
            if mo.sum() >= 3:
                b_.append(np.median(1.0 / (np.sqrt(np.maximum(G[mo], 1.0)) * L.LN10)))
        sd_f.append(np.mean(a_) if a_ else np.nan)
        sd_o.append(np.mean(b_) if b_ else np.nan)
        kh.append(np.mean(c_) if c_ else np.nan)
    ax.plot(mid, sd_f, "o-", color="#D55E00", lw=2, label="Poisson MEASUREMENT kernel, $1/(\\sqrt{g}\\ln 10)$")
    ax.plot(mid, kh, "s-", color="#0072B2", lw=2, label=f"kNN POPULATION resolution (scale {KNN})")
    ax.axhline(L.GRID_H, color="0.2", ls="--", lw=1.4)
    ax.text(0.4, L.GRID_H * 1.15, "GRID_H — below this a node is a DELTA", fontsize=8, color="0.2",
            ha="right")
    ax.set_yscale("log")
    ax.set_xlabel("log$_{10}$ $\\rho_g$", fontsize=9)
    ax.set_ylabel("kernel width (decades)", fontsize=9)
    ax.set_title("the mechanism: the measurement width collapses as $\\rho^{-1/2}$", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    # ── bottom-middle: roughness by band ───────────────────────────────────────────────────────────
    ax = axes[1][1]
    curves = [("reference, PRODUCTION", "#D55E00", lambda s, mk: L.oracle_landscape(s)),
              ("reference, CORRECTED", "0.12", lambda s, mk: L.oracle_landscape(s, knn_scale=KNN)),
              ("fit h=0", "#E69F00",
               lambda s, mk: L.recipe(s, sel=region_substrate(s, mk),
                                      w=L.recipe_weights(s, region_substrate(s, mk), mk))),
              ("fit kNN 0.5", "#0072B2",
               lambda s, mk: L.recipe(s, sel=region_substrate(s, mk),
                                      w=L.recipe_weights(s, region_substrate(s, mk), mk), knn_scale=KNN))]
    rb = [(-5, -2), (-2, -1), (-1, 0), (0, 1), (1, 2.5)]
    x = np.arange(len(rb))
    for i, (nm, col, fn) in enumerate(curves):
        v = [np.nanmean([roughness(fn(s, L.masks(s)), lo, hi) for s in ss]) for lo, hi in rb]
        ax.bar(x + (i - 1.5) * 0.2, v, 0.19, color=col, label=nm)
    ax.axhspan(0, 4, color="#009E73", alpha=0.15, zorder=0)
    ax.text(4.4, 4.6, "smooth-bump band", fontsize=8, color="#006B4F", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([f"({lo:+.0f},{hi:+.1f}]" for lo, hi in rb], fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("band of log$_{10}$ $\\rho_g$", fontsize=9)
    ax.set_ylabel("roughness: TV of log P per decade", fontsize=9)
    ax.set_title("the comb, quantified — height-free", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    # ── bottom-right: shape-based bandwidth selection ──────────────────────────────────────────────
    ax = axes[1][2]
    t = np.array([truth_stats(s, 50.0) for s in ss])
    t_lo, t_hi = np.nanmean(t[:, 3]) - 0.03, np.nanmean(t[:, 3]) + 0.03
    ref, fitw = [], []
    for sc in SCALES:
        r, f = [], []
        for s in ss:
            orc = L.oracle_landscape(s, knn_scale=sc)
            mk = L.masks(s)
            sel = region_substrate(s, mk)
            P = L.recipe(s, sel=sel, w=L.recipe_weights(s, sel, mk), knn_scale=sc)
            if orc is None:
                continue
            tco = L.two_component(orc)
            r.append(tco["enr_width"])
            if P is not None:   # score the fit at the REFERENCE's split, the standing convention
                f.append(L.two_component(P, split=tco["split"])["enr_width"])
        ref.append(np.nanmean(r))
        fitw.append(np.nanmean(f))
    xs = [max(sc, 0.12) for sc in SCALES]
    ax.axhspan(t_lo, t_hi, color="#009E73", alpha=0.2, zorder=0,
               label=f"TRUTH: enriched sd {np.nanmean(t[:, 3]):.2f} dec")
    ax.plot(xs, ref, "o-", color="0.12", lw=2, label="reference, rendered at scale")
    ax.plot(xs, fitw, "s-", color="#0072B2", lw=2, label="fit, rendered at scale")
    ax.axvline(KNN, color="#D55E00", ls="--", lw=1.5)
    ax.text(KNN * 1.1, 0.155, "kNN 0.5\n(W2's choice, now\nvalidated on truth)", fontsize=8,
            color="#D55E00")
    ax.set_ylim(0.12, 0.48)
    ax.set_xscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(["h=0"] + [str(s) for s in SCALES[1:]], fontsize=8)
    ax.set_xlabel("kNN population-resolution scale", fontsize=9)
    ax.set_ylabel("rendered enriched-mode width (decades)", fontsize=9)
    ax.set_title("selecting the bandwidth by SHAPE, not by EMD", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    fig.suptitle("N5 — the landscape's resolution was a MEASUREMENT width (~$1/\\sqrt{g}$), so it collapsed "
                 "exactly where the enriched mode lives. The FIT was already protected; the REFERENCE was not.",
                 fontsize=12.5)
    fig.tight_layout(rect=(0.005, 0.005, 1, 0.955))
    fig.savefig(OUT, dpi=110)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
