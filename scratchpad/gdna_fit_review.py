"""PER-SCENARIO REVIEW of the gDNA prior fit — every scenario, fit vs ORACLE, plus a mode/bandwidth audit.

Two questions drive this (owner, 2026-07-26):
  1. IS THE FIT ROBUST — does it capture ALL the modes, at the right bandwidth, WITHOUT overfitting?
  2. DOES THE PROJECTION pull a node's observed density to the right gDNA level?

This file answers (1). It renders, for all 32 conditions, three curves on one log10-rho_g axis:
  * ORACLE      — `fit_poisson` on the TRUE per-node gDNA counts (the same estimator as the fit, so the
                  comparison isolates the INPUT counts + weights, not the smoother)
  * PRODUCTION  — the shipped `DensityNPMLE` (`calibrate._fit_gdna_hyperprior`), cached by gdna_cache_build
  * RECIPE      — the 2026-07-21 exploration landscape (`gdna_explore_lib.recipe`)

and audits each fit for: EMD to oracle, the enriched-region mass ratio, MODE recall/precision (are the
oracle's modes present? does the fit invent modes that are not there?) and the spread ratio (sd_fit/sd_oracle
< 1 = too sharp = overfitting the training counts; > 1 = over-smoothed).

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_fit_review.py [--suite ambig] [--out FIG.png]
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_explore_lib import modes, spread  # noqa: E402  — shared metrics; were duplicated here
from gdna_resurrect import load_recipe  # noqa: E402

GRID = L.GRID
ENR_THR = 0.0          # decades; above this is the capture-enriched gDNA level
MODE_TOL = 0.30        # decades; a fitted mode "matches" an oracle mode within this


def support(P, tol=1e-12):
    """The fit's own FITTED SUPPORT, recovered from the rendered curve.

    `DensityNPMLE` is fitted on its own grid and read back with `np.interp(..., left=logP[0],
    right=logP[-1])` — so outside that grid the solver sees a FLAT CLAMPED SHELF, not a real tail. Those
    shelves are exactly-constant runs at the ends of the rendered density, so they are recoverable here.
    Returns (lo, hi) in decades, or (GRID[0], GRID[-1]) when nothing is clamped."""
    if P is None:
        return GRID[0], GRID[-1]
    lo = 0
    while lo + 1 < len(P) and abs(P[lo + 1] - P[0]) <= tol * max(P.max(), 1e-30):
        lo += 1
    hi = len(P) - 1
    while hi - 1 > lo and abs(P[hi - 1] - P[-1]) <= tol * max(P.max(), 1e-30):
        hi -= 1
    return float(GRID[lo]), float(GRID[hi])


def audit(P, orc):
    """Compare a fit against the oracle landscape: EMD, enriched mass ratio, mode recall/precision, spread."""
    if P is None or orc is None:
        return dict(emd=np.nan, enr=np.nan, recall=np.nan, spurious=np.nan, sd_ratio=np.nan, n_modes=0)
    hi = GRID > ENR_THR
    om, fm = float(orc[hi].sum()), float(P[hi].sum())
    mo, mf = modes(orc), modes(P)
    # recall: oracle modes with a fitted mode within MODE_TOL. spurious: fitted modes with no oracle mode.
    rec = np.mean([any(abs(a - b) <= MODE_TOL for b, _ in mf) for a, _ in mo]) if mo else np.nan
    spu = np.mean([not any(abs(b - a) <= MODE_TOL for a, _ in mo) for b, _ in mf]) if mf else np.nan
    sdo = spread(orc)
    return dict(emd=L.emd(P, orc), enr=(fm / om if om > 1e-4 else np.nan), recall=rec, spurious=spu,
                sd_ratio=(spread(P) / sdo if sdo > 0 else np.nan), n_modes=len(mf), n_oracle_modes=len(mo))


def review(suite="ambig", out=None):
    rec_fn = load_recipe()
    scen = L.load_scenarios(suite)
    scen = sorted(scen, key=lambda s: (s["group"][0], s["group"][2], s["group"][1], s["group"][3]))
    rows = []
    for s in scen:
        orc = L.oracle_landscape(s)
        a_prod = audit(s.get("prod_P"), orc)
        a_rec = audit(rec_fn(s), orc)
        rows.append((s, orc, a_prod, a_rec))

    # ---------------- the table ----------------
    print(f"=== {suite}: gDNA prior fit vs ORACLE, per scenario "
          f"(PROD = shipped DensityNPMLE, REC = exploration recipe) ===")
    hdr = (f"{'condition':50s} | {'EMD':>5s} {'enr':>5s} {'mrec':>5s} {'spur':>5s} {'sd/o':>5s} {'nM':>3s}"
           f" | {'EMD':>5s} {'enr':>5s} {'mrec':>5s} {'spur':>5s} {'sd/o':>5s} {'nM':>3s} | {'oM':>3s}")
    print(hdr)
    print("-" * len(hdr))

    def f(v, w=5, p=2):
        return f"{v:{w}.{p}f}" if np.isfinite(v) else f"{'--':>{w}s}"

    for s, orc, ap, ar in rows:
        print(f"{s['cond']:50s} | {f(ap['emd'])} {f(ap['enr'])} {f(ap['recall'])} {f(ap['spurious'])} "
              f"{f(ap['sd_ratio'])} {ap['n_modes']:3d} | "
              f"{f(ar['emd'])} {f(ar['enr'])} {f(ar['recall'])} {f(ar['spurious'])} "
              f"{f(ar['sd_ratio'])} {ar['n_modes']:3d} | {ap['n_oracle_modes']:3d}")

    # the production fit's own SUPPORT vs where the oracle actually puts its mass — a coverage failure here
    # means the solver reads a flat clamped shelf, not a prior, over part of the density axis.
    print("\n  PRODUCTION fitted support vs the ORACLE's mass range (decades of log10 rho_g):")
    print(f"  {'condition':50s} {'fit lo':>7s} {'fit hi':>7s} | {'oracle p1':>9s} {'oracle p99':>10s}"
          f" | {'oracle mass BELOW fit lo':>25s}")
    for s, orc, _ap, _ar in rows:
        P = s.get("prod_P")
        if P is None:
            continue
        slo, shi = support(P)
        c = np.cumsum(orc)
        p1, p99 = float(GRID[np.searchsorted(c, 0.01)]), float(GRID[np.searchsorted(c, 0.99)])
        below = float(orc[GRID < slo].sum())
        flag = "  <== UNCOVERED" if below > 0.10 else ""
        print(f"  {s['cond']:50s} {slo:7.2f} {shi:7.2f} | {p1:9.2f} {p99:10.2f} | {below:25.3f}{flag}")

    for label, key in (("PRODUCTION", 2), ("RECIPE", 3)):
        vals = {k: np.array([r[key][k] for r in rows], float) for k in
                ("emd", "enr", "recall", "spurious", "sd_ratio")}
        print(f"\n  {label:10s} mean EMD {np.nanmean(vals['emd']):.3f} | enriched-mass ratio "
              f"{np.nanmean(vals['enr']):.2f} | mode recall {np.nanmean(vals['recall']):.2f} | "
              f"spurious-mode rate {np.nanmean(vals['spurious']):.2f} | sd/oracle "
              f"{np.nanmean(vals['sd_ratio']):.2f}")
    print("\n  mrec = fraction of ORACLE modes recovered (want 1.00) · spur = fraction of FITTED modes that are"
          "\n  not in the oracle (want 0.00 — this is the overfitting canary) · sd/o < 1 = too sharp")

    # ---------------- the figure ----------------
    if out is None:
        return rows
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(rows)
    ncol, nrow = 4, int(np.ceil(n / 4))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.1 * nrow), sharex=True)
    for ax, (s, orc, ap, ar) in zip(axes.ravel(), rows):
        P_prod, P_rec = s.get("prod_P"), rec_fn(s)
        ax.fill_between(GRID, 0, orc, color="0.78", lw=0, label="ORACLE")
        if P_prod is not None:
            slo, shi = support(P_prod)
            ax.axvspan(GRID[0], slo, color="C0", alpha=.07, lw=0)   # clamped shelf = outside the fitted grid
            ax.axvspan(shi, GRID[-1], color="C0", alpha=.07, lw=0)
            ax.plot(GRID, P_prod, "C0", lw=1.5, label="PRODUCTION (NPMLE)")
        if P_rec is not None:
            ax.plot(GRID, P_rec, "C3", lw=1.5, label="RECIPE")
        for loc, _h in modes(orc):
            ax.axvline(loc, color="0.45", ls=":", lw=0.9)
        cap, dna, ss, nr = s["group"]
        ax.set_title(f"{cap} · {dna} · ss{ss} · nrna={nr}\n"
                     f"EMD prod {ap['emd']:.3f} / rec {ar['emd']:.3f}   "
                     f"enr {ap['enr']:.2f} / {ar['enr']:.2f}", fontsize=8)
        ax.set_xlim(-5, 2.0)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=7)
    fig.suptitle(f"gDNA prior fit vs ORACLE — every {suite} scenario "
                 f"(dotted = oracle mode locations)", fontsize=13)
    fig.supxlabel("log10 rho_g   (gDNA density)")
    fig.tight_layout(rect=(0, 0.01, 1, 0.985))
    fig.savefig(out, dpi=100)
    print(f"\nwrote {out}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="ambig")
    ap.add_argument("--out", default="/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_fit_review.png")
    a = ap.parse_args()
    review(a.suite, a.out)
