"""N5b — QUANTIFY the comb, and find what drives it (fit AND oracle).

The band medians (gdna_n5_bandwidth_diag.py) show the FIT's kernel collapsing 41x across the axis and going
sub-grid above +1.0. They do NOT explain the ORACLE's spikiness, whose median widths look adequate. So the
median is the wrong statistic: the landscape's HEIGHT at x is a sum of unit-mass kernels, and a node whose
kernel is narrower than the grid puts all of its mass in ONE CELL. Such a node contributes a spike towering
over the smooth background no matter how few of them there are.

Measured here per band:
  delta_frac   share of NODES whose kernel sd is below GRID_H — they cannot be rendered as anything but
               a delta, whatever the estimator
  delta_mass   the share of the band's MASS those nodes carry
  roughness    total variation of log P per decade, restricted to the band — a height-free measure of
               "spiky", directly comparable between curves

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n5_roughness.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds  # noqa: E402

_EPS = 1e-12
BANDS = [(-5, -2), (-2, -1), (-1, 0), (0, 1), (1, 2.5)]


def region_substrate(s, mk):
    """The substrate with BOUNDARIES EXCLUDED (owner decision, 2026-07-27)."""
    return (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]


def roughness(P, lo, hi):
    """Total variation of log P per decade inside a band — height-free, so a tall smooth bump scores 0
    and a comb of the same height scores large."""
    m = (L.GRID > lo) & (L.GRID <= hi)
    if m.sum() < 3:
        return np.nan
    lp = np.log(np.maximum(P[m], 1e-30))
    return float(np.abs(np.diff(lp)).sum() / (hi - lo))


def delta_stats(g, E, lo, hi):
    """Share of nodes (and of mass) in the band whose Poisson kernel is narrower than the grid can render."""
    a = np.log10(np.maximum(g, 1.0)) - np.log10(np.maximum(E, _EPS))
    m = (a > lo) & (a <= hi)
    if m.sum() < 1:
        return np.nan, np.nan, 0
    sd = 1.0 / (np.sqrt(np.maximum(g[m], 1.0)) * L.LN10)
    d = sd < L.GRID_H
    return float(d.mean()), float(g[m][d].sum() / max(g[m].sum(), _EPS)), int(m.sum())


def main():
    ss = conds()
    print("=== A. Who is rendered as a DELTA (kernel sd < GRID_H = %.4f dec) ===\n" % L.GRID_H)
    print(f"{'band':>14s} | {'FIT: n':>7s} {'delta_frac':>11s} {'delta_mass':>11s} | "
          f"{'ORACLE: n':>10s} {'delta_frac':>11s} {'delta_mass':>11s}")
    for lo, hi in BANDS:
        f, o = [], []
        for s in ss:
            mk = L.masks(s)
            sel = region_substrate(s, mk)
            f.append(delta_stats(s["g_hat"][sel], s["eff"][sel], lo, hi))
            osel = L.live_region(s)
            o.append(delta_stats(s["G"][osel], s["eff"][osel], lo, hi))
        fa, oa = np.array(f, float), np.array(o, float)
        print(f"{f'({lo:+.0f}, {hi:+.1f}]':>14s} | {np.nanmean(fa[:, 2]):7.0f} "
              f"{np.nanmean(fa[:, 0]):11.3f} {np.nanmean(fa[:, 1]):11.3f} | "
              f"{np.nanmean(oa[:, 2]):10.0f} {np.nanmean(oa[:, 0]):11.3f} {np.nanmean(oa[:, 1]):11.3f}")

    print("\n=== B. ROUGHNESS — total variation of log P per decade (height-free) ===\n")
    curves = [
        ("ORACLE (production reference)", lambda s, mk: L.oracle_landscape(s)),
        ("fit h=0 (no resolution term)", lambda s, mk: L.recipe(s, sel=region_substrate(s, mk),
                                                                w=L.recipe_weights(s, region_substrate(s, mk), mk))),
        ("fit kNN 0.5 (production)", lambda s, mk: L.recipe(s, sel=region_substrate(s, mk),
                                                            w=L.recipe_weights(s, region_substrate(s, mk), mk),
                                                            knn_scale=0.5)),
        ("fit kNN 2.0", lambda s, mk: L.recipe(s, sel=region_substrate(s, mk),
                                               w=L.recipe_weights(s, region_substrate(s, mk), mk),
                                               knn_scale=2.0)),
    ]
    hdr = " ".join(f"{f'({lo:+.0f},{hi:+.1f}]':>13s}" for lo, hi in BANDS)
    print(f"{'curve':32s} {hdr}")
    for name, fn in curves:
        vals = []
        for lo, hi in BANDS:
            r = []
            for s in ss:
                P = fn(s, L.masks(s))
                if P is not None:
                    r.append(roughness(P, lo, hi))
            vals.append(np.nanmean(r))
        print(f"{name:32s} " + " ".join(f"{v:13.1f}" for v in vals))
    print("\n  A smooth unimodal bump scores ~2-4 per decade. Anything above ~20 is a comb.")


if __name__ == "__main__":
    main()
