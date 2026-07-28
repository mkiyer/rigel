"""N5 — WHY THE LANDSCAPE IS SMOOTH BELOW log10 rho = 0 AND A COMB ABOVE IT.

Owner observation (2026-07-27): on a log-y axis the landscape is beautifully smooth through the depleted
bulk and becomes violently spiky above 0, worse the higher you go — AND THE ORACLE DOES IT TOO. Hypothesis:
a bandwidth that is effectively LINEAR-scale while the axis is log.

This measures the kernel width every node actually contributes, as a function of where it sits.

THE ARITHMETIC BEING TESTED.  The per-node kernel is the Poisson likelihood of `g` given rate `rho*E`,
rendered on x = log10 rho.  Its curvature at the mode is

    d2/dx2 [ g*(x*ln10 + ln E) - E*10^x ]  =  -g*(ln10)^2        =>   sd_x = 1/(sqrt(g)*ln10)  decades

so the width **shrinks as rho^(-1/2)** at fixed E: a node ten times denser carries ten times the counts and
so a sqrt(10)-times narrower kernel.  That is a measurement precision masquerading as a population
resolution, and it collapses exactly where the enriched mode lives.

Three scales are compared at each x:
    sd_poisson  the kernel the recipe actually uses
    GRID_H      the axis's own resolution — below this a node is a DELTA and cannot be rendered at all
    h_amise     1.06 * sd_hat * n^(-1/5) over the LOCAL population — what a sample of that size supports

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n5_bandwidth_diag.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds  # noqa: E402

_EPS = 1e-12
LN10 = L.LN10
BANDS = [(-5, -2), (-2, -1), (-1, 0), (0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.5)]


def knn_h(a, scale=0.5):
    """The production kNN width, verbatim from `gdna_explore_lib.recipe`."""
    n = a.size
    kk = max(int(round(np.sqrt(n))), 2)
    srt = np.sort(a)
    idx = np.searchsorted(srt, a)
    lo, hi = np.maximum(idx - kk, 0), np.minimum(idx + kk, n - 1)
    return scale * np.maximum(np.maximum(a - srt[lo], srt[hi] - a), L.GRID_H)


def main():
    ss = conds()
    print("=== Per-node kernel width vs position on the log10-rho axis ===")
    print(f"    GRID_H (the axis's own resolution) = {L.GRID_H:.4f} decades. Below it, a node is a DELTA.\n")
    print(f"{'band of log10 rho':>18s} {'n':>7s} {'med count g':>12s} {'sd_poisson':>11s} "
          f"{'/GRID_H':>8s} {'kNN h (0.5)':>12s} {'h_amise':>9s} {'amise/poisson':>14s}")
    for lo, hi in BANDS:
        N, G, SD, KH, AM = [], [], [], [], []
        for s in ss:
            mk = L.masks(s)
            # DECISION 2026-07-27 (owner): boundaries are EXCLUDED from the training substrate.
            sel = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
            g, E = s["g_hat"][sel], s["eff"][sel]
            a = np.log10(np.maximum(g, 1.0)) - np.log10(np.maximum(E, _EPS))
            h = knn_h(a)
            m = (a > lo) & (a <= hi)
            if m.sum() < 3:
                continue
            N.append(int(m.sum()))
            G.append(float(np.median(g[m])))
            SD.append(float(np.median(1.0 / (np.sqrt(np.maximum(g[m], 1.0)) * LN10))))
            KH.append(float(np.median(h[m])))
            # what a sample of THIS many points, with THIS spread, can actually resolve
            AM.append(1.06 * float(np.std(a[m])) * max(m.sum(), 2) ** (-0.2))
        if not N:
            continue
        sd, am = np.mean(SD), np.mean(AM)
        print(f"{f'({lo:+.1f}, {hi:+.1f}]':>18s} {np.mean(N):7.0f} {np.mean(G):12.1f} {sd:11.4f} "
              f"{sd / L.GRID_H:8.2f} {np.mean(KH):12.4f} {am:9.4f} {am / max(sd, 1e-9):14.1f}")

    print("\n=== The SAME measurement on the ORACLE reference (true counts G) ===")
    print(f"{'band of log10 rho':>18s} {'n':>7s} {'med count G':>12s} {'sd_poisson':>11s} "
          f"{'/GRID_H':>8s} {'h_amise':>9s} {'amise/poisson':>14s}")
    for lo, hi in BANDS:
        N, G, SD, AM = [], [], [], []
        for s in ss:
            sel = L.live_region(s)
            Gc, E = s["G"][sel], s["eff"][sel]
            a = np.log10(np.maximum(Gc, 1.0)) - np.log10(np.maximum(E, _EPS))
            m = (a > lo) & (a <= hi)
            if m.sum() < 3:
                continue
            N.append(int(m.sum()))
            G.append(float(np.median(Gc[m])))
            SD.append(float(np.median(1.0 / (np.sqrt(np.maximum(Gc[m], 1.0)) * LN10))))
            AM.append(1.06 * float(np.std(a[m])) * max(m.sum(), 2) ** (-0.2))
        if not N:
            continue
        sd, am = np.mean(SD), np.mean(AM)
        print(f"{f'({lo:+.1f}, {hi:+.1f}]':>18s} {np.mean(N):7.0f} {np.mean(G):12.1f} {sd:11.4f} "
              f"{sd / L.GRID_H:8.2f} {am:9.4f} {am / max(sd, 1e-9):14.1f}")

    print("\n  sd_poisson = 1/(sqrt(g)*ln10): the kernel SHRINKS as rho^(-1/2), so the axis's right-hand end")
    print("  is rendered with a width the axis cannot even represent. The ORACLE has the same defect, and it")
    print("  is worse there: G is the TRUTH, not a measurement, so its kernel is not a posterior at all.")


if __name__ == "__main__":
    main()
