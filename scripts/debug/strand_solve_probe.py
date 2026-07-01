"""Isolate the strand-solve at the UNIT level (no messages, no global, no dataset, no noise).

A pure-gDNA single-strand node has EXACT 50/50 per-strand counts → +frac=0.5 deterministically. The
strand mixture (gDNA mean 0.5, RNA mean κ) should call it f_g≈1 for any κ away from 0.5 (0.5 ≠ κ). Any
f_g<1 here is a pure MODEL bug, not noise/dataset. Sweep κ × od_g × od_r × N. Also a pure-RNA control
(+frac=κ → want f_g≈0) and a known-mixture control.
"""
from __future__ import annotations
import numpy as np
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all


def solve(u_pos, u_neg, kappa, od_g, od_r, n_grid=60):
    z = np.zeros(1)
    n = float(u_pos + u_neg)
    d = _solve_nodes_logodds_all(
        np.array([float(u_pos)]), np.array([float(u_neg)]), z, z,
        np.array([True]), np.array([False]), np.array([True]),
        np.array([n]), z, kappa=kappa, od_g=od_g, od_r=od_r, n_grid=n_grid)
    return float(d.gdna_frac[0])


def pure_gdna(N, kappa, od_g, od_r):  # exact 50/50 → should be f_g≈1
    return solve(N / 2.0, N / 2.0, kappa, od_g, od_r)


def pure_rna(N, kappa, od_g, od_r):   # +frac=κ exactly → should be f_g≈0
    return solve(round(N * kappa), round(N * (1 - kappa)), kappa, od_g, od_r)


if __name__ == "__main__":
    print("=== PURE gDNA node (exact 50/50; WANT f_g≈1.0) — f_g vs κ, od_g (N=1000, od_r=0.05) ===")
    print(f"{'kappa':>6} | " + " ".join(f"od_g={o:<4}" for o in (0.0, 0.05, 0.1, 0.2)))
    for kappa in (0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.1, 0.05, 0.01):
        row = [f"{pure_gdna(1000, kappa, o, 0.05):.3f}    " for o in (0.0, 0.05, 0.1, 0.2)]
        print(f"{kappa:>6} | " + " ".join(row))
    print("\n=== PURE gDNA, effect of N (od_g=0.1, od_r=0.05) ===")
    print(f"{'kappa':>6} | " + " ".join(f"N={n:<5}" for n in (50, 200, 1000, 5000)))
    for kappa in (0.5, 0.4, 0.3, 0.2, 0.1, 0.01):
        row = [f"{pure_gdna(n, kappa, 0.1, 0.05):.3f}" for n in (50, 200, 1000, 5000)]
        print(f"{kappa:>6} | " + "   ".join(row))
    print("\n=== PURE RNA control (+frac=κ exact; WANT f_g≈0.0) (N=1000, od_g=0.1, od_r=0.05) ===")
    for kappa in (0.5, 0.4, 0.3, 0.2, 0.1, 0.01):
        print(f"  κ={kappa}: f_g={pure_rna(1000, kappa, 0.1, 0.05):.3f}")
