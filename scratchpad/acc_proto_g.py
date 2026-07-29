"""EXPERIMENT G — which second moment? sum(L) or sum(1/L)?

At a genomic POINT the crossing opportunity for a length-w fragment is exactly w
start positions, so depositing 1/w cancels it EXACTLY:

    E[ sum 1/L ]  =  sum_c sum_w rho_c * w * f_c(w) * (1/w)  =  sum_c rho_c  =  rho_total

-> a DIVISOR-FREE, FL-MODEL-FREE estimate of total molecular density. No E_c anywhere.

The count gives the second equation:      E[N] = rho_g E_g[w] + rho_r E_r[w]
so (N, sum 1/L) identifies (rho_g, rho_r) iff the two MEAN lengths differ.

Compare against v3's (N, sum L), whose second row is E_c[w^2].

At a NODE the opportunity is max(0, B-w+1), which 1/w does NOT cancel — so the
identity is an EDGE property. Tested below.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260801)


def fl(mu, sd, n=1600):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


W = np.arange(1600, dtype=np.float64)


def edge_rows(pmf):
    """(E[w], E[w^2], 1.0) — the count row, the sum-L row, and the sum-1/L row."""
    return np.sum(pmf * W), np.sum(pmf * W * W), 1.0


def node_rows(pmf, B):
    fit = np.maximum(0.0, B - W + 1.0)
    inv = np.divide(1.0, W, out=np.zeros_like(W), where=W > 0)
    return np.sum(pmf * fit), np.sum(pmf * W * fit), np.sum(pmf * fit * inv)


def crossing_pmf(pmf, mode, B=None):
    q = pmf * W if mode == "edge" else pmf * np.maximum(0.0, B - W + 1.0)
    return q / q.sum()


def trial(pg, pr, rho, scale, mode, B=None, Mg=None, Mr=None, mis=None):
    """Draw one realisation; return (N, sumL, sum1/L) and the true f_g."""
    lam_g = rho[0] * (Mg[0] if mode == "edge" else Mg[0]) * scale
    lam_r = rho[1] * (Mr[0] if mode == "edge" else Mr[0]) * scale
    ng, nr = RNG.poisson(lam_g), RNG.poisson(lam_r)
    wg = RNG.choice(pg.size, size=ng, p=crossing_pmf(pg, mode, B))
    wr = RNG.choice(pr.size, size=nr, p=crossing_pmf(pr, mode, B))
    w = np.concatenate([wg, wr]).astype(np.float64)
    n = w.size
    return np.array([n, w.sum(), np.sum(1.0 / w) if n else 0.0]) / scale


def solve(obs2, M2):
    try:
        return np.maximum(np.linalg.solve(M2, obs2), 0.0)
    except np.linalg.LinAlgError:
        return np.array([np.nan, np.nan])


def compare(gmu, rmu, fg, mode, B=None, scale=500.0, reps=800, fl_err=0.0):
    pg, pr = fl(gmu, gmu * 0.25), fl(rmu, rmu * 0.25)
    rows = edge_rows if mode == "edge" else (lambda p: node_rows(p, B))
    Mg, Mr = np.array(rows(pg)), np.array(rows(pr))
    # the design the SOLVER uses — possibly built on a mis-specified gDNA FL
    pg_used = fl(gmu * (1 + fl_err), gmu * 0.25 * (1 + fl_err)) if fl_err else pg
    Mgu = np.array(rows(pg_used))
    rho = np.array([fg, 1 - fg]) * 0.02

    A = np.column_stack([Mgu[[0, 1]], Mr[[0, 1]]])   # (N, sum L)
    Bm = np.column_stack([Mgu[[0, 2]], Mr[[0, 2]]])  # (N, sum 1/L)
    ea, eb, tot_direct, tot_a = [], [], [], []
    for _ in range(int(reps)):
        o = trial(pg, pr, rho, scale, mode, B, Mg, Mr)
        ra = solve(o[[0, 1]], A)
        rb = solve(o[[0, 2]], Bm)
        if np.isfinite(ra).all() and ra.sum() > 0:
            ea.append(ra[0] / ra.sum())
            tot_a.append(ra.sum())
        if np.isfinite(rb).all() and rb.sum() > 0:
            eb.append(rb[0] / rb.sum())
        tot_direct.append(o[2])              # sum(1/L) alone, no model at all
    ea, eb = np.asarray(ea), np.asarray(eb)
    td, ta = np.asarray(tot_direct), np.asarray(tot_a)
    rho_tot = rho.sum()
    tag = f"{mode}" + (f" B={B}" if B else "")
    err = f" [gDNA FL +{100*fl_err:.0f}%]" if fl_err else ""
    print(f"  {tag:11s} gDNA {gmu:3.0f}/RNA {rmu:3.0f} f_g {fg:.2f}{err}")
    print(f"      (N, sum L)    f_g {ea.mean():.3f} +/- {ea.std():.3f}  bias {ea.mean()-fg:+.4f}"
          f"   cond {np.linalg.cond(A):8.1f}")
    print(f"      (N, sum 1/L)  f_g {eb.mean():.3f} +/- {eb.std():.3f}  bias {eb.mean()-fg:+.4f}"
          f"   cond {np.linalg.cond(Bm):8.1f}")
    if mode == "edge":
        print(f"      rho_total: sum(1/L) alone {td.mean()/rho_tot:.4f}x truth "
              f"(+/-{td.std()/rho_tot:.4f})   vs (N,sumL) model {ta.mean()/rho_tot:.4f}x")


if __name__ == "__main__":
    print("=" * 84)
    print("EDGE — does sum(1/L) beat sum(L) as the second moment?")
    print("=" * 84)
    for rmu in (200, 150, 100, 70):
        compare(50, rmu, 0.15, "edge")
    print()
    print("=" * 84)
    print("ROBUSTNESS — the FL model is MIS-SPECIFIED by 10 %")
    print("=" * 84)
    compare(50, 200, 0.15, "edge", fl_err=0.10)
    compare(50, 200, 0.50, "edge", fl_err=0.10)
    print()
    print("=" * 84)
    print("NODE — the 1/L identity is an EDGE property; does it still work contained?")
    print("=" * 84)
    for B in (300, 1000, 5000):
        compare(50, 200, 0.15, "node", B=B)
