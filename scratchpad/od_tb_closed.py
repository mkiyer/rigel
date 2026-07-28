"""(a) validate the closed-form bias bound; (b) the exact truncated-moment correction and its
conditioning; (c) the contamination hazard of correcting; (d) implementation realities."""
from __future__ import annotations

import numpy as np
import od_tb_core as T
from od_tb_agg import D, estimator_terms, seed_view, size_hist, pooled_bias
from od_tb_bias import ALPHA_BONF, min_n_rejectable

CEIL = 0.2


def part_a():
    print("=" * 100)
    print("(a) CLOSED-FORM BOUND   1 - r  <=  q(1-od)/(od(1-q))   [q = P(reject | truth)]")
    print("    exact vs bound, per-seed, a_null = 3")
    print(f"    {'n':>6s} {'od':>7s} {'alpha':>9s} {'q':>10s} {'1-r exact':>11s} "
          f"{'bound':>11s} {'ratio':>7s}")
    for alpha in (ALPHA_BONF, 1e-3, 0.01, 0.05):
        for n in (300, 1523):
            for od in (0.05, 0.1, 0.2):
                pk, e_all, e_kept, r = T.seed_moments(n, od, 3.0, alpha)
                q = 1.0 - pk
                bound = q * (1 - od) / (od * (1 - q)) if q > 0 else 0.0
                exact = 1.0 - r
                print(f"    {n:6d} {od:7.3f} {alpha:9.3g} {q:10.3e} {exact:11.3e} "
                      f"{bound:11.3e} {exact/bound if bound>0 else np.nan:7.3f}")
    print()


def g_of_od(u, c, a_null, alpha, od):
    """E[ sum over ALL seeds of excess*1{keep} ] under true dispersion od -- the estimating fn."""
    nmin = min_n_rejectable(a_null, alpha)
    tot = 0.0
    for n, cnt in zip(u, c):
        n = int(n)
        if n < 2:
            continue
        scale = n * (n - 1) / 4.0
        if nmin < 0 or n < nmin:
            tot += cnt * scale * od
        else:
            tot += cnt * T.seed_moments(n, od, a_null, alpha)[2]
    return tot


def part_b():
    print("=" * 100)
    print("(b) THE EXACT CORRECTION: solve  sum_s E[excess_s 1{keep_s}](od) = sum_{kept} excess_s")
    print("    Conditioning = d g / d od, relative to the uncorrected slope sum(scale).")
    print("    a_null = 3.  'slope ratio' < 1 => the correction AMPLIFIES noise by 1/ratio.")
    print(f"    {'sample':10s} {'alpha':>9s} {'od':>7s} {'g/(S*od)':>10s} {'slope ratio':>12s}")
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        _, sc, nc = estimator_terms(s, t, w, k)
        u, c = size_hist(nc)
        S = float(sum(int(cn) * int(n) * (int(n) - 1) / 4.0 for n, cn in zip(u, c) if n >= 2))
        for alpha in (ALPHA_BONF, 0.01, 0.05):
            for od in (0.05, 0.15):
                h = 1e-4
                g0 = g_of_od(u, c, 3.0, alpha, od)
                gp = g_of_od(u, c, 3.0, alpha, od + h)
                gm = g_of_od(u, c, 3.0, alpha, od - h)
                print(f"    {name:10s} {alpha:9.3g} {od:7.3f} {g0/(S*od):10.6f} "
                      f"{(gp-gm)/(2*h)/S:12.6f}")
    print()


def part_c():
    print("=" * 100)
    print("(c) THE HAZARD OF CORRECTING: contaminated seeds are rejected with prob ~1, but the")
    print("    correction credits them the null's keep probability -> od is pulled DOWN.")
    print("    Simulation: N_legit seeds ~ BetaBinom(n, od_true), plus n_bad contaminated seeds")
    print("    (all-antisense, size n_bad_size), a_null = 3, alpha = Bonferroni.")
    rng = np.random.default_rng(0)
    a_null, alpha = 3.0, ALPHA_BONF
    n_legit, n_size, od_true = 2000, 400, 0.10
    for n_bad in (0, 5, 20):
        a_t = T.a_of_od(od_true)
        p = rng.beta(a_t, a_t, n_legit)
        kk = rng.binomial(n_size, p)
        sizes = np.full(n_legit, n_size)
        if n_bad:
            sizes = np.concatenate([sizes, np.full(n_bad, n_size)])
            kk = np.concatenate([kk, np.zeros(n_bad, int)])
        scale = sizes * (sizes - 1) / 4.0
        excess = (kk - sizes / 2.0) ** 2 - sizes / 4.0
        tp = T.two_sided_p(n_size, a_null)
        keep = tp[np.clip(kk, 0, n_size)] >= alpha
        raw = excess.sum() / scale.sum()
        scr = excess[keep].sum() / scale[keep].sum()
        u, c = np.array([n_size]), np.array([sizes.size])
        lo, hi = 1e-6, 0.999
        rhs = excess[keep].sum()
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g_of_od(u, c, a_null, alpha, mid) < rhs:
                lo = mid
            else:
                hi = mid
        print(f"    contaminants={n_bad:3d}  raw={raw:8.5f}  screened={scr:8.5f}  "
              f"trunc-corrected={0.5*(lo+hi):8.5f}   (truth {od_true})")
    print()


def part_d():
    print("=" * 100)
    print("(d) IMPLEMENTATION REALITIES on the real seed arrays")
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        nc = w * t
        nonint_n = float(np.mean(np.abs(t - np.rint(t)) > 1e-9))
        nonint_k = float(np.mean(np.abs(s - np.rint(s)) > 1e-9))
        print(f"    {name:10s} w: min {w.min():.4f} p1 {np.percentile(w,1):.4f} "
              f"med {np.median(w):.4f} mean {w.mean():.4f} | frac w<0.99 {np.mean(w<0.99):.4%} "
              f"| non-integer total {nonint_n:.2%} sense {nonint_k:.2%} "
              f"| max n_c {nc.max():.1f}")
    print()


if __name__ == "__main__":
    part_a()
    part_d()
    part_c()
    part_b()
