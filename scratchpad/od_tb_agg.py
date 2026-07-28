"""Aggregate truncation bias on the REAL cfRNA seed-size distributions, and the
bias-vs-contamination trade, exactly."""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import od_tb_core as T
from od_tb_bias import ALPHA_BONF, min_n_rejectable

D = pickle.load(open(Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/5a99d9c9-dd1f-4d8a-b7a3-a9d62be385dc/scratchpad/seeds.pkl"), "rb"))
CEIL = 0.2


def seed_view(d):
    s, t, w, k = d["sense"], d["total"], d["weight"], d["kappa"]
    ok = t > 0
    return s[ok], t[ok], w[ok], k


def estimator_terms(s, t, w, k):
    node_mean = 0.5 * w + k * (1.0 - w)
    excess = (s - t * node_mean) ** 2 - t * node_mean * (1.0 - node_mean)
    nc = w * t
    scale = np.maximum(nc * (nc - 1.0), 0.0) * 0.25
    return excess, scale, nc


def size_hist(t):
    n = np.rint(t).astype(np.int64)
    u, c = np.unique(n, return_counts=True)
    return u, c


def pooled_bias(u, c, a_null, alpha, od):
    """Exact limit of the truncated estimator on this size multiset, both-sums version."""
    nmin = min_n_rejectable(a_null, alpha)
    num = 0.0
    den = 0.0
    num0 = 0.0
    qw = 0.0
    for n, cnt in zip(u, c):
        n = int(n)
        if n < 2:
            continue
        scale = n * (n - 1) / 4.0
        if nmin < 0 or n < nmin:  # cannot be rejected at all -> exact
            num += cnt * scale * od
            den += cnt * scale
            num0 += cnt * scale * od
            continue
        pk, e_all, e_kept, _ = T.seed_moments(n, od, a_null, alpha)
        num += cnt * e_kept
        den += cnt * scale * pk
        num0 += cnt * e_all
        qw += cnt * scale * (1.0 - pk)
    den_full = num0 / od if od > 0 else sum(int(cnt) * int(n) * (int(n) - 1) / 4.0
                                            for n, cnt in zip(u, c) if n >= 2)
    return dict(od_hat=num / den if den > 0 else np.nan,
                od_hat_numonly=num / den_full,
                q_scale_weighted=qw / den_full,
                frac_eligible_scale=sum(int(cnt) * int(n) * (int(n) - 1) / 4.0
                                        for n, cnt in zip(u, c)
                                        if nmin > 0 and n >= max(nmin, 2)) / den_full)


def main():
    alphas = [ALPHA_BONF, 1e-4, 1e-3, 1e-2, 5e-2]
    ods = [0.01, 0.0345, 0.1, 0.143, 0.2]

    print("=== REAL SEED SIZE DISTRIBUTIONS (n_c = gdna_weight * total, rounded) ===")
    print(f"{'sample':10s} {'seeds':>9s} {'n>=2':>9s} {'sum pairs':>13s} {'med':>5s} "
          f"{'p99':>7s} {'max':>7s} {'#n>=264':>8s} {'#n>=1385':>9s} {'sc>=264':>8s} {'sc>=1385':>9s}")
    views = {}
    for name, d in D.items():
        s, t, w, k = seed_view(d)
        ex, sc, nc = estimator_terms(s, t, w, k)
        views[name] = (s, t, w, k, ex, sc, nc)
        n = np.rint(nc).astype(np.int64)
        pairs = (n * (n - 1) / 2.0)
        tot_sc = sc.sum()
        print(f"{name:10s} {n.size:9,d} {int((n>=2).sum()):9,d} {pairs.sum():13,.0f} "
              f"{np.median(n):5.1f} {np.percentile(n,99):7.0f} {n.max():7d} "
              f"{int((n>=264).sum()):8,d} {int((n>=1385).sum()):9,d} "
              f"{sc[n>=264].sum()/tot_sc:8.3%} {sc[n>=1385].sum()/tot_sc:9.3%}")

    for a_null in (2.0, 3.0):
        print(f"\n{'='*110}\n=== a_null = {a_null:.0f}  (screening ceiling od = {T.od_of_a(a_null):.4f}) ===")
        for alpha in alphas:
            nmin = min_n_rejectable(a_null, alpha)
            print(f"\n  alpha = {alpha:.3g}  (smallest rejectable n = {nmin})")
            print(f"    {'sample':10s} {'sc-elig':>8s} | " +
                  " ".join(f"{'od='+format(od,'.4g'):>16s}" for od in ods))
            print(f"    {'':10s} {'':8s} | " + " ".join(f"{'bias(both) q':>16s}" for od in ods))
            for name in D:
                s, t, w, k, ex, sc, nc = views[name]
                u, c = size_hist(nc)
                cells = []
                elig = None
                for od in ods:
                    r = pooled_bias(u, c, a_null, alpha, od)
                    elig = r["frac_eligible_scale"]
                    cells.append(f"{r['od_hat']/od:8.6f} {r['q_scale_weighted']:6.1e}")
                print(f"    {name:10s} {elig:8.3%} | " + " ".join(cells))


if __name__ == "__main__":
    main()
