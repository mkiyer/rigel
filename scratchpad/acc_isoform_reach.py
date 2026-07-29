"""Q/D2 — how much do isoforms disagree about `reach`?

`reach_M(s)` = exonic bases remaining downstream of genomic position `s` within a
transcript. A position usually belongs to several isoforms with different reaches, and
the honest opportunity is sum_t theta_t * F(reach_t) with theta unknown.

Candidate conventions:
  (a) MAXIMAL   F(max_t reach_t)      -- "open if ANY isoform supports it"   (recommended)
  (c) MEAN      mean_t F(reach_t)     -- assumes equal expression

If (a) and (c) agree on nearly all exonic mass, the choice is immaterial and D2 closes.
Measured on the real human index, sampling genes with >= 2 real transcripts.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"
FL_MU, FL_SD = 200.0, 50.0
N_GENES = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
RNG = np.random.default_rng(11)


def cdf(mu, sd, n=1600):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return np.cumsum(p / p.sum())


F = cdf(FL_MU, FL_SD)


def main():
    iv = pd.read_feather(D + "intervals.feather")
    tx = pd.read_feather(D + "transcripts.feather")
    real = (~tx["is_synthetic"].astype(bool)) & (~tx["is_nrna"].astype(bool))
    tx = tx[real]
    ex = iv[iv["interval_type"] == 0]
    ex = ex[ex["t_index"].isin(set(tx.index.astype(int)))]

    gene_of = dict(zip(tx.index.astype(int), tx["g_index"].astype(int)))
    strand_of = dict(zip(tx.index.astype(int), tx["strand"].astype(int)))
    ex = ex.sort_values(["t_index", "start"])
    by_t = {t: g[["start", "end"]].to_numpy(np.int64) for t, g in ex.groupby("t_index")}

    genes: dict[int, list[int]] = {}
    for t in by_t:
        genes.setdefault(gene_of[int(t)], []).append(int(t))
    multi = [g for g, ts in genes.items() if len(ts) >= 2]
    print(f"{len(by_t):,} real transcripts in {len(genes):,} genes; "
          f"{len(multi):,} genes have >= 2 isoforms")
    pick = RNG.choice(len(multi), size=min(N_GENES, len(multi)), replace=False)
    sample = [multi[i] for i in pick]

    tot_bp = 0
    zone_bp = 0
    agree = 0
    diffs = []
    for g in sample:
        acc: dict[int, list[float]] = {}
        for t in genes[g]:
            blocks = by_t[t]
            pos = np.concatenate([np.arange(a, b) for a, b in blocks])
            L = pos.size
            if L == 0:
                continue
            # transcript coordinate u: 5'->3'.  + strand ascends genomically, - descends.
            u = np.arange(L) if strand_of[t] == 1 else np.arange(L)[::-1]
            reach = L - u                      # exonic bases remaining, inclusive
            fv = F[np.clip(reach, 0, F.size - 1)]
            for p_, f_ in zip(pos, fv):
                acc.setdefault(int(p_), []).append(float(f_))
        for _p, vals in acc.items():
            tot_bp += 1
            v = np.asarray(vals)
            mx, mn = v.max(), v.mean()
            if mx < 0.99:                       # this position is in SOMEONE's taper zone
                zone_bp += 1
            if abs(mx - mn) < 1e-9:
                agree += 1
            else:
                diffs.append(mx - mn)

    d = np.asarray(diffs) if diffs else np.zeros(1)
    print(f"\nsampled {len(sample):,} multi-isoform genes -> {tot_bp:,} exonic positions")
    print(f"  positions where MAXIMAL == MEAN exactly (all isoforms agree): "
          f"{100*agree/max(tot_bp,1):.1f} %")
    print(f"  positions in SOMEONE's taper zone (max F < 0.99):            "
          f"{100*zone_bp/max(tot_bp,1):.1f} %")
    print("\n  where they disagree, |F_max - F_mean| :")
    for q in (50, 75, 90, 95, 99):
        print(f"     p{q:<3d} {np.percentile(d, q):.4f}")
    print(f"     mean over ALL exonic positions: {d.sum()/max(tot_bp,1):.4f}")
    print(f"\n  => using MAXIMAL instead of MEAN changes the mature opportunity by "
          f"{100*d.sum()/max(tot_bp,1):.2f} % of exonic bp on average")


if __name__ == "__main__":
    main()
