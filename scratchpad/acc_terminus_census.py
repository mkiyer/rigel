"""How much of the human transcriptome sits inside the zone where the admissible-start
effective length departs from |A| — i.e. within one fragment length of a transcript
terminus, measured in TRANSCRIPT coordinates (the constraint is a transcript-length
constraint, not a genomic one).

Reports, per transcript, the share of exonic positions whose remaining transcript
length is short enough that F(remaining) < 1, weighted by exonic bp.
"""

import numpy as np
import pandas as pd

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"


def cdf(mu, sd, n=1200):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return np.cumsum(p / p.sum())


iv = pd.read_feather(D + "intervals.feather")
tx = pd.read_feather(D + "transcripts.feather")
ex = iv[iv["interval_type"] == 0]
real = set(map(int, tx.index[(~tx["is_synthetic"].astype(bool)) & (~tx["is_nrna"].astype(bool))]))
lens = ex.groupby("t_index")["end"].sum() - ex.groupby("t_index")["start"].sum()
lens = lens[[int(i) in real for i in lens.index]].to_numpy(np.int64)

print(f"{lens.size:,} real transcripts; median transcript length {int(np.median(lens))} bp\n")
print(f"{'FL pmf':22s} {'exonic bp with':>16s} {'mean E_r/|A|':>13s} {'transcripts':>12s}")
print(f"{'':22s} {'F(rem) < 0.99':>16s} {'over-all exons':>13s} {'wholly in zone':>12s}")
print("-" * 68)
for mu, sd in ((60, 15), (180, 60), (200, 50), (300, 90)):
    F = cdf(mu, sd)
    tot_bp = 0.0
    zone_bp = 0.0
    eff_sum = 0.0
    whole = 0
    for L in lens:
        rem = np.arange(L, 0, -1)  # remaining transcript length from each position
        f = F[np.clip(rem, 0, F.size - 1)]
        tot_bp += L
        zone_bp += float((f < 0.99).sum())
        eff_sum += float(f.sum())
        if f[0] < 0.99:
            whole += 1
    print(f"mu={mu:3d} sd={sd:3d}          {100*zone_bp/tot_bp:15.1f}% "
          f"{eff_sum/tot_bp:13.3f} {100*whole/lens.size:11.1f}%")
