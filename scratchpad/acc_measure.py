"""Decisive measurements for the accumulator redesign, on the REAL human index.

M1  region-size distribution + the "tiny node" census (region_eff_length ~ 0)
M2  TSS/TES positions that fall strictly INSIDE a region (invisible to a
    boundary-keyed flag; would become NEW partition events under the redesign)
M3  the crossing census: what fraction of fragment start positions produce a
    fragment that crosses >=1, >=2, >=3 region boundaries at a given FL
"""

import sys

import numpy as np
import pandas as pd

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"
FL = int(sys.argv[1]) if len(sys.argv) > 1 else 200

reg = pd.read_feather(D + "regions.feather")
iv = pd.read_feather(D + "intervals.feather")
tx = pd.read_feather(D + "transcripts.feather")

print(f"=== index: {len(reg):,} regions, {len(tx):,} transcripts, FL = {FL} ===\n")

# ---------------------------------------------------------------- M1
L = reg["length"].to_numpy(np.int64)
print("M1  REGION SIZE DISTRIBUTION")
qs = [1, 5, 10, 25, 50, 75, 90, 99]
print("    percentiles bp:", {q: int(np.percentile(L, q)) for q in qs})
for thr in (FL, 2 * FL):
    n = int((L < thr).sum())
    print(f"    regions shorter than {thr:5d} bp: {n:7,} ({100*n/len(L):5.1f} %)"
          f"  covering {100*L[L<thr].sum()/L.sum():5.2f} % of bp")
# contained eff-length E[max(0, L - l + 1)] at a delta FL
eff_contained = np.maximum(L - FL + 1, 0)
zero = int((eff_contained <= 0).sum())
print(f"    CONTAINED eff-length == 0 (cannot hold a contained fragment):"
      f" {zero:,} regions ({100*zero/len(L):.1f} %)")
print(f"    ... those regions hold {100*L[eff_contained<=0].sum()/L.sum():.2f} % of all bp\n")

# ---------------------------------------------------------------- M2
# TSS/TES from real (non-synthetic) transcripts; exon intervals are type 0.
ex = iv[iv["interval_type"] == 0]
real = tx.index[(~tx["is_synthetic"].astype(bool)) & (~tx["is_nrna"].astype(bool))]
real_set = set(map(int, real))
g = ex.groupby("t_index").agg(ref=("ref", "first"), start=("start", "min"),
                              end=("end", "max"), strand=("strand", "first"))
g = g[[int(i) in real_set for i in g.index]]
print(f"M2  TSS/TES CENSUS  ({len(g):,} real transcripts)")

term = {}   # ref -> set of terminus positions
for ref, s, e in zip(g["ref"].to_numpy(), g["start"].to_numpy(), g["end"].to_numpy()):
    term.setdefault(str(ref), set()).update((int(s), int(e)))

# region interfaces per ref
iface = {ref: np.union1d(gr["start"].to_numpy(np.int64), gr["end"].to_numpy(np.int64))
         for ref, gr in reg.groupby("ref_name", sort=False)}

n_term = n_on_iface = n_interior = 0
interior_by_ref = {}
for ref, pos in term.items():
    arr = np.fromiter(pos, np.int64, len(pos))
    n_term += arr.size
    ifa = iface.get(ref)
    if ifa is None:
        n_interior += arr.size
        continue
    hit = np.isin(arr, ifa)
    n_on_iface += int(hit.sum())
    n_interior += int((~hit).sum())
    interior_by_ref[ref] = arr[~hit]
print(f"    distinct terminus positions           : {n_term:,}")
print(f"    ... that ARE a region interface today : {n_on_iface:,} ({100*n_on_iface/n_term:.1f} %)")
print(f"    ... STRICTLY INSIDE a region (invisible): {n_interior:,} ({100*n_interior/n_term:.1f} %)")
print(f"    -> making TSS/TES partition events adds {n_interior:,} regions "
      f"({100*n_interior/len(reg):+.1f} % on {len(reg):,})\n")

# ---------------------------------------------------------------- M3
print("M3  CROSSING CENSUS (delta FL, per-bp start positions)")
tot_bp = 0
cross = np.zeros(6, np.int64)  # index k = # of boundaries crossed, capped at 5
for ref, gr in reg.groupby("ref_name", sort=False):
    ends = np.sort(gr["end"].to_numpy(np.int64))
    starts = np.sort(gr["start"].to_numpy(np.int64))
    lo, hi = int(starts[0]), int(ends[-1])
    if hi - lo <= FL:
        continue
    # for each interior boundary position b, a fragment starting in [b-FL, b) crosses it.
    # count, per start position s in [lo, hi-FL), how many boundaries lie in (s, s+FL].
    b = ends[:-1]  # interior interfaces
    s = np.arange(lo, hi - FL, dtype=np.int64)
    # number of boundaries in (s, s+FL] = searchsorted(b, s+FL, 'right') - searchsorted(b, s, 'right')
    k = np.searchsorted(b, s + FL, side="right") - np.searchsorted(b, s, side="right")
    tot_bp += s.size
    kk = np.minimum(k, 5)
    cross += np.bincount(kk, minlength=6)
print(f"    start positions sampled: {tot_bp:,}")
for k in range(6):
    lbl = f"{k}" if k < 5 else ">=5"
    print(f"    crosses {lbl:>3} boundaries: {cross[k]:12,} ({100*cross[k]/tot_bp:5.2f} %)")
n_cross = tot_bp - cross[0]
print(f"    -> CROSSING at all : {n_cross:,} ({100*n_cross/tot_bp:.2f} %)")
print(f"    -> of the crossers, {100*(n_cross-cross[1])/max(n_cross,1):.2f} % cross MORE THAN ONE"
      f" (today: split fractionally across >=2 boundaries)")
