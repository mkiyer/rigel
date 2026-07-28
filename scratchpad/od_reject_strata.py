"""STRATIFY the dispersion by seed size. Is the excess a property of BIG seeds (contamination)
or of ALL seeds (real dispersion)?  The n=2 stratum is an assumption-free readout:
under BetaBinom(2, 1/2, od),  P(both same strand) = (1 + od)/2   =>   od = 2*P(same) - 1.
"""

from __future__ import annotations

import pickle
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"
_EPS = 1e-12


def stratify(tag, sense, total, weight, kappa):
    nc = weight * total
    ni = np.rint(nc).astype(np.int64)
    ki = np.clip(np.rint(sense * np.where(total > 0, nc / np.maximum(total, _EPS), 1.0)), 0, ni)
    ki = ki.astype(np.int64)
    print(f"\n{tag}")
    print(
        f"   {'n stratum':>14} {'#seeds':>9} {'%pairs':>8} {'od_hat(n)':>10}"
        f"  {'(n=2: %same-strand)':>22}"
    )
    edges = [(2, 2), (3, 3), (4, 4), (5, 6), (7, 10), (11, 20), (21, 50), (51, 200),
             (201, 1000), (1001, 10**9)]
    allpairs = np.sum(ni * (ni - 1) / 2.0)
    for lo, hi in edges:
        m = (ni >= lo) & (ni <= hi)
        if not m.any():
            continue
        n, k = ni[m].astype(float), ki[m].astype(float)
        x2 = np.mean((k - n / 2.0) ** 2)
        # per-stratum MoM: E[(K-n/2)^2] = n/4 (1 + (n-1) od)  (use mean n within the stratum)
        nb = np.mean(n)
        od = (4.0 * x2 / nb - 1.0) / max(nb - 1.0, _EPS)
        pairs = np.sum(n * (n - 1) / 2.0) / allpairs
        extra = ""
        if lo == 2 and hi == 2:
            extra = f"{np.mean(k != 1):22.4%}"
        print(f"   {f'{lo}-{hi}':>14} {int(m.sum()):9d} {pairs:8.2%} {od:10.4f}  {extra:>22}")


def main():
    print("=" * 96)
    print("SIZE-STRATIFIED OVERDISPERSION  (flat => real dispersion; rising => size-linked contamination)")
    print("=" * 96)
    real = pickle.load(open(REAL, "rb"))
    for k in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[k]["gdna"]
        stratify("REAL  " + k, s, t, w, kap)
    syn = pickle.load(open(SYNTH, "rb"))
    for k, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        stratify("SYNTH " + k, s, t, w, kap)


if __name__ == "__main__":
    main()
