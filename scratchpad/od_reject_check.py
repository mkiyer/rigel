"""Sanity-check the stratified readout: are the seed masses integral, and is gdna_weight ~ 1?"""

from __future__ import annotations

import pickle
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"


def check(tag, sense, total, weight, kappa):
    ok = total > 0
    s, t, w = sense[ok], total[ok], weight[ok]
    integral_t = np.mean(np.abs(t - np.rint(t)) < 1e-6)
    integral_s = np.mean(np.abs(s - np.rint(s)) < 1e-6)
    # EXACT n=2, pure-gDNA-weight seeds only: no rounding, no mixture
    m2 = (np.abs(t - 2.0) < 1e-9) & (w > 0.999)
    same = np.abs(s - 1.0) > 1e-9  # sense is 0 or 2 => same strand
    p_same = float(np.mean(same[m2])) if m2.any() else np.nan
    m3 = (np.abs(t - 3.0) < 1e-9) & (w > 0.999)
    # n=3: E[(K-1.5)^2] = 3/4 (1+2od)
    od3 = (4.0 * np.mean((s[m3] - 1.5) ** 2) / 3.0 - 1.0) / 2.0 if m3.any() else np.nan
    print(
        f"{tag:46s} w>=0.999: {np.mean(w > 0.999):6.2%}  median w={np.median(w):.4f}  "
        f"integral N/K: {integral_t:6.2%}/{integral_s:6.2%} | "
        f"EXACT n=2,w=1: m={int(m2.sum()):7d} P(same)={p_same:7.4f} od={2 * p_same - 1:8.4f}"
        f" | n=3 od={od3:8.4f} (m={int(m3.sum())})"
    )


def main():
    real = pickle.load(open(REAL, "rb"))
    for k in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[k]["gdna"]
        check("REAL  " + k, s, t, w, kap)
    print()
    syn = pickle.load(open(SYNTH, "rb"))
    for k, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        check("SYNTH " + k[:36], s, t, w, kap)


if __name__ == "__main__":
    main()
