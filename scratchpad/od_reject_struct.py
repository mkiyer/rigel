"""Seed-size / information STRUCTURE: real vs synthetic. Where does the evidence live?"""

from __future__ import annotations

import pickle
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from od_reject_frames import terms  # noqa: E402

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"


def row(tag, sense, total, weight, kappa):
    excess, scale, nc, _ = terms(sense, total, weight, kappa)
    ok = (total > 0) & (nc > 0)
    excess, scale, nc = excess[ok], scale[ok], nc[ok]
    P = np.sum(scale)
    o = np.argsort(-scale)
    top1 = scale[o[0]] / P
    top10 = np.sum(scale[o[:10]]) / P
    # how many seeds hold 50% / 90% of the pairs?
    cw = np.cumsum(scale[o]) / P
    n50 = int(np.searchsorted(cw, 0.5) + 1)
    n90 = int(np.searchsorted(cw, 0.9) + 1)
    q = np.percentile(nc, [50, 90, 99, 99.9])
    print(
        f"{tag:52s} m={nc.size:7d} n_c q50/q90/q99/q99.9="
        f"{q[0]:7.0f}{q[1]:8.0f}{q[2]:8.0f}{q[3]:9.0f} max={nc.max():9.0f} | "
        f"pairs: top1={top1:7.2%} top10={top10:7.2%} n(50%)={n50:6d} n(90%)={n90:7d}"
    )


def main():
    print("SEED-SIZE AND PAIR-CONCENTRATION STRUCTURE\n")
    syn = pickle.load(open(SYNTH, "rb"))
    for k, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        row("SYNTH " + k.replace("gdna_gdna", "g").replace("_capture", "_cap"), s, t, w, kap)
    print()
    real = pickle.load(open(REAL, "rb"))
    for k in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[k]["gdna"]
        row("REAL  " + k, s, t, w, kap)


if __name__ == "__main__":
    main()
