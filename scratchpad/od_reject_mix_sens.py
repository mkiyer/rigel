"""Sensitivity of the two-component EM to its initialisation and to od_r."""

from __future__ import annotations

import pickle
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from od_reject_mix import em_mixture  # noqa: E402

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"


def main():
    real = pickle.load(open(REAL, "rb"))
    syn = pickle.load(open(SYNTH, "rb"))
    cases = [("vcap", real["vcap"]["gdna"]), ("LBX0588", real["LBX0588"]["gdna"])]
    cases += [
        (k[9:29], v["gdna"])
        for k, v in syn.items()
        if "present" in k
    ]
    print(f"{'case':24s} {'od_r':>6} {'init od0':>9} {'init pi0':>9} {'-> od_g':>9} {'1-pi':>9}")
    for nm, (s, t, w, kap, _) in cases:
        for od_r in (0.01, 0.05, 0.15):
            for od0, pi0 in ((0.001, 0.99), (0.05, 0.9), (0.19, 0.5)):
                g, pi, r, n, k = em_mixture(s, t, kap, od_r=od_r, od0=od0, pi0=pi0)
                print(f"{nm:24s} {od_r:6.2f} {od0:9.3f} {pi0:9.2f} {g:9.4f} {1 - pi:9.3%}")
        print()


if __name__ == "__main__":
    main()
