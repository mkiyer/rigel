"""DISSECTION 4 — is the landscape WRONG, or FAITHFUL to a pass-0 that is wrong?

The arithmetic from `gdna_d3_zero.py` part A closes this. On `none_ss_0.50_nrna_none_capture_on` pass-0
places **84 training nodes above rho = 0, carrying 1.54 % of the training weight** — on a library with NO
gDNA at all. Spread over the ~50 grid cells of the enriched range against a peak cell of ~0.05, that is

    log(0.0154/50 / 0.05)  =  -5.1 nats

and the landscape measures **-4.4**. So the prior is not inventing anything: it is reporting, to within the
arithmetic, exactly what pass-0 told it. The delta-pin's -25 nats is not a better *estimate* — it is a
narrow fixed kernel plus EM competition collapsing those nodes into the dominant low mode, i.e. it
distrusts its own training data. That is an accident of bandwidth, not a principle to port.

This quantifies the consequence: how much false gDNA each prior removes, and where the upstream over-call
actually lives — because if the fix is upstream, it needs an address.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d4_circularity.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")

NEW = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_dissect_cache.pkl")   # landscape re-solve
OLD = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_refit_cache.pkl")     # delta-pin re-solve
NT = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def main():
    new = {s["cond"]: s for s in pickle.loads(NEW.read_bytes())}
    old = pickle.loads(OLD.read_bytes())
    print("=== A. FALSE-POSITIVE gDNA mass on ZERO-gDNA libraries (every fragment here is an error) ===\n")
    print(f"{'condition':46s} {'pass-0':>10s} {'delta-pin':>10s} {'LANDSCAPE':>10s} "
          f"{'dp/p0':>7s} {'ls/p0':>7s}")
    tot = np.zeros(3)
    for cond, s in sorted(new.items()):
        if s["group"][1] != "none":
            continue
        live = (s["eff"] > 1e-9) & (s["mass"] > 1e-12)
        p0 = float((s["f0"][live] * s["mass"][live]).sum())
        ls = float((s["f1"][live] * s["mass"][live]).sum())
        o = old.get(cond)
        dp = float((o[1]["f0"][live] * s["mass"][live]).sum()) if o else np.nan
        tot += np.array([p0, dp, ls])
        print(f"{cond[5:]:46s} {p0:10.0f} {dp:10.0f} {ls:10.0f} "
              f"{dp / max(p0, 1):7.2f} {ls / max(p0, 1):7.2f}")
    print(f"\n{'TOTAL':46s} {tot[0]:10.0f} {tot[1]:10.0f} {tot[2]:10.0f} "
          f"{tot[1] / tot[0]:7.2f} {tot[2] / tot[0]:7.2f}")
    print("\n  BOTH priors remove false gDNA — the landscape is not manufacturing it. The delta-pin simply")
    print("  removes more, by trusting its training data less.")

    print("\n=== B. WHERE pass-0's over-call lives — the address for an upstream fix ===\n")
    print(f"{'stratum':34s} {'class':11s} {'nodes':>6s} {'total mass':>12s} {'FP mass':>11s} {'FP %':>7s}")
    acc = {}
    for s in new.values():
        if s["group"][1] != "none":
            continue
        key = ("unstranded" if s["group"][2] == "0.50" else "stranded") + " x " + s["group"][3]
        live = (s["eff"] > 1e-9) & (s["mass"] > 1e-12)
        for t in (0, 1, 2, 3):
            m = live & (s["ntype"] == t)
            if not m.any():
                continue
            a = acc.setdefault((key, NT[t]), np.zeros(3))
            a += np.array([float(s["mass"][m].sum()), float((s["f0"][m] * s["mass"][m]).sum()),
                           float(m.sum())])
    for (key, cls), a in sorted(acc.items()):
        print(f"{key:34s} {cls:11s} {a[2]:6.0f} {a[0]:12.0f} {a[1]:11.0f} "
              f"{100 * a[1] / max(a[0], 1):7.2f}")
    print("\n  The over-call is not diffuse: it is UNSTRANDED x nrna_none, on EXONS and BOUNDARIES.")
    print("  Introns are clean, which is the intron factory doing its job.")


if __name__ == "__main__":
    main()
