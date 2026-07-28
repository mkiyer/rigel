"""P-2 RESIDUAL, PART D — split P-2's move into the pin's TWO BRANCHES.

`pin_derivation.md` §4/§8 says the pin is **correct and BP-admissible for a COMPLETE claim** and that "the
entire defect is the partial-claim branch", where the budget `S` borrows the DESTINATION's own density for a
component the message did not supply.  P-2 removed the pin from the delivered message on BOTH branches.

So split every node by whether the pin's budget was CONTAMINATED — i.e. whether any component with no
message precision had a NONZERO own density to lend the budget:

    contaminated := Σ_c [ prec_c == 0 and own_c > 0 ]  > 0

(a structurally dead strand has own density 0, so it lends nothing and does NOT contaminate — the docstring's
own carve-out.)  If P-2's GAIN is on the contaminated branch and its LOSS on the clean branch, then gating
the pin on `clean` keeps every win and repays the residual.

PREDICTION (registered before running): the `gdna_none` and capture-ON gains are contaminated-branch; the
capture-OFF × gDNA-bearing loss is clean-branch.  FALSIFIER: if the loss is also on the contaminated branch,
the gate cannot recover it and the diagnosis is wrong.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_d_branch.py [--refit 0]
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np

_EPS = 1e-9
CLASSES = ("intergenic", "intron", "exon", "boundary")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refit", type=int, default=0)
    a = ap.parse_args()
    A = pickle.loads(Path(f"scratchpad/p2r_dump_pin_r{a.refit}.pkl").read_bytes())
    B = pickle.loads(Path(f"scratchpad/p2r_dump_p2_r{a.refit}.pkl").read_bytes())

    print(f"=== D1. P-2's per-condition move, SPLIT BY BRANCH   [refit={a.refit}] ===")
    print("   contaminated = the pin's budget borrowed the destination's own density for an")
    print("   unsupplied component (the BP violation).  clean = it did not (M and E only).\n")
    print(f"{'condition':50s} {'delta':>8s} | {'CLEAN':>8s} {'nodes':>7s} {'mass%':>6s} | "
          f"{'CONTAM':>8s} {'nodes':>7s} {'mass%':>6s}")
    for c in A:
        p, q = A[c], B[c]
        m = np.isfinite(p["fo"]) & (p["mass"] > _EPS)
        w = p["mass"]
        # the pin's own supplied-test, on the fused mode-fusion precisions
        contam = (((q["cpg"] <= 0) & (q["og"] > 0))
                  | ((q["cpp"] <= 0) & (q["op"] > 0))
                  | ((q["cpn"] <= 0) & (q["on"] > 0)))
        d = (np.abs(q["fg"] - q["fo"]) - np.abs(p["fg"] - p["fo"])) * w
        tw = w[m].sum()
        cl, ct = m & ~contam, m & contam
        print(f"{c[5:]:50s} {d[m].sum() / tw:+8.4f} | {d[cl].sum() / tw:+8.4f} {int(cl.sum()):7d} "
              f"{100 * w[cl].sum() / tw:5.1f}% | {d[ct].sum() / tw:+8.4f} {int(ct.sum()):7d} "
              f"{100 * w[ct].sum() / tw:5.1f}%")

    print(f"\n\n=== D2. the CLEAN branch: what does the pin's budget say there?   [refit={a.refit}] ===")
    print("   S/M = the delivered claim's implied fragment count over the node's observed mass.")
    print("   S/M > 1 ⇒ the level claims over-account, and the pin used to shrink them by M/S.\n")
    reg = [c for c in A if "ss_0.50" in c and "capture_off" in c and not c[5:].startswith("none_")]
    grp = [("REGRESSING unstr x capOFF x gDNA", reg),
           ("control    none_* capOFF", [c for c in A if c[5:].startswith("none_")]),
           ("control    capture ON", [c for c in A if "capture_on" in c and "ss_0.50" in c]),
           ("control    stranded capOFF", [c for c in A if "ss_0.99" in c])]
    print(f"{'group':34s} {'n clean':>8s} {'med S/M':>8s} {'p90 S/M':>8s} {'%S>M':>6s} "
          f"{'d(fg) P-2-pin':>14s} {'oracle':>7s} {'net delta':>10s}")
    for gname, cs in grp:
        sm, dz, fo, dd, ww = [], [], [], [], []
        for c in cs:
            p, q = A[c], B[c]
            m = np.isfinite(p["fo"]) & (p["mass"] > _EPS)
            contam = (((q["cpg"] <= 0) & (q["og"] > 0))
                      | ((q["cpp"] <= 0) & (q["op"] > 0))
                      | ((q["cpn"] <= 0) & (q["on"] > 0)))
            sel = m & ~contam & ((q["cpg"] + q["cpp"] + q["cpn"]) > 0)
            S = q["cg"] * q["E_g"] + (q["cp"] + q["cn"]) * q["E_r"]
            sm.append((S / np.maximum(q["mass"], _EPS))[sel])
            dz.append((q["fg"] - p["fg"])[sel])
            fo.append(q["fo"][sel])
            dd.append(((np.abs(q["fg"] - q["fo"]) - np.abs(p["fg"] - p["fo"])) * q["mass"])[sel])
            ww.append(q["mass"][sel])
        sm, dz, fo, dd, ww = (np.concatenate(x) for x in (sm, dz, fo, dd, ww))
        if not ww.size:
            continue
        print(f"{gname:34s} {ww.size:8d} {np.median(sm):8.3f} {np.quantile(sm, 0.9):8.3f} "
              f"{100 * np.average(sm > 1.0, weights=ww):5.1f}% {np.average(dz, weights=ww):+14.3f} "
              f"{np.average(fo, weights=ww):7.3f} {dd.sum() / ww.sum():+10.4f}")

    print(f"\n\n=== D3. CLEAN-branch nodes in the regressing stratum, by class   [refit={a.refit}] ===")
    print(f"{'class':12s} {'n':>7s} {'mass':>12s} {'med S/M':>8s} {'oracle':>7s} {'fg pin':>7s} "
          f"{'fg P-2':>7s} {'net delta':>10s}")
    for ci, cn in enumerate(CLASSES):
        sm, fo, fa, fb, dd, ww = [], [], [], [], [], []
        for c in reg:
            p, q = A[c], B[c]
            contam = (((q["cpg"] <= 0) & (q["og"] > 0))
                      | ((q["cpp"] <= 0) & (q["op"] > 0))
                      | ((q["cpn"] <= 0) & (q["on"] > 0)))
            sel = (np.isfinite(p["fo"]) & (p["mass"] > _EPS) & ~contam & (p["cls"] == ci)
                   & ((q["cpg"] + q["cpp"] + q["cpn"]) > 0))
            S = q["cg"] * q["E_g"] + (q["cp"] + q["cn"]) * q["E_r"]
            sm.append((S / np.maximum(q["mass"], _EPS))[sel])
            fo.append(q["fo"][sel]); fa.append(p["fg"][sel]); fb.append(q["fg"][sel])
            dd.append(((np.abs(q["fg"] - q["fo"]) - np.abs(p["fg"] - p["fo"])) * q["mass"])[sel])
            ww.append(q["mass"][sel])
        sm, fo, fa, fb, dd, ww = (np.concatenate(x) for x in (sm, fo, fa, fb, dd, ww))
        if not ww.size:
            continue
        tw = sum(A[c]["mass"][np.isfinite(A[c]["fo"]) & (A[c]["mass"] > _EPS)].sum() for c in reg)
        print(f"{cn:12s} {ww.size:7d} {ww.sum():12,.0f} {np.median(sm):8.3f} "
              f"{np.average(fo, weights=ww):7.3f} {np.average(fa, weights=ww):7.3f} "
              f"{np.average(fb, weights=ww):7.3f} {dd.sum() / tw:+10.4f}")


if __name__ == "__main__":
    main()
