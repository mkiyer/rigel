"""P-2 RESIDUAL, PART C — diff the pinned tree against the unpinned one, node by node.

Part A falsified the "off-simplex claim" reading by looking at ONE tree.  The residual is a DIFFERENCE, so
attribute it as one: which nodes got worse when `_pin_v` stopped rewriting the delivered message, in which
direction, of which class, and what changed in the message packet they received.

Uses the bench's own node mask (`isfinite(oracle) & mass > 0`, no solvable gate) so the numbers reconcile
with `pass0_oracle_bench.py`'s headline.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_c_diff.py [--refit 0|1]
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
    A = pickle.loads(Path(f"scratchpad/p2r_dump_pin_r{a.refit}.pkl").read_bytes())   # pinned (pre-P-2)
    B = pickle.loads(Path(f"scratchpad/p2r_dump_p2_r{a.refit}.pkl").read_bytes())    # unpinned (P-2)

    print(f"=== C1. the mwae move, bench mask, per condition   [refit={a.refit}] ===")
    print(f"{'condition':50s} {'pin':>8s} {'P-2':>8s} {'delta':>8s} {'worse mass':>11s} "
          f"{'of which':>9s} {'of which':>9s}")
    print(f"{'':50s} {'':>8s} {'':>8s} {'':>8s} {'share':>11s} {'UNDER':>9s} {'OVER':>9s}")
    for c in A:
        p, q = A[c], B[c]
        m = np.isfinite(p["fo"]) & (p["mass"] > _EPS)
        w = p["mass"][m]
        ea = np.abs(p["fg"] - p["fo"])[m]
        eb = np.abs(q["fg"] - q["fo"])[m]
        d = (eb - ea) * w
        worse = d > 0
        # direction of the SOLVED move on the nodes that got worse
        dz = (q["fg"] - p["fg"])[m]
        und = worse & (dz < 0)   # solved f_g went DOWN  -> gDNA under-called harder
        ovr = worse & (dz > 0)   # solved f_g went UP    -> gDNA over-called harder
        print(f"{c[5:]:50s} {np.average(ea, weights=w):8.4f} {np.average(eb, weights=w):8.4f} "
              f"{np.average(eb - ea, weights=w):+8.4f} {d[worse].sum() / max(w.sum(), 1):11.4f} "
              f"{d[und].sum() / max(d[worse].sum(), _EPS):8.1%} "
              f"{d[ovr].sum() / max(d[worse].sum(), _EPS):8.1%}")

    print(f"\n\n=== C2. the REGRESSING stratum, by node class   [refit={a.refit}] ===")
    reg = [c for c in A if "ss_0.50" in c and "capture_off" in c and not c[5:].startswith("none_")]
    print(f"{'class':12s} {'n':>6s} {'n worse':>8s} {'net delta':>11s} {'mean fg pin':>12s} "
          f"{'mean fg P-2':>12s} {'mean oracle':>12s}")
    for ci, cn in enumerate(CLASSES):
        tot = np.zeros(4)
        fga, fgb, fob, ws = [], [], [], []
        nw = 0
        for c in reg:
            p, q = A[c], B[c]
            m = np.isfinite(p["fo"]) & (p["mass"] > _EPS) & (p["cls"] == ci)
            if not m.any():
                continue
            w = p["mass"][m]
            d = (np.abs(q["fg"] - q["fo"])[m] - np.abs(p["fg"] - p["fo"])[m]) * w
            tot += np.array([m.sum(), w.sum(), d.sum(), (d > 0).sum()])
            nw += int((d > 0).sum())
            fga.append(p["fg"][m]); fgb.append(q["fg"][m]); fob.append(p["fo"][m]); ws.append(w)
        if tot[0] == 0:
            continue
        fga, fgb, fob, ws = (np.concatenate(x) for x in (fga, fgb, fob, ws))
        print(f"{cn:12s} {int(tot[0]):6d} {nw:8d} {tot[2] / max(tot[1], 1):+11.4f} "
              f"{np.average(fga, weights=ws):12.3f} {np.average(fgb, weights=ws):12.3f} "
              f"{np.average(fob, weights=ws):12.3f}")

    print(f"\n\n=== C3. THE WORSE NODES: what claim did they receive?   [refit={a.refit}] ===")
    print("   partial := the incoming packet supplies RNA precision but NO gDNA precision.\n")
    print(f"{'population':34s} {'n':>7s} {'wmass':>12s} {'%partial':>9s} {'d(fg)':>8s} "
          f"{'oracle':>8s} {'fg pin':>8s} {'fg P-2':>8s}")
    for tag, want_worse in (("WORSE (P-2 hurt)", True), ("BETTER/flat", False)):
        nn = 0
        acc = {k: [] for k in ("w", "part", "dz", "fo", "fa", "fb")}
        for c in reg:
            p, q = A[c], B[c]
            m = np.isfinite(p["fo"]) & (p["mass"] > _EPS)
            d = np.abs(q["fg"] - q["fo"]) - np.abs(p["fg"] - p["fo"])
            sel = m & ((d > 1e-9) if want_worse else (d <= 1e-9))
            nn += int(sel.sum())
            acc["w"].append(p["mass"][sel])
            acc["part"].append(((q["prec_p"] + q["prec_n"] > 0) & (q["prec_g"] <= 0))[sel])
            acc["dz"].append((q["fg"] - p["fg"])[sel])
            acc["fo"].append(p["fo"][sel]); acc["fa"].append(p["fg"][sel]); acc["fb"].append(q["fg"][sel])
        A_ = {k: np.concatenate(v) for k, v in acc.items()}
        w = A_["w"]
        print(f"{tag:34s} {nn:7d} {w.sum():12,.0f} {np.average(A_['part'], weights=w):8.1%} "
              f"{np.average(A_['dz'], weights=w):+8.3f} {np.average(A_['fo'], weights=w):8.3f} "
              f"{np.average(A_['fa'], weights=w):8.3f} {np.average(A_['fb'], weights=w):8.3f}")

    print(f"\n\n=== C4. TOP 12 nodes by regression mass (all regressing conditions)   [refit={a.refit}] ===")
    rows = []
    for c in reg:
        p, q = A[c], B[c]
        m = np.isfinite(p["fo"]) & (p["mass"] > _EPS)
        d = (np.abs(q["fg"] - q["fo"]) - np.abs(p["fg"] - p["fo"])) * p["mass"]
        for i in np.argsort(-np.where(m, d, -np.inf))[:12]:
            i = int(i)
            rows.append((d[i], c, i, CLASSES[int(p["cls"][i])], p["mass"][i], p["fo"][i],
                         p["fg"][i], q["fg"][i], q["mode_g"][i], q["prec_g"][i],
                         q["mode_p"][i], q["prec_p"][i], q["mode_n"][i], q["prec_n"][i],
                         p["mode_p"][i], p["prec_p"][i], p["og"][i] / max(p["rho0"][i], _EPS)))
    rows.sort(key=lambda r: -r[0])
    print(f"{'dmass':>10s} {'node':>6s} {'cls':9s} {'mass':>10s} {'orc':>6s} {'pin':>6s} {'P-2':>6s} "
          f"{'e^moG':>7s} {'pG':>7s} {'e^moP':>7s} {'pP':>7s} {'e^moN':>7s} {'pN':>7s} "
          f"{'e^moP(pin)':>10s} {'fg_own':>7s}  cond")
    for r in rows[:12]:
        print(f"{r[0]:10,.0f} {r[2]:6d} {r[3]:9s} {r[4]:10,.0f} {r[5]:6.3f} {r[6]:6.3f} {r[7]:6.3f} "
              f"{np.exp(r[8]):7.3f} {r[9]:7.3g} {np.exp(r[10]):7.3f} {r[11]:7.3g} "
              f"{np.exp(r[12]):7.3f} {r[13]:7.3g} {np.exp(r[14]):10.3f} {r[16]:7.3f}  {r[1][5:]}")


if __name__ == "__main__":
    main()
