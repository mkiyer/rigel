"""W1 — THE SUBSTRATE STUDY: which nodes should train the gDNA prior?

Approached from the taxonomy, not from a threshold (production plan §2.2). Four DISTINCT reasons a node
might be treated differently, which the current code conflates:

    circularity     -> STRUCTURAL exclusion   (AMBIG: the two-root ambiguity the prior exists to resolve)
    identifiability -> STRUCTURAL inclusion   (intergenic / zero-count structural = the depleted anchor)
    precision       -> CONTINUOUS weight, never a cutoff
    geometry        -> matched reference      (boundaries CROSS, regions CONTAIN)

The headline question for this run: the `tau_prior` boundary gate is 3 bare constants and a cliff, it
clips to its lower bound on 19/32 conditions, and it discards boundaries carrying 12-16% of boundary
weight. **What does simply removing it cost?**

Boundaries are held INCLUDED throughout (owner, 2026-07-27); include/exclude is a later ablation.

Every arm uses the SAME kernel and the SAME weights unless the arm is explicitly about weights, so a
substrate result is never confounded with a weighting change.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w1_substrate.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402

MODE_TOL = 0.30
ENR_THR = 0.0


# ---------------- substrate variants: (s, mk) -> boolean node mask ----------------
def _region_core(s, mk):
    """Region half, unchanged across arms: non-AMBIG live regions + the zero-count structural anchor."""
    return (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]


def _bnd_live(s, mk):
    return mk["boundary"] & (s["eff"] > 1e-9) & (s["mass"] > 1e-12)


def sub_nogate(s, mk):
    """THE DEFAULT (landed 2026-07-27): every live single/gonly boundary is admitted and the continuous
    weight decides how much it counts. The `tau_prior` admission cliff it replaced is gone — measured FREE
    to remove (EMD −0.0034 ambig / −0.0255 quick), so that arm no longer exists to compare against."""
    return L.recipe_substrate(s, mk)


def sub_nogate_ambig(s, mk):
    """+ AMBIG admitted — tests the CIRCULARITY exclusion (expected to be the one that must stay)."""
    live_r = (s["eff"] > 1e-9) & s["is_region"] & (s["mass"] > 1e-12)
    return sub_nogate(s, mk) | (live_r & mk["ambig"]) | (_bnd_live(s, mk) & mk["ambig"])


def sub_nogate_no_structzero(s, mk):
    """- the zero-count structural anchor: tests whether 'gDNA is absent here' is load-bearing."""
    return sub_nogate(s, mk) & ~mk["struct_zero"]


def sub_nogate_no_gonly(s, mk):
    """- structural-gDNA (intergenic) nodes: tests the IDENTIFIABILITY inclusion."""
    return sub_nogate(s, mk) & ~mk["gonly"]


ARMS = [
    ("DEFAULT (gate removed)", sub_nogate, "w"),
    ("+ AMBIG in", sub_nogate_ambig, "w"),
    ("- struct_zero", sub_nogate_no_structzero, "w"),
    ("- gonly", sub_nogate_no_gonly, "w"),
    ("FLAT weights", sub_nogate, "flat"),
    ("precision CUTOFF", sub_nogate, "cut"),
]


def weights_for(s, sel, mk, kind):
    if kind == "w":
        return None                                   # the production continuous weight
    if kind == "flat":
        return np.ones(int(sel.sum()))                # precision ignored entirely
    # a hard median-variance cutoff — the bad-practice control the gate is an instance of
    v = np.maximum(np.nan_to_num(s["var"][sel], nan=1e9, posinf=1e9), 0.0)
    return (v <= np.median(v)).astype(float)


def score(P, s):
    """Against the production target (region oracle) plus the mode-structure metrics."""
    orc = L.oracle_landscape(s)
    if P is None or orc is None:
        return None
    hi = L.GRID > ENR_THR
    om, fm = float(orc[hi].sum()), float(P[hi].sum())
    mo, mf = L.modes(orc), L.modes(P)
    rec = np.mean([any(abs(a - b) <= MODE_TOL for b, _ in mf) for a, _ in mo]) if mo else np.nan
    spu = np.mean([not any(abs(b - a) <= MODE_TOL for a, _ in mo) for b, _ in mf]) if mf else np.nan
    return dict(emd=L.emd(P, orc), recall=rec, spurious=spu,
                enr=(fm / om if om > 1e-4 else np.nan))


def main():
    print("W1 SUBSTRATE STUDY — boundaries INCLUDED throughout; kernel + weights fixed unless stated\n")
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        print(f"=== {su} ===")
        print(f"{'arm':28s} {'n_train':>8s} {'EMD':>7s} {'d(EMD)':>8s} {'recall':>7s} {'spurious':>9s}"
              f" {'enr mass':>9s}")
        base = None
        for name, fn, wkind in ARMS:
            rows, ntr = [], []
            for s in scen:
                mk = L.masks(s)
                sel = fn(s, mk)
                if not sel.any():
                    continue
                P = L.recipe(s, sel=sel, w=weights_for(s, sel, mk, wkind))
                r = score(P, s)
                if r:
                    rows.append((r["emd"], r["recall"], r["spurious"], r["enr"]))
                    ntr.append(int(sel.sum()))
            if not rows:
                continue
            a = np.array(rows, dtype=float)
            emd = float(np.nanmean(a[:, 0]))
            if base is None:
                base = emd
            print(f"{name:28s} {int(np.mean(ntr)):8d} {emd:7.3f} {emd - base:+8.4f} "
                  f"{np.nanmean(a[:, 1]):7.2f} {np.nanmean(a[:, 2]):9.2f} {np.nanmean(a[:, 3]):9.2f}")
        print()
    print("d(EMD) is vs the 'today' arm; + = worse. recall/spurious are the mode-structure metrics")
    print("(recall want 1.00, spurious want 0.00). enr mass = fit/oracle mass above log10 rho_g = 0.")


if __name__ == "__main__":
    main()
