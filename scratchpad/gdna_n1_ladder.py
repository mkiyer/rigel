"""N1 — the SUBSTRATE ladder: exactly which step drops the truly-enriched weight share.

`gdna_n1_massflow.py` reduced the enriched-mass deficit to ONE number: the share of training weight sitting
on truly-enriched nodes is 0.0612 against the oracle's 0.1474 (ratio 0.415), while PLACEMENT is 0.873.
This walks the substrate from the oracle's set to the production one, ONE change per rung, and reports the
truly-enriched share at each. Every rung is a count share except the last, which adds the weight.

    w_enr = sum_{i in S, true > split} w_i / sum_{i in S} w_i

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_ladder.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds, truth  # noqa: E402

_EPS = 1e-12


def rungs(s, mk):
    """(name, selection) from the oracle's substrate down to production's, one change per rung."""
    live = s["eff"] > 1e-9
    isr = mk["region"]
    havem = mk["havemass"]
    reg_prod = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
    bnd_prod = mk["boundary"] & live & havem & (mk["single"] | mk["gonly"])
    yield "0  ORACLE substrate: live REGION nodes", live & isr
    yield "1  - AMBIG regions (circularity)", live & isr & ~mk["ambig"]
    yield "2  - regions with no mass (+ struct_zero back)", reg_prod
    yield "3  + BOUNDARY nodes  = production substrate", reg_prod | bnd_prod


def main():
    ss = conds()
    print("=== N1: the substrate ladder — truly-enriched share of the training set ===")
    print("    (13 gDNA-bearing capture-ON/VSTRONG conditions; share = truly-enriched / all, at the "
          "oracle's split)\n")
    print(f"{'rung':46s} {'n_train':>8s} {'n_enr':>7s} {'share':>7s} {'vs oracle':>10s}")
    acc = {}
    for s in ss:
        split = L.two_component(L.oracle_landscape(s))["split"]
        tr = truth(s)
        mk = L.masks(s)
        for name, sel in rungs(s, mk):
            enr = sel & (tr > split)
            acc.setdefault(name, []).append((int(sel.sum()), int(enr.sum()),
                                             enr.sum() / max(sel.sum(), 1)))
        # the final rung, with the production weights
        sel = L.recipe_substrate(s, mk)
        w = L.recipe_weights(s, sel, mk)
        enr = (tr[sel] > split)
        acc.setdefault("4  + reliability WEIGHT = production fit", []).append(
            (int(sel.sum()), int(enr.sum()), float(w[enr].sum() / max(w.sum(), _EPS))))
    base = float(np.mean([r[2] for r in acc["0  ORACLE substrate: live REGION nodes"]]))
    for name, rs in acc.items():
        n = float(np.mean([r[0] for r in rs]))
        ne = float(np.mean([r[1] for r in rs]))
        sh = float(np.mean([r[2] for r in rs]))
        print(f"{name:46s} {n:8.0f} {ne:7.0f} {sh:7.4f} {sh / base:10.3f}")

    # --- who the boundary nodes are, since rung 3 is the big one --------------------------------------
    print("\n=== rung 3 detail: the BOUNDARY population that gets added ===")
    print(f"{'':46s} {'n':>8s} {'n_enr':>7s} {'share':>7s} {'med w':>7s}")
    for lbl, pick in (("regions in production substrate", "reg"), ("boundaries added", "bnd")):
        N, NE, W = [], [], []
        for s in ss:
            split = L.two_component(L.oracle_landscape(s))["split"]
            tr, mk = truth(s), L.masks(s)
            live, havem = s["eff"] > 1e-9, mk["havemass"]
            reg = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
            bnd = mk["boundary"] & live & havem & (mk["single"] | mk["gonly"])
            sel = reg if pick == "reg" else bnd
            full = L.recipe_substrate(s, mk)
            w = np.zeros(sel.shape)
            w[full] = L.recipe_weights(s, full, mk)
            N.append(int(sel.sum()))
            NE.append(int((sel & (tr > split)).sum()))
            W.append(float(np.median(w[sel])) if sel.any() else np.nan)
        print(f"{lbl:46s} {np.mean(N):8.0f} {np.mean(NE):7.0f} "
              f"{np.mean(NE) / max(np.mean(N), 1):7.4f} {np.nanmean(W):7.3f}")


if __name__ == "__main__":
    main()
