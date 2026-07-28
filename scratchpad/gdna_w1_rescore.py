"""W1 RE-SCORE — the substrate table, against the CORRECTED reference (N5).

W1's table was decided by EMD against `oracle_landscape` at its default, which N5 showed is a COMB above
log10 rho = 0 (roughness 40.7 where a smooth bump is 2-4): it renders the TRUTH with a Poisson MEASUREMENT
kernel whose width collapses as rho^(-1/2). EMD to a comb rewards a fit that is also a comb, so every
substrate verdict taken on it is suspect until re-scored. This is that re-score, and it is the last
EMD-decided result still resting on the old instrument.

Two changes from the original table, both deliberate and both stated per column so neither hides the other:
  * the BASE substrate is now REGIONS ONLY (owner decision, 2026-07-27), so `+ boundaries` appears as an arm
    rather than as the default;
  * every arm is scored against BOTH references, so the effect of correcting the instrument is visible
    rather than assumed.

EMD is reported because it is what the original table used, but N5 measured it as nearly flat in the kernel
and monotone in smoothing — so the SHAPE columns (enriched mass and roughness) carry the verdict.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w1_rescore.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n5_roughness import region_substrate, roughness  # noqa: E402

KNN = 0.5
_EPS = 1e-12


def _bnd_live(s, mk):
    return mk["boundary"] & (s["eff"] > 1e-9) & (s["mass"] > _EPS)


def arms(s, mk):
    """(name, selection, weight-kind). BASE = regions only; each arm changes exactly one decision."""
    base = region_substrate(s, mk)
    live_r = (s["eff"] > 1e-9) & s["is_region"] & (s["mass"] > _EPS)
    yield "BASE (regions only)", base, "w"
    yield "+ AMBIG in", base | (live_r & mk["ambig"]), "w"
    yield "+ boundaries in", base | (_bnd_live(s, mk) & (mk["single"] | mk["gonly"])), "w"
    yield "- struct_zero anchor", base & ~mk["struct_zero"], "w"
    yield "- gonly (intergenic)", base & ~mk["gonly"], "w"
    yield "FLAT weights", base, "flat"
    yield "precision CUTOFF", base, "cut"


def weights_for(s, sel, mk, kind):
    if kind == "w":
        return L.recipe_weights(s, sel, mk)
    if kind == "flat":
        return np.ones(int(sel.sum()))
    v = np.maximum(np.nan_to_num(s["var"][sel], nan=1e9, posinf=1e9), 0.0)
    return (v <= np.median(v)).astype(float)


def run(suite, keep, title):
    scen = [s for s in L.load_scenarios(suite) if keep(s)]
    acc = {}
    for s in scen:
        comb = L.oracle_landscape(s)                      # the OLD instrument
        corr = L.oracle_landscape(s, knn_scale=KNN)       # the CORRECTED instrument
        if comb is None or corr is None:
            continue
        split = L.two_component(corr)["split"]
        hi = L.GRID > split
        ref_mass = float(corr[hi].sum())
        mk = L.masks(s)
        for name, sel, wk in arms(s, mk):
            if not sel.any():
                continue
            P = L.recipe(s, sel=sel, w=weights_for(s, sel, mk, wk), knn_scale=KNN)
            if P is None:
                continue
            acc.setdefault(name, []).append((
                L.emd(P, comb), L.emd(P, corr),
                float(P[hi].sum()) / max(ref_mass, 1e-9), roughness(P, 1.0, 2.5), int(sel.sum())))
    if not acc:
        return
    print(f"\n=== {title} (n={len(scen)}) ===")
    print(f"{'arm':24s} {'n_train':>8s} {'EMD/comb':>9s} {'d':>8s} | {'EMD/CORRECTED':>14s} {'d':>8s} | "
          f"{'enr mass':>9s} {'d':>8s} | {'rough>+1':>9s}")
    b = np.array(acc["BASE (regions only)"], float).mean(0)
    for name, rs in acc.items():
        a = np.array(rs, float).mean(0)
        print(f"{name:24s} {a[4]:8.0f} {a[0]:9.4f} {a[0] - b[0]:+8.4f} | {a[1]:14.4f} {a[1] - b[1]:+8.4f} | "
              f"{a[2]:9.3f} {a[2] - b[2]:+8.3f} | {a[3]:9.1f}")


def main():
    for su in ("ambig", "quick"):
        run(su, lambda s: True, f"{su}: ALL conditions")
        run(su, lambda s: s["group"][1] != "none" and s["group"][0] in ("ON", "VSTRONG"),
            f"{su}: gDNA-bearing capture-ON/VSTRONG — where the enriched mode exists")
        run(su, lambda s: s["group"][1] == "none", f"{su}: ZERO-gDNA — the fabrication guard")
        run(su, lambda s: s["group"][1] in ("gdna1", "gdna5"),
            f"{su}: LOW gDNA (gdna1 / gdna5) — the hardest regime")
    print("\n  'enr mass' is the fit's enriched mass as a fraction of the CORRECTED reference's "
          "(1.00 = right).")
    print("  'rough>+1' = TV of log P per decade above +1; a smooth unimodal bump scores 2-4.")


if __name__ == "__main__":
    main()
