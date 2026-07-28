"""N2 (RE-RUN after W4) — the ITERATIVE AMBIG resolution, now judged against the LANDSCAPE.

    pass-0 -> prior #1 (AMBIG excluded) -> RE-SOLVE -> prior #2 (re-solved AMBIG INCLUDED)

is non-circular by construction: the AMBIG estimates prior #2 trains on were produced by a prior that never
saw them.  The first run (`gdna_n2_iterative.py`, 2026-07-27) CONFIRMED the premise (Q1: a re-solve cuts
AMBIG density error 42 %) but its verdict was **contaminated** — prior #1 was the shipped δ-pin
`DensityNPMLE`, the object this workstream retires, and it degrades the non-AMBIG nodes it touches.

THREE things differ here, and each is a correction the record demands:

1. **prior #1 is now the `GdnaLandscape`** (W4, in the tree). The belief cache was rebuilt at this tree, so
   both iterations come from the shipped object.
2. **the reference is rendered at POPULATION resolution** (`knn_scale=0.5`). N5: the default `knn_scale=0`
   makes `oracle_landscape` a COMB above `log10 ρ ≈ 0` (roughness 40.7), and EMD to a comb rewards a comb.
   Every W1/W2/W3 number taken before the fix was scored against it.
3. **the substrate is REGIONS ONLY** (owner decision D1, 2026-07-27), matching the shipped
   `_fit_gdna_hyperprior`. `gdna_explore_lib.recipe_substrate` still includes boundaries — it is the
   pre-decision version, deliberately not edited here so the other tools keep their meaning.

Scoring is at the **ORACLE's own split** (N1: scoring each landscape at its own split reads 91 % recovery
and hides the entire defect), and enriched mass is read alongside the ZERO-gDNA fabrication guard, because
broadening a kernel raises "enriched mass" for free.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n2b_iterative.py   (needs scratchpad/gdna_refit_cache.pkl)
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402

_EPS = 1e-12
_KNN = 0.5  # the one resolution law, applied to the fit AND the reference (N5)
REFIT = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_refit_cache.pkl")


def joined(suite="ambig"):
    rf = pickle.loads(REFIT.read_bytes())
    for s in L.load_scenarios(suite):
        r = rf.get(s["cond"])
        if r is None:
            continue
        assert r[0]["f0"].shape == s["f0"].shape and np.allclose(r[0]["f0"], s["f0"]), s["cond"]
        yield s, r


def production_substrate(s, mk):
    """Exactly what the shipped `_fit_gdna_hyperprior` selects: REGION nodes, AMBIG out, boundaries out,
    plus the zero-count structural anchor on non-exon regions."""
    return (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]


def q1():
    print("=== Q1: does the RE-SOLVE improve the AMBIG nodes' gDNA density?  [prior #1 = LANDSCAPE] ===")
    print("    mean |log10 rho_hat - log10 rho_true| over live AMBIG REGION nodes, in decades\n")
    print(f"{'stratum':22s} {'n':>7s} {'pass-0':>8s} {'re-solve':>9s} {'delta':>8s} {'better':>8s}")
    rows = {}
    for s, r in joined():
        amb = (s["fp"] & s["fn"]) & s["is_region"] & (s["eff"] > 1e-9) & (s["mass"] > _EPS)
        if not amb.any():
            continue
        E = np.maximum(s["eff"][amb], _EPS)
        tru = np.log10(np.maximum(s["G"][amb], 1.0) / E)
        e0 = np.abs(np.log10(np.maximum(r[0]["f0"][amb] * s["mass"][amb], 1.0) / E) - tru)
        e1 = np.abs(np.log10(np.maximum(r[1]["f0"][amb] * s["mass"][amb], 1.0) / E) - tru)
        for k in ("all", f"capture {s['group'][0]}",
                  "unstranded" if s["group"][2] == "0.50" else "stranded"):
            a = rows.setdefault(k, [0, [], []])
            a[0] += int(amb.sum())
            a[1].append(float(e0.mean()))
            a[2].append(float(e1.mean()))
    for k, (n, a, b) in rows.items():
        a, b = np.array(a), np.array(b)
        print(f"{k:22s} {n:7d} {a.mean():8.4f} {b.mean():9.4f} {b.mean() - a.mean():+8.4f} "
              f"{int((b < a).sum()):>4d}/{len(a)}")


ARMS = [
    ("#1  pass-0, AMBIG out (ships)", 0, False),
    ("    pass-0, AMBIG in (circular)", 0, True),
    ("#2  re-solved, AMBIG out", 1, False),
    ("#2  re-solved, AMBIG IN (the test)", 1, True),
]


def q2(keep, title, guard=False):
    print(f"\n=== Q2 [{title}] ===\n")
    acc = {a[0]: [] for a in ARMS}
    n_cond = 0
    for s, r in joined():
        if not keep(s):
            continue
        n_cond += 1
        mk = L.masks(s)
        base_sel = production_substrate(s, mk)
        orc = L.oracle_landscape(s, sel=L.live_region(s), knn_scale=_KNN)
        split = L.two_component(orc)["split"]
        hi = L.GRID > split
        nonamb = L.live_region(s) & ~(s["fp"] & s["fn"])
        orc_na = L.oracle_landscape(s, sel=nonamb, knn_scale=_KNN)
        om = float(orc[hi].sum())
        fab = L.GRID > 0.0                      # the fixed fabrication band (N5's gamut convention)
        for name, it, with_amb in ARMS:
            s2 = dict(s)
            s2["f0"] = r[it]["f0"]
            s2["g_hat"] = r[it]["f0"] * s["mass"]
            s2["var"] = r[it]["var"]
            sel = base_sel
            if with_amb:
                sel = sel | ((s["eff"] > 1e-9) & mk["region"] & mk["havemass"] & mk["ambig"])
            P = L.recipe(s2, sel=sel, w=L.recipe_weights(s2, sel, mk), knn_scale=_KNN)
            if P is None:
                continue
            acc[name].append((float(P[hi].sum()) / max(om, 1e-9), L.emd(P, orc), L.emd(P, orc_na),
                              float(P[fab].sum()), float(orc[fab].sum()), int(sel.sum())))
    print(f"{'arm':36s} {'enr recovery':>13s} {'EMD vs ALL':>11s} {'EMD vs NON-AMBIG':>17s} "
          f"{'n_train':>8s}" + ("   fabrication (fit/oracle)" if guard else ""))
    for name, _, _ in ARMS:
        a = np.array(acc[name])
        if not a.size:
            continue
        extra = ""
        if guard:
            fit_m, orc_m = a[:, 3].mean(), a[:, 4].mean()
            extra = f"   {fit_m:.4f} / {orc_m:.4f}  = {fit_m / max(orc_m, 1e-9):5.1f}x"
        print(f"{name:36s} {a[:, 0].mean():13.3f} {a[:, 1].mean():11.4f} {a[:, 2].mean():17.4f} "
              f"{a[:, 5].mean():8.0f}{extra}")
    print(f"\n  n = {n_cond} conditions.  'EMD vs NON-AMBIG' is W1's criterion: the prior is APPLIED to the")
    print("  non-AMBIG population, so scoring against an oracle that contains AMBIG rewards substrate-")
    print("  matching rather than accuracy.  Enriched recovery is at the ORACLE's split (N1).")


if __name__ == "__main__":
    q1()
    q2(lambda s: s["group"][1] != "none" and s["group"][0] in ("ON", "VSTRONG"),
       "gDNA-bearing capture-ON/VSTRONG, n=13 — where the enriched mode exists")
    q2(lambda s: s["group"][1] == "none",
       "ZERO-gDNA, n=9 — THE FABRICATION GUARD (any enriched mass here is false)", guard=True)
    q2(lambda s: True, "ALL 32 conditions — W1's own stratum")
