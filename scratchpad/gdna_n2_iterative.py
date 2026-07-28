"""N2 — the ITERATIVE resolution of the AMBIG circularity-vs-enriched-mass conflict.

HANDOFF_17 §0: AMBIG must be excluded from the prior's training set (it is the two-root ambiguity the prior
exists to resolve, and W1 measured that admitting it LOSES against the non-AMBIG population the prior is
applied to), yet AMBIG is where a third of the truly-enriched gDNA lives and excluding it costs ~20 % of the
enriched census (N1). The proposed way out uses the architecture that already ships:

    pass-0 -> prior #1 (AMBIG excluded) -> RE-SOLVE -> prior #2 (re-solved AMBIG INCLUDED)

which is non-circular by construction: the AMBIG estimates prior #2 trains on were produced by a prior that
never saw them.

Two questions, in order:
  Q1  does the re-solve actually IMPROVE the AMBIG nodes' density estimates? (if not, nothing else matters)
  Q2  does prior #2 recover enriched mass WITHOUT losing on the non-AMBIG population — the exact criterion
      that decided against admitting AMBIG in the first place?

⚠ Prior #1 here is the SHIPPED `DensityNPMLE`, the object this workstream is retiring. This measures the
ARCHITECTURE, not the final prior; repeat after W4.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n2_iterative.py   (needs scratchpad/gdna_refit_cache.pkl)
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
REFIT = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_refit_cache.pkl")


def joined(suite="ambig"):
    """The substrate cache (oracle truth + node types) joined to the paired belief cache, on condition."""
    rf = pickle.loads(REFIT.read_bytes())
    for s in L.load_scenarios(suite):
        r = rf.get(s["cond"])
        if r is None:
            continue
        assert r[0]["f0"].shape == s["f0"].shape and np.allclose(r[0]["f0"], s["f0"]), s["cond"]
        yield s, r


def q1():
    print("=== Q1: does the RE-SOLVE improve the AMBIG nodes' gDNA density? ===")
    print("    (mean |log10 rho_hat - log10 rho_true| over live AMBIG REGION nodes, in decades)\n")
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
        for k in ("all", f"capture {s['group'][0]}", "unstranded" if s["group"][2] == "0.50" else "stranded"):
            a = rows.setdefault(k, [0, [], []])
            a[0] += int(amb.sum())
            a[1].append(float(e0.mean()))
            a[2].append(float(e1.mean()))
    for k, (n, a, b) in rows.items():
        a, b = np.array(a), np.array(b)
        print(f"{k:22s} {n:7d} {a.mean():8.4f} {b.mean():9.4f} {b.mean() - a.mean():+8.4f} "
              f"{int((b < a).sum()):>4d}/{len(a)}")


def q2(keep, title):
    print(f"\n=== Q2 [{title}] — enriched mass, and the non-AMBIG EMD that decided against AMBIG ===\n")
    arms = [
        ("#1  pass-0, AMBIG out (today)", 0, False),
        ("    pass-0, AMBIG in (circular)", 0, True),
        ("#2  re-solved, AMBIG out", 1, False),
        ("#2  re-solved, AMBIG IN (the test)", 1, True),
    ]
    acc = {a[0]: [] for a in arms}
    for s, r in joined():
        if not keep(s):
            continue
        orc = L.oracle_landscape(s)
        split = L.two_component(orc)["split"]
        hi = L.GRID > split
        nonamb = L.live_region(s) & ~(s["fp"] & s["fn"])
        orc_na = L.oracle_landscape(s, sel=nonamb)          # the population the prior is APPLIED to
        om = float(orc[hi].sum())
        for name, it, with_amb in arms:
            # rebuild the scenario view with this iteration's beliefs, changing NOTHING else
            s2 = dict(s)
            s2["f0"] = r[it]["f0"]
            s2["g_hat"] = r[it]["f0"] * s["mass"]
            s2["var"] = r[it]["var"]
            mk = L.masks(s2)
            sel = L.recipe_substrate(s2, mk)
            if with_amb:
                live = s["eff"] > 1e-9
                sel = sel | (live & mk["havemass"] & mk["ambig"])
            P = L.recipe(s2, sel=sel, w=L.recipe_weights(s2, sel, mk), knn_scale=0.5)
            if P is None:
                continue
            acc[name].append((float(P[hi].sum()) / max(om, 1e-9), L.emd(P, orc), L.emd(P, orc_na)))
    print(f"{'arm':36s} {'enr recovery':>13s} {'EMD vs ALL':>11s} {'EMD vs NON-AMBIG':>17s}")
    for name, _, _ in arms:
        a = np.array(acc[name])
        print(f"{name:36s} {a[:, 0].mean():13.3f} {a[:, 1].mean():11.4f} {a[:, 2].mean():17.4f}")
    print("\n  'EMD vs NON-AMBIG' is W1's criterion: the prior is applied to the non-AMBIG population, so")
    print("  scoring against an oracle that contains AMBIG rewards substrate-matching rather than accuracy.")


if __name__ == "__main__":
    q1()
    q2(lambda s: s["group"][1] != "none" and s["group"][0] in ("ON", "VSTRONG"),
       "gDNA-bearing capture-ON/VSTRONG, n=13 — where the enriched mode exists")
    q2(lambda s: True, "ALL 32 conditions — W1's own stratum")
