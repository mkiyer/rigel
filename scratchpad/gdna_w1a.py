"""W1a — FIX THE TWO SCORING CONFOUNDS before any substrate or bandwidth tuning is read.

Both confounds bias the very comparisons W1 (substrate) and W2 (kernel) are decided on.

  A. THE SMEAR.  `oracle_landscape` renders the TRUE counts through the SAME zero-native Poisson smoother
     the fits use, so fit-vs-oracle cancels the smoother — which is what makes it a fair test of the input
     counts, and also what makes it BLIND to how much the smoother inflates the width. `oracle_pointmass`
     is the population itself (unit mass per node at its true density, smoothed only by the grid's own
     resolution). The gap between them is the smear, and it decides deconvolve-vs-accept.

  B. THE BOUNDARY CONFOUND.  The recipe trains partly on BOUNDARY nodes; `oracle_landscape` defaults to
     REGION-ONLY. So "include boundaries" is currently scored against a population it was never fitting.
     Measured three ways here:
       B1  do boundary nodes have a DIFFERENT true density distribution than region nodes?  (if yes,
           including them shifts the landscape away from the region population — a real bias, independent
           of any solver error)
       B2  how much of the recipe's apparent error is just this mismatch?  (fit scored vs region-only
           oracle vs vs a MATCHED region+boundary oracle)
       B3  what does dropping boundaries from the fit actually cost, scored fairly?

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w1a.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402

SUITES = ("ambig", "quick")


def live(s):
    return s["eff"] > 1e-9


def _fmt(v, p=3):
    return f"{v:.{p}f}" if np.isfinite(v) else "  --"


def part_a():
    """THE SMEAR: how much wider is the Poisson-smoothed oracle than the population it represents?"""
    print("=" * 100)
    print("A. THE SMEAR — oracle_landscape (Poisson-smoothed) vs oracle_pointmass (the population itself)")
    print("=" * 100)
    print(f"{'suite':6s} {'stratum':22s} {'sd pointmass':>13s} {'sd landscape':>13s} {'inflation':>10s}"
          f" {'EMD':>8s}")
    for su in SUITES:
        scen = L.load_scenarios(su)
        for name, keep in (("all", lambda s: True),
                           ("capture OFF", lambda s: s["group"][0] == "OFF"),
                           ("capture ON", lambda s: s["group"][0] == "ON"),
                           ("zero-gDNA (none_*)", lambda s: s["group"][1] == "none")):
            rows = []
            for s in scen:
                if not keep(s):
                    continue
                pm, ls = L.oracle_pointmass(s), L.oracle_landscape(s)
                if pm is None or ls is None:
                    continue
                rows.append((L.spread(pm), L.spread(ls), L.emd(ls, pm)))
            if rows:
                a = np.array(rows)
                print(f"{su:6s} {name:22s} {_fmt(a[:, 0].mean()):>13s} {_fmt(a[:, 1].mean()):>13s} "
                      f"{_fmt(a[:, 1].mean() / max(a[:, 0].mean(), 1e-9), 2):>10s} "
                      f"{_fmt(a[:, 2].mean()):>8s}")
    print("\n  inflation = sd(landscape)/sd(pointmass). 1.00 = the smoother adds nothing (the population")
    print("  is already wider than the per-node counting error). >1 = measurement error is leaking into")
    print("  the reference, i.e. the smear is real and a deconvolution is warranted.")


def part_b():
    """THE BOUNDARY CONFOUND."""
    print()
    print("=" * 100)
    print("B1. Do BOUNDARY nodes have a different TRUE density distribution than REGION nodes?")
    print("=" * 100)
    print(f"{'suite':6s} {'stratum':22s} {'region med':>11s} {'bnd med':>9s} {'shift':>7s}"
          f" {'region sd':>10s} {'bnd sd':>7s} {'EMD(r,b)':>9s}")

    def med(P):
        return float(L.GRID[np.searchsorted(np.cumsum(P), 0.5)]) if P is not None else np.nan

    for su in SUITES:
        scen = L.load_scenarios(su)
        for name, keep in (("capture OFF", lambda s: s["group"][0] == "OFF"),
                           ("capture ON", lambda s: s["group"][0] == "ON"),
                           ("zero-gDNA (none_*)", lambda s: s["group"][1] == "none")):
            rows = []
            for s in scen:
                if not keep(s):
                    continue
                r = L.oracle_pointmass(s, live(s) & s["is_region"])
                b = L.oracle_pointmass(s, live(s) & ~s["is_region"])
                if r is None or b is None:
                    continue
                rows.append((med(r), med(b), L.spread(r), L.spread(b), L.emd(b, r)))
            if rows:
                a = np.array(rows)
                print(f"{su:6s} {name:22s} {_fmt(a[:, 0].mean(), 2):>11s} {_fmt(a[:, 1].mean(), 2):>9s} "
                      f"{_fmt(a[:, 1].mean() - a[:, 0].mean(), 2):>7s} {_fmt(a[:, 2].mean(), 2):>10s} "
                      f"{_fmt(a[:, 3].mean(), 2):>7s} {_fmt(a[:, 4].mean()):>9s}")
    print("\n  A large |shift| means the two populations genuinely differ, so training on boundaries to")
    print("  predict regions imports a bias — the owner's 'partial probe overlap underestimates' concern,")
    print("  measured on ORACLE truth with no solver involved.")

    print()
    print("=" * 100)
    print("B2/B3. How much of the recipe's error is the MISMATCH, and what does dropping boundaries cost?")
    print("=" * 100)
    print("  (WEIGHTING HELD FIXED across arms — only the substrate differs. Dropping a node class and the")
    print("   reliability weight at the same time is not a substrate experiment.)\n")
    print(f"{'suite':6s} {'arm':26s} {'vs REGION':>10s} {'vs MATCHED':>11s} {'vs POINTMASS':>13s}")
    for su in SUITES:
        scen = L.load_scenarios(su)
        out = {}
        for arm, drop_bnd in (("boundaries IN (today)", False), ("boundaries OUT", True)):
            reg, mat, pm = [], [], []
            for s in scen:
                mk = L.masks(s)
                train = L.recipe_substrate(s, mk)
                if drop_bnd:
                    train = train & s["is_region"]
                P = L.recipe(s, sel=train)          # SAME weights, only the substrate changes
                if P is None:
                    continue
                reg.append(L.emd(P, L.oracle_landscape(s)))                    # production target
                mat.append(L.emd(P, L.oracle_landscape(s, live(s) & train)))   # the nodes it FIT on
                pm.append(L.emd(P, L.oracle_pointmass(s)))                     # the population itself
            if reg:
                out[arm] = tuple(float(np.mean(x)) for x in (reg, mat, pm))
                print(f"{su:6s} {arm:26s} {_fmt(out[arm][0]):>10s} {_fmt(out[arm][1]):>11s} "
                      f"{_fmt(out[arm][2]):>13s}")
        if len(out) == 2:
            a, b = out["boundaries IN (today)"], out["boundaries OUT"]
            print(f"{'':6s} {'  -> cost of dropping':26s} {_fmt(b[0] - a[0], 4):>10s} "
                  f"{_fmt(b[1] - a[1], 4):>11s} {_fmt(b[2] - a[2], 4):>13s}   (+ = dropping is WORSE)")
    print("\n  'vs MATCHED' scores each arm against the true density of the nodes it actually fit (the")
    print("  estimator's own bias); 'vs REGION' is the production-relevant target; 'vs POINTMASS' is the")
    print("  population without the smoother. A large gap between columns IS the confound.")


if __name__ == "__main__":
    part_a()
    part_b()
