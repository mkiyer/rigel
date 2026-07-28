"""ENRICHED-MODE CENSUS — does the fitted landscape carry the enriched mode at all, and where?

At the 2026-07-21 pause the landscape was recorded as ENRICHED-BLIND: "the enriched mode has ~no mass"
(the reliability weight silenced the high-variance unstranded enriched nodes), so the projection had
nothing to pull toward. This measures, per capture-ON condition, the mass and LOCATION of the enriched
mode in the fit vs the oracle.

  enriched region = log10 rho_g > THR (the capture-enriched gDNA level; the depleted bulk sits < -1)
  mass  = integrated probability there        location = argmax within it

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_enriched_mode.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_resurrect import load_recipe  # noqa: E402

THR = 0.0   # decades; the enriched mode sits at ~+0.7..+1.3 under capture, the depleted bulk at < -1


def census(suite="ambig"):
    rec = load_recipe()
    hi = L.GRID > THR
    print(f"=== {suite}: enriched-mode census (mass and location above log10 rho_g = {THR:+.1f}) ===")
    print(f"{'condition':52s} {'oracle mass':>11s} {'fit mass':>9s} {'ratio':>7s} | "
          f"{'oracle loc':>10s} {'fit loc':>8s} {'shift':>7s}")
    rows = []
    for s in L.load_scenarios(suite):
        if s["group"][0] == "OFF":
            continue                       # the enriched mode only exists under capture
        orc, fit = L.oracle_landscape(s), rec(s)
        if orc is None or fit is None:
            continue
        om, fm = float(orc[hi].sum()), float(fit[hi].sum())
        if om < 1e-4:
            continue                       # no enriched mode in the truth -> nothing to recover
        ol = float(L.GRID[hi][np.argmax(orc[hi])])
        fl = float(L.GRID[hi][np.argmax(fit[hi])])
        rows.append((om, fm))
        print(f"{s['cond']:52s} {om:11.4f} {fm:9.4f} {fm / om:7.2f} | "
              f"{ol:+10.2f} {fl:+8.2f} {fl - ol:+7.2f}")
    if rows:
        om = np.array([r[0] for r in rows])
        fm = np.array([r[1] for r in rows])
        print(f"\n  MEAN enriched-mass recovery ratio (fit/oracle): {np.mean(fm / om):.2f}   "
              f"[1.00 = the mode is carried at its true mass; ~0 = 'enriched-blind', the 07-21 state]")


if __name__ == "__main__":
    for su in ("ambig", "quick"):
        census(su)
        print()
