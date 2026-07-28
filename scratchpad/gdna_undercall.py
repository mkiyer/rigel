"""Re-measure the ENRICHED-gDNA UNDER-CALL table — `archive/enrichment_sensitivity_worklog.md` §1.

That table is the reason the hyperprior track was PAUSED on 2026-07-21: pass-0 systematically under-called
gDNA density at enriched, unstranded, single-strand nodes (-0.13..-0.30 decades), so the fitted landscape's
enriched mode was under-massed and the projection could not pull anything up. The owner's call was to fix
pass-0 at the source and come back. Pass-0 has since been rebuilt. This re-runs the identical measurement.

Substrate: single-strand REGION nodes (NOT ambig - those are what the prior must predict), G > 0.
  obs = log10(g_hat/E)   tru = log10(G/E)   under-call = obs - tru  (negative = under-called)

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_undercall.py [ambig|quick]
"""
from __future__ import annotations

import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import gdna_explore_lib as L  # noqa: E402

# The 2026-07-21 measurement, for the side-by-side (worklog §1): (all single-strand, enriched single-strand).
# ⚠ Those numbers were taken with `capture_verystrong` mislabelled as OFF (fixed in gdna_cache_build._group),
# so the OFF rows are not a clean like-for-like — VSTRONG has its own row here and had none there.
WAS = {("OFF", "0.50"): (+0.415, -0.135), ("OFF", "0.99"): (-0.026, 0.0),
       ("ON", "0.50"): (+0.097, -0.169), ("ON", "0.99"): (-0.048, -0.012)}


def main(suite="ambig"):
    scen = L.load_scenarios(suite)
    cell = defaultdict(lambda: ([], []))   # (cap, ss) -> (all-single deltas, enriched-single deltas)
    per_scen = []
    for s in scen:
        mk = L.masks(s)
        live = (s["eff"] > 1e-9) & (s["mass"] > 1e-12)
        sel = live & mk["region"] & mk["single"] & (s["G"] > 0)
        if sel.sum() < 5:
            continue
        obs = np.log10(np.maximum(s["g_hat"], 1e-9) / s["eff"])
        tru = np.log10(np.maximum(s["G"], 1e-9) / s["eff"])
        d = (obs - tru)[sel]
        enr = (tru > 0.0)[sel]
        cap, _dna, ss, _nr = s["group"]
        cell[(cap, ss)][0].append(float(d.mean()))
        if enr.sum() >= 5:
            cell[(cap, ss)][1].append(float(d[enr].mean()))
            per_scen.append((s["cond"], int(enr.sum()), float(d[enr].mean())))

    print(f"=== {suite}: single-strand REGION under-call (obs - oracle), decades ===")
    print(f"{'cap':>8} {'ss':>5} | {'all single':>11} {'(was)':>8} | {'ENRICHED':>10} {'(was)':>8}  n_scen")
    for k in sorted(cell):
        a, e = cell[k]
        wa, we = WAS.get(k, (np.nan, np.nan))
        es = f"{np.mean(e):+.3f}" if e else "    n/a"
        ws = f"{wa:+.3f}" if np.isfinite(wa) else "   --"
        wse = f"{we:+.3f}" if np.isfinite(we) else "   --"
        print(f"{k[0]:>8} {k[1]:>5} | {np.mean(a):+11.3f} {ws:>8} | {es:>10} {wse:>8}  {len(a)}")

    per_scen.sort(key=lambda r: r[2])
    print("\n  worst ENRICHED under-callers (the 2026-07-21 dissection targets):")
    for cond, n, dd in per_scen[:6]:
        print(f"    {dd:+.3f}  n_enr={n:4d}  {cond}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ambig")
