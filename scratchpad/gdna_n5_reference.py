"""N5c — FIX THE SCORING INSTRUMENT, THEN RE-SELECT THE BANDWIDTH.

`gdna_n5_roughness.py` found the comb is in the ORACLE, not the fit: above +1 decade the reference's
roughness is **40.7** (a smooth bump scores 2-4) while the production kNN 0.5 fit scores **5.1**. The cause
is that `oracle_landscape` renders the TRUTH with a MEASUREMENT kernel — the Poisson likelihood of `G`,
width `1/(sqrt(G)*ln10)`, which shrinks as rho^(-1/2) and goes sub-grid exactly where the enriched mode is.
`G` is not an observation of `rho`; it IS `rho*E`. There is no posterior width to render.

That matters far beyond the picture: **every bandwidth and substrate decision in W1/W2/W3 was scored by EMD
against this comb**, and EMD to a comb rewards a fit that is also a comb. It is the mechanical reason
"EMD prefers h=0" and why kNN 2.0 was recorded as "over-widening the enriched component 2.2x".

Part A validates the corrected reference against ground truth that is known independently (owner + the
oracle census): the depleted level is a point mass at -1.28 (IQR 0.02-0.04 dec) and the enriched exon
population sits at +1.41 with spread 0.23-0.32 dec. A reference that renders those two facts is right; one
that renders the enriched mode as 40 spikes is not.

Part B re-runs the kernel sweep against BOTH references and reports whether the preference moves.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n5_reference.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds  # noqa: E402
from gdna_n5_roughness import region_substrate, roughness  # noqa: E402

SCALES = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)


def truth_stats(s, gmin=1.0):
    """The population as it actually is, from the oracle counts — no landscape, no kernel, no smoother.

    `gmin=1` is the LANDSCAPE's own convention (every live node, floored at the one-count wall), which is
    what a landscape must reproduce. `gmin=50` is the DENSITY-STRUCTURE convention the production plan used
    for the two-mode census, where `G/E` is precise; both are printed so neither is mistaken for the other."""
    live = L.live_region(s) & (s["G"] >= gmin if gmin > 1 else np.ones_like(s["G"], bool))
    x = np.log10(np.maximum(s["G"][live], 1.0)) - np.log10(np.maximum(s["eff"][live], 1e-12))
    dep = x[(s["ntype"][live] < 2)]
    enr = x[x > 0.0]
    return (float(np.median(dep)) if dep.size else np.nan,
            float(np.subtract(*np.percentile(dep, [75, 25]))) if dep.size else np.nan,
            float(np.median(enr)) if enr.size else np.nan,
            float(np.std(enr)) if enr.size > 2 else np.nan)


def part_a():
    ss = conds()
    print("=== A. Does the reference render the population that is actually there? ===\n")
    print("  GROUND TRUTH, straight off the oracle counts (no landscape involved):")
    for gmin, lbl in ((1.0, "all live nodes (the LANDSCAPE's convention)"),
                      (50.0, "G > 50 (the density-structure census)")):
        t = np.array([truth_stats(s, gmin) for s in ss])
        print(f"    {lbl:44s} depleted {np.nanmean(t[:, 0]):+.3f} dec (IQR {np.nanmean(t[:, 1]):.3f})"
              f" | enriched {np.nanmean(t[:, 2]):+.3f} dec (sd {np.nanmean(t[:, 3]):.3f})")
    t = np.array([truth_stats(s) for s in ss])
    print()
    print(f"{'oracle rendering':26s} {'dep loc':>8s} {'dep wid':>8s} {'enr loc':>8s} {'enr wid':>8s} "
          f"{'enr mass':>9s} {'roughness >+1':>14s}")
    for sc in SCALES:
        tc, rg = [], []
        for s in ss:
            P = L.oracle_landscape(s, knn_scale=sc)
            if P is None:
                continue
            tc.append(L.two_component(P))
            rg.append(roughness(P, 1.0, 2.5))
        name = "PRODUCTION (measurement)" if sc == 0 else f"population kNN {sc}"
        def m(k):
            return float(np.nanmean([x[k] for x in tc]))
        print(f"{name:26s} {m('dep_loc'):+8.3f} {m('dep_width'):8.3f} {m('enr_loc'):+8.3f} "
              f"{m('enr_width'):8.3f} {m('enr_mass'):9.4f} {np.nanmean(rg):14.1f}")
    print("\n  The enriched WIDTH is the discriminator: truth says ~0.23-0.32 dec. A reference reporting")
    print("  much less is resolving spikes it has no sample to support.")


def part_b():
    ss = conds()
    print("\n=== B. Select the fit's kernel by SHAPE, against the corrected reference ===")
    print("    (substrate = REGIONS ONLY, per the owner's 2026-07-27 decision. EMD is reported but is")
    print("     nearly flat in the kernel and monotone in smoothing — it does not discriminate here.)\n")
    ref = [L.two_component(L.oracle_landscape(s, knn_scale=0.5)) for s in ss]

    def rm(rs, k):
        return float(np.nanmean([x[k] for x in rs]))
    print(f"{'CORRECTED REFERENCE':22s} {rm(ref, 'dep_loc'):+8.3f} {rm(ref, 'dep_width'):8.3f} "
          f"{rm(ref, 'enr_loc'):+8.3f} {rm(ref, 'enr_width'):8.3f} {rm(ref, 'enr_mass'):9.4f}")
    print(f"\n{'fit kernel':22s} {'dep loc':>8s} {'dep wid':>8s} {'enr loc':>8s} {'enr wid':>8s} "
          f"{'enr mass':>9s} {'EMD':>8s} {'roughness >+1':>14s}")
    for sc in SCALES:
        tc, rg, em = [], [], []
        for s, r in zip(ss, ref):
            mk = L.masks(s)
            sel = region_substrate(s, mk)
            P = L.recipe(s, sel=sel, w=L.recipe_weights(s, sel, mk), knn_scale=sc)
            if P is None:
                continue
            tc.append(L.two_component(P, split=r["split"]))
            rg.append(roughness(P, 1.0, 2.5))
            em.append(L.emd(P, L.oracle_landscape(s, knn_scale=0.5)))
        name = "h=0 (no resolution)" if sc == 0 else f"kNN {sc}"
        print(f"{name:22s} {rm(tc, 'dep_loc'):+8.3f} {rm(tc, 'dep_width'):8.3f} "
              f"{rm(tc, 'enr_loc'):+8.3f} {rm(tc, 'enr_width'):8.3f} {rm(tc, 'enr_mass'):9.4f} "
              f"{np.nanmean(em):8.4f} {np.nanmean(rg):14.1f}")
    print("\n  The gate is the SHAPE: match the reference's enriched width and mass without smearing the")
    print("  depleted component, and get roughness into the smooth-bump band (~2-4 per decade).")


if __name__ == "__main__":
    part_a()
    part_b()
