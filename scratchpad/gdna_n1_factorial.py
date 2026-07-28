"""N1 — the FACTORIAL: substrate x weights x counts, against BOTH references.

The ladder (gdna_n1_ladder.py) attributes the truly-enriched weight-share loss one rung at a time, which
assumes the rungs compose.  This measures them jointly, and adds the question the ladder cannot answer:

    is the boundary dilution a DEFECT, or a MISMATCHED REFERENCE?

`oracle_landscape`'s default substrate is REGION-ONLY, so a fit trained on regions+boundaries is scored
against a population it was never asked to describe (production plan §W1a's own warning).  Every arm is
therefore scored twice:

    REGION reference  — oracle over live region nodes            (what §W3's 0.140 is)
    MATCHED reference — oracle over the SAME nodes the arm used  (isolates counts+weights from substrate)

Both use true counts G and uniform weights, so the reference differs from the arm only in what is being
tested.  Threshold: the REGION oracle's own split, held fixed across all arms of a condition.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_factorial.py
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds, kernels  # noqa: E402

_EPS = 1e-12


def substrates(s, mk):
    """The substrate axis. `prod` is production; the others change exactly one structural decision."""
    live = s["eff"] > 1e-9
    reg = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
    bnd = mk["boundary"] & live & mk["havemass"] & (mk["single"] | mk["gonly"])
    reg_a = (live & mk["region"] & mk["havemass"]) | mk["struct_zero"]          # regions incl. AMBIG
    bnd_a = mk["boundary"] & live & mk["havemass"]                              # boundaries incl. AMBIG
    return {
        "prod (reg+bnd, no AMBIG)": reg | bnd,
        "region-only, no AMBIG": reg,
        "prod + AMBIG": reg_a | bnd_a,
        "region-only + AMBIG": reg_a,
    }


def mass_above(sel, g, E, w, split):
    """The exact weighted kernel mass above `split` for a node set."""
    if not sel.any():
        return np.nan
    pn = kernels(g[sel], E[sel])
    wt = w / max(float(w.sum()), _EPS)
    return float((wt * pn[:, L.GRID > split].sum(1)).sum())


def main():
    ss = conds()
    axes = list(itertools.product(("prod (reg+bnd, no AMBIG)", "region-only, no AMBIG",
                                   "prod + AMBIG", "region-only + AMBIG"),
                                  ("recipe w", "flat w"), ("pass-0 g", "oracle G")))
    out = {a: [] for a in axes}
    ref_region, ref_matched = [], {a: [] for a in axes}
    for s in ss:
        split = L.two_component(L.oracle_landscape(s))["split"]
        mk = L.masks(s)
        subs = substrates(s, mk)
        E = s["eff"]
        prod_sel = L.recipe_substrate(s, mk)
        ref_region.append(mass_above(L.live_region(s), s["G"], E,
                                     np.ones(int(L.live_region(s).sum())), split))
        for a in axes:
            sname, wname, cname = a
            sel = subs[sname]
            n = int(sel.sum())
            if wname == "recipe w":
                # the production weight, evaluated on whichever substrate this arm uses
                w = L.recipe_weights(s, sel, mk)
                if sname != "prod (reg+bnd, no AMBIG)":
                    pass  # recipe_weights is per-node; substrate only selects which nodes
            else:
                w = np.ones(n)
            g = s["g_hat"] if cname == "pass-0 g" else s["G"]
            out[a].append(mass_above(sel, g, E, w, split))
            ref_matched[a].append(mass_above(sel, s["G"], E, np.ones(n), split))
        assert prod_sel.sum() == subs["prod (reg+bnd, no AMBIG)"].sum()

    r_reg = float(np.nanmean(ref_region))
    print("=== N1 factorial: enriched mass above the REGION-oracle split (13 conditions) ===\n")
    print(f"REGION reference (oracle, live regions, true counts, flat): {r_reg:.4f}\n")
    print(f"{'substrate':26s} {'weights':9s} {'counts':9s} {'mass':>7s} {'/REGION':>8s} "
          f"{'matched':>8s} {'/MATCHED':>9s}")
    for a in axes:
        m = float(np.nanmean(out[a]))
        rm = float(np.nanmean(ref_matched[a]))
        print(f"{a[0]:26s} {a[1]:9s} {a[2]:9s} {m:7.4f} {m / r_reg:8.3f} {rm:8.4f} {m / rm:9.3f}")

    print("\n  /REGION  = recovery against the §W3 headline reference (region-only oracle) — the 0.40 number")
    print("  /MATCHED = recovery against an oracle on the SAME nodes — what is left once the substrate's")
    print("             own composition is removed from the comparison.")


if __name__ == "__main__":
    main()
