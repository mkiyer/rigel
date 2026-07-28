"""W4 — does the PRODUCTION fitter reproduce the validated exploration recipe?

Everything the landscape's design rests on (N1, N5, W1-W3) was measured with
`scripts/debug/gdna_explore_lib.recipe` on a FIXED grid, `linspace(-5, 2.5, 260)`. The production
`gdna_landscape.fit_gdna_landscape` derives its grid from the data's own support instead, which removes an
asserted constant from the ledger but is a behaviour change and must not ride along unmeasured.

So: same substrate, same weights, same kernel, one difference at a time.

    A  exploration recipe, fixed grid                       (the validated object)
    B  production fitter, derived grid                      (what will ship)

Both are rendered onto the exploration grid for comparison, so the metric sees the landscape, not the axis.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w4_verify.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n5_roughness import region_substrate, roughness  # noqa: E402

from rigel.calibration.gdna_landscape import fit_gdna_landscape  # noqa: E402

LN10 = np.log(10.0)


def onto_explore_grid(landscape):
    """Render a GdnaLandscape onto the exploration grid as a normalised density, so A and B are compared
    as landscapes rather than as griddings."""
    lp = np.interp(L.GRID * LN10, landscape.log_rho, landscape.logP,
                   left=landscape.logP[0], right=landscape.logP[-1])
    d = np.exp(lp - lp.max())
    return d / d.sum()


def main():
    for suite in ("ambig", "quick"):
        rows = []
        for s in L.load_scenarios(suite):
            mk = L.masks(s)
            sel = region_substrate(s, mk)
            if int(sel.sum()) < 5:
                continue
            a = L.recipe(s, sel=sel, w=L.recipe_weights(s, sel, mk), knn_scale=0.5)
            lb = fit_gdna_landscape(
                s["g_hat"][sel], s["mass"][sel], s["eff"][sel], s["var"][sel],
                anchor=mk["struct_zero"][sel], strength=1.0,
            )
            if a is None or lb is None:
                continue
            b = onto_explore_grid(lb)
            orc = L.oracle_landscape(s, knn_scale=0.5)
            ta, tb = L.two_component(a), L.two_component(b)
            rows.append((
                L.emd(a, b), L.emd(a, orc), L.emd(b, orc),
                ta["enr_mass"], tb["enr_mass"], ta["enr_width"], tb["enr_width"],
                roughness(a, 1.0, 2.5), roughness(b, 1.0, 2.5),
                float(lb.log_rho[0] / LN10), float(lb.log_rho[-1] / LN10),
            ))
        r = np.array(rows).mean(0)
        print(f"=== {suite} (n={len(rows)}) ===")
        print(f"  EMD(A exploration, B production)   {r[0]:.4f}   <- the whole cost of the derived grid")
        print(f"  EMD to the corrected oracle        A {r[1]:.4f}   B {r[2]:.4f}   ({r[2] - r[1]:+.4f})")
        print(f"  enriched mass                      A {r[3]:.4f}   B {r[4]:.4f}   ({r[4] - r[3]:+.4f})")
        print(f"  enriched width                     A {r[5]:.4f}   B {r[6]:.4f}   ({r[6] - r[5]:+.4f})")
        print(f"  roughness above +1                 A {r[7]:.1f}      B {r[8]:.1f}")
        print(f"  derived grid span                  [{r[9]:+.2f}, {r[10]:+.2f}]  "
              f"(exploration: [-5.00, +2.50])\n")


if __name__ == "__main__":
    main()
