"""DISSECTION 2 — the FLOOR: what the retiring delta-pin does that the landscape does not.

Owner recollection (2026-07-27): *"in the production version we had a much higher uniform floor added to
the landscape... that seemed to help quite a bit in controlling the problem of snapping to the wrong mode."*
Confirmed in `npmle.py`, and it is not merely higher — it is **SCOPED**:

    uni = where(in_support, floor_eps/(bhi-blo) * step, 0.0)      floor_eps = 0.02, support = [q0.5%, q99.5%]
    # "interior uniform floor (fills −∞ valleys between modes) bounded to the fitted support;
    #  OUTSIDE the support the real Gaussian tails are left intact (the false-positive suppression)"

Two jobs, and they pull in OPPOSITE directions:
  * INSIDE the support a floor is protective — a node genuinely sitting in the valley between the two modes
    must not meet −∞ and get snapped to whichever mode is nearer.
  * OUTSIDE it a floor is harmful — the absence of any training node out at high rho IS information, carried
    by the kernels' own decaying tails, and flattening it is exactly a false-positive channel.

The landscape's Laplace floor is uniform over the WHOLE grid, so it does the first and breaks the second.
This measures whether that is what the `gdna_none` regression is made of — in `log P`, which is what psi
consumes, not in mass, which is what an earlier probe looked at and which cannot see a 9-nat difference
between two negligible numbers.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d2_floor.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n5_roughness import region_substrate  # noqa: E402

from rigel.calibration.gdna_landscape import fit_gdna_landscape  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402

LN10 = np.log(10.0)
PROBES = (-3.0, -2.0, -1.0, 0.0, 1.0)


def at(obj, x_dec):
    """log P at a given log10 density, RELATIVE TO THE PEAK — the pull psi actually feels."""
    return float(np.interp(x_dec * LN10, obj.log_rho, obj.logP,
                           left=obj.logP[0], right=obj.logP[-1]) - obj.logP.max())


def main():
    print("=== log P relative to peak (nats) at log10 rho_g = ... — LESS NEGATIVE = weaker suppression ===")
    print("    A false positive needs the prior to permit high density. On a ZERO-gDNA library every nat")
    print("    of permission at rho > 0 is a nat the strand/message evidence no longer has to overcome.\n")
    hdr = " ".join(f"{f'{p:+.0f}':>8s}" for p in PROBES)
    print(f"{'condition':44s} {'object':12s} {hdr}")
    for s in L.load_scenarios("ambig"):
        if s["group"][1] != "none" or s["group"][2] != "0.50":
            continue
        mk = L.masks(s)
        sel = region_substrate(s, mk)
        objs = [
            ("LANDSCAPE", fit_gdna_landscape(s["g_hat"][sel], s["mass"][sel], s["eff"][sel],
                                             s["var"][sel], anchor=mk["struct_zero"][sel])),
            ("delta-pin", DensityNPMLE.fit(s["g_hat"][sel], s["eff"][sel], var_g=s["var"][sel],
                                           bandwidth=0.15)),
        ]
        for i, (name, o) in enumerate(objs):
            lbl = s["cond"][5:] if i == 0 else ""
            print(f"{lbl:44s} {name:12s} " + " ".join(f"{at(o, p):8.1f}" for p in PROBES))
    print("\n  The delta-pin's tails are the REAL kernel tails outside its fitted support; the landscape's")
    print("  are lifted to a uniform floor everywhere. If the landscape is systematically LESS negative at")
    print("  rho > 0, that is the false-positive channel, and the fix is to SCOPE the floor, not shrink it.")


if __name__ == "__main__":
    main()
