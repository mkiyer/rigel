"""How many NODES (and how much MASS) sit in the N band where raising od turns on the
kappa=1/2 spurious-variance artifact?  The artifact saturates for N >~ 1/od.

  od = 0.0345 (a=14, today)  -> on above N ~ 29
  od = 0.10   (a=4.5, C/2)   -> on above N ~ 10
  od = 0.20   (a=2,  ceiling)-> on above N ~ 5

So the band N in (5, 29] is where shrink-to-ceiling ADDS spurious strand confidence
relative to today.
"""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")


def main():
    for name in ["LBX0190", "LBX0588", "MO_3021", "vcap"]:
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        pl = d["payload"]
        rc = np.asarray(pl.region_contained, float)
        reg = rc[:, 0] + rc[:, 1]                     # unspliced pos + neg per region
        fl = np.asarray(pl.boundary_flux_left, float)
        fr = np.asarray(pl.boundary_flux_right, float)
        bl = fl[:, 0] + fl[:, 1]
        br = fr[:, 0] + fr[:, 1]
        allN = np.concatenate([reg, bl, br])
        print(f"=== {name} ===  nodes {allN.size:,}  unspliced mass {allN.sum():,.0f}")
        bands = [(0.0, 1.0), (1.0, 5.0), (5.0, 10.0), (10.0, 29.0), (29.0, 100.0), (100.0, 1e18)]
        print(f"   {'band N':>14s} {'nodes':>12s} {'% nodes':>9s} {'mass':>14s} {'% mass':>9s}")
        for lo, hi in bands:
            m = (allN > lo) & (allN <= hi)
            print(f"   {f'({lo:g},{hi:g}]':>14s} {int(m.sum()):12,} {100*m.mean():8.2f}% "
                  f"{allN[m].sum():14,.0f} {100*allN[m].sum()/max(allN.sum(),1):8.2f}%")
        band = (allN > 5) & (allN <= 29)
        nz = allN > 0
        print(f"   -> band (5,29]: {int(band.sum()):,} nodes = {100*band.mean():.2f}% of all, "
              f"{100*band.sum()/max(nz.sum(),1):.2f}% of NONZERO nodes, "
              f"{100*allN[band].sum()/max(allN.sum(),1):.2f}% of mass")
        big = allN > 29
        print(f"   -> N>29 (artifact already saturated today): {int(big.sum()):,} nodes, "
              f"{100*allN[big].sum()/max(allN.sum(),1):.2f}% of mass")
        print()


if __name__ == "__main__":
    main()
