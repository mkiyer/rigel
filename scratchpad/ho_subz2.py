"""SUBSTRATE ACCURACY *and* TRUST, per arm, from the `subacc_dump.py` npz files.

Splits every node into the two populations the hyperprior question is about:
  SUBSTRATE  = live & REGION & (single|gonly)            <- what `_fit_gdna_hyperprior` trains on
  EXCLUDED   = boundary | (REGION & AMBIG)               <- what it never sees
and reports, for each, the mass-weighted |f_g - oracle| (ACCURACY, the metric that matters if the refit
RESETS) and z2|Q1 (TRUST, the metric that matters if the refit ITERATES).

z2|Q1 uses a FIXED Var(f_g) threshold taken from the FIRST arm, so it is apples-to-apples across arms
(`p4_trust.py`'s own caveat).

    python scratchpad/ho_subz2.py Shd0 Samb0 Sexcl0
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
ARMS = sys.argv[1:]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}


def masks(d):
    isr = d["isr"].astype(bool)
    live = d["live"].astype(bool)
    amb = d["ambig"].astype(bool)
    ok = np.isfinite(d["fo"]) & (d["mass"] > _EPS)
    sub = live & isr & ~amb & ok
    exc = ((~isr) | (isr & amb)) & ok
    return ok, sub, exc, isr, amb


THR = None
print(f"{'arm':<10}{'population':<28}{'nodes':>9}{'mass':>15}{'mwae':>9}{'ERRmass':>12}"
      f"{'CW(fix)':>11}{'z2|Q1':>9}{'z2 all':>9}")
print("-" * 112)
for a in ARMS:
    d = D[a]
    ok, sub, exc, isr, amb = masks(d)
    mass, fg, fo, var = d["mass"], d["fg"], d["fo"], d["var"]
    err = np.abs(fg - fo)
    emass = err * mass
    fin = np.isfinite(var) & ok
    if THR is None:
        THR = float(np.quantile(var[fin], 0.25))
    conf = fin & (var <= THR)

    def z2(m):
        k = m & fin
        den = float(np.sum(mass[k] * var[k]))
        return float(np.sum(mass[k] * err[k] ** 2)) / den if den > 0 else float("nan")

    rt = d["rt"]
    rows = [
        ("SUBSTRATE (fit set)", sub),
        ("  x exon", sub & (rt == 2)),
        ("  x intron", sub & (rt == 1)),
        ("  x intergenic", sub & (rt == 0)),
        ("EXCLUDED (ambig|bnd)", exc),
        ("  x AMBIG region", exc & isr),
        ("  x boundary", exc & ~isr),
        ("ALL", ok),
    ]
    for lab, m in rows:
        if not m.any():
            continue
        print(f"{a:<10}{lab:<28}{int(m.sum()):>9,}{mass[m].sum():>15,.0f}"
              f"{np.average(err[m], weights=mass[m]):>9.4f}{emass[m].sum():>12,.0f}"
              f"{emass[m & conf].sum():>11,.0f}{z2(m & conf):>9.2f}{z2(m):>9.2f}")
    print("-" * 112)

if len(ARMS) > 1:
    print(f"\nΔ vs {ARMS[0]}  (mwae; negative = better)")
    d0 = D[ARMS[0]]
    ok0, sub0, exc0, isr0, amb0 = masks(d0)
    print(f"  {'population':<28}" + "".join(f"{a:>13}" for a in ARMS[1:]))
    for lab, sel in (("SUBSTRATE", "sub"), ("  x exon", "subex"), ("  x intron", "subin"),
                     ("EXCLUDED", "exc"), ("ALL", "ok")):
        row = f"  {lab:<28}"
        for a in ARMS[1:]:
            d = D[a]
            ok, sub, exc, isr, amb = masks(d)
            m = {"sub": sub, "subex": sub & (d["rt"] == 2), "subin": sub & (d["rt"] == 1),
                 "exc": exc, "ok": ok}[sel]
            m0 = {"sub": sub0, "subex": sub0 & (d0["rt"] == 2), "subin": sub0 & (d0["rt"] == 1),
                  "exc": exc0, "ok": ok0}[sel]
            e = np.average(np.abs(d["fg"][m] - d["fo"][m]), weights=d["mass"][m])
            e0 = np.average(np.abs(d0["fg"][m0] - d0["fo"][m0]), weights=d0["mass"][m0])
            row += f"{e - e0:>+13.4f}"
        print(row)
