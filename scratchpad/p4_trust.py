"""z2 RESTRICTED TO THE CONFIDENT QUARTILE (`z2|Q1`) — the statistic the docs quote, which
`pass0_error_table.py` prints only over ALL quartiles.

⚠ `pass0_error_table.py` defines "confident" by the run's OWN 25th percentile of Var(log f_g). That makes
CWRONG **incomparable across arms** whenever an arm moves the variance distribution — which is exactly what a
precision change does. So this prints both: the per-arm quartile (`self`) and a FIXED threshold taken from the
first file on the command line (`fix`), which is the apples-to-apples number.

    python scratchpad/p4_trust.py BASE.npz [ARM.npz ...]
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
paths = sys.argv[1:] or ["/tmp/pass0_state.npz"]
THR = None
for p in paths:
    d = np.load(p, allow_pickle=True)
    mass, err, var = d["mass"], d["err"], d["var"]
    amb, cls = d["amb"].astype(bool), d["cls"]
    raw = np.where(mass > _EPS, err / np.maximum(mass, _EPS), 0.0)
    fin = np.isfinite(var)
    q1 = float(np.quantile(var[fin], 0.25))
    if THR is None:
        THR = q1
    selfconf = fin & (var <= q1)
    fixconf = fin & (var <= THR)

    def z2(m, _fin=fin, _mass=mass, _raw=raw, _var=var):
        k = m & _fin
        den = float(np.sum(_mass[k] * _var[k]))
        return float(np.sum(_mass[k] * _raw[k] ** 2)) / den if den > 0 else float("nan")

    print(f"\n=== {p}   ERR={err.sum():,.0f}   own q1={q1:.5g}  (fixed thr={THR:.5g})"
          f"   med Var(log f_g)={np.median(var[fin]):.5g}")
    print(f"  CWRONG  self-quartile {err[selfconf].sum():>10,.0f}   FIXED-threshold "
          f"{err[fixconf].sum():>10,.0f}   nodes@fixed {fixconf.mean():.1%}")
    print(f"  {'population':<26}{'ERRreads':>12}{'CW fix':>11}{'z2|Q1 fix':>11}{'z2|Q1 self':>12}")
    rows = [("ALL", np.ones(mass.shape, bool))]
    for c in ("exon", "boundary", "intron", "intergenic"):
        for lab, m2 in ((" single", ~amb), (" AMBIG", amb)):
            rows.append((c + lab, (cls == c) & m2))
    rows.append(("boundary both", cls == "boundary"))
    for lab, m in rows:
        if not m.any():
            continue
        print(f"  {lab:<26}{err[m].sum():>12,.0f}{err[m & fixconf].sum():>11,.0f}"
              f"{z2(m & fixconf):>11.2f}{z2(m & selfconf):>12.2f}")
