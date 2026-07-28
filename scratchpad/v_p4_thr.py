"""Robustness of P4's trust conclusion: sweep the confidence threshold convention.

Three conventions, all applied to the SAME base+arm pair:
  self   each arm's own 25th percentile of Var(log f_g)   (what pass0_error_table does)
  fixed  the BASE arm's 25th percentile, applied to both  (P4's claim)
  rank   the most-confident N nodes of each arm, N = 25 % of nodes (equal-count, arm-relative)
plus a sweep of the fixed threshold over base quantiles 10/25/50 %.
"""
from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
base_p, arm_p = sys.argv[1], sys.argv[2]
B = np.load(base_p, allow_pickle=True)
A = np.load(arm_p, allow_pickle=True)
assert bool((B["cond"] == A["cond"]).all())
mass, cls, amb = B["mass"], B["cls"], B["amb"].astype(bool)


def z2(d, m):
    mass_, err, var = d["mass"], d["err"], d["var"]
    raw = np.where(mass_ > _EPS, err / np.maximum(mass_, _EPS), 0.0)
    k = m & np.isfinite(var)
    den = float(np.sum(mass_[k] * var[k]))
    return float(np.sum(mass_[k] * raw[k] ** 2)) / den if den > 0 else float("nan")


POPS = [("ALL", np.ones(mass.shape, bool)), ("boundary single", (cls == "boundary") & ~amb),
        ("boundary AMBIG", (cls == "boundary") & amb), ("exon single", (cls == "exon") & ~amb),
        ("intron single", (cls == "intron") & ~amb)]

print(f"base={base_p}  arm={arm_p}")
for qq in (0.10, 0.25, 0.50):
    tb = float(np.quantile(B["var"][np.isfinite(B["var"])], qq))
    ta = float(np.quantile(A["var"][np.isfinite(A["var"])], qq))
    selB_self, selA_self = B["var"] <= tb, A["var"] <= ta
    selA_fix = A["var"] <= tb
    # equal-count rank selection
    nsel = int(round(qq * mass.size))
    rk = np.argsort(A["var"])[:nsel]
    selA_rank = np.zeros(mass.shape, bool)
    selA_rank[rk] = True
    print(f"\n  q={qq:.0%}   base thr={tb:.5g}  arm thr={ta:.5g}  "
          f"arm nodes under base thr={selA_fix.mean():.1%}")
    print(f"    {'population':<20}{'base':>9}{'arm self':>10}{'arm fixed':>11}{'arm rank':>10}")
    for lab, m in POPS:
        print(f"    {lab:<20}{z2(B, m & selB_self):>9.2f}{z2(A, m & selA_self):>10.2f}"
              f"{z2(A, m & selA_fix):>11.2f}{z2(A, m & selA_rank):>10.2f}")
