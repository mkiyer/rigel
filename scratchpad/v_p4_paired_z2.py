"""THE SELECTION-FREE TEST of P4's headline #7 ("`_far` is not the cause of boundary z2 = 5.58").

Both the self-quartile and the fixed-threshold views re-select the population per arm. Hold the node set
FIXED (the BASE arm's confident quartile) and let each arm bring its OWN error and its OWN variance:

    z2_arm = sum(mass * raw_arm^2) / sum(mass * var_arm)     over the SAME nodes

That is the paired calibration question with no selection at all.
"""
from __future__ import annotations

import sys

import numpy as np

import sys as _sys
_sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from z2 import lin_var  # noqa: E402  -- THE single z2 denominator (log Var -> linear Var)

_EPS = 1e-9
B = np.load(sys.argv[1], allow_pickle=True)
arms = sys.argv[2:]
mass, cls, amb = B["mass"], B["cls"], B["amb"].astype(bool)
vb = B["var"]
q1 = float(np.quantile(vb[np.isfinite(vb)], 0.25))
POPS = [("ALL", np.ones(mass.shape, bool)), ("exon single", (cls == "exon") & ~amb),
        ("exon AMBIG", (cls == "exon") & amb), ("boundary single", (cls == "boundary") & ~amb),
        ("boundary AMBIG", (cls == "boundary") & amb), ("intron single", (cls == "intron") & ~amb)]


def parts(d, m):
    ms, err, var = d["mass"], d["err"], d["var"]
    raw = np.where(ms > _EPS, err / np.maximum(ms, _EPS), 0.0)
    k = m & np.isfinite(var)
    return (float(np.sum(ms[k] * raw[k] ** 2)),
            float(np.sum(ms[k] * lin_var(var[k], d["fg"][k]))))


print(f"node set held FIXED = base's confident quartile (var <= q1 = {q1:.5g}); "
      f"each arm uses its OWN err and var on those SAME nodes")
hdr = f"  {'population':<18}{'n':>7}" + "".join(f"{a.split('_')[-1][:9]:>10}" for a in [sys.argv[1]] + arms)
print(hdr)
for lab, m in POPS:
    k = m & np.isfinite(vb) & (vb <= q1)
    if k.sum() < 10:
        continue
    out = []
    for p in [sys.argv[1]] + arms:
        d = np.load(p, allow_pickle=True)
        num, den = parts(d, k)
        out.append(num / den if den > 0 else float("nan"))
    print(f"  {lab:<18}{int(k.sum()):>7,}" + "".join(f"{v:>10.2f}" for v in out))

print("\n  and the two halves separately (numerator = mass-weighted MSE, denominator = declared var):")
for lab, m in POPS:
    k = m & np.isfinite(vb) & (vb <= q1)
    if k.sum() < 10:
        continue
    row = []
    for p in [sys.argv[1]] + arms:
        d = np.load(p, allow_pickle=True)
        num, den = parts(d, k)
        row.append((num, den))
    print(f"  {lab:<18}" + "".join(f"  MSE={n:>11,.0f} VAR={v:>11,.0f}" for n, v in row))
