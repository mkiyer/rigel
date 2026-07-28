"""P4b — compare an arbitrary set of captured arms on the same held-fixed node set.

    python scratchpad/p4b_arms.py rho2_r0 rho1_r0 fix1.0_r0 fix0.5_r0 fix0.0_r0

The FIRST arm is the reference: it defines the confident-quartile node set used for every arm's trust view.
"""

from __future__ import annotations

import sys
import numpy as np

_EPS = 1e-9
names = sys.argv[1:]
D = {n: dict(np.load(f"/tmp/p4b_{n}.npz", allow_pickle=True)) for n in names}
ref = D[names[0]]
for n, v in D.items():
    assert (v["cond"] == ref["cond"]).all() and (v["nid"] == ref["nid"]).all(), n

cond, mass, fo = ref["cond"], ref["mass"], ref["fo"]
cls, amb = ref["cls"], ref["amb"].astype(bool)
solv, fit = ref["solv"].astype(bool), ref["fit"].astype(bool)
conds = sorted(set(cond.tolist()))
STRATA = {
    "ALL 32": lambda c: True,
    "stranded ss_0.99": lambda c: "ss_0.99" in c,
    "unstranded x capON": lambda c: "ss_0.50" in c and "capture_on" in c,
    "capture OFF": lambda c: "capture_off" in c,
    "capture ON": lambda c: "capture_on" in c,
    "verystrong": lambda c: "verystrong" in c,
    "gdna_none": lambda c: c.startswith("gdna_none"),
}


def cm(f):
    return np.array([f(c) for c in cond])


print(f"{'stratum':<24}" + "".join(f"{n:>13}" for n in names))
print("-" * (24 + 13 * len(names)))
for lab, f in STRATA.items():
    m = cm(f)
    print(f"{lab:<24}" + "".join(
        f"{float(np.average(np.abs(D[n]['fg'][m] - fo[m]), weights=mass[m])):>13.4f}" for n in names))
m = fit & solv
print(f"{'FIT SUBSTRATE':<24}" + "".join(
    f"{float(np.average(np.abs(D[n]['fg'][m] - fo[m]), weights=mass[m])):>13.4f}" for n in names))

# per-condition better/worse vs the reference
print()
for n in names[1:]:
    b = w = fl = 0
    for c in conds:
        k = cond == c
        x = float(np.average(np.abs(D[names[0]]["fg"][k] - fo[k]), weights=mass[k]))
        y = float(np.average(np.abs(D[n]["fg"][k] - fo[k]), weights=mass[k]))
        if abs(y - x) < 5e-5:
            fl += 1
        elif y < x:
            b += 1
        else:
            w += 1
    print(f"  {n:<14} vs {names[0]}:  better {b} / worse {w} / flat {fl}")

# trust view on the reference arm's confident quartile
v0 = D[names[0]]["var"]
fin = np.isfinite(v0) & solv
q1 = float(np.quantile(v0[fin], 0.25))
conf = fin & (v0 <= q1)
print(f"\nTRUST (held-fixed node set = {names[0]}'s confident quartile, thr {q1:.5g}, {conf.sum():,} nodes)")
rows = [("ALL (solvable)", solv), ("FIT SUBSTRATE", fit & solv)]
for c in ("exon", "intron", "boundary"):
    for lab, mm in ((" single", ~amb), (" AMBIG", amb)):
        rows.append((c + lab, (cls == c) & mm & solv))
print(f"  {'population':<20}" + "".join(f"{n[:11]:>12}" for n in names) + "   |" +
      "".join(f"{'z2 ' + n[:8]:>12}" for n in names))
for lab, m in rows:
    k = m & conf
    if not k.any():
        continue
    cw, z2 = [], []
    for n in names:
        d = D[n]
        e = np.abs(d["fg"] - fo) * mass
        cw.append(float(e[k].sum()))
        den = float(np.sum(mass[k] * d["var"][k]))
        z2.append(float(np.sum(mass[k] * (e[k] / np.maximum(mass[k], _EPS)) ** 2)) / den if den > 0 else np.nan)
    print(f"  {lab:<20}" + "".join(f"{x:>12,.0f}" for x in cw) + "   |" + "".join(f"{x:>12.2f}" for x in z2))
