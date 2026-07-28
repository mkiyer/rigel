"""DOES THE DECLARED PRECISION STILL MATTER UNDER THE RESET STRATEGY?

`_fit_gdna_hyperprior` passes `var_g = belief.var_gdna[sel]` into `DensityNPMLE.fit`, where it is the
per-cell belief width tau in the Poisson-lognormal. So the substrate's DECLARED precision is an input to the
prior even if pass-0's beliefs are afterwards thrown away. Four fits per condition are on file:

    belief       arm's f_g, arm's widths        <- what ships
    beliefsharp  arm's f_g, ZERO widths         <- the width contributes nothing
    orc          oracle f_g, arm's widths
    orcsharp     oracle f_g, ZERO widths        <- the perfect prior

If TV(beliefsharp, orcsharp) > TV(belief, orcsharp) the declared widths are EARNING their place.

    python scratchpad/subacc_prec.py pk2 noown2
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
ARMS = sys.argv[1:] or ["noown2"]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}
CONDS = sorted(set(D[ARMS[0]]["cond"].tolist()))


def tv(d, ca, cb):
    ka, kb = f"{ca}|log_rho", f"{cb}|log_rho"
    if ka not in d or kb not in d:
        return None
    lra, lrb = d[ka], d[kb]
    g = np.linspace(min(lra[0], lrb[0]), max(lra[-1], lrb[-1]), 600)
    pa = np.exp(np.interp(g, lra, d[f"{ca}|logP"], left=-np.inf, right=-np.inf))
    pb = np.exp(np.interp(g, lrb, d[f"{cb}|logP"], left=-np.inf, right=-np.inf))
    pa, pb = pa / max(pa.sum(), _EPS), pb / max(pb.sum(), _EPS)
    return 0.5 * float(np.abs(pa - pb).sum())


print("=" * 104)
print("Distance to the PERFECT prior (orcsharp), with and without the solver's declared belief widths")
print("=" * 104)
print(f"{'arm':<10}{'TV(belief,perfect)':>20}{'TV(beliefsharp,perfect)':>26}{'widths help?':>16}"
      f"{'TV(orc,perfect)':>18}")
for a in ARMS:
    d = D[a]
    t1 = np.array([tv(d, f"{c}|belief", f"{c}|orcsharp") for c in CONDS if tv(d, f"{c}|belief", f"{c}|orcsharp")])
    t2 = np.array([tv(d, f"{c}|beliefsharp", f"{c}|orcsharp") for c in CONDS
                   if tv(d, f"{c}|beliefsharp", f"{c}|orcsharp")])
    t3 = np.array([tv(d, f"{c}|orc", f"{c}|orcsharp") for c in CONDS if tv(d, f"{c}|orc", f"{c}|orcsharp")])
    better = int(np.sum(t1 < t2))
    print(f"{a:<10}{t1.mean():>20.4f}{t2.mean():>26.4f}"
          f"{f'{better}/{t1.size} yes':>16}{t3.mean():>18.4f}")

print("\n  per-condition, HEAD arm: TV(belief,perfect) vs TV(beliefsharp,perfect)")
d = D[ARMS[-1]]
print(f"    {'condition':<48}{'with widths':>13}{'no widths':>12}{'delta':>10}")
rows = []
for c in CONDS:
    x, y = tv(d, f"{c}|belief", f"{c}|orcsharp"), tv(d, f"{c}|beliefsharp", f"{c}|orcsharp")
    if x is None or y is None:
        continue
    rows.append((y - x, c, x, y))
rows.sort(reverse=True)
for dl, c, x, y in rows[:6] + rows[-4:]:
    print(f"    {c[5:]:<48}{x:>13.3f}{y:>12.3f}{dl:>+10.3f}")
