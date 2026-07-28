"""Verify P4's '14 % manufactured precision' claim: PAIRED per-node Var(log f_g), HEAD vs arms."""
from __future__ import annotations

import sys

import numpy as np

base = np.load(sys.argv[1], allow_pickle=True)
arms = sys.argv[2:]
cb, mb, vb = base["cond"], base["mass"], base["var"]
ab, clb = base["amb"].astype(bool), base["cls"]
print(f"base nodes = {cb.size:,}")
for p in arms:
    d = np.load(p, allow_pickle=True)
    aligned = (d["cond"].size == cb.size) and bool((d["cond"] == cb).all()) and bool(
        np.allclose(d["mass"], mb)
    )
    v2 = d["var"]
    print(f"\n{p}   node sets aligned (cond + mass identical) = {aligned}   n={d['cond'].size:,}")
    fin = np.isfinite(vb) & np.isfinite(v2)
    rows = [("ALL", np.ones(cb.shape, bool))]
    for c in ("boundary", "exon", "intron", "intergenic"):
        for lab, m2 in ((" single", ~ab), (" AMBIG", ab)):
            rows.append((c + lab, (clb == c) & m2))
    print(f"  {'population':<22}{'n':>9}{'med HEAD':>11}{'med arm':>11}{'ratio':>8}"
          f"{'mean ratio':>12}{'frac up':>9}")
    for lab, m in rows:
        k = m & fin
        if k.sum() < 10:
            continue
        r = v2[k] / np.maximum(vb[k], 1e-300)
        print(f"  {lab:<22}{k.sum():>9,}{np.median(vb[k]):>11.5f}{np.median(v2[k]):>11.5f}"
              f"{np.median(v2[k]) / np.median(vb[k]):>8.3f}{np.median(r):>12.3f}{np.mean(r > 1.0):>9.1%}")
