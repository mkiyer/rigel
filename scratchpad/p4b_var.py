"""P4b — (a) the refit=1 half of the decomposition, (b) is the belief-free frame's z2 win a UNIFORM
variance inflation (mush) or a TARGETED one (honesty)? Per-class median Var(log f_g) ratios."""

from __future__ import annotations

import numpy as np

_EPS = 1e-9
ARMS = ["rho2_r0", "rho1_r0", "fix0.0_r0", "rho2_r1", "rho1_r1", "fix0.0_r1"]
D = {n: dict(np.load(f"/tmp/p4b_{n}.npz", allow_pickle=True)) for n in ARMS}
ref = D["rho2_r0"]
cond, mass, fo = ref["cond"], ref["mass"], ref["fo"]
cls, amb, solv, fit = ref["cls"], ref["amb"].astype(bool), ref["solv"].astype(bool), ref["fit"].astype(bool)

print("(2b) refit=1 — where does iteration 2's effect land?  (ERR reads = Σ mass·|f_g−oracle|)")
e2 = np.abs(D["rho2_r1"]["fg"] - fo) * mass
e1 = np.abs(D["rho1_r1"]["fg"] - fo) * mass
for lab, m in (("fit substrate", fit & solv), ("EXCLUDED (AMBIG+bnd)", ~fit & solv)):
    print(f"    {lab:<24} ERR rho2 {e2[m].sum():>11,.0f}  ERR rho1 {e1[m].sum():>11,.0f}"
          f"  Δ {e1[m].sum() - e2[m].sum():>+11,.0f}  ({(e1[m].sum() / max(e2[m].sum(), 1) - 1) * 100:+.1f}%)")

print("\n(3b) VARIANCE PROFILE — median stated Var(log f_g) per class, and the ratio vs HEAD.")
rows = [("ALL solvable", solv), ("fit substrate", fit & solv)]
for c in ("exon", "intron", "boundary"):
    for lab, mm in ((" single", ~amb), (" AMBIG", amb)):
        rows.append((c + lab, (cls == c) & mm & solv))
hdr = f"  {'population':<20}{'n':>8}" + "".join(f"{n:>13}" for n in ARMS)
print(hdr)
print("  " + "-" * (len(hdr) - 2))
for lab, m in rows:
    k = m & np.isfinite(D["rho2_r0"]["var"])
    if not k.any():
        continue
    print(f"  {lab:<20}{k.sum():>8,}" + "".join(f"{np.median(D[n]['var'][k]):>13.5g}" for n in ARMS))
print("\n  ratio to HEAD (same refit):")
for lab, m in rows:
    k = m & np.isfinite(D["rho2_r0"]["var"])
    if not k.any():
        continue
    r0 = np.median(D["fix0.0_r0"]["var"][k]) / max(np.median(D["rho2_r0"]["var"][k]), _EPS)
    r1 = np.median(D["fix0.0_r1"]["var"][k]) / max(np.median(D["rho2_r1"]["var"][k]), _EPS)
    print(f"  {lab:<20}  fix0.0/HEAD  r0 {r0:>7.2f}x   r1 {r1:>7.2f}x")
