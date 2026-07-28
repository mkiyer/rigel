"""W-derivation, part 2: verify the sandwich, price the RNA blow-up, and ask what W would be NEEDED.

Run: OMP_NUM_THREADS=1 python scratchpad/od_prior_weight2.py
"""
from __future__ import annotations

import pickle

import numpy as np

from od_prior_weight import EPS, SEEDS, i_max, terms

rng = np.random.default_rng(1)

print("=" * 100)
print("D1 — does the SANDWICH recover the nominal pair count when the model is TRUE?")
print("=" * 100)
print("   I_sand = (SUM scale)^2 / SUM (e_s - od*scale_s)^2    [Huber-White / GMM ratio variance]")
print(f"\n   {'setup':>34s} {'I=pairs':>12s} {'E[I_sand]':>12s} {'ratio':>8s}")
for label, sizes, od_true in (
    ("500 seeds n=10, od=0 (model true)", np.full(500, 10), 0.0),
    ("500 seeds n=50, od=0 (model true)", np.full(500, 50), 0.0),
    ("500 seeds n=50, od=0.1 (model true)", np.full(500, 50), 0.1),
    ("500 n=50 + 5 CONTAMINANTS", np.full(500, 50), 0.0),
):
    vals, moms = [], []
    for _ in range(300):
        n = sizes.astype(float)
        if od_true <= 0:
            k = rng.binomial(sizes, 0.5).astype(float)
        else:
            a = 0.5 * (1 - od_true) / od_true
            p = rng.beta(a, a, size=sizes.size)
            k = rng.binomial(sizes, p).astype(float)
        if "CONTAMIN" in label:  # 5 seeds are transcripts: 1500 fragments, 1% sense
            n = np.concatenate([n, np.full(5, 1500.0)])
            k = np.concatenate([k, np.full(5, 15.0)])
        e = (k - n / 2) ** 2 - n / 4
        sc = n * (n - 1) / 4
        mom = e.sum() / sc.sum()
        vals.append(sc.sum() ** 2 / max(((e - mom * sc) ** 2).sum(), EPS))
        moms.append(mom)
    n = sizes.astype(float)
    pairs = (n * (n - 1) / 2).sum()
    if "CONTAMIN" in label:
        pairs += 5 * 1500 * 1499 / 2
    print(
        f"   {label:>34s} {pairs:12,.0f} {np.mean(vals):12,.0f} {np.mean(vals) / pairs:8.4f}"
        f"   od_mom={np.mean(moms):.4f}"
    )
print("\n   Rows 1-3: ratio ~ 1 => the sandwich reproduces the nominal information when the model holds.")
print("   Row 4: 5 contaminated seeds out of 505 collapse it => it prices misspecification, free.")

print("\n" + "=" * 100)
print("D2 — the MoM is INADMISSIBLE: od is an intraclass correlation, so od <= 1 structurally")
print("=" * 100)
d = pickle.load(open(SEEDS, "rb"))
print(f"   {'sample':>40s} {'comp':>5s} {'kappa':>7s} {'od_mom':>10s} {'algebraic cap':>14s}")
for name, g in d.items():
    for comp in ("gdna", "rna"):
        if comp not in g:
            continue
        s, t, w, kappa, _ = g[comp]
        mu = 0.5 if comp == "gdna" else kappa
        ex, sc, nc = terms(s, t, w if comp == "gdna" else None, kappa, mu)
        m = sc > 0
        if not m.any():
            continue
        mom = float(ex[m].sum() / sc[m].sum())
        cap = max(mu, 1 - mu) / min(mu, 1 - mu)  # (mu_s - mu)^2 / mu(1-mu) at mu_s in {0,1}
        print(f"   {name[:40]:>40s} {comp:>5s} {kappa:7.4f} {mom:10.4f} {cap:14.1f}")
print("\n   'algebraic cap' = max_{mu_s} (mu_s-mu)^2 / (mu(1-mu)) — what a pure MEAN DISPLACEMENT")
print("   contributes to the excess/scale ratio.  The RNA fits sit AT that cap => they are measuring")
print("   antisense-dominant boundary sides, i.e. a wrong MEAN, not a dispersion.")

print("\n" + "=" * 100)
print("D3 — what W would be NEEDED to bring each library under the ceiling?")
print("=" * 100)
print("   solve (I*od_mom + W*od0)/(I+W) = od_max  for W, with od0 = od_max/2 = 0.1")
print(f"\n   {'sample':>40s} {'I=pairs':>15s} {'od_mom':>9s} {'W needed':>16s} {'W/I':>8s} {'tau needed':>11s}")
for name, g in d.items():
    if "gdna" not in g:
        continue
    s, t, w, kappa, _ = g["gdna"]
    ex, sc, nc = terms(s, t, w, kappa, 0.5)
    m = sc > 0
    mom = float(ex[m].sum() / sc[m].sum())
    I = float((nc[m] * (nc[m] - 1) / 2).sum())
    od_max, od0 = 0.2, 0.1
    if mom <= od_max:
        print(f"   {name[:40]:>40s} {I:15,.0f} {mom:9.4f} {'(already under)':>16s}")
        continue
    W = I * (mom - od_max) / (od_max - od0)
    print(
        f"   {name[:40]:>40s} {I:15,.0f} {mom:9.4f} {W:16,.0f} {W / I:8.2f} {1 / np.sqrt(W):11.6f}"
    )
print("\n   'tau needed' = the prior SD on od that W implies.  A prior claiming od is known to +-0.0005")
print("   is not assertable.  => NO defensible W rescues a biased MoM.  W is not the lever.")

print("\n" + "=" * 100)
print("D4 — empirical Bayes on the ADMISSIBLE (clipped) fits, and the closed forms")
print("=" * 100)
v = np.array([0.2, 0.0031, 0.2, 0.0923])
print(f"   clipped gDNA fits {v} -> tau^2={v.var(ddof=1):.5f}  W_eb={1 / v.var(ddof=1):.0f} pairs")
print("   (4 libraries, and each value is itself contaminated => not a defensible estimate)\n")
print(f"   {'a_ceil':>7s} {'od_max':>8s} {'od0=od_max/2':>13s} {'W=12/od_max^2':>14s} {'Beta(a,a) of od0':>17s}")
for a in (2.0, 3.0, 4.0):
    od_max = 1.0 / (2 * a + 1)
    od0 = od_max / 2
    print(
        f"   {a:7.0f} {od_max:8.4f} {od0:13.4f} {12 / od_max**2:14.0f} "
        f"{0.5 * (1 - od0) / od0:17.2f}"
    )
print("\n   closed form:  od0 = 1/(2(2a+1)),   W = 12(2a+1)^2  pairs.   a=2 -> od0=0.1, W=300.")
print(f"   for reference i_max at those od0: " + ", ".join(f"{i_max(1 / (2 * (2 * a + 1))):.0f}" for a in (2, 3, 4)))
