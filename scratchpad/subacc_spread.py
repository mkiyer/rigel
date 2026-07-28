"""(5) THE SUBSTRATE'S ACCURACY SPREAD + the controls the headline needs.

  * where the >0.1 / >0.25 substrate error lives (condition x class), and the signed bias
  * a scale for (4)'s TV: how far apart are two DIFFERENT conditions' oracle priors?
  * floor sensitivity for the density metric

    python scratchpad/subacc_spread.py pk2 p1dp p4land noown2 p1e1
"""

from __future__ import annotations

import sys

import numpy as np

_EPS = 1e-9
_LN10 = np.log(10.0)
ARMS = sys.argv[1:] or ["noown2"]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}
HEAD = ARMS[-1]


def sel_em(d):
    return (d["live"].astype(bool) & d["isr"].astype(bool)
            & (d["single"].astype(bool) | d["gonly"].astype(bool)) & ~d["intergenic"].astype(bool)
            & np.isfinite(d["fo"]))


print("=" * 122)
print("(5a) SIGNED bias on the substrate: mass-weighted mean (f_g - oracle).  >0 = the prior will be fitted "
      "TOO gDNA-rich")
print("=" * 122)
print(f"{'population':<40}" + "".join(f"{a:>13}" for a in ARMS))
for lab, key in (("SUBSTRATE", "all"), ("  exon", "exon"), ("  intron", "intron"),
                 ("  capture OFF", "coff"), ("  capture ON/verystrong", "con"),
                 ("  stranded ss_0.99", "s99"), ("  unstranded ss_0.50", "s50")):
    row = f"{lab:<40}"
    for a in ARMS:
        d = D[a]
        s = sel_em(d)
        c = d["cond"]
        m = {"all": s, "exon": s & (d["rt"] == 2), "intron": s & (d["rt"] == 1),
             "coff": s & np.array(["capture_off" in x for x in c]),
             "con": s & np.array(["capture_off" not in x for x in c]),
             "s99": s & np.array(["ss_0.99" in x for x in c]),
             "s50": s & np.array(["ss_0.50" in x for x in c])}[key]
        row += f"{np.average(d['fg'][m] - d['fo'][m], weights=d['mass'][m]):>13.4f}"
    print(row)

print("\n" + "=" * 122)
print("(5b) WHERE the substrate's >0.1 error mass lives (HEAD arm = " + HEAD + ")")
print("=" * 122)
d = D[HEAD]
s = sel_em(d)
e = np.abs(d["fg"] - d["fo"])
w = d["mass"]
tot = w[s].sum()
for thr in (0.1, 0.25):
    bad = s & (e > thr)
    print(f"\n  |err| > {thr}:  {bad.sum():,} nodes, {w[bad].sum():,.0f} mass "
          f"({w[bad].sum() / tot:.1%} of substrate mass)")
    print(f"    {'condition':<48}{'bad mass':>13}{'share of bad':>14}{'of its own sub mass':>21}")
    rows = sorted({c for c in d["cond"]}, key=lambda c: -w[bad & (d["cond"] == c)].sum())
    cum = 0.0
    for c in rows[:8]:
        m = bad & (d["cond"] == c)
        own = w[s & (d["cond"] == c)].sum()
        cum += w[m].sum()
        print(f"    {c[5:]:<48}{w[m].sum():>13,.0f}{w[m].sum() / w[bad].sum():>14.1%}"
              f"{w[m].sum() / max(own, _EPS):>21.1%}")
    print(f"    {'(top 8 cumulative)':<48}{cum:>13,.0f}{cum / w[bad].sum():>14.1%}")
    for lab, m2 in ((" exon", d["rt"] == 2), (" intron", d["rt"] == 1)):
        m = bad & m2
        print(f"    class{lab:<43}{w[m].sum():>13,.0f}{w[m].sum() / w[bad].sum():>14.1%}"
              f"{w[m].sum() / max(w[s & m2].sum(), _EPS):>21.1%}")

print("\n" + "=" * 122)
print("(5c) CONTROL for the TV scale in (4): TV between DIFFERENT conditions' ORACLE priors")
print("=" * 122)
CONDS = sorted(set(d["cond"].tolist()))


def cmp_curves(lra, lpa, lrb, lpb):
    lo, hi = min(lra[0], lrb[0]), max(lra[-1], lrb[-1])
    g = np.linspace(lo, hi, 600)
    pa = np.exp(np.interp(g, lra, lpa, left=-np.inf, right=-np.inf))
    pb = np.exp(np.interp(g, lrb, lpb, left=-np.inf, right=-np.inf))
    pa, pb = pa / max(pa.sum(), _EPS), pb / max(pb.sum(), _EPS)
    return 0.5 * float(np.abs(pa - pb).sum())


tvs = []
for i, ci in enumerate(CONDS):
    for cj in CONDS[i + 1:]:
        ka, kb = f"{ci}|orcsharp|", f"{cj}|orcsharp|"
        if ka + "logP" in d and kb + "logP" in d:
            tvs.append(cmp_curves(d[ka + "log_rho"], d[ka + "logP"], d[kb + "log_rho"], d[kb + "logP"]))
tvs = np.array(tvs)
print(f"  cross-condition oracle-vs-oracle TV: mean {tvs.mean():.3f}  median {np.median(tvs):.3f}  "
      f"p10 {np.percentile(tvs, 10):.3f}  p90 {np.percentile(tvs, 90):.3f}   (n={tvs.size} pairs)")
print("  => the fitted-vs-oracle TV of ~0.14 must be read against THIS as the 'no information' scale.")

print("\n" + "=" * 122)
print("(5d) FLOOR SENSITIVITY for the density metric dlog(rho_g) = |log(f_g*M) - log(fo*M)|")
print("=" * 122)
print(f"{'floor on g_hat':<40}" + "".join(f"{a:>13}" for a in ARMS))
for lab, fl in (("1 count (the KDE resolution wall)", 1.0), ("0.1 count", 0.1), ("1e-6 (the _collapse g_eps)", 1e-6)):
    row = f"{lab:<40}"
    for a in ARMS:
        da = D[a]
        sa = sel_em(da)
        M = da["mass"][sa]
        gh = np.maximum(da["fg"][sa] * M, fl)
        go = np.maximum(da["fo"][sa] * M, fl)
        row += f"{np.average(np.abs(np.log(gh) - np.log(go)), weights=M):>13.4f}"
    print(row)
print("\n  NOTE M and E_g are belief-free, so rho_g error IS f_g error in log space: the only difference "
      "between\n  the two metrics is the weighting of small f_g.")
