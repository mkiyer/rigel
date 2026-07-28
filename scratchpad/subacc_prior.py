"""(4) HOW GOOD COULD THE PRIOR BE — the achievable gDNA hyperprior vs the ORACLE-fitted one.

`subacc_dump.py` fitted `_fit_gdna_hyperprior` three ways per condition: from the arm's belief, from the
ORACLE f_g on the same substrate nodes with the arm's own belief widths (`orc`), and from the oracle with
zero width (`orcsharp` — the perfect prior). This scores the gap: both P(log rho) are rendered on a common
grid and compared by total variation / Hellinger, and by where their mass sits (decades of log10 rho).

    python scratchpad/subacc_prior.py pk2 p1dp p4land noown2 p1e1
"""

from __future__ import annotations

import sys

import numpy as np

_LN10 = np.log(10.0)
ARMS = sys.argv[1:] or ["noown2"]
D = {a: dict(np.load(f"/tmp/subacc_{a}.npz", allow_pickle=True)) for a in ARMS}
d0 = D[ARMS[0]]
CONDS = sorted(set(d0["cond"].tolist()))


def curve(d, cond, tag, grid):
    k = f"{cond}|{tag}|"
    if k + "logP" not in d:
        return None
    lr, lp = d[k + "log_rho"], d[k + "logP"]
    p = np.exp(np.interp(grid, lr, lp, left=-np.inf, right=-np.inf))
    s = p.sum()
    return p / s if s > 0 else None


def stats(p, grid):
    m = float((p * grid).sum()) / _LN10
    v = float((p * (grid - m * _LN10) ** 2).sum()) / _LN10**2
    c = np.cumsum(p)
    q = [float(grid[np.searchsorted(c, t)]) / _LN10 for t in (0.1, 0.5, 0.9)]
    return m, np.sqrt(max(v, 0.0)), q


def compare(d, cond, tag_a, tag_b):
    ka, kb = f"{cond}|{tag_a}|log_rho", f"{cond}|{tag_b}|log_rho"
    if ka not in d or kb not in d:
        return None
    lo = min(d[ka][0], d[kb][0])
    hi = max(d[ka][-1], d[kb][-1])
    grid = np.linspace(lo, hi, 600)
    pa, pb = curve(d, cond, tag_a, grid), curve(d, cond, tag_b, grid)
    if pa is None or pb is None:
        return None
    tv = 0.5 * float(np.abs(pa - pb).sum())
    hel = float(np.sqrt(max(0.0, 1.0 - np.sqrt(pa * pb).sum())))
    sa, sb = stats(pa, grid), stats(pb, grid)
    return tv, hel, sa, sb


print("=" * 126)
print("(4) THE FITTED HYPERPRIOR: solver-substrate belief vs ORACLE on the same nodes")
print("    TV = total-variation distance between the two P(log rho) (0 = identical, 1 = disjoint)")
print("    d(median) / d(mean) in DECADES of log10 rho_g;  width = sd of log10 rho under each")
print("=" * 126)
for tag, lab in (("orc", "ORACLE f_g, arm's own belief widths"), ("orcsharp", "ORACLE f_g, zero width (PERFECT)")):
    print(f"\n--- vs {lab} ---")
    print(f"{'arm':<10}{'TV mean':>9}{'TV med':>9}{'TV p90':>9}{'Hel mean':>10}"
          f"{'|d median|':>12}{'|d mean|':>10}{'width fit':>11}{'width orc':>11}{'cells fit':>11}{'cells orc':>11}")
    for a in ARMS:
        d = D[a]
        tvs, hels, dmed, dmu, wf, wo, cf, co = [], [], [], [], [], [], [], []
        for c in CONDS:
            r = compare(d, c, "belief", tag)
            if r is None:
                continue
            tv, hel, sa, sb = r
            tvs.append(tv)
            hels.append(hel)
            dmed.append(abs(sa[2][1] - sb[2][1]))
            dmu.append(abs(sa[0] - sb[0]))
            wf.append(sa[1])
            wo.append(sb[1])
            cf.append(float(d[f"{c}|belief|ncell"][0]))
            co.append(float(d[f"{c}|{tag}|ncell"][0]))
        print(f"{a:<10}{np.mean(tvs):>9.3f}{np.median(tvs):>9.3f}{np.percentile(tvs, 90):>9.3f}"
              f"{np.mean(hels):>10.3f}{np.mean(dmed):>12.3f}{np.mean(dmu):>10.3f}"
              f"{np.mean(wf):>11.3f}{np.mean(wo):>11.3f}{np.mean(cf):>11.0f}{np.mean(co):>11.0f}")

print("\n--- per-condition TV(belief, ORACLE-sharp), and the two priors' medians (log10 rho_g) ---")
print(f"{'condition':<48}" + "".join(f"{a:>9}" for a in ARMS) + f"{'med(fit)':>10}{'med(orc)':>10}")
rows = []
for c in CONDS:
    r = compare(D[ARMS[-1]], c, "belief", "orcsharp")
    rows.append((r[0] if r else -1.0, c, r))
rows.sort(reverse=True)
for tv0, c, r in rows:
    line = f"{c[5:]:<48}"
    for a in ARMS:
        rr = compare(D[a], c, "belief", "orcsharp")
        line += f"{rr[0]:>9.3f}" if rr else f"{'-':>9}"
    line += f"{r[2][2][1]:>10.2f}{r[3][2][1]:>10.2f}" if r else ""
    print(line)

print("\n--- ORACLE-vs-ORACLE control: does the belief WIDTH alone move the prior? TV(orc, orcsharp) ---")
for a in ARMS:
    d = D[a]
    t = [compare(d, c, "orc", "orcsharp")[0] for c in CONDS if compare(d, c, "orc", "orcsharp")]
    print(f"{a:<10} TV mean {np.mean(t):.3f}  median {np.median(t):.3f}  max {np.max(t):.3f}")
