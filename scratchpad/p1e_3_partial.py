"""P1e (3b) — IS z2 BETTER THAN |delta| ALONE?  The partial-rank test.

z2 = delta^2/(aSa + 1/n_dst) is |delta| divided by a DECLARED variance. The only question that matters is
whether the denominator adds information. So: stratify by |delta| quintile and ask whether z2 still separates
the error inside each stratum (and the mirror). If it does not, the simpler statistic wins.

    OMP_NUM_THREADS=1 python scratchpad/p1e_3_partial.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

ORDER = [
    "gdna300_ss0.99_present_capOFF",
    "gdna300_ss0.50_present_capOFF",
    "gdna300_ss0.99_present_capON",
    "gdna100_ss0.50_present_VERYSTRONG",
    "none_ss0.50_present_capOFF",
    "gdna100_ss0.50_none_capOFF",
    "gdna300_ss0.50_none_capOFF",
]
P = {k: [] for k in ("err", "z2", "absd", "co", "M", "cls", "lam", "cond", "denom")}
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, _ = L.message_table(inp, dbg)
    live = (t["nsup"] > 0) & np.isfinite(t["fg_msg"]) & np.isfinite(t["fo"])
    P["err"].append(np.abs(t["fg_msg"] - t["fo"])[live])
    P["z2"].append(t["z2"][live])
    P["absd"].append(np.abs(t["delta"])[live])
    P["co"].append((t["S"] / np.maximum(t["M"], 1e-9))[live])
    P["M"].append(t["M"][live])
    P["cls"].append(t["cls"][live])
    P["lam"].append(t["lam_emit"][live])
    P["denom"].append((t["aSa"] + 1.0 / np.maximum(t["n_dst"], 1e-9))[live])
    P["cond"].append(np.full(int(live.sum()), name))
P = {k: np.concatenate(v) for k, v in P.items()}


def wspearman(a, b, w):
    """mass-weighted rank correlation (ranks unweighted; the correlation weighted)."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b, w = a[m], b[m], w[m]
    if a.size < 5:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ma, mb = np.average(ra, weights=w), np.average(rb, weights=w)
    ca = np.average((ra - ma) * (rb - mb), weights=w)
    va = np.average((ra - ma) ** 2, weights=w)
    vb = np.average((rb - mb) ** 2, weights=w)
    return float(ca / np.sqrt(va * vb)) if va > 0 and vb > 0 else float("nan")


def block(tag, sel):
    e, z, d, w = P["err"][sel], P["z2"][sel], P["absd"][sel], P["M"][sel]
    if e.size < 100:
        return
    print(f"\n  {tag}  (n = {e.size}, {w.sum():,.0f} reads)")
    print(f"    plain Spearman vs |err| :  z2 {L.spearman(z, e):+.3f}   |delta| {L.spearman(d, e):+.3f}")
    print(f"    MASS-WEIGHTED Spearman  :  z2 {wspearman(z, e, w):+.3f}   |delta| {wspearman(d, e, w):+.3f}")
    dq = np.digitize(d, np.quantile(d, [0.2, 0.4, 0.6, 0.8]))
    zq = np.digitize(z, np.quantile(z, [0.2, 0.4, 0.6, 0.8]))
    print("    PARTIAL: within each |delta| quintile, Spearman(z2, |err|)  "
          "[if ~0, the declared variance adds nothing]")
    row = []
    for i in range(5):
        m = dq == i
        row.append(f"Q{i + 1}={L.spearman(z[m], e[m]):+.3f}(n={int(m.sum())})")
    print("      " + "  ".join(row))
    print("    PARTIAL: within each z2 quintile, Spearman(|delta|, |err|)")
    row = []
    for i in range(5):
        m = zq == i
        row.append(f"Q{i + 1}={L.spearman(d[m], e[m]):+.3f}(n={int(m.sum())})")
    print("      " + "  ".join(row))
    print("    mean |err| on the 5x5 grid  (rows = |delta| quintile, cols = z2 quintile)")
    print("        " + "".join(f"{'z2 Q' + str(j + 1):>12}" for j in range(5)))
    for i in range(5):
        cells = []
        for j in range(5):
            m = (dq == i) & (zq == j)
            cells.append(f"{e[m].mean():>8.3f}({m.sum():>3d})" if m.sum() else f"{'-':>12}")
        print(f"    d Q{i + 1} " + "".join(f"{c:>12}" for c in cells))


print("=" * 118)
print("POOLED over 7 conditions — is the DECLARED VARIANCE (the denominator) carrying information?")
print("=" * 118)
block("ALL live messages", np.ones(P["err"].size, bool))
block("lambda-EMITTING", P["lam"])
block("lambda-EMITTING x exon (the graft)", P["lam"] & (P["cls"] == "exon"))
block("lambda-EMITTING x boundary", P["lam"] & (P["cls"] == "boundary"))
block("lambda-EMITTING x intron", P["lam"] & (P["cls"] == "intron"))
for c in ORDER:
    block(f"lambda-EMITTING @ {c}", P["lam"] & (P["cond"] == c))
