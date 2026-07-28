"""P1e (3c) — the prediction test with the RIGHT error metric per message TYPE.

A one-component gDNA-only claim does not deliver a composition at all — it delivers a LEVEL, ``mo_g``, which
literally asserts a gDNA mass SHARE  a_g = rho_g*E_g/M.  So score it on that: |clip(a_g,0,1) - oracle f_g|.
A lambda-emitting claim delivers a composition; score it on that.  Reported separately, plus the
"count-inadmissible" share (a_g > 1: more gDNA fragments claimed than the node sequenced in total).

    OMP_NUM_THREADS=1 python scratchpad/p1e_6_clean.py
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
P = {k: [] for k in ("err", "ag", "z2", "absd", "co", "M", "cls", "lam", "cond", "fo")}
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    live = (t["nsup"] > 0) & np.isfinite(t["fo"])
    d = t["dst"]
    ag = t["alpha_g"] * t["S"] / np.maximum(t["M"], 1e-9)  # exp(mo_g): the claimed gDNA MASS SHARE
    err = np.where(
        t["lam_emit"],
        np.abs(t["fg_msg"] - t["fo"]),
        np.abs(np.clip(ag, 0.0, 1.0) - t["fo"]),
    )
    for k, v in (("err", err), ("ag", ag), ("z2", t["z2"]), ("absd", np.abs(t["delta"])),
                 ("co", t["S"] / np.maximum(t["M"], 1e-9)), ("M", t["M"]), ("cls", t["cls"]),
                 ("lam", t["lam_emit"]), ("fo", t["fo"])):
        P[k].append(v[live])
    P["cond"].append(np.full(int(live.sum()), name))
P = {k: np.concatenate(v) for k, v in P.items()}
notig = P["cls"] != "intergenic"


def blk(tag, sel):
    e, z, dd, c, w = P["err"][sel], P["z2"][sel], P["absd"][sel], P["co"][sel], P["M"][sel]
    if e.size < 50:
        return
    qs = np.quantile(z, [0.2, 0.4, 0.6, 0.8])
    b = np.digitize(z, qs)
    print(f"\n  {tag}   (n = {e.size}, {w.sum():,.0f} reads)")
    print(f"    {'z2 quintile':<14}{'n':>7}{'z2 med':>11}{'|err| MEAN':>12}{'|err| med':>11}"
          f"{'mass-wtd |err|':>16}{'reads':>14}{'a_g>1':>9}")
    for i in range(5):
        m = b == i
        if not m.sum():
            continue
        print(f"    Q{i + 1:<13}{int(m.sum()):>7}{np.median(z[m]):>11.4f}{e[m].mean():>12.4f}"
              f"{np.median(e[m]):>11.4f}{np.average(e[m], weights=np.maximum(w[m], 1e-9)):>16.4f}"
              f"{w[m].sum():>14,.0f}{100 * np.mean(P['ag'][sel][m] > 1.0):>8.1f}%")
    print(f"    SPEARMAN vs |err|:  z2 {L.spearman(z, e):+.3f} | |delta| {L.spearman(dd, e):+.3f} "
          f"| claim/obs {L.spearman(c, e):+.3f}")


print("=" * 122)
print("POOLED, 7 conditions — the message's DELIVERED claim vs the destination's ORACLE f_g")
print("=" * 122)
blk("ALL live", np.ones(P["err"].size, bool))
blk("ALL live, EXCLUDING intergenic destinations", notig)
blk("lambda-EMITTING (composition delivered)", notig & P["lam"])
blk("ONE-COMPONENT claims (gDNA LEVEL delivered)", notig & ~P["lam"])
blk("ONE-COMPONENT x exon dst", notig & ~P["lam"] & (P["cls"] == "exon"))
blk("ONE-COMPONENT x boundary dst", notig & ~P["lam"] & (P["cls"] == "boundary"))
blk("ONE-COMPONENT x intron dst", notig & ~P["lam"] & (P["cls"] == "intron"))
print("\n  --- and the intergenic destinations on their own (the trivially-right gDNA anchors) ---")
blk("intergenic dst", ~notig)

print("\n" + "=" * 122)
print("THE COUNT-INADMISSIBLE CLAIM: a_g = rho_g*E_g/M > 1  ('more gDNA fragments than the node has')")
print("=" * 122)
print(f"  {'stratum':<44}{'n':>8}{'%a_g>1':>9}{'median a_g|a_g>1':>19}{'p99 a_g':>11}"
      f"{'reads on a_g>1':>17}{'mean|err| there':>17}")
for tag, sel in (("ALL live (excl. intergenic)", notig),
                 ("  one-component", notig & ~P["lam"]),
                 ("  lambda-emitting", notig & P["lam"]),
                 ("intergenic dst", ~notig)):
    a, e, w = P["ag"][sel], P["err"][sel], P["M"][sel]
    hi = a > 1.0
    print(f"  {tag:<44}{a.size:>8}{100 * hi.mean():>8.1f}%"
          f"{(np.median(a[hi]) if hi.any() else np.nan):>19.3f}{np.percentile(a, 99):>11.2f}"
          f"{w[hi].sum():>17,.0f}{(e[hi].mean() if hi.any() else np.nan):>17.4f}")
