"""P1e (3d) — the UNSATURATED target: the gDNA LEVEL error in log space.

The clip-to-[0,1] composition metric saturates on the one-component claims (a_g > 1 on 100 % of them at
exon destinations), so |err| there is a constant of the node and neither statistic can be scored on it. The
honest, unsaturated target for a level claim is

    L_g = |log( rho_g^claim / rho_g^oracle )|        (and the same for the RNA arm)

    OMP_NUM_THREADS=1 python scratchpad/p1e_7_level.py
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
_EPS = 1e-9
P = {k: [] for k in ("lg", "lr", "z2", "absd", "M", "cls", "lam", "cond", "supg", "supr")}
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    d = t["dst"]
    M, E_g, E_r, fo = nf["M"], nf["E_g"], nf["E_r"], nf["fo"]
    og_t, or_t = fo * M / np.maximum(E_g, _EPS), (1 - fo) * M / np.maximum(E_r, _EPS)
    cg = t["alpha_g"] * t["S"] / np.maximum(E_g[d], _EPS)
    cr = (t["alpha_p"] + t["alpha_n"]) * t["S"] / np.maximum(E_r[d], _EPS)
    live = (t["nsup"] > 0) & np.isfinite(fo[d])
    lg = np.where(t["sup_g"] & (og_t[d] > _EPS) & (cg > _EPS),
                  np.abs(np.log(np.maximum(cg, _EPS) / np.maximum(og_t[d], _EPS))), np.nan)
    lr = np.where(t["sup_r"] & (or_t[d] > _EPS) & (cr > _EPS),
                  np.abs(np.log(np.maximum(cr, _EPS) / np.maximum(or_t[d], _EPS))), np.nan)
    for k, v in (("lg", lg), ("lr", lr), ("z2", t["z2"]), ("absd", np.abs(t["delta"])),
                 ("M", t["M"]), ("cls", t["cls"]), ("lam", t["lam_emit"]),
                 ("supg", t["sup_g"]), ("supr", t["sup_r"])):
        P[k].append(v[live])
    P["cond"].append(np.full(int(live.sum()), name))
P = {k: np.concatenate(v) for k, v in P.items()}


def blk(tag, sel, tgt):
    e, z, dd, w = P[tgt][sel], P["z2"][sel], P["absd"][sel], P["M"][sel]
    m = np.isfinite(e)
    e, z, dd, w = e[m], z[m], dd[m], w[m]
    if e.size < 50:
        return
    b = np.digitize(z, np.quantile(z, [0.2, 0.4, 0.6, 0.8]))
    print(f"\n  {tag}  [target = {tgt}]  (n = {e.size}, {w.sum():,.0f} reads)")
    print(f"    {'z2 quintile':<14}{'n':>7}{'z2 med':>11}{'MEAN |log err|':>16}{'median':>10}"
          f"{'mass-wtd':>12}{'x fold':>10}")
    for i in range(5):
        k = b == i
        if not k.sum():
            continue
        print(f"    Q{i + 1:<13}{int(k.sum()):>7}{np.median(z[k]):>11.4f}{e[k].mean():>16.4f}"
              f"{np.median(e[k]):>10.4f}{np.average(e[k], weights=np.maximum(w[k], 1e-9)):>12.4f}"
              f"{np.exp(np.average(e[k], weights=np.maximum(w[k], 1e-9))):>9.2f}x")
    print(f"    SPEARMAN:  z2 {L.spearman(z, e):+.3f}   |delta| {L.spearman(dd, e):+.3f}"
          f"   [n.b. for one-component claims a_g is a deterministic function of delta]")


print("=" * 118)
print("POOLED, 7 conditions — the gDNA LEVEL error |log(claimed/oracle)| and the RNA arm's")
print("=" * 118)
notig = P["cls"] != "intergenic"
for tgt in ("lg", "lr"):
    blk("ALL supplying that arm (excl. intergenic)", notig, tgt)
    blk("  lambda-emitting", notig & P["lam"], tgt)
    blk("  one-component", notig & ~P["lam"], tgt)
    blk("  exon destinations", notig & (P["cls"] == "exon"), tgt)
    blk("  boundary destinations", notig & (P["cls"] == "boundary"), tgt)
    blk("  intron destinations", notig & (P["cls"] == "intron"), tgt)
