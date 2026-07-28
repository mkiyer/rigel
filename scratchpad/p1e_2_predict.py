"""P1e (3) — DOES z2 PREDICT THE MESSAGE'S ACTUAL ERROR?  (needs the oracle; this is the case for P1e)

For every live message: its DELIVERED composition  f_g^msg = rho_g E_g / (rho_g E_g + rho_R E_r)  (exactly
the sigmoid of the lambda it emits, and invariant under the common-factor pin) against the destination
node's ORACLE f_g. Binned by z2 quintile; Spearman against the three candidate statistics:

    z2 = delta^2/(aSa + 1/n_dst)          the surprise
    |delta| = |log(M/S)|                  the raw violation  (monotone-equivalent to |log(claim/obs)|)
    claim/obs = S/M                       PASS0_FINISH_PLAN §P1b's simpler statistic (SIGNED)

    OMP_NUM_THREADS=1 python scratchpad/p1e_2_predict.py
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

POOL = {k: [] for k in ("err", "z2", "absd", "co", "M", "cls", "lam", "cond")}

for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    live = (t["nsup"] > 0) & np.isfinite(t["fg_msg"]) & np.isfinite(t["fo"])
    err = np.abs(t["fg_msg"] - t["fo"])
    co = t["S"] / np.maximum(t["M"], 1e-9)
    absd = np.abs(t["delta"])
    print(f"\n{'=' * 118}\n{name}\n{'=' * 118}")

    for tag, sel in (
        ("ALL live messages", live),
        ("lambda-EMITTING only (both components supplied)", live & t["lam_emit"]),
    ):
        e, z, d, c, w = err[sel], t["z2"][sel], absd[sel], co[sel], t["M"][sel]
        if e.size < 20:
            continue
        qs = np.quantile(z, [0.2, 0.4, 0.6, 0.8])
        b = np.digitize(z, qs)
        print(f"\n  {tag}   (n = {e.size})")
        print(f"    {'z2 quintile':<16}{'n':>6}{'z2 med':>10}{'z2 max':>11}"
              f"{'|err| MEAN':>12}{'|err| med':>11}{'mass-wtd |err|':>16}{'reads':>14}")
        for i in range(5):
            m = b == i
            if m.sum() == 0:
                continue
            print(f"    Q{i + 1:<15}{int(m.sum()):>6}{np.median(z[m]):>10.4f}{z[m].max():>11.1f}"
                  f"{e[m].mean():>12.4f}{np.median(e[m]):>11.4f}"
                  f"{np.average(e[m], weights=np.maximum(w[m], 1e-9)):>16.4f}{w[m].sum():>14,.0f}")
        print(f"    SPEARMAN vs |err|:   z2 = {L.spearman(z, e):+.3f}   "
              f"|delta| = {L.spearman(d, e):+.3f}   claim/obs (signed) = {L.spearman(c, e):+.3f}   "
              f"|log claim/obs| = {L.spearman(np.abs(np.log(np.maximum(c, 1e-12))), e):+.3f}")
        # and the same for claim/obs bins, to reproduce P1b's table shape
        cq = np.quantile(c, [0.2, 0.4, 0.6, 0.8])
        cb = np.digitize(c, cq)
        print(f"    {'claim/obs bin':<16}{'n':>6}{'c/o med':>10}{'':>11}"
              f"{'|err| MEAN':>12}{'|err| med':>11}{'mass-wtd |err|':>16}")
        for i in range(5):
            m = cb == i
            if m.sum() == 0:
                continue
            print(f"    Q{i + 1:<15}{int(m.sum()):>6}{np.median(c[m]):>10.4f}{'':>11}"
                  f"{e[m].mean():>12.4f}{np.median(e[m]):>11.4f}"
                  f"{np.average(e[m], weights=np.maximum(w[m], 1e-9)):>16.4f}")

    for k, v in (("err", err[live]), ("z2", t["z2"][live]), ("absd", absd[live]),
                 ("co", co[live]), ("M", t["M"][live]), ("cls", t["cls"][live]),
                 ("lam", t["lam_emit"][live])):
        POOL[k].append(v)
    POOL["cond"].append(np.full(int(live.sum()), name))

P = {k: np.concatenate(v) for k, v in POOL.items()}
print(f"\n{'=' * 118}\nPOOLED over the 7 conditions   (n = {P['err'].size})\n{'=' * 118}")
for tag, sel in (("ALL live", np.ones(P["err"].size, bool)),
                 ("lambda-emitting", P["lam"]),
                 ("exon (= graft) only", P["cls"] == "exon"),
                 ("boundary only", P["cls"] == "boundary"),
                 ("intron only", P["cls"] == "intron")):
    e, z, d, c, w = (P["err"][sel], P["z2"][sel], P["absd"][sel], P["co"][sel], P["M"][sel])
    if e.size < 20:
        continue
    qs = np.quantile(z, [0.2, 0.4, 0.6, 0.8])
    b = np.digitize(z, qs)
    print(f"\n  {tag}  (n = {e.size})")
    print(f"    {'z2 quintile':<16}{'n':>7}{'z2 med':>10}{'|err| MEAN':>12}{'mass-wtd |err|':>16}")
    for i in range(5):
        m = b == i
        print(f"    Q{i + 1:<15}{int(m.sum()):>7}{np.median(z[m]):>10.4f}{e[m].mean():>12.4f}"
              f"{np.average(e[m], weights=np.maximum(w[m], 1e-9)):>16.4f}")
    print(f"    SPEARMAN vs |err|:   z2 = {L.spearman(z, e):+.3f}   |delta| = {L.spearman(d, e):+.3f}   "
          f"claim/obs = {L.spearman(c, e):+.3f}")
