"""SINGLE-STRAND x CAPTURE, step 2 — inside the EXONS: mode or precision, and which channel?

Step 1 located it: exons carry 77-92 % of the capture-ON single-strand error and the driver is the MESSAGES
(Δmsg +0.0005 off -> +0.052 on). In the nrna_none pair the exon SELF-solve is bit-identical off/on
(0.0049 -> 0.0049) while SOLVED degrades 4x. So this is a pure message defect on a clean control.

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_2_exons.py
"""

from __future__ import annotations

import numpy as np

d = np.load("/tmp/suite_nodes.npz", allow_pickle=True)
cond = d["cond"]
PAIRS = [
    ("gdna_gdna300_ss_0.99_nrna_none_capture_off", "gdna_gdna300_ss_0.99_nrna_none_capture_on"),
    ("gdna_gdna300_ss_0.99_nrna_present_capture_off", "gdna_gdna300_ss_0.99_nrna_present_capture_on"),
    ("gdna_gdna100_ss_0.99_nrna_present_capture_off", "gdna_gdna100_ss_0.99_nrna_present_capture_on"),
]


def sub(c, cls="exon"):
    m = (cond == c) & (d["dof"] != "ambig") & (d["cls"] == cls)
    return {k: d[k][m] for k in d.files}


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print("EXONS, single-strand.  err = |solved-oracle|.  DIRECTION: signed (solved-oracle) mass-weighted.")
print(f"{'condition':<42}{'oracle':>8}{'self':>8}{'solved':>8}{'signed':>9}{'tau_own':>9}"
      f"{'c_tau':>9}{'lamfg':>8}{'cm_g':>9}{'mo_g':>8}{'cm_p':>9}{'mo_p':>8}{'cm_n':>9}{'mo_n':>8}")
for coff, con in PAIRS:
    for c in (coff, con):
        s = sub(c)
        w = s["mass"]
        print(f"  {c[5:]:<40}{mw(s['oracle'], w):>8.4f}{mw(s['self'], w):>8.4f}{mw(s['solved'], w):>8.4f}"
              f"{mw(s['solved'] - s['oracle'], w):>+9.4f}{mw(s['tau_own'], w):>9.2f}{mw(s['c_tau'], w):>9.2f}"
              f"{mw(s['lam_fg'], w):>8.4f}{mw(s['cm_g'], w):>9.2f}{mw(np.clip(s['mo_g'], 0, 2), w):>8.4f}"
              f"{mw(s['cm_p'], w):>9.2f}{mw(np.clip(s['mo_p'], 0, 2), w):>8.4f}"
              f"{mw(s['cm_n'], w):>9.2f}{mw(np.clip(s['mo_n'], 0, 2), w):>8.4f}")
    print()

# ── the worst exons under capture, per pair: what happened to each channel ────────────────────────────────
print("=" * 130)
print("WORST EXONS by err-mass (capture ON), with their OFF twin")
print("=" * 130)
for coff, con in PAIRS[:2]:
    a, b = sub(coff), sub(con)
    ia = {int(n): i for i, n in enumerate(a["node"])}
    err = np.abs(b["solved"] - b["oracle"])
    em = err * b["mass"]
    top = np.argsort(-em)[:12]
    print(f"\n{con[5:]}")
    print(f"  {'node':>6}{'mass':>10}{'orc':>7}{'self':>7}{'solv':>7}|{'OFFsolv':>8}{'OFForc':>7}|"
          f"{'tau_own':>8}{'c_tau':>8}{'lamfg':>7}|{'cm_g':>8}{'mo_g':>7}|{'cm_p':>8}{'mo_p':>7}|"
          f"{'cm_n':>8}{'mo_n':>7}|{'nL':>6}{'nR':>6}")
    for i in top:
        nd = int(b["node"][i])
        j = ia.get(nd, -1)
        offs = f"{a['solved'][j]:>8.3f}{a['oracle'][j]:>7.3f}" if j >= 0 else f"{'-':>8}{'-':>7}"
        print(f"  {nd:>6}{b['mass'][i]:>10,.0f}{b['oracle'][i]:>7.3f}{b['self'][i]:>7.3f}"
              f"{b['solved'][i]:>7.3f}|{offs}|{b['tau_own'][i]:>8.2f}{b['c_tau'][i]:>8.2f}"
              f"{b['lam_fg'][i]:>7.3f}|{b['cm_g'][i]:>8.2f}{min(b['mo_g'][i], 9.99):>7.3f}|"
              f"{b['cm_p'][i]:>8.2f}{min(b['mo_p'][i], 9.99):>7.3f}|{b['cm_n'][i]:>8.2f}"
              f"{min(b['mo_n'][i], 9.99):>7.3f}|{b['nl_cls'][i][:5]:>6}{b['nr_cls'][i][:5]:>6}")
