"""SINGLE-STRAND x CAPTURE, step 1 — decompose the 10x degradation: self-solve? message MODE? message PRECISION?

HANDOFF_7 §4-5. Matched capture-off/on pairs, restricted to single-strand nodes (dof != 'ambig'), joined
PER NODE (the region partition is identical across conditions, so node ids match).

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_1_decompose.py
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


def sub(c):
    m = (cond == c) & (d["dof"] != "ambig")
    return {k: d[k][m] for k in d.files}


def mw(x, w):
    return np.average(x, weights=w) if w.sum() > 0 else np.nan


print(f"{'pair / stratum':<46}{'nodes':>7}{'mass':>12}{'self':>8}{'solved':>8}{'Δmsg':>8}"
      f"{'mo_g err':>10}{'lam err':>9}{'cm_g':>10}{'c_tau':>10}")
for coff, con in PAIRS:
    for lab, c in (("OFF", coff), ("ON ", con)):
        s = sub(c)
        w = s["mass"]
        se = np.abs(s["self"] - s["oracle"])
        de = np.abs(s["solved"] - s["oracle"])
        # message MODE error: the gDNA claim delivered to psi (exp(mo_g) is a fraction) and the lambda claim
        mge = np.abs(np.clip(s["mo_g"], 0, 1) - s["oracle"])
        lme = np.abs(s["lam_fg"] - s["oracle"])
        # precisions: only where the channel is live
        cg, ct = s["cm_g"], s["c_tau"]
        print(f"  {c[5:]:<44}{len(w):>7}{w.sum():>12,.0f}{mw(se, w):>8.4f}{mw(de, w):>8.4f}"
              f"{mw(de, w) - mw(se, w):>+8.4f}{mw(mge, w):>10.4f}{mw(lme, w):>9.4f}"
              f"{mw(cg, w):>10.1f}{mw(ct, w):>10.1f}")
    print()

# ── per-node paired view: which nodes degrade, and does their SELF degrade too? ────────────────────────────
print("\n" + "=" * 118)
print("PER-NODE PAIRED (same node id, capture off vs on) — single-strand only")
print("=" * 118)
for coff, con in PAIRS:
    a, b = sub(coff), sub(con)
    ia = {int(n): i for i, n in enumerate(a["node"])}
    common = np.array([int(n) for n in b["node"] if int(n) in ia])
    bi = {int(n): i for i, n in enumerate(b["node"])}
    ai = np.array([ia[n] for n in common])
    bj = np.array([bi[n] for n in common])
    wa, wb = a["mass"][ai], b["mass"][bj]
    sea, seb = np.abs(a["self"][ai] - a["oracle"][ai]), np.abs(b["self"][bj] - b["oracle"][bj])
    dea, deb = np.abs(a["solved"][ai] - a["oracle"][ai]), np.abs(b["solved"][bj] - b["oracle"][bj])
    print(f"\n{con[5:]}   ({len(common)} common single-strand nodes)")
    print(f"  self   OFF {mw(sea, wa):.4f}  ON {mw(seb, wb):.4f}    "
          f"solved OFF {mw(dea, wa):.4f}  ON {mw(deb, wb):.4f}")
    # split by region class (coarse) — where does the degradation land?
    print(f"    {'cls':<12}{'n':>7}{'massON':>12}{'selfOFF':>9}{'selfON':>9}{'solvOFF':>9}{'solvON':>9}"
          f"{'errmassON':>12}{'share':>7}")
    emb = deb * wb
    for cl in ("exon", "intron", "boundary", "intergenic"):
        m = b["cls"][bj] == cl
        if not m.any():
            continue
        print(f"    {cl:<12}{int(m.sum()):>7}{wb[m].sum():>12,.0f}{mw(sea[m], wa[m]):>9.4f}"
              f"{mw(seb[m], wb[m]):>9.4f}{mw(dea[m], wa[m]):>9.4f}{mw(deb[m], wb[m]):>9.4f}"
              f"{emb[m].sum():>12,.0f}{emb[m].sum() / emb.sum():>7.1%}")
