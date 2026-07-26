"""Part 2 of the adversarial verification. Three decisive tests x3 did not run.

  A. VARIANCE COMPONENTS.  x3's omega_hat = Var(pooled rows) - E[per-row floor] is a law-of-total-variance
     estimator.  Its premise is testable directly: split the pooled variance into BETWEEN-edge and
     WITHIN-edge (across conditions) parts.  The WITHIN part is pure counting noise (different simulated
     libraries, disjoint fragment sets), so it must equal E[floor].  If E[floor] > Var_within the floor
     is over-stated; if <, under-stated.  And Var_between - E[floor]/n_cond is the honest omega.
  B. GTF gene grouping: x3 groups isoforms by `transcript_id.split(".")[0]`.  Check against gene_id.
  C. The 2-D cross-tab cells x3 reads as "expression raises the variance": what is in them?
"""

from __future__ import annotations

import collections
import dataclasses
import importlib
import re
import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
GTF = SUITE / "reference" / "genes.gtf"
_EPS = 1e-9

# ══ B. GTF — independent check of the alternative-splicing census ═════════════════════════════════
g2t, t2g, tx = collections.defaultdict(set), {}, collections.defaultdict(list)
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    gid = re.search(r'gene_id "([^"]+)"', f[8])
    tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    if gid:
        g2t[gid.group(1)].add(tid)
        t2g[tid] = gid.group(1)
    tx[tid].append((int(f[3]) - 1, int(f[4])))
TX = {t: sorted(v) for t, v in tx.items()}
print("=" * 112)
print("B. GTF CENSUS — independent, using the real gene_id (x3 used transcript_id.split('.')[0])")
print("=" * 112)
iso = collections.Counter(len(v) for v in g2t.values())
print(f"  by gene_id:   genes={len(g2t)}  transcripts={len(TX)}  isoforms/gene: "
      + ", ".join(f"{k}:{v}" for k, v in sorted(iso.items())))
iso2 = collections.Counter()
_p = collections.defaultdict(set)
for t in TX:
    _p[t.split(".")[0]].add(t)
for v in _p.values():
    iso2[len(v)] += 1
print(f"  by x3's rule: genes={len(_p)}  isoforms/gene: "
      + ", ".join(f"{k}:{v}" for k, v in sorted(iso2.items())))
# do multi-isoform genes actually differ in structure?
n_same, n_diff = 0, 0
for g, ts in g2t.items():
    if len(ts) < 2:
        continue
    structs = {tuple(TX[t]) for t in ts}
    intr = {tuple((a[1], b[0]) for a, b in zip(TX[t], TX[t][1:])) for t in ts}
    n_same += len(structs) == 1
    n_diff += len(intr) > 1
print(f"  multi-isoform genes: {sum(1 for v in g2t.values() if len(v) > 1)};  "
      f"identical exon structure: {n_same};  genes whose isoforms differ in their INTRON SET: {n_diff}")
JP = collections.Counter()
for t, ex in TX.items():
    for a, b in zip(ex, ex[1:]):
        JP[a[1]] += 1
        JP[b[0]] += 1
# a coordinate shared by 2 tx OF THE SAME GENE = real alternative splicing; across genes = overlap
JG = collections.defaultdict(set)
for t, ex in TX.items():
    for a, b in zip(ex, ex[1:]):
        JG[a[1]].add(t2g.get(t, t))
        JG[b[0]].add(t2g.get(t, t))
print(f"  distinct junction coordinates = {len(JP)}; shared by >=2 transcripts = "
      f"{sum(1 for v in JP.values() if v > 1)}; of those, spanning >1 GENE = "
      f"{sum(1 for k, v in JP.items() if v > 1 and len(JG[k]) > 1)}")
print()

# ══ the measurement (same definitions as x3/x5) ══════════════════════════════════════════════════
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
R_START = np.asarray(index.region_df["start"], np.int64)
R_END = np.asarray(index.region_df["end"], np.int64)
N_REG = R_START.shape[0]
B_POS = np.array([int(R_START[b]) if b < N_REG else int(R_END[N_REG - 1]) for b in range(N_REG + 1)])
TT = collections.Counter()
for t, ex in TX.items():
    TT[ex[0][0]] += 1
    TT[ex[-1][1]] += 1
B_NTERM = np.array([TT.get(int(p), 0) for p in B_POS], np.int64)
B_NJUNC = np.array([JP.get(int(p), 0) for p in B_POS], np.int64)

CONDS = ["gdna_gdna300_ss_0.99_nrna_present_capture_off",
         "gdna_gdna300_ss_0.50_nrna_present_capture_off",
         "gdna_gdna100_ss_0.99_nrna_present_capture_off",
         "gdna_gdna100_ss_0.50_nrna_present_capture_off",
         "gdna_gdna300_ss_0.99_nrna_none_capture_off",
         "gdna_gdna300_ss_0.50_nrna_none_capture_off",
         "gdna_gdna100_ss_0.99_nrna_none_capture_off",
         "gdna_gdna100_ss_0.50_nrna_none_capture_off"]
F = ["ci", "ex", "bn", "face", "phi", "s", "s_nu", "n_R", "n_spl", "k_nu",
     "G_b", "G_e", "nterm", "njunc"]


def measure(cond, ci):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                          SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind, idx = np.asarray(chain.kind), np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = (pool("mat_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_pos") + pool("nas_uns_neg"))
    Rs = pool("mat_spl") + pool("nas_spl")
    E_g, E_r = us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    NSP = (us["spl_n_pos_l"] + us["spl_n_neg_l"], us["spl_n_pos_r"] + us["spl_n_neg_r"])
    NUN = (us["n_unspl_l"], us["n_unspl_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)
    fR = np.where((Ru + G) > _EPS, Ru / np.maximum(Ru + G, _EPS), 0.0)
    rows = []
    for face, nbr in ((1, li), (0, ri)):
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not (mu > _EPS) or not np.isfinite(rho_R[i]) or rho_R[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS):
                continue
            den = rho_R[i] * (rho_g[b] / rho_g[i])
            rows.append((ci, i, b, face, (rho_nu[b] + mu) / den, mu / den, rho_nu[b] / den,
                         Ru[i] + Rs[i], NSP[face][b], (NUN[0][b] + NUN[1][b]) * fR[b],
                         G[b], G[i], B_NTERM[idx[b]], B_NJUNC[idx[b]]))
    return rows


rows = []
for ci, c in enumerate(CONDS):
    rows += measure(c, ci)
A = np.array(rows, float)
C = {f: A[:, k] for k, f in enumerate(F)}
ok = np.isfinite(C["phi"]) & (C["phi"] > 0) & (C["s"] > 0)
C = {f: C[f][ok] for f in F}
lp = np.log(C["phi"])
w_mu = C["s"] / np.maximum(C["phi"], _EPS)
w_nu = C["s_nu"] / np.maximum(C["phi"], _EPS)
fl3 = (w_nu ** 2 * polygamma(1, np.maximum(C["k_nu"], 0.5))
       + w_mu ** 2 * polygamma(1, np.maximum(C["n_spl"], 0.5))
       + polygamma(1, np.maximum(C["n_R"], 0.5)))
flS = fl3 + polygamma(1, np.maximum(C["G_b"], 0.5)) + polygamma(1, np.maximum(C["G_e"], 0.5))
term = C["nterm"] > 0

print("=" * 112)
print("A. VARIANCE COMPONENTS on the 8 clean conditions (gdna100+gdna300, capture-OFF).")
print("   Var_within = the SAME edge across independent simulated libraries = PURE counting noise.")
print("   If the floor is right, E[floor] == Var_within.  omega = Var_between - E[floor]/n_cond.")
print("=" * 112)
key = (C["ex"] * 1e7 + C["bn"] * 10 + C["face"]).astype(np.int64)
print(f"  {'stratum':<22}{'rows':>6}{'edges':>7}{'Var_pool':>10}{'Var_betw':>10}{'Var_with':>10}"
      f"{'E[fl_x3]':>10}{'E[fl+step]':>12}{'omega_x3':>10}{'omega_VC':>10}")
for lab, m0 in (("ALL", np.ones_like(term, bool)), ("TERMINUS", term), ("JUNCTION-only", ~term)):
    uk, inv = np.unique(key[m0], return_inverse=True)
    lm = lp[m0]
    means = np.array([lm[inv == k].mean() for k in range(uk.shape[0])])
    nk = np.array([int((inv == k).sum()) for k in range(uk.shape[0])])
    wv = np.array([lm[inv == k].var(ddof=1) if (inv == k).sum() > 1 else np.nan
                   for k in range(uk.shape[0])])
    vb = float(np.var(means))
    vw = float(np.nanmean(wv))
    nc = float(np.mean(nk))
    f3, fs = float(np.mean(fl3[m0])), float(np.mean(flS[m0]))
    print(f"  {lab:<22}{int(m0.sum()):>6}{uk.shape[0]:>7}{np.var(lm):>10.4f}{vb:>10.4f}"
          f"{vw:>10.4f}{f3:>10.4f}{fs:>12.4f}{max(0, np.var(lm) - f3):>10.4f}"
          f"{max(0.0, vb - fs / nc):>10.4f}")
print()
print("  READ: Var_within is what the counting noise ACTUALLY is (it is measured, not modelled).")
print("  Compare it to E[floor]. Then omega_VC = Var_between - E[floor]/n_cond is the honest premise")
print("  variance, and it does not depend on the floor being right except through a /n_cond term.")
print()
print("  the SAME, using the MEASURED Var_within in place of the modelled floor:")
for lab, m0 in (("ALL", np.ones_like(term, bool)), ("TERMINUS", term), ("JUNCTION-only", ~term)):
    uk, inv = np.unique(key[m0], return_inverse=True)
    lm = lp[m0]
    means = np.array([lm[inv == k].mean() for k in range(uk.shape[0])])
    nk = np.array([int((inv == k).sum()) for k in range(uk.shape[0])])
    wv = np.array([lm[inv == k].var(ddof=1) if (inv == k).sum() > 1 else np.nan
                   for k in range(uk.shape[0])])
    vw = float(np.nanmean(wv))
    om = max(0.0, float(np.var(means)) - vw / float(np.mean(nk)))
    print(f"    {lab:<22} omega (measured-noise version) = {om:.4f}")
print()

# ══ C. what is inside the cross-tab cells x3 reads as an expression effect? ═══════════════════════
print("=" * 112)
print("C. THE 2-D CROSS-TAB CELLS x3 READS AS 'EXPRESSION RAISES THE VARIANCE'")
print("=" * 112)


def qb(x, q):
    e = np.unique(np.quantile(x, np.linspace(0, 1, q + 1)))
    return np.clip(np.searchsorted(e, x, "right") - 1, 0, e.shape[0] - 2)


bs, br = qb(C["n_spl"], 4), qb(C["n_R"], 4)
print(f"  {'cell':<14}{'n':>6}{'Var(log phi)':>14}{'med log phi':>13}{'% TERMINUS':>12}"
      f"{'med n_spl':>11}{'med n_R':>10}")
for i in range(4):
    for j in range(4):
        m = (bs == i) & (br == j)
        if int(m.sum()) < 8:
            continue
        print(f"  nspl Q{i+1} nR Q{j+1}".ljust(14)
              + f"{int(m.sum()):>6}{np.var(lp[m]):>14.3f}{np.median(lp[m]):>13.3f}"
              f"{term[m].mean():>12.1%}{np.median(C['n_spl'][m]):>11.0f}"
              f"{np.median(C['n_R'][m]):>10.0f}")
print()
print("  the SAME cross-tab restricted to JUNCTION-only boundaries (the terminus class removed):")
print("  " + "n_spl\\n_R".ljust(14) + "".join(f"{f'nR Q{j+1}':>14}" for j in range(4)) + f"{'ROW':>14}")
for i in range(4):
    row = f"  Q{i+1}".ljust(14)
    for j in range(4):
        m = (bs == i) & (br == j) & (~term)
        row += (f"{np.var(lp[m]):>8.3f}/{int(m.sum()):<5d}" if int(m.sum()) >= 8 else f"{'-':>14}")
    m = (bs == i) & (~term)
    row += f"{np.var(lp[m]):>8.3f}/{int(m.sum()):<5d}"
    print(row)

# ══ D. THE EXCESS TABLES, RECOMPUTED WITH THE FLOOR THAT MATCHES THE MEASURED NOISE ═══════════════
print()
print("=" * 112)
print("D. x3's 'irreducible premise error' AND 'count-flatness', with the step term restored.")
print("   (E[fl+step] was just shown to match the MEASURED counting noise Var_within to <6%.)")
print("=" * 112)
print(f"  {'stratum / n_spl bin':<30}{'n':>6}{'Var':>9}{'fl_x3':>9}{'excess_x3':>11}"
      f"{'fl+step':>10}{'excess_corr':>13}")
for lab, m0 in (("JUNCTION-only", ~term), ("TERMINUS", term)):
    for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
        m = m0 & (C["n_spl"] >= lo) & (C["n_spl"] < hi)
        if int(m.sum()) < 15:
            continue
        v = float(np.var(lp[m]))
        f3, fs = float(np.mean(fl3[m])), float(np.mean(flS[m]))
        print(f"  {lab + ' [' + str(lo) + ',' + (str(hi) if hi < 1e11 else 'inf') + ')':<30}"
              f"{int(m.sum()):>6}{v:>9.3f}{f3:>9.4f}{v - f3:>11.3f}{fs:>10.4f}{v - fs:>13.3f}")
# class D equivalent: junction-only AND single-isoform exon region
print()
print("  and the 'structurally perfect' subset (junction-only boundary, n_spl bins):")
print(f"  the pooled JUNCTION-only omega:  x3 = {np.var(lp[~term]) - np.mean(fl3[~term]):.4f}   "
      f"corrected = {max(0.0, np.var(lp[~term]) - np.mean(flS[~term])):.4f}   "
      f"variance-components = 0.0610")
