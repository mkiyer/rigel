"""ADVERSARIAL VERIFICATION of scratchpad/x3_share.py's "share-vs-count" claim set.

Re-measures phi from scratch (same definitions as x3) and adds the three things x3's floor omits:

  1. THE ORACLE STEP'S OWN NOISE.  phi's denominator carries `step = rho_g(bnd)/rho_g(exon)`, a ratio of
     two COUNTED gDNA pools.  x3's "trigamma floor" has NO gDNA term at all -- only k_nu, k_mu, k_R.
     Var(log step) ~= psi'(G_bnd) + psi'(G_exon).  Adding it changes every omega_hat.
  2. THE DECOMPOSITION OF THE FLOOR.  x3 compares its 3-term trigamma floor against a 1-term `1/n_spl`
     and attributes the gap to "trigamma vs 1/n".  Decompose it.
  3. LOCATION vs SCALE in the dispersion regression.  y = log|log phi - GLOBAL median| conflates a
     group MEAN shift with a group SPREAD.  Re-run centred within cells.

    OMP_NUM_THREADS=1 python scratchpad/x5_verify.py
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

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna5_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_off",
]
N_PRESENT = 8

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
rdf = index.region_df
R_START = np.asarray(rdf["start"], np.int64)
R_END = np.asarray(rdf["end"], np.int64)
R_LEN = R_END - R_START
N_REG = R_START.shape[0]

# ── GTF ──
tx = collections.defaultdict(list)
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    tx[tid].append((int(f[3]) - 1, int(f[4])))
TX = {t: sorted(v) for t, v in tx.items()}
JUNC_POS, TSS_TES = collections.Counter(), collections.Counter()
EXON_IV = []
for t, ex in TX.items():
    TSS_TES[ex[0][0]] += 1
    TSS_TES[ex[-1][1]] += 1
    for a, b in zip(ex, ex[1:]):
        JUNC_POS[a[1]] += 1
        JUNC_POS[b[0]] += 1
    for e in ex:
        EXON_IV.append((e[0], e[1], t))
_ex_s = np.array([e[0] for e in EXON_IV])
_ex_e = np.array([e[1] for e in EXON_IV])
_ex_t = [e[2] for e in EXON_IV]
reg_n_tx = np.zeros(N_REG, np.int64)
for r in range(N_REG):
    hit = (_ex_s < R_END[r]) & (_ex_e > R_START[r])
    reg_n_tx[r] = len({_ex_t[i] for i in np.flatnonzero(hit)})
B_POS = np.array([int(R_START[b]) if b < N_REG else int(R_END[N_REG - 1]) for b in range(N_REG + 1)])
B_NJUNC = np.array([JUNC_POS.get(int(p), 0) for p in B_POS], np.int64)
B_NTERM = np.array([TSS_TES.get(int(p), 0) for p in B_POS], np.int64)

FIELDS = ["cond_i", "exon_node", "bnd_node", "face", "reg_idx", "bnd_idx",
          "phi", "s", "s_nu", "n_R", "n_spl", "k_nu", "k_R",
          "G_bnd", "G_exon", "n_uns_bnd", "step", "b_njunc", "b_nterm", "reg_n_tx", "nrna",
          "Ru_bnd", "Rs_ex", "Ru_ex"]


def measure(cond, ci):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                          SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
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
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)
    f_R_bnd = np.where((Ru + G) > _EPS, Ru / np.maximum(Ru + G, _EPS), 0.0)
    nrna = 1.0 if "nrna_present" in cond else 0.0

    rows = []
    for face, nbr in ((1, li), (0, ri)):
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not (mu > _EPS) or not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS):
                continue
            step = rho_g[b] / rho_g[i]
            den = rho_R_ex[i] * step
            rows.append((ci, i, b, face, idx[i], idx[b],
                         (rho_nu_b[b] + mu) / den, mu / den, rho_nu_b[b] / den,
                         Ru[i] + Rs[i], NSP[face][b],
                         (NUN[0][b] + NUN[1][b]) * f_R_bnd[b], Ru[i] + Rs[i],
                         G[b], G[i], NUN[0][b] + NUN[1][b], step,
                         B_NJUNC[idx[b]], B_NTERM[idx[b]], reg_n_tx[idx[i]], nrna,
                         Ru[b], Rs[i], Ru[i]))
    return rows


all_rows = []
for ci, cond in enumerate(CONDS):
    all_rows += measure(cond, ci)
A = np.array(all_rows, float)
C = {f: A[:, k] for k, f in enumerate(FIELDS)}
ok = np.isfinite(C["phi"]) & (C["phi"] > 0) & (C["s"] > 0)
C = {f: C[f][ok] for f in FIELDS}
lp = np.log(C["phi"])
N = lp.shape[0]
print(f"REPRODUCTION: total graft edge rows = {N}   (x3 reported 5162)")
print(f"  pooled Var(log phi) = {np.var(lp):.4f}   (x3: 0.718)")
print(f"  median phi          = {np.median(C['phi']):.4f}  (x3: 1.053 on the K subset)")
print()

# ══ x3's floor, verbatim ═════════════════════════════════════════════════════════════════════════
w_mu = C["s"] / np.maximum(C["phi"], _EPS)
w_nu = C["s_nu"] / np.maximum(C["phi"], _EPS)
t_nu = w_nu ** 2 * polygamma(1, np.maximum(C["k_nu"], 0.5))
t_mu = w_mu ** 2 * polygamma(1, np.maximum(C["n_spl"], 0.5))
t_R = polygamma(1, np.maximum(C["n_R"], 0.5))
floor_x3 = t_nu + t_mu + t_R
# ══ the MISSING term: the oracle capture step's own counting noise ════════════════════════════════
t_step = polygamma(1, np.maximum(C["G_bnd"], 0.5)) + polygamma(1, np.maximum(C["G_exon"], 0.5))
floor_full = floor_x3 + t_step

print("=" * 112)
print("1. THE FLOOR x3 USES HAS NO gDNA TERM, YET phi'S DENOMINATOR IS A RATIO OF TWO COUNTED gDNA POOLS")
print("=" * 112)
print(f"  {'condition':<44}{'n':>6}{'med G_bnd':>11}{'med G_exon':>12}"
      f"{'E[t_step]':>11}{'E[floor_x3]':>13}")
for ci, cond in enumerate(CONDS):
    m = C["cond_i"] == ci
    if int(m.sum()) < 5:
        continue
    print(f"  {cond[5:]:<44}{int(m.sum()):>6}{np.median(C['G_bnd'][m]):>11.1f}"
          f"{np.median(C['G_exon'][m]):>12.1f}{np.mean(t_step[m]):>11.4f}"
          f"{np.mean(floor_x3[m]):>13.4f}")
print()

_term = C["b_nterm"] > 0
_lowg = np.isin(C["cond_i"], [i for i, c in enumerate(CONDS) if ("gdna1_" in c or "gdna5_" in c)])
_g300 = np.isin(C["cond_i"], [i for i, c in enumerate(CONDS) if "gdna300" in c])


def omega(m, fl):
    v = float(np.var(lp[m]))
    f = float(np.mean(fl[m]))
    return int(m.sum()), v, f, max(0.0, v - f)


print("=" * 112)
print("2. EVERY HEADLINE omega_hat, RECOMPUTED WITH THE STEP TERM PUT BACK")
print("=" * 112)
print(f"  {'stratum':<46}{'n':>6}{'Var':>8}{'x3 floor':>10}{'x3 om':>8}"
      f"{'+step floor':>13}{'CORRECTED om':>14}")
for lab, m in (
    ("K: ALL graft edges (drop gdna1/gdna5)", ~_lowg),
    ("K: TERMINUS boundary", (~_lowg) & _term),
    ("K: JUNCTION-only boundary", (~_lowg) & (~_term)),
    ("gdna300 ONLY: ALL", _g300),
    ("gdna300 ONLY: TERMINUS", _g300 & _term),
    ("gdna300 ONLY: JUNCTION-only", _g300 & (~_term)),
    ("gdna100 ONLY: JUNCTION-only",
     np.isin(C["cond_i"], [i for i, c in enumerate(CONDS) if "gdna100" in c]) & (~_term)),
):
    n, v, f0, o0 = omega(m, floor_x3)
    _, _, f1, o1 = omega(m, floor_full)
    print(f"  {lab:<46}{n:>6}{v:>8.3f}{f0:>10.4f}{o0:>8.3f}{f1:>13.4f}{o1:>14.3f}")
print()
n1, _, _, oT = omega((~_lowg) & _term, floor_full)
n2, _, _, oJ = omega((~_lowg) & (~_term), floor_full)
print(f"  headline ratio with the step term:  {oT:.3f} / {oJ:.3f} = "
      f"{(oT / oJ if oJ > 1e-6 else float('inf')):.1f}x   (x3 claimed 1.748/0.086 = 20.3x)")
n3, _, _, oT3 = omega(_g300 & _term, floor_full)
n4, _, _, oJ3 = omega(_g300 & (~_term), floor_full)
print(f"  gdna300-only, step term in:         {oT3:.3f} / {oJ3:.3f} = "
      f"{(oT3 / oJ3 if oJ3 > 1e-6 else float('inf')):.1f}x")
print()

# ── the per-condition MEDIAN phi drift = the Jensen bias of the oracle step ──
print("=" * 112)
print("3. IS THE ORACLE FRAME BIASED AT gdna100 TOO?  E[1/step] ~ (1 + 1/G_bnd) -> phi biased UP")
print("=" * 112)
print(f"  {'condition':<44}{'n':>6}{'med phi':>9}{'med G_bnd':>11}{'pred bias 1+1/G_b':>19}")
for ci, cond in enumerate(CONDS):
    m = C["cond_i"] == ci
    if int(m.sum()) < 5:
        continue
    gb = np.median(C["G_bnd"][m])
    print(f"  {cond[5:]:<44}{int(m.sum()):>6}{np.median(C['phi'][m]):>9.3f}{gb:>11.1f}"
          f"{1.0 + 1.0 / max(gb, 1e-9):>19.3f}")
print()

# ══ 4. DECOMPOSE the "trigamma vs 1/n" claim ══════════════════════════════════════════════════════
print("=" * 112)
print("4. x3's CLAIM 'trigamma is 1.74x the 1/n floor (0.3889 vs 0.2240)' -- DECOMPOSED.")
print("   0.3889 is a THREE-term floor; 0.2240 is E[1/n_spl], ONE term. Apples to oranges.")
print("=" * 112)
for lo, hi, lab in ((0, 30, "[0,30)"), (30, 100, "[30,100)"), (1000, 1e12, "[1000,inf)")):
    m = (C["n_spl"] >= lo) & (C["n_spl"] < hi)
    if int(m.sum()) < 5:
        continue
    same_1n = (w_nu[m] ** 2 / np.maximum(C["k_nu"][m], 0.5)
               + w_mu[m] ** 2 / np.maximum(C["n_spl"][m], 0.5)
               + 1.0 / np.maximum(C["n_R"][m], 0.5))
    print(f"  n_spl {lab:<12} n={int(m.sum()):>5}")
    print(f"      E[full trigamma floor]        = {np.mean(floor_x3[m]):.4f}   (x3's column)")
    print(f"      E[same 3 terms with 1/k]      = {np.mean(same_1n):.4f}   <- the HONEST "
          f"trigamma-vs-1/n comparison: ratio {np.mean(floor_x3[m]) / max(np.mean(same_1n), 1e-12):.2f}x")
    print(f"      E[1/n_spl] alone              = "
          f"{np.mean(1.0 / np.maximum(C['n_spl'][m], 0.5)):.4f}   (x3's '1/n' column)")
    print(f"      term split: nu={np.mean(t_nu[m]):.4f}  mu={np.mean(t_mu[m]):.4f}  "
          f"R={np.mean(t_R[m]):.4f}   step(omitted)={np.mean(t_step[m]):.4f}")
print()

# ══ 5. LOCATION vs SCALE in the dispersion regression ═════════════════════════════════════════════
print("=" * 112)
print("5. THE DISPERSION REGRESSION: y = log|log phi - GLOBAL median| MIXES A GROUP MEAN SHIFT INTO")
print("   THE SPREAD. Re-centre within (n_spl quartile x n_R quartile) cells and re-run.")
print("=" * 112)


def qb(x, q):
    e = np.unique(np.quantile(x, np.linspace(0, 1, q + 1)))
    return np.clip(np.searchsorted(e, x, "right") - 1, 0, e.shape[0] - 2)


def ols(y, X, names):
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ (X.T @ y)
    r = y - X @ beta
    n, k = X.shape
    S = (X * r[:, None]).T @ (X * r[:, None])
    V = XtXi @ S @ XtXi * (n / max(n - k, 1))
    return beta, np.sqrt(np.diag(V)), names


ls = np.log(np.maximum(C["n_spl"], 0.5))
lr = np.log(np.maximum(C["n_R"], 0.5))
one = np.ones_like(lp)
bs, br = qb(C["n_spl"], 4), qb(C["n_R"], 4)
print(f"  corr(log n_spl, log n_R) = {np.corrcoef(ls, lr)[0, 1]:.3f}   "
      f"(collinearity: joint coefficients are unstable if this is near 1)")
print()
print("  cell MEDIANS of log phi (the location surface x3's global centring ignores):")
print("  " + "n_spl\\n_R".ljust(14) + "".join(f"{f'nR Q{j+1}':>12}" for j in range(4)))
for i in range(4):
    row = f"  Q{i+1}".ljust(14)
    for j in range(4):
        m = (bs == i) & (br == j)
        row += (f"{np.median(lp[m]):>12.3f}" if int(m.sum()) >= 8 else f"{'-':>12}")
    print(row)
print()
for lab, ctr in (("GLOBAL median (x3's y)", np.full(N, np.median(lp))),
                 ("cell median (n_spl x n_R)", None),
                 ("class x cell median", None)):
    if lab.startswith("cell"):
        ctr = np.zeros(N)
        for i in range(4):
            for j in range(4):
                m = (bs == i) & (br == j)
                if m.sum() > 0:
                    ctr[m] = np.median(lp[m])
    elif lab.startswith("class"):
        ctr = np.zeros(N)
        for t in (0, 1):
            for i in range(4):
                for j in range(4):
                    m = (bs == i) & (br == j) & (_term == bool(t))
                    if m.sum() >= 8:
                        ctr[m] = np.median(lp[m])
                    elif m.sum() > 0:
                        ctr[m] = np.median(lp[(bs == i) & (_term == bool(t))])
    dev = np.abs(lp - ctr)
    y = np.log(np.maximum(dev, 1e-6))
    k = dev > 1e-9
    b, se, nm = ols(y[k], np.column_stack([one[k], ls[k], lr[k]]),
                    ("const", "log n_spl", "log n_R"))
    print(f"  centring = {lab:<30} n={int(k.sum())}  "
          + "   ".join(f"{a}={x:+.3f}({s:.3f})" for a, x, s in zip(nm, b, se)))
print()
print("  and the same, WITHIN the junction-only stratum only (x3's I2 'class C/D' claim):")
for lab, m0 in (("JUNCTION-only", ~_term), ("TERMINUS", _term)):
    ctr = np.zeros(N)
    for i in range(4):
        mm = (bs == i) & m0
        if mm.sum() > 0:
            ctr[mm] = np.median(lp[mm])
    dev = np.abs(lp - ctr)
    y = np.log(np.maximum(dev, 1e-6))
    k = (dev > 1e-9) & m0
    b, se, nm = ols(y[k], np.column_stack([one[k], ls[k], lr[k]]),
                    ("const", "log n_spl", "log n_R"))
    print(f"    {lab:<16} n={int(k.sum())}  "
          + "   ".join(f"{a}={x:+.3f}({s:.3f})" for a, x, s in zip(nm, b, se)))
print()

# ══ 6. PSEUDO-REPLICATION: collapse to one row per DISTINCT edge, then redo the headline ══════════
print("=" * 112)
print("6. COLLAPSE THE PSEUDO-REPLICATION: one row per DISTINCT edge (median log phi over conditions)")
print("=" * 112)
key = (C["exon_node"] * 1e7 + C["bnd_node"] * 10 + C["face"]).astype(np.int64)
sel = ~_lowg
uk, inv = np.unique(key[sel], return_inverse=True)
lpc = np.array([np.median(lp[sel][inv == k]) for k in range(uk.shape[0])])
flc = np.array([np.mean(floor_full[sel][inv == k]) for k in range(uk.shape[0])])
fl3 = np.array([np.mean(floor_x3[sel][inv == k]) for k in range(uk.shape[0])])
tmc = np.array([C["b_nterm"][sel][inv == k][0] > 0 for k in range(uk.shape[0])])
nc = np.array([np.median(C["n_spl"][sel][inv == k]) for k in range(uk.shape[0])])
print(f"  distinct edges: {uk.shape[0]}  (x3 reported 612)")
print(f"  {'stratum':<34}{'n':>6}{'Var':>9}{'x3 floor':>11}{'x3 om':>9}"
      f"{'+step floor':>13}{'CORRECTED om':>14}")
for lab, m in (("ALL", np.ones_like(tmc, bool)), ("TERMINUS", tmc), ("JUNCTION-only", ~tmc)):
    # a median of ~4-6 conditions shrinks the per-row count noise by ~1/n_cond -- but the conditions
    # are 0.80-0.97 correlated, so the shrink is far less than 1/n_cond; report the raw floor as an
    # UPPER bound and the un-shrunk one as the LOWER bound on omega.
    v = float(np.var(lpc[m]))
    print(f"  {lab:<34}{int(m.sum()):>6}{v:>9.3f}{np.mean(fl3[m]):>11.4f}"
          f"{max(0.0, v - np.mean(fl3[m])):>9.3f}{np.mean(flc[m]):>13.4f}"
          f"{max(0.0, v - np.mean(flc[m])):>14.3f}")
print()
print("  count-flatness inside JUNCTION-only, on DISTINCT edges (excess vs the +step floor):")
for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
    m = (~tmc) & (nc >= lo) & (nc < hi)
    if int(m.sum()) < 10:
        continue
    print(f"    n_spl [{lo},{hi if hi < 1e11 else 'inf'})".ljust(26)
          + f"n={int(m.sum()):>4}  Var={np.var(lpc[m]):>7.3f}  floor={np.mean(flc[m]):>7.4f}  "
            f"excess={np.var(lpc[m]) - np.mean(flc[m]):>7.3f}")
print()

# ══ 7. IS THE REGION-CONTAINED TRUTH MASS ACTUALLY INTEGER (caveat 8)? ═══════════════════════════
print("=" * 112)
print("7. caveat-8 check: are the counts used in the trigamma floor really COUNTS?")
print("=" * 112)
fr = np.abs(C["n_R"] - np.round(C["n_R"]))
print(f"  n_R  (exon contained RNA mass): max |x - round(x)| = {fr.max():.4f}; "
      f"frac non-integer = {(fr > 1e-6).mean():.3%}")
fg = np.abs(C["G_exon"] - np.round(C["G_exon"]))
print(f"  G_exon                        : max |x - round(x)| = {fg.max():.4f}; "
      f"frac non-integer = {(fg > 1e-6).mean():.3%}")
fb = np.abs(C["G_bnd"] - np.round(C["G_bnd"]))
print(f"  G_bnd  (boundary l+r mass)    : max |x - round(x)| = {fb.max():.4f}; "
      f"frac non-integer = {(fb > 1e-6).mean():.3%}")
fk = np.abs(C["k_nu"] - np.round(C["k_nu"]))
print(f"  k_nu = n_unspl*f_R (x3's arm)  : max |x - round(x)| = {fk.max():.4f}; "
      f"frac non-integer = {(fk > 1e-6).mean():.3%}")
print()

# ══ 8. DOES THE EXON'S SPLICED MASS DOUBLE-COUNT THE JUNCTION'S? (floor correlation) ═════════════
print("=" * 112)
print("8. INDEPENDENCE of numerator and denominator (the law-of-total-variance premise).")
print("   A spliced fragment at the junction has its BLOCKS in the exon: is it in BOTH pools?")
print("=" * 112)
print(f"  median  Rs_exon / n_spl        = {np.median(C['Rs_ex'] / np.maximum(C['n_spl'], 1)):.3f}")
print(f"  median  Rs_exon / (Ru+Rs)_exon = "
      f"{np.median(C['Rs_ex'] / np.maximum(C['n_R'], _EPS)):.3f}   "
      f"<- the exon's RNA that is SPLICED mass")
print(f"  corr(log n_spl, log Rs_exon)   = "
      f"{np.corrcoef(np.log(np.maximum(C['n_spl'], .5)), np.log(np.maximum(C['Rs_ex'], .5)))[0, 1]:.3f}")
print("  If the same fragments feed both, Cov>0 and the floor OVERSTATES -> excess understated.")
