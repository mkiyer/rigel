"""P1d SHAPE — is the graft's premise error driven by the junction's SHARE of the exon's RNA, or merely by
EXPRESSION (counting noise)?  Measurement only; nothing here is landed.

Starts from scratchpad/x1_graft.py (the phi measurement) and adds, per GRAFT EDGE:

    phi   = [rho_nu(bnd) + rho_mu(bnd)] / (rho_R(exon) * oracle_step)      the graft share   (x1)
    s     =  rho_mu(bnd)               / (rho_R(exon) * oracle_step)      the junction's SHARE of exon RNA
    s_nu  =  rho_nu(bnd)               / (rho_R(exon) * oracle_step)      the continue arm   (phi = s_nu + s)
    n_R   = the exon's TRUE RNA count (contained unspliced RNA + contained spliced)  -- the EXPRESSION proxy
    n_spl = the junction's integer spliced FLUX on the grafting face   (the count the graft is charged 1/n on)
    n_spl_other / mu_other = the OTHER flank's spliced flux / density  (the prior-free stand-in)
    exon length, boundary eff-spl length, region/junction annotation from the GTF.

The GTF annotation is the mechanism test: does this suite even HAVE minor-isoform junctions, and if so do
they explain phi != 1?

    OMP_NUM_THREADS=1 python scratchpad/x3_share.py
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
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
GTF = SUITE / "reference" / "genes.gtf"
_EPS = 1e-9

CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna5_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_present_capture_off",
    # nascent-FREE twins: the direct test of the "unspliced nascent RNA" mechanism
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_off",
]
N_PRESENT = 8  # the first N_PRESENT entries of CONDS are nrna_present

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
rdf = index.region_df
R_START = np.asarray(rdf["start"], np.int64)
R_END = np.asarray(rdf["end"], np.int64)
R_LEN = R_END - R_START
N_REG = R_START.shape[0]


# ── GTF structure: does the suite have alternative splicing at all? ──────────────────────────────────
def load_gtf():
    tx = collections.defaultdict(list)
    for line in open(GTF):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[tid].append((int(f[3]) - 1, int(f[4])))  # 0-based half-open
    return {t: sorted(v) for t, v in tx.items()}


TX = load_gtf()
JUNC_POS = collections.Counter()  # genomic coordinate -> # transcripts with a junction edge there
TSS_TES = collections.Counter()  # genomic coordinate -> # transcript 5'/3' termini there
EXON_IV = []
for t, ex in TX.items():
    TSS_TES[ex[0][0]] += 1
    TSS_TES[ex[-1][1]] += 1
    for a, b in zip(ex, ex[1:]):
        JUNC_POS[a[1]] += 1
        JUNC_POS[b[0]] += 1
    for e in ex:
        EXON_IV.append((e[0], e[1], t))
EXON_IV.sort()
_ex_s = np.array([e[0] for e in EXON_IV])
_ex_e = np.array([e[1] for e in EXON_IV])
_ex_t = [e[2] for e in EXON_IV]

reg_n_tx = np.zeros(N_REG, np.int64)  # # transcripts with an exon overlapping the region
reg_n_term = np.zeros(N_REG, np.int64)  # # transcript termini STRICTLY inside the region
for r in range(N_REG):
    s, e = R_START[r], R_END[r]
    hit = (_ex_s < e) & (_ex_e > s)
    reg_n_tx[r] = len({_ex_t[i] for i in np.flatnonzero(hit)})
for p, c in TSS_TES.items():
    j = int(np.searchsorted(R_START, p, "right")) - 1
    if 0 <= j < N_REG and R_START[j] < p < R_END[j]:
        reg_n_term[j] += c


# per-boundary annotations. boundary bi sits at the START of region bi (end of region bi-1); the last
# boundary (bi == N_REG) sits at the end of the last region.
def bnd_pos(bi):
    return int(R_START[bi]) if bi < N_REG else int(R_END[N_REG - 1])


B_POS = np.array([bnd_pos(b) for b in range(N_REG + 1)], np.int64)
B_NJUNC = np.array([JUNC_POS.get(int(p), 0) for p in B_POS], np.int64)
B_NTERM = np.array([TSS_TES.get(int(p), 0) for p in B_POS], np.int64)


def gtf_summary():
    g2t = collections.defaultdict(list)
    for t in TX:
        g2t[t.split(".")[0]].append(t)
    iso = collections.Counter(len(v) for v in g2t.values())
    njx = collections.Counter(JUNC_POS.values())
    print("=" * 122)
    print("A. DOES THE SIMULATOR MAKE ALTERNATIVE SPLICING?  (reference/genes.gtf)")
    print("=" * 122)
    print(f"  genes={len(g2t)}  transcripts={len(TX)}  isoforms/gene: "
          + ", ".join(f"{k}:{v}" for k, v in sorted(iso.items())))
    print(f"  distinct junction COORDINATES = {len(JUNC_POS)}; "
          "#transcripts sharing a coordinate: "
          + ", ".join(f"{k}:{v}" for k, v in sorted(njx.items())))
    print(f"  regions={N_REG}  annotated regions with >1 overlapping transcript = "
          f"{int((reg_n_tx > 1).sum())} / {int((reg_n_tx > 0).sum())}")
    print(f"  regions containing a transcript TERMINUS strictly inside = {int((reg_n_term > 0).sum())}")
    print(f"  boundaries carrying a junction coordinate = {int((B_NJUNC > 0).sum())} / {N_REG + 1}; "
          f"carrying a transcript terminus = {int((B_NTERM > 0).sum())}")
    print()


# ── the per-edge measurement ─────────────────────────────────────────────────────────────────────────
FIELDS = [
    "cond_i", "exon_node", "bnd_node", "face", "reg_idx", "bnd_idx",
    "phi", "s", "s_nu", "n_R", "n_spl", "n_spl_other", "mu", "mu_other",
    "exon_len", "e_spl", "step", "n_ct_exon", "n_ct_bnd", "f_R_bnd",
    "nas_frac_exon", "M_exon", "E_r_exon", "reg_n_tx", "reg_n_term", "b_njunc", "b_nterm",
    "k_nu", "k_mu", "k_R", "nrna",
]


def measure(cond, ci):
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
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
    Rmat = pool("mat_uns_pos") + pool("mat_uns_neg")
    Rnas = pool("nas_uns_pos") + pool("nas_uns_neg")
    Ru = Rmat + Rnas
    Rs = pool("mat_spl") + pool("nas_spl")
    E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    NSP = (us["spl_n_pos_l"] + us["spl_n_neg_l"], us["spl_n_pos_r"] + us["spl_n_neg_r"])
    NUN = (us["n_unspl_l"], us["n_unspl_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    _node_region_type(chain, ra)  # parity with x1

    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)
    f_R_bnd = np.where((Ru + G) > _EPS, Ru / np.maximum(Ru + G, _EPS), 0.0)
    nas_frac = np.where(Ru > _EPS, Rnas / np.maximum(Ru, _EPS), 0.0)
    nrna = 1.0 if "nrna_present" in cond else 0.0

    rows = []
    for face, nbr, oface, onbr in ((1, li, 0, ri), (0, ri, 1, li)):
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
            ob = onbr[i]
            if ob >= 0 and is_bnd[ob]:
                mu_o = (SPf[oface][ob] + SNf[oface][ob]) / max(ESP[oface][ob], _EPS)
                n_o = NSP[oface][ob]
            else:
                mu_o, n_o = np.nan, np.nan
            k_mu = NSP[face][b]
            k_nu = (NUN[0][b] + NUN[1][b]) * f_R_bnd[b]
            k_R = Ru[i] + Rs[i]
            rows.append((
                ci, i, b, face, idx[i], idx[b],
                (rho_nu_b[b] + mu) / den, mu / den, rho_nu_b[b] / den,
                Ru[i] + Rs[i], k_mu, n_o, mu, mu_o,
                R_LEN[idx[i]], ESP[face][b], step, NUN[0][i], NUN[0][b] + NUN[1][b],
                f_R_bnd[b], nas_frac[i], M[i], E_r[i],
                reg_n_tx[idx[i]], reg_n_term[idx[i]], B_NJUNC[idx[b]], B_NTERM[idx[b]],
                k_nu, k_mu, k_R, nrna,
            ))
    return rows


# ── helpers ──────────────────────────────────────────────────────────────────────────────────────────
def qbins(x, q):
    e = np.unique(np.quantile(x, np.linspace(0, 1, q + 1)))
    b = np.clip(np.searchsorted(e, x, "right") - 1, 0, e.shape[0] - 2)
    return e, b


def ols(y, X, names):
    """OLS with HC1 robust SEs. X already includes the intercept column."""
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ (X.T @ y)
    r = y - X @ beta
    n, k = X.shape
    S = (X * r[:, None]).T @ (X * r[:, None])
    V = XtXi @ S @ XtXi * (n / max(n - k, 1))
    return beta, np.sqrt(np.diag(V)), names


def sp_corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 5:
        return np.nan, np.nan
    ra = np.argsort(np.argsort(a[m])).astype(float)
    rb = np.argsort(np.argsort(b[m])).astype(float)
    return float(np.corrcoef(a[m], b[m])[0, 1]), float(np.corrcoef(ra, rb)[0, 1])


# ── main ─────────────────────────────────────────────────────────────────────────────────────────────
gtf_summary()

all_rows = []
for ci, cond in enumerate(CONDS):
    all_rows += measure(cond, ci)
A = np.array(all_rows, float)
C = {f: A[:, k] for k, f in enumerate(FIELDS)}
ok = np.isfinite(C["phi"]) & (C["phi"] > 0) & (C["s"] > 0)
C = {f: C[f][ok] for f in FIELDS}
lp = np.log(C["phi"])
print(f"total graft edges (capture-OFF, {len(CONDS)} conditions): {lp.shape[0]}")
print()

# ── B. per-condition, and the pseudo-replication check ───────────────────────────────────────────────
print("=" * 122)
print("B. PER-CONDITION phi, AND WHETHER POOLING CONDITIONS IS PSEUDO-REPLICATION")
print("=" * 122)
print(f"{'condition':<46}{'n':>6}{'med phi':>9}{'Var(log phi)':>14}{'sd':>8}{'nascent?':>10}")
for ci, cond in enumerate(CONDS):
    m = C["cond_i"] == ci
    if int(m.sum()) < 5:
        continue
    print(f"{cond[5:]:<46}{int(m.sum()):>6}{np.median(C['phi'][m]):>9.3f}"
          f"{np.var(lp[m]):>14.3f}{np.std(lp[m]):>8.3f}"
          f"{'present' if C['nrna'][m][0] else 'NONE':>10}")
key = C["exon_node"] * 1e7 + C["bnd_node"] * 10 + C["face"]
for a, b in ((0, 1), (0, 2), (0, 8), (8, 9)):
    ma, mb = C["cond_i"] == a, C["cond_i"] == b
    ka, kb = key[ma], key[mb]
    _, ia, ib = np.intersect1d(ka, kb, return_indices=True)
    if ia.size > 20:
        r = np.corrcoef(lp[ma][ia], lp[mb][ib])[0, 1]
        print(f"  edge-matched corr(log phi)  {CONDS[a][5:]:<40} vs {CONDS[b][5:]:<40} "
              f"n={ia.size:>4}  r={r:.3f}")
print()

# ── C. THE MECHANISM ─────────────────────────────────────────────────────────────────────────────────
print("=" * 122)
print("C. THE MECHANISM — what actually breaks the graft premise here")
print("=" * 122)
sub = C["cond_i"] < N_PRESENT
sub0 = ~sub
for lab, m in (("nrna PRESENT", sub), ("nrna NONE   ", sub0)):
    print(f"  {lab}: n={int(m.sum()):>5}  med phi={np.median(C['phi'][m]):.3f}  "
          f"Var(log phi)={np.var(lp[m]):.3f}   med s={np.median(C['s'][m]):.3f}  "
          f"med s_nu={np.median(C['s_nu'][m]):.3f}")
print()
print("  C1. by the exon region's ISOFORM MULTIPLICITY (# transcripts with an exon over the region)")
print(f"  {'n_tx over exon region':<28}{'n':>7}{'med phi':>9}{'Var(log phi)':>14}{'med s':>9}{'med s_nu':>10}")
for lo, hi, lab in ((1, 1, "1  (single isoform)"), (2, 2, "2"), (3, 99, ">=3")):
    m = (C["reg_n_tx"] >= lo) & (C["reg_n_tx"] <= hi)
    if int(m.sum()) < 5:
        continue
    print(f"  {lab:<28}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{np.var(lp[m]):>14.3f}"
          f"{np.median(C['s'][m]):>9.3f}{np.median(C['s_nu'][m]):>10.3f}")
print()
print("  C2. by transcript TERMINI and junction multiplicity")
print(f"  {'subset':<42}{'n':>7}{'med phi':>9}{'Var(log phi)':>14}{'med s':>9}")
for lab, m in (
    ("terminus INSIDE exon region", C["reg_n_term"] > 0),
    ("no terminus inside", C["reg_n_term"] == 0),
    ("boundary IS a junction coord (>=1 tx)", C["b_njunc"] > 0),
    ("boundary junction shared by >=2 tx", C["b_njunc"] > 1),
    ("boundary NOT a junction coord", C["b_njunc"] == 0),
    ("boundary carries a tx terminus", C["b_nterm"] > 0),
):
    if int(m.sum()) < 5:
        continue
    print(f"  {lab:<42}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{np.var(lp[m]):>14.3f}"
          f"{np.median(C['s'][m]):>9.3f}")
print()
print("  C3. THE CLEAN SUBSET — single-isoform exon region, no terminus inside, boundary a genuine")
print("      junction used by exactly ONE transcript. If the premise is exact anywhere it is here.")
clean = (C["reg_n_tx"] == 1) & (C["reg_n_term"] == 0) & (C["b_njunc"] == 1) & (C["b_nterm"] == 0)
for lab, m in (("clean & nrna PRESENT", clean & sub), ("clean & nrna NONE", clean & sub0),
               ("dirty & nrna NONE", (~clean) & sub0)):
    if int(m.sum()) < 5:
        print(f"      {lab:<24} n={int(m.sum())}  (too few)")
        continue
    print(f"      {lab:<24} n={int(m.sum()):>5}  med phi={np.median(C['phi'][m]):.3f}  "
          f"Var(log phi)={np.var(lp[m]):.3f}  med s={np.median(C['s'][m]):.3f}  "
          f"med s_nu={np.median(C['s_nu'][m]):.3f}")
print()
print("  C4. by EXON-REGION LENGTH vs the fragment length (frag_mean=200, frag_std=0). The contained")
print("      eff-length E_r = E[max(0, L-l)] collapses on short regions — a GEOMETRIC mechanism.")
e_len, bl = qbins(C["exon_len"], 5)
print(f"  {'exon_len bin':<24}{'n':>7}{'med E_r':>10}{'med phi':>9}{'Var(log phi)':>14}{'med n_R':>10}"
      f"{'med s':>9}")
for k in range(len(e_len) - 1):
    m = bl == k
    if int(m.sum()) < 5:
        continue
    print(f"  [{e_len[k]:>7.0f},{e_len[k + 1]:>7.0f})".ljust(24)
          + f"{int(m.sum()):>7}{np.median(C['E_r_exon'][m]):>10.1f}{np.median(C['phi'][m]):>9.3f}"
          f"{np.var(lp[m]):>14.3f}{np.median(C['n_R'][m]):>10.1f}{np.median(C['s'][m]):>9.3f}")
print()

# ── 1. THE 2-D CROSS-TAB ─────────────────────────────────────────────────────────────────────────────
print("=" * 122)
print("1. 2-D CROSS-TAB  Var(log phi)/n   rows = junction spliced count n_spl (quartiles), "
      "cols = exon RNA count n_R (quartiles)")
print("=" * 122)
qs, bs = qbins(C["n_spl"], 4)
qr, br = qbins(C["n_R"], 4)
print("   n_spl edges:", np.round(qs, 1), "   n_R edges:", np.round(qr, 1))
print("  " + "n_spl \\ n_R".ljust(24)
      + "".join(f"{f'nR Q{j + 1}':>16}" for j in range(len(qr) - 1)) + f"{'ROW':>16}")
for i in range(len(qs) - 1):
    row = f"  Q{i + 1} [{qs[i]:.0f},{qs[i + 1]:.0f})".ljust(24)
    for j in range(len(qr) - 1):
        m = (bs == i) & (br == j)
        row += (f"{np.var(lp[m]):>9.3f}/{int(m.sum()):<6d}" if int(m.sum()) >= 8 else f"{'-':>16}")
    m = bs == i
    row += f"{np.var(lp[m]):>9.3f}/{int(m.sum()):<6d}"
    print(row)
row = "  COL".ljust(24)
for j in range(len(qr) - 1):
    m = br == j
    row += f"{np.var(lp[m]):>9.3f}/{int(m.sum()):<6d}"
print(row + f"{np.var(lp):>9.3f}/{lp.shape[0]:<6d}")
print()
print("   the same cross-tab with the TRUE SHARE s on the rows instead of the count:")
qsh, bsh = qbins(C["s"], 4)
print("   s edges:", np.round(qsh, 3))
print("  " + "s \\ n_R".ljust(24)
      + "".join(f"{f'nR Q{j + 1}':>16}" for j in range(len(qr) - 1)) + f"{'ROW':>16}")
for i in range(len(qsh) - 1):
    row = f"  Q{i + 1} [{qsh[i]:.3f},{qsh[i + 1]:.3f})".ljust(24)
    for j in range(len(qr) - 1):
        m = (bsh == i) & (br == j)
        row += (f"{np.var(lp[m]):>9.3f}/{int(m.sum()):<6d}" if int(m.sum()) >= 8 else f"{'-':>16}")
    m = bsh == i
    row += f"{np.var(lp[m]):>9.3f}/{int(m.sum()):<6d}"
    print(row)
print()

# ── 2. Var(log phi) by the SHARE s, with the exact trigamma Poisson floor ────────────────────────────
print("=" * 122)
print("2. Var(log phi) BY THE TRUE SHARE s, quintiles.  Poisson floor uses the EXACT trigamma psi'(k):")
print("   floor = w_nu^2 psi'(k_nu) + w_mu^2 psi'(k_mu) + psi'(k_R),  w = the arm's share of phi")
print("=" * 122)
w_mu = C["s"] / np.maximum(C["phi"], _EPS)
w_nu = C["s_nu"] / np.maximum(C["phi"], _EPS)
floor = (w_nu ** 2 * polygamma(1, np.maximum(C["k_nu"], 0.5))
         + w_mu ** 2 * polygamma(1, np.maximum(C["k_mu"], 0.5))
         + polygamma(1, np.maximum(C["k_R"], 0.5)))
q5, b5 = qbins(C["s"], 5)
print(f"  {'s bin':<22}{'n':>6}{'med s':>8}{'med phi':>9}{'Var(log phi)':>14}{'E[floor]':>10}"
      f"{'excess':>10}{'med n_spl':>11}{'med n_R':>9}")
for k in range(len(q5) - 1):
    m = b5 == k
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  [{q5[k]:.3f},{q5[k + 1]:.3f})".ljust(22)
          + f"{int(m.sum()):>6}{np.median(C['s'][m]):>8.3f}{np.median(C['phi'][m]):>9.3f}"
          f"{v:>14.3f}{fl:>10.4f}{v - fl:>10.3f}{np.median(C['n_spl'][m]):>11.0f}"
          f"{np.median(C['n_R'][m]):>9.0f}")
print()
print("  the ledger's own axis (junction spliced count bins), recomputed with the exact trigamma floor:")
print(f"  {'n_spl bin':<22}{'n':>6}{'Var(log phi)':>14}{'E[trigamma floor]':>20}{'1/n_spl':>10}"
      f"{'excess':>10}{'med s':>9}")
for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
    m = (C["n_spl"] >= lo) & (C["n_spl"] < hi)
    if int(m.sum()) < 5:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  [{lo},{hi if hi < 1e11 else 'inf'})".ljust(22)
          + f"{int(m.sum()):>6}{v:>14.3f}{fl:>20.4f}"
          f"{np.mean(1.0 / np.maximum(C['n_spl'][m], 0.5)):>10.4f}{v - fl:>10.3f}"
          f"{np.median(C['s'][m]):>9.3f}")
print()

# ── 3. the regression ────────────────────────────────────────────────────────────────────────────────
print("=" * 122)
print("3. REGRESSION.  y = log|log phi - med(log phi)|.  For any scale family E[log|X|] = log sigma + c,")
print("   so a coefficient b is d log sigma / d log x  and  2b is d log Var / d log x.  HC1 robust SEs.")
print("=" * 122)
dev = np.abs(lp - np.median(lp))
y = np.log(np.maximum(dev, 1e-6))
keep = dev > 1e-6
ls = np.log(np.maximum(C["n_spl"], 0.5))
lr = np.log(np.maximum(C["n_R"], 0.5))
lsh = np.log(C["s"])
one = np.ones_like(y)
for lab, cols, nm in (
    ("log n_spl alone", (one, ls), ("const", "log n_spl")),
    ("log n_R alone", (one, lr), ("const", "log n_R")),
    ("BOTH jointly", (one, ls, lr), ("const", "log n_spl", "log n_R")),
    ("log s alone", (one, lsh), ("const", "log s")),
    ("log s + log n_R", (one, lsh, lr), ("const", "log s", "log n_R")),
    ("log s + log n_spl + log n_R", (one, lsh, ls, lr), ("const", "log s", "log n_spl", "log n_R")),
):
    X = np.column_stack([c[keep] for c in cols])
    beta, se, nms = ols(y[keep], X, nm)
    print(f"  {lab:<30}n={int(keep.sum())}  "
          + "   ".join(f"{n}={b:+.3f}({s:.3f})" for n, b, s in zip(nms, beta, se)))
print()

# ── 4. IS THE SHARE ESTIMABLE PRIOR-FREE? ────────────────────────────────────────────────────────────
print("=" * 122)
print("4. THE PRIOR-FREE QUESTION — can s be estimated at the graft site without the oracle?")
print("=" * 122)
cand = {
    "(a)  mu/(mu+mu_other)   two-junction split": C["mu"] / np.maximum(C["mu"] + C["mu_other"], _EPS),
    "(a') n_spl/(n_spl+n_spl_other)": C["n_spl"] / np.maximum(C["n_spl"] + C["n_spl_other"], _EPS),
    "(b)  mu/(M_exon/E_r_exon)  obs-mass bound":
        C["mu"] / np.maximum(C["M_exon"] / np.maximum(C["E_r_exon"], _EPS), _EPS),
    "(c)  n_spl alone": C["n_spl"],
    "(c') log n_spl": ls,
    "(d)  n_spl/(n_spl + n_unspl_contained_exon)": C["n_spl"] / np.maximum(C["n_spl"] + C["n_ct_exon"], _EPS),
}
print(f"  {'candidate':<46}{'n':>6}{'Pearson r (raw)':>18}{'Pearson r (logs)':>19}{'Spearman':>11}")
for lab, v in cand.items():
    m = np.isfinite(v)
    r_raw, r_sp = sp_corr(v[m], C["s"][m])
    ml = m & (v > 0)
    r_log, _ = sp_corr(np.log(v[ml]), np.log(C["s"][ml]))
    print(f"  {lab:<46}{int(m.sum()):>6}{r_raw:>18.3f}{r_log:>19.3f}{r_sp:>11.3f}")
print()
print("  does any candidate reproduce the Var(log phi) GRADIENT the true s produces?")
for best in ("(a)  mu/(mu+mu_other)   two-junction split", "(b)  mu/(M_exon/E_r_exon)  obs-mass bound",
             "(c') log n_spl"):
    v = cand[best]
    m = np.isfinite(v)
    qc, bc = qbins(v[m], 5)
    lpm, sm = lp[m], C["s"][m]
    print(f"  quintiles of {best}   (n={int(m.sum())})")
    print(f"  {'candidate bin':<26}{'n':>6}{'med TRUE s':>12}{'Var(log phi)':>14}")
    for k in range(len(qc) - 1):
        mm = bc == k
        print(f"  [{qc[k]:.4f},{qc[k + 1]:.4f})".ljust(26)
              + f"{int(mm.sum()):>6}{np.median(sm[mm]):>12.3f}{np.var(lpm[mm]):>14.3f}")
    print()
print("  pooled method-of-moments omega (the fallback if the per-junction shape is not estimable):")
for lab, m in (("nrna PRESENT", sub), ("nrna NONE", sub0), ("ALL", np.ones_like(sub, bool))):
    print(f"    {lab:<16} n={int(m.sum()):>5}  Var={np.var(lp[m]):.3f}  "
          f"E[trigamma floor]={np.mean(floor[m]):.4f}  "
          f"omega_hat={max(0.0, np.var(lp[m]) - float(np.mean(floor[m]))):.3f}")

# ── D. THE STRUCTURAL SHARE — the hypothesis stated in ANNOTATION units, no measured density anywhere ──
# s_struct = (abundance of transcripts that SPLICE INTO this exon region through THIS boundary)
#          / (abundance of all RNA present in the exon region: every transcript with an exon over it,
#             plus, when nascent is on, every transcript whose SPAN covers it).
# This is the hypothesis's own quantity — "the junction's share of the exon's transcripts" — and it is
# built from the GTF + the simulator's TRUE per-transcript abundances only. Because it never touches phi's
# numerator or denominator, binning Var(log phi) by it is NOT tautological, unlike binning by s itself
# (phi = s_nu + s and s_nu is a ~6% arm, so phi ~ s almost identically).
import pandas as pd  # noqa: E402

_ab = pd.read_csv(SUITE / "truth_abundances_nrna_present.tsv", sep="\t")
AB = dict(zip(_ab["transcript_id"], _ab["mrna_abundance"].astype(float)))
NAB = dict(zip(_ab["transcript_id"], _ab["nrna_abundance"].astype(float)))

A_mat = np.zeros(N_REG)     # mature abundance present in the region (exonic)
A_nas = np.zeros(N_REG)     # nascent abundance present in the region (span-covering)
ACC_AB = np.zeros(N_REG + 1)  # abundance SPLICING IN from the left at this boundary (acceptor)
DON_AB = np.zeros(N_REG + 1)  # abundance SPLICING OUT to the right at this boundary (donor)
_atomic_viol = 0
for t, ex in TX.items():
    a, nn = AB.get(t, 0.0), NAB.get(t, 0.0)
    for k, (s0, e0) in enumerate(ex):
        lo = int(np.searchsorted(R_START, s0, "right")) - 1
        hi = int(np.searchsorted(R_START, e0, "left"))
        for r in range(max(lo, 0), min(hi, N_REG)):
            if R_END[r] <= s0 or R_START[r] >= e0:
                continue
            if R_START[r] < s0 or R_END[r] > e0:
                _atomic_viol += 1
            A_mat[r] += a
        if k > 0:  # acceptor at s0
            j = int(np.searchsorted(R_START, s0))
            if j <= N_REG and (j == N_REG or R_START[min(j, N_REG - 1)] == s0):
                ACC_AB[j] += a
        if k < len(ex) - 1:  # donor at e0
            j = int(np.searchsorted(R_START, e0))
            if j <= N_REG and (j == N_REG or R_START[min(j, N_REG - 1)] == e0):
                DON_AB[j] += a
    sp0, sp1 = ex[0][0], ex[-1][1]
    lo = int(np.searchsorted(R_START, sp0, "right")) - 1
    hi = int(np.searchsorted(R_START, sp1, "left"))
    for r in range(max(lo, 0), min(hi, N_REG)):
        if R_END[r] <= sp0 or R_START[r] >= sp1:
            continue
        A_nas[r] += nn

ri_ = C["reg_idx"].astype(int)
bi_ = C["bnd_idx"].astype(int)
face1 = C["face"] == 1  # boundary on the LEFT of the exon -> the exon is the ACCEPTOR side
num_ab = np.where(face1, ACC_AB[np.clip(bi_, 0, N_REG)], DON_AB[np.clip(bi_, 0, N_REG)])
den_ab = A_mat[ri_] + np.where(C["nrna"] > 0, A_nas[ri_], 0.0)
s_struct = np.where(den_ab > 0, num_ab / np.maximum(den_ab, _EPS), np.nan)
s_other = np.where(C["mu"] > _EPS, C["s"] * C["mu_other"] / np.maximum(C["mu"], _EPS), np.nan)

print()
print("=" * 122)
print("D. THE STRUCTURAL SHARE (annotation + TRUE per-transcript abundances; touches neither side of phi)")
print("=" * 122)
print(f"  region/exon atomicity violations (an exon partially covering a region): {_atomic_viol}")
print(f"  s_struct finite on {int(np.isfinite(s_struct).sum())}/{s_struct.shape[0]} edges; "
      f"exactly 1.0 (the whole exon splices in here) on {int(np.isclose(s_struct, 1.0).sum())}")
r_raw, r_sp = sp_corr(s_struct, C["s"])
print(f"  corr(s_struct, TRUE s):  Pearson {r_raw:.3f}   Spearman {r_sp:.3f}")
r_raw2, r_sp2 = sp_corr(s_struct, C["phi"])
print(f"  corr(s_struct, phi):     Pearson {r_raw2:.3f}   Spearman {r_sp2:.3f}")
print()
print("  D1. Var(log phi) BINNED BY s_struct — the NON-TAUTOLOGICAL version of table 2")
print(f"  {'s_struct bin':<24}{'n':>7}{'med s_struct':>14}{'med TRUE s':>12}{'med phi':>9}"
      f"{'Var(log phi)':>14}{'med n_spl':>11}{'med n_R':>10}")
for lab, m in (
    ("== 0  (no junction)", np.isclose(s_struct, 0.0)),
    ("(0, 0.25)", (s_struct > 0) & (s_struct < 0.25)),
    ("[0.25, 0.60)", (s_struct >= 0.25) & (s_struct < 0.60)),
    ("[0.60, 0.95)", (s_struct >= 0.60) & (s_struct < 0.95)),
    ("[0.95, 1.0]  (sole)", s_struct >= 0.95),
):
    m = m & np.isfinite(s_struct)
    if int(m.sum()) < 5:
        continue
    print(f"  {lab:<24}{int(m.sum()):>7}{np.median(s_struct[m]):>14.3f}{np.median(C['s'][m]):>12.3f}"
          f"{np.median(C['phi'][m]):>9.3f}{np.var(lp[m]):>14.3f}{np.median(C['n_spl'][m]):>11.0f}"
          f"{np.median(C['n_R'][m]):>10.0f}")
print()
print("  D2. the same, restricted to nrna NONE (no nascent arm at all, so phi is purely the mature graft)")
print(f"  {'s_struct bin':<24}{'n':>7}{'med TRUE s':>12}{'med phi':>9}{'Var(log phi)':>14}")
for lab, mm in (
    ("== 0  (no junction)", np.isclose(s_struct, 0.0)),
    ("(0, 0.60)", (s_struct > 0) & (s_struct < 0.60)),
    ("[0.60, 0.95)", (s_struct >= 0.60) & (s_struct < 0.95)),
    ("[0.95, 1.0]  (sole)", s_struct >= 0.95),
):
    m = mm & np.isfinite(s_struct) & sub0
    if int(m.sum()) < 5:
        continue
    print(f"  {lab:<24}{int(m.sum()):>7}{np.median(C['s'][m]):>12.3f}{np.median(C['phi'][m]):>9.3f}"
          f"{np.var(lp[m]):>14.3f}")
print()
print("  D3. Var(log phi) binned by the OTHER flank's true share s_other (also non-tautological: s_other")
print("      appears nowhere in phi). The hypothesis predicts: the more the OTHER junction carries, the")
print("      less THIS one represents.")
mo = np.isfinite(s_other) & (s_other > 0)
qo, bo = qbins(s_other[mo], 5)
lpo, so_s = lp[mo], C["s"][mo]
print(f"  {'s_other bin':<24}{'n':>7}{'med TRUE s(this)':>18}{'Var(log phi)':>14}")
for k in range(len(qo) - 1):
    m = bo == k
    print(f"  [{qo[k]:.3f},{qo[k + 1]:.3f})".ljust(24)
          + f"{int(m.sum()):>7}{np.median(so_s[m]):>18.3f}{np.var(lpo[m]):>14.3f}")
print()

# ── E. VARIANCE DECOMPOSITION BY STRUCTURAL CLASS ────────────────────────────────────────────────────
print("=" * 122)
print("E. WHERE THE VARIANCE LIVES — a disjoint structural partition of every graft edge")
print("=" * 122)
cls_noj = np.isclose(s_struct, 0.0) | ~np.isfinite(s_struct)
cls_term = (~cls_noj) & (C["b_nterm"] > 0)
cls_minor = (~cls_noj) & (~cls_term) & (s_struct < 0.95)
cls_sole = (~cls_noj) & (~cls_term) & (s_struct >= 0.95)
tot_v = np.var(lp)
gm = np.mean(lp)
print(f"  {'class':<52}{'n':>7}{'frac':>7}{'med phi':>9}{'Var(log phi)':>14}"
      f"{'share of TOTAL sq-dev':>24}")
for lab, m in (
    ("A: boundary carries NO junction of this exon", cls_noj),
    ("B: boundary is a transcript TERMINUS (TSS/TES)", cls_term),
    ("C: MINOR-ISOFORM junction (s_struct < 0.95)", cls_minor),
    ("D: SOLE junction into the exon (s_struct >= 0.95)", cls_sole),
):
    n = int(m.sum())
    if n < 5:
        continue
    sq = float(np.sum((lp[m] - gm) ** 2))
    print(f"  {lab:<52}{n:>7}{n / lp.shape[0]:>7.1%}{np.median(C['phi'][m]):>9.3f}"
          f"{np.var(lp[m]):>14.3f}{sq / np.sum((lp - gm) ** 2):>24.1%}")
print(f"  TOTAL Var(log phi) = {tot_v:.3f}")
print()
print("  E1. class D (sole junction) alone, broken down by n_spl — is the residual pure counting noise?")
print(f"  {'n_spl bin':<20}{'n':>7}{'med phi':>9}{'Var(log phi)':>14}{'E[trigamma floor]':>20}{'excess':>10}")
for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
    m = cls_sole & (C["n_spl"] >= lo) & (C["n_spl"] < hi)
    if int(m.sum()) < 8:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  [{lo},{hi if hi < 1e11 else 'inf'})".ljust(20)
          + f"{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{v:>14.3f}{fl:>20.4f}{v - fl:>10.3f}")
print()

# ── F. COVERAGE / SELECTION CENSUS ───────────────────────────────────────────────────────────────────
print("=" * 122)
print("F. SELECTION — which exon nodes this measurement can even see")
print("=" * 122)
_c = CONDS[0]
_inp = _scan_and_truth(SUITE, _c, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
_dbg: dict = {}
calmod.calibrate(_inp["payload"], ra, _inp["strand_model"], np.asarray(_inp["gdna_fl_pmf"]),
                 np.asarray(_inp["rna_fl_pmf"]),
                 dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=_dbg)
_us = _dbg["capture"]["_uni_static"]
_ex = _us["is_exon"]
_idx = np.asarray(_dbg["chain"].ref_idx, np.int64)
_Er = _us["E_r"]
print(f"  exon NODES in the chain: {int(_ex.sum())}")
print(f"    with a positive contained RNA eff-length E_r: {int((_ex & (_Er > 1e-9)).sum())}")
_L = R_LEN[np.clip(_idx, 0, N_REG - 1)]
print(f"    exon regions shorter than the fixed fragment length (200 bp): "
      f"{int((_ex & (_L < 200)).sum())}  -> E_r ~ 0, NO contained mass, phi UNMEASURABLE here")
print(f"    exon-region length: min={int(_L[_ex].min())} p25={int(np.percentile(_L[_ex], 25))} "
      f"med={int(np.median(_L[_ex]))} p75={int(np.percentile(_L[_ex], 75))} max={int(_L[_ex].max())}")
print("  NOTE: the shipped graft still fires on those short exons; this measurement cannot see them, so")
print("  every Var(log phi) here is conditioned on the exon being long enough to hold a fragment.")

# ── G. THE CIRCULARITY CONTROL, AND THE ONLY HONEST PRIOR-FREE TEST ──────────────────────────────────
# Every candidate whose numerator is mu inherits phi's own numerator: log(b) = log(phi) + log(f_R,exon) +
# log(step) + o(1), so binning Var(log phi) by it is partly MECHANICAL, exactly as binning by s is
# (phi = s_nu + s with s_nu a ~6% arm). The two controls that are NOT circular:
#   G1  corr(candidate, s_struct) — s_struct comes from the GTF + true abundances and touches no density.
#   G2  the candidate's Var(log phi) gradient measured WITHIN class D (sole junction, s_struct >= 0.95),
#       where the true structural share is ~1 for every edge. Any gradient surviving there is NOT share
#       information; it is the count/mechanical part.
print()
print("=" * 122)
print("G. CIRCULARITY CONTROL — does any prior-free candidate carry SHARE information, or only count?")
print("=" * 122)
_dist = np.unique(C["exon_node"] * 1e7 + C["bnd_node"] * 10 + C["face"]).shape[0]
print(f"  effective n: {lp.shape[0]} edge-rows but only {_dist} DISTINCT graft edges "
      f"(14 conditions re-measure the same edges; edge-matched corr(log phi) = 0.80-0.97),")
print("  so every n in this report is ~8.6x inflated. Treat ~600 as the independent sample size.")
print()
print(f"  G1  {'candidate':<44}{'r vs TRUE s':>13}{'r vs s_struct':>15}{'rho vs s_struct':>17}"
      f"{'rho vs n_spl':>14}")
for lab, v in cand.items():
    m = np.isfinite(v) & np.isfinite(s_struct)
    r1, _ = sp_corr(v[m], C["s"][m])
    r2, rs2 = sp_corr(v[m], s_struct[m])
    _, rn = sp_corr(v[m], C["n_spl"][m])
    print(f"      {lab:<44}{r1:>13.3f}{r2:>15.3f}{rs2:>17.3f}{rn:>14.3f}")
print(f"      {'(reference) the TRUE s itself':<44}"
      f"{1.0:>13.3f}{sp_corr(C['s'], s_struct)[0]:>15.3f}{sp_corr(C['s'], s_struct)[1]:>17.3f}"
      f"{sp_corr(C['s'], C['n_spl'])[1]:>14.3f}")
print()
print("  G2  the Var(log phi) gradient WITHIN class D (sole junction, s_struct>=0.95): no share variation")
print("      is left there, so a surviving gradient is mechanical/count, not share.")
for lab in ("(b)  mu/(M_exon/E_r_exon)  obs-mass bound", "(c') log n_spl",
            "(a)  mu/(mu+mu_other)   two-junction split"):
    v = cand[lab]
    m = cls_sole & np.isfinite(v)
    qc, bc = qbins(v[m], 4)
    lpm = lp[m]
    line = "      ".join(
        f"Q{k + 1}:{np.var(lpm[bc == k]):.3f}/{int((bc == k).sum())}" for k in range(len(qc) - 1)
    )
    vall = cand[lab]
    mall = np.isfinite(vall)
    qa, ba = qbins(vall[mall], 4)
    lpa = lp[mall]
    line2 = "      ".join(
        f"Q{k + 1}:{np.var(lpa[ba == k]):.3f}/{int((ba == k).sum())}" for k in range(len(qa) - 1)
    )
    print(f"      {lab}")
    print(f"        ALL edges     {line2}")
    print(f"        class D only  {line}")
print()
print("  G3  a structure-aware omega vs a pooled one (method of moments, trigamma floor subtracted):")
print(f"      {'stratum':<52}{'n':>7}{'Var':>9}{'floor':>9}{'omega_hat':>11}")
for lab, m in (
    ("POOLED (one omega for every graft)", np.ones_like(cls_sole, bool)),
    ("A: no junction of this exon at the boundary", cls_noj),
    ("B: boundary is a transcript terminus", cls_term),
    ("C: minor-isoform junction", cls_minor),
    ("D: sole junction", cls_sole),
):
    if int(m.sum()) < 5:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"      {lab:<52}{int(m.sum()):>7}{v:>9.3f}{fl:>9.4f}{max(0.0, v - fl):>11.3f}")
print()
print("  G4  and the STRUCTURAL classes A/B/C/D are themselves prior-free (they are ANNOTATION facts):")
print("      class membership needs only the GTF + the region partition, not any abundance. Check that")
print("      the annotation-only version of the split still separates:")
_ann_sole = (C["b_njunc"] > 0) & (C["b_nterm"] == 0) & (C["reg_n_tx"] == 1)
_ann_term = C["b_nterm"] > 0
_ann_nojx = C["b_njunc"] == 0
_ann_multi = (C["b_njunc"] > 0) & (C["b_nterm"] == 0) & (C["reg_n_tx"] > 1)
print(f"      {'annotation-only class':<52}{'n':>7}{'med phi':>9}{'Var':>9}{'floor':>9}{'omega_hat':>11}")
for lab, m in (
    ("no junction coordinate at the boundary", _ann_nojx),
    ("boundary carries a transcript terminus", _ann_term),
    ("junction, but exon region has >1 isoform", _ann_multi),
    ("junction, exon region single-isoform", _ann_sole),
):
    if int(m.sum()) < 5:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"      {lab:<52}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{v:>9.3f}{fl:>9.4f}"
          f"{max(0.0, v - fl):>11.3f}")

# ── H. THE KILL SHOT — is the ledger's 40x COUNT gradient just class composition? ─────────────────────
print()
print("=" * 122)
print("H. IS THE LEDGER'S COUNT GRADIENT A COUNT EFFECT, OR CLASS COMPOSITION?")
print("=" * 122)
_qn, _bn = qbins(C["n_spl"], 4)
print("  H1. structural-class COMPOSITION by junction spliced-count quartile")
print(f"  {'n_spl quartile':<24}{'n':>7}{'A no-junc':>11}{'B terminus':>12}{'C minor':>10}{'D sole':>9}")
for k in range(len(_qn) - 1):
    m = _bn == k
    n = int(m.sum())
    print(f"  Q{k + 1} [{_qn[k]:.0f},{_qn[k + 1]:.0f})".ljust(24)
          + f"{n:>7}{(cls_noj & m).sum() / n:>11.1%}{(cls_term & m).sum() / n:>12.1%}"
          f"{(cls_minor & m).sum() / n:>10.1%}{(cls_sole & m).sum() / n:>9.1%}")
print()
print("  H2. the EXCESS (Var - trigamma floor) by n_spl bin, WITHIN each structural class.")
print("      If the count carried premise information, the excess would fall with n INSIDE a class.")
print(f"  {'class':<22}" + "".join(f"{f'n_spl<{h}' if h < 1e11 else 'n_spl>=1000':>16}"
                                   for h in (30, 100, 300, 1000, 1e12)))
for lab, cm in (("A no-junction", cls_noj), ("B terminus", cls_term),
                ("C minor-isoform", cls_minor), ("D sole junction", cls_sole)):
    row = f"  {lab:<22}"
    for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
        m = cm & (C["n_spl"] >= lo) & (C["n_spl"] < hi)
        row += (f"{np.var(lp[m]) - float(np.mean(floor[m])):>9.3f}/{int(m.sum()):<6d}"
                if int(m.sum()) >= 15 else f"{'-':>16}")
    print(row)
print()
print("  H3. the boundary-terminus class, split by whether the coordinate is ALSO a junction")
print(f"  {'subset':<50}{'n':>7}{'med phi':>9}{'med s':>9}{'Var':>9}{'omega_hat':>11}")
for lab, m in (
    ("pure TERMINUS (no junction coord there)", (C["b_nterm"] > 0) & (C["b_njunc"] == 0)),
    ("terminus AND junction at the same coord", (C["b_nterm"] > 0) & (C["b_njunc"] > 0)),
    ("junction only, no terminus", (C["b_nterm"] == 0) & (C["b_njunc"] > 0)),
    ("neither (internal partition seam)", (C["b_nterm"] == 0) & (C["b_njunc"] == 0)),
):
    if int(m.sum()) < 5:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  {lab:<50}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{np.median(C['s'][m]):>9.3f}"
          f"{v:>9.3f}{max(0.0, v - fl):>11.3f}")
print()
print("  H4. class-A mechanism check — is a 'no-junction' boundary an INTERNAL SEAM inside an exon?")
print("      (the accumulator gives a SPLICED fragment's contiguous crossings the spliced channel too,")
print("       so ρ_μ = S/E_spl at such a seam is a CONTIGUOUS crossing misread as a splice event)")
_int_seam = np.zeros(N_REG + 1, bool)
for t, ex in TX.items():
    for s0, e0 in ex:
        lo = int(np.searchsorted(R_START, s0, "right"))
        hi = int(np.searchsorted(R_START, e0, "left"))
        for b in range(lo, hi):  # boundaries strictly inside this exon
            _int_seam[b] = True
_is_seam = _int_seam[np.clip(bi_, 0, N_REG)]
mA = cls_noj
print(f"      class-A edges = {int(mA.sum())}; of those, the boundary is strictly INSIDE an annotated "
      f"exon: {int((mA & _is_seam).sum())} ({(mA & _is_seam).sum() / max(int(mA.sum()), 1):.0%})")
print(f"      Var(log phi) on internal-seam edges = {np.var(lp[mA & _is_seam]):.3f}, "
      f"med phi = {np.median(C['phi'][mA & _is_seam]):.3f}, med s = {np.median(C['s'][mA & _is_seam]):.3f}")

# ── I. IS THE CLASS MIX AN ARTIFACT OF THIS MEASUREMENT'S FILTERS? + within-class regressions ─────────
print()
print("=" * 122)
print("I. THE SOLVER'S OWN GRAFT CENSUS (no oracle filters at all) + within-class regressions")
print("=" * 122)
_us2 = _us
_ex2, _bnd2 = _us2["is_exon"], _us2["is_bnd"]
_li2, _ri2 = _us2["left"], _us2["right"]
_SP2 = (_us2["SP_l"], _us2["SP_r"])
_SN2 = (_us2["SN_l"], _us2["SN_r"])
_idx2 = np.asarray(_dbg["chain"].ref_idx, np.int64)
tot_edges = 0
tot_mass = 0.0
cnt = collections.Counter()
mass = collections.Counter()
for face, nbr in ((1, _li2), (0, _ri2)):
    for i in np.flatnonzero(_ex2):
        b = nbr[i]
        if b < 0 or not _bnd2[b]:
            continue
        S = _SP2[face][b] + _SN2[face][b]
        if not (S > 0):
            continue
        bb = int(_idx2[b])
        rr = int(_idx2[i])
        nj, nt = B_NJUNC[min(bb, N_REG)], B_NTERM[min(bb, N_REG)]
        k = ("no-junction seam" if nj == 0 and nt == 0 else
             "terminus boundary" if nt > 0 else
             "junction, multi-isoform exon" if reg_n_tx[rr] > 1 else "junction, single-isoform exon")
        cnt[k] += 1
        mass[k] += float(S)
        tot_edges += 1
        tot_mass += float(S)
print(f"  graft edges the SOLVER actually forms (exon <- boundary with spliced mass on the exon's face), "
      f"cond {CONDS[0][5:]}: {tot_edges}")
print(f"  {'class':<36}{'edges':>8}{'% edges':>9}{'spliced MASS grafted':>23}{'% mass':>9}")
for k in ("no-junction seam", "terminus boundary", "junction, multi-isoform exon",
          "junction, single-isoform exon"):
    print(f"  {k:<36}{cnt[k]:>8}{cnt[k] / tot_edges:>9.1%}{mass[k]:>23.1f}"
          f"{mass[k] / tot_mass:>9.1%}")
print(f"  (this measurement resolves {int((C['cond_i'] == 0).sum())} of those {tot_edges} edges — the rest")
print("   lack contained exon mass (exon region shorter than the 200 bp fragment) or oracle gDNA.)")
print()
print("  I2. the deliverable-3 regression, run WITHIN each structural class")
print(f"  {'class':<24}{'n':>7}   coefficients of  y = log|log phi - median|  (HC1 SEs)")
for lab, m in (("ALL", np.ones_like(cls_sole, bool)), ("A no-junction", cls_noj),
               ("B terminus", cls_term), ("C minor-isoform", cls_minor), ("D sole junction", cls_sole)):
    mm = m & keep
    if int(mm.sum()) < 40:
        continue
    X = np.column_stack([one[mm], ls[mm], lr[mm]])
    beta, se, _ = ols(y[mm], X, ("const", "log n_spl", "log n_R"))
    print(f"  {lab:<24}{int(mm.sum()):>7}   " + "   ".join(
        f"{n}={b:+.3f}({s:.3f})" for n, b, s in zip(("const", "log n_spl", "log n_R"), beta, se)))

# ── J. THE HEADLINE — the simplest prior-free split that survives every control ───────────────────────
print()
print("=" * 122)
print("J. HEADLINE — every region boundary is a JUNCTION coordinate, a TERMINUS coordinate, or both.")
print("   That single annotation bit splits the graft's premise error 14x, with no fitting at all.")
print("=" * 122)
_term = C["b_nterm"] > 0
print(f"  {'stratum':<44}{'n':>7}{'med phi':>9}{'Var(log phi)':>14}{'floor':>9}{'omega_hat':>11}")
for lab, m in (
    ("boundary carries a TERMINUS  (all conds)", _term),
    ("   ... nrna PRESENT", _term & sub),
    ("   ... nrna NONE", _term & sub0),
    ("JUNCTION-only boundary       (all conds)", ~_term),
    ("   ... nrna PRESENT", (~_term) & sub),
    ("   ... nrna NONE", (~_term) & sub0),
):
    if int(m.sum()) < 5:
        continue
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  {lab:<44}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{v:>14.3f}{fl:>9.4f}"
          f"{max(0.0, v - fl):>11.3f}")
print()
print("  and the count-flatness inside each stratum (excess = Var - trigamma floor):")
print(f"  {'stratum':<24}" + "".join(f"{f'n_spl<{h}' if h < 1e11 else 'n_spl>=1000':>16}"
                                     for h in (30, 100, 300, 1000, 1e12)))
for lab, cm in (("TERMINUS boundary", _term), ("JUNCTION-only", ~_term)):
    row = f"  {lab:<24}"
    for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e12)):
        m = cm & (C["n_spl"] >= lo) & (C["n_spl"] < hi)
        row += (f"{np.var(lp[m]) - float(np.mean(floor[m])):>9.3f}/{int(m.sum()):<6d}"
                if int(m.sum()) >= 15 else f"{'-':>16}")
    print(row)

# ── K. ROBUSTNESS — drop the low-gDNA conditions (gdna1/gdna5), where the ORACLE CAPTURE STEP itself is
#      count-starved (rho_g is estimated from very few gDNA fragments, so `step` is noisy and inflates phi).
_lowg = np.isin(C["cond_i"], [i for i, c in enumerate(CONDS) if ("gdna1_" in c or "gdna5_" in c)])
print()
print("=" * 122)
print("K. ROBUSTNESS — the headline with the count-starved-gDNA conditions (gdna1, gdna5) DROPPED")
print("=" * 122)
print(f"  {'stratum':<44}{'n':>7}{'med phi':>9}{'Var(log phi)':>14}{'floor':>9}{'omega_hat':>11}")
for lab, m in (
    ("ALL graft edges", ~_lowg),
    ("boundary carries a TERMINUS", (~_lowg) & (C["b_nterm"] > 0)),
    ("JUNCTION-only boundary", (~_lowg) & (C["b_nterm"] == 0)),
    ("   ... nrna PRESENT", (~_lowg) & (C["b_nterm"] == 0) & sub),
    ("   ... nrna NONE", (~_lowg) & (C["b_nterm"] == 0) & sub0),
):
    v, fl = np.var(lp[m]), float(np.mean(floor[m]))
    print(f"  {lab:<44}{int(m.sum()):>7}{np.median(C['phi'][m]):>9.3f}{v:>14.3f}{fl:>9.4f}"
          f"{max(0.0, v - fl):>11.3f}")
