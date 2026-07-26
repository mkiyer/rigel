"""ADVERSARIAL VERIFICATION of the 'graft-peel-asymmetry' claim set.

Checks, in order:
  V1. How many times does the combine's `_transport` actually fire per sweep?  (len(_pin), len(_uni))
      -> the T5 "50% reach" number.
  V2. THE FRAME OF phi's DENOMINATOR. The solver's rho_R at a REGION node is defined by
      rho_R * E_r = the RNA share of the node's UNSPLICED contained mass (M = substrate.contained
      .mass_unspliced; the combine forms mo_p = log(cp*E_r/M)). The report's denominator is
      (Ru + Rs)/E_r -- it ADDS the exon's SPLICED contained mass. Measure both.
  V3. Spliced-vs-unspliced contained mass at exons (how big is the V2 error?).
  V4. T6's over-confidence uses the FUSED c_tau (both messages) against a PER-MESSAGE premise
      variance. Recompute with the per-message ttau.
  V5. w_prec > 1 / w_dens > 1 sanity (strand gates zero tpp but not spl_prec).
  V6. Does the suite have ALTERNATIVE SPLICING at all?

    OMP_NUM_THREADS=1 python scratchpad/v1_verify_graft.py
"""

from __future__ import annotations

import dataclasses
import importlib
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
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

# ── V6: alternative splicing in the annotation ────────────────────────────────────────────────
td = index.t_df
print("=" * 100)
print("V6. ALTERNATIVE SPLICING IN THE SUITE ANNOTATION")
print(f"   t_df cols: {list(td.columns)[:14]}")
real = td[~td["is_nrna"]] if "is_nrna" in td.columns else td
gcol = "gene_id" if "gene_id" in td.columns else td.columns[1]
per_gene = real.groupby(gcol).size()
print(f"   genes={per_gene.size}  transcripts(mRNA rows)={len(real)}  "
      f"mean tx/gene={per_gene.mean():.2f}  max={per_gene.max()}  "
      f"frac genes with >1 tx = {float((per_gene > 1).mean()):.3f}")


def qbin(x, k=5):
    q = np.quantile(x, np.linspace(0.0, 1.0, k + 1))
    q[0] -= 1e-12
    q[-1] += 1e-12
    return np.clip(np.digitize(x, q[1:-1]), 0, k - 1), q


for cond in CONDS:
    print("=" * 100)
    print(f"CONDITION {cond}")
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
    n = kind.shape[0]

    # ── V1 ────────────────────────────────────────────────────────────────────────────────────
    print(f"  V1. len(_pin) = {len(cap['_pin'])}   len(_uni) = {len(cap['_uni'])}   "
          f"len(_lvl) = {len(cap['_lvl'])}")

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    NSP, NSN = (us["spl_n_pos_l"], us["spl_n_pos_r"]), (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    pin_all = cap["_pin"]
    pin = pin_all[-2:]
    lvl = cap["_lvl"][-4:]
    uni = cap["_uni"][-1]

    # ── V3: does M (the solver's mass) equal Ru+G, i.e. is the spliced mass EXCLUDED? ──────────
    ex = np.flatnonzero(is_exon)
    resid = (M[ex] - (G[ex] + Ru[ex]))
    print(f"  V3. at EXON nodes: max|M - (G+Ru)| = {np.abs(resid).max():.3e}   "
          f"max|M - (G+Ru+Rs)| = {np.abs(M[ex] - (G[ex]+Ru[ex]+Rs[ex])).max():.3e}")
    okr = Ru[ex] > _EPS
    print(f"      median Rs/Ru at exons = {np.median((Rs[ex]/np.maximum(Ru[ex],_EPS))[okr]):.4f}"
          f"   p75={np.percentile((Rs[ex]/np.maximum(Ru[ex],_EPS))[okr],75):.4f}"
          f"   frac exons with Rs>0 = {float(np.mean(Rs[ex] > 0)):.3f}")

    # ── V2: phi with BOTH denominators, on the SAME edge set ──────────────────────────────────
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R_rep = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)  # the REPORT's
    rho_R_uns = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)  # the SOLVER's currency
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rows = []
    for df, nbr in ((0, li), (1, ri)):
        face = 1 - df
        pe = pin[df]
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SP[face][b] + SN[face][b]) / max(ESP[face][b], _EPS)
            nspl = NSP[face][b] + NSN[face][b]
            smass = SP[face][b] + SN[face][b]
            if not (mu > _EPS) or not np.isfinite(rho_R_rep[i]) or rho_R_rep[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS) or not (rho_nu_b[b] >= 0):
                continue
            if not (rho_R_uns[i] > _EPS):
                continue
            step = rho_g[b] / rho_g[i]
            num = rho_nu_b[b] + mu
            rows.append((num / (rho_R_rep[i] * step), num / (rho_R_uns[i] * step), nspl, smass, i))
    A = np.asarray(rows, float)
    ok = (A[:, 0] > 0) & (A[:, 1] > 0) & np.isfinite(A[:, 0]) & np.isfinite(A[:, 1])
    A = A[ok]
    lp_rep, lp_uns = np.log(A[:, 0]), np.log(A[:, 1])
    nspl, smass = A[:, 2], A[:, 3]
    print(f"  V2. n={A.shape[0]}   Var(log phi) REPORT[(Ru+Rs)/E_r] = {np.var(lp_rep):.4f}"
          f"   |   CORRECT[Ru/E_r] = {np.var(lp_uns):.4f}"
          f"   med log phi {np.median(lp_rep):+.3f} vs {np.median(lp_uns):+.3f}")
    vsh = 1.0 / np.maximum(smass, _EPS)
    print(f"      over-confidence (Var/median v_charged): REPORT {np.var(lp_rep)/np.median(vsh):.1f}x"
          f"   CORRECT {np.var(lp_uns)/np.median(vsh):.1f}x")
    bi, q = qbin(nspl)
    print(f"      {'n_spl bin':<22}{'edges':>7}{'Var REPORT':>12}{'Var CORRECT':>13}"
          f"{'med v':>10}{'oc REP':>9}{'oc COR':>9}")
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        vr, vc, mv = np.var(lp_rep[m]), np.var(lp_uns[m]), float(np.median(vsh[m]))
        print(f"      [{q[k]:>8.1f},{q[k+1]:>8.1f}]{'':<3}{int(m.sum()):>7d}{vr:>12.3f}{vc:>13.3f}"
              f"{mv:>10.5f}{vr/mv:>8.0f}x{vc/mv:>8.0f}x")

    # ── V4: T6 per-message ttau vs fused c_tau ────────────────────────────────────────────────
    ct = np.asarray(uni["c_tau"], float)
    atau, btau = np.asarray(uni.get("atau", np.zeros(n)), float), None
    exg = np.zeros(n, bool)
    for pe in pin:
        exg |= pe["graft"]
    live = exg & (ct > _EPS)
    print(f"  V4. fused c_tau at graft-fed exons: med 1/c_tau = "
          f"{np.median(1.0/np.maximum(ct[live],_EPS)):.4f}  n={int(live.sum())}"
          f"   (frac of graft-fed exons with c_tau==0, i.e. NO lambda claim at all: "
          f"{float(np.mean(ct[exg] <= _EPS)):.3f})")

    # ── V5: w_prec / w_dens out-of-range ──────────────────────────────────────────────────────
    wp_bad = wd_bad = tot = 0
    for df in (0, 1):
        pe = pin[df]
        gr = pe["graft"]
        pp = pe["tpp"] + pe["tpn"]
        m = gr & (pp > 0)
        wpv = pe["spl_prec"][m] / np.maximum(pp[m], _EPS)
        wp_bad += int((wpv > 1.0 + 1e-9).sum())
        tot += int(m.sum())
    print(f"  V5. w_prec > 1 on {wp_bad}/{tot} graft edges "
          f"(strand gate zeroes tpp but not spl_prec)")
    del wd_bad, atau, btau, polygamma
