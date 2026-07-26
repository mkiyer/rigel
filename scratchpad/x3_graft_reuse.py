"""Two closing checks for the P1d gap report.

(1) REUSE: is `_peel_share` already inert-but-correct at a GRAFT destination? It is called unconditionally
    at `slice(None)` in the vectorized combine, so if `SP/SN == 0` on REGION nodes then `rho_mu = 0` there,
    `peel_continue_share` hits its structural `w = 1` limit, and the whole machinery (including the
    destination's fused RNA LEVEL `nu`, the denominator of phi) is already computed for free at every exon.

(2) SHAPE: the per-edge estimator vs a pooled omega vs max(pooled, per-edge), z2 by quintile.

    OMP_NUM_THREADS=1 python scratchpad/x3_graft_reuse.py
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
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

for cond in CONDS:
    print("=" * 100)
    print(cond)
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
    isr = kind == REGION
    is_bnd, is_exon, li, ri = us["is_bnd"], us["is_exon"], us["left"], us["right"]
    n = kind.shape[0]
    SP = (us["SP_l"], us["SP_r"])
    SN = (us["SN_l"], us["SN_r"])
    NSP = (us["spl_n_pos_l"], us["spl_n_pos_r"])
    NSN = (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))

    # ── (1) REUSE ────────────────────────────────────────────────────────────────────────────────
    smax = float(max(np.abs(SP[k][isr]).max() for k in (0, 1)))
    nmax = float(max(np.abs(SN[k][isr]).max() for k in (0, 1)))
    lvl = cap["_lvl"][-4:]
    wr = [float(np.median(e["w"][isr])) for e in lvl]
    we = [float(np.median(e["w"][is_exon])) for e in lvl]
    print(f"  max |SP| over REGION nodes = {smax:.3e}   max |SN| = {nmax:.3e}")
    print(f"  _peel_share's w at REGION nodes (median, 4 face/strand slots) = {wr}")
    print(f"  _peel_share's w at EXON   nodes (median)                      = {we}")
    nu_ok = [float(np.mean(np.isfinite(e['v_nu'][is_exon]) & (e['nu'][is_exon] > 0))) for e in lvl]
    print(f"  fraction of EXON nodes where _peel_share already yields a LIVE fused RNA level nu: {nu_ok}")

    # ── (2) SHAPE ────────────────────────────────────────────────────────────────────────────────
    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        idx = np.asarray(chain.ref_idx, np.int64)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    E_g, E_r = us["E_g"], us["E_r"]
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)
    L = []
    for df, nbr in ((0, li), (1, ri)):
        face, oth = 1 - df, (ri if df == 0 else li)
        for i in np.flatnonzero(is_exon):
            b, ob = nbr[i], oth[i]
            if b < 0 or ob < 0 or not is_bnd[b] or not is_bnd[ob]:
                continue
            fo = df
            mu = (SP[face][b] + SN[face][b]) / max(ESP[face][b], _EPS)
            mo = (SP[fo][ob] + SN[fo][ob]) / max(ESP[fo][ob], _EPS)
            nt = NSP[face][b] + NSN[face][b]
            no = NSP[fo][ob] + NSN[fo][ob]
            if not (mu > _EPS and mo > _EPS and nt > 0 and no > 0):
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS and rho_R[i] > _EPS):
                continue
            ph = (rho_nu[b] + mu) / (rho_R[i] * rho_g[b] / rho_g[i])
            if np.isfinite(ph) and ph > 0:
                L.append((np.log(ph), np.log(mu / mo), nt, no))
    A = np.asarray(L, float)
    lp, d, nt, no = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
    pois = polygamma(1, np.maximum(nt, 1.0)) + polygamma(1, np.maximum(no, 1.0))
    v_edge = np.maximum(d * d - pois, 0.0) / 2.0
    om = max(0.0, float(np.var(d) - pois.mean())) / 2.0
    e2 = lp * lp
    print(f"\n  n={A.shape[0]}  E[(log phi)^2]={e2.mean():.4f}  pooled omega (MoM, trigamma)={om:.4f}")
    q = np.quantile(v_edge, np.linspace(0, 1, 6))
    q[0] -= 1e-12
    q[-1] += 1e-12
    bi = np.clip(np.digitize(v_edge, q[1:-1]), 0, 4)
    print(f"  {'v_edge quintile':<24}{'n':>6}{'E[dsq]':>9}"
          f"{'z2 per-edge':>13}{'z2 pooled':>11}{'z2 max()':>10}{'z2 sum':>9}")
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        a, bb = e2[m].mean(), v_edge[m].mean()
        mx = np.maximum(v_edge[m], om).mean()
        sm = (v_edge[m] + om).mean()
        print(f"  [{q[k]:>9.4f},{q[k+1]:>9.4f}]{'':<2}{int(m.sum()):>6d}{a:>9.3f}"
              f"{a/max(bb,1e-9):>13.2f}{a/om:>11.2f}{a/mx:>10.2f}{a/sm:>9.2f}")
    print(f"  {'ALL':<24}{A.shape[0]:>6d}{e2.mean():>9.3f}"
          f"{e2.mean()/max(v_edge.mean(),1e-9):>13.2f}{e2.mean()/om:>11.2f}"
          f"{e2.mean()/np.maximum(v_edge, om).mean():>10.2f}"
          f"{e2.mean()/(v_edge+om).mean():>9.2f}")
