"""RC2 — the single-hop COUNTERFACTUAL: rebuild the combine with the graft handled differently, re-run the
exact ψ replay, and measure.

The relay is held FIXED (it is the shipped one from the capture); only the COMBINE's `_transport` is
recomputed, so this is a one-hop counterfactual — a lower bound on the effect, since the same defect also
runs inside the relay.

Variants
    shipped        exactly what ships
    graft_no_r     the boundary's MEASURED mature joins AFTER the reframe:  tp = rp[src]·r + gp
    graft_s2t      the graft still ×r, but σ²_transfer is CHARGED on it instead of being set to 0
    both

    OMP_NUM_THREADS=1 python scratchpad/rc2_c_cf.py [COND ...]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.enrichment_frame import mismatch_deflate, mismatch_gap, transfer_logvar  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.node_init import own_composition_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex import _mixture_strand_loglik  # noqa: E402
from rigel.calibration.simplex_logodds import (  # noqa: E402
    _log1m_fg, _log_fg, _logodds_grid, _lse, _posterior_median_fg, _regrid_global,
)
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9
BINS = [(0.0, 0.30), (0.30, 0.60), (0.60, 0.90), (0.90, 0.99), (0.99, 1.01)]
VARIANTS = ("shipped", "graft_no_r", "graft_s2t", "both")

ap = argparse.ArgumentParser()
ap.add_argument("conds", nargs="*", default=["gdna_gdna300_ss_0.99_nrna_none_capture_on",
                                             "gdna_gdna300_ss_0.99_nrna_present_capture_on",
                                             "gdna_gdna100_ss_0.99_nrna_present_capture_on",
                                             "gdna_gdna300_ss_0.50_nrna_present_capture_on"])
a = ap.parse_args()
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
OUT: dict[str, list] = {}


def run(cond):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap, st, geo = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    n = np.asarray(cap["f_g"]).shape[0]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    og, op, on = us["og"], us["op"], us["on"]
    is_bnd, ex_a = us["is_bnd"].astype(bool), us["is_exon"].astype(bool)
    fp_a, fn_a = us["fp"].astype(bool), us["fn"].astype(bool)
    is_amb = fp_a & fn_a
    li, ri = us["left"], us["right"]
    logvar_tot, tau_own = us["logvar_tot"], us["tau_own"]
    struct = us["struct_lock"].astype(bool)
    ni_fg = np.asarray(cap["fg_loc"])
    v_own_g, v_own_r = own_composition_logvar(ni_fg, tau_own, struct)
    v_own_lam = np.where(struct, 0.0, np.where(tau_own > _EPS, 1.0 / np.maximum(tau_own, _EPS), np.inf))
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    spf = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    snf = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    vl, vr = li >= 0, ri >= 0
    sl, sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)
    fwd = tuple(us[f"fwd_{x}"] for x in ("g", "p", "n", "pg", "pp", "pn", "mg", "mp", "mn", "tau"))
    bwd = tuple(us[f"bwd_{x}"] for x in ("g", "p", "n", "pg", "pp", "pn", "mg", "mp", "mn", "tau"))

    def transport(src, valid, df, sf, arrs, dstf, srcf, *, graft_no_r, graft_s2t):
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
        framed = valid & (srcf[src] > _EPS) & (dstf > _EPS)
        r = np.where(framed, dstf / np.maximum(srcf[src], _EPS), np.where(valid, 1.0, 0.0))
        graft = ex_a & is_bnd[src] & valid
        gp = np.where(graft, spf[sf][src], 0.0)
        gn = np.where(graft, snf[sf][src], 0.0)
        tg = rg[src] * r
        if graft_no_r:  # the MEASURED mature is already in the destination's (exonic) capture frame
            tp, tn = rp[src] * r + gp, rn[src] * r + gn
        else:
            tp, tn = (rp[src] + gp) * r, (rn[src] + gn) * r
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)
        s2t_g = (logvar_tot + logvar_tot[src]) if graft_s2t else s2t

        def _dv(p_, s2_):
            return np.where(valid & (p_[src] > 0.0), 1.0 / (1.0 / np.maximum(p_[src], _EPS) + s2_), 0.0)

        tpg, tpp, tpn = _dv(pg, s2t), _dv(pp, s2t), _dv(pn, s2t)
        tmg, tmp, tmn = _dv(mg, s2t), _dv(mp, s2t), _dv(mn, s2t)
        ttau = _dv(tau, s2t)
        _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
        _s2 = np.where(np.isfinite(s2t_g), s2t_g, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = is_bnd & ex_a[src] & valid
        tp = np.where(peel, np.maximum(tp - spf[df], 0.0), tp)
        tn = np.where(peel, np.maximum(tn - snf[df], 0.0), tn)
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        sg = np.where(tpg > 0.0, tg, og)
        sp_ = np.where(tpp > 0.0, tp, op)
        sn_ = np.where(tpn > 0.0, tn, on)
        s = sg * E_g + (sp_ + sn_) * E_r
        k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
        tg, tp, tn = tg * k, tp * k, tn * k
        g_g, c_g = mismatch_gap(tg, og)
        g_p, c_p = mismatch_gap(tp, op)
        g_n, c_n = mismatch_gap(tn, on)
        tpg, tpp, tpn = (mismatch_deflate(tpg, g_g, c_g, v_own_g),
                         mismatch_deflate(tpp, g_p, c_p, v_own_r),
                         mismatch_deflate(tpn, g_n, c_n, v_own_r))
        tmg, tmp, tmn = (mismatch_deflate(tmg, g_g, c_g, v_own_g),
                         mismatch_deflate(tmp, g_p, c_p, v_own_r),
                         mismatch_deflate(tmn, g_n, c_n, v_own_r))
        ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
        g_R, c_R = mismatch_gap(tp + tn, op + on)
        ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, v_own_lam)
        return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau

    def combine(**kw):
        A = transport(sl, vl, 0, 1, fwd, uni["rho_lf"], uni["rho_rf"], **kw)
        B = transport(sr, vr, 1, 0, bwd, uni["rho_rf"], uni["rho_lf"], **kw)
        cpg, cpp, cpn = A[3] + B[3], A[4] + B[4], A[5] + B[5]

        def fz(x, y, pa, pb, p):
            return np.where(p > _EPS, (pa * x + pb * y) / np.maximum(p, _EPS), 0.0)

        cg = fz(A[0], B[0], A[3], B[3], cpg)
        cp = fz(A[1], B[1], A[4], B[4], cpp)
        cn = fz(A[2], B[2], A[5], B[5], cpn)
        cm_g, cm_p, cm_n = A[6] + B[6], A[7] + B[7], A[8] + B[8]
        c_tau = A[9] + B[9]
        mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
        mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
        mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
        cR = cp + cn
        mo_R = np.log(np.maximum(cR * E_r / np.maximum(M, _EPS), _EPS))
        lam = mo_g - mo_R
        c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
        return dict(mo_g=mo_g, cm_g=cm_g, mo_p=mo_p, cm_p=cm_p, mo_n=mo_n, cm_n=cm_n,
                    lam=lam, c_tau=c_tau)

    # ── the ψ replay (single-strand 1-D, exact) ──
    K = int(cc.sweep_n_grid_single_strand or cc.sweep_n_grid)
    L = float(cc.sweep_logodds_window)
    lam_g, fgrid = _logodds_grid(K, L)
    LG, LA = _log_fg(lam_g), _log1m_fg(lam_g)
    GLP = _regrid_global(cap["global_lp"], cc.sweep_n_grid, K, L) if cap["global_lp"] is not None else None
    IP = _regrid_global(cap["intron_prior"], cc.sweep_n_grid, K, L) if cap["intron_prior"] is not None else None
    prs = dbg["calibration_priors"]
    KAP = float(dbg["rna_sense_frac"])
    ODG, ODR = float(prs.gdna_strand_overdispersion), float(prs.rna_strand_overdispersion)
    up, un = np.asarray(st.u_pos, float), np.asarray(st.u_neg, float)
    single = fp_a ^ fn_a

    def psolve(idx, msg):
        out = np.zeros(idx.size)
        for j, i in enumerate(idx):
            i = int(i)
            pl = bool(fp_a[i]) and not bool(fn_a[i])
            nl = bool(fn_a[i]) and not bool(fp_a[i])
            fa = 1.0 - fgrid
            fpg = fa if pl else np.zeros_like(fa)
            fng = fa if nl else np.zeros_like(fa)
            psi = _mixture_strand_loglik(
                np.array([[up[i]]]), np.array([[up[i] + un[i]]]), fgrid[None, :], fpg[None, :],
                fng[None, :], KAP, ODG, ODR, np.array([[cap["fg_init"][i]]]),
                np.array([[cap["fpos_init"][i]]]), np.array([[cap["fneg_init"][i]]]))[0]
            psi = psi + (0.5 * LG if GLP is None else GLP[i]) + 0.5 * LA
            if IP is not None:
                psi = psi + IP[i]
            psi = psi - 0.5 * msg["cm_g"][i] * (LG - msg["mo_g"][i]) ** 2
            psi = psi - 0.5 * msg["cm_p"][i] * (LA - msg["mo_p"][i]) ** 2
            psi = psi - 0.5 * msg["cm_n"][i] * (LA - msg["mo_n"][i]) ** 2
            psi = psi - 0.5 * msg["c_tau"][i] * (lam_g - msg["lam"][i]) ** 2
            post = np.exp(psi - _lse(psi[None, :], axis=1, keepdims=True))
            out[j] = _posterior_median_fg(post, fgrid)[0]
        return out

    Gp_, Gn_, Rp_, Rn_ = _oracle_per_node(inp, chain)
    fo = np.where(Gp_ + Gn_ + Rp_ + Rn_ > _EPS,
                  (Gp_ + Gn_) / np.maximum(Gp_ + Gn_ + Rp_ + Rn_, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(n)])
    mass = np.asarray(cap["mass_global"])
    ok = np.flatnonzero(np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable_mask"], bool)
                        & single & (tau_own > 0))
    res = {}
    for v in VARIANTS:
        msg = combine(graft_no_r=v in ("graft_no_r", "both"), graft_s2t=v in ("graft_s2t", "both"))
        res[v] = psolve(ok, msg)
    fid = np.max(np.abs(res["shipped"] - np.asarray(cap["f_g"])[ok]))
    print(f"### {cond}   replay fidelity(shipped) = {fid:.2e}  {'OK' if fid < 1e-9 else '*** MISMATCH ***'}")
    for k, v in dict(oracle=fo[ok], mass=mass[ok], self=np.asarray(cap["fg_loc"])[ok],
                     cls=cls[ok], **{f"v_{v}": res[v] for v in VARIANTS}).items():
        OUT.setdefault(k, []).append(np.asarray(v))
    _ = is_amb


for c in a.conds:
    run(c)

D = {k: np.concatenate(v) for k, v in OUT.items()}
w = D["mass"].astype(float)
o = D["oracle"].astype(float)
print(f"\n{'=' * 130}\nSINGLE-HOP COUNTERFACTUAL — mass-weighted mwae over full-rank single-strand nodes "
      f"(n={w.size:,})\n{'=' * 130}")
print(f"{'stratum':<26}{'n':>7}{'self':>9}" + "".join(f"{v:>13}" for v in VARIANTS))
for lab, m in (("ALL full-rank", np.ones(w.size, bool)),
               ("  exon", D["cls"] == "exon"),
               ("  boundary", D["cls"] == "boundary"),
               ("  intron", D["cls"] == "intron"),
               ("oracle >= 0.90", o >= 0.90),
               ("oracle < 0.60", o < 0.60)):
    if not m.any():
        continue
    ww = w[m]
    print(f"{lab:<26}{int(m.sum()):>7}"
          f"{np.average(np.abs(D['self'][m].astype(float) - o[m]), weights=ww):>9.4f}"
          + "".join(f"{np.average(np.abs(D['v_' + v][m].astype(float) - o[m]), weights=ww):>13.4f}"
                   for v in VARIANTS))
print(f"\n{'oracle bin':<13}{'n':>7}{'oracle':>9}{'self':>9}" + "".join(f"{v:>13}" for v in VARIANTS))
for lo, hi in BINS:
    m = (o >= lo) & (o < hi)
    if not m.any():
        continue
    ww = w[m]
    print(f"[{lo:.2f},{hi:.2f})  {int(m.sum()):>7}{np.average(o[m], weights=ww):>9.3f}"
          f"{np.average(D['self'][m].astype(float), weights=ww):>9.3f}"
          + "".join(f"{np.average(D['v_' + v][m].astype(float), weights=ww):>13.3f}" for v in VARIANTS))
