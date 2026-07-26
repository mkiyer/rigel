"""RC3 — offline, bit-faithful replay of ``bp_solver._unified_solve``'s COMBINE (`_transport` + the fuse +
the three-stream handoff to ψ), reconstructed from the ``_capture`` dict.

Nothing here edits ``src/``. Everything ``_transport`` computes is recoverable from
``_capture['_uni_static']`` (the relay outputs + the per-face geometry) + ``_capture['_uni'][-1]``
(the last ρ-iteration's frames), so the whole combine can be re-run offline and every intermediate
inspected. ``validate()`` proves the replay is bit-identical to the shipped combine before any
counter is believed.

Used by ``rc3_hunt.py`` (the bug-class firing scan) and ``rc3_fix_ab.py`` (the ψ-replay A/B).
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.enrichment_frame import (  # noqa: E402
    mismatch_deflate,
    mismatch_gap,
    transfer_logvar,
)
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.node_init import own_composition_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1.0e-9
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}

_INDEX = None
_RA = None


def load(cond, refit=0):
    global _INDEX, _RA
    if _INDEX is None:
        _INDEX = TranscriptIndex.load(str(SUITE / "rigel_index"))
        _RA = RegionArrays.from_region_df(_INDEX.region_df, _INDEX.ref_name_to_id)
    cfg = PipelineConfig()
    inp = _scan_and_truth(
        SUITE, cond, _INDEX, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=refit)
    calmod.calibrate(
        inp["payload"], _RA, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    return dict(cond=cond, index=_INDEX, ra=_RA, inp=inp, dbg=dbg, cfg=cfg, cc=cc)


def oracle(ctx):
    Gp, Gn, Rp, Rn = _oracle_per_node(ctx["inp"], ctx["dbg"]["chain"])
    G, R = Gp + Gn, Rp + Rn
    return np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)


def build(ctx):
    """Reconstruct every static `_unified_solve` needs, then replay both `_transport` calls."""
    dbg = ctx["dbg"]
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    n = int(chain.n_nodes)

    ESP = (np.asarray(geo.eff_spl_left, np.float64), np.asarray(geo.eff_spl_right, np.float64))
    SP = (np.asarray(geo.spliced_pos_left, np.float64), np.asarray(geo.spliced_pos_right, np.float64))
    SN = (np.asarray(geo.spliced_neg_left, np.float64), np.asarray(geo.spliced_neg_right, np.float64))
    SPN = (np.asarray(geo.spliced_n_pos_left, np.float64), np.asarray(geo.spliced_n_pos_right, np.float64))
    SNN = (np.asarray(geo.spliced_n_neg_left, np.float64), np.asarray(geo.spliced_n_neg_right, np.float64))
    spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))

    M = np.asarray(us["M"], np.float64)
    E_g = np.asarray(us["E_g"], np.float64)
    E_r = np.asarray(us["E_r"], np.float64)
    og, op, on = us["og"], us["op"], us["on"]
    ex_a = np.asarray(us["is_exon"], bool)
    is_bnd_a = np.asarray(us["is_bnd"], bool)
    fp_a, fn_a = np.asarray(us["fp"], bool), np.asarray(us["fn"], bool)
    is_amb = fp_a & fn_a
    logvar_tot = np.asarray(us["logvar_tot"], np.float64)
    tau_own = np.asarray(us["tau_own"], np.float64)
    struct = np.asarray(us["struct_lock"], bool)
    fg_loc = np.asarray(cap["fg_loc"], np.float64)
    v_own_g, v_own_r = own_composition_logvar(fg_loc, tau_own, struct)
    v_own_lam = np.where(struct, 0.0, np.where(tau_own > _EPS, 1.0 / np.maximum(tau_own, _EPS), np.inf))

    fwd = tuple(us[k] for k in ("fwd_g", "fwd_p", "fwd_n", "fwd_pg", "fwd_pp", "fwd_pn",
                                "fwd_mg", "fwd_mp", "fwd_mn", "fwd_tau"))
    bwd = tuple(us[k] for k in ("bwd_g", "bwd_p", "bwd_n", "bwd_pg", "bwd_pp", "bwd_pn",
                                "bwd_mg", "bwd_mp", "bwd_mn", "bwd_tau"))
    li, ri = np.asarray(us["left"]), np.asarray(us["right"])
    vl, vr = li >= 0, ri >= 0
    sl, sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)
    rho_lf, rho_rf = np.asarray(uni["rho_lf"], np.float64), np.asarray(uni["rho_rf"], np.float64)

    def _pin_v(g, p, nn, pg_, pp_, pn_):
        sg = np.where(pg_ > 0.0, g, og)
        sp = np.where(pp_ > 0.0, p, op)
        sn = np.where(pn_ > 0.0, nn, on)
        s = sg * E_g + (sp + sn) * E_r
        k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
        return g * k, p * k, nn * k, k

    def _transport(src, valid, df, sf, arrs, dst_face_v, src_face_v):
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
        T = {}
        framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
        r = np.where(framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0))
        graft = ex_a & is_bnd_a[src] & valid
        gp = np.where(graft, spl_p_f[sf][src], 0.0)
        gn = np.where(graft, spl_n_f[sf][src], 0.0)
        tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

        def _dv(p, s2=s2t):
            return np.where(valid & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)
        ttau = _dv(tau, s2t)
        _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
        _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = is_bnd_a & ex_a[src] & valid
        T["peel"] = peel
        T["graft"] = graft
        T["r"] = r
        T["s2t"] = s2t
        T["tp_pre_peel"], T["tn_pre_peel"] = tp.copy(), tn.copy()
        T["mature_p"], T["mature_n"] = spl_p_f[df], spl_n_f[df]
        tp = np.where(peel, np.maximum(tp - spl_p_f[df], 0.0), tp)
        tn = np.where(peel, np.maximum(tn - spl_n_f[df], 0.0), tn)
        T["tp_post_peel"], T["tn_post_peel"] = tp.copy(), tn.copy()
        T["tpp_pre_gate"], T["tpn_pre_gate"] = tpp.copy(), tpn.copy()
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        T["tg_pre_pin"], T["tp_pre_pin"], T["tn_pre_pin"] = tg.copy(), tp.copy(), tn.copy()
        tg, tp, tn, kpin = _pin_v(tg, tp, tn, tpg, tpp, tpn)
        T["k_pin"] = kpin
        T["pin_sub_g"] = (tpg <= 0.0) & valid
        T["pin_sub_p"] = (tpp <= 0.0) & valid
        T["pin_sub_n"] = (tpn <= 0.0) & valid
        ttau_pre_gate = ttau.copy()
        ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
        T["lam_gate_killed"] = (ttau_pre_gate > 0.0) & (ttau <= 0.0)
        g_g, c_g = mismatch_gap(tg, og)
        g_p, c_p = mismatch_gap(tp, op)
        g_n, c_n = mismatch_gap(tn, on)
        T["contra_g"], T["contra_p"], T["contra_n"] = c_g, c_p, c_n
        T["tpg_pre_dl"], T["tmp_pre_dl"], T["tmn_pre_dl"] = tpg.copy(), tmp.copy(), tmn.copy()
        T["tmg_pre_dl"] = tmg.copy()
        tpg = mismatch_deflate(tpg, g_g, c_g, v_own_g)
        tpp = mismatch_deflate(tpp, g_p, c_p, v_own_r)
        tpn = mismatch_deflate(tpn, g_n, c_n, v_own_r)
        tmg = mismatch_deflate(tmg, g_g, c_g, v_own_g)
        tmp = mismatch_deflate(tmp, g_p, c_p, v_own_r)
        tmn = mismatch_deflate(tmn, g_n, c_n, v_own_r)
        g_R, c_R = mismatch_gap(tp + tn, op + on)
        ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, v_own_lam)
        T.update(tg=tg, tp=tp, tn=tn, tpg=tpg, tpp=tpp, tpn=tpn, tmg=tmg, tmp=tmp, tmn=tmn, ttau=ttau,
                 spc=_spc, snc=_snc, sp_mass=_sp, sn_mass=_sn,
                 sp_count=np.where(graft, SPN[sf][src], 0.0), sn_count=np.where(graft, SNN[sf][src], 0.0))
        return T

    A = _transport(sl, vl, 0, 1, fwd, rho_lf, rho_rf)
    B = _transport(sr, vr, 1, 0, bwd, rho_rf, rho_lf)

    def _fuse_v(a, pa, b, pb):
        p = pa + pb
        return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

    cg, cpg = _fuse_v(A["tg"], A["tpg"], B["tg"], B["tpg"])
    cp, cpp = _fuse_v(A["tp"], A["tpp"], B["tp"], B["tpp"])
    cn, cpn = _fuse_v(A["tn"], A["tpn"], B["tn"], B["tpn"])
    cm_g, cm_p, cm_n = A["tmg"] + B["tmg"], A["tmp"] + B["tmp"], A["tmn"] + B["tmn"]
    c_tau = A["ttau"] + B["ttau"]
    mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
    mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
    mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
    cR = cp + cn
    mo_R = np.log(np.maximum(cR * E_r / np.maximum(M, _EPS), _EPS))
    lam_msg = mo_g - mo_R
    c_tau_pre = c_tau.copy()
    c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
    tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1.0, 1.0)
    th_msg = np.arcsin(tau_tilt)
    th_prec = np.where(is_amb, cm_p + cm_n, 0.0)

    rt, _ = _node_region_type(chain, ctx["ra"])
    kind = np.asarray(chain.kind)
    cls = np.where(kind != 0, -1, rt)

    return dict(
        n=n, chain=chain, cap=cap, geo=geo, uni=uni, us=us, A=A, B=B,
        M=M, E_g=E_g, E_r=E_r, og=og, op=op, on=on, ex_a=ex_a, is_bnd_a=is_bnd_a,
        fp_a=fp_a, fn_a=fn_a, is_amb=is_amb, tau_own=tau_own, struct=struct,
        v_own_g=v_own_g, v_own_r=v_own_r, v_own_lam=v_own_lam,
        li=li, ri=ri, sl=sl, sr=sr, vl=vl, vr=vr, cls=cls, rt=rt, kind=kind,
        cg=cg, cp=cp, cn=cn, cR=cR, cpg=cpg, cpp=cpp, cpn=cpn,
        cm_g=cm_g, cm_p=cm_p, cm_n=cm_n, c_tau=c_tau, c_tau_pre=c_tau_pre,
        mo_g=mo_g, mo_p=mo_p, mo_n=mo_n, mo_R=mo_R, lam_msg=lam_msg,
        th_msg=th_msg, th_prec=th_prec, tau_tilt=tau_tilt,
        spl_p_f=spl_p_f, spl_n_f=spl_n_f, SP=SP, SN=SN, SPN=SPN, SNN=SNN, ESP=ESP,
        fg_loc=fg_loc, mass=np.asarray(cap["mass_global"], np.float64),
        solved=np.asarray(cap["f_g"], np.float64), solvable=np.asarray(cap["solvable"], bool),
    )


def validate(ctx, S):
    uni = S["uni"]
    out = {}
    for k in ("cg", "cp", "cn", "cm_g", "cm_p", "cm_n", "c_tau", "mo_g", "mo_p", "mo_n", "lam_msg"):
        a, b = np.asarray(S[k], np.float64), np.asarray(uni[k], np.float64)
        d = np.max(np.abs(a - b)) if a.size else 0.0
        out[k] = float(d)
    return out


if __name__ == "__main__":
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_present_capture_on"
    ctx = load(cond)
    S = build(ctx)
    print(f"# {cond}")
    for k, v in validate(ctx, S).items():
        print(f"  max|replay-shipped| {k:>8} = {v:.3e}")
