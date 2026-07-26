"""RC3 — a FULL offline re-implementation of ``bp_solver._unified_solve`` (relay + combine + ψ + write-back),
driven entirely off the ``_capture`` dict, with each candidate fix behind a switch.

Validated bit-for-bit against the shipped ``f_g`` with all switches OFF, so any delta is attributable to the
switch alone. Nothing in ``src/`` is touched.

Switches (all structural / presence tests — no constants):
  ``peel_zero_gate``  a fully-consumed peel emits NOTHING (mode AND precision), instead of "ρ_R = 0" at live
                      precision (M3: ``ρ_ν ≤ 0`` is a prior TRUNCATION, not an emission).
  ``mature_split``    carry the measured MATURE density as its own relay component: the graft ADDS to it, the
                      peel REMOVES ALL of it (mature either splices or terminates at the exon edge — none of
                      it continues), instead of subtracting only the near junction's measured flux.
  ``pin_partial``     ``_pin_v`` may not renormalize a message onto a component it does not carry AND the
                      destination has no own density for: the message keeps its own scale (k = 1) there
                      rather than being pushed to the ``f_c = 1`` vertex.
  ``dl_absence``      the DL contradiction mask applies even where ``v_own = ∞``: a message asserting
                      ``ρ_c = 0`` against a destination that has the component live is killed regardless of
                      whether the destination has composition evidence.
  ``theta_gate``      the θ tilt message is emitted only when BOTH strand arms reached the node (else the
                      tilt is structurally undetermined, not "pinned at the wall").
  ``tau_break``       the composition τ stream is RESET (not relayed) through a node with no free RNA strand
                      — the same structural break that already zeroes the ``mp``/``mn`` measurement streams.
  ``spl_count``       the graft's measurement precision is the spliced COUNT (``spliced_n_*``), not the
                      fractional spliced MASS.
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from rc3_replay import build, load, oracle  # noqa: E402

from rigel.calibration.enrichment_frame import (  # noqa: E402
    mismatch_deflate,
    mismatch_gap,
    transfer_logvar,
)
from rigel.calibration.node_geometry import node_total_density  # noqa: E402
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402

_EPS = 1e-9
_RHO_ITERS = 2

SWITCHES = ("peel_zero_gate", "mature_split", "pin_partial", "dl_absence", "theta_gate",
            "tau_break", "spl_count", "mature_reach", "mu_reach", "frame_missing",
            "peel_m3", "graft_m2")


def solve(ctx, S, **sw):
    for k in sw:
        assert k in SWITCHES, k
    on = {k: bool(sw.get(k, False)) for k in SWITCHES}
    dbg = ctx["dbg"]
    chain, cap, geo, st = dbg["chain"], dbg["capture"], dbg["geometry"], dbg["statics"]
    cc = ctx["cc"]
    us = S["us"]
    n = S["n"]
    M, E_g, E_r = S["M"], S["E_g"], S["E_r"]
    og, op, onn = S["og"], S["op"], S["on"]
    pg_own = np.asarray(us["pg_own"], np.float64)
    pp_own = np.asarray(us["pp_own"], np.float64)
    pn_own = np.asarray(us["pn_own"], np.float64)
    mg_own = np.asarray(us["mg_own"], np.float64)
    tau_own = S["tau_own"]
    struct = S["struct"]
    ex_a, is_bnd_a, fp_a, fn_a = S["ex_a"], S["is_bnd_a"], S["fp_a"], S["fn_a"]
    is_amb = S["is_amb"]
    logvar_tot = np.asarray(us["logvar_tot"], np.float64)
    v_own_g, v_own_r, v_own_lam = S["v_own_g"], S["v_own_r"], S["v_own_lam"]
    spl_p_f, spl_n_f = S["spl_p_f"], S["spl_n_f"]
    SP, SN, SPN, SNN = S["SP"], S["SN"], S["SPN"], S["SNN"]
    order = [int(x) for x in np.asarray(us["order"])]
    li, ri = S["li"], S["ri"]
    sl, sr = S["sl"], S["sr"]
    vl, vr = S["vl"], S["vr"]
    accept_l = (SP[0] + SN[0]) > _EPS
    accept_r = (SP[1] + SN[1]) > _EPS
    # the graft's measurement precision source
    GP = (SPN if on["spl_count"] else SP)
    GN = (SNN if on["spl_count"] else SN)
    # the MATURE continuity gate — already computed by `build_node_statics`, currently used only to pick the
    # node prior: mature RNA is contiguous-exon-only, so a MATURE-derived claim has no authority at a node
    # where the mature molecule structurally does not exist (an intron, a splice-junction seam).
    mract_p = np.asarray(st.mrna_active_pos, bool)
    mract_n = np.asarray(st.mrna_active_neg, bool)

    def _rho_faces(fgc):
        ru, rw = node_total_density(chain, geo, fgc)
        rs = rw - ru
        return ru, ru + np.where(accept_l, rs, 0.0), ru + np.where(accept_r, rs, 0.0)

    fg_in = np.asarray(cap["fg_init"], np.float64)
    _, rho_l0, rho_r0 = _rho_faces(fg_in)

    def _damp(p, s2t):
        return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

    def _m2(rho_nu, rho_mu, p_nu, n_spl):
        """M2 (`enrichment_frame.graft_rna_logvar`): Var(log ρ_R) = w_ν²·v_ν + w_μ²·v_μ on a SUM."""
        tot = rho_nu + rho_mu
        if tot <= _EPS:
            return 0.0
        wn, wm = rho_nu / tot, rho_mu / tot
        v_mu = (1.0 / n_spl) if n_spl > 0.0 else np.inf
        if p_nu <= 0.0:
            v = np.inf if wn > 0.0 else wm * wm * v_mu
        else:
            v = wn * wn / p_nu + (wm * wm * v_mu if wm > 0.0 else 0.0)
        return 0.0 if (not np.isfinite(v) or v <= 0.0) else 1.0 / v

    def _m2v(rho_nu, rho_mu, p_nu, n_spl):
        tot = rho_nu + rho_mu
        wn = np.where(tot > _EPS, rho_nu / np.maximum(tot, _EPS), 0.0)
        wm = np.where(tot > _EPS, rho_mu / np.maximum(tot, _EPS), 0.0)
        v_mu = np.where(n_spl > 0.0, 1.0 / np.maximum(n_spl, _EPS), np.inf)
        v_nu = np.where(p_nu > 0.0, 1.0 / np.maximum(p_nu, _EPS), np.inf)
        a = np.where(wn > 0.0, wn * wn * np.where(np.isfinite(v_nu), v_nu, 0.0), 0.0)
        a = np.where((wn > 0.0) & ~np.isfinite(v_nu), np.inf, a)
        b = np.where(wm > 0.0, wm * wm * np.where(np.isfinite(v_mu), v_mu, 0.0), 0.0)
        b = np.where((wm > 0.0) & ~np.isfinite(v_mu), np.inf, b)
        v = a + b
        return np.where(np.isfinite(v) & (v > 0.0), 1.0 / np.maximum(v, _EPS), 0.0)

    def _m3v(T, rho_nu, p_T, n_spl):
        u = np.where(rho_nu > 0.0, T / np.maximum(rho_nu, _EPS), 0.0)
        v_mu = np.where(n_spl > 0.0, 1.0 / np.maximum(n_spl, _EPS), 0.0)
        v = np.where(p_T > 0.0, u * u / np.maximum(p_T, _EPS), np.inf) + (u - 1.0) ** 2 * v_mu
        return np.where(np.isfinite(v) & (v > 0.0), 1.0 / np.maximum(v, _EPS), 0.0)

    def _m3(T, rho_nu, p_T, n_spl):
        """M3 (`enrichment_frame.peel_rna_logvar`): Var(log ρ_ν) = u²·Var(log T) + (u−1)²·v_μ on a DIFFERENCE."""
        if p_T <= 0.0 or rho_nu <= 0.0:
            return 0.0
        u = T / rho_nu
        v_mu = (1.0 / n_spl) if n_spl > 0.0 else 0.0
        v = u * u / p_T + (u - 1.0) ** 2 * v_mu
        return 0.0 if v <= 0.0 else 1.0 / v

    def _fuse(a, pa, b, pb):
        p = pa + pb
        return ((pa * a + pb * b) / p, p) if p > _EPS else (a, 0.0)

    def _relay(seq, nbr, dst_face, src_face, df, sf):
        rg, rp, rn = og.copy(), op.copy(), onn.copy()
        mu_p = np.zeros(n)   # the MATURE (measured-spliced) sub-component of rp / rn
        mu_n = np.zeros(n)
        pg, pp, pn = pg_own.copy(), pp_own.copy(), pn_own.copy()
        mg, mp, mn = mg_own.copy(), np.zeros(n), np.zeros(n)
        tau = tau_own.copy()
        for i in seq:
            s = nbr[i]
            if s < 0:
                continue
            rho_src = src_face[s]
            _framed = rho_src > _EPS and dst_face[i] > _EPS
            r = (dst_face[i] / rho_src) if _framed else 1.0
            _gr = ex_a[i] and is_bnd_a[s]
            s2t = 0.0 if _gr else (logvar_tot[i] + logvar_tot[s])
            gp = spl_p_f[sf][s] if _gr else 0.0
            gn = spl_n_f[sf][s] if _gr else 0.0
            tg, tp, tn = rg[s] * r, (rp[s] + gp) * r, (rn[s] + gn) * r
            tmu_p, tmu_n = (mu_p[s] + gp) * r, (mu_n[s] + gn) * r
            tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)
            tmg, tmp, tmn = _damp(mg[s], s2t), _damp(mp[s], s2t), _damp(mn[s], s2t)
            ttau = _damp(tau[s], s2t)
            if on["frame_missing"] and not _framed:
                # no ρ_tot on one side ⇒ the reframe is UNDEFINED, not "identity". A density delivered in the
                # source's own frame is a claim about a different node's enrichment.
                tpg = tpp = tpn = tmg = tmp = tmn = ttau = 0.0
            if _gr:
                _spc = GP[sf][s] / (1.0 + GP[sf][s] * s2t) if GP[sf][s] > _EPS else 0.0
                _snc = GN[sf][s] / (1.0 + GN[sf][s] * s2t) if GN[sf][s] > _EPS else 0.0
                if on["graft_m2"]:  # M2: the GRAFT is a SUM ⇒ share-weighted (convex), not an inverse-var fuse
                    tpp = _m2(rp[s] * r, gp * r, tpp, SPN[sf][s])
                    tpn = _m2(rn[s] * r, gn * r, tpn, SNN[sf][s])
                else:
                    tpp += _spc
                    tpn += _snc
                tmp += _spc
                tmn += _snc
            if is_bnd_a[i] and ex_a[s]:  # EXON → boundary: PEEL the departing mature
                if on["mature_split"]:
                    sub_p, sub_n = tmu_p, tmu_n
                else:
                    sub_p, sub_n = spl_p_f[df][i], spl_n_f[df][i]
                np_, nn_ = tp - sub_p, tn - sub_n
                if on["peel_m3"]:  # M3: a DIFFERENCE ⇒ u-weighted (≥1); ρ_ν ≤ 0 is a PRIOR TRUNCATION
                    if np_ <= 0.0:
                        tp, tpp, tmp = 0.0, 0.0, 0.0
                    else:
                        tpp = _m3(tp, np_, tpp, SPN[df][i])
                        tp = np_
                    if nn_ <= 0.0:
                        tn, tpn, tmn = 0.0, 0.0, 0.0
                    else:
                        tpn = _m3(tn, nn_, tpn, SNN[df][i])
                        tn = nn_
                elif on["peel_zero_gate"]:
                    if np_ <= 0.0:
                        tp, tpp, tmp = 0.0, 0.0, 0.0
                    else:
                        tp = np_
                    if nn_ <= 0.0:
                        tn, tpn, tmn = 0.0, 0.0, 0.0
                    else:
                        tn = nn_
                else:
                    tp, tn = max(np_, 0.0), max(nn_, 0.0)
                tmu_p, tmu_n = 0.0, 0.0
            if not fp_a[i]:
                tp, tpp, tmp, tmu_p = 0.0, 0.0, 0.0, 0.0
            if not fn_a[i]:
                tn, tpn, tmn, tmu_n = 0.0, 0.0, 0.0, 0.0
            if on["mature_reach"]:  # the mp/mn stream is ENTIRELY mature-derived (mp_own ≡ 0)
                if not mract_p[i]:
                    tmp = 0.0
                if not mract_n[i]:
                    tmn = 0.0
            if on["mu_reach"]:  # the grafted mature DENSITY may not survive into a non-mature node either
                if not mract_p[i]:
                    tp, tmu_p = max(tp - tmu_p, 0.0), 0.0
                if not mract_n[i]:
                    tn, tmu_n = max(tn - tmu_n, 0.0), 0.0
            if on["tau_break"] and not (fp_a[i] or fn_a[i]):
                ttau = 0.0
            rg[i], pg[i] = _fuse(og[i], pg_own[i], tg, tpg)
            rp[i], pp[i] = _fuse(op[i], pp_own[i], tp, tpp)
            rn[i], pn[i] = _fuse(onn[i], pn_own[i], tn, tpn)
            # the mature sub-component rides the SAME fusion weights as its parent (it is part of tp/tn)
            wp = tpp / (pp_own[i] + tpp) if (pp_own[i] + tpp) > _EPS else 0.0
            wn = tpn / (pn_own[i] + tpn) if (pn_own[i] + tpn) > _EPS else 0.0
            mu_p[i] = wp * tmu_p
            mu_n[i] = wn * tmu_n
            mg[i] = mg_own[i] + tmg
            mp[i] = 0.0 + tmp
            mn[i] = 0.0 + tmn
            tau[i] = tau_own[i] + ttau
        return rg, rp, rn, pg, pp, pn, mg, mp, mn, tau, mu_p, mu_n

    fwd = _relay(order, li, rho_l0, rho_r0, 0, 1)
    bwd = _relay(order[::-1], ri, rho_r0, rho_l0, 1, 0)

    def _pin_v(g, p, nn_, pg_, pp_, pn_):
        sg = np.where(pg_ > 0.0, g, og)
        sp = np.where(pp_ > 0.0, p, op)
        sn = np.where(pn_ > 0.0, nn_, onn)
        s = sg * E_g + (sp + sn) * E_r
        k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
        if on["pin_partial"]:
            # a component the message does NOT carry and the destination has NO own density for cannot be
            # supplied by anybody — renormalizing as if it were ABSENT is the vertex assertion the pin exists
            # to avoid. Structural presence test, no constant.
            miss = ((pg_ <= 0.0) & (og <= _EPS)) | ((pp_ <= 0.0) & (op <= _EPS)) \
                | ((pn_ <= 0.0) & (onn <= _EPS))
            k = np.where(miss, 1.0, k)
        return g * k, p * k, nn_ * k

    def _transport(src, valid, df, sf, arrs, dst_face_v, src_face_v):
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau, mu_p, mu_n = arrs
        framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
        r = np.where(framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0))
        graft = ex_a & is_bnd_a[src] & valid
        gp = np.where(graft, spl_p_f[sf][src], 0.0)
        gn = np.where(graft, spl_n_f[sf][src], 0.0)
        tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
        tmu_p, tmu_n = (mu_p[src] + gp) * r, (mu_n[src] + gn) * r
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

        gate = (valid & framed) if on["frame_missing"] else valid

        def _dv(p, s2=s2t):
            return np.where(gate & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)
        ttau = _dv(tau, s2t)
        _sp = np.where(graft, GP[sf][src], 0.0)
        _sn = np.where(graft, GN[sf][src], 0.0)
        _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
        if on["graft_m2"]:
            tpp = np.where(graft, _m2v(rp[src] * r, gp * r, tpp, np.where(graft, SPN[sf][src], 0.0)), tpp)
            tpn = np.where(graft, _m2v(rn[src] * r, gn * r, tpn, np.where(graft, SNN[sf][src], 0.0)), tpn)
        else:
            tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = is_bnd_a & ex_a[src] & valid
        sub_p = tmu_p if on["mature_split"] else spl_p_f[df]
        sub_n = tmu_n if on["mature_split"] else spl_n_f[df]
        np_ = np.where(peel, tp - sub_p, tp)
        nn_ = np.where(peel, tn - sub_n, tn)
        if on["peel_m3"]:
            dead_p, dead_n = peel & (np_ <= 0.0), peel & (nn_ <= 0.0)
            tpp = np.where(peel & ~dead_p, _m3v(tp, np_, tpp, SPN[df]), tpp)
            tpn = np.where(peel & ~dead_n, _m3v(tn, nn_, tpn, SNN[df]), tpn)
            tp = np.where(dead_p, 0.0, np.where(peel, np_, tp))
            tn = np.where(dead_n, 0.0, np.where(peel, nn_, tn))
            tpp, tmp = np.where(dead_p, 0.0, tpp), np.where(dead_p, 0.0, tmp)
            tpn, tmn = np.where(dead_n, 0.0, tpn), np.where(dead_n, 0.0, tmn)
        elif on["peel_zero_gate"]:
            dead_p, dead_n = peel & (np_ <= 0.0), peel & (nn_ <= 0.0)
            tp = np.where(dead_p, 0.0, np.maximum(np_, 0.0))
            tn = np.where(dead_n, 0.0, np.maximum(nn_, 0.0))
            tpp, tmp = np.where(dead_p, 0.0, tpp), np.where(dead_p, 0.0, tmp)
            tpn, tmn = np.where(dead_n, 0.0, tpn), np.where(dead_n, 0.0, tmn)
        else:
            tp, tn = np.maximum(np_, 0.0), np.maximum(nn_, 0.0)
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        if on["mature_reach"]:
            tmp, tmn = np.where(mract_p, tmp, 0.0), np.where(mract_n, tmn, 0.0)
        if on["mu_reach"]:
            tp = np.where(mract_p, tp, np.maximum(tp - tmu_p, 0.0))
            tn = np.where(mract_n, tn, np.maximum(tn - tmu_n, 0.0))
        tg, tp, tn = _pin_v(tg, tp, tn, tpg, tpp, tpn)
        ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
        g_g, c_g = mismatch_gap(tg, og)
        g_p, c_p = mismatch_gap(tp, op)
        g_n, c_n = mismatch_gap(tn, onn)
        vg_, vr_, vl_ = v_own_g, v_own_r, v_own_lam
        if on["dl_absence"]:
            # a message asserting ρ_c = 0 where the destination HAS the component live is the b̂²→∞ limit
            # regardless of how well the destination knows its own composition.
            fin = np.zeros_like(v_own_g)
            vg_ = np.where(c_g, fin, v_own_g)
            vr_p = np.where(c_p, fin, v_own_r)
            vr_n = np.where(c_n, fin, v_own_r)
        else:
            vr_p = vr_n = v_own_r
        tpg = mismatch_deflate(tpg, g_g, c_g, vg_)
        tpp = mismatch_deflate(tpp, g_p, c_p, vr_p)
        tpn = mismatch_deflate(tpn, g_n, c_n, vr_n)
        tmg = mismatch_deflate(tmg, g_g, c_g, vg_)
        tmp = mismatch_deflate(tmp, g_p, c_p, vr_p)
        tmn = mismatch_deflate(tmn, g_n, c_n, vr_n)
        g_R, c_R = mismatch_gap(tp + tn, op + onn)
        if on["dl_absence"]:
            vl_ = np.where(c_g | c_R, np.zeros_like(v_own_lam), v_own_lam)
        ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, vl_)
        return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau

    def _fuse_v(a, pa, b, pb):
        p = pa + pb
        return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

    solvable = np.asarray(cap["solvable"], bool)
    f_cur = fg_in.copy()
    dc = None
    for _ in range(_RHO_ITERS):
        _, rho_lf, rho_rf = _rho_faces(f_cur)
        A = _transport(sl, vl, 0, 1, fwd, rho_lf, rho_rf)
        B = _transport(sr, vr, 1, 0, bwd, rho_rf, rho_lf)
        cg, cpg = _fuse_v(A[0], A[3], B[0], B[3])
        cp, cpp = _fuse_v(A[1], A[4], B[1], B[4])
        cn, cpn = _fuse_v(A[2], A[5], B[2], B[5])
        cm_g, cm_p, cm_n = A[6] + B[6], A[7] + B[7], A[8] + B[8]
        c_tau = A[9] + B[9]
        mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
        mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
        mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
        cR = cp + cn
        mo_R = np.log(np.maximum(cR * E_r / np.maximum(M, _EPS), _EPS))
        lam_msg = mo_g - mo_R
        c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
        tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1.0, 1.0)
        th_msg = np.arcsin(tau_tilt)
        th_prec = np.where(is_amb, cm_p + cm_n, 0.0)
        if on["theta_gate"]:
            th_prec = np.where((cp > _EPS) & (cn > _EPS), th_prec, 0.0)
        dc = _solve_nodes_logodds_all(
            st.u_pos, st.u_neg, st.free_pos, st.free_neg, st.mass_unspliced, st.mass_spliced,
            kappa=float(dbg["rna_sense_frac"]),
            od_g=float(dbg["calibration_priors"].gdna_strand_overdispersion),
            od_r=float(dbg["calibration_priors"].rna_strand_overdispersion),
            n_grid=cc.sweep_n_grid, L=cc.sweep_logodds_window,
            n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
            global_logprior=cap["global_lp"],
            gdna_imp_mode=mo_g, gdna_imp_prec=cm_g,
            rna_imp_mode=(mo_p, mo_n), rna_imp_prec=(cm_p, cm_n),
            lam_logprior=cap["intron_prior"],
            lam_imp_mode=lam_msg, lam_imp_prec=c_tau,
            theta_imp_mode=th_msg, theta_imp_prec=th_prec,
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
        )
        f_cur = np.clip(np.asarray(dc.gdna_frac, np.float64), 0.0, 1.0)
    return np.where(solvable, f_cur, fg_in)


if __name__ == "__main__":
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    ctx = load(cond)
    S = build(ctx)
    base = solve(ctx, S)
    print(f"# {cond}   max|offline - shipped| = {np.max(np.abs(base - S['solved'])):.3e}")
