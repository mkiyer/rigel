"""RC1 — offline EXACT replay of bp_solver's combine (`_transport`) + psi solve, with counterfactual knobs.

Reads nothing from src/ except pure functions; rebuilds the combine from the captured relay output so the
arithmetic can be re-run with `r` forced to 1, the graft removed, etc.  Validated bit-for-bit against the
shipped f_g before any counterfactual is applied.
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
from rigel.calibration.node_geometry import (  # noqa: E402
    _node_region_type,
    node_total_density,
)
from rigel.calibration.node_init import own_composition_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1.0e-9
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_RHO_ITERS = 2

_index = None
_ra = None


def load(cond):
    global _index, _ra
    if _index is None:
        _index = TranscriptIndex.load(str(SUITE / "rigel_index"))
        _ra = RegionArrays.from_region_df(_index.region_df, _index.ref_name_to_id)
    cfg = PipelineConfig()
    inp = _scan_and_truth(
        SUITE, cond, _index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], _ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    return _index, _ra, inp, dbg, cfg


class Ctx:
    """Everything the combine needs, pulled out of the capture + geometry."""

    def __init__(self, cond):
        index, ra, inp, dbg, cfg = load(cond)
        self.cond, self.index, self.ra, self.inp, self.dbg, self.cfg = cond, index, ra, inp, dbg, cfg
        cap = dbg["capture"]
        self.cap = cap
        self.us = us = cap["_uni_static"]
        self.uni = cap["_uni"]
        self.chain = chain = dbg["chain"]
        self.geom = g = dbg["geometry"]
        self.st = dbg["statics"]
        self.cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)

        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        self.Gp, self.Gn, self.Rp, self.Rn = Gp, Gn, Rp, Rn
        self.G, self.R = Gp + Gn, Rp + Rn
        self.fo = np.where(self.G + self.R > _EPS, self.G / np.maximum(self.G + self.R, _EPS), np.nan)

        self.M = np.asarray(us["M"], float)
        self.E_g = np.asarray(us["E_g"], float)
        self.E_r = np.asarray(us["E_r"], float)
        self.og, self.op, self.on = us["og"], us["op"], us["on"]
        self._pg_own = np.asarray(us["pg_own"], float)
        self._pp_own = np.asarray(us["pp_own"], float)
        self._pn_own = np.asarray(us["pn_own"], float)
        self.ex_a = np.asarray(us["is_exon"], bool)
        self.is_bnd = np.asarray(us["is_bnd"], bool)
        self.fp_a = np.asarray(us["fp"], bool)
        self.fn_a = np.asarray(us["fn"], bool)
        self.is_amb = self.fp_a & self.fn_a
        self.left = np.asarray(us["left"], np.int64)
        self.right = np.asarray(us["right"], np.int64)
        self.logvar_tot = np.asarray(us["logvar_tot"], float)
        self.tau_own = np.asarray(us["tau_own"], float)
        self.struct = np.asarray(us["struct_lock"], bool)
        self.mass = np.asarray(cap["mass_global"], float)
        self.fg_loc = np.asarray(cap["fg_loc"], float)
        self.solvable = np.asarray(cap["solvable"], bool)
        rt, _ = _node_region_type(chain, ra)
        self.rt = rt
        self.kind = np.asarray(chain.kind)
        self.cls = np.array(
            [CLS[int(rt[i])] if self.kind[i] == 0 else "boundary" for i in range(self.M.size)]
        )

        # per-face spliced mass + one-sided spliced eff-length + the derived mature density
        self.SP = (np.asarray(us["SP_l"], float), np.asarray(us["SP_r"], float))
        self.SN = (np.asarray(us["SN_l"], float), np.asarray(us["SN_r"], float))
        self.ESP = (np.asarray(g.eff_spl_left, float), np.asarray(g.eff_spl_right, float))
        self.spl_p_f = tuple(
            np.where(self.SP[k] > _EPS, self.SP[k] / np.maximum(self.ESP[k], _EPS), 0.0)
            for k in (0, 1)
        )
        self.spl_n_f = tuple(
            np.where(self.SN[k] > _EPS, self.SN[k] / np.maximum(self.ESP[k], _EPS), 0.0)
            for k in (0, 1)
        )
        self.accept_l = (self.SP[0] + self.SN[0]) > _EPS
        self.accept_r = (self.SP[1] + self.SN[1]) > _EPS

        self.fwd = (us["fwd_g"], us["fwd_p"], us["fwd_n"], us["fwd_pg"], us["fwd_pp"],
                    us["fwd_pn"], us["fwd_mg"], us["fwd_mp"], us["fwd_mn"], us["fwd_tau"])
        self.bwd = (us["bwd_g"], us["bwd_p"], us["bwd_n"], us["bwd_pg"], us["bwd_pp"],
                    us["bwd_pn"], us["bwd_mg"], us["bwd_mp"], us["bwd_mn"], us["bwd_tau"])

        self.v_own_g, self.v_own_r = own_composition_logvar(self.fg_loc, self.tau_own, self.struct)
        self.v_own_lam = np.where(
            self.struct, 0.0,
            np.where(self.tau_own > _EPS, 1.0 / np.maximum(self.tau_own, _EPS), np.inf),
        )
        n = self.M.size
        li, ri = self.left, self.right
        self.vl, self.vr = li >= 0, ri >= 0
        self.sl, self.sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)

        prs = dbg["calibration_priors"]
        self.psi_kw = dict(
            kappa=float(dbg["rna_sense_frac"]),
            od_g=float(prs.gdna_strand_overdispersion),
            od_r=float(prs.rna_strand_overdispersion),
            n_grid=self.cc.sweep_n_grid, L=self.cc.sweep_logodds_window,
            n_tilt=self.cc.sweep_n_tilt, n_grid_ss=self.cc.sweep_n_grid_single_strand,
            fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
            lam_logprior=cap["intron_prior"], global_logprior=cap["global_lp"],
        )
        self.psi_pos = (self.st.u_pos, self.st.u_neg, self.st.free_pos, self.st.free_neg,
                        self.st.mass_unspliced, self.st.mass_spliced)
        self.fg_init = np.asarray(cap["fg_init"], float)

    # ── the reframe frames, lazily from the current belief ──
    def rho_faces(self, fgc):
        ru, rw = node_total_density(self.chain, self.geom, fgc)
        rs = rw - ru
        return (ru, ru + np.where(self.accept_l, rs, 0.0), ru + np.where(self.accept_r, rs, 0.0))

    def _pin_v(self, g, p, n_, pg_, pp_, pn_):
        sg = np.where(pg_ > 0.0, g, self.og)
        sp = np.where(pp_ > 0.0, p, self.op)
        sn = np.where(pn_ > 0.0, n_, self.on)
        s = sg * self.E_g + (sp + sn) * self.E_r
        k = np.where((s > _EPS) & (self.M > _EPS), self.M / np.maximum(s, _EPS), 1.0)
        return g * k, p * k, n_ * k, k

    # ── the RELAY, offline (bp_solver._relay), with the same graft knob ──────────────────────────────
    def relay(self, seq, nbr, dst_face, src_face, df, sf, graft_mode="src", src_face_unspl=None,
              peel_mode="dst", meas_local=False):
        og, op, on = self.og, self.op, self.on
        pgo = self._pg_own
        ppo = self._pp_own
        pno = self._pn_own
        mgo = np.where(self.struct, pgo, 0.0)
        rg, rp, rn = og.copy(), op.copy(), on.copy()
        pg, pp, pn = pgo.copy(), ppo.copy(), pno.copy()
        mg = mgo.copy()
        mp = np.zeros_like(ppo)
        mn = np.zeros_like(pno)
        tau = self.tau_own.copy()
        lv = self.logvar_tot

        def _damp(p, s2):
            return 1.0 / (1.0 / p + s2) if p > 0.0 else 0.0

        def _fuse(a, pa, b, pb):
            p = pa + pb
            return ((pa * a + pb * b) / p, p) if p > _EPS else (a, 0.0)

        for i in seq:
            s = nbr[i]
            if s < 0:
                continue
            _gr = bool(self.ex_a[i]) and bool(self.is_bnd[s])
            rho_src = src_face_unspl[s] if (_gr and src_face_unspl is not None) else src_face[s]
            r = (dst_face[i] / rho_src) if (rho_src > _EPS and dst_face[i] > _EPS) else 1.0
            s2t = 0.0 if _gr else (lv[i] + lv[s])
            gp = self.spl_p_f[sf][s] if _gr else 0.0
            gn = self.spl_n_f[sf][s] if _gr else 0.0
            if graft_mode == "off":
                gp = gn = 0.0
            _pl = bool(self.is_bnd[i]) and bool(self.ex_a[s])
            rps, rns = rp[s], rn[s]
            if _pl and peel_mode == "src":
                rps = max(rps - self.spl_p_f[df][i], 0.0)
                rns = max(rns - self.spl_n_f[df][i], 0.0)
            if graft_mode == "src":
                tg, tp, tn = rg[s] * r, (rps + gp) * r, (rns + gn) * r
            else:  # 'dst' / 'post_pin' — the relay has no pin, so both mean "add after the reframe"
                tg, tp, tn = rg[s] * r, rps * r + gp, rns * r + gn
            tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)
            tmg, tmp, tmn = _damp(mg[s], s2t), _damp(mp[s], s2t), _damp(mn[s], s2t)
            if meas_local:  # a junction count measures its OWN flanking exon, not every node downstream
                tmp = tmn = 0.0
            ttau = _damp(tau[s], s2t)
            if _gr and graft_mode != "off":
                sp_, sn_ = self.SP[sf][s], self.SN[sf][s]
                _spc = sp_ / (1.0 + sp_ * s2t) if sp_ > _EPS else 0.0
                _snc = sn_ / (1.0 + sn_ * s2t) if sn_ > _EPS else 0.0
                tpp += _spc
                tpn += _snc
                tmp += _spc
                tmn += _snc
            if _pl and peel_mode != "src":
                tp = max(tp - self.spl_p_f[df][i], 0.0)
                tn = max(tn - self.spl_n_f[df][i], 0.0)
            if not self.fp_a[i]:
                tp, tpp, tmp = 0.0, 0.0, 0.0
            if not self.fn_a[i]:
                tn, tpn, tmn = 0.0, 0.0, 0.0
            rg[i], pg[i] = _fuse(og[i], pgo[i], tg, tpg)
            rp[i], pp[i] = _fuse(op[i], ppo[i], tp, tpp)
            rn[i], pn[i] = _fuse(on[i], pno[i], tn, tpn)
            mg[i] = mgo[i] + tmg
            mp[i] = 0.0 + tmp
            mn[i] = 0.0 + tmn
            tau[i] = self.tau_own[i] + ttau
        return rg, rp, rn, pg, pp, pn, mg, mp, mn, tau

    def build_relays(self, graft_mode="src", unspl_frame=False, peel_mode="dst", meas_local=False):
        order = [int(x) for x in np.asarray(self.us["order"])]
        rl0, rr0 = self.us["rho_l0"], self.us["rho_r0"]
        ru0 = self.us["rho_node0"] if unspl_frame else None
        fwd = self.relay(order, self.left, rl0, rr0, 0, 1, graft_mode, ru0, peel_mode, meas_local)
        bwd = self.relay(order[::-1], self.right, rr0, rl0, 1, 0, graft_mode, ru0, peel_mode,
                         meas_local)
        return fwd, bwd

    def transport(self, src, valid, df, sf, arrs, dst_face_v, src_face_v,
                  force_r1=False, no_graft=False, no_dl=False, trace=None,
                  graft_mode="src", src_face_unspl=None, peel_mode="dst"):
        """graft_mode: 'src' = shipped (mature joins the source RNA BEFORE the reframe, so it is pinned
        against the SOURCE's cliff-depleted gDNA); 'dst' = added after the reframe, still pinned;
        'post_pin' = an ABSOLUTE claim about the destination's RNA density, added after `_pin_v` and
        excluded from the pin normalizer."""
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
        graft0 = self.ex_a & self.is_bnd[src] & valid
        sfv = src_face_v[src]
        if src_face_unspl is not None:  # variant C: the mature is not carried by r, so drop it from the frame
            sfv = np.where(graft0, src_face_unspl[src], sfv)
        framed = valid & (sfv > _EPS) & (dst_face_v > _EPS)
        r = np.where(framed, dst_face_v / np.maximum(sfv, _EPS),
                     np.where(valid, 1.0, 0.0))
        r_true = r.copy()
        if force_r1:
            r = np.where(valid, 1.0, 0.0)
        graft = self.ex_a & self.is_bnd[src] & valid
        gp = np.where(graft, self.spl_p_f[sf][src], 0.0)
        gn = np.where(graft, self.spl_n_f[sf][src], 0.0)
        if no_graft or graft_mode == "off":
            gp = np.zeros_like(gp)
            gn = np.zeros_like(gn)
        peel0 = self.is_bnd & self.ex_a[src] & valid
        rp_s, rn_s = rp[src], rn[src]
        if peel_mode == "src":  # the junction flux lives in the EXON frame -> peel BEFORE the reframe
            rp_s = np.where(peel0, np.maximum(rp_s - self.spl_p_f[df], 0.0), rp_s)
            rn_s = np.where(peel0, np.maximum(rn_s - self.spl_n_f[df], 0.0), rn_s)
        if graft_mode == "src":
            tg, tp, tn = rg[src] * r, (rp_s + gp) * r, (rn_s + gn) * r
        elif graft_mode == "dst":
            tg, tp, tn = rg[src] * r, rp_s * r + gp, rn_s * r + gn
        else:  # post_pin / off
            tg, tp, tn = rg[src] * r, rp_s * r, rn_s * r
        s2t = transfer_logvar(self.logvar_tot, self.logvar_tot[src], graft)

        def _dv(p, s2=s2t):
            return np.where(valid & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)
        ttau = _dv(tau, s2t)
        _sp = np.where(graft, self.SP[sf][src], 0.0)
        _sn = np.where(graft, self.SN[sf][src], 0.0)
        if no_graft:
            _sp = np.zeros_like(_sp)
            _sn = np.zeros_like(_sn)
        _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = peel0
        if peel_mode != "src":
            tp = np.where(peel, np.maximum(tp - self.spl_p_f[df], 0.0), tp)
            tn = np.where(peel, np.maximum(tn - self.spl_n_f[df], 0.0), tn)
        tp, tpp, tmp = (np.where(self.fp_a, tp, 0.0), np.where(self.fp_a, tpp, 0.0),
                        np.where(self.fp_a, tmp, 0.0))
        tn, tpn, tmn = (np.where(self.fn_a, tn, 0.0), np.where(self.fn_a, tpn, 0.0),
                        np.where(self.fn_a, tmn, 0.0))
        pre_pin = (tg.copy(), tp.copy(), tn.copy())
        tg, tp, tn, pin_k = self._pin_v(tg, tp, tn, tpg, tpp, tpn)
        if graft_mode == "post_pin":
            # the mature is a MEASURED absolute density in the destination's own frame: add it AFTER the
            # pin, so it is never normalized against the source's (cliff-depleted) gDNA leg.
            tp = tp + gp
            tn = tn + gn
        ttau = np.where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0.0)
        if not no_dl:
            g_g, c_g = mismatch_gap(tg, self.og)
            g_p, c_p = mismatch_gap(tp, self.op)
            g_n, c_n = mismatch_gap(tn, self.on)
            tpg = mismatch_deflate(tpg, g_g, c_g, self.v_own_g)
            tpp = mismatch_deflate(tpp, g_p, c_p, self.v_own_r)
            tpn = mismatch_deflate(tpn, g_n, c_n, self.v_own_r)
            tmg = mismatch_deflate(tmg, g_g, c_g, self.v_own_g)
            tmp = mismatch_deflate(tmp, g_p, c_p, self.v_own_r)
            tmn = mismatch_deflate(tmn, g_n, c_n, self.v_own_r)
            g_R, c_R = mismatch_gap(tp + tn, self.op + self.on)
            ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, self.v_own_lam)
        if trace is not None:
            trace.update(r=r, r_true=r_true, graft=graft, gp=gp, gn=gn, s2t=s2t, pin=pin_k,
                         pre_pin=pre_pin, sp=_sp, sn=_sn, spc=_spc, snc=_snc, peel=peel,
                         src=src, sf=sf, df=df, rg_s=rg[src], rp_s=rp[src], rn_s=rn[src],
                         pp_s=pp[src], pg_s=pg[src], pn_s=pn[src], mp_s=mp[src], mn_s=mn[src],
                         dst_face=dst_face_v, src_face=src_face_v[src])
        return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau

    def run(self, force_r1=False, no_graft=False, no_dl=False, kill=None, trace_node=None,
            graft_mode="src", relay_graft=None, unspl_frame=False, peel_mode="dst",
            meas_local=False):
        """Replay the 2-iteration combine + psi. `kill` is a set of channel names to zero:
        {'lam','gdna','rna_p','rna_n','theta'}."""
        kill = kill or set()
        fwd, bwd = self.fwd, self.bwd
        if relay_graft is not None:
            fwd, bwd = self.build_relays(relay_graft, unspl_frame, peel_mode, meas_local)
        f_cur = self.fg_init.copy()
        tr = {}
        out = None
        for it in range(_RHO_ITERS):
            ru_f, rho_lf, rho_rf = self.rho_faces(f_cur)
            ruf = ru_f if unspl_frame else None
            ta = {} if trace_node is not None else None
            tb = {} if trace_node is not None else None
            A = self.transport(self.sl, self.vl, 0, 1, fwd, rho_lf, rho_rf,
                               force_r1, no_graft, no_dl, ta, graft_mode, ruf, peel_mode)
            B = self.transport(self.sr, self.vr, 1, 0, bwd, rho_rf, rho_lf,
                               force_r1, no_graft, no_dl, tb, graft_mode, ruf, peel_mode)
            ag, ap, an, apg, app, apn, amg, amp, amn, atau = A
            bg, bp, bn, bpg, bpp, bpn, bmg, bmp, bmn, btau = B

            def _fv(a, pa, b, pb):
                p = pa + pb
                return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

            cg, cpg = _fv(ag, apg, bg, bpg)
            cp, cpp = _fv(ap, app, bp, bpp)
            cn, cpn = _fv(an, apn, bn, bpn)
            cm_g, cm_p, cm_n = amg + bmg, amp + bmp, amn + bmn
            c_tau = atau + btau
            mo_g = np.log(np.maximum(cg * self.E_g / np.maximum(self.M, _EPS), _EPS))
            mo_p = np.log(np.maximum(cp * self.E_r / np.maximum(self.M, _EPS), _EPS))
            mo_n = np.log(np.maximum(cn * self.E_r / np.maximum(self.M, _EPS), _EPS))
            cR = cp + cn
            mo_R = np.log(np.maximum(cR * self.E_r / np.maximum(self.M, _EPS), _EPS))
            lam_msg = mo_g - mo_R
            c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
            tau_tilt = np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1, 1)
            th_msg = np.arcsin(tau_tilt)
            th_prec = np.where(self.is_amb, cm_p + cm_n, 0.0)
            if "lam" in kill:
                c_tau = np.zeros_like(c_tau)
            if "gdna" in kill:
                cm_g = np.zeros_like(cm_g)
            if "rna_p" in kill:
                cm_p = np.zeros_like(cm_p)
            if "rna_n" in kill:
                cm_n = np.zeros_like(cm_n)
            if "theta" in kill:
                th_prec = np.zeros_like(th_prec)
            dc = _solve_nodes_logodds_all(
                *self.psi_pos, **self.psi_kw,
                gdna_imp_mode=mo_g, gdna_imp_prec=cm_g,
                rna_imp_mode=(mo_p, mo_n), rna_imp_prec=(cm_p, cm_n),
                lam_imp_mode=lam_msg, lam_imp_prec=c_tau,
                theta_imp_mode=th_msg, theta_imp_prec=th_prec,
            )
            f_cur = np.clip(np.asarray(dc.gdna_frac, float), 0.0, 1.0)
            out = dict(cg=cg, cp=cp, cn=cn, cpg=cpg, cpp=cpp, cpn=cpn, cm_g=cm_g, cm_p=cm_p,
                       cm_n=cm_n, c_tau=c_tau, mo_g=mo_g, mo_p=mo_p, mo_n=mo_n, lam_msg=lam_msg,
                       th_msg=th_msg, th_prec=th_prec, fg=f_cur, ag=ag, ap=ap, an=an,
                       apg=apg, app=app, apn=apn, amp=amp, amn=amn, bg=bg, bp=bp, bn=bn,
                       bpg=bpg, bpp=bpp, bpn=bpn, bmp=bmp, bmn=bmn, rho_lf=rho_lf, rho_rf=rho_rf)
            if trace_node is not None:
                tr[it] = (ta, tb)
        fg = np.where(self.solvable, f_cur, self.fg_init)
        out["fg_final"] = fg
        out["trace"] = tr
        return out


def fidelity(ctx):
    base = ctx.run()
    ship = np.asarray(ctx.cap["f_g"], float)
    d = np.max(np.abs(base["fg_final"] - ship))
    u = ctx.uni[-1]
    dd = {k: float(np.nanmax(np.abs(np.asarray(base[k]) - np.asarray(u[k]))))
          for k in ("cp", "cn", "cg", "cm_p", "cm_n", "cm_g", "c_tau", "mo_g", "mo_p", "mo_n")}
    return d, dd, base
