"""DERIVATION STEP 2 — replay the relay offline, exactly, then test PINNING AT EVERY HOP.

Step 1 established: (a) the reframe estimates the true capture step well (slope ~1, corr 0.92-0.96), so the
single-r channel mixing is the SMALLER effect; (b) the imputation premise itself fails by x2-3 (x3.22 even with
capture OFF); (c) analytically, one hop's mass violation is bounded by the eff-length ratio
(k = [Sum a_c^dst E_c]/[Sum a_c^src E_c], so ~1.04x on a contained region, 1.5x at a boundary crossing) --
therefore the measured p99 of 31-288x is ACCUMULATED drift, not per-hop error.

This script replays `_relay` bit-exactly from the captured statics (validated against the shipped fwd/bwd
arrays), then measures:
  A. the PER-HOP pin factor k_hop = M_i / Sum_c (incoming, reframed, routed context)_c E_c  -- is it bounded
     as the analysis predicts, and where does it exceed the bound (graft/peel edges break the identity)?
  B. what pinning at EVERY hop does to the accumulated overshoot, and to the relayed densities.

No solver edits; pure offline numpy.

    OMP_NUM_THREADS=1 python scratchpad/derive_2_relay_pin.py
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

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1.0e-9

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def relay(us, geom, *, pin: str = "none", forward: bool = True, r_true=None):
    """Bit-exact replay of `_relay`'s DENSITY arithmetic (the precisions do not affect the densities except
    through `_fuse`'s weights, which we reproduce). ``pin=True`` adds the candidate fix: after fusing, rescale
    the node's running belief so that Sum_c rho_c E_c == M -- the same operator `_pin_v` applies at the combine.
    pin: "none" | "fused" (scale the fused belief) | "context" (apply `_pin_v`'s OWN semantics to the
    incoming context before fusing -- substitute the node's own density for any component the context does not
    supply, so a PARTIAL claim stays partial, which is the whole point of the operator).
    Returns (rg, rp, rn, per-hop k, per-hop edge kind)."""
    og, op, on = us["og"].copy(), us["op"].copy(), us["on"].copy()
    pg_own, pp_own, pn_own = us["pg_own"], us["pp_own"], us["pn_own"]
    M, Eg, Er = us["M"], us["E_g"], us["E_r"]
    rho_l0, rho_r0 = us["rho_l0"], us["rho_r0"]
    ex_a, is_bnd = us["is_exon"].astype(bool), us["is_bnd"].astype(bool)
    fp_a, fn_a = us["fp"].astype(bool), us["fn"].astype(bool)
    lv = us["logvar_tot"]
    SP = (us["SP_l"], us["SP_r"])
    SN = (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geom.eff_spl_left, np.float64), np.asarray(geom.eff_spl_right, np.float64))
    spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    order = us["order"].tolist()
    nbr = us["left"] if forward else us["right"]
    dst_face, src_face = (rho_l0, rho_r0) if forward else (rho_r0, rho_l0)
    df, sf = (0, 1) if forward else (1, 0)
    seq = order if forward else order[::-1]

    rg, rp, rn = og.copy(), op.copy(), on.copy()
    pg, pp, pn = pg_own.copy(), pp_own.copy(), pn_own.copy()
    khop, kind_hop = [], []

    def _damp(p, s2):
        return 1.0 / (1.0 / p + s2) if p > 0.0 else 0.0

    def _fuse(a, pa, b, pb):
        p = pa + pb
        return ((pa * a + pb * b) / p, p) if p > _EPS else (a, 0.0)

    for i in seq:
        s = int(nbr[i])
        if s < 0:
            continue
        rho_src = src_face[s]
        r = (dst_face[i] / rho_src) if (rho_src > _EPS and dst_face[i] > _EPS) else 1.0
        if r_true is not None and r_true[i] > _EPS and r_true[s] > _EPS:
            r = r_true[i] / r_true[s]          # the ORACLE capture step (gDNA is uniform => pure capture)
        _gr = bool(ex_a[i] and is_bnd[s])
        s2t = 0.0 if _gr else (lv[i] + lv[s])
        gp = spl_p_f[sf][s] if _gr else 0.0
        gn = spl_n_f[sf][s] if _gr else 0.0
        tg, tp, tn = rg[s] * r, (rp[s] + gp) * r, (rn[s] + gn) * r
        tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)
        if _gr:
            _spc = SP[sf][s] / (1.0 + SP[sf][s] * s2t) if SP[sf][s] > _EPS else 0.0
            _snc = SN[sf][s] / (1.0 + SN[sf][s] * s2t) if SN[sf][s] > _EPS else 0.0
            tpp += _spc
            tpn += _snc
        if is_bnd[i] and ex_a[s]:
            tp = max(tp - spl_p_f[df][i], 0.0)
            tn = max(tn - spl_n_f[df][i], 0.0)
        if not fp_a[i]:
            tp, tpp = 0.0, 0.0
        if not fn_a[i]:
            tn, tpn = 0.0, 0.0
        # (A) the per-hop pin factor the INCOMING context would need
        ssum = tg * Eg[i] + (tp + tn) * Er[i]
        khop.append(M[i] / ssum if (ssum > _EPS and M[i] > _EPS) else np.nan)
        kind_hop.append("graft" if _gr else ("peel" if (is_bnd[i] and ex_a[s]) else "plain"))
        if pin == "context":   # `_pin_v` semantics, applied to the incoming context
            sg = tg if tpg > 0.0 else og[i]
            sp = tp if tpp > 0.0 else op[i]
            sn = tn if tpn > 0.0 else on[i]
            sv = sg * Eg[i] + (sp + sn) * Er[i]
            if sv > _EPS and M[i] > _EPS:
                kk_ = M[i] / sv
                tg, tp, tn = tg * kk_, tp * kk_, tn * kk_
        rg[i], pg[i] = _fuse(og[i], pg_own[i], tg, tpg)
        rp[i], pp[i] = _fuse(op[i], pp_own[i], tp, tpp)
        rn[i], pn[i] = _fuse(on[i], pn_own[i], tn, tpn)
        if pin == "fused":  # anchor the FUSED belief to this node's OBSERVED mass
            s2 = rg[i] * Eg[i] + (rp[i] + rn[i]) * Er[i]
            if s2 > _EPS and M[i] > _EPS:
                k = M[i] / s2
                rg[i], rp[i], rn[i] = rg[i] * k, rp[i] * k, rn[i] * k
    return rg, rp, rn, np.asarray(khop, np.float64), np.asarray(kind_hop)


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    us, geom = dbg["capture"]["_uni_static"], dbg["geometry"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, dbg["chain"])
    G, R = Gp + Gn, Rp + Rn
    M, Eg, Er = us["M"], us["E_g"], us["E_r"]

    rg0, rp0, rn0, k0, kh = relay(us, geom, pin="none")
    err = max(
        np.nanmax(np.abs(rg0 - us["fwd_g"])), np.nanmax(np.abs(rp0 - us["fwd_p"])),
        np.nanmax(np.abs(rn0 - us["fwd_n"])),
    )
    rg1, rp1, rn1, _, _ = relay(us, geom, pin="fused")
    rg2, rp2, rn2, _, _ = relay(us, geom, pin="context")
    ok = M > _EPS
    ov0 = ((rg0 * Eg + (rp0 + rn0) * Er) / np.maximum(M, _EPS))[ok]
    ov1 = ((rg1 * Eg + (rp1 + rn1) * Er) / np.maximum(M, _EPS))[ok]
    ov2 = ((rg2 * Eg + (rp2 + rn2) * Er) / np.maximum(M, _EPS))[ok]
    kk = k0[np.isfinite(k0) & (k0 > 0)]
    lk = np.abs(np.log(kk))

    print(f"\n{'=' * 98}\n{cond}     [replay fidelity vs shipped relay: max|Δρ| = {err:.3e}]\n{'=' * 98}")
    print(f"  A. PER-HOP pin factor |log k_hop|:  median {np.median(lk):.3f} (x{np.exp(np.median(lk)):.2f})   "
          f"p90 {np.percentile(lk, 90):.3f} (x{np.exp(np.percentile(lk, 90)):.2f})   "
          f"p99 {np.percentile(lk, 99):.3f} (x{np.exp(np.percentile(lk, 99)):.1f})   "
          f"max x{np.exp(lk.max()):.3g}")
    print(f"     fraction of hops needing >1.5x rescale: {(lk > np.log(1.5)).mean():.1%}   "
          f">10x: {(lk > np.log(10)).mean():.2%}")
    print(f"  B. ACCUMULATED overshoot Σ_c ρ_c E_c / M over nodes:")
    print(f"     {'':<12}{'median':>9}{'p90':>9}{'p99':>10}{'max':>12}{'% over 1':>10}{'% over 2':>10}")
    for lab, ov in (("as shipped", ov0), ("pin=fused", ov1), ("pin=context", ov2)):
        print(f"     {lab:<12}{np.median(ov):>9.3f}{np.percentile(ov, 90):>9.3f}"
              f"{np.percentile(ov, 99):>10.2f}{ov.max():>12.4g}{(ov > 1.001).mean():>9.1%}{(ov > 2).mean():>9.1%}")
    # does pinning move the relayed RNA toward the truth?
    sel = ok & (R + G > 0) & (Er > _EPS)
    tr = R[sel] / Er[sel]
    lk_all = np.abs(np.log(np.where(np.isfinite(k0) & (k0 > 0), k0, 1.0)))
    print("     per-hop |log k| by EDGE KIND (which hops break the mass identity?):")
    for t in ("plain", "graft", "peel"):
        m = (kh == t) & np.isfinite(k0) & (k0 > 0)
        if m.sum() < 5: continue
        print(f"        {t:<7} n={int(m.sum()):>6}  median x{np.exp(np.median(lk_all[m])):>6.2f}"
              f"   p90 x{np.exp(np.percentile(lk_all[m],90)):>8.2f}   >1.5x: {(lk_all[m]>np.log(1.5)).mean():>6.1%}")
    for lab, rp_, rn_ in (("as shipped", rp0, rn0), ("pin=fused", rp1, rn1), ("pin=context", rp2, rn2)):
        cl = (rp_ + rn_)[sel]
        print(f"     relayed RNA density vs oracle [{lab:<10}]: median claim/truth = "
              f"{np.median((cl + 1e-12) / (tr + 1e-12)):>8.2f}x   "
              f"mass-wt |Δ| = {np.average(np.abs(cl - tr), weights=M[sel]):.4f}")
