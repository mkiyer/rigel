"""ADVERSARY #3 — the DerSimonian-Laird composition-mismatch plan, measured ON REAL SOLVER STATE.

Runs the production calibrate() with the diagnostic capture on cached ambig_dense_10mb conditions, then
reconstructs the PROPOSED DL deflation exactly as HANDOFF_5 §4.3 specifies it, using the captured
post-`_pin_v` message densities (ag/ap/an, bg/bp/bn), the captured per-message precisions (apg/app,
bpg/bpp) and the node's OWN self-solve (`_ni.rho_*` = og/op/on, `_ni.f_g`, `_ni.tau_lam`,
`_ni.struct_lock`).

    G_c      = log(t_c / o_c)
    v_own,g  = (1-fg_own)^2/tau_own ; v_own,R = fg_own^2/tau_own ; 0 if struct_lock ; +inf if tau_own<=EPS
    excess_c = max(0, G_c^2 - v_msg,c - v_own,c)      v_msg,c = 1/p_msg,c
    p_new    = 1/(1/p_old + excess_c)

We run with RIGEL_S2T_OFF=1 so the captured precisions are the UNDAMPED relay precisions (the shipped
(log r)^2 removed), then optionally re-apply the M5 Var(log r) = logvar_tot[dst]+logvar_tot[src] (0 on the
graft) analytically -- that is the state HANDOFF_5 step 1 restores.

Reports, per condition:
  R0  reach: which nodes DL can even touch (tau_own regimes; struct_lock n solvable)
  R1  the damping distribution + how many messages are killed
  R2  REGRET: for damped messages, was the message BETTER than the own belief (vs the oracle)?
  R3  the b=0 floor: damping applied where message and own AGREE
  R4  left/right asymmetry (both damped? opposite directions?)
  R5  iteration-2 feedback: f_cur at UNSOLVABLE nodes, and the rho_face swing it causes
  R6  partial messages (the _pin_v r-dependence hazard) + zero-density components (log(0) hazard)

    OMP_NUM_THREADS=1 python scratchpad/dlwf_adv_3_real.py [--conds a,b,c]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
SCRATCH = Path("/Users/mkiyer/proj/rigel/scratchpad/_dlwf_work")
_EPS = 1.0e-9

DEFAULT_CONDS = (
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",  # stranded (the regression arm)
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",  # unstranded x capture (the error-mass arm)
    "gdna_gdna5_ss_0.50_nrna_present_capture_verystrong",  # verystrong (the over-damped arm)
)


def _load(cond, index, cfg):
    from selfsolve_diag import _scan_and_truth

    return _scan_and_truth(SUITE, cond, index, cfg, SCRATCH, SUITE / "_selfsolve_cache")


def _run(inp, ra, cfg, s2t_off: bool):
    calmod = importlib.import_module("rigel.calibration.calibrate")
    if s2t_off:
        os.environ["RIGEL_S2T_OFF"] = "1"
    else:
        os.environ.pop("RIGEL_S2T_OFF", None)
    dbg = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"],
        ra,
        inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        cc,
        _debug=dbg,
    )
    return dbg


def _pct(x, w=None):
    x = np.asarray(x, float)
    if x.size == 0:
        return "n=0"
    q = np.percentile(x, [5, 25, 50, 75, 95])
    return f"p5={q[0]:.3g} p25={q[1]:.3g} p50={q[2]:.3g} p75={q[3]:.3g} p95={q[4]:.3g}"


def analyse(cond, index, cfg, ra):
    from flagship_interrogate import _oracle_per_node
    from rigel.calibration.node_geometry import _node_region_type

    inp = _load(cond, index, cfg)
    dbg = _run(inp, ra, cfg, s2t_off=True)
    chain = dbg["chain"]
    cap = dbg["capture"]
    st = cap["_uni_static"]
    uni = cap["_uni"]
    n_it = len(uni)

    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    Gt, Rt = Gp + Gn, Rp + Rn
    oracle_fg = np.where(Gt + Rt > _EPS, Gt / np.maximum(Gt + Rt, _EPS), np.nan)

    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    from rigel.calibration.node_chain import REGION

    is_reg = kind == REGION
    cls = np.where(is_reg, rt, -1)  # 0 intergenic 1 intron 2 exon -1 boundary
    CLSN = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary", 3: "other"}

    M = np.asarray(cap["mass_global"], float)
    E_g = np.asarray(cap["eff_global"], float)
    E_r = np.where(
        is_reg,
        np.asarray(cap["eff_rna_l"], float),
        np.asarray(cap["eff_rna_l"], float) + np.asarray(cap["eff_rna_r"], float),
    )
    og, op, on = st["og"], st["op"], st["on"]
    fg_own = np.asarray(cap["fg_loc"], float)
    tau_own = np.asarray(cap["_tau0_lam"], float)
    lock = np.asarray(st["struct_lock"], bool)
    solvable = np.asarray(cap["solvable"], bool)
    left, right = np.asarray(st["left"]), np.asarray(st["right"])
    logvar_tot = np.asarray(st["logvar_tot"], float)
    is_bnd = np.asarray(st["is_bnd"], bool)
    ex_a = np.asarray(st["is_exon"], bool)
    n = M.shape[0]

    print(f"\n{'='*100}\nCONDITION {cond}   nodes={n}")
    # ── R0. REACH ──────────────────────────────────────────────────────────────────────────────────────
    v_own_g = np.where(
        lock, 0.0, np.where(tau_own > _EPS, (1.0 - fg_own) ** 2 / np.maximum(tau_own, _EPS), np.inf)
    )
    v_own_R = np.where(
        lock, 0.0, np.where(tau_own > _EPS, fg_own**2 / np.maximum(tau_own, _EPS), np.inf)
    )
    has_own = np.isfinite(v_own_g)
    print("\n[R0] DL REACH (which destinations have a finite v_own -> DL can damp at all)")
    print(f"   solvable nodes                : {solvable.sum():7d}")
    print(f"   struct_lock nodes             : {lock.sum():7d}   (struct_lock & solvable = "
          f"{int((lock & solvable).sum())}  <-- the v_own=0 'full G^2' regime)")
    print(f"   solvable & tau_own>0 (DL live): {int((solvable & has_own).sum()):7d} "
          f"({100.0*(solvable & has_own).sum()/max(solvable.sum(),1):.1f}% of solvable)")
    print(f"   solvable & tau_own=0 (DL OFF) : {int((solvable & ~has_own).sum()):7d}")
    for c in (0, 1, 2, -1):
        m = solvable & (cls == c)
        if m.sum() == 0:
            continue
        live = m & has_own
        print(f"      {CLSN[c]:<11} solvable={m.sum():6d}  DL-live={live.sum():6d} "
              f"({100.0*live.sum()/m.sum():5.1f}%)  mass_share_DLlive="
              f"{100.0*M[live].sum()/max(M[m].sum(),_EPS):5.1f}%")

    # ── per-message DL reconstruction (iteration index -1 = the final combine) ────────────────────────
    out = {}
    for tag, src, td, tp_, tn_, ppg, ppp in (
        ("L", left, "ag", "ap", "an", "apg", "app"),
        ("R", right, "bg", "bp", "bn", "bpg", "bpp"),
    ):
        u = uni[-1]
        valid = src >= 0
        srcc = np.clip(src, 0, n - 1)
        graft = ex_a & is_bnd[srcc] & valid
        s2t_m5 = np.where(graft, 0.0, logvar_tot + logvar_tot[srcc])  # HANDOFF step-1 M5 restore
        rec = {}
        for cname, tk, pk, o_c, v_own in (
            ("g", td, ppg, og, v_own_g),
            ("+", tp_, ppp, op, v_own_R),
        ):
            t = np.asarray(u[tk], float)
            p_raw = np.asarray(u[pk], float)  # captured with RIGEL_S2T_OFF -> undamped
            p = np.where(p_raw > 0.0, 1.0 / (1.0 / np.maximum(p_raw, _EPS) + s2t_m5), 0.0)
            live = valid & (p > 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                G = np.log(np.maximum(t, _EPS) / np.maximum(o_c, _EPS))
            v_msg = np.where(p > 0.0, 1.0 / np.maximum(p, _EPS), np.inf)
            exc = np.maximum(0.0, G * G - v_msg - v_own)
            exc = np.where(np.isfinite(v_own), exc, 0.0)  # tau=0 -> no damping
            p_new = np.where(p > 0.0, 1.0 / (1.0 / np.maximum(p, _EPS) + exc), 0.0)
            damp = np.where(p > 0.0, p_new / np.maximum(p, _EPS), 1.0)
            rec[cname] = dict(
                t=t, p=p, G=G, v_msg=v_msg, exc=exc, damp=damp, live=live,
                zero_t=(t <= _EPS) & valid, zero_o=(o_c <= _EPS) & valid,
            )
        out[tag] = rec

    # ── R1. damping distribution ─────────────────────────────────────────────────────────────────────
    print("\n[R1] DAMPING of the mode-fusion stream  (p_new/p_old; 1.0 = untouched, ->0 = message killed)")
    for tag in ("L", "R"):
        for cname in ("g", "+"):
            r = out[tag][cname]
            m = r["live"] & solvable & has_own
            if m.sum() == 0:
                continue
            d = r["damp"][m]
            print(f"   msg {tag} comp {cname}: n={m.sum():6d}  {_pct(d)}")
            print(f"        killed(<0.1)={100.0*np.mean(d<0.1):5.1f}%   "
                  f"halved(<0.5)={100.0*np.mean(d<0.5):5.1f}%   untouched(=1)={100.0*np.mean(d>=0.999):5.1f}%"
                  f"   mass-weighted mean damp={np.average(d, weights=M[m]+_EPS):.3f}")

    # ── R2. REGRET: was the damped message better than the own belief? ───────────────────────────────
    def _logit(x):
        x = np.clip(x, 1e-4, 1 - 1e-4)
        return np.log(x / (1 - x))

    print("\n[R2] REGRET — among DAMPED gDNA messages, who was closer to the ORACLE?")
    for tag in ("L", "R"):
        r = out[tag]["g"]
        fg_msg = np.clip(r["t"] * E_g / np.maximum(M, _EPS), 1e-6, 1 - 1e-6)
        ok = r["live"] & solvable & has_own & np.isfinite(oracle_fg) & (M > 0)
        dm = ok & (r["damp"] < 0.5)
        if dm.sum() == 0:
            print(f"   msg {tag}: no damped messages")
            continue
        e_msg = np.abs(_logit(fg_msg) - _logit(oracle_fg))
        e_own = np.abs(_logit(fg_own) - _logit(oracle_fg))
        better = e_msg < e_own
        print(f"   msg {tag}: damped n={dm.sum():6d}  message-was-BETTER-than-own "
              f"{100.0*np.mean(better[dm]):5.1f}% of them "
              f"({100.0*M[dm & better].sum()/max(M[dm].sum(),_EPS):5.1f}% by mass)")
        print(f"        median |logit err|  msg={np.median(e_msg[dm]):.2f}  own={np.median(e_own[dm]):.2f}")
        kp = ok & (r["damp"] >= 0.5)
        if kp.sum():
            print(f"        (kept n={kp.sum():6d}: message-better {100.0*np.mean(better[kp]):5.1f}%)")

    # ── R3. b=0 floor: damping where message and own AGREE closely ──────────────────────────────────
    print("\n[R3] b=0 FLOOR — damping applied where the message and the own belief AGREE (|G| small)")
    for tag in ("L", "R"):
        r = out[tag]["g"]
        m = r["live"] & solvable & has_own
        for lo, hi in ((0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 1e9)):
            s = m & (np.abs(r["G"]) >= lo) & (np.abs(r["G"]) < hi)
            if s.sum() == 0:
                continue
            print(f"   msg {tag} |G| in [{lo},{hi}): n={s.sum():6d}  median damp={np.median(r['damp'][s]):.3f}"
                  f"  frac damped(<0.5)={100.0*np.mean(r['damp'][s]<0.5):5.1f}%"
                  f"  median v_own={np.median(np.where(np.isfinite(v_own_g[s]), v_own_g[s], np.nan)):.3g}"
                  f"  median v_msg={np.median(r['v_msg'][s]):.3g}")
        break  # L only, R is symmetric

    # ── R4. left/right asymmetry ────────────────────────────────────────────────────────────────────
    print("\n[R4] LEFT/RIGHT — do the two messages get damped together / in opposite directions?")
    rl, rr = out["L"]["g"], out["R"]["g"]
    both = rl["live"] & rr["live"] & solvable & has_own
    if both.sum():
        dl, dr = rl["damp"][both], rr["damp"][both]
        gl, gr = rl["G"][both], rr["G"][both]
        agree = np.sign(gl) == np.sign(gr)
        print(f"   both-live n={both.sum()}: BOTH damped(<0.5) {100.0*np.mean((dl<0.5)&(dr<0.5)):5.1f}%  "
              f"exactly one damped {100.0*np.mean((dl<0.5)^(dr<0.5)):5.1f}%  neither {100.0*np.mean((dl>=0.5)&(dr>=0.5)):5.1f}%")
        print(f"   G_L, G_R SAME sign (both disagree with own the same way): {100.0*np.mean(agree):5.1f}%")
        print(f"      of those, both damped: {100.0*np.mean((dl<0.5)&(dr<0.5)&agree)/max(np.mean(agree),1e-9):5.1f}%"
              f"  | opposite-sign both damped: {100.0*np.mean((dl<0.5)&(dr<0.5)&~agree)/max(np.mean(~agree),1e-9):5.1f}%")
        print(f"   corr(G_L,G_R)={np.corrcoef(gl, gr)[0,1]:+.3f}")

    # ── R5. iteration feedback ─────────────────────────────────────────────────────────────────────
    print("\n[R5] _RHO_ITERS FEEDBACK — f_cur at UNSOLVABLE nodes and the rho-face swing")
    fg_out_last = np.asarray(uni[-1]["fg_out"], float)
    uns = ~solvable
    print(f"   unsolvable nodes: n={uns.sum()}   f_cur(iter1 out) at unsolvable: "
          f"mean={fg_out_last[uns].mean():.4f} max={fg_out_last[uns].max():.4f} "
          f"(init belief f_g there = 1.0 for intergenic/seam)")
    print(f"   -> _rho_faces(f_cur) on iteration 2 treats those anchors as f_g={fg_out_last[uns].mean():.3f}")
    if n_it >= 2:
        a, b = np.asarray(uni[0]["rho_lf"], float), np.asarray(uni[1]["rho_lf"], float)
        ok = (a > _EPS) & (b > _EPS)
        sw = np.abs(np.log(b[ok] / a[ok]))
        print(f"   |log rho_lf(it2/it1)| over all nodes: {_pct(sw)}")
        ok_u = ok & uns
        if ok_u.sum():
            print(f"   ... restricted to UNSOLVABLE nodes: {_pct(np.abs(np.log(b[ok_u]/a[ok_u])))}")
        # what that does to G on edges whose SOURCE is unsolvable
        d1 = np.asarray(uni[0]["ag"], float)
        d2 = np.asarray(uni[1]["ag"], float)
        okg = (d1 > _EPS) & (d2 > _EPS) & (left >= 0)
        src_uns = uns[np.clip(left, 0, n - 1)]
        if (okg & src_uns).sum():
            print(f"   |ΔG_g| between iters, edges with UNSOLVABLE source: "
                  f"{_pct(np.abs(np.log(d2[okg & src_uns] / d1[okg & src_uns])))}")
        if (okg & ~src_uns).sum():
            print(f"   |ΔG_g| between iters, edges with SOLVABLE source:   "
                  f"{_pct(np.abs(np.log(d2[okg & ~src_uns] / d1[okg & ~src_uns])))}")

    # ── R6. partial messages + zero-density hazards ────────────────────────────────────────────────
    print("\n[R6] PARTIAL messages (the _pin_v r-dependence) + zero-density components")
    for tag in ("L", "R"):
        rg_, rp_ = out[tag]["g"], out[tag]["+"]
        valid = rg_["live"] | rp_["live"]
        partial = valid & ~(rg_["live"] & rp_["live"])
        print(f"   msg {tag}: valid n={valid.sum():6d}  PARTIAL (a component has 0 precision) "
              f"{100.0*partial.sum()/max(valid.sum(),1):5.1f}%")
        z = rg_["zero_t"] | rp_["zero_t"]
        print(f"        message component with density t_c == 0 (log(0) -> G^2 set by _EPS): "
              f"{int((z & solvable).sum()):6d} nodes "
              f"({100.0*(z & solvable & has_own).sum()/max((solvable&has_own).sum(),1):5.2f}% of DL-live)")
        zo = rg_["zero_o"] | rp_["zero_o"]
        print(f"        OWN density o_c == 0 at a valid destination: {int((zo & solvable).sum()):6d} "
              f"(G = +inf hazard; DL-live {(zo & solvable & has_own).sum()})")

    # ── R7. WHY was it damped? attribute the kill to the _EPS floor vs a genuine gap ────────────────
    print("\n[R7] KILL ATTRIBUTION — is the damping driven by the _EPS log-floor or by a real gap?")
    for tag in ("L", "R"):
        for cname, o_c in (("g", og), ("+", op)):
            r = out[tag][cname]
            m = r["live"] & solvable & has_own & (r["damp"] < 0.1)
            if m.sum() == 0:
                continue
            zt = r["zero_t"] & m
            zo = (np.asarray(o_c, float) <= _EPS) & m
            gen = m & ~zt & ~zo
            print(f"   msg {tag} comp {cname}: killed n={m.sum():6d} -> t_c==0 (EPS floor) {100.0*zt.sum()/m.sum():5.1f}%"
                  f" | o_c==0 (EPS floor) {100.0*zo.sum()/m.sum():5.1f}% | genuine gap {100.0*gen.sum()/m.sum():5.1f}%")
            if gen.sum():
                print(f"        genuine-gap |G| {_pct(np.abs(r['G'][gen]))}")

    # ── R8. per-strand v_own (the plan uses the RNA-TOTAL law on each strand) ───────────────────────
    print("\n[R8] PER-STRAND v_own — plan uses v_own,R = fg^2/tau for EACH strand; the true own")
    print("     Var(log f_pos) includes the TILT uncertainty. Ratio true/plan > 1 => plan OVER-damps.")
    vp_loc = np.asarray(cap["vp_loc"], float)
    vg_loc = np.asarray(cap["vg_loc"], float)
    fp_a = np.asarray(st["fp"], bool)
    fn_a = np.asarray(st["fn"], bool)
    amb = fp_a & fn_a
    for name, msk in (("AMBIG", amb & solvable & has_own), ("single-strand", (fp_a ^ fn_a) & solvable & has_own)):
        if msk.sum() == 0:
            continue
        ratio = np.where(v_own_R > _EPS, vp_loc / np.maximum(v_own_R, _EPS), np.nan)[msk & (v_own_R > _EPS)]
        rg_ = np.where(v_own_g > _EPS, vg_loc / np.maximum(v_own_g, _EPS), np.nan)[msk & (v_own_g > _EPS)]
        print(f"   {name:<14} n={msk.sum():6d}  Var(log f_pos)_solver / v_own,R(plan): {_pct(ratio[np.isfinite(ratio)])}")
        print(f"   {'':<14}          Var(log f_g)_solver  / v_own,g(plan): {_pct(rg_[np.isfinite(rg_)])}")

    # ── R9. DOUBLE USE of the own belief in the fusion (lambda-space Gaussian proxy) ────────────────
    print("\n[R9] DOUBLE USE — DL deflates the message with the own belief, then psi fuses the SAME own")
    print("     belief. Compare: (a) plan (own full weight, msg deflated); (b) honest random-effects")
    print("     (BOTH inflated by b^2); (c) no-DL. Gaussian proxy on lambda = logit f_g.")

    def _lam(x):
        return np.log(np.clip(x, 1e-6, 1 - 1e-6) / (1 - np.clip(x, 1e-6, 1 - 1e-6)))

    r = out["L"]["g"]
    fg_msg = np.clip(r["t"] * E_g / np.maximum(M, _EPS), 1e-6, 1 - 1e-6)
    m = r["live"] & solvable & has_own & np.isfinite(oracle_fg) & (M > 0) & (tau_own > _EPS)
    if m.sum():
        lam_o, lam_m, lam_t = _lam(fg_own[m]), _lam(fg_msg[m]), _lam(oracle_fg[m])
        # message lambda-precision: convert the per-component log-f precision to lambda via the Jacobian
        # d log f_g/d lambda = (1 - f_g); p_lambda = p_logf * (1-f_g)^2
        p_m = r["p"][m] * (1.0 - fg_msg[m]) ** 2
        b2 = r["exc"][m]
        t_o = tau_own[m]
        p_m_dl = 1.0 / (1.0 / np.maximum(p_m, _EPS) + b2)
        p_o_re = 1.0 / (1.0 / np.maximum(t_o, _EPS) + b2)  # honest RE inflates the own arm too
        for nm, po, pm in (("no-DL   ", t_o, p_m), ("plan-DL ", t_o, p_m_dl), ("honestRE", p_o_re, p_m_dl)):
            lam_post = (po * lam_o + pm * lam_m) / np.maximum(po + pm, _EPS)
            fgp = 1.0 / (1.0 + np.exp(-lam_post))
            err = np.abs(fgp - oracle_fg[m])
            print(f"   {nm}: post prec median={np.median(po+pm):9.4g}  |f_g err| mean={err.mean():.4f} "
                  f"median={np.median(err):.4f}  mass-wt={np.average(err, weights=M[m]):.4f}")
        oc = (t_o + p_m_dl) / np.maximum(p_o_re + p_m_dl, _EPS)
        print(f"   OVER-CONFIDENCE of the plan's fused precision vs honest random-effects: {_pct(oc)}")

    # ── R10. damping vs tau_own strength (risk 2: weak-but-nonzero own throttles good messages) ────
    print("\n[R10] DAMPING vs tau_own (v_own = (1-fg)^2/tau). Weak own belief = large v_own.")
    r = out["L"]["g"]
    base = r["live"] & solvable & has_own
    qs = np.percentile(tau_own[base], [0, 20, 40, 60, 80, 100])
    for i in range(5):
        s = base & (tau_own >= qs[i]) & (tau_own <= qs[i + 1])
        if s.sum() == 0:
            continue
        print(f"   tau_own in [{qs[i]:.3g},{qs[i+1]:.3g}]: n={s.sum():5d} median v_own={np.median(v_own_g[s]):.4g}"
              f" median |G|={np.median(np.abs(r['G'][s])):.3f} median damp={np.median(r['damp'][s]):.4f}"
              f" killed={100.0*np.mean(r['damp'][s]<0.1):5.1f}%")
    return dict(cond=cond, out=out, solvable=solvable, has_own=has_own, M=M, cls=cls,
                oracle=oracle_fg, fg_own=fg_own, E_g=E_g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conds", default=",".join(DEFAULT_CONDS))
    a = ap.parse_args()
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    SCRATCH.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    for cond in a.conds.split(","):
        analyse(cond.strip(), index, cfg, ra)


if __name__ == "__main__":
    main()
