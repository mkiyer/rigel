"""DERIVATION #1 — THE τ-STREAM DL FORM.  Verdict: (B-λ), DL in λ space.  MC + exact algebra.

HANDOFF_5 §4.4 step 4 leaves open which gap the DerSimonian–Laird excess for the COMPOSITION (τ) stream
should be built from:

  (A-comp)  a PER-COMPONENT gap  G_c = log(t_c/o_c)          (c ∈ {g, +, −}; "which component?" undefined)
  (B-λ)     the λ-space gap      G_λ = (mo_g − mo_R) − (log f_g^own − log f_R^own) = G_g − G_R

THE ANSWER IS (B-λ), and the reason is not a preference — (A-comp) is a category error with a definite,
bounded, ALWAYS-UNSAFE sign:

  1. DL is a moment identity ON ONE ESTIMATOR.  The τ stream delivers exactly one estimator to ψ:
     ``lam_msg = mo_g − mo_R`` with precision ``c_tau`` (bp_solver ~line 591; simplex_logodds writes it as
     ``−½·lp_·(λ − lam_msg)²`` — a Gaussian on λ ITSELF, not on log f_c).  DL prices the bias OF THAT
     estimator, so the gap must be that estimator's residual against an independent estimate of the SAME
     quantity: G_λ vs λ_own.  G_c is the residual of a DIFFERENT estimator (the mode-fusion / measurement
     message on log f_c) and shares no moment identity with c_tau.
  2. UNITS.  ``v_msg,λ = 1/c_tau`` and ``v_own,λ = 1/τ_own`` are λ-space variances; ``G_c²`` is a log-f_c-space
     squared gap.  (A-comp) subtracts λ-space variances from a component-space gap — dimensionally incoherent.
  3. THE SIGN IS ALWAYS UNSAFE, and the factor is bounded in [1,4].  After ``_pin_v`` BOTH the message and the
     own belief are normalized compositions, so ``G_g`` and ``G_R`` have OPPOSITE signs and
     ``|G_g| + |G_R| = |G_λ|`` EXACTLY (§1 below proves it).  Hence
          G_λ²/4  ≤  max_c G_c²  ≤  G_λ²          ⇒   (A-comp) under-damps the τ stream by 1×–4×,
     with the worst case (4×) exactly when the composition gap straddles f_g = ½ — i.e. on the BIG mismatches
     the term exists to kill.  Under-damping = over-confident message = the §1 FORBIDDEN direction.

Also answered (the second question): ``v_msg,λ = 1/ttau`` is CORRECT and ``Var(log k) = v_g + v_R`` is WRONG
(§4).  The Schur τ IS a λ-space precision by construction, and mixing does NOT break the identity — because
DL is self-limiting: the delivered variance is ``max(v_msg, G² − v_own)``, so ``v_msg`` only matters where the
message already agrees with the own belief.  ``v_g + v_R`` double-counts the shared count (2/n, which cancels
in the ratio k) and drops the exact −1 rank-1 correlation.

    OMP_NUM_THREADS=1 python scratchpad/dlwf_tau_lambda.py [--draws 200000]

NOTE on the predecessor MC (``scratchpad/cliff_mc_3.py``): CONFIRMED BUGGY (§0).  Its ``transport_case``
anchors the destination's OWN mode at ``log(s_d_own·E_c/M_d)`` but the truth is
``log f_c^dst = log(ρ_c^dst·E_c/M_d) = log(s_c·ρ_tot_d·E_c/M_d)`` — it DROPS ``ρ_tot_d``, inflating every G by
``log ρ_tot_d`` (= +3.53 in cases A–D).  Its case-C "matched" gDNA row therefore reports σ²_DL = 12.3 against
a true MSE of 0.00023 — pure artifact.  Every anchor in THIS file is ``o_c = log f_c^own`` (a log-FRACTION of
the destination's own mass, exactly ``_ni.rho_c·E_c/M`` — which is what bp_solver's ``mo_c`` is compared to).
"""

from __future__ import annotations

import argparse

import numpy as np

from rigel.calibration.enrichment_frame import composition_logvar
from rigel.calibration.node_init import own_composition_logvar, strand_evidence
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

_EPS = 1.0e-9
rng = np.random.default_rng(20260725)

KAPPA = 0.99  # stranded library (the arm the anchor/stranded regression lives on)
OD_G = 0.0  # Poisson (ω = 0 by construction — the synthetic suite)
OD_R = 0.0
N_GRID = 60
N_GRID_SS = 601
L_WIN = 10.0


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §0  the predecessor's bug, confirmed by construction
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def sec0_confirm_bug():
    print("\n" + "=" * 108)
    print("§0  CONFIRM the cliff_mc_3.py own-anchor bug (before trusting anything it printed)")
    print("=" * 108)
    for name, fg_own, rho_tot_d, E in (("case A/C dst", 0.953, 34.0, 1.0),):
        M_d = rho_tot_d * E  # E_g=E_r=1 ⇒ M_d = ρ_tot_d
        buggy = np.log(fg_own * E / M_d)  # what cliff_mc_3 computes
        correct = np.log(fg_own * rho_tot_d * E / M_d)  # = log f_g^own
        print(f"  {name}: own mode  buggy={buggy:+.4f}   correct={correct:+.4f}   "
              f"Δ={buggy - correct:+.4f}  (= −log ρ_tot_d = {-np.log(rho_tot_d):+.4f})")
    print("  ⇒ G inflated by +log ρ_tot_d = +3.53 in cases A–D.  Its case-C 'matched' row reports")
    print("    σ²_DL(gDNA)=12.3 vs true MSE 0.00023 — a pure artifact of the dropped ρ_tot_d.")
    print("    Cases B/C per-component rows and case A's magnitudes are UNUSABLE; case E (b² sweep) is fine.")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §1  the exact post-pin geometry:  G_g·G_R < 0,  |G_g|+|G_R| = |G_λ|  ⇒  max_c G_c² ∈ [G_λ²/4, G_λ²]
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def sec1_gap_geometry():
    print("\n" + "=" * 108)
    print("§1  EXACT geometry of the gaps after ÷M (`_pin_v`): the per-component gaps are the RANK-1")
    print("    shadows of the ONE λ gap — max_c G_c² UNDER-states G_λ² by 1×–4×, always.")
    print("=" * 108)
    print("    (both message and own are normalized compositions ⇒ G_g and G_R have OPPOSITE signs and")
    print("     |G_g|+|G_R| = |G_λ| exactly; the mean Jacobian J̄_c = G_c/G_λ → ½ as the gap grows)")
    print(f"  {'f_g^own':>9} {'f_g^msg':>9} | {'G_g':>8} {'G_R':>8} {'G_λ':>8} | {'J̄_g':>6} {'J̄_R':>6} |"
          f" {'max G_c²':>9} {'G_λ²':>9} | {'under-damp':>10}")
    for fo, fm in ((0.953, 0.073), (0.953, 0.36), (0.985, 0.50), (0.50, 0.05),
                   (0.20, 0.80), (0.95, 0.90), (0.60, 0.61), (0.999, 0.001)):
        Gg = np.log(fm / fo)
        GR = np.log((1 - fm) / (1 - fo))
        Gl = Gg - GR
        mx = max(Gg * Gg, GR * GR)
        print(f"  {fo:9.3f} {fm:9.3f} | {Gg:+8.3f} {GR:+8.3f} {Gl:+8.3f} | {Gg / Gl:6.3f} {-GR / Gl:6.3f} |"
              f" {mx:9.4f} {Gl * Gl:9.4f} | {Gl * Gl / max(mx, 1e-12):9.2f}×")
    print("  ⇒ the ratio G_λ²/max_c G_c² is in [1,4] ALWAYS (=4 when the gap straddles f_g=½ — the big")
    print("    mismatches the term exists for).  (A-comp) on the τ stream is therefore over-confident by")
    print("    up to 4× exactly where over-confidence is most damaging.")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §2  the DL moment identity, per space (draws)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def sec2_moment_identity(N):
    print("\n" + "=" * 108)
    print("§2  DL MOMENT IDENTITY, per space.  E[G²] = b² + v_msg + v_own holds in the space the estimator")
    print("    lives in.  b̂²_λ recovers the λ-space bias; b̂²_c recovers only its J̄²-shrunk shadow.")
    print("=" * 108)
    tau_own, tau_msg = 1.6, 2.0
    v_own_l, v_msg_l = 1.0 / tau_own, 1.0 / tau_msg
    print(f"  τ_own={tau_own}  (v_own,λ={v_own_l:.4f})   τ_msg={tau_msg} (v_msg,λ={v_msg_l:.4f})")
    print(f"  {'λ_true':>7} {'b_λ':>6} | {'b̂²_λ (B-λ)':>12} {'b_λ²':>8} {'rel':>7} |"
          f" {'b̂²_g':>9} {'b̂²_R':>9} {'max':>9} | {'shortfall':>9}")
    for lam_true, b in ((3.0, -5.5), (3.0, -2.0), (0.0, -3.0), (3.0, 0.0), (-2.0, 4.0)):
        lam_own = lam_true + rng.normal(0.0, np.sqrt(v_own_l), N)
        lam_msg = lam_true + b + rng.normal(0.0, np.sqrt(v_msg_l), N)
        # the per-component gaps are formed from the COMPOSITIONS the two λ's imply (exactly as the solver
        # forms mo_c and o_c) — no linearization anywhere.
        fo, fm = _sig(lam_own), _sig(lam_msg)
        Gg, GR = np.log(fm / fo), np.log((1 - fm) / (1 - fo))
        Gl = lam_msg - lam_own
        # v_own,c / v_msg,c are the FOUNDATION Jacobians of the same τ (node_init.own_composition_logvar)
        vog, vor = (1 - fo) ** 2 * v_own_l, fo**2 * v_own_l
        vmg, vmr = (1 - fm) ** 2 * v_msg_l, fm**2 * v_msg_l
        bl = float(np.mean(np.maximum(0.0, Gl**2 - v_msg_l - v_own_l)))
        bg = float(np.mean(np.maximum(0.0, Gg**2 - vmg - vog)))
        br = float(np.mean(np.maximum(0.0, GR**2 - vmr - vor)))
        mx = max(bg, br)
        rel = "  n/a " if b == 0 else f"{abs(bl - b * b) / (b * b):6.1%}"
        sf = "  n/a " if b == 0 else f"{b * b / max(mx, 1e-9):7.2f}×"
        print(f"  {lam_true:7.2f} {b:+6.2f} | {bl:12.4f} {b * b:8.4f} {rel:>7} |"
              f" {bg:9.4f} {br:9.4f} {mx:9.4f} | {sf:>9}")
    print("  ⇒ (B-λ) is UNBIASED for b_λ² (rel ≤ ~1% at real mismatches; the b=0 row is the known SAFE")
    print("    positive truncation bias).  (A-comp)'s best component recovers only a 1–4× shrunk b², so")
    print("    feeding it to the λ-precision leaves the τ message that many times too confident.")


def _sig(x):
    return 1.0 / (1.0 + np.exp(-np.asarray(x, np.float64)))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §3  the real mechanics: reframe → graft → ÷M pin → three streams → DL → the REAL ψ solve
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def _logvar_tot(f_g, E_g, E_r, tau, n, lock):
    """bp_solver's logvar_tot (M5 σ²_transfer input), verbatim."""
    fg = float(np.clip(f_g, 0.0, 1.0))
    fgfr = fg * (1.0 - fg)
    if lock:
        var_fg = 0.0
    elif tau > _EPS:
        var_fg = min(fgfr * fgfr / tau, fgfr)
    else:
        var_fg = fgfr
    return float(composition_logvar(fg, E_g, E_r, var_fg, n))


def transport(dst, src, *, graft, s2t_law="M5"):
    """ONE edge of bp_solver._transport (single neighbour, so the fuse is a pass-through).

    Mirrors: reframe by r → graft the source's measured mature (a density AT the source) → damp every stream
    by σ²_transfer → add the graft's measurement count → ÷M pin → the three delivered objects
    (mo_g/mo_p/mo_n modes, cm_* measurement precisions, c_tau composition precision, lam_msg)."""
    r = dst["rho_tot_face"] / src["rho_tot_face"]
    lv_d = _logvar_tot(dst["fg_own"], dst["E_g"], dst["E_r"], dst["tau_own"], dst["n"], dst["lock"])
    lv_s = _logvar_tot(src["fg"], src["E_g"], src["E_r"], src["tau"], src["n"], src["lock"])
    if s2t_law == "M5":  # HANDOFF_5 §4.3 step 1: honest Var(log r), 0 on the graft
        s2t = 0.0 if graft else (lv_d + lv_s)
    elif s2t_law == "logr2":  # the shipped proxy being replaced
        s2t = float(np.log(max(r, _EPS)) ** 2)
    else:
        s2t = 0.0

    def dv(p):
        return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

    gp = src["spl_p"] if graft else 0.0
    gn = src["spl_n"] if graft else 0.0
    tg, tp, tn = src["rho_g"] * r, (src["rho_p"] + gp) * r, (src["rho_n"] + gn) * r
    tpg, tpp, tpn = dv(src["pg"]), dv(src["pp"]), dv(src["pn"])  # mode-fusion stream
    tmg, tmp, tmn = dv(src["mg"]), dv(src["mp"]), dv(src["mn"])  # measurement stream
    ttau = dv(src["tau"])  # composition stream
    if graft:  # the grafted mature is a COUNT — never τ-gated
        SP, SN = src["SP"], src["SN"]
        spc = SP / (1.0 + SP * s2t) if SP > _EPS else 0.0
        snc = SN / (1.0 + SN * s2t) if SN > _EPS else 0.0
        tpp, tpn = tpp + spc, tpn + snc
        tmp, tmn = tmp + spc, tmn + snc
    if not dst["free_p"]:
        tp, tpp, tmp = 0.0, 0.0, 0.0
    if not dst["free_n"]:
        tn, tpn, tmn = 0.0, 0.0, 0.0
    # ── _pin_v: the message is a claim about THIS node's mass (uninformed components filled from own) ──
    og, op, on = dst["og"], dst["op"], dst["on"]
    sg = tg if tpg > 0.0 else og
    sp = tp if tpp > 0.0 else op
    sn = tn if tpn > 0.0 else on
    s = sg * dst["E_g"] + (sp + sn) * dst["E_r"]
    k = dst["M"] / s if (s > _EPS and dst["M"] > _EPS) else 1.0
    tg, tp, tn = tg * k, tp * k, tn * k
    mo_g = np.log(max(tg * dst["E_g"] / dst["M"], _EPS))
    mo_p = np.log(max(tp * dst["E_r"] / dst["M"], _EPS))
    mo_n = np.log(max(tn * dst["E_r"] / dst["M"], _EPS))
    mo_R = np.log(max((tp + tn) * dst["E_r"] / dst["M"], _EPS))
    return dict(mo_g=mo_g, mo_p=mo_p, mo_n=mo_n, mo_R=mo_R, lam_msg=mo_g - mo_R,
                tg=tg, tp=tp, tn=tn, cm_g=tmg, cm_p=tmp, cm_n=tmn, c_tau=ttau,
                cpg=tpg, cpp=tpp, cpn=tpn, r=r, s2t=s2t)


def dl_deflate(msg, dst, law):
    """The HANDOFF_5 §4.3 step-2 DL deflation.  The measurement streams are per-component under EVERY law
    (that is settled — they ARE log-f_c estimators); the laws differ ONLY in the τ stream."""
    og, op, on = dst["og"], dst["op"], dst["on"]
    oR = op + on
    fg_o, tau_o, lock = dst["fg_own"], dst["tau_own"], dst["lock"]
    v_fg, v_fr = own_composition_logvar(np.array([fg_o]), np.array([tau_o]), np.array([lock]))
    v_own_g, v_own_R = float(v_fg[0]), float(v_fr[0])
    v_own_l = 0.0 if lock else (np.inf if tau_o <= _EPS else 1.0 / tau_o)
    # per-component gaps (E_c/M cancel ⇒ log(t_c/o_c) IS mo_c − o_c); guard both sides at _EPS
    def _G(t, o):
        return np.log(max(t, _EPS) / max(o, _EPS))
    G_g = _G(msg["tg"], og)
    G_R = _G(msg["tp"] + msg["tn"], oR)
    G_p = _G(msg["tp"], op) if op > _EPS else G_R
    G_n = _G(msg["tn"], on) if on > _EPS else G_R
    G_l = G_g - G_R

    def _ex(G2, v_msg, v_own):
        if not np.isfinite(v_own):
            return 0.0
        return max(0.0, G2 - v_msg - v_own)

    out = dict(msg)
    if law == "none":
        return out, dict(G_g=G_g, G_R=G_R, G_l=G_l, ex_tau=0.0)
    # measurement streams — per component, identical under all DL laws
    e_g = _ex(G_g**2, 1.0 / msg["cm_g"] if msg["cm_g"] > 0 else 0.0, v_own_g)
    e_p = _ex(G_p**2, 1.0 / msg["cm_p"] if msg["cm_p"] > 0 else 0.0, v_own_R)
    e_n = _ex(G_n**2, 1.0 / msg["cm_n"] if msg["cm_n"] > 0 else 0.0, v_own_R)
    out["cm_g"] = _damp(msg["cm_g"], e_g)
    out["cm_p"] = _damp(msg["cm_p"], e_p)
    out["cm_n"] = _damp(msg["cm_n"], e_n)
    # ── THE QUESTION: the τ (composition) stream ──
    v_msg_l = 1.0 / msg["c_tau"] if msg["c_tau"] > 0 else 0.0
    if law == "B-lam":
        ex_tau = _ex(G_l**2, v_msg_l, v_own_l)
    elif law == "A-g":
        ex_tau = _ex(G_g**2, v_msg_l, v_own_g)
    elif law == "A-R":
        ex_tau = _ex(G_R**2, v_msg_l, v_own_R)
    elif law == "A-max":
        ex_tau = max(_ex(G_g**2, v_msg_l, v_own_g), _ex(G_R**2, v_msg_l, v_own_R))
    else:
        raise ValueError(law)
    out["c_tau"] = _damp(msg["c_tau"], ex_tau)
    return out, dict(G_g=G_g, G_R=G_R, G_l=G_l, ex_tau=ex_tau)


def _damp(p, excess):
    """p_new = 1/(1/p + excess) — the §4.3 step-2 deflation.  excess = +inf ⇒ the stream is switched off."""
    if p <= 0.0:
        return 0.0
    return 1.0 / (1.0 / p + excess) if np.isfinite(excess) else 0.0


def psi_solve(dst, msgs):
    """THE REAL ψ (simplex_logodds), batched over the laws: strand BB + both Jeffreys arms + the gDNA/RNA
    MEASUREMENT messages + the SINGLE-λ composition message.  Returns f_g (posterior median) per law."""
    m = len(msgs)
    one = np.ones(m)
    fp = np.full(m, dst["free_p"])
    fn = np.full(m, dst["free_n"])
    dc = _solve_nodes_logodds_all(
        one * dst["u_pos"], one * dst["u_neg"], fp, fn,
        one * dst["M"], one * dst["S"],
        kappa=KAPPA, od_g=OD_G, od_r=OD_R, n_grid=N_GRID, L=L_WIN,
        n_tilt=N_GRID, n_grid_ss=N_GRID_SS,
        global_logprior=None,
        gdna_imp_mode=np.array([x["mo_g"] for x in msgs]),
        gdna_imp_prec=np.array([x["cm_g"] for x in msgs]),
        rna_imp_mode=(np.array([x["mo_p"] for x in msgs]), np.array([x["mo_n"] for x in msgs])),
        rna_imp_prec=(np.array([x["cm_p"] for x in msgs]), np.array([x["cm_n"] for x in msgs])),
        lam_imp_mode=np.array([x["lam_msg"] for x in msgs]),
        lam_imp_prec=np.array([x["c_tau"] for x in msgs]),
        theta_imp_mode=np.array([x.get("th_msg", 0.0) for x in msgs]),
        theta_imp_prec=np.array([x.get("th_prec", 0.0) for x in msgs]),
        fg_ref=one * dst["fg_ref"], fpos_ref=one * dst["fpos_ref"], fneg_ref=one * dst["fneg_ref"],
    )
    return np.asarray(dc.gdna_frac, np.float64)


def self_solve(dst):
    """The message-free ψ = the node's OWN belief (what `_ni` produces)."""
    z = dict(mo_g=0.0, mo_p=0.0, mo_n=0.0, mo_R=0.0, lam_msg=0.0,
             cm_g=0.0, cm_p=0.0, cm_n=0.0, c_tau=0.0)
    return float(psi_solve(dst, [z])[0])


LAWS = ("none", "A-g", "A-R", "A-max", "B-lam")


def make_dst(*, truth, n, fg_target, ambig=False, M=500.0, E_g=1000.0, E_r=900.0, cliff=407.0):
    """Build a destination node whose OWN belief and τ_own are SELF-CONSISTENT with its own counts:
    the counts are chosen so the message-free ψ lands at ``fg_target``, then ``fg_own`` is taken from the
    ACTUAL self-solve and ``τ_own`` from `node_init.strand_evidence` AT that f_g (single-strand only —
    AMBIG gets τ_own = 0 by the Schur gate).  This is exactly `build_node_init`'s own arithmetic."""
    p = KAPPA + fg_target * (0.5 - KAPPA)
    if ambig:
        u_pos, u_neg = 0.5 * n, 0.5 * n
    else:
        u_pos, u_neg = p * n, n - p * n
    dst = dict(M=M, S=0.0, E_g=E_g, E_r=E_r, n=n, u_pos=u_pos, u_neg=u_neg,
               free_p=True, free_n=bool(ambig), truth=truth, lock=False,
               fg_ref=fg_target, fpos_ref=(1 - fg_target) * (0.5 if ambig else 1.0),
               fneg_ref=(1 - fg_target) * 0.5 if ambig else 0.0)
    fg_self = self_solve(dst)
    dst["fg_own"] = fg_self
    if ambig:  # both strands live ⇒ the tilt is free ⇒ the strand cancels out of λ (Schur) ⇒ τ = 0
        dst["fp_own"] = dst["fn_own"] = 0.5 * (1 - fg_self)
        dst["tau_own"] = 0.0
    else:
        dst["fp_own"], dst["fn_own"] = 1 - fg_self, 0.0
        i_str, _ = strand_evidence(np.array([u_pos]), np.array([u_neg]), np.array([fg_self]),
                                   kappa=KAPPA, od_g=OD_G, od_r=OD_R,
                                   n_gdna_obs=1e6, n_rna_obs=1e6,
                                   is_region=np.array([True]), locked=np.array([False]))
        dst["tau_own"] = float(i_str[0])
    dst["og"] = dst["fg_own"] * M / E_g
    dst["op"] = dst["fp_own"] * M / E_r
    dst["on"] = dst["fn_own"] * M / E_r
    dst["rho_tot_face"] = M * (dst["fg_own"] / E_g + (1 - dst["fg_own"]) / E_r)
    dst["cliff"] = cliff
    return dst


def make_src(dst, *, fg, tau, SP=0.0, n=14.0, E=200.0, E_spl=150.0):
    """A source neighbour at ``dst.cliff``× lower total density (the enrichment cliff), composition ``fg``,
    optionally carrying a measured spliced junction of ``SP`` counts (the GRAFT)."""
    return dict(rho_g=fg * n / E, rho_p=(1 - fg) * n / E, rho_n=0.0,
                spl_p=SP / E_spl, spl_n=0.0, SP=SP, SN=0.0,
                pg=2.0, pp=2.0, pn=0.0, mg=0.0, mp=0.0, mn=0.0, tau=tau,
                n=n, fg=fg, E_g=E, E_r=E, lock=False,
                rho_tot_face=dst["rho_tot_face"] / dst["cliff"])


def run_scenario(title, dst, src, *, graft, note=""):
    print("\n" + "-" * 108)
    print(f"  {title}")
    if note:
        print(f"  {note}")
    fg_self = dst["fg_own"]
    print(f"  dst: TRUTH f_g={dst['truth']:.4f} | own self-solve f_g={fg_self:.4f}  τ_own={dst['tau_own']:.4g}"
          f"  n={dst['n']:g}   AMBIG={dst['free_p'] and dst['free_n']}")
    base = transport(dst, src, graft=graft, s2t_law="M5")
    shipped = transport(dst, src, graft=graft, s2t_law="logr2")
    print(f"  msg: cliff r={base['r']:.1f}  σ²_transfer(M5)={base['s2t']:.4g}  [(log r)²={shipped['s2t']:.3g}]"
          f"  delivered f_g^msg={np.exp(base['mo_g']):.4f}  λ_msg={base['lam_msg']:+.3f}  "
          f"streams: c_tau={base['c_tau']:.3g} cm_p={base['cm_p']:.3g}")
    rows, msgs, labels = [], [], []
    for law in LAWS:
        out, dg = dl_deflate(base, dst, law)
        msgs.append(out)
        labels.append("M5 only (no DL)" if law == "none" else f"M5 + DL[{law}]")
        rows.append((law, out, dg))
    msgs.append(shipped)
    labels.append("shipped (log r)²")
    rows.append(("logr2", shipped, dl_deflate(base, dst, "none")[1]))
    fgs = psi_solve(dst, msgs)
    d0 = rows[0][2]
    print(f"  gaps: G_g={d0['G_g']:+.3f}  G_R={d0['G_R']:+.3f}  G_λ={d0['G_l']:+.3f}   |   "
          f"G_g²={d0['G_g'] ** 2:.3f}  G_R²={d0['G_R'] ** 2:.3f}  G_λ²={d0['G_l'] ** 2:.3f}"
          f"  (G_λ²/max G_c² = {d0['G_l'] ** 2 / max(d0['G_g'] ** 2, d0['G_R'] ** 2, 1e-12):.2f}×)")
    print(f"  {'law':<18} {'excess_τ':>9} {'c_tau':>9} {'cm_p':>8} | {'fused f_g':>9} {'|err|':>7}"
          f" {'msg pull':>9}")
    res = {}
    for (law, out, dg), lab, fg in zip(rows, labels, fgs):
        err = abs(fg - dst["truth"])
        print(f"  {lab:<18} {dg['ex_tau']:9.3f} {out['c_tau']:9.4f} {out['cm_p']:8.4f} |"
              f" {fg:9.4f} {err:7.4f} {fg_self - fg:+9.4f}")
        res[law] = (fg, err, fg_self - fg)
    return res


def sec3_mechanics():
    print("\n" + "=" * 108)
    print("§3  THE REAL MECHANICS + THE REAL ψ.  reframe → graft → ÷M pin → 3 streams → DL → simplex_logodds.")
    print("    'msg pull' = own self-solve − fused: the DAMAGE the message does.  τ_own and the own belief are")
    print("    derived from the node's OWN counts (strand_evidence), so v_own = 1/τ_own is not a free knob.")
    print("=" * 108)

    print("\n### (i) ANCHOR — weak-but-CORRECT own belief vs a strong WRONG message across a ~400× cliff")
    print("###     dst EXON truth 0.985 ; src BOUNDARY f_g=0.36 + a 47-count junction (GRAFT ⇒ σ²_transfer=0)")
    summary = []
    for n_dst in (828.0, 311.0, 155.0):
        for tau_src in (2.0, 10.0):
            dst = make_dst(truth=0.985, n=n_dst, fg_target=0.953)
            src = make_src(dst, fg=0.36, tau=tau_src, SP=47.0)
            r = run_scenario(f"(i) anchor   n_dst={n_dst:g} (τ_own={dst['tau_own']:.3g})   τ_src={tau_src:g}",
                             dst, src, graft=True)
            summary.append(("anchor n=%g τs=%g" % (n_dst, tau_src), r))

    print("\n### (ii) MATCHED composition across a big cliff — the cliff ALONE must not damp")
    dst2 = make_dst(truth=0.97, n=828.0, fg_target=0.97)
    src2 = make_src(dst2, fg=dst2["fg_own"] - 0.005, tau=2.0, SP=0.0)
    r2 = run_scenario("(ii) matched composition, r≈407", dst2, src2, graft=False,
                      note="src composition ≈ dst own ⇒ G_λ ≈ 0 ⇒ every DL law must leave the message intact")
    summary.append(("matched-cliff", r2))

    print("\n### (iii) AMBIG destination (τ_own = 0 ⇒ v_own = ∞) — DL must be INERT")
    dst3 = make_dst(truth=0.30, n=200.0, fg_target=0.5, ambig=True)
    src3 = make_src(dst3, fg=0.36, tau=2.0, SP=47.0)
    r3 = run_scenario("(iii) AMBIG dst, τ_own=0", dst3, src3, graft=True,
                      note="the nodes where messages are the ONLY information — every DL law must be inert")
    summary.append(("AMBIG", r3))

    print("\n### (iv) STRESS — the message is RIGHT and the own belief is WRONG (the §4.5 worry).")
    print("###      Honest cost accounting: B-λ damps a CORRECT message harder than A-comp does.")
    dst4 = make_dst(truth=0.20, n=828.0, fg_target=0.90)  # own strand solve badly wrong (0.90 vs 0.20)
    src4 = make_src(dst4, fg=0.20, tau=2.0, SP=0.0)
    r4 = run_scenario("(iv) own WRONG (0.90) / msg RIGHT (0.20)", dst4, src4, graft=False,
                      note="a confidently-wrong own belief is NOT rescued by pass-0 under either law — by design")
    summary.append(("own-wrong", r4))

    print("\n  ── §3 SUMMARY (|err| vs truth) " + "─" * 74)
    print(f"  {'scenario':<24}" + "".join(f"{law:>11}" for law in LAWS + ("logr2",)))
    for name, r in summary:
        print(f"  {name:<24}" + "".join(f"{r[law][1]:11.4f}" for law in LAWS + ("logr2",)))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §4  v_msg,λ:  1/ttau (the relayed Schur τ)  vs  Var(log k) = v_g + v_R  — and the self-limiting cap
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def sec4_vmsg(N):
    print("\n" + "=" * 108)
    print("§4  IS v_msg,λ = 1/ttau LEGITIMATE?  (the Schur τ was not built as v_g+v_R)")
    print("=" * 108)
    print("  (a) The RANK-1 truth.  ρ_g and ρ_R from ONE source share M (⇒ the count cancels in the ratio) and")
    print("      move on ONE λ: d log f_g = (1−f_g)dλ, d log f_R = −f_g dλ ⇒ corr = −1 EXACTLY.")
    print(f"      {'f_g':>6} {'Var(log k) TRUE':>16} {'v_g+v_R':>10} {'ratio':>8}   (τ_src=2, n_src=14)")
    tau_s, n_s = 2.0, 14.0
    for f in (0.05, 0.25, 0.5, 0.75, 0.95, 0.999):
        v_g = (1 - f) ** 2 / tau_s + 1.0 / n_s
        v_R = f**2 / tau_s + 1.0 / n_s
        true = 1.0 / tau_s  # ((1−f)+f)²/τ ; the 1/n cancels in the ratio
        print(f"      {f:6.3f} {true:16.4f} {v_g + v_R:10.4f} {(v_g + v_R) / true:8.2f}×")
    print("      ⇒ v_g+v_R equals Var(log k) at NO interior f_g.  Its composition part ((1−f)²+f²)/τ understates")
    print("        the truth 1/τ by up to 2× (it drops the −1 covariance, worst at f=½), and it DOUBLE-COUNTS")
    print("        the shared count (+2/n, which cancels exactly in the ratio k).  On a thin source the +2/n")
    print("        over-states it instead.  Neither error is principled: 1/ttau IS the λ variance, use it.")
    print()
    print("  (b) Does mixing a λ-space gap with the Schur-τ precision break the DL identity?  NO — and here is")
    print("      the exact reason: the DELIVERED variance is self-limiting,")
    print("            v_deliv = v_msg + max(0, G² − v_msg − v_own) = MAX(v_msg, G² − v_own),")
    print("      which is MONOTONE in v_msg and has the v_msg-FREE floor G² − v_own.  So a mis-stated v_msg")
    print("      cannot make the message more confident than the observed gap allows.")
    print(f"      {'G²':>7} {'v_own':>7} | {'v_msg=1/ttau':>13} {'v_deliv':>9} | {'v_msg=v_g+v_R':>14} {'v_deliv':>9}"
          f" | {'cap 1/(G²−v_own)':>17}")
    v_own = 1.0 / 1.6
    for G2 in (30.8, 6.6, 2.0, 0.5, 0.05):
        for vm1, vm2 in ((0.5, 1.07),):
            d1 = max(vm1, G2 - v_own)
            d2 = max(vm2, G2 - v_own)
            cap = f"{1.0 / (G2 - v_own):17.4f}" if G2 > v_own else "     (no excess)"
            print(f"      {G2:7.3f} {v_own:7.4f} | {vm1:13.3f} {d1:9.4f} | {vm2:14.3f} {d2:9.4f} |"
                  f" {cap:>17}")
    print("      ⇒ identical wherever the excess is live (G² > v_msg+v_own).  v_msg only matters where the")
    print("        message AGREES with the own belief — precisely where it is harmless.  The only way to break")
    print("        the identity is to make the CAP wrong, i.e. to use a gap from a different space: (A-comp)")
    print("        caps at 1/(G_c²−v_own,c) = 1/(J̄²·(G_λ²−v_own,λ)) — 1–4× too confident (§1).")
    print()
    print("  (c) Residual bias of b̂² when the stated v_msg is wrong:  E[b̂²] = b² + (v_msg,true − v_msg,stated).")
    lam_t, b = 3.0, -5.5
    v_own_l, v_true = 1.0 / 1.6, 1.0 / 2.0
    for stated, tag in ((0.5, "1/ttau (correct)"), (1.07, "v_g+v_R (n=14)"), (0.05, "10× over-confident τ")):
        lam_own = lam_t + rng.normal(0.0, np.sqrt(v_own_l), N)
        lam_msg = lam_t + b + rng.normal(0.0, np.sqrt(v_true), N)
        G = lam_msg - lam_own
        bh = float(np.mean(np.maximum(0.0, G**2 - stated - v_own_l)))
        dl = float(np.mean(stated + np.maximum(0.0, G**2 - stated - v_own_l)))
        print(f"      v_msg stated={stated:5.3f} [{tag:>20}]  b̂²={bh:8.4f} (b²={b * b:.3f})  "
              f"v_deliv={dl:8.4f}   → precision {1.0 / dl:7.4f}")
    print("      ⇒ the delivered precision is ~invariant to the stated v_msg (the cap), and any UNDER-stated")
    print("        v_msg errs toward MORE damping — the §1-safe direction.  Use 1/ttau.")


# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
# §5  THE PIN-SAFETY INVARIANT — the governing principle, made exact
# ══════════════════════════════════════════════════════════════════════════════════════════════════════════
def sec5_pin_safety():
    print("\n" + "=" * 108)
    print("§5  THE PIN-SAFETY INVARIANT (the §1 governing principle as an exact inequality).")
    print("=" * 108)
    print("  A message OUT-WEIGHS the own belief iff v_deliv < v_own.  With v_deliv = max(v_msg, G²−v_own):")
    print("        message dominates  ⟺  v_msg < v_own  AND  G² < 2·v_own  ⟺  |G| < √2·σ_own.")
    print("  (B-λ) applies this in λ — the space the composition DOF lives in — so:")
    print("        ▶ NO message can dominate a node's own composition belief unless it agrees with it to")
    print("          within √2 own-σ.  A message that disagrees by more can never PIN the node.  That is")
    print("          exactly 'pass-0 must be weak and correctable', with no constant.")
    print("  (A-comp) applies it in log f_c, where G_c = J̄_c·G_λ and v_own,c = J_c²·v_own,λ, so the same")
    print("  crossover lands at  |G_λ| < σ_own,λ·√((1+J_c²)/J̄_c²)  ≥  √2·σ_own,λ  — ALWAYS looser:")
    print(f"      {'J̄ (mean)':>10} {'J (own)':>8} | {'B-λ threshold':>14} {'A-comp threshold':>17} {'looser by':>10}")
    for jb, j in ((0.5, 0.5), (0.5, 0.95), (0.557, 0.964), (0.27, 0.95), (0.9, 0.9), (1.0, 1.0)):
        tb = np.sqrt(2.0)
        ta = np.sqrt((1.0 + j * j) / (jb * jb))
        print(f"      {jb:10.3f} {j:8.3f} | {tb:11.3f}·σ {ta:14.3f}·σ {ta / tb:9.2f}×")
    print("  ⇒ (A-comp) lets a message that disagrees with the own belief by up to ~2.5 own-σ still PIN the")
    print("    node.  (B-λ) caps it at 1.41 own-σ, for every node, at every f_g.  Over-confidence is the")
    print("    one direction the owner forbids ⇒ (B-λ).")
    print()
    print("  EMPIRICAL CONFIRMATION on the real ψ — the τ stream ALONE (no junction ⇒ cm_*=0), sweeping the")
    print("  source composition.  'moved' = (λ_fused − λ_own)/(λ_msg − λ_own): 1 = the message wins, 0 = the")
    print("  own belief wins.  The 0.5 crossing is the pinning threshold.")
    dst = make_dst(truth=0.985, n=828.0, fg_target=0.953)
    sig = 1.0 / np.sqrt(dst["tau_own"])
    lam_own = np.log(dst["fg_own"] / (1 - dst["fg_own"]))
    print(f"  dst own f_g={dst['fg_own']:.4f}  τ_own={dst['tau_own']:.4g}  σ_own={sig:.3f}   "
          f"(B-λ threshold |G_λ|={np.sqrt(2) * sig:.3f})")
    print(f"  {'f_g^msg':>8} {'G_λ':>8} {'|G_λ|/σ':>8} |" + "".join(f"{law:>12}" for law in LAWS))
    for fgm in (0.94, 0.90, 0.80, 0.60, 0.40, 0.20, 0.05):
        src = make_src(dst, fg=fgm, tau=2.0, SP=0.0)
        base = transport(dst, src, graft=False, s2t_law="M5")
        msgs, dgs = [], []
        for law in LAWS:
            o, d = dl_deflate(base, dst, law)
            msgs.append(o)
            dgs.append(d)
        fgs = psi_solve(dst, msgs)
        moved = [(np.log(max(f, 1e-9) / max(1 - f, 1e-9)) - lam_own) / (base["lam_msg"] - lam_own)
                 for f in fgs]
        Gl = dgs[0]["G_l"]
        print(f"  {fgm:8.3f} {Gl:+8.3f} {abs(Gl) / sig:8.2f} |" + "".join(f"{m:12.3f}" for m in moved))
    print("  ⇒ (B-λ) crosses 'moved = ½' between 1.0 and 1.75·σ_own, as the √2·σ derivation predicts.  A-g is")
    print("    still at ~0.44 at 1.75·σ and hands over 36% of the node at 2.7·σ (B-λ: 0.28 / 0.13), and the")
    print("    A-variants DISAGREE WITH EACH OTHER — the same law in the wrong space, with an arbitrary choice.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=200_000)
    a = ap.parse_args()
    print("# DERIVATION #1 — the τ-stream DL form.   VERDICT: (B-λ).   draws=%d" % a.draws)
    sec0_confirm_bug()
    sec1_gap_geometry()
    sec2_moment_identity(a.draws)
    sec3_mechanics()
    sec4_vmsg(a.draws)
    sec5_pin_safety()


if __name__ == "__main__":
    main()
