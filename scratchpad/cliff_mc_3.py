"""CROSS-CLIFF PRECISION — derivation agent #3 MC.

QUESTION (my assigned angle). Is `σ²_transfer = Var(log r)` the RIGHT object for the cross-cliff precision
damping, or are the graft-exemption + the 1/n-smallness the two bugs? Should the PRECISION always pay
Var(log r), and should Var(log r) carry a composition-mismatch term beyond the counting 1/n?

MY ANSWER (derived here, MC-validated below).
  (i)  The MODE reframe cancellation is exact and the graft MODE exemption is correct — KEEP it.
  (ii) `Var(log r)` (the enrichment-ratio *measurement* noise, ~1/n) genuinely CANCELS in the composition on a
       matched reframe, so making the precision "pay Var(log r)" on the graft is NOT the fix — the doc is right
       that it cancels. The real missing object is a SEPARATE, PER-COMPONENT composition-MISMATCH variance.
  (iii)The exact delivered-message mode error for component c is
              ε_c = mo_c − log f_c^dst,true = log(s_c^src) − log(s_c^dst,true) ,   s_c = ρ_c/ρ_tot  (the SHARE)
       i.e. the whole error is the source-vs-dest SHARE mismatch of component c. ρ_tot^dst and M_dst cancel
       (confirms M4). The reframe r does NOT appear — it is fully absorbed into the shares. So the honest
       delivered variance is
              σ²_c = Var(log s_c^src)          [source sampling, = the doc's v_source]
                   + b_c² ,  b_c = log s_c^src − log s_c^dst,true   [the composition-mismatch BIAS, per component]
       This is MSE = variance + bias². The bug is that the shipped law uses only the first term.
  (iv) b_c is a PRIOR quantity (physical between-region composition variability) — not knowable from counts in
       pass-0. But it is ESTIMABLE at solve time from the conflict between the message and the destination's OWN
       belief, by the DerSimonian–Laird between-source (random-effects) estimator — NO magic number:
              b_c² ≈ max(0, G_c² − v_src,c − v_own,c) ,   G_c = mo_c^msg − mo_c^own  (the observed share gap)
       → σ²_c,delivered = v_src,c + max(0, G_c² − v_src,c − v_own,c).

WHY THIS SATISFIES EVERY CONSTRAINT.
  * PER-COMPONENT and majority-dominated-aware: across a cliff the MAJORITY share barely moves (small G) while
    the MINORITY share collapses (huge G) — so the minority message is damped hard and the majority is spared.
    A single scalar Var(log r) applied to all components CANNOT express this asymmetry; that is why the doc's
    object structurally could not fix the minority inflation.
  * Matched transport (same composition, any cliff) → G≈0 → excess≈0 → precision PRESERVED (full).
  * UNSTRANDED / AMBIG destination has NO own composition belief → v_own = ∞ → excess = max(0, G²−∞) = 0 →
    the message propagates at full v_src → the M5 unstranded win is UNTOUCHED (no over-damping).
  * No magic number: the bias² is a method-of-moments (DerSimonian–Laird) between-source variance.

This MC constructs the actual reframe/÷M mechanics (densities → r → mode) and shows:
  A) the anchor case (small mismatched RNA source across a ~425× cliff) — the shipped law is confidently WRONG
     (overconfidence factor ~10²–10³); the DL law delivers precision ~0.
  B) a matched same-scale transport — precision preserved (DL ≈ shipped, both = v_src).
  C) a matched-composition BIG cliff — precision preserved (the cliff alone does NOT damp).
  D) the AMBIG destination (no own belief) — excess = 0, message propagates.
  E) DL estimator recovers the true bias² (the load-bearing claim).

    OMP_NUM_THREADS=1 python scratchpad/cliff_mc_3.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260724)


def _lognormal(mean, var_log, size):
    s = np.sqrt(max(float(var_log), 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def _fmt(x):
    return f"{x:12.5g}"


def transport_case(
    name,
    N,
    *,
    # destination TRUTH (the region we are imputing into)
    fg_d_true,
    rho_tot_d,
    # destination OWN belief (its self-solve; None ⇒ AMBIG, no own composition info ⇒ v_own = ∞)
    fg_d_own,
    v_own_g,
    v_own_R,
    # source belief (the neighbour that emits the message)
    fg_s,
    rho_tot_s,
    v_src_g,  # source gDNA-share log-variance (its sampling precision)
    v_src_R,  # source RNA-share log-variance (e.g. 1/n_spl for a measured junction)
    tau_own_lam=None,  # dst own belief precision on λ=logit f_g (for the fusion punchline)
    E_g=1.0,
    E_r=1.0,
):
    """One source→dest transport. E_g=E_r=1 ⇒ share s_c = fraction f_c (the eff-lengths cancel in ε_c anyway;
    the derivation is exact for any E). Returns nothing; prints the per-component verdict."""
    # ── truth & source densities (shares × total). With E_g=E_r these shares equal the f_c fractions. ──
    s_g_d, s_R_d = fg_d_true, 1.0 - fg_d_true
    s_g_s, s_R_s = fg_s, 1.0 - fg_s
    rho_g_d, rho_R_d = s_g_d * rho_tot_d, s_R_d * rho_tot_d
    rho_g_s0, rho_R_s0 = s_g_s * rho_tot_s, s_R_s * rho_tot_s
    r0 = rho_tot_d / rho_tot_s
    # dst mass identity M_d = Σ ρ_c^dst·E_c (so f_c^dst are true fractions). M_d CANCELS in the mode error
    # (verified: A–E are invariant to M_d), but the fusion punchline needs the absolute modes to be real
    # fractions, so we use the honest mass rather than an arbitrary constant.
    M_d = rho_g_d * E_g + rho_R_d * E_r

    print(f"\n── {name}")
    print(f"   dst truth f_g={fg_d_true:.4f}  ρ_tot(dst)={rho_tot_d:g}    "
          f"src f_g={fg_s:.4f}  ρ_tot(src)={rho_tot_s:g}    cliff r≈{r0:.1f}")

    # ── DRAWS: source sampling noise on each component's share (its belief), + dst own-belief noise. ──
    rho_g_s = _lognormal(rho_g_s0, v_src_g, N)
    rho_R_s = _lognormal(rho_R_s0, v_src_R, N)
    rho_tot_s_d = rho_g_s + rho_R_s
    r = rho_tot_d / rho_tot_s_d  # the measured reframe (dst total is well-counted ⇒ treated fixed)

    # the delivered ψ modes: mo_c = log(ρ_c^src · r · E_c / M_d)
    mo_g = np.log(rho_g_s * r * E_g / M_d)
    mo_R = np.log(rho_R_s * r * E_r / M_d)
    # the TRUTH the message is estimating: log f_c^dst
    tru_g = np.log(rho_g_d * E_g / M_d)
    tru_R = np.log(rho_R_d * E_r / M_d)

    for cname, mo, tru, v_src, s_s, s_d in (
        ("gDNA", mo_g, tru_g, v_src_g, s_g_s, s_g_d),
        ("RNA ", mo_R, tru_R, v_src_R, s_R_s, s_R_d),
    ):
        err = mo - tru
        emp_var = float(np.var(err))          # sampling spread of the delivered mode
        emp_bias = float(np.mean(err))        # the systematic mismatch = b_c
        emp_mse = float(np.mean(err**2))      # the TRUE honest variance = var + bias²
        b_true = np.log(s_s) - np.log(s_d)    # analytic bias = log(s_c^src / s_c^dst,true)

        # SHIPPED law: precision = 1/v_src (only the source sampling; NO mismatch term)
        p_shipped = 1.0 / v_src
        overconf = emp_mse / v_src            # how confidently-wrong the shipped precision is

        # DERIVED (DL) law: fold the estimated squared bias from the conflict with the dst's OWN belief.
        if fg_d_own is None:                  # AMBIG: no own composition belief ⇒ v_own = ∞ ⇒ excess = 0
            excess = np.zeros(N)
            own_note = "AMBIG (v_own=∞)"
        else:
            s_d_own = fg_d_own if cname == "gDNA" else 1.0 - fg_d_own
            v_own = v_own_g if cname == "gDNA" else v_own_R
            o = np.log(s_d_own * (E_g if cname == "gDNA" else E_r) / M_d) + _lognormal(1.0, v_own, N) * 0 \
                + rng.normal(0.0, np.sqrt(v_own), N)   # own belief mode + its own sampling noise
            G = mo - o                        # observed share gap message-vs-own
            excess = np.maximum(0.0, G**2 - v_src - v_own)
            own_note = f"own f_g={fg_d_own:.3f} v_own={v_own:.3g}"
        sig2_DL = v_src + excess
        p_DL = float(np.mean(1.0 / sig2_DL))
        sig2_DL_mean = float(np.mean(sig2_DL))

        sig2_DL_med = float(np.median(sig2_DL))
        print(f"   [{cname}] b_true={b_true:+.3f}  emp_bias={emp_bias:+.3f}  "
              f"emp_var={emp_var:.4g}  true MSE={emp_mse:.4g}")
        print(f"          SHIPPED prec=1/v_src={p_shipped:9.4g}  → OVERCONFIDENCE {overconf:8.1f}×   "
              f"(claims {p_shipped:.3g}, earns {1.0/emp_mse:.3g})")
        print(f"          DERIVED(DL) σ²(median)={sig2_DL_med:9.4g}  σ²(mean)={sig2_DL_mean:9.4g}  "
              f"[{own_note}]   honest MSE={emp_mse:.4g}")

    # ── FUSION punchline: the RNA message (the corrupting channel) vs the dst's own belief in λ=logit f_g.
    # This is what actually decides f_g. p_msg (shipped or DL) competes with τ_own. Show the fused f_g. ──
    if tau_own_lam is not None:
        def _sigm(x):
            return 1.0 / (1.0 + np.exp(-x))

        lam_own = np.log(fg_d_own / (1.0 - fg_d_own))
        # the RNA message implies an f_g (1 − exp(log f_R)); use its mean mode as the message's λ pull.
        fR_msg = min(float(np.exp(np.mean(mo_R))), 0.999)
        lam_msg = np.log(max(1.0 - fR_msg, 1e-6) / max(fR_msg, 1e-6))
        p_ship = 1.0 / v_src_R
        o_R = np.log(1.0 - fg_d_own) + rng.normal(0.0, np.sqrt(v_own_R), N)  # own belief mode on log f_R
        p_dl = float(np.median(1.0 / (v_src_R + np.maximum(
            0.0, (mo_R - o_R) ** 2 - v_src_R - v_own_R))))
        fg_ship = _sigm((tau_own_lam * lam_own + p_ship * lam_msg) / (tau_own_lam + p_ship))
        fg_dl = _sigm((tau_own_lam * lam_own + p_dl * lam_msg) / (tau_own_lam + p_dl))
        print(f"   FUSE f_g (RNA msg vs own τ={tau_own_lam:g}): truth={fg_d_true:.3f}  own={fg_d_own:.3f}  "
              f"SHIPPED(p={p_ship:.1f})→{fg_ship:.3f}   DERIVED-DL(p={p_dl:.3g})→{fg_dl:.3f}")


def dl_recovers_bias(N):
    """E — the load-bearing claim: the DerSimonian–Laird excess recovers the TRUE squared bias b², so
    σ²_DL → the true MSE, when the dst own belief is unbiased for the truth. Swept over the mismatch size."""
    print("\n\n═══ E. DL excess recovers the true bias²  (b̂² = G² − v_src − v_own → b²) ═══")
    v_src, v_own = 0.04, 0.30
    print(f"   v_src={v_src}  v_own={v_own}   (own belief UNBIASED for the truth)")
    for b in (0.0, 0.5, 1.5, 2.6, 3.75):
        # message mode = truth + b + src-noise ; own = truth + own-noise ; both estimate 'truth'
        truth = 0.0
        msg = truth + b + rng.normal(0.0, np.sqrt(v_src), N)
        own = truth + rng.normal(0.0, np.sqrt(v_own), N)
        G = msg - own
        bhat2 = float(np.mean(np.maximum(0.0, G**2 - v_src - v_own)))
        sig2_dl = v_src + bhat2
        mse = float(np.mean((msg - truth) ** 2))
        print(f"   b={b:4.2f}  b²={b*b:7.3f}   b̂²(DL)={bhat2:7.3f}   "
              f"σ²_DL={sig2_dl:7.3f}   true MSE={mse:7.3f}   rel(σ²_DL vs MSE) {abs(sig2_dl-mse)/mse:6.1%}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# CROSS-CLIFF precision MC (agent #3)   draws={N:,}")
    print("# ε_c = mo_c − log f_c^dst,true = log(s_c^src) − log(s_c^dst,true)  (share mismatch; r,M_dst cancel)")
    print("# honest σ²_c = v_src,c + b_c² ;  b_c² estimated by DerSimonian–Laird excess (no magic number)")

    print("\n\n═══ A. ANCHOR: small mismatched RNA source across a ~425× cliff (node-1909-like) ═══")
    # dst = high-mass mostly-gDNA exon (oracle 0.985); its own weak self-solve says 0.953 (τ≈1.6).
    # src = tiny high-RNA boundary (64% RNA) with a real 47-count spliced junction (v_src,R≈1/26).
    transport_case(
        "anchor (mismatched: dst gDNA / src RNA)",
        N,
        fg_d_true=0.985, rho_tot_d=34.0,
        fg_d_own=0.953, v_own_g=0.05, v_own_R=0.56,   # weak own belief (τ≈1.6)
        fg_s=0.36, rho_tot_s=0.08,
        v_src_g=0.04, v_src_R=1.0 / 26.0,             # 47-count junction ⇒ strong RNA share (prec≈26)
        tau_own_lam=1.6,
    )

    print("\n\n═══ B. MATCHED same-scale transport (no cliff, same composition) — precision PRESERVED ═══")
    transport_case(
        "matched same-scale",
        N,
        fg_d_true=0.60, rho_tot_d=10.0,
        fg_d_own=0.60, v_own_g=0.10, v_own_R=0.10,
        fg_s=0.60, rho_tot_s=10.0,
        v_src_g=0.04, v_src_R=0.04,
        tau_own_lam=1.6,
    )

    print("\n\n═══ C. MATCHED-COMPOSITION BIG cliff (dense gDNA next to dense gDNA) — cliff alone must NOT damp ═══")
    transport_case(
        "matched-composition, r≈425",
        N,
        fg_d_true=0.97, rho_tot_d=34.0,
        fg_d_own=0.97, v_own_g=0.05, v_own_R=0.30,
        fg_s=0.965, rho_tot_s=0.08,                  # same composition, 425× less dense
        v_src_g=0.04, v_src_R=0.10,
        tau_own_lam=1.6,
    )

    print("\n\n═══ D. AMBIG destination (no own composition belief) across the cliff — message PROPAGATES ═══")
    transport_case(
        "AMBIG dst (v_own=∞)",
        N,
        fg_d_true=0.50, rho_tot_d=34.0,
        fg_d_own=None, v_own_g=0.0, v_own_R=0.0,      # None ⇒ excess forced to 0
        fg_s=0.36, rho_tot_s=0.08,
        v_src_g=0.04, v_src_R=1.0 / 26.0,
    )

    dl_recovers_bias(N)


if __name__ == "__main__":
    main()
