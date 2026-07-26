"""ADVERSARIAL verifier #2 — attack the COMBINE resolution (single-λ message, resolution 1).

Mandate: pressure-test whether the "single-λ message" fix actually
  (A) reproduces the joint destination-composition posterior on the matched graft,
  (B) stays FINITE at the pure-gDNA anchor (f_g→1), and
  (C) handles a PARTIAL (one-component) message correctly,
and whether it hides a singularity or a magic number (e.g. a regime switch between the
λ-message and the ÷M density mode, or a latent spliced double-count).

Everything is a pure-numpy delta/exact MC. Composition is a SINGLE latent λ=logit f_g
(rank-1) so the gDNA/RNA arms are deterministic multiples — the whole point of M6.

    OMP_NUM_THREADS=1 python scratchpad/mc_adv_2.py
"""

from __future__ import annotations

import numpy as np

rng = np.random.default_rng(20260724)


def _ln(mean, var_log, size):
    s = np.sqrt(max(var_log, 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def _draw_source(N, *, f_g0, tau_lam, n, n_spl, E_g, E_r, E_spl, M0, S0, sigma2_transfer):
    """Draw a source belief with ONE composition latent λ (rank-1), shared count M, indep spliced, indep r."""
    lam0 = np.log(f_g0 / (1.0 - f_g0))
    dlam = rng.normal(0.0, np.sqrt(max(tau_lam, 1e-300)) ** -1 if tau_lam > 0 else 0.0, N)
    lam = lam0 + dlam
    fg = 1.0 / (1.0 + np.exp(-lam))                       # rank-1 composition
    M = _ln(M0, 1.0 / n, N)                               # shared unspliced count
    S = _ln(S0, 1.0 / n_spl, N) if S0 > 0 else np.zeros(N)  # independent spliced count
    r = _ln(1.0, sigma2_transfer, N) if sigma2_transfer > 0 else np.ones(N)
    rho_g = fg * M / E_g
    rho_nu = (1.0 - fg) * M / E_r
    rho_mu = S / E_spl
    rho_R = rho_nu + rho_mu
    return fg, rho_g, rho_R, rho_nu, rho_mu, r


def joint_fg_dst(rho_g, rho_R, r, E_g_dst, E_r_dst):
    """The TRUTH ψ should recover: the matched-graft destination fraction (r cancels in the ratio)."""
    cg, cR = rho_g * r, rho_R * r
    fgx = cg * E_g_dst / (cg * E_g_dst + cR * E_r_dst)
    return fgx


def two_msg_var(fgd0, v_g_msg, v_R_msg):
    """ψ combining TWO independent Gaussians on (log f_g, log f_R): the code's current behaviour.
    λ-Fisher = (1−f)²/v_g + f²/v_R ; report Var(log f_g^dst) = (1−f)²/λFisher."""
    lam_fisher = (1.0 - fgd0) ** 2 / v_g_msg + fgd0 ** 2 / v_R_msg
    return (1.0 - fgd0) ** 2 / lam_fisher


SEP = "─" * 100


def scenario(tag, *, f_g0, tau_lam, n, n_spl, E_g, E_r, E_spl, M0, S0, sigma2_transfer,
             Egx, Erx, N=1_500_000):
    fg, rho_g, rho_R, rho_nu, rho_mu, r = _draw_source(
        N, f_g0=f_g0, tau_lam=tau_lam, n=n, n_spl=n_spl, E_g=E_g, E_r=E_r, E_spl=E_spl,
        M0=M0, S0=S0, sigma2_transfer=sigma2_transfer)

    fgx = joint_fg_dst(rho_g, rho_R, r, Egx, Erx)
    var_joint = float(np.var(np.log(np.maximum(fgx, 1e-300))))
    fgd0 = float(np.mean(fgx))

    # empirical joint composition variance in the k-parameterization (r cancels in k too)
    logk = np.log(np.maximum(rho_g, 1e-300)) - np.log(np.maximum(rho_R, 1e-300))
    var_logk_full = float(np.var(logk))          # INCLUDES spliced noise (via rho_R)
    # composition-only Var(log k): remove the spliced by using rho_nu instead of rho_R
    logk_comp = np.log(np.maximum(rho_g, 1e-300)) - np.log(np.maximum(rho_nu, 1e-300))
    var_logk_comp = float(np.var(logk_comp))

    # operating-point per-component message variances (what the code feeds ψ)
    var_log_fg = (1.0 - f_g0) ** 2 / tau_lam if tau_lam > 0 else 0.0
    var_log_fR = f_g0 ** 2 / tau_lam if tau_lam > 0 else 0.0
    rg0, rn0, rm0 = f_g0 * M0 / E_g, (1.0 - f_g0) * M0 / E_r, (S0 / E_spl if S0 > 0 else 0.0)
    rR0 = rn0 + rm0
    w_nu = rn0 / rR0 if rR0 > 0 else 0.0
    w_mu = rm0 / rR0 if rR0 > 0 else 0.0
    v_g_msg = var_log_fg + 1.0 / n + sigma2_transfer
    v_R_msg = w_nu ** 2 * (var_log_fR + 1.0 / n) + w_mu ** 2 / n_spl + sigma2_transfer

    # ── the candidate combine schemes, each → predicted Var(log f_g^dst) ──
    va = two_msg_var(fgd0, v_g_msg, v_R_msg)                     # (a) code today: 2 indep msgs
    vd_g_alone = v_g_msg                                        # (d) single gDNA density msg ALONE
    v_c_full = (1.0 - fgd0) ** 2 * var_logk_full               # (c1) single-λ, FULL joint Var(log k)
    v_c_comp = (1.0 - fgd0) ** 2 * var_logk_comp               # (c-comp) single-λ, comp-only Var(log k)
    # (c2) resolution-1 LITERAL: single-λ(comp-only) + a SEPARATE spliced message on log f_R.
    #      spliced msg λ-info = f² / (w_mu²/n_spl). Does adding it to comp-only λ reproduce joint,
    #      or does using FULL Var(log k) + separate spliced DOUBLE-count?
    lamfish_comp = 1.0 / var_logk_comp
    lamfish_spl = fgd0 ** 2 / (w_mu ** 2 / n_spl) if w_mu > 0 else 0.0
    v_res1_correct = (1.0 - fgd0) ** 2 / (lamfish_comp + lamfish_spl)   # comp-only + spliced
    lamfish_full = 1.0 / var_logk_full
    v_res1_double = (1.0 - fgd0) ** 2 / (lamfish_full + lamfish_spl)    # FULL(k) + spliced = DOUBLE

    def row(name, pred):
        rel = abs(pred - var_joint) / max(var_joint, 1e-300)
        ratio = pred / max(var_joint, 1e-300)
        flag = "OK " if 0.90 <= ratio <= 1.11 else ("OVER" if ratio < 1 else "UNDR")
        print(f"    {flag}  {name:<48} pred {pred:11.6g}  ratio {ratio:6.3f}  rel {rel:6.1%}")

    print(f"\n{SEP}\n{tag}")
    print(f"    f_g^src={f_g0}  τ_λ={tau_lam}  n={n}  n_spl={n_spl}  w_μ={w_mu:.3f}  "
          f"f_g^dst={fgd0:.4f}  σ²t={sigma2_transfer}")
    print(f"    JOINT TRUTH Var(log f_g^dst) = {var_joint:.6g}   "
          f"[Var(log k) full={var_logk_full:.5g}  comp-only={var_logk_comp:.5g}]")
    row("(c1) single-λ  FULL Var(log k)      [correct?]", v_c_full)
    row("(c2a) res-1: comp-Var(log k) + spliced msg", v_res1_correct)
    row("(c2b) res-1 MISREAD: FULL Var(k)+spliced", v_res1_double)
    row("(a)  TWO indep msgs (code today)", va)
    row("(d)  single gDNA density msg ALONE", vd_g_alone)
    return dict(var_joint=var_joint, va=va, vd=vd_g_alone, vc=v_c_full, fgd0=fgd0,
                var_logk_full=var_logk_full, w_mu=w_mu)


def anchor_sweep(N=1_500_000):
    """(B) the PURE-gDNA ANCHOR attack: sweep f_g^src → 1 and watch the single-λ Var(log k) form.
    Does τ_T = 1/Var(log k) → 0 (message dies), and does the destination Jacobian (1−f)² rescue the
    joint? At exactly f_g=1 there is NO k (k=∞) — show the density mode is the only finite mechanism."""
    print(f"\n{SEP}\nANCHOR ATTACK (B): single-λ Var(log k) vs ÷M density mode as f_g^src → 1")
    print(f"    {'f_g^src':>8}  {'Var(logk)':>11}  {'τ_T=1/Vk':>10}  {'joint Vfgd':>11}  "
          f"{'(1-f)²·Vk':>11}  {'dens-mode':>10}  {'match?':>7}")
    E_g, E_r, E_spl = 110.0, 200.0, 100.0
    n, n_spl, M0, S0 = 800, 10 ** 9, 900.0, 0.0   # no spliced: RNA is pure composition
    tau_lam = 60.0
    for f_g0 in (0.4, 0.7, 0.9, 0.97, 0.995, 0.9995):
        fg, rho_g, rho_R, rho_nu, rho_mu, r = _draw_source(
            N, f_g0=f_g0, tau_lam=tau_lam, n=n, n_spl=n_spl, E_g=E_g, E_r=E_r, E_spl=E_spl,
            M0=M0, S0=S0, sigma2_transfer=0.0)
        fgx = joint_fg_dst(rho_g, rho_R, r, 380.0, 290.0)
        var_joint = float(np.var(np.log(np.maximum(fgx, 1e-300))))
        fgd0 = float(np.mean(fgx))
        logk = np.log(np.maximum(rho_g, 1e-300)) - np.log(np.maximum(rho_R, 1e-300))
        var_logk = float(np.var(logk))
        tau_T = 1.0 / var_logk
        v_klam = (1.0 - fgd0) ** 2 * var_logk       # single-λ delivered Var(log f_g^dst)
        # ÷M gDNA density mode: Var(log f_g^dst)=Var(log f_g^src)+1/n  (M4 anchor law)
        v_dens = (1.0 - f_g0) ** 2 / tau_lam + 1.0 / n
        ratio = v_klam / max(var_joint, 1e-300)
        match = "OK" if 0.9 <= ratio <= 1.11 else f"{ratio:.2f}x"
        print(f"    {f_g0:8.4f}  {var_logk:11.5g}  {tau_T:10.5g}  {var_joint:11.6g}  "
              f"{v_klam:11.6g}  {v_dens:10.5g}  {match:>7}")
    print("    NOTE the density mode (last col) is what ψ SHOULD receive on the FRACTION; the single-λ")
    print("    Var(log k) form tracks the JOINT (they agree — the (1−f)² Jacobian rescues the blow-up)")
    print("    but at exactly f_g=1 Var(log k)=∞ AND (1−f)²=0 → 0·∞: undefined. The two mechanisms")
    print("    deliver DIFFERENT numbers (joint→0 fraction-var vs density-mode→1/n): the handoff is real.")


def partial_msg():
    """(C) PARTIAL message: source is structurally pure gDNA (ρ_R ≡ 0). There is NO k, NO λ-message.
    Only the ÷M gDNA density mode is defined and finite. Confirm single-λ cannot represent it."""
    print(f"\n{SEP}\nPARTIAL-MESSAGE ATTACK (C): structurally pure-gDNA seam (ρ_R ≡ 0)")
    N = 1_500_000
    n, E_g = 800, 110.0
    M = _ln(900.0, 1.0 / n, N)
    rho_g = M / E_g            # f_g ≡ 1 structurally
    rho_R = np.zeros(N)
    logk = np.log(rho_g) - np.log(np.maximum(rho_R, 1e-300))
    print(f"    log k mean = {np.mean(logk):.3g} (→ +inf as ρ_R→0); Var(log k) = "
          f"{np.var(logk):.3g}  →  single-λ Var(log k) is UNDEFINED / degenerate")
    v_dens = float(np.var(np.log(rho_g)))   # the density mode: pure count 1/n
    print(f"    ÷M density-mode Var(log ρ_g) = {v_dens:.5g}  (pred 1/n = {1.0/n:.5g})  FINITE ✓")
    print("    ⇒ resolution 1 MUST fall back to the ÷M density mode here; 'single-λ message' is a")
    print("      separate mechanism from the anchor/partial mechanism — it is a TWO-regime scheme.")


def main():
    print("# ADVERSARIAL MC — attacking the single-λ (resolution 1) combine fix\n")

    print("=" * 100)
    print("(A) MATCHED GRAFT — does single-λ reproduce the joint, and where do the schemes land?")
    # Case A: substantial spliced flux (w_μ~0.85) — BOTH composition & spliced live
    scenario("[A1] substantial spliced (w_μ≈0.85), both components carry composition",
             f_g0=0.40, tau_lam=60.0, n=800, n_spl=600, E_g=110.0, E_r=200.0, E_spl=100.0,
             M0=900.0, S0=1500.0, sigma2_transfer=0.0, Egx=380.0, Erx=290.0)
    # Case B: NO spliced (w_μ=0) — RNA is pure composition (the cleanest double-count)
    scenario("[A2] NO spliced (w_μ=0), RNA is pure composition",
             f_g0=0.40, tau_lam=60.0, n=800, n_spl=10 ** 9, E_g=110.0, E_r=200.0, E_spl=100.0,
             M0=900.0, S0=1e-9, sigma2_transfer=0.0, Egx=380.0, Erx=290.0)
    # Case C: high f_g (0.7) no spliced — closer to anchor
    scenario("[A3] f_g=0.70, no spliced",
             f_g0=0.70, tau_lam=60.0, n=800, n_spl=10 ** 9, E_g=110.0, E_r=200.0, E_spl=100.0,
             M0=900.0, S0=1e-9, sigma2_transfer=0.0, Egx=380.0, Erx=290.0)
    # Case D: spliced-DOMINATED (w_μ~0.85) — does FULL Var(log k) already contain the spliced?
    scenario("[A4] spliced-dominated (w_μ≈0.85) — spliced double-count probe",
             f_g0=0.25, tau_lam=40.0, n=500, n_spl=300, E_g=110.0, E_r=200.0, E_spl=100.0,
             M0=900.0, S0=2500.0, sigma2_transfer=0.0, Egx=380.0, Erx=290.0)

    anchor_sweep()
    partial_msg()


if __name__ == "__main__":
    main()
