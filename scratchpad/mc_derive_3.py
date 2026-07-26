"""INDEPENDENT DERIVATION #3 — per-component message-variance model for the composition (unified) solver.

Framed from COUNT-ZERO-INFORMATION first principles. A fragment count carries no intrinsic gDNA/RNA
information; it only sets precision. So the derivation strategy is: write the message mode mo_c EXACTLY as
the code builds it (reframe by r=rho_tot(dst)/rho_tot(src), then divide by the destination's OWN observed
mass M_dst), then identify every common-mode factor (M_src, M_dst, the capture efficiency e) that cancels
in the composition, and keep only the noise that SURVIVES.

The message mode the psi solve consumes is, per component c on log f_c^dst:

    mo_c = log( rho_c^src * r * E_c^dst / M_dst ),    r = rho_tot(dst)/rho_tot(src)

and it enters psi as  -1/2 * p_c * (log f_c^dst - mo_c)^2, i.e. p_c = 1 / Var(mo_c).

KEY ALGEBRAIC CANCELLATIONS (derived, MC-verified below by ablation):

  rho_tot(dst) = M_dst * B_dst   (B_dst = f_g^dst/E_g + (1-f_g^dst)/E_r, the composition-weighted 1/E blend),
  so  r/M_dst = B_dst / rho_tot(src),  and

    mo_c = log rho_c^src + log B_dst - log rho_tot(src) + log E_c^dst.

  * M_dst CANCELS EXACTLY: it builds rho_tot(dst) (r's numerator) and the ÷M_dst normalizer. It contributes
    ZERO message variance. (count-zero-info, direction-INDEPENDENT).
  * rho_c^src = f_c^src * M_src/E_c^src and rho_tot(src) = rho_g+rho_nu+rho_mu shares M_src -> M_src cancels
    in the RATIO rho_c^src/rho_tot(src) = phi_c^src, EXCEPT against the spliced density rho_mu (independent
    count n_spl). So the source count 1/n survives ONLY share-weighted by w_mu (the spliced weight).
  * the capture efficiency e(dst)/e(src) is a fixed constant -> cancels identically (T1); its only stochastic
    residual is Var(log B_dst) (+ Var(log B_src)).

Thus  mo_g = log phi_g^src + log(B_dst E_g),  phi_g^src = rho_g/(rho_g+rho_nu+rho_mu).

GRAFT (boundary->exon, RNA is the SUM rho_R=rho_nu+rho_mu): mo built from phi. Share-weighted, convex.
PEEL  (exon->boundary, RNA-continue is the DIFFERENCE rho_nu = rho_R^x * r - rho_mu): u-weighted, r load-bearing.

METHOD: draw the physical primitives, compute mo_c EXACTLY as the code would (including building rho_tot(dst)
from an explicit Poisson M_dst so the cancellation is REAL, not assumed), measure Var(mo_c) empirically, and
compare against the derived closed form to <8% rel. Ablate each primitive's variance to confirm which sources
survive.  OMP_NUM_THREADS=1 python scratchpad/mc_derive_3.py
"""

from __future__ import annotations

import numpy as np

rng = np.random.default_rng(31337)
NDR = 1_500_000


def beta_draws(mean, var, size):
    m = float(mean)
    v = min(float(var), m * (1.0 - m) * 0.98)
    if v <= 0:
        return np.full(size, m)
    c = m * (1.0 - m) / v - 1.0
    return rng.beta(m * c, (1.0 - m) * c, size=size)


def lognm(mean, var_log, size):
    s = np.sqrt(max(float(var_log), 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<50} pred {pred:11.6g}  emp {emp:11.6g}  rel {rel:6.2%}")
    return rel < tol, rel


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# GRAFT — boundary -> exon.  Full physical model with EXPLICIT M_src and M_dst (Poisson), so the count-zero
# cancellations are emergent, not assumed.
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
def graft(name, f_g, var_fg, n, n_spl, E_g, E_r, E_spl, M0, S0, E_gd, E_rd, var_fg_dst,
          e_true=7.3, ablate=None):
    """boundary(src) -> exon(dst) GRAFT.  Returns empirical & predicted Var(mo_g), Var(mo_R)."""
    ab = ablate or {}
    # --- source belief + counts (the boundary) ---
    fg_s = beta_draws(f_g, 0.0 if ab.get("fg") else var_fg, NDR)          # source composition belief
    M_s = lognm(M0, 0.0 if ab.get("M") else 1.0 / n, NDR)                  # source unspliced count
    S_s = lognm(S0, 0.0 if ab.get("S") else 1.0 / n_spl, NDR)             # boundary spliced count (indep)
    rho_g = fg_s * M_s / E_g
    rho_nu = (1.0 - fg_s) * M_s / E_r
    rho_mu = S_s / E_spl
    rho_tot_s = rho_g + rho_nu + rho_mu                                    # WITH spliced (matched set)
    rho_R_s = rho_nu + rho_mu                                             # graft RNA = SUM

    # --- destination (the exon): explicit Poisson mass + its own composition belief for B_dst ---
    fg_d = beta_draws(f_g, 0.0 if ab.get("Bdst") else var_fg_dst, NDR)     # dst composition belief -> B_dst
    B_dst_g = fg_d / E_gd + (1.0 - fg_d) / E_rd
    # dst observed mass: Poisson around the true dst total.  M_dst noise MUST cancel -> include it explicitly.
    lam_dst = e_true * rho_tot_s.mean() * ((E_gd + E_rd) * 0.5)            # arbitrary true dst mean mass
    n_dst = lam_dst
    M_d = lognm(lam_dst, 0.0 if ab.get("Mdst") else 1.0 / n_dst, NDR)
    rho_tot_d = M_d * B_dst_g                                              # rho_tot(dst) = M_dst * B_dst

    r = rho_tot_d / rho_tot_s                                             # enrichment reframe (carries M_dst)
    # message modes exactly as the code: reframe every component by r, ÷ M_dst.
    mo_g = np.log(rho_g * r * E_gd / M_d)                                  # -> log f_g^dst
    mo_R = np.log(rho_R_s * r * E_rd / M_d)                                # -> log f_R^dst

    emp_g = float(np.var(mo_g))
    emp_R = float(np.var(mo_R))

    # ---- PREDICTED (closed form, count-zero-info) ----
    a, b, c = f_g * M0 / E_g, (1.0 - f_g) * M0 / E_r, S0 / E_spl
    tot = a + b + c
    w_g, w_nu, w_mu = a / tot, b / tot, c / tot
    B_d0 = f_g / E_gd + (1.0 - f_g) / E_rd
    var_logB = ((1.0 / E_gd - 1.0 / E_rd) / B_d0) ** 2 * var_fg_dst        # Var(log B_dst), delta method
    # gDNA:  d log phi_g = A_g df_g + w_mu dlogM - w_mu dlogS ;  + d log B_dst
    A_g = (w_nu + w_mu) / f_g + w_nu / (1.0 - f_g)
    pred_g = A_g**2 * var_fg + w_mu**2 * (1.0 / n + 1.0 / n_spl) + var_logB
    # RNA-total:  phi_R = (b+c)/tot.  d log phi_R = coefficients below (delta method, exact)
    #   d log(b+c) = om_nu dlogb + om_mu dlogc ;  d log tot = w_g dloga + w_nu dlogb + w_mu dlogc
    om_nu, om_mu = b / (b + c), c / (b + c)
    #   dloga = df/f + dlogM ; dlogb = -df/(1-f)+dlogM ; dlogc = dlogS
    cf = (om_nu - w_nu) * (-1.0 / (1.0 - f_g)) + (-w_g) * (1.0 / f_g)      # coeff of df_g
    cM = (om_nu - w_nu) + (-w_g)                                          # coeff of dlogM
    cS = (om_mu - w_mu)                                                   # coeff of dlogS
    pred_R = cf**2 * var_fg + cM**2 * (1.0 / n) + cS**2 * (1.0 / n_spl) + var_logB

    print(f"\n GRAFT [{name}]  f_g={f_g}  w_mu={w_mu:.3f}  A_g={A_g:.3f}")
    report("Var(mo_g)  gDNA graft", pred_g, emp_g)
    report("Var(mo_R)  RNA-sum graft", pred_R, emp_R)
    return dict(emp_g=emp_g, pred_g=pred_g, emp_R=emp_R, pred_R=pred_R,
                w_mu=w_mu, A_g=A_g, var_logB=var_logB)


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# PEEL — exon -> boundary.  RNA-continue is a DIFFERENCE; enrichment is load-bearing; ÷M_b adds (u-1)^2/n_b.
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
def peel(name, f_g_x, var_fg_x, n_x, E_g, E_r, M0_x, target_u, E_spl, n_b, n_spl,
         var_fg_b, M0_b=3000.0, ablate=None):
    """exon(src) -> boundary(dst) PEEL. rho_nu = rho_R^x * r - rho_mu ; enrichment r=rho_tot_b/rho_tot_x is
    built from the observed M_b (so its cancellation with ÷M_b is REAL). S0_b is chosen to hit target_u."""
    ab = ablate or {}
    # operating point: phi_R^x, B_b, T0, pick rho_mu0 to hit target u.
    a0, b0 = f_g_x / E_g, (1.0 - f_g_x) / E_r
    phiR0, wg0 = b0 / (a0 + b0), a0 / (a0 + b0)
    B_b0 = f_g_x / E_g + (1.0 - f_g_x) / E_r
    T0 = phiR0 * M0_b * B_b0                                              # T = phi_R^x * M_b * B_b
    rho_mu0 = T0 * (1.0 - 1.0 / target_u)
    S0_b = rho_mu0 * E_spl
    # source exon
    fg_x = beta_draws(f_g_x, 0.0 if ab.get("fgx") else var_fg_x, NDR)
    M_x = lognm(M0_x, 0.0 if ab.get("Mx") else 1.0 / n_x, NDR)
    rho_R_x = (1.0 - fg_x) * M_x / E_r
    rho_tot_x = fg_x * M_x / E_g + rho_R_x
    # destination boundary
    S_b = lognm(S0_b, 0.0 if ab.get("S") else 1.0 / n_spl, NDR)
    rho_mu = S_b / E_spl
    fg_b = beta_draws(f_g_x, 0.0 if ab.get("Bb") else var_fg_b, NDR)
    B_b = fg_b / E_g + (1.0 - fg_b) / E_r
    M_b = lognm(M0_b, 0.0 if ab.get("Mb") else 1.0 / n_b, NDR)
    rho_tot_b = M_b * B_b
    r = rho_tot_b / rho_tot_x                                             # carries M_b (its ÷M_b twin cancels)
    T = rho_R_x * r
    rho_nu = T - rho_mu
    keep = rho_nu > 0
    frac = keep.mean()
    mo_nu = np.log(np.maximum(rho_nu, 1e-300) * E_r / M_b)                 # -> log f_nu^dst (boundary)
    emp = float(np.var(mo_nu[keep]))

    # ---- predicted (u-weighted difference, incl ÷M_b conversion) ----
    u = target_u
    # Var(log phi_R^x): M_x cancels (exon has no independent spliced count) -> pure composition.
    var_logPhiR = (wg0 / (f_g_x * (1.0 - f_g_x))) ** 2 * var_fg_x
    # transfer Var(log B_b) (dst composition) — B_x CANCELS (rho_tot_x is the reframe denom of phi_R^x).
    var_logBb = ((1.0 / E_g - 1.0 / E_r) / B_b0) ** 2 * var_fg_b
    pred = u**2 * (var_logPhiR + var_logBb) + (u - 1.0) ** 2 * (1.0 / n_b + 1.0 / n_spl)
    print(f"\n PEEL [{name}]  u={u:.2f} (emp {np.mean(T[keep])/np.mean(rho_nu[keep]):.2f})  kept {frac:.1%}")
    report("Var(mo_nu)  RNA-continue peel", pred, emp, tol=0.12)
    return dict(emp=emp, pred=pred, u=u, var_logBb=var_logBb, var_logPhiR=var_logPhiR)


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# COMBINE — psi gets a gDNA message on log f_g AND an RNA message on log f_R FROM ONE SOURCE, as INDEPENDENT
# Gaussians on the same 1-DOF composition (lambda).  Does that reproduce the exact joint dst-composition var?
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────
def combine(name, f_g, var_lambda, n, n_spl, E_g, E_r, E_spl, M0, S0, E_gd, E_rd):
    """One source belief -> two messages (mo_g on log f_g, mo_R on log f_R).  Compare the psi posterior var on
    lambda from treating them INDEPENDENT vs the true single-scalar (joint) transport variance."""
    # source belief is ONE latent lambda = logit f_g ; draw it, no separate counts here (isolate the
    # composition-correlation effect the combine mishandles).  Include the graft structure (spliced sum).
    var_fg = var_lambda * (f_g * (1.0 - f_g)) ** 2
    fg_s = beta_draws(f_g, var_fg, NDR)
    M_s = lognm(M0, 1.0 / n, NDR)
    S_s = lognm(S0, 1.0 / n_spl, NDR)
    rho_g = fg_s * M_s / E_g
    rho_nu = (1.0 - fg_s) * M_s / E_r
    rho_mu = S_s / E_spl
    rho_tot = rho_g + rho_nu + rho_mu
    rho_R = rho_nu + rho_mu
    # transport to a dst frame (fixed B_dst -> isolate the source-correlation, no transfer noise)
    B_d = f_g / E_gd + (1.0 - f_g) / E_rd
    mo_g = np.log(rho_g / rho_tot * B_d * E_gd)
    mo_R = np.log(rho_R / rho_tot * B_d * E_rd)
    Vg = float(np.var(mo_g))
    VR = float(np.var(mo_R))
    corr = float(np.corrcoef(mo_g, mo_R)[0, 1])

    # The dst composition is 1-DOF: lambda_d = logit f_g^dst.  psi consumes the TWO messages INDEPENDENTLY:
    #   psi(lambda) = -1/2 p_g (log f_g(lam) - mo_g)^2  - 1/2 p_R (log f_R(lam) - mo_R)^2 ,  p_g=1/Vg, p_R=1/VR.
    # Linearize at the operating point: d log f_g/dlam = 1-f_g = J_g ; d log f_R/dlam = -f_g -> |J_R|=f_g.
    # A single message alone pins lambda at  lam_g = lam0 + (mo_g - logf_g0)/J_g  (var Vg/J_g^2), and the RNA
    # message alone at  lam_R = lam0 + (logf_R0 - mo_R)/f_g  (var VR/f_g^2).  The independent-combine argmax is
    # the precision-weighted mean, and psi REPORTS posterior var = 1/(p_g J_g^2 + p_R f_g^2).
    Jg, Jf = (1.0 - f_g), f_g
    lam0 = np.log(f_g / (1.0 - f_g))
    logfg0, logfr0 = np.log(f_g), np.log(1.0 - f_g)
    lam_g = lam0 + (mo_g - logfg0) / Jg          # lambda estimate from the gDNA message
    lam_R = lam0 + (logfr0 - mo_R) / Jf          # lambda estimate from the RNA message (opposite sign)
    p_g, p_R = 1.0 / Vg, 1.0 / VR
    wg = p_g * Jg**2 / (p_g * Jg**2 + p_R * Jf**2)
    wR = 1.0 - wg
    lam_hat = wg * lam_g + wR * lam_R            # the independent-combine point estimate
    var_reported = 1.0 / (p_g * Jg**2 + p_R * Jf**2)   # what psi CLAIMS (two independent messages)
    var_actual = float(np.var(lam_hat))               # the naive estimator's TRUE spread
    var_single_g = float(np.var(lam_g))               # one honest scalar transport (gDNA arm)
    var_single_R = float(np.var(lam_R))               # one honest scalar transport (RNA arm)
    # The physically-correct JOINT-TRUTH is the estimator's actual delivered spread (var_actual): what the
    # independent-combine point estimate ACTUALLY varies by, given the true joint distribution of the two
    # messages. The GLS 1/(1^T Sigma_eps^-1 1) is DEGENERATE here (shown below): the two lambda-observations
    # are ~collinear (corr ~-0.99) because both are deterministic functions of ONE source latent, so GLS
    # spuriously "cancels" a difference-mode that carries no real information -> an artifactual tiny variance.
    Se = np.cov(lam_g, lam_R)
    ones = np.array([1.0, 1.0])
    var_gls = 1.0 / float(ones @ np.linalg.solve(Se, ones))
    ratio = var_reported / var_actual                 # reported vs the honest delivered spread

    print(f"\n COMBINE [{name}]  f_g={f_g}  corr(mo_g,mo_R)={corr:+.3f}")
    print(f"     single-msg gDNA={var_single_g:.5g}  single-msg RNA={var_single_R:.5g}")
    print(f"     two-msg REPORTED={var_reported:.5g}  joint-truth(delivered spread)={var_actual:.5g}"
          f"  [GLS degenerate={var_gls:.2g}]")
    print(f"     ratio reported/joint-truth = {ratio:.3f}   (1.0 = calibrated; <1 = OVERCONFIDENT)")
    return dict(corr=corr, var_reported=var_reported, var_actual=var_actual, var_gls=var_gls,
                var_single_g=var_single_g, var_single_R=var_single_R, ratio=ratio)


def main():
    print(f"# MC derivation #3 — per-component message variance   draws={NDR:,}")
    # E_g != E_r everywhere so the transfer term Var(log B) is real.
    # (A) substantial spliced flux -> counting does NOT fully cancel (w_mu ~ 0.4)
    graft("substantial-spliced", f_g=0.40, var_fg=0.004, n=800, n_spl=600,
          E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1500.0,
          E_gd=120.0, E_rd=260.0, var_fg_dst=0.004)
    # (B) NO spliced -> count must CANCEL (w_mu=0); Var(mo_g) collapses to (A_g)^2 Var(f_g)+Var(logB)
    graft("no-spliced", f_g=0.40, var_fg=0.004, n=800, n_spl=10**9,
          E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1e-9,
          E_gd=120.0, E_rd=260.0, var_fg_dst=0.004)
    # (C) ANCHOR limit: near-pure gDNA source, structural (Var(f_g)->0), no spliced.  k=rho_g/rho_R singular;
    #     the density-mode mo_g stays FINITE = Var(log B_dst).
    graft("anchor f_g->1", f_g=0.985, var_fg=2e-5, n=1200, n_spl=10**9,
          E_g=110.0, E_r=200.0, E_spl=100.0, M0=2000.0, S0=1e-9,
          E_gd=120.0, E_rd=260.0, var_fg_dst=1e-5)

    # ABLATIONS — confirm which sources survive on the graft (A).
    print("\n --- GRAFT ablations (scenario A), empirical Var(mo_g) ---")
    base = graft("A/ablate-none", 0.40, 0.004, 800, 600, 110.0, 200.0, 100.0, 900.0, 1500.0,
                 120.0, 260.0, 0.004)["emp_g"]
    for key, lab in [("Mdst", "M_dst (must CANCEL -> ~0 change)"),
                     ("M", "M_src (survives only w_mu-weighted)"),
                     ("S", "S_spl (survives w_mu-weighted)"),
                     ("fg", "f_g^src composition (dominant)"),
                     ("Bdst", "B_dst transfer (small)")]:
        v = graft(f"A/ablate-{key}", 0.40, 0.004, 800, 600, 110.0, 200.0, 100.0, 900.0, 1500.0,
                  120.0, 260.0, 0.004, ablate={key: True})["emp_g"]
        print(f"     ablate {lab:<42}  Var {v:10.6g}  (Δ from base {100*(v-base)/base:+6.1f}%)")

    # PEEL — enrichment load-bearing; u-weighted difference.
    peel("comfortable u~2", f_g_x=0.05, var_fg_x=0.0008, n_x=5000, E_g=110.0, E_r=200.0,
         M0_x=4000.0, target_u=2.0, E_spl=100.0, n_b=800, n_spl=1500, var_fg_b=0.004)
    peel("tight u~6", f_g_x=0.05, var_fg_x=0.0008, n_x=5000, E_g=110.0, E_r=200.0,
         M0_x=4000.0, target_u=6.0, E_spl=100.0, n_b=800, n_spl=1500, var_fg_b=0.004)
    print("\n --- PEEL ablation: enrichment B_b (u^2-amplified) is load-bearing ---")
    p_on = peel("peel transfer ON (u~4)", 0.05, 0.0008, 5000, 110.0, 200.0, 4000.0, 4.0, 100.0,
                800, 1500, 0.004)["emp"]
    p_off = peel("peel transfer OFF (B_b ablated)", 0.05, 0.0008, 5000, 110.0, 200.0, 4000.0, 4.0,
                 100.0, 800, 1500, 0.004, ablate={"Bb": True})["emp"]
    print(f"     Var(mo_nu): transfer ON {p_on:.5g}  vs transfer OFF {p_off:.5g}"
          f"   (B_b transfer carries {100*(p_on-p_off)/p_on:.0f}% )")

    # COMBINE — the independence question.
    combine("balanced f_g=0.4", f_g=0.40, var_lambda=0.05, n=800, n_spl=600,
            E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1500.0, E_gd=120.0, E_rd=260.0)
    combine("gDNA-lean f_g=0.15", f_g=0.15, var_lambda=0.05, n=800, n_spl=600,
            E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1500.0, E_gd=120.0, E_rd=260.0)
    combine("gDNA-rich f_g=0.7", f_g=0.70, var_lambda=0.05, n=800, n_spl=600,
            E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1500.0, E_gd=120.0, E_rd=260.0)
    combine("no-spliced f_g=0.4", f_g=0.40, var_lambda=0.05, n=800, n_spl=10**9,
            E_g=110.0, E_r=200.0, E_spl=100.0, M0=900.0, S0=1e-9, E_gd=120.0, E_rd=260.0)


if __name__ == "__main__":
    main()
