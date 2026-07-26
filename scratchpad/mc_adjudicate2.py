"""ADJUDICATION MC v2 — CORRECTED spliced model (composition-certain measured count, independent of f_g).

The bug in v1: rho_mu was tied to (1-fg), making the spliced carry the lambda latent. In truth rho_mu = S/E_spl
is a MEASURED pure-RNA density, independent of the gDNA composition (v_mu = 1/n_spl, NO composition term).
"""
import numpy as np

rng = np.random.default_rng(20260724)
N = 2_000_000


def beta_draws(mean, var, size):
    m = float(mean)
    v = min(float(var), m * (1 - m) * 0.98)
    c = m * (1 - m) / v - 1.0
    return rng.beta(m * c, (1 - m) * c, size=size)


def lognormal(var_log, size):
    s = np.sqrt(max(var_log, 0.0))
    return np.exp(rng.normal(-0.5 * s * s, s, size=size))


def combine_case(name, f_g, tau_lam, n_b, n_s, w_mu, E_g=150.0, E_r=200.0):
    var_fg = (f_g * (1 - f_g)) ** 2 / tau_lam
    fg = beta_draws(f_g, var_fg, N)
    Mb = lognormal(1.0 / n_b, N)          # node count noise (shared g & nu)
    Sfac = lognormal(1.0 / n_s, N)        # independent spliced count noise

    # operating-point densities defining shares
    rho_g0 = f_g / E_g
    rho_nu0 = (1 - f_g) / E_r
    # w_mu = rho_mu0 / (rho_nu0 + rho_mu0)  =>  rho_mu0 = w_mu/(1-w_mu) * rho_nu0
    rho_mu0 = (w_mu / (1 - w_mu)) * rho_nu0 if w_mu < 1 else 0.0

    rho_g = fg * Mb / E_g
    rho_nu = (1 - fg) * Mb / E_r
    rho_mu = rho_mu0 * Sfac               # MEASURED, composition-INDEPENDENT
    rho_R = rho_nu + rho_mu

    ag = rho_g * E_g
    aR = rho_R * E_r
    lam_dst = np.log(ag) - np.log(aR)     # true delivered composition constraint (matched-frame => M cancels)
    true_var = float(np.var(lam_dst))

    # ---- two-message per-arm predicted variances ----
    v_g = (1 - f_g) ** 2 / tau_lam
    # RNA-sum share-weighted: continue share (1-w_mu) carries composition f_g^2/tau; spliced share w_mu carries 1/n_s (+node 1/n_b via nu? nu shares Mb)
    # continue var (log rho_nu) = f_g^2/tau + 1/n_b ; spliced var (log rho_mu)=1/n_s ; sum share-weighted:
    w_nu = 1 - w_mu
    v_R = w_nu ** 2 * (f_g ** 2 / tau_lam + 1.0 / n_b) + w_mu ** 2 * (1.0 / n_s)
    p_g, p_R = 1.0 / v_g, 1.0 / v_R
    stated_var_two = 1.0 / (p_g * (1 - f_g) ** 2 + p_R * f_g ** 2)
    ratio_two = stated_var_two / true_var

    # ---- single-lambda via Var(log k) ----
    logk = np.log(rho_g / rho_R)
    emp_var_logk = float(np.var(logk))
    var_lam_src = 1.0 / tau_lam
    # Var(log k): d log k = d log rho_g - d log rho_R. rho_g: fg,Mb. rho_R: nu(fg,Mb)+mu(Sfac).
    # share-weighted: w_nu*(dlog nu) + w_mu*(dlog mu). Mb cancels between rho_g and w_nu part? partial.
    # Use the established T2 form: Var(log k)= w_mu^2*(1/n_b+1/n_s) + (1 - w_mu*f_g)^2 * var_lam_src  ... but that
    # assumed Mb shared across g and the WHOLE R. Here mu doesn't share Mb, so redo the count part:
    # d log k = [dlog rho_g] - [w_nu dlog rho_nu + w_mu dlog rho_mu]
    #   rho_g  : (1-fg)dlam_fg? ; work in lambda. dlog fg=(1-fg)dlam ; dlog(1-fg)=-fg dlam. Mb: +dlogMb in g and nu.
    #   dlog rho_g = (1-fg)dlam + dMb
    #   dlog rho_nu= -fg dlam + dMb ; dlog rho_mu = dS
    #   dlog rho_R = w_nu(-fg dlam + dMb) + w_mu dS
    #   dlog k = (1-fg)dlam + dMb - w_nu(-fg dlam+dMb) - w_mu dS
    #          = [(1-fg)+w_nu fg]dlam + (1-w_nu)dMb - w_mu dS
    #          = [1 - fg + w_nu fg]dlam + w_mu dMb - w_mu dS
    coef_lam = 1 - f_g + w_nu * f_g   # = 1 - w_mu*f_g
    pred_var_logk = coef_lam ** 2 * var_lam_src + w_mu ** 2 * (1.0 / n_b) + w_mu ** 2 * (1.0 / n_s)
    ratio_single = pred_var_logk / true_var

    print(f"\n[{name}] f_g={f_g} tau={tau_lam} n_b={n_b} n_s={n_s} w_mu={w_mu}")
    print(f"  true Var(lam)={true_var:.5g}  emp Var(log k)={emp_var_logk:.5g}  pred Var(log k)={pred_var_logk:.5g} ({abs(pred_var_logk-emp_var_logk)/emp_var_logk:.1%})")
    print(f"  Q1 two-message ratio  = {ratio_two:.3f}  ({1/ratio_two:.2f}x over-confident)")
    print(f"  Q2 single-lambda(full k) ratio = {ratio_single:.3f}  ({'OK <8%' if abs(ratio_single-1)<0.08 else 'OFF'})")


print("Q1/Q2  COMBINE over-confidence + single-lambda-via-Var(log k), CORRECT spliced model")
print("=" * 90)
combine_case("pure-comp symmetric", 0.5, 60.0, 10**9, 10**9, 0.0, 150.0, 150.0)
combine_case("spliced-free fg=0.4", 0.40, 60.0, 10**9, 10**9, 0.0)
combine_case("light-spliced   fg=0.4", 0.40, 60.0, 800, 40, 0.15)
combine_case("subst-spliced   fg=0.4", 0.40, 60.0, 800, 600, 0.50)
combine_case("heavy-spliced   fg=0.4", 0.40, 60.0, 800, 200, 0.90)
combine_case("heavy-spliced n_s=big", 0.40, 60.0, 800, 5000, 0.90)
