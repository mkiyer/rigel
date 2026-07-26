"""ADJUDICATION MC — resolve the four disputed points, not re-derive everything.

Q1  Combine over-confidence: measure ratio (two-message stated var / joint-truth var) as a function of the
    independent-count content, to reconcile 2x (der1/2), 1.06-1.96x (adv2), 1.2-3.6x (adv1).
Q2  Does single-lambda-via-Var(log k) reproduce the joint? Across spliced regimes incl heavy-spliced.
Q3  adv-2's spliced fork: (a) full Var(log k) alone, (b) full Var(log k)+separate spliced msg,
    (c) comp-only Var(log k)+separate spliced. Which reproduces the joint?
Q4  Anchor: is Var(log k) singular at f_g->1 while the divM gDNA density mode is finite?
Q5  adv-1's M1 near-wall claim: Var(log f_g)=var_fg/f_g^2 for a small+uncertain minority component.
"""
import numpy as np

rng = np.random.default_rng(20260724)
N = 2_000_000


def beta_draws(mean, var, size):
    m = float(mean)
    v = min(float(var), m * (1 - m) * 0.98)
    c = m * (1 - m) / v - 1.0
    return rng.beta(m * c, (1 - m) * c, size=size)


def lognormal(mean, var_log, size):
    s = np.sqrt(max(var_log, 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def combine_case(name, f_g, tau_lam, n_b, n_s, w_mu, E_g=150.0, E_r=200.0):
    """Faithful ψ two-message combine vs joint truth, on a matched-set graft.

    Source composition dof lambda drives BOTH arms. gDNA arm: log f_g. RNA arm: log f_R = log(1-f_g), where
    f_R = f_nu + f_mu (imputed continue + measured spliced). w_mu = f_mu/f_R is the spliced share.
    """
    var_fg = (f_g * (1 - f_g)) ** 2 / tau_lam            # foundation: Var(f_g) from tau_lam
    fg = beta_draws(f_g, var_fg, N)
    Mb = lognormal(1.0, 1.0 / n_b, N)                    # node count noise (shared by g and nu)
    # RNA total = continue + spliced. continue ~ (1-fg) (composition), spliced ~ measured w_mu share w/ own count
    Sfac = lognormal(1.0, 1.0 / n_s, N)                  # independent spliced count noise
    rho_g = fg * Mb / E_g
    rho_nu = (1 - fg) * (1 - w_mu) * Mb / E_r
    rho_mu = (1 - fg) * w_mu * Sfac / E_r                # spliced carries the SAME composition (1-fg) but own count
    rho_R = rho_nu + rho_mu
    # destination log-fractions delivered (matched frame, M cancels in the fraction => use content ratio)
    # f_g^dst propto rho_g*E_g ; f_R^dst propto rho_R*E_r ; normalize
    ag = rho_g * E_g
    aR = rho_R * E_r
    logfg = np.log(ag / (ag + aR))
    logfR = np.log(aR / (ag + aR))
    lam_dst = logfg - logfR                              # the true delivered composition constraint

    true_var_lam = float(np.var(lam_dst))

    # ---- predicted per-arm message variances (the derivation's laws, at the operating point) ----
    v_g = (1 - f_g) ** 2 / tau_lam                       # Var(log f_g) transport seed (matched graft, M cancels)
    # RNA-sum share-weighted: w_nu^2 * f_g^2/tau + w_mu^2*(1/n_b+1/n_s)  [comp part uses f_g^2/tau for continue]
    v_R = (1 - w_mu) ** 2 * (f_g ** 2 / tau_lam) + w_mu ** 2 * (1.0 / n_b + 1.0 / n_s)
    p_g = 1.0 / v_g
    p_R = 1.0 / v_R
    # two independent messages: lambda-Fisher = p_g*(dlogfg/dlam)^2 + p_R*(dlogfR/dlam)^2 = p_g*(1-fg)^2+p_R*fg^2
    stated_prec_two = p_g * (1 - f_g) ** 2 + p_R * f_g ** 2
    stated_var_two = 1.0 / stated_prec_two
    ratio_two = stated_var_two / true_var_lam

    # single-lambda via Var(log k): k = rho_g/rho_R, lam = log k (matched frame E_g=E_r; here E differ so lam=logfg-logfR)
    logk = np.log(rho_g / rho_R)
    var_logk_emp = float(np.var(logk))
    # predicted Var(log k) = w_mu^2*(1/n_b+1/n_s) + (1 - w_mu*f_g)^2 * var_lam_src, var_lam_src=var_fg/(fg(1-fg))^2=1/tau
    var_lam_src = var_fg / (f_g * (1 - f_g)) ** 2
    pred_var_logk = w_mu ** 2 * (1.0 / n_b + 1.0 / n_s) + (1 - w_mu * f_g) ** 2 * var_lam_src
    # single message stated var (on lam) = pred_var_logk ; ratio to joint truth
    ratio_single = pred_var_logk / true_var_lam

    # adv-2 fork (c): comp-only Var(log k) [drop the w_mu count term] + a SEPARATE spliced message 1/n_s on RNA
    comp_only_k = (1 - w_mu * f_g) ** 2 * var_lam_src
    # separate spliced as independent lambda constraint: precision w_mu^2/(1/n_s) mapped... treat as extra lam prec
    sep_spliced_prec = 1.0 / (w_mu ** 2 * (1.0 / n_b + 1.0 / n_s)) if w_mu > 0 else 0.0
    prec_c = 1.0 / comp_only_k + sep_spliced_prec
    ratio_c = (1.0 / prec_c) / true_var_lam
    # adv-2 fork (b): FULL Var(log k) + separate spliced (double counts spliced)
    prec_b = 1.0 / pred_var_logk + sep_spliced_prec
    ratio_b = (1.0 / prec_b) / true_var_lam

    print(f"\n[{name}] f_g={f_g} tau={tau_lam} n_b={n_b} n_s={n_s} w_mu={w_mu}")
    print(f"  true Var(lam_dst)={true_var_lam:.5g}   emp Var(log k)={var_logk_emp:.5g}  pred Var(log k)={pred_var_logk:.5g} ({abs(pred_var_logk-var_logk_emp)/var_logk_emp:.1%})")
    print(f"  Q1 two-message ratio  = {ratio_two:.3f}   ({1/ratio_two:.2f}x over-confident)")
    print(f"  Q2 single-lambda ratio= {ratio_single:.3f}   ({'OK' if abs(ratio_single-1)<0.08 else 'OFF'})")
    print(f"  Q3 fork(b) full-k+sep = {ratio_b:.3f} (double-count spliced)   fork(c) comp-k+sep = {ratio_c:.3f} (drop cov)")
    return ratio_two, ratio_single


print("=" * 100)
print("Q1/Q2/Q3  COMBINE over-confidence + single-lambda fix + adv-2 spliced fork")
print("=" * 100)
# pure-composition symmetric control (der1/2 claim exactly 2x)
combine_case("pure-comp symmetric", f_g=0.5, tau_lam=60.0, n_b=10**9, n_s=10**9, w_mu=0.0, E_g=150.0, E_r=150.0)
# spliced-free, asymmetric f_g
combine_case("spliced-free fg=0.4", f_g=0.40, tau_lam=60.0, n_b=10**9, n_s=10**9, w_mu=0.0)
# light spliced (adv-1's 3.6x regime?)
combine_case("light-spliced fg=0.4", f_g=0.40, tau_lam=60.0, n_b=800, n_s=40, w_mu=0.15)
# substantial spliced
combine_case("subst-spliced fg=0.4", f_g=0.40, tau_lam=60.0, n_b=800, n_s=600, w_mu=0.50)
# heavy spliced (adv-1: single-lambda still holds?)
combine_case("heavy-spliced fg=0.4", f_g=0.40, tau_lam=60.0, n_b=800, n_s=200, w_mu=0.90)

print("\n" + "=" * 100)
print("Q4  ANCHOR: Var(log k) singular at f_g->1 while divM gDNA density mode is finite")
print("=" * 100)
for f_g in (0.9, 0.99, 0.999):
    tau = 60.0
    n = 400
    var_fg = (f_g * (1 - f_g)) ** 2 / tau
    fg = beta_draws(f_g, var_fg, N)
    Mb = lognormal(1.0, 1.0 / n, N)
    E_g, E_r = 150.0, 200.0
    rho_g = fg * Mb / E_g
    rho_R = (1 - fg) * Mb / E_r
    logk = np.log(rho_g / np.maximum(rho_R, 1e-300))
    # divM gDNA density mode: message on log f_g = log(rho_g * E_g / M). here mass_dst = Mb (its own).
    logfg_msg = np.log(rho_g * E_g / Mb)   # = log fg exactly => var = Var(log f_g)
    var_logk = float(np.var(logk))
    var_gmode = float(np.var(logfg_msg))
    pred_gmode = (1 - f_g) ** 2 / tau  # + 1/n would be there if M didn't cancel; here matched-own so M cancels
    print(f"  f_g={f_g}: Var(log k)={var_logk:.4g} (BLOWS UP)   Var(divM gDNA mode)={var_gmode:.4g}  pred (1-fg)^2/tau={pred_gmode:.4g}")

print("\n" + "=" * 100)
print("Q5  adv-1 M1 near-wall: Var(log f_g)=var_fg/f_g^2 for a small+uncertain MINORITY component")
print("=" * 100)
for f_g, var_fg in ((0.40, 0.004), (0.08, 0.004), (0.04, 0.01), (0.95, 0.002)):
    fg = beta_draws(f_g, var_fg, N)
    emp = float(np.var(np.log(fg)))
    pred = var_fg / f_g ** 2
    print(f"  f_g={f_g} var_fg={var_fg} (CV={np.sqrt(var_fg)/f_g:.2f}): pred Var(log f_g)={pred:.4g}  emp={emp:.4g}  rel {abs(pred-emp)/emp:.1%}")
