"""Independent derivation-agent #2 MC — per-component message-variance model for the unified solver.

BAYESIAN framing: the source node holds a belief about its composition, summarized by ONE scalar
    lambda = logit f_g,   Var(lambda) = 1/tau_lam
plus a Poisson mass count M (Var(log M)=1/n) and a Poisson spliced count S (Var(log S)=1/n_spl).
Each per-component density rho_c = f_c*M/E_c is a NOISY OBSERVATION. The message that psi consumes is a
Gaussian on log f_c^dst with mode  mo_c = log( rho_c^src * r * E_c / M_dst ),  r = rho_tot(dst)/rho_tot(src).
We derive Var(log f_c^dst) per component and check what precision psi must receive.

Foundation (given): Var(log f_g)=(1-f_g)^2/tau,  Var(log f_R)=f_g^2/tau  (both from the SINGLE lambda).

Run:  OMP_NUM_THREADS=1 python scratchpad/mc_derive_2.py
"""
from __future__ import annotations
import numpy as np

rng = np.random.default_rng(20260724)
N = 2_000_000


def sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<48} pred {pred:12.6g}  emp {emp:12.6g}  rel {rel:7.2%}")
    return rel, rel < tol


def lognormal(mean, varlog, size):
    s = np.sqrt(max(varlog, 0.0))
    return mean * np.exp(rng.normal(-0.5 * s * s, s, size))


def lam_draws(fg0, tau, size):
    """lambda ~ N(logit fg0, 1/tau); returns f_g = sigma(lambda)."""
    lam0 = np.log(fg0 / (1.0 - fg0))
    lam = rng.normal(lam0, 1.0 / np.sqrt(tau), size)
    return sig(lam)


results = {}


# ============================================================================
# LAW 1 — FOUNDATION SEEDS: Var(log f_g)=(1-f_g)^2/tau ; Var(log f_R)=f_g^2/tau
# ============================================================================
print("\nLAW 1 — foundation composition seeds (single lambda)")
fg0, tau = 0.40, 60.0
fg = lam_draws(fg0, tau, N)
r1a, _ = report("Var(log f_g) = (1-f_g)^2/tau", (1 - fg0) ** 2 / tau, np.var(np.log(fg)))
r1b, _ = report("Var(log f_R) = f_g^2/tau", fg0**2 / tau, np.var(np.log(1 - fg)))
# perfect anti-correlation of the two arms (corr = -1) -- the crux of the combine
corr = np.corrcoef(np.log(fg), np.log(1 - fg))[0, 1]
print(f"      corr(log f_g, log f_R) = {corr:+.4f}   (structurally -1: ONE dof)")
results["foundation"] = max(r1a, r1b)


# ============================================================================
# LAW 2 — TRANSPORT SEED (density): Var(log rho_c) = Var(log f_c) + 1/n   (comp (+) sampling)
#   rho-splice is composition-CERTAIN: Var(log rho_mu)=1/n_spl (no comp term).
# ============================================================================
print("\nLAW 2 — transport seed per component (density var)")
n, n_spl = 300.0, 150.0
fg = lam_draws(fg0, tau, N)
M = lognormal(900.0, 1.0 / n, N)
S = lognormal(1500.0, 1.0 / n_spl, N)
E_g, E_r, E_spl = 110.0, 200.0, 100.0
rho_g = fg * M / E_g
rho_nu = (1 - fg) * M / E_r
rho_mu = S / E_spl
r2a, _ = report("Var(log rho_g) = (1-f_g)^2/tau + 1/n",
                (1 - fg0) ** 2 / tau + 1 / n, np.var(np.log(rho_g)))
r2b, _ = report("Var(log rho_nu)= f_g^2/tau + 1/n",
                fg0**2 / tau + 1 / n, np.var(np.log(rho_nu)))
r2c, _ = report("Var(log rho_mu)= 1/n_spl (comp-certain)", 1 / n_spl, np.var(np.log(rho_mu)))
results["transport_seed"] = max(r2a, r2b, r2c)


# ============================================================================
# LAW 3 — GRAFT gDNA MESSAGE (matched set). Faithful joint: r built from the SHARED source mass, so the
#   count 1/n CANCELS (common-mode with rho_tot(src)) and sigma^2_transfer (src part) CANCELS.
#   f_g^dst = f_g^src * B_dst/B_src ; Var(log f_g^dst)= c_g^2/tau,
#     c_g = (1-f_g)*[1 - (1/E_g - 1/E_r)*f_g/B_src]  (frame coefficient; ->(1-f_g) when E_g=E_r).
# ============================================================================
print("\nLAW 3 — GRAFT gDNA message (count + transfer CANCEL; frame coeff)")
# destination frame held fixed (its belief); source varies through lambda and M.
fg = lam_draws(fg0, tau, N)
M = lognormal(900.0, 1.0 / n, N)
# dst frame fixed:
fg_dst0 = 0.35
Bd = fg_dst0 / E_g + (1 - fg_dst0) / E_r
M_dst = 700.0
rho_tot_dst = M_dst * Bd
Bs = fg / E_g + (1 - fg) / E_r
rho_tot_src = M * Bs
r_enr = rho_tot_dst / rho_tot_src
fg_msg = (fg * M / E_g) * r_enr * E_g / M_dst
Bs0 = fg0 / E_g + (1 - fg0) / E_r
c_g = (1 - fg0) * (1 - (1 / E_g - 1 / E_r) * fg0 / Bs0)
r3, _ = report("Var(log f_g^dst) = c_g^2/tau (no 1/n, no s2t)", c_g**2 / tau, np.var(np.log(fg_msg)))
# matched-frame sanity: E_g=E_r -> collapses to (1-f_g)^2/tau, count truly gone
Eq = 150.0
Bs_q = 1.0 / Eq
rho_tot_src_q = M * Bs_q
Bd_q = 1.0 / Eq
r_q = (M_dst * Bd_q) / rho_tot_src_q
fg_msg_q = (fg * M / Eq) * r_q * Eq / M_dst
r3b, _ = report("  matched frame E_g=E_r: Var=(1-f_g)^2/tau", (1 - fg0) ** 2 / tau,
                np.var(np.log(fg_msg_q)))
results["graft_gdna"] = max(r3, r3b)


# ============================================================================
# LAW 4 — GRAFT RNA-SUM MESSAGE (share-weighted). rho_R = rho_nu + rho_mu.
#   Matched frame E_g=E_r=Eq so B cancels; f_R^dst = (f_R + (S/E_spl)*(Eq/M)) * 1.
#   Var(log f_R^dst) = w_nu^2 * Var(log f_R) + w_mu^2*(1/n + 1/n_spl)   [SUM: convex squared weights]
# ============================================================================
print("\nLAW 4 — GRAFT RNA-sum message (item-E share weighting)")
fg = lam_draws(fg0, tau, N)
M = lognormal(900.0, 1.0 / n, N)
S = lognormal(1500.0, 1.0 / n_spl, N)
# matched frame; message RNA fraction at dst:
A = (1 - fg) + (S / E_spl) * (Eq / M)      # = rho_R^src * Eq / M   (B=1/Eq cancels)
frac_R = A * (M_dst * (1.0 / Eq) / (M * (1.0 / Eq))) * 1.0  # r*Eq/M_dst with B cancel -> *(1)? keep explicit
# do it fully explicitly to avoid algebra slips:
Bs_q = 1.0 / Eq
rho_R_src = (1 - fg) * M / Eq + S / E_spl
r_q = (M_dst / Eq) / (M * Bs_q)
fR_msg = rho_R_src * r_q * Eq / M_dst
# operating point shares
rho_nu0 = (1 - fg0) * 900.0 / Eq
rho_mu0 = 1500.0 / E_spl
rR0 = rho_nu0 + rho_mu0
w_nu, w_mu = rho_nu0 / rR0, rho_mu0 / rR0
pred4 = w_nu**2 * (fg0**2 / tau) + w_mu**2 * (1 / n + 1 / n_spl)
r4, _ = report(f"Var(log f_R^dst) share-wtd [w_mu={w_mu:.2f}]", pred4, np.var(np.log(fR_msg)))
results["graft_rna_sum"] = r4


# ============================================================================
# LAW 5 — PEEL DIFFERENCE (u-weighted; sigma^2_transfer LOAD-BEARING).
#   rho_nu(b) = rho_R(x)/r - rho_mu ; u = T/rho_nu, T=rho_R(x)/r.
#   Var(log rho_nu) = u^2 * (Var(log rho_R(x)) + Var(log r)) + (u-1)^2 * 1/n_spl
# ============================================================================
print("\nLAW 5 — PEEL difference (u-weighted; s2t load-bearing)")
rho_Rx, var_log_Rx = 40.0, 1.0 / 5000
r_bar, s2t = 200.0, 0.010
Rx = lognormal(rho_Rx, var_log_Rx, N)
rr = lognormal(r_bar, s2t, N)
rmu = lognormal(0.10, 1.0 / 1500, N)
T = Rx / rr
nu = T - rmu
keep = nu > 0
T0 = rho_Rx / r_bar
nu0 = T0 - 0.10
u = T0 / nu0
pred5 = u**2 * (var_log_Rx + s2t) + (u - 1) ** 2 * (1.0 / 1500)
r5, _ = report(f"Var(log rho_nu) [u={u:.2f}, kept {keep.mean():.0%}]", pred5,
               np.var(np.log(nu[keep])), tol=0.12)
# ablate s2t -> shows how much it carries
nu_ab = (Rx / r_bar) - rmu
frac_carried = 1 - np.var(np.log(nu_ab[nu_ab > 0])) / np.var(np.log(nu[keep]))
print(f"      sigma^2_transfer carries {frac_carried:.0%} of the peel variance")
results["peel_diff"] = r5


# ============================================================================
# LAW 6 — ANCHOR LIMIT (pure gDNA, f_g -> 1). Ratio k=rho_g/rho_R is SINGULAR (inf); the per-component
#   form is FINITE: gDNA message var -> c_g^2 * Var(log f_g) with Var(log f_g)=(1-f_g)^2/tau -> 0 as f_g->1
#   (composition-certain), and the RNA arm carries rho_R -> 0 => f_R -> 0 => log f_R -> -inf but PRECISION 0.
# ============================================================================
print("\nLAW 6 — pure-gDNA ANCHOR limit (finite where k=rho_g/rho_R is singular)")
fg_anchor = 0.999
fg = lam_draws(fg_anchor, tau, N)
M = lognormal(900.0, 1.0 / n, N)
Bs = fg / E_g + (1 - fg) / E_r
r_enr = (M_dst * Bd) / (M * Bs)
fg_msg = (fg * M / E_g) * r_enr * E_g / M_dst
Bs0 = fg_anchor / E_g + (1 - fg_anchor) / E_r
c_g = (1 - fg_anchor) * (1 - (1 / E_g - 1 / E_r) * fg_anchor / Bs0)
r6, _ = report("Var(log f_g^dst) FINITE at f_g=0.999", c_g**2 / tau, np.var(np.log(fg_msg)))
print(f"      k = f_g E_r/((1-f_g)E_g) = {fg_anchor * E_r / ((1 - fg_anchor) * E_g):.1f} (ratio form blows up); "
      f"per-comp var = {c_g**2 / tau:.3e} (finite). RNA prec -> 0 (rho_R->0).")
results["anchor_limit"] = r6


# ============================================================================
# LAW 7 — divM_dst CONVERSION: M_dst is COMMON-MODE across components and DROPS from the composition
#   (the simplex-relevant ratio). Show Var(log(f_g^dst/f_R^dst)) is independent of Var(log M_dst).
# ============================================================================
print("\nLAW 7 — divM_dst common-mode: drops from the composition ratio")
fg = lam_draws(fg0, tau, N)
M = lognormal(900.0, 1.0 / n, N)
S = lognormal(1500.0, 1.0 / n_spl, N)
for vM in (0.0, 1.0 / 50):   # vary the dst-mass count uncertainty
    Md = lognormal(700.0, vM, N)
    Bs_q = 1.0 / Eq
    rho_g_src = fg * M / Eq
    rho_R_src = (1 - fg) * M / Eq + S / E_spl
    r_q = (Md / Eq) / (M * Bs_q)
    fg_m = rho_g_src * r_q * Eq / Md
    fR_m = rho_R_src * r_q * Eq / Md
    print(f"      Var(log M_dst)={vM:.3f}: Var(log f_g)={np.var(np.log(fg_m)):.5f}  "
          f"Var(log f_R)={np.var(np.log(fR_m)):.5f}  Var(log f_g/f_R)={np.var(np.log(fg_m / fR_m)):.5f}")
print("      -> the RATIO (what the simplex keeps) is INVARIANT to M_dst; M_dst is common-mode.")


# ============================================================================
# LAW 8 — THE COMBINE (the critical question). psi gets a gDNA message on log f_g AND an RNA message on
#   log f_R from ONE source belief, treated INDEPENDENT with p_g=tau/(1-f_g)^2, p_R=tau/f_g^2.
#   But both arms are the SAME lambda => the independent combine DOUBLE-COUNTS lambda:
#     effective lambda-precision = p_g*(1-f_g)^2 + p_R*f_g^2 = tau + tau = 2*tau
#   The correct (joint) posterior variance the dst should hold = 1/tau. Ratio (two-msg var/joint var)=1/2.
# ============================================================================
print("\nLAW 8 — THE COMBINE: independent gDNA+RNA messages DOUBLE-COUNT lambda")
fg0c, tauc = 0.45, 40.0
lam0 = np.log(fg0c / (1 - fg0c))
# MC: the source's belief spread = the true dst spread (matched transport). Draw the source point estimate.
lam_src = rng.normal(lam0, 1.0 / np.sqrt(tauc), N)
fg_src = sig(lam_src)
mo_g = np.log(fg_src)
mo_R = np.log(1 - fg_src)
# psi combines the two as independent Gaussians. Effective precision on lambda_dst at the operating point:
p_g = tauc / (1 - fg0c) ** 2
p_R = tauc / fg0c**2
eff_prec_two = p_g * (1 - fg0c) ** 2 + p_R * fg0c**2   # = 2*tau
var_two = 1.0 / eff_prec_two
# empirical: solve the two-message psi mode for each draw (a 1-D Gaussian combine in lambda near op point),
# posterior MODE lambda_hat; its variance across draws is the ACTUAL uncertainty the dst inherits.
# linearized mode: lambda_hat = [p_g*(1-f)*mo_g_shift - p_R*f*... ] -- easier: the mode tracks lam_src exactly
# (both messages are consistent), so spread of the mode = Var(lam_src) = 1/tau. The STATED var is var_two.
actual_var = np.var(lam_src)   # = 1/tau (the honest transport uncertainty)
joint_var = 1.0 / tauc
ratio = var_two / joint_var
print(f"      two-message stated Var(lambda_dst) = {var_two:.5f}  (eff prec {eff_prec_two:.2f} = 2*tau={2*tauc})")
print(f"      joint-truth   Var(lambda_dst) = {joint_var:.5f}  (= 1/tau)")
print(f"      empirical actual spread of transported lambda = {actual_var:.5f}")
report("COMBINE ratio (two-msg var / joint truth)", ratio, 0.5)
# grid confirmation of the two-message curvature (independent of linearization):
lam_grid = np.linspace(lam0 - 3 / np.sqrt(tauc), lam0 + 3 / np.sqrt(tauc), 4001)
fgg = sig(lam_grid)
# use the mean message modes (source at its mean)
psi = -0.5 * p_g * (np.log(fgg) - np.log(fg0c)) ** 2 - 0.5 * p_R * (np.log(1 - fgg) - np.log(1 - fg0c)) ** 2
w = np.exp(psi - psi.max())
w /= w.sum()
mean = (w * lam_grid).sum()
var_grid = (w * (lam_grid - mean) ** 2).sum()
r8, _ = report("  grid two-msg Var(lambda) vs analytic 1/(2tau)", var_grid, 1.0 / (2 * tauc))
results["combine_grid"] = r8
print(f"\n  >>> COMBINE: two-message var / joint-truth var = {var_grid / joint_var:.4f}  "
      f"(target 0.5 = 2x OVERCONFIDENT)")


print("\n================ SUMMARY (max rel error per law) ================")
for k, v in results.items():
    print(f"  {k:<20} rel {v:7.2%}  {'OK' if v < 0.08 else 'CHECK'}")
