"""AMBIG study, step 8 — DERIVATION + MC ARBITER: the AMBIG lambda-likelihood, non-degenerate.

OWNER'S CORRECTION (2026-07-25). The plug-in d_hat = (p_obs-1/2)/(kappa-1/2) is degenerate: kappa is FITTED
(posterior Beta(n_same+1, n_opp+1)) and on unstranded data |kappa_hat - 1/2| is SMALLER THAN ITS OWN SD, yet
it sits in a denominator. Measured: |k-.5|/sd(k) = 0.0-1.3 on ss_0.50 and med|d_hat| up to 118, 99 % of mass
with |d_hat| > 1. The constraint is real but it is STATISTICS, not arithmetic, and it has a distribution.

THE DERIVATION. Never form d_hat. The node's strand count is
    u_pos ~ BetaBin(n, p, omega),   p = 1/2 + (kappa-1/2)*d,   d = (1-f_g)*tau,   tau in [-1,1]
so the likelihood depends on (f_g, kappa) ONLY through the single scale

    A(f_g, kappa) = (1 - f_g) * |kappa - 1/2|        <-- the achievable strand-excess amplitude

(the sign of kappa-1/2 is absorbed by tau's symmetry). The lambda-factor is then the tilt- and kappa-marginal

    L(f_g) = INT INT BetaBin(u_pos; n, 1/2 + A*tau, omega) * pi(tau) * pi(kappa) d tau d kappa

with pi(tau) the solver's existing arcsine reference (uniform in theta = arcsin tau) and pi(kappa) the Beta
posterior the strand model already fits. NO DIVISION ANYWHERE, so there is no degeneracy: as kappa -> 1/2,
A -> 0 for EVERY f_g, and L becomes flat -- "no information" is recovered as a LIMIT, not as a guard.

In the sharp limit (sigma_e << A, kappa known) this has the closed form
    p(f_g) ~ 1 / sqrt( (1-f_g)^2 - d^2 )        on f_g in [0, 1-|d|]
-- an integrable spike at the bound, with median 1 - |d|*cosh( arccosh(1/|d|) / 2 ).

Claims to arbitrate:
  M9a  the likelihood depends on (f_g, kappa) only through A            (exact identity)
  M9b  sharp-limit closed form + its median                             (vs direct numerical marginalization)
  M9c  the UNSTRANDED limit is FLAT once kappa is marginalized          (the owner's point; no guard)
  M9d  frequentist calibration of the posterior median over resampling
  M9e  the plug-in d_hat estimator BLOWS UP where the marginal does not (the defect, quantified)

    OMP_NUM_THREADS=1 python scratchpad/ambig_8_mc.py
"""

from __future__ import annotations

import numpy as np
from scipy.special import betaln, logsumexp

rng = np.random.default_rng(20260725)
_EPS = 1e-12


def _report(name, pred, emp, tol=0.05):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    print(f"  {'OK ' if rel < tol else '***'} {name:<58} pred {pred:12.6g}  emp {emp:12.6g}  rel {rel:7.2%}")
    return rel < tol


def bb_logpmf(k, n, p, omega):
    """Beta-binomial log-pmf with mean p and overdispersion omega (omega -> 0 gives the Binomial)."""
    k = np.asarray(k, float)
    p = np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    if omega <= 0:
        return k * np.log(p) + (n - k) * np.log1p(-p)
    s = (1.0 - omega) / omega
    a, b = p * s, (1.0 - p) * s
    return betaln(k + a, n - k + b) - betaln(a, b)


def tilt_grid(K):
    """The solver's tilt reference: uniform in theta = arcsin(tau) (the Berger-Bernardo conditional)."""
    th = np.linspace(-np.pi / 2, np.pi / 2, K)
    return np.sin(th)


def lam_logfactor(u_pos, n, fg_grid, *, kappa_draws, omega, K_t=2001):
    """log L(f_g): the tilt- and kappa-marginal strand log-likelihood on an f_g grid. Fully in log space —
    the counts here reach 2e6, so any linear accumulation underflows. No division by (kappa-1/2) anywhere."""
    tau = tilt_grid(K_t)[None, :]
    acc = []
    for kap in np.atleast_1d(kappa_draws):
        A = ((1.0 - fg_grid)[:, None]) * abs(kap - 0.5)
        lp = bb_logpmf(u_pos, n, 0.5 + A * tau, omega)  # (F, K_t)
        acc.append(logsumexp(lp, axis=1) - np.log(K_t))
    return logsumexp(np.asarray(acc), axis=0) - np.log(len(acc))


def lam_factor(u_pos, n, fg_grid, *, kappa_draws, omega, K_t=2001):
    """The normalized posterior density on the f_g grid (flat reference)."""
    lg = lam_logfactor(u_pos, n, fg_grid, kappa_draws=kappa_draws, omega=omega, K_t=K_t)
    L = np.exp(lg - lg.max())
    return L / max(np.trapezoid(L, fg_grid), _EPS)


def _median(L, grid):
    c = np.cumsum(L)
    if not np.isfinite(c[-1]) or c[-1] <= 0:
        return np.nan
    return grid[np.searchsorted(c / c[-1], 0.5)]


def simulate(n, f_g, tau_true, kappa, omega=0.0, size=1):
    """Draw u_pos for a node with composition (f_g, tau) under library strand-specificity kappa."""
    d = (1.0 - f_g) * tau_true
    p = 0.5 + (kappa - 0.5) * d
    if omega <= 0:
        return rng.binomial(n, p, size=size).astype(float)
    s = (1.0 - omega) / omega
    return rng.binomial(n, rng.beta(p * s, (1 - p) * s, size=size)).astype(float)


FG = np.linspace(1e-4, 1 - 1e-4, 1200)
ok_all = True

print("# AMBIG lambda-likelihood — MC arbiter\n")

# ── M9a: the likelihood depends on (f_g, kappa) ONLY through A = (1-f_g)|kappa-1/2| ──────────────────────
print("═══ M9a  reduction to the single scale A = (1-f_g)*|kappa-1/2| ═══\n")
n, u = 4000.0, 2300.0
for kA, kB in ((0.01, 0.99), (0.20, 0.80), (0.35, 0.72)):
    # pick f_g's so that A matches exactly:  (1-fa)|kA-.5| == (1-fb)|kB-.5|
    fa = 0.30
    fb = 1.0 - (1.0 - fa) * abs(kA - 0.5) / abs(kB - 0.5)
    if not (0 < fb < 1):
        continue
    La = lam_factor(u, n, np.array([fa]), kappa_draws=[kA], omega=0.0)
    Lb = lam_factor(u, n, np.array([fb]), kappa_draws=[kB], omega=0.0)
    # compare the raw (unnormalized) values via the tilt-marginal directly
    tau = tilt_grid(2001)
    va = np.exp(bb_logpmf(u, n, 0.5 + (1 - fa) * abs(kA - 0.5) * tau, 0.0) - bb_logpmf(u, n, u / n, 0.0)).mean()
    vb = np.exp(bb_logpmf(u, n, 0.5 + (1 - fb) * abs(kB - 0.5) * tau, 0.0) - bb_logpmf(u, n, u / n, 0.0)).mean()
    ok_all &= _report(f"M9a kappa={kA:.2f} (f_g={fa:.3f}) == kappa={kB:.2f} (f_g={fb:.3f})", va, vb, tol=1e-9)

# ── M9b: the sharp-limit closed form and its median ──────────────────────────────────────────────────────
print("\n═══ M9b  sharp limit:  p(f_g) ∝ 1/sqrt((1-f_g)^2 - d^2),  median = 1 - |d|cosh(arccosh(1/|d|)/2) ═══\n")
for d_true, kappa in ((0.084, 0.01), (0.206, 0.01), (0.400, 0.05)):
    n = 2_000_000.0  # sharp limit: sigma_e -> 0
    p = 0.5 + (kappa - 0.5) * d_true
    u = p * n
    L = lam_factor(u, n, FG, kappa_draws=[kappa], omega=0.0)
    med_num = _median(L, FG)
    m = abs(d_true)
    med_cf = 1.0 - m * np.cosh(np.arccosh(1.0 / m) / 2.0)
    ok_all &= _report(f"M9b median  |d|={m:.3f}  (bound B={1 - m:.3f})", med_cf, med_num, tol=0.02)

# ── M9c: THE OWNER'S POINT — the unstranded limit is FLAT once kappa is marginalized ─────────────────────
print("\n═══ M9c  unstranded: kappa marginalized over its Beta posterior ⇒ the factor is FLAT ═══\n")
n_obs = 321_205.0  # spliced observations backing kappa (measured on the suite)
for kap_hat, lab in ((0.49969, "ss_0.50 capON  (|k-.5|/sd = 0.4)"), (0.50004, "ss_0.50 capOFF (|k-.5|/sd = 0.0)")):
    sd_k = np.sqrt(kap_hat * (1 - kap_hat) / (n_obs + 3.0))
    n, f_g_true = 1000.0, 0.60
    u = simulate(n, f_g_true, 0.9, kap_hat, size=1)[0]
    a = kap_hat * (n_obs + 2.0)
    b = (1.0 - kap_hat) * (n_obs + 2.0)
    kdraw = rng.beta(a, b, size=400)
    L_plug = lam_factor(u, n, FG, kappa_draws=[kap_hat], omega=0.0)
    L_marg = lam_factor(u, n, FG, kappa_draws=kdraw, omega=0.0)
    rng_plug = L_plug.max() / max(L_plug.min(), _EPS)
    rng_marg = L_marg.max() / max(L_marg.min(), _EPS)
    d_plug = (u / n - 0.5) / (kap_hat - 0.5)
    print(f"      {lab}:  sd(kappa)={sd_k:.5f}   plug-in |d_hat| = {abs(d_plug):8.2f}  "
          f"(B = {max(0.0, 1 - abs(d_plug)):.3f})")
    print(f"        L range (max/min): plug-in kappa {rng_plug:9.3f}   kappa-marginalized {rng_marg:9.3f}")
    ok_all &= _report(f"M9c FLAT after marginalizing kappa [{lab[:7]}]", 1.0, rng_marg, tol=0.05)

# ── M9d: calibration of the posterior median over repeated sampling ──────────────────────────────────────
print("\n═══ M9d  frequentist behaviour of the marginal posterior median (stranded, kappa well determined) ═══\n")
kappa = 0.00989
n_obs = 320_676.0
a, b = kappa * (n_obs + 2.0), (1 - kappa) * (n_obs + 2.0)
for f_g_true, tau_true, n in ((0.79, 1.00, 1500.0), (0.79, 0.60, 1500.0), (0.44, 1.00, 800.0)):
    meds, bnds = [], []
    for _ in range(120):
        u = simulate(n, f_g_true, tau_true, kappa, size=1)[0]
        L = lam_factor(u, n, FG, kappa_draws=rng.beta(a, b, size=8), omega=0.0)
        meds.append(_median(L, FG))
        bnds.append(max(0.0, 1.0 - abs((u / n - 0.5) / (kappa - 0.5))))
    meds, bnds = np.array(meds), np.array(bnds)
    print(f"      truth f_g={f_g_true:.2f} tau={tau_true:.2f} n={n:.0f} | marginal median "
          f"{meds.mean():.4f} ± {meds.std():.4f}  (bias {meds.mean() - f_g_true:+.4f}) | "
          f"plug-in bound {bnds.mean():.4f} (bias {bnds.mean() - f_g_true:+.4f})")

# ── M9e: where the plug-in blows up and the marginal does not ────────────────────────────────────────────
print("\n═══ M9e  the defect, quantified: plug-in vs marginal as kappa -> 1/2 ═══\n")
print(f"      {'kappa':>9}{'|k-.5|/sd':>11}{'plug |d|':>10}{'plug B':>8}{'marg med':>10}{'truth':>8}")
n_obs = 321_205.0
for kap_hat in (0.010, 0.200, 0.400, 0.470, 0.490, 0.49969, 0.50004):
    sd_k = np.sqrt(kap_hat * (1 - kap_hat) / (n_obs + 3.0))
    a, b = kap_hat * (n_obs + 2.0), (1 - kap_hat) * (n_obs + 2.0)
    n, f_g_true = 1200.0, 0.70
    u = simulate(n, f_g_true, 1.0, kap_hat, size=1)[0]
    L = lam_factor(u, n, FG, kappa_draws=rng.beta(a, b, size=200), omega=0.0)
    dpl = abs((u / n - 0.5) / (kap_hat - 0.5))
    print(f"      {kap_hat:>9.5f}{abs(kap_hat - 0.5) / sd_k:>11.1f}{dpl:>10.3f}"
          f"{max(0.0, 1 - dpl):>8.3f}{_median(L, FG):>10.3f}{f_g_true:>8.3f}")

print(f"\n{'ALL CHECKS PASSED' if ok_all else '*** SOME CHECKS FAILED'}")
