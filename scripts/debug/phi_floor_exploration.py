"""Explore the numerical-stability landscape of phi to pick CalibrationConfig.phi_floor.

phi (NB dispersion = Gamma exposure-prior spread) enters the calibrator only as
r = 1/phi, in three places:

  (A) exposure posterior:  alpha = 1/phi + M_g_tot,  beta = 1/phi + rho0*Leff
                           omega = alpha/beta,  log_omega_var = 1/alpha
  (B) count channel:       NB(n_u; r=1/phi, p=r/(r+mu))   -> Poisson(mu) as phi->0
  (C) phi Newton M-step:   gradient uses digamma(n+r) - digamma(r), etc.

As phi -> 0 (r -> inf) the model -> Poisson. The "floor" is the smallest phi for
which all three are still computed accurately. Below it, decreasing phi changes
nothing computable and only invites precision noise / data being rounded away.

We characterize each failure boundary empirically and pick the floor with margin.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import nbinom, poisson
from scipy.special import digamma, gammaln

EPS = np.finfo(np.float64).eps  # 2.22e-16
print(f"float64 eps = {EPS:.3e}, tiny = {np.finfo(np.float64).tiny:.3e}\n")

phis = np.array([10**-k for k in range(0, 16)], dtype=np.float64)  # 1e0 .. 1e-15

# Representative gDNA means and observed counts (rho0~1e-3, Leff 1e2-1e4, omega 0.1-10):
mus = np.array([0.1, 1.0, 10.0, 100.0, 1000.0])
ns = np.array([0, 1, 10, 100, 1000])

# ---- (A) data-retention: does (1/phi + M) still 'see' the data M? -------------
print("=" * 78)
print("(A) Exposure posterior: smallest M_g_tot still resolvable in alpha=1/phi+M")
print("    (M is lost when 1/phi + M == 1/phi in float64)")
print(f"{'phi':>10} {'1/phi':>12} {'M=0.5 kept?':>12} {'M=1 kept?':>10} {'M=10 kept?':>11}")
for phi in phis:
    r = 1.0 / phi
    keeps = [(r + M) != r for M in (0.5, 1.0, 10.0)]
    print(f"{phi:>10.0e} {r:>12.3e} {str(keeps[0]):>12} {str(keeps[1]):>10} {str(keeps[2]):>11}")

# ---- (B) NB count channel accuracy vs the Poisson limit -----------------------
print("\n" + "=" * 78)
print("(B) NB(n; r=1/phi, p) vs Poisson(n; mu): max |logpmf diff| over (n,mu) grid")
print("    As phi->0 this should DECREASE monotonically; precision noise makes it")
print("    plateau/erratic. The plateau floor is where smaller phi buys nothing.")
print(f"{'phi':>10} {'1/phi':>12} {'max|NB-Poisson|':>18} {'max rel':>12}")
prev = None
for phi in phis:
    r = 1.0 / phi
    max_abs = 0.0
    max_rel = 0.0
    for mu in mus:
        p = r / (r + mu)
        nb = nbinom.logpmf(ns, r, p)
        po = poisson.logpmf(ns, mu)
        d = np.abs(nb - po)
        max_abs = max(max_abs, float(np.nanmax(d)))
        rel = d / np.maximum(np.abs(po), 1e-300)
        max_rel = max(max_rel, float(np.nanmax(rel)))
    trend = ""
    if prev is not None:
        trend = "  (down)" if max_abs < prev else "  <-- not decreasing (precision noise)"
    prev = max_abs
    print(f"{phi:>10.0e} {r:>12.3e} {max_abs:>18.3e} {max_rel:>12.3e}{trend}")

# ---- (C) phi Newton gradient term: digamma(n+r) - digamma(r) cancellation -----
print("\n" + "=" * 78)
print("(C) Newton term digamma(n+r)-digamma(r) vs accurate value; large r cancels")
print("    accurate(n,r) = sum_{j=0}^{n-1} 1/(r+j); compare to digamma difference")
print(f"{'phi':>10} {'1/phi':>12} {'max rel err (n<=1000)':>22}")
for phi in phis:
    r = 1.0 / phi
    worst = 0.0
    for n in ns:
        if n == 0:
            continue
        approx = digamma(n + r) - digamma(r)
        exact = float(np.sum(1.0 / (r + np.arange(n))))
        worst = max(worst, abs(approx - exact) / max(abs(exact), 1e-300))
    print(f"{phi:>10.0e} {r:>12.3e} {worst:>22.3e}")

# ---- (D) gammaln cancellation magnitude (informational) -----------------------
print("\n" + "=" * 78)
print("(D) gammaln(r) magnitude (digits lost in gammaln(n+r)-gammaln(r) ~ log10)")
print(f"{'phi':>10} {'1/phi':>12} {'gammaln(1/phi)':>16} {'~digits lost':>13}")
for phi in phis:
    r = 1.0 / phi
    g = gammaln(r)
    digits_lost = max(0.0, np.log10(max(g, 1.0)) - np.log10(max(100 * np.log(max(r, 2)), 1.0)))
    print(f"{phi:>10.0e} {r:>12.3e} {g:>16.3e} {digits_lost:>13.1f}")
