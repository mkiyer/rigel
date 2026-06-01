"""M-step — fit the library hyperparameters (doc 03 §5).

Each outer iteration, after the E-step soft-allocation and the exposure update,
the M-step re-estimates the library-wide hyperparameters from the soft-allocated
masses. ``ρ_0`` and ``ε_s`` are closed forms; ``φ`` and ``ρ_d_bb`` are single
global scalars found by a bounded 1-D minimiser on their negative log-likelihood
(``scipy.optimize.minimize_scalar``). ``κ_rna`` and ``ρ_r_bb`` are **not** fit
here — they come from the clean spliced channel (PR 3) and do not change with the
EM (PR05 §III.1); ``κ_d = 0.5`` is biology.

Constants (PR05 §II.1, Q6): ``_PHI_MAX`` (NB-overdispersion ceiling),
``_EPS_S_PRIOR`` (the Beta(1,1) unit pseudo-count), ``_PI_FLOOR`` (the
no-divide-by-zero floor in the π_g-prior update). The φ lower bound and the BB
bounds reuse ``config.phi_floor`` (PR 2) and ``_BB_FLOOR`` (PR 3).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import betabinom, nbinom

from .estep import _PI_CLIP
from .strand_balance import _BB_FLOOR

# NB-overdispersion ceiling for the φ minimiser (doc 03 §8 `(1e-6, 1e2)`).
_PHI_MAX = 100.0
# Beta(1,1) unit pseudo-count for the ε_s failsafe (doc 03 §5.5).
_EPS_S_PRIOR = 1.0
# No-divide-by-zero floor inside the π_g-prior update (doc 03 §5.6 / §8).
_PI_FLOOR = 1.0e-9
# gDNA strand dispersion when too few regions to fit (matches the doc 03 §7 init).
_RHO_D_BB_FALLBACK = 0.01
_MIN_OBS_FOR_BB = 2


def update_rho_0(m_g_tot: np.ndarray, omega: np.ndarray, l_eff: np.ndarray) -> float:
    """ρ_0 = Σ M_g_tot / Σ(ω·L_eff) (doc 03 §5.1, closed form, physical length)."""
    denom = float(np.sum(omega * l_eff))
    total = float(np.sum(m_g_tot))
    if denom <= 0.0:
        return 1.0 / max(float(np.sum(l_eff)), 1.0)
    return total / denom


def update_eps_s(pi_g_contained: np.ndarray, n_spliced: np.ndarray) -> float:
    """ε_s = (1 + Σ π_g·n_s) / (1 + Σ n_s) — Beta(1,1) failsafe (doc 03 §5.5)."""
    n_s = np.asarray(n_spliced, dtype=np.float64)
    num = _EPS_S_PRIOR + float(np.sum(pi_g_contained * n_s))
    den = _EPS_S_PRIOR + float(np.sum(n_s))
    return min(max(num / den, _BB_FLOOR), 1.0 - _BB_FLOOR)


def fit_phi(n_u: np.ndarray, mu_g: np.ndarray, m_d_unspl: np.ndarray, *, phi_floor: float) -> float:
    """φ via a bounded 1-D minimiser on the NB NLL (doc 03 §5.4).

    Treats the mixture mean μ = μ_g + μ_d as known (plug-in); a single global
    scalar on ``(phi_floor, _PHI_MAX)``.
    """
    n = np.asarray(n_u, dtype=np.float64)
    mu = np.maximum(
        np.asarray(mu_g, dtype=np.float64) + np.asarray(m_d_unspl, dtype=np.float64), 1e-12
    )
    if float(n.sum()) <= 0.0:
        return phi_floor

    def nll(phi: float) -> float:
        inv = 1.0 / phi
        return -float(np.sum(nbinom.logpmf(n, inv, inv / (inv + mu))))

    res = minimize_scalar(nll, bounds=(phi_floor, _PHI_MAX), method="bounded")
    return float(np.clip(res.x, phi_floor, _PHI_MAX))


def fit_rho_d_bb(k_plus_g: np.ndarray, n_g: np.ndarray) -> float:
    """ρ_d_bb via a bounded 1-D minimiser on the gDNA BB NLL (κ_d=0.5; doc 03 §5.2).

    Consumes the soft-allocated gDNA sense / total counts, rounded to integers
    (Beta-Binomial support). Falls back to the init when fewer than two regions
    carry gDNA.
    """
    k = np.rint(np.maximum(np.asarray(k_plus_g, dtype=np.float64), 0.0)).astype(np.int64)
    n = np.rint(np.maximum(np.asarray(n_g, dtype=np.float64), 0.0)).astype(np.int64)
    mask = n > 0
    k = np.minimum(k[mask], n[mask])
    n = n[mask]
    if k.size < _MIN_OBS_FOR_BB:
        return _RHO_D_BB_FALLBACK

    def nll(rho: float) -> float:
        a = 0.5 * (1.0 - rho) / rho  # symmetric about κ_d = 0.5
        return -float(np.sum(betabinom.logpmf(k, n, a, a)))

    res = minimize_scalar(nll, bounds=(_BB_FLOOR, 1.0 - _BB_FLOOR), method="bounded")
    return float(np.clip(res.x, _BB_FLOOR, 1.0 - _BB_FLOOR))


def update_pi_g_prior(
    omega: np.ndarray, rho_0: float, l_eff_gdna: np.ndarray, n_u: np.ndarray
) -> np.ndarray:
    """Data-driven π_g prior for the next iteration (doc 03 §5.6).

    ``μ_g = ω ρ_0 L_eff_gdna`` (the FL-corrected contained gDNA mean); ``μ_d`` is
    the count excess over it. The ``_PI_FLOOR`` `maximum` is the only soft guard.
    """
    mu_g = omega * rho_0 * np.asarray(l_eff_gdna, dtype=np.float64)
    mu_d = np.maximum(np.asarray(n_u, dtype=np.float64) - mu_g, _PI_FLOOR)
    prior = mu_g / (mu_g + mu_d)
    return np.clip(prior, _PI_CLIP, 1.0 - _PI_CLIP)


__all__ = ["update_rho_0", "update_eps_s", "fit_phi", "fit_rho_d_bb", "update_pi_g_prior"]
