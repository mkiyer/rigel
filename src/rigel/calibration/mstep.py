"""M-step — fit the library hyperparameters (doc 03 §5).

Each outer iteration, after the E-step soft-allocation and the exposure update,
the M-step re-estimates the library-wide hyperparameters. ``ρ_0`` is a closed
form; the **exposure dispersion** (the NB count overdispersion ≡ the variance of
the per-region gDNA exposure ω) is the proper EM update for the Gamma
exposure-prior shape, solved from the exposure posteriors
(:func:`update_exposure_dispersion`); ``ρ_d_bb`` is a bounded 1-D minimiser on
its negative log-likelihood. ``κ_rna`` and ``ρ_r_bb`` are **not** fit here — they
come from the clean spliced channel (PR 3) and do not change with the EM (PR05
§III.1); ``κ_d = 0.5`` is biology. (The gDNA splice-artifact rate ``ε_s`` was
removed in PR 7 — spliced fragments are RNA; artifacts are handled upstream by
the ``alignable`` splice blacklist.)

Constants (PR05 §II.1, Q6): ``_EXPOSURE_DISPERSION_MAX`` (the exposure-dispersion
ceiling), ``_PI_FLOOR`` (the no-divide-by-zero floor in the π_g-prior update).
The exposure-dispersion lower bound and the BB bounds reuse
``config.exposure_dispersion_floor`` (PR 2) and ``_BB_FLOOR`` (PR 3).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import digamma
from scipy.stats import betabinom

from .estep import _PI_CLIP
from .signature import TS_AMBIG

# Upper bound on the exposure dispersion (the NB count overdispersion ≡ variance
# of the per-region gDNA exposure ω; doc 03 §8). A numerical ceiling only: the
# dispersion enters elsewhere as 1/dispersion, so capping it keeps that channel
# well-conditioned even when the per-region gDNA exposure is extremely
# heterogeneous. (Was ``_PHI_MAX``; see the no-greek-in-code naming preference.)
_EXPOSURE_DISPERSION_MAX = 100.0
# No-divide-by-zero floor inside the π_g-prior update (doc 03 §5.6 / §8).
_PI_FLOOR = 1.0e-9
# gDNA strand dispersion when too few regions to fit (matches the doc 03 §7 init).
_RHO_D_BB_FALLBACK = 0.01
_MIN_OBS_FOR_BB = 2
# Beta-Binomial shape floor keeping the κ_d=0.5 gDNA BB shape params finite and
# positive at the bounds (doc 03 §8). Used by the ρ_d_bb fit below. (The RNA
# strand BB no longer needs a floor — PR 9 made it the posterior-predictive.)
_BB_FLOOR = 1.0e-6


def decodable_node_masks(
    ts_class: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-node decodability masks for the density/dispersion M-steps (PR 8).

    A *region* is decodable iff its transcript-strand class is not AMBIG (NONE /
    POS / NEG — its gDNA/RNA split is recoverable by the strand channel, or it is
    intergenic where gDNA is unstranded). An internal *boundary* is decodable iff
    it joins two same-reference regions and **at least one** of them is decodable:
    a boundary-crossing fragment can be deconvolved whenever one side carries
    strand context (the premise of using crossing reads — and the calibration
    signal of hybrid-capture libraries lives at boundaries, so we keep them).

    Returns ``(decodable_region, decodable_left_boundary, decodable_right_boundary)``
    as ``bool[R]``, indexed by the region the boundary side is attributed to. A
    boundary side attributed to region ``r`` is decodable iff ``r`` itself is
    decodable **or** the same-reference neighbour on the far side is decodable. A
    decodable region's own boundary sides are therefore always kept (so a library
    with no AMBIG region is scored identically to before — terminals included);
    only an AMBIG region's sides depend on its decodable neighbour.
    """
    ts = np.asarray(ts_class)
    ref = np.asarray(ref_id)
    r = ts.shape[0]
    dec_region = ts != TS_AMBIG
    has_internal_left = np.zeros(r, dtype=bool)  # same-ref neighbour to the left
    has_internal_right = np.zeros(r, dtype=bool)  # same-ref neighbour to the right
    if r > 1:
        same = ref[:-1] == ref[1:]  # internal seam between r and r+1
        has_internal_right[:-1] = same
        has_internal_left[1:] = same
    dec_left_bnd = dec_region.copy()
    dec_left_bnd[1:] |= has_internal_left[1:] & dec_region[:-1]  # far (left) neighbour decodable
    dec_right_bnd = dec_region.copy()
    dec_right_bnd[:-1] |= has_internal_right[:-1] & dec_region[1:]  # far (right) neighbour decodable
    return dec_region, dec_left_bnd, dec_right_bnd


def update_rho_0(
    m_g_contained: np.ndarray,
    m_g_left: np.ndarray,
    m_g_right: np.ndarray,
    omega: np.ndarray,
    l_eff: np.ndarray,
    *,
    dec_region: np.ndarray,
    dec_left_bnd: np.ndarray,
    dec_right_bnd: np.ndarray,
) -> float:
    """ρ_0 = Σ_decodable M_g / Σ_decodable(ω·L_eff) (doc 03 §5.1; PR 8 decodable-node).

    Only **decodable** nodes contribute. An AMBIG region's *contained* gDNA is
    withheld (undecodable — no strand split, no crossing), but its *boundary* mass
    is kept whenever the boundary has a decodable neighbour. A region's exposure
    capacity ``ω·L_eff`` enters the denominator iff *any* of its nodes is decodable,
    so boundary mass is never counted in the numerator without matching capacity in
    the denominator. Physical length ``L_eff`` (consistent with the exposure
    posterior, doc 03 §4); the fully FL-corrected node form is deferred (PR 8 doc).
    """
    num = float(
        np.sum(
            np.where(dec_region, m_g_contained, 0.0)
            + np.where(dec_left_bnd, m_g_left, 0.0)
            + np.where(dec_right_bnd, m_g_right, 0.0)
        )
    )
    any_decodable = dec_region | dec_left_bnd | dec_right_bnd
    denom = float(np.sum(np.where(any_decodable, omega * l_eff, 0.0)))
    if denom <= 0.0:
        return 1.0 / max(float(np.sum(l_eff)), 1.0)
    return num / denom


def update_exposure_dispersion(
    omega: np.ndarray,
    log_omega_var: np.ndarray,
    *,
    floor: float,
    ceil: float = _EXPOSURE_DISPERSION_MAX,
) -> float:
    """EM M-step for the exposure dispersion (doc 03 §5.4; oscillation diagnosis).

    The "exposure dispersion" is the NB count overdispersion ≡ the variance of
    the per-region gDNA exposure ``ω`` (``ω ~ Gamma(s, s)`` with shape/rate
    ``s = 1/dispersion``, so ``E[ω]=1``, ``Var[ω]=dispersion``). The **proper EM
    M-step** maximises the expected complete-data log-likelihood, whose only
    dispersion term is the exposure prior ``Σ_r E[log Gamma(ω_r | s, s)]``. Its
    stationary condition is the standard Gamma-shape MLE::

        log s − ψ(s) = mean_r( E[ω_r] − E[log ω_r] ) − 1

    using the per-region exposure posteriors ``Gamma(α_r, β_r)`` that
    :func:`exposure_posterior` already produced — recovered here from its outputs
    (``α_r = 1/log_omega_var_r``, ``β_r = α_r/ω_r``), so ``E[ω_r] = ω_r`` and
    ``E[log ω_r] = ψ(α_r) − log β_r``. Solved by monotone bisection on ``s`` (the
    LHS is positive and strictly decreasing); the dispersion is ``1/s`` clipped to
    ``[floor, ceil]``.

    Unlike the former count-NB fit, this uses **only** the exposure posteriors —
    never the raw counts — so empty (zero-count) regions cannot drive the
    count-misfit feedback that produced the period-2 limit cycle (see
    ``docs/acc_caljointmodel/calibration_oscillation_diagnosis.md``). EM
    guarantees a monotone increase of the objective, hence convergence.
    """
    omega = np.asarray(omega, dtype=np.float64)
    log_omega_var = np.asarray(log_omega_var, dtype=np.float64)
    if omega.size == 0:
        return float(np.clip(1.0, floor, ceil))

    alpha = 1.0 / log_omega_var  # α_post  (since log_omega_var = 1/α_post)
    beta = alpha / omega  # β_post  (since ω = α_post/β_post)
    e_log_omega = digamma(alpha) - np.log(beta)
    # c ≥ 0 always: x − log x ≥ 1 ⇒ mean(E[ω] − E[log ω]) ≥ 1 (Jensen).
    c = float(np.mean(omega - e_log_omega)) - 1.0
    if not np.isfinite(c) or c <= 0.0:
        return floor  # s → ∞ ⇒ dispersion → 0

    # Solve log s − ψ(s) = c by geometric bisection over s ∈ [1/ceil, 1/floor];
    # LHS is strictly decreasing, so the bracket also enforces the [floor, ceil] clip.
    s_lo, s_hi = 1.0 / ceil, 1.0 / floor
    for _ in range(100):
        s_mid = (s_lo * s_hi) ** 0.5
        if np.log(s_mid) - float(digamma(s_mid)) > c:
            s_lo = s_mid  # LHS too large ⇒ need larger s
        else:
            s_hi = s_mid
    s = (s_lo * s_hi) ** 0.5
    return float(np.clip(1.0 / s, floor, ceil))


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


__all__ = [
    "decodable_node_masks",
    "update_rho_0",
    "update_exposure_dispersion",
    "fit_rho_d_bb",
    "update_pi_g_prior",
]
