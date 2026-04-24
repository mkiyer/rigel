"""Mappability-aware EM mixture deconvolution for gDNA calibration (v5).

Two-component mixture over the **full genome partition** produced by
``index.build_region_table`` — intergenic, intronic, and exonic atomic
bins alike.  Every region with mappable bp above ``mappable_floor`` is
eligible.  Spliced-read evidence (only possible in exonic regions) is a
*physical* RNA anchor (spliced fragments cannot arise from gDNA).

Latent classes:
  * G (not-expressed; gDNA-only):  k^s_i = 0 (hard).
  * R (expressed; gDNA + RNA):      any k^s_i.

Signals fused as additive LLRs on logit γ_i (unit weights, no rails):

  1. Count LLR — Poisson(λ_G · E_i) vs Poisson-LogNormal marginal
     over an expressed-rate prior μ ~ LogNormal(μ_R, σ_R²).  Uses
     21-node Gauss-HermiteE quadrature on log μ.

  2. Strand LLR — Beta-Binomial under G (symmetric around p=0.5, shape
     κ_G) vs Beta-Binomial under R integrated over the RNA-fraction
     prior ρ ~ Uniform(0, 1), with mean p(ρ) = ρ·SS + (1-ρ)·0.5 and
     shape κ_R.  κ_G and κ_R are fit independently each M-step by
     γ-weighted (resp. (1-γ)-weighted) soft MLE on the full data.
     17-node Gauss-Legendre quadrature on ρ.

  3. Fragment-length LLR — shape-normalised histogram ratio
     (γ-weighted gDNA FL vs (1-γ)-weighted RNA FL), from iter 1.

See ``docs/calibration/strand_channel_theory_and_redesign.md`` and
``docs/calibration/strand_channel_implementation_plan.md``.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.special import gammaln as _vlgamma

from ..frag_length_model import FragmentLengthModel
from ._fl_model import build_gdna_fl_model


logger = logging.getLogger(__name__)

_EPS: float = 1e-12
_LOG_2PI: float = math.log(2.0 * math.pi)

# Gauss-HermiteE nodes/weights: ∫ f(x) exp(-x²/2) dx ≈ Σ w_i f(x_i).
_GH_NODES, _GH_WEIGHTS = np.polynomial.hermite_e.hermegauss(21)
# Gaussian expectation helpers: E_{x~N(0,1)}[f(x)] = Σ (w_i/√(2π)) f(x_i).
_GH_LOGW = np.log(_GH_WEIGHTS) - 0.5 * _LOG_2PI

# 9-node Gauss-Legendre on [0, 1] for the ρ integral under R-class.
# The integrand BB(k; n, κ·p(ρ), κ·(1-p(ρ))) is smooth in ρ over
# [0.5, SS] with SS ≤ 1 — 9 Gauss-Legendre nodes integrate polynomials
# up to degree 17 exactly and are overkill for the per-region LLR
# accuracy needed by the downstream E-step (sigmoid saturates well
# before 0.01 nats per-region error).  Halving node count roughly
# halves the dominant R-class BB log-likelihood cost.
# np.polynomial.legendre.leggauss returns nodes on [-1, 1] with weights
# summing to 2; we remap to ρ = (x+1)/2 with w' = w/2 so Σ w' = 1.
_RHO_NODES_X, _RHO_WEIGHTS_X = np.polynomial.legendre.leggauss(9)
_RHO_NODES: np.ndarray = 0.5 * (_RHO_NODES_X + 1.0)          # shape (9,)
_RHO_LOGW: np.ndarray = np.log(0.5 * _RHO_WEIGHTS_X)         # shape (9,)

# Beta-Binomial dispersion bracket.  Lower bound: near-U Beta(0.05, 0.05),
# extreme overdispersion.  Upper bound: effectively Binomial (BB variance
# differs from Binomial by < 0.1% at κ = 1000).
_KAPPA_LO: float = 0.1
_KAPPA_HI: float = 1000.0


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EMFit:
    """EM fit results for the two-component gDNA mixture.

    Capture-class (composite-Poisson) mode: when
    ``capture_class_mode`` is True (an ``(E_on, E_off)`` attribution
    was supplied), ``lam_G_on`` and ``lam_G_off`` are fit
    independently via a Poisson-superposition nested EM and
    ``lam_G`` reports the effective coverage-weighted rate
    ``(λ_on · ΣE_on + λ_off · ΣE_off) / ΣE``.  In single-rate mode the
    three rate fields are identical.

    The R-class expressed-rate prior is kept global and unscaled: a
    single ``(μ_R, σ_R²)`` LogNormal applies to on- and off-target
    regions alike — the prior's cross-transcript variance absorbs any
    residual probe-enrichment in RNA counts (Phase 3a Option (a)).
    """

    # Effective / legacy global gDNA rate (fragments per base).  In
    # capture mode this is the coverage-weighted average of
    # ``lam_G_on`` and ``lam_G_off``; in non-capture mode they all
    # coincide.  Kept for backward compatibility with downstream
    # callers that read a single rate.
    lam_G: float
    # LogNormal(μ_R, σ_R²) parameters of the expressed-rate prior
    mu_R: float
    sigma_R: float
    # Mixing proportion π = P(z=G)
    pi: float
    # π conditional on k^s = 0 (the soft prior for the E-step)
    pi_soft: float
    # Per-region posterior γ_i = P(z=G | data)
    gamma: np.ndarray
    # Per-region strand-only posterior γ̃^strand_i = sigmoid(logit(0.5) + LLR_strand_i).
    # Diagnostic: what does strand alone say about this region?
    region_gamma_strand: np.ndarray
    # Beta-Binomial dispersions (separate for each class).
    kappa_G: float = 0.0
    kappa_R: float = 0.0
    # Convergence
    n_iter: int = 0
    converged: bool = False
    # Effective counts
    n_eligible: int = 0
    n_soft: int = 0
    n_spliced_hard: int = 0
    # Iteration trace (optional diagnostics)
    history: list[dict] = field(default_factory=list)

    # -- Composite-Poisson capture-class rates --
    # In non-capture mode both equal ``lam_G``.  In capture mode
    # ``lam_G_on`` is the rate per on-target mappable bp and
    # ``lam_G_off`` the rate per off-target mappable bp.
    lam_G_on: float = 0.0
    lam_G_off: float = 0.0
    capture_class_mode: bool = False


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _poisson_logpmf(k: np.ndarray, rate: np.ndarray) -> np.ndarray:
    """Log-pmf of Poisson(k | rate).

    Numerically stable; handles rate==0 (returns 0 when k==0, -inf else)
    and non-integer k (treats k as real via lgamma(k+1)).
    """
    k = np.asarray(k, dtype=np.float64)
    rate = np.asarray(rate, dtype=np.float64)
    out = np.full(np.broadcast_shapes(k.shape, rate.shape), -np.inf, dtype=np.float64)
    safe = rate > 0
    # log(Poisson) = k*log(rate) - rate - lgamma(k+1)
    lograte = np.zeros_like(rate)
    lograte[safe] = np.log(rate[safe])
    vals = k * lograte - rate - _vlgamma(k + 1.0)
    vals = np.asarray(vals, dtype=np.float64)
    np.copyto(out, vals, where=np.broadcast_to(safe, out.shape))
    # Degenerate branch: rate == 0 → logPMF = 0 if k==0 else -inf
    zero_rate = ~safe & (k == 0.0)
    out = np.where(np.broadcast_to(zero_rate, out.shape), 0.0, out)
    return out


def _logsumexp(a: np.ndarray, axis: int = 0) -> np.ndarray:
    # Pure-numpy idiom — benchmarked ~2× faster than scipy's
    # ``scipy.special.logsumexp`` for (K, N) arrays with K ≪ N
    # because scipy's version introduces array-API ``.astype`` and
    # ``where``-mask copies that dominate for moderate N.
    m = np.max(a, axis=axis, keepdims=True)
    finite = np.isfinite(m)
    m_safe = np.where(finite, m, 0.0)
    s = np.log(np.sum(np.exp(a - m_safe), axis=axis, keepdims=True))
    return np.squeeze(m_safe + s, axis=axis)


# ---------------------------------------------------------------------------
# Count LLR (Poisson vs Poisson-LogNormal)
# ---------------------------------------------------------------------------


def _count_llr_poisson_ln(
    k_u: np.ndarray,
    E: np.ndarray,
    lam_G: float,
    mu_R: float,
    sigma_R: float,
) -> np.ndarray:
    """Per-region count LLR: log P(k|G) − log P(k|R).

    P(k|G) = Poisson(k | λ_G · E).
    P(k|R) = ∫ Poisson(k | (λ_G + μ) · E) · LogNormal(μ; μ_R, σ_R²) dμ,
             evaluated by 21-node Gauss-HermiteE quadrature on log μ.
    """
    k = np.asarray(k_u, dtype=np.float64)
    E = np.asarray(E, dtype=np.float64)

    # lgamma(k+1) is constant across both classes and all 21 quadrature nodes;
    # compute once and reuse (≈10× cheaper via scipy.gammaln than frompyfunc,
    # and eliminates the 21-node Python for-loop entirely).
    lgamma_k1 = _vlgamma(k + 1.0)

    # G-class log-likelihood (vector over regions). In practice lam_G > 0 and
    # E > 0 for every soft region (mappable-bp floor is applied upstream in
    # _stats.compute_region_stats), so we can skip the safe-mask branch.
    rate_G = lam_G * E
    ll_G = k * np.log(np.maximum(rate_G, _EPS)) - rate_G - lgamma_k1

    # R-class: quadrature nodes on log μ → μ_j = exp(μ_R + σ_R · x_j).
    # Single broadcast over (K, N) avoids the per-node Python loop and the
    # associated (K,) copies from _poisson_logpmf.
    mu_j = np.exp(mu_R + sigma_R * _GH_NODES)  # (K,)
    rate_R = (lam_G + mu_j[:, None]) * E[None, :]  # (K, N)
    log_p_R_nodes = k[None, :] * np.log(rate_R) - rate_R - lgamma_k1[None, :]
    # logsumexp with weights
    log_p_R = _logsumexp(log_p_R_nodes + _GH_LOGW[:, None], axis=0)

    return ll_G - log_p_R


def _count_llr_poisson_ln_composite(
    k_u: np.ndarray,
    E_on: np.ndarray,
    E_off: np.ndarray,
    lam_G_on: float,
    lam_G_off: float,
    mu_R: float,
    sigma_R: float,
) -> np.ndarray:
    """Composite-Poisson per-region count LLR: log P(k|G) − log P(k|R).

    P(k|G) = Poisson(k | λ_on · E_on + λ_off · E_off).
    P(k|R) = ∫ Poisson(k | λ_on·E_on + λ_off·E_off + μ · (E_on+E_off))
             · LogNormal(μ; μ_R, σ_R²) dμ, evaluated by 21-node
    Gauss-HermiteE quadrature on log μ.

    When ``E_off`` is identically zero and ``lam_G_off`` is passed as
    anything (it drops out), this reduces to :func:`_count_llr_poisson_ln`
    with ``lam_G = lam_G_on`` and ``E = E_on`` — preserving the
    single-rate result.
    """
    k = np.asarray(k_u, dtype=np.float64)
    E_on = np.asarray(E_on, dtype=np.float64)
    E_off = np.asarray(E_off, dtype=np.float64)
    lgamma_k1 = _vlgamma(k + 1.0)

    # G-class: composite Poisson rate (sum of two independent
    # channels).  Floor with _EPS to avoid log(0) on the rare region
    # where both attributions are zero (e.g. fully-unmapped block);
    # such regions are filtered out upstream by the mappable_bp floor.
    rate_G = lam_G_on * E_on + lam_G_off * E_off
    ll_G = k * np.log(np.maximum(rate_G, _EPS)) - rate_G - lgamma_k1

    # R-class: μ_R prior on a single expressed-rate acting on total
    # mappable_bp (E_on + E_off).  Option (a): keep μ global /
    # unscaled across capture classes; the LogNormal σ_R absorbs any
    # residual probe enrichment in RNA counts.
    E_tot = E_on + E_off
    mu_j = np.exp(mu_R + sigma_R * _GH_NODES)  # (K,)
    rate_R = rate_G[None, :] + mu_j[:, None] * E_tot[None, :]  # (K, N)
    log_p_R_nodes = k[None, :] * np.log(np.maximum(rate_R, _EPS)) - rate_R - lgamma_k1[None, :]
    log_p_R = _logsumexp(log_p_R_nodes + _GH_LOGW[:, None], axis=0)

    return ll_G - log_p_R


# ---------------------------------------------------------------------------
# Strand LLR: Beta-Binomial (G) vs ρ-integrated Beta-Binomial (R)
# ---------------------------------------------------------------------------


def _compute_k_sense_valid(
    stats: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (valid_mask, k_sense[valid], n_unspliced[valid]).

    Valid ⇔ tx_strand != 0 and n_unspliced ≥ 1.  A single oriented
    fragment contributes a real, bounded BB log-likelihood; there is
    no theoretical reason to discard n = 1 regions.  ``k_sense`` is
    the R1-antisense-convention sense-read count.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    tx_strand = stats["tx_strand"]
    valid = (n_unspliced >= 1) & (tx_strand != 0)
    if not valid.any():
        return valid, np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
    k_sense = np.zeros(n_unspliced.shape[0], dtype=np.float64)
    k_sense[tx_strand == 1] = n_unspliced[tx_strand == 1] - n_pos[tx_strand == 1]
    k_sense[tx_strand == -1] = n_pos[tx_strand == -1]
    return valid, k_sense[valid].astype(np.float64), n_unspliced[valid].astype(np.float64)


# ---------------------------------------------------------------------------
# Beta-Binomial strand channel — sufficiency-grouped fast path
# ---------------------------------------------------------------------------
#
# The BB strand log-likelihood depends only on (k_sense, n_unspliced).
# Empirically ~97% of valid regions share a (k, n) pattern with another
# region — VCaP dna01m: 258K valid regions → 6.6K unique pairs.  Grouping
# turns every EM step from O(N) into O(U) where U ≪ N.


class _StrandBBCache:
    """Pre-computed invariants for Beta-Binomial strand log-likelihoods.

    Groups the N valid (k, n) pairs into U unique patterns with an
    ``inverse`` index back to the full N-vector.  All per-κ kernel
    work happens on the U-long arrays.  The full per-region output
    is recovered by a single gather ``out_u[inverse]``.
    """

    __slots__ = (
        "valid", "inverse", "uk", "un", "un_minus_uk", "log_binom",
        "p_rho", "one_minus_p_rho",
    )

    def __init__(self, stats: dict[str, np.ndarray]):
        self.valid, k_all, n_all = _compute_k_sense_valid(stats)
        if k_all.size == 0:
            self.inverse = np.empty(0, dtype=np.intp)
            self.uk = np.empty(0, dtype=np.float64)
            self.un = np.empty(0, dtype=np.float64)
            self.un_minus_uk = np.empty(0, dtype=np.float64)
            self.log_binom = np.empty(0, dtype=np.float64)
        else:
            # Combine (k, n) into one int64 key; n fits well under 2^31 on
            # any realistic RNA-seq region so this is exact.
            n_stride = int(n_all.max()) + 2
            key = k_all.astype(np.int64) * n_stride + n_all.astype(np.int64)
            _, first_idx, inverse = np.unique(
                key, return_index=True, return_inverse=True,
            )
            self.inverse = inverse.astype(np.intp, copy=False)
            self.uk = k_all[first_idx].astype(np.float64, copy=False)
            self.un = n_all[first_idx].astype(np.float64, copy=False)
            self.un_minus_uk = self.un - self.uk
            # log C(n, k) — invariant in κ and ρ; cache once.
            self.log_binom = (
                _vlgamma(self.un + 1.0)
                - _vlgamma(self.uk + 1.0)
                - _vlgamma(self.un_minus_uk + 1.0)
            )
        # ρ grid depends on SS; filled by ``set_ss``.
        self.p_rho = None
        self.one_minus_p_rho = None

    def set_ss(self, strand_specificity: float) -> None:
        ss = float(strand_specificity)
        p = _RHO_NODES * ss + (1.0 - _RHO_NODES) * 0.5
        self.p_rho = p
        self.one_minus_p_rho = 1.0 - p

    def aggregate_weights(self, w_valid: np.ndarray) -> np.ndarray:
        """Sum per-valid-region weights into per-unique-pair weights.

        ``w_valid`` has shape (N_valid,).  Result has shape (U,).
        """
        if self.uk.size == 0:
            return np.zeros(0, dtype=np.float64)
        return np.bincount(
            self.inverse, weights=w_valid, minlength=self.uk.shape[0],
        ).astype(np.float64, copy=False)


def _bb_g_loglik_unique(cache: _StrandBBCache, kappa_G: float) -> np.ndarray:
    """G-class BB log-likelihood on the U unique (k, n) pairs.

    Symmetric BB with α = β = κ_G / 2.  The log-Beta denominator
    collapses to a scalar ``2·lgamma(α) − lgamma(2α)``.
    """
    a = 0.5 * float(kappa_G)
    log_B_num = (
        _vlgamma(cache.uk + a)
        + _vlgamma(cache.un_minus_uk + a)
        - _vlgamma(cache.un + 2.0 * a)
    )
    log_B_den = 2.0 * math.lgamma(a) - math.lgamma(2.0 * a)
    return cache.log_binom + log_B_num - log_B_den


def _bb_r_loglik_unique(cache: _StrandBBCache, kappa_R: float) -> np.ndarray:
    """R-class ρ-integrated BB log-likelihood on the U unique pairs."""
    if cache.p_rho is None:
        raise RuntimeError("_StrandBBCache.set_ss() must be called first")
    kR = float(kappa_R)
    alpha_j = (kR * cache.p_rho)[:, None]              # (K_rho, 1)
    beta_j = (kR * cache.one_minus_p_rho)[:, None]     # (K_rho, 1)
    # lgamma(n + κ) is ρ-independent: compute as (U,) once and broadcast.
    lgamma_n_plus_k = _vlgamma(cache.un + kR)          # (U,)
    log_B_num = (
        _vlgamma(cache.uk[None, :] + alpha_j)
        + _vlgamma(cache.un_minus_uk[None, :] + beta_j)
        - lgamma_n_plus_k[None, :]
    )                                                  # (K_rho, U)
    log_B_den = (
        _vlgamma(kR * cache.p_rho)
        + _vlgamma(kR * cache.one_minus_p_rho)
        - math.lgamma(kR)
    )                                                  # (K_rho,)
    log_p_nodes = cache.log_binom[None, :] + log_B_num - log_B_den[:, None]
    return _logsumexp(log_p_nodes + _RHO_LOGW[:, None], axis=0)


def _bb_logpmf(
    k: np.ndarray,
    n: np.ndarray,
    alpha: np.ndarray | float,
    beta: np.ndarray | float,
) -> np.ndarray:
    """Log PMF of BetaBinomial(k; n, α, β).  Reference implementation.

    Kept for tests and external callers.  The EM hot path uses the
    grouped-unique-pairs code above which is ~30-100× cheaper.
    """
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    log_binom = _vlgamma(n + 1.0) - _vlgamma(k + 1.0) - _vlgamma(n - k + 1.0)
    log_B_num = _vlgamma(k + alpha) + _vlgamma(n - k + beta) - _vlgamma(n + alpha + beta)
    log_B_den = _vlgamma(alpha) + _vlgamma(beta) - _vlgamma(alpha + beta)
    return log_binom + log_B_num - log_B_den


def _strand_llr_bb_rho(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    kappa_G: float,
    kappa_R: float,
    cache: "_StrandBBCache | None" = None,
) -> np.ndarray:
    """Per-region strand LLR = log P(k|n,G,κ_G) − log P(k|n,R,κ_R).

    Zero for invalid regions (n=0 or tx_strand=0).  Computes the LLR
    once on the U unique (k, n) pairs, then scatters to the N valid
    regions via ``cache.inverse``.  At SS → 0.5 the R integrand
    collapses onto G and the LLR vanishes automatically.
    """
    n_regions = stats["n_unspliced"].shape[0]
    llr = np.zeros(n_regions, dtype=np.float64)
    if cache is None:
        cache = _StrandBBCache(stats)
        cache.set_ss(strand_specificity)
    if cache.uk.size == 0:
        return llr
    ll_g_u = _bb_g_loglik_unique(cache, kappa_G)
    ll_r_u = _bb_r_loglik_unique(cache, kappa_R)
    llr_u = ll_g_u - ll_r_u
    # Gather to full valid-region vector, then place in the N-long output.
    llr[cache.valid] = llr_u[cache.inverse]
    return llr


def _strand_posterior(llr: np.ndarray, pi_0: float = 0.5) -> np.ndarray:
    """Strand-only γ̃_i = sigmoid(logit(π₀) + LLR_strand_i).  Diagnostic."""
    pi_safe = float(np.clip(pi_0, _EPS, 1.0 - _EPS))
    logit = math.log(pi_safe / (1.0 - pi_safe))
    return 1.0 / (1.0 + np.exp(-np.clip(logit + llr, -500.0, 500.0)))


def _golden_section_max_log(
    obj,
    a: float = _KAPPA_LO,
    b: float = _KAPPA_HI,
    log_tol: float = 0.03,
    max_iter: int = 40,
) -> float:
    """Argmax of ``κ ↦ obj(κ)`` on [a, b] searched in log-κ space.

    κ spans four orders of magnitude (0.1 → 1000).  A log-space
    golden-section search with relative tolerance ~3% needs ~9
    evaluations and gives more than enough precision for the E-step
    posterior (BB log-likelihood is slowly varying in κ near the
    optimum).
    """
    la, lb = math.log(a), math.log(b)
    gr = (math.sqrt(5.0) + 1.0) / 2.0
    lc = lb - (lb - la) / gr
    ld = la + (lb - la) / gr
    fc, fd = obj(math.exp(lc)), obj(math.exp(ld))
    for _ in range(max_iter):
        if abs(lb - la) < log_tol:
            break
        if fc > fd:
            lb, ld, fd = ld, lc, fc
            lc = lb - (lb - la) / gr
            fc = obj(math.exp(lc))
        else:
            la, lc, fc = lc, ld, fd
            ld = la + (lb - la) / gr
            fd = obj(math.exp(ld))
    return math.exp(0.5 * (la + lb))


def _estimate_kappa_G(
    cache: _StrandBBCache,
    gamma_valid: np.ndarray,
    min_eff_regions: float = 2.0,
    kappa_lo: float = _KAPPA_LO,
) -> float:
    """γ-weighted MLE of κ_G on the G-class symmetric Beta-Binomial.

    Operates on the U unique (k, n) pairs with aggregated weights.
    Falls back to the Binomial limit when effective G-weight is too
    small to constrain overdispersion.

    ``kappa_lo`` clips the search lower bound; caller typically sets
    it to ``CalibrationConfig.kappa_gdna_min`` (default 3.0) to keep the
    G-class strictly unimodal (α = β = κ/2 > 1) and prevent the
    U-shape pathology at low gDNA loads.
    """
    if cache.uk.size < 2:
        return _KAPPA_HI
    if float(gamma_valid.sum()) < min_eff_regions:
        return _KAPPA_HI
    w_u = cache.aggregate_weights(gamma_valid)
    return _golden_section_max_log(
        lambda kG: float(np.dot(w_u, _bb_g_loglik_unique(cache, kG))),
        a=kappa_lo,
    )


def _estimate_kappa_R(
    cache: _StrandBBCache,
    gamma_valid: np.ndarray,
    min_eff_regions: float = 2.0,
) -> float:
    """(1-γ)-weighted MLE of κ_R on the ρ-integrated R-class BB.

    Operates on the U unique (k, n) pairs.
    """
    if cache.uk.size < 2:
        return _KAPPA_HI
    w = 1.0 - gamma_valid
    if float(w.sum()) < min_eff_regions:
        return _KAPPA_HI
    w_u = cache.aggregate_weights(w)
    return _golden_section_max_log(
        lambda kR: float(np.dot(w_u, _bb_r_loglik_unique(cache, kR))),
    )


# ---------------------------------------------------------------------------
# Fragment-length LLR
# ---------------------------------------------------------------------------


def _fl_llr(
    fl_region_ids: np.ndarray,
    fl_frag_lens: np.ndarray,
    gdna_fl_model: FragmentLengthModel,
    rna_fl_model: FragmentLengthModel,
    n_regions: int,
) -> np.ndarray:
    """Per-region Σ log[q_G(f) / q_R(f)] over fragments in the region.

    Both histograms are normalised and smoothed with symmetric
    Dirichlet (α = 1/M per bin) so only shape drives the LLR.
    """
    llr = np.zeros(n_regions, dtype=np.float64)
    if len(fl_region_ids) == 0:
        return llr
    if gdna_fl_model.total_weight <= 0 or rna_fl_model.total_weight <= 0:
        return llr

    max_size = gdna_fl_model.max_size
    M = max_size + 1
    T_g = max(float(gdna_fl_model.total_weight), 1.0)
    T_r = max(float(rna_fl_model.total_weight), 1.0)
    alpha = 1.0 / M
    q_g = gdna_fl_model.counts[:M] / T_g + alpha
    q_r = rna_fl_model.counts[:M] / T_r + alpha
    log_ratio = np.log(q_g) - np.log(q_r)

    rid = np.asarray(fl_region_ids, dtype=np.intp)
    flen = np.asarray(fl_frag_lens, dtype=np.intp)
    # Bounds validated once up-front in run_em; see bounds-check block there.
    # Per-iter loop avoids re-masking 15M+ items.
    contrib = log_ratio[flen]
    # Scatter-add per region.
    np.add.at(llr, rid, contrib)
    return llr


# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------


def _robust_initial_lam_G(
    n_u_soft: np.ndarray,
    E_soft: np.ndarray,
) -> float:
    """Pooled Poisson rate over soft eligible regions.

    This is the one-shot MLE

        λ̂ = Σ k^u_i / Σ E_i    (over soft regions)

    which is guaranteed positive whenever any soft region carries a
    fragment.  It is biased **high** by expressed contamination — the
    soft pool contains low-μ expressed regions that escaped the
    hard-splice anchor — but EM's first E-step uses this as a starting
    point and the count LLR then pulls λ̂ down as expressed regions
    migrate toward γ=0.  The only failure mode is "all soft counts are
    zero", in which case λ̂ = 0 is the correct answer.

    Chosen over a median/percentile split because the median of
    ``k^u / E`` is identically 0 whenever at least half the regions
    have zero counts (a common case when λ_G·Ē < 1), which collapses
    the numerator to 0 and kills the EM.  The pooled mean is
    free of that failure mode and has no tunable threshold.
    """
    if n_u_soft.size == 0:
        return 0.0
    num = float(np.asarray(n_u_soft, dtype=np.float64).sum())
    den = float(np.asarray(E_soft, dtype=np.float64).sum())
    return num / den if den > _EPS else 0.0


def _seed_gamma(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
) -> tuple[np.ndarray, float, int]:
    """Anchor-only initialisation.

    * Hard anchor from splice evidence: (k^s > 0) & eligible → γ = 0.
      Spliced evidence is *definitionally* RNA — no ambiguity.
    * All other (eligible or ineligible) regions start at γ = π_init
      = 0.5, i.e. a truly uninformative prior.  The first E-step then
      partitions them via the count, strand, and FL LLRs.

    Notes
    -----
    This is deliberately **symmetry-breaking from one side only**.
    The old percentile-based gDNA seed injected a magic threshold and
    mis-classified low-μ expressed regions with rate ≈ λ_G in small
    samples.  We rely on the count-LLR E-step to discover the gDNA
    class from its Poisson signature rather than planting it by fiat.
    """
    n = stats["n_unspliced"].shape[0]
    n_s = stats["n_spliced"]

    expressed_seed = (n_s > 0) & eligible
    pi_init = 0.5
    gamma = np.full(n, pi_init, dtype=np.float64)
    gamma[expressed_seed] = 0.0

    return gamma, pi_init, int(expressed_seed.sum())


# ---------------------------------------------------------------------------
# E-step
# ---------------------------------------------------------------------------


def _e_step(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
    lam_G_on: float,
    lam_G_off: float,
    mu_R: float,
    sigma_R: float,
    pi_soft: float,
    strand_specificity: float,
    kappa_G: float,
    kappa_R: float,
    llr_fl: np.ndarray | None,
    E_on: np.ndarray,
    E_off: np.ndarray,
    cache: "_StrandBBCache | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute γ_i = P(z=G | data_i) and the per-region strand LLR.

    Hard:
      * k^s > 0 → γ = 0
      * ineligible → γ = π_soft

    Count channel: composite-Poisson (``lam_G_on · E_on + lam_G_off ·
    E_off``) under G; single global LogNormal(μ_R, σ_R²) prior over
    an additive expressed rate under R.  When no capture BED was
    supplied, ``E_on`` is identically zero and the composite rate
    collapses to ``lam_G_off · E_off = lam_G · mappable_bp`` —
    bit-identical to the pre-composite single-rate pathway.

    Returns (gamma, llr_strand_full) where ``llr_strand_full`` covers
    all regions (zero on invalid rows) — kept for diagnostics and the
    final γ̃^strand computation.
    """
    n = stats["n_unspliced"].shape[0]
    n_s = stats["n_spliced"]
    gamma = np.full(n, pi_soft, dtype=np.float64)

    hard_expr = (n_s > 0) & eligible
    gamma[hard_expr] = 0.0

    llr_strand_full = _strand_llr_bb_rho(
        stats, strand_specificity, kappa_G, kappa_R, cache=cache,
    )

    soft = eligible & ~hard_expr
    if not soft.any():
        return gamma, llr_strand_full

    pi_safe = float(np.clip(pi_soft, _EPS, 1.0 - _EPS))
    log_prior_odds = math.log(pi_safe / (1.0 - pi_safe))

    if E_on is None or not bool(np.any(E_on[soft] > 0)):
        # Non-capture fast path: ``E_on ≡ 0`` → the composite model
        # reduces algebraically to the classical single-rate form, but
        # the factored ``rate_R = (λ + μ) · E`` computation preserves
        # bit-identical FP with the pre-composite code path.
        llr_count = _count_llr_poisson_ln(
            stats["n_unspliced"][soft],
            stats["mappable_bp"][soft],
            lam_G_off,
            mu_R,
            sigma_R,
        )
    else:
        llr_count = _count_llr_poisson_ln_composite(
            stats["n_unspliced"][soft],
            E_on[soft],
            E_off[soft],
            lam_G_on,
            lam_G_off,
            mu_R,
            sigma_R,
        )
    llr_strand = llr_strand_full[soft]

    if llr_fl is not None:
        llr_fl_soft = llr_fl[soft]
    else:
        llr_fl_soft = 0.0

    log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft
    log_odds = np.clip(log_odds, -500.0, 500.0)
    gamma[soft] = 1.0 / (1.0 + np.exp(-log_odds))
    return gamma, llr_strand_full


# ---------------------------------------------------------------------------
# M-step
# ---------------------------------------------------------------------------


def _m_step(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
    gamma: np.ndarray,
    E_on: np.ndarray,
    E_off: np.ndarray,
    lam_G_on_prev: float,
    lam_G_off_prev: float,
    *,
    capture_class_mode: bool,
    inner_iters: int = 2,
) -> tuple[float, float, float, float, float, float]:
    """Update (λ_G_on, λ_G_off, μ_R, σ_R, π, π_soft) from posteriors.

    Composite-Poisson M-step:

    * In capture-class mode (``capture_class_mode=True``) the on/off
      rates are updated by a 2-iteration nested EM that analytically
      unmixes the two independent Poisson channels.  The Poisson
      superposition theorem gives

          E[k_on_i | k_i, λ_on, λ_off] = k_i · (λ_on · E_on_i) / rate_G_i,

      and the γ-weighted M-step is
      ``λ_on = Σ γ_i · k̂_on_i / Σ γ_i · E_on_i`` (symmetric off).
      Two inner iterations suffice — the objective is strictly
      convex in (λ_on, λ_off) and the fixed point is well inside
      the convergence basin after one inner pass.

    * Without capture-class (``E_on`` identically zero), the update
      collapses to the classical closed-form
      ``λ = Σ γ·k / Σ γ·E_off`` with ``λ_on = λ_off`` — bit-identical
      to the pre-composite single-rate code path.

    * (μ_R, σ_R²) use the anchor-only (hard-spliced) update operating
      on the per-region effective rate
      ``rate_G_i / (E_on_i + E_off_i)``.  Anchor regions are
      definitively RNA; their raw per-bp rate minus the G-channel
      offset gives the expressed μ̂_i.

    * π, π_soft use the usual γ-means over eligible regions.
    """
    n_u = stats["n_unspliced"]
    n_s = stats["n_spliced"]

    if not eligible.any():
        return 0.0, 0.0, 0.0, 1.0, 0.5, 0.5

    g = gamma[eligible]
    ku = n_u[eligible].astype(np.float64)
    ns = n_s[eligible]
    eon = E_on[eligible]
    eoff = E_off[eligible]
    E_tot = eon + eoff

    # -- λ_G_on, λ_G_off: nested inner EM over the Poisson channels --
    if capture_class_mode:
        lam_on = float(lam_G_on_prev)
        lam_off = float(lam_G_off_prev)
        # Seed from the coverage-weighted pooled rate if the previous
        # iteration supplied zeros (iter 1 before any gDNA was
        # attributed).  Keeps the inner update away from the 0/0
        # fixed-point of both rates being zero.
        if lam_on <= 0.0 and lam_off <= 0.0:
            den = float(np.sum(g * E_tot))
            lam_seed = float(np.sum(g * ku)) / den if den > _EPS else 0.0
            lam_on = lam_seed
            lam_off = lam_seed
        for _ in range(int(inner_iters)):
            # Inner E-step: unmix k_i into (k_on, k_off) under the
            # current (λ_on, λ_off).  The 1e-12 floor in the
            # denominator guards against regions where both rate
            # components are exactly zero (e.g. fully-off-target
            # region with λ_off driven to zero transiently).
            rate_on = lam_on * eon
            rate_off = lam_off * eoff
            denom = rate_on + rate_off + 1e-12
            k_on_hat = g * ku * (rate_on / denom)
            k_off_hat = g * ku - k_on_hat
            # Inner M-step: γ-weighted Poisson MLE per channel.
            den_on = float(np.sum(g * eon))
            den_off = float(np.sum(g * eoff))
            lam_on = float(np.sum(k_on_hat)) / den_on if den_on > _EPS else 0.0
            lam_off = float(np.sum(k_off_hat)) / den_off if den_off > _EPS else 0.0
    else:
        # Non-capture path: eon is identically zero; rate_G =
        # lam_off · eoff; single closed-form update.  Use the
        # original-dtype ``n_u`` / ``mappable_bp`` arrays (typically
        # float32) for bit-identical FP with the pre-composite single-
        # rate code path.
        _ku_raw = n_u[eligible]
        _E_raw = stats["mappable_bp"][eligible]
        num = float(np.sum(g * _ku_raw))
        den = float(np.sum(g * _E_raw))
        lam_off = num / den if den > _EPS else 0.0
        lam_on = lam_off  # keep fields consistent; unused downstream

    # Per-region effective G rate for the anchor μ_R update.  In
    # non-capture mode this is exactly ``lam_off`` (eon ≡ 0); broadcast
    # the scalar to preserve bit-identical FP with the pre-composite
    # single-rate code path.
    if capture_class_mode:
        with np.errstate(divide="ignore", invalid="ignore"):
            lam_eff = np.where(
                E_tot > 0,
                (lam_on * eon + lam_off * eoff) / E_tot,
                0.0,
            )
    else:
        lam_eff = np.full_like(E_tot, lam_off)

    # -- μ_R, σ_R: anchor-only update on the eligible subset --
    anchor_local = ns > 0
    if anchor_local.any():
        if capture_class_mode:
            ku_a = ku[anchor_local]
            E_a = E_tot[anchor_local]
            lam_a = lam_eff[anchor_local]
            rate_a = np.where(E_a > 0, ku_a / E_a, 0.0)
        else:
            # Preserve pre-composite FP: work in the original stats
            # dtypes (typically float32) for bit-identical ``rate_a``
            # and therefore identical ``μ_R``, ``σ_R``.  Under the
            # non-capture path ``E_tot ≡ mappable_bp`` and
            # ``lam_eff ≡ lam_off``.
            anchor_mask_full = eligible.copy()
            anchor_mask_full[eligible] = anchor_local
            ku_a = n_u[anchor_mask_full]
            ee_a = stats["mappable_bp"][anchor_mask_full]
            rate_a = np.where(ee_a > 0, ku_a / ee_a, 0.0)
            lam_a = lam_off
            E_a = ee_a
        med_E = float(np.median(E_a[E_a > 0])) if np.any(E_a > 0) else 1.0
        eps_mu = 1.0 / max(med_E, 1.0)
        mu_hat = np.maximum(rate_a - lam_a, eps_mu)
        log_mu = np.log(mu_hat)
        if log_mu.size >= 2:
            mu_R = float(np.mean(log_mu))
            var_R = float(np.var(log_mu))
            sigma_R = math.sqrt(max(var_R, 1e-4))
        else:
            mu_R = float(log_mu.item())
            sigma_R = 1.5
    else:
        lam_scalar = lam_off if lam_off > 0 else max(lam_on, 1e-12)
        mu_R = math.log(max(10.0 * lam_scalar, 1e-6))
        sigma_R = 2.0

    # -- Mixing proportions --
    pi = float(np.mean(g))
    soft_mask = ns == 0
    if soft_mask.any():
        pi_soft = float(np.mean(g[soft_mask]))
    else:
        pi_soft = pi
    pi = float(np.clip(pi, _EPS, 1.0 - _EPS))
    pi_soft = float(np.clip(pi_soft, _EPS, 1.0 - _EPS))

    return lam_on, lam_off, mu_R, sigma_R, pi, pi_soft


# ---------------------------------------------------------------------------
# FL model construction helper (for per-iteration RNA/gDNA FL models)
# ---------------------------------------------------------------------------


def _build_fl_histogram(
    fl_region_ids: np.ndarray,
    fl_frag_lens: np.ndarray,
    weights: np.ndarray,
    n_regions: int,
    max_size: int = 1000,
) -> FragmentLengthModel:
    """Build a plain weighted FL histogram (no prior blending)."""
    model = FragmentLengthModel(max_size=max_size)
    if len(fl_region_ids) == 0:
        model.finalize()
        return model
    rid = np.asarray(fl_region_ids, dtype=np.intp)
    flen = np.asarray(fl_frag_lens, dtype=np.intp)
    # Upstream (EM driver) guarantees rid ∈ [0, n_regions) and
    # flen ∈ [0, max_size]; we only need to drop the flen==0 bin (no signal
    # there) and let bincount handle everything else. Zero-weight entries are
    # valid for np.bincount and cost less than explicit masking on 15M fragments.
    valid = flen >= 1
    if not valid.any():
        model.finalize()
        return model
    rid_v = rid[valid]
    flen_v = flen[valid]
    w_v = weights[rid_v]
    counts = np.bincount(flen_v, weights=w_v, minlength=max_size + 1).astype(np.float64)
    model.counts[: max_size + 1] = counts[: max_size + 1]
    model._total_weight = float(w_v.sum())
    model.finalize()
    return model


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run_em(
    stats: dict[str, np.ndarray],
    fl_table: pd.DataFrame | None,
    strand_specificity: float,
    *,
    mean_frag_len: float,
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    kappa_gdna_min: float = 3.0,
    diagnostics: bool = False,
    capture_e: tuple[np.ndarray, np.ndarray] | None = None,
) -> EMFit:
    """Run the two-component mixture EM and return an :class:`EMFit`.

    Parameters
    ----------
    stats
        Output of :func:`rigel.calibration._stats.compute_region_stats`.
        Must contain ``n_unspliced``, ``n_spliced``, ``n_pos``,
        ``tx_strand``, ``mappable_bp``, and ``ref``.
    fl_table
        Per-fragment (region_id, frag_len) table; may be ``None`` or
        empty.  When available, the FL LLR channel is enabled from
        iteration 1.
    strand_specificity
        Library SS ∈ [0.5, 1.0].  Feeds the R-class ρ-integrated
        Beta-Binomial; at SS = 0.5 the strand LLR vanishes identically
        without any explicit gating.
    mean_frag_len
        Library mean fragment length — eligibility floor for
        ``mappable_bp_i``.
    kappa_gdna_min
        Lower bound on the G-class Beta-Binomial dispersion.  Default
        3.0 keeps α = β = κ_G/2 > 1 — the G-class stays unimodal and
        peaked at 0.5, the correct shape for real gDNA.  Values below
        2.0 admit the biologically-nonsense U-shape that absorbs
        authentic-RNA tails at low contamination.
    capture_e
        Optional ``(E_on, E_off)`` tuple from
        :func:`rigel.calibration._annotate.annotate_capture_class`.
        When provided, the gDNA density channel fits two rates
        (``λ_on``, ``λ_off``) via a composite-Poisson superposition
        model with per-bp attribution.  When ``None``, a single global
        rate is fit over ``mappable_bp`` (classical behavior).
    """
    n_u = stats["n_unspliced"]
    n_s = stats["n_spliced"]
    E = stats["mappable_bp"]
    n_regions = n_u.shape[0]

    # Composite-Poisson capture-class attribution.  When no BED is
    # supplied we emulate the single-rate pathway by setting
    # ``E_on ≡ 0`` and ``E_off ≡ mappable_bp`` — the composite count
    # LLR and M-step then reduce to the classical single-rate update
    # for ``lam_G_off`` (and ``lam_G_on`` is inert).
    capture_class_mode = capture_e is not None
    if capture_e is None:
        E_on = np.zeros(n_regions, dtype=np.float64)
        E_off = E.astype(np.float64, copy=False)
    else:
        _e_on, _e_off = capture_e
        E_on = np.asarray(_e_on, dtype=np.float64)
        E_off = np.asarray(_e_off, dtype=np.float64)
        if E_on.shape != (n_regions,) or E_off.shape != (n_regions,):
            raise ValueError(
                f"capture_e shapes {E_on.shape}/{E_off.shape} do not match "
                f"n_regions {n_regions}"
            )

    # Eligibility: region must carry ≥ one fragment's worth of mappable bp.
    mappable_floor = max(1.0, float(mean_frag_len))
    eligible = E >= mappable_floor

    if not eligible.any():
        logger.warning(
            "EM calibration: no eligible regions (mappable_bp < %.1f everywhere)",
            mappable_floor,
        )
        return EMFit(
            lam_G=0.0,
            mu_R=0.0,
            sigma_R=1.0,
            pi=1.0,
            pi_soft=1.0,
            gamma=np.ones(n_regions, dtype=np.float64),
            region_gamma_strand=np.full(n_regions, 0.5, dtype=np.float64),
            kappa_G=0.0,
            kappa_R=0.0,
            n_iter=0,
            converged=True,
            n_eligible=0,
            n_soft=0,
            n_spliced_hard=0,
            lam_G_on=0.0,
            lam_G_off=0.0,
            capture_class_mode=capture_class_mode,
        )

    fl_region_ids = np.empty(0, dtype=np.intp)
    fl_frag_lens = np.empty(0, dtype=np.intp)
    if fl_table is not None and len(fl_table) > 0:
        fl_region_ids = fl_table["region_id"].to_numpy(dtype=np.intp)
        fl_frag_lens = fl_table["frag_len"].to_numpy(dtype=np.intp)
        _FL_MAX_SIZE = 1000
        _ok = (
            (fl_region_ids >= 0)
            & (fl_region_ids < n_regions)
            & (fl_frag_lens >= 0)
            & (fl_frag_lens <= _FL_MAX_SIZE)
        )
        if not _ok.all():
            fl_region_ids = fl_region_ids[_ok]
            fl_frag_lens = fl_frag_lens[_ok]

    # --- Anchor-only seed ---
    gamma, pi_init, n_spliced_hard = _seed_gamma(stats, eligible)

    # Initial λ_G: pooled rate over soft eligible regions.
    soft_mask = eligible & (n_s == 0)
    if soft_mask.any():
        lam_G_init = _robust_initial_lam_G(n_u[soft_mask], E[soft_mask])
    else:
        lam_G_init = 0.0

    # Seed both capture-class rates with the same pooled init; they
    # diverge after the first composite M-step.  In single-rate mode
    # ``lam_G_on`` is inert (E_on ≡ 0) so only ``lam_G_off`` drives
    # the subsequent E-step.
    lam_G_on = lam_G_init
    lam_G_off = lam_G_init

    # Initial μ_R, σ_R from spliced-seed regions (definitively RNA).
    expressed_seed_mask = (n_s > 0) & eligible
    if expressed_seed_mask.any():
        ee = E[expressed_seed_mask]
        rates = np.where(ee > 0, n_u[expressed_seed_mask] / ee, 0.0)
        med_E = float(np.median(ee[ee > 0])) if np.any(ee > 0) else 1.0
        eps_mu = 1.0 / max(med_E, 1.0)
        mu_hat = np.maximum(rates - lam_G_init, eps_mu)
        log_mu = np.log(mu_hat)
        mu_R = float(np.mean(log_mu))
        sigma_R = max(float(np.std(log_mu)), 0.5)
    else:
        mu_R = math.log(max(10.0 * lam_G_init, 1e-6))
        sigma_R = 2.0

    # --- Initial (κ_G, κ_R): Binomial limit.  γ-weighted M-step
    # adapts κ as the classification firms up.
    kappa_G = _KAPPA_HI
    kappa_R = _KAPPA_HI

    n_soft = int(soft_mask.sum())
    n_elig = int(eligible.sum())
    pi = pi_init
    pi_soft = pi_init

    if capture_class_mode:
        n_regions_on = int((E_on > 0).sum())
        logger.info(
            "EM calibration init (composite-Poisson, capture-class): "
            "n_elig=%d, n_soft=%d, n_spliced_hard=%d, π=%.3f, π_soft=%.3f, "
            "λ_G(seed)=%.3e, μ_R=%.2f, σ_R=%.2f, κ_G=%.2f, κ_R=%.2f  "
            "(ΣE_on/ΣE=%.3f, n_regions_on=%d/%d, SS=%.3f)",
            n_elig, n_soft, n_spliced_hard,
            pi, pi_soft, lam_G_init, mu_R, sigma_R, kappa_G, kappa_R,
            float(E_on.sum()) / max(float(E_on.sum() + E_off.sum()), 1.0),
            n_regions_on, n_regions, float(strand_specificity),
        )
    else:
        logger.info(
            "EM calibration init: n_elig=%d, n_soft=%d, n_spliced_hard=%d, "
            "π=%.3f, π_soft=%.3f, λ_G=%.3e, μ_R=%.2f, σ_R=%.2f, "
            "κ_G=%.2f, κ_R=%.2f  (SS=%.3f)",
            n_elig, n_soft, n_spliced_hard,
            pi, pi_soft, lam_G_init, mu_R, sigma_R,
            kappa_G, kappa_R, float(strand_specificity),
        )

    history: list[dict] = []
    converged = False
    n_iter = 0
    gdna_fl_model = None

    prev_pi_soft = pi_soft
    llr_strand_full = np.zeros(n_regions, dtype=np.float64)

    bb_cache = _StrandBBCache(stats)
    bb_cache.set_ss(strand_specificity)

    for iteration in range(max_iterations):
        n_iter = iteration + 1

        # --- FL LLR channel (after the first iteration) ---
        llr_fl = None
        if (
            gdna_fl_model is not None
            and fl_region_ids.size > 0
            and gdna_fl_model.total_weight > 0
        ):
            rna_fl_model = _build_fl_histogram(
                fl_region_ids, fl_frag_lens, 1.0 - gamma, n_regions,
            )
            if rna_fl_model.total_weight > 0:
                llr_fl = _fl_llr(
                    fl_region_ids, fl_frag_lens,
                    gdna_fl_model, rna_fl_model, n_regions,
                )

        # --- E-step ---
        gamma, llr_strand_full = _e_step(
            stats, eligible, lam_G_on, lam_G_off, mu_R, sigma_R,
            pi_soft, strand_specificity, kappa_G, kappa_R, llr_fl,
            E_on, E_off, cache=bb_cache,
        )

        # --- M-step (composite-Poisson nested EM) ---
        (lam_G_on, lam_G_off, mu_R, sigma_R, pi, pi_soft) = _m_step(
            stats, eligible, gamma, E_on, E_off,
            lam_G_on, lam_G_off,
            capture_class_mode=capture_class_mode,
        )

        # --- κ updates (γ-weighted MLEs on unique (k, n) pairs) ---
        gamma_valid = gamma[bb_cache.valid]
        kappa_G = _estimate_kappa_G(
            bb_cache, gamma_valid, kappa_lo=float(kappa_gdna_min),
        )
        kappa_R = _estimate_kappa_R(bb_cache, gamma_valid)

        # --- FL model update (γ-weighted gDNA FL) ---
        if fl_region_ids.size > 0:
            gdna_fl_model = _build_fl_histogram(
                fl_region_ids, fl_frag_lens, gamma, n_regions,
            )

        delta = abs(pi_soft - prev_pi_soft)
        if diagnostics:
            hist_entry = {
                "iter": n_iter,
                "lam_G_on": lam_G_on,
                "lam_G_off": lam_G_off,
                "mu_R": mu_R,
                "sigma_R": sigma_R,
                "pi": pi,
                "pi_soft": pi_soft,
                "kappa_G": kappa_G,
                "kappa_R": kappa_R,
                "delta_pi_soft": delta,
                "mean_gamma_soft": float(np.mean(gamma[eligible & (n_s == 0)]))
                    if (eligible & (n_s == 0)).any() else float("nan"),
            }
            history.append(hist_entry)
        if capture_class_mode:
            logger.debug(
                "EM iter %d: λ_on/λ_off=%.3e/%.3e μ_R=%.2f σ_R=%.2f π=%.3f "
                "π_soft=%.3f κ_G=%.2f κ_R=%.2f Δ=%.2e",
                n_iter, lam_G_on, lam_G_off, mu_R, sigma_R, pi, pi_soft,
                kappa_G, kappa_R, delta,
            )
        else:
            logger.debug(
                "EM iter %d: λ_G=%.3e μ_R=%.2f σ_R=%.2f π=%.3f π_soft=%.3f "
                "κ_G=%.2f κ_R=%.2f Δ=%.2e",
                n_iter, lam_G_off, mu_R, sigma_R, pi, pi_soft,
                kappa_G, kappa_R, delta,
            )

        if delta < convergence_tol:
            converged = True
            break
        prev_pi_soft = pi_soft

    # Effective / legacy single-rate λ_G: coverage-weighted average
    # of the two channels.  In non-capture mode this equals
    # ``lam_G_off`` exactly (E_on ≡ 0) — short-circuit to avoid the
    # 1-ULP drift that ``(0 + lam_off*E)/E`` can introduce.
    if not capture_class_mode:
        lam_G_eff = float(lam_G_off)
    else:
        sum_on = float(E_on.sum())
        sum_off = float(E_off.sum())
        sum_tot = sum_on + sum_off
        lam_G_eff = (
            (lam_G_on * sum_on + lam_G_off * sum_off) / sum_tot
            if sum_tot > _EPS
            else float(lam_G_off)
        )

    if capture_class_mode:
        enrichment = lam_G_on / max(lam_G_off, _EPS)
        logger.info(
            "EM calibration %s in %d iters (composite-Poisson): "
            "λ_on=%.3e, λ_off=%.3e (λ_on/λ_off=%.1fx), λ_eff=%.3e, "
            "μ_R=%.2f, σ_R=%.2f, π=%.3f, κ_G=%.2f, κ_R=%.2f",
            "converged" if converged else "did NOT converge",
            n_iter, lam_G_on, lam_G_off, enrichment, lam_G_eff,
            mu_R, sigma_R, pi, kappa_G, kappa_R,
        )
    else:
        logger.info(
            "EM calibration %s in %d iters: λ_G=%.3e, μ_R=%.2f, σ_R=%.2f, π=%.3f, "
            "κ_G=%.2f, κ_R=%.2f",
            "converged" if converged else "did NOT converge",
            n_iter, lam_G_eff, mu_R, sigma_R, pi, kappa_G, kappa_R,
        )

    # Strand-only diagnostic posterior (what does strand alone say?).
    region_gamma_strand = _strand_posterior(llr_strand_full, pi_0=0.5)

    return EMFit(
        lam_G=float(lam_G_eff),
        mu_R=float(mu_R),
        sigma_R=float(sigma_R),
        pi=float(pi),
        pi_soft=float(pi_soft),
        gamma=gamma,
        region_gamma_strand=region_gamma_strand,
        kappa_G=float(kappa_G),
        kappa_R=float(kappa_R),
        n_iter=n_iter,
        converged=converged,
        n_eligible=n_elig,
        n_soft=n_soft,
        n_spliced_hard=n_spliced_hard,
        history=history,
        lam_G_on=float(lam_G_on),
        lam_G_off=float(lam_G_off),
        capture_class_mode=capture_class_mode,
    )
