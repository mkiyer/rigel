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
# halves the dominant ``_bb_r_loglik_cached`` cost.
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
    """EM fit results for the two-component gDNA mixture."""

    # Global gDNA rate (fragments per base)
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


def _bb_logpmf(
    k: np.ndarray,
    n: np.ndarray,
    alpha: np.ndarray | float,
    beta: np.ndarray | float,
) -> np.ndarray:
    """Log PMF of BetaBinomial(k; n, α, β), vectorized and numerically safe.

    Uses the ``log C(n,k) + log B(k+α, n-k+β) - log B(α, β)`` form, which
    has no ``gammaln(0)`` traps at k=0 or k=n (the log-Beta terms stay
    finite as long as α, β > 0).  All inputs broadcast.

    Kept as the canonical reference; the hot path in the EM now routes
    through ``_StrandBBCache`` which pre-computes the k/n-only terms
    (``log_binom``) once per EM run, avoiding ~3 of the 6 gammaln calls
    per evaluation.
    """
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    log_binom = _vlgamma(n + 1.0) - _vlgamma(k + 1.0) - _vlgamma(n - k + 1.0)
    log_B_num = _vlgamma(k + alpha) + _vlgamma(n - k + beta) - _vlgamma(n + alpha + beta)
    log_B_den = _vlgamma(alpha) + _vlgamma(beta) - _vlgamma(alpha + beta)
    return log_binom + log_B_num - log_B_den


class _StrandBBCache:
    """Precomputed quantities for Beta-Binomial strand log-likelihoods.

    The EM calls ``_bb_g_loglik`` / ``_bb_r_loglik`` many times per
    iteration (golden-section search over κ), but the region-level
    arrays (k, n) and the binomial coefficient ``log C(n,k)`` never
    change during an EM run.  This cache computes them once, plus
    the shifted counts and the 17-node ρ grid used by the R channel,
    so each κ evaluation costs ~3 vector ``gammaln`` calls instead
    of ~6 for G and ~19 for R (was ~34 for R before).
    """

    __slots__ = (
        "valid", "k", "n", "n_minus_k", "log_binom",
        "p_rho", "one_minus_p_rho",
    )

    def __init__(self, stats: dict[str, np.ndarray]):
        self.valid, k, n = _compute_k_sense_valid(stats)
        self.k = k
        self.n = n
        self.n_minus_k = n - k
        # log C(n, k) — depends only on (k, n); cache once.
        if k.size:
            self.log_binom = (
                _vlgamma(n + 1.0) - _vlgamma(k + 1.0) - _vlgamma(n - k + 1.0)
            )
        else:
            self.log_binom = np.empty(0, dtype=np.float64)
        # ρ quadrature grid is fixed by the module constants and SS,
        # but SS is needed at call time; populate on first use via setter.
        self.p_rho = None
        self.one_minus_p_rho = None

    def set_ss(self, strand_specificity: float) -> None:
        ss = float(strand_specificity)
        p = _RHO_NODES * ss + (1.0 - _RHO_NODES) * 0.5  # (17,)
        self.p_rho = p
        self.one_minus_p_rho = 1.0 - p


def _bb_g_loglik_cached(cache: _StrandBBCache, kappa_G: float) -> np.ndarray:
    """Cached version of G-class BB log-likelihood (symmetric, α=β=κ/2).

    Avoids recomputing the (k, n)-only ``log_binom`` term and collapses
    the denominator to a single scalar ``2·lgamma(α) - lgamma(2α)``.
    """
    a = 0.5 * float(kappa_G)
    # Numerator: lgamma(k+a) + lgamma(n-k+a) - lgamma(n+2a)
    log_B_num = (
        _vlgamma(cache.k + a)
        + _vlgamma(cache.n_minus_k + a)
        - _vlgamma(cache.n + 2.0 * a)
    )
    # Denominator (scalar): 2*lgamma(a) - lgamma(2a)
    log_B_den = 2.0 * math.lgamma(a) - math.lgamma(2.0 * a)
    return cache.log_binom + log_B_num - log_B_den


def _bb_r_loglik_cached(cache: _StrandBBCache, kappa_R: float) -> np.ndarray:
    """Cached version of R-class ρ-integrated BB log-likelihood.

    Optimisations vs. the reference implementation:
      * ``log_binom`` (k/n-only) is cached once per EM run.
      * ``lgamma(n + κ_R)`` depends on n and κ only (not on ρ), so it
        is a vector of size N computed once per κ — not a (17, N)
        tensor as it would be if fused with the α_j, β_j broadcast.
      * ``lgamma(α_j) + lgamma(β_j) - lgamma(α_j + β_j)`` has shape (17,)
        (scalars over regions).
    """
    if cache.p_rho is None:
        raise RuntimeError("_StrandBBCache.set_ss() must be called first")
    kR = float(kappa_R)
    alpha_j = (kR * cache.p_rho)[:, None]              # (17, 1)
    beta_j = (kR * cache.one_minus_p_rho)[:, None]     # (17, 1)
    # Numerator per node: lgamma(k + α_j) + lgamma(n-k + β_j) - lgamma(n + κ)
    # The third term does NOT depend on ρ; compute as (N,) once and broadcast.
    lgamma_n_plus_k = _vlgamma(cache.n + kR)           # (N,)
    log_B_num = (
        _vlgamma(cache.k[None, :] + alpha_j)
        + _vlgamma(cache.n_minus_k[None, :] + beta_j)
        - lgamma_n_plus_k[None, :]
    )                                                  # (17, N)
    # Denominator (17,): lgamma(α_j) + lgamma(β_j) - lgamma(κ)  [since α+β=κ]
    log_B_den = (
        _vlgamma(kR * cache.p_rho)
        + _vlgamma(kR * cache.one_minus_p_rho)
        - math.lgamma(kR)
    )                                                  # (17,)
    log_p_nodes = cache.log_binom[None, :] + log_B_num - log_B_den[:, None]
    return _logsumexp(log_p_nodes + _RHO_LOGW[:, None], axis=0)


def _bb_g_loglik(k: np.ndarray, n: np.ndarray, kappa_G: float) -> np.ndarray:
    """Per-region log P(k | n, G, κ_G) under BB(κ_G/2, κ_G/2).  Shape (N,).

    Uncached convenience wrapper (kept for tests / external callers).
    """
    a = 0.5 * float(kappa_G)
    return _bb_logpmf(k, n, a, a)


def _bb_r_loglik(
    k: np.ndarray,
    n: np.ndarray,
    kappa_R: float,
    strand_specificity: float,
) -> np.ndarray:
    """Per-region log P(k | n, R, κ_R) marginalising ρ ~ Uniform(0, 1).

    Evaluates
        ∫₀¹ BB(k; n, κ_R · p(ρ), κ_R · (1-p(ρ))) dρ
    with p(ρ) = ρ · SS + (1-ρ) · 0.5, via 17-node Gauss-Legendre on [0, 1].
    Returns log-values of shape (N,).  Uncached convenience wrapper.
    """
    ss = float(strand_specificity)
    # p(ρ) spans [0.5, SS].  At SS → 0.5 the R and G integrands coincide
    # and the LLR collapses to 0 automatically — no gating required.
    p_rho = _RHO_NODES * ss + (1.0 - _RHO_NODES) * 0.5     # (17,)
    alpha = (float(kappa_R) * p_rho)[:, None]              # (17, 1)
    beta = (float(kappa_R) * (1.0 - p_rho))[:, None]       # (17, 1)
    log_p_nodes = _bb_logpmf(k[None, :], n[None, :], alpha, beta)  # (17, N)
    return _logsumexp(log_p_nodes + _RHO_LOGW[:, None], axis=0)


def _strand_llr_bb_rho(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    kappa_G: float,
    kappa_R: float,
    cache: "_StrandBBCache | None" = None,
) -> np.ndarray:
    """Per-region strand LLR = log P(k|n,G,κ_G) − log P(k|n,R,κ_R).

    Zero for invalid regions (n=0 or tx_strand=0).  No library-level
    gate: SS near 0.5 makes G and R integrands coincide and the LLR
    vanishes naturally.  If ``cache`` is supplied, uses the fast
    path; otherwise builds a one-shot cache internally.
    """
    n_regions = stats["n_unspliced"].shape[0]
    llr = np.zeros(n_regions, dtype=np.float64)
    if cache is None:
        cache = _StrandBBCache(stats)
        cache.set_ss(strand_specificity)
    if cache.k.size == 0:
        return llr
    ll_g = _bb_g_loglik_cached(cache, kappa_G)
    ll_r = _bb_r_loglik_cached(cache, kappa_R)
    llr[cache.valid] = ll_g - ll_r
    return llr


def _strand_posterior(llr: np.ndarray, pi_0: float = 0.5) -> np.ndarray:
    """Strand-only γ̃_i = sigmoid(logit(π₀) + LLR_strand_i).  Diagnostic."""
    pi_safe = float(np.clip(pi_0, _EPS, 1.0 - _EPS))
    logit = math.log(pi_safe / (1.0 - pi_safe))
    return 1.0 / (1.0 + np.exp(-np.clip(logit + llr, -500.0, 500.0)))


def _golden_section_max(
    f,
    a: float,
    b: float,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> float:
    """Argmax of *f* on [*a*, *b*] via golden-section search."""
    gr = (math.sqrt(5.0) + 1.0) / 2.0
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    fc, fd = f(c), f(d)
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - (b - a) / gr
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + (b - a) / gr
            fd = f(d)
    return (a + b) / 2.0


def _golden_section_max_logx(
    f_on_x,
    a: float,
    b: float,
    log_tol: float = 0.03,
    max_iter: int = 40,
) -> float:
    """Argmax of ``x ↦ f_on_x(x)`` for x ∈ [a, b] searched in log space.

    κ spans four orders of magnitude (0.1 → 1000), so a linear-space
    search with tol=1e-3 wastes ~28 evals.  Searching on ln(κ) with a
    relative tolerance (default ``log_tol`` = 0.03 ↔ ~3% in κ) needs
    only ~9 evals while preserving convergence quality: BB log-
    likelihoods are slowly varying in κ around the optimum, and the
    E-step only cares about the integrated posterior which is
    insensitive to κ at the 1% level.
    """
    la, lb = math.log(a), math.log(b)
    gr = (math.sqrt(5.0) + 1.0) / 2.0
    lc = lb - (lb - la) / gr
    ld = la + (lb - la) / gr
    fc, fd = f_on_x(math.exp(lc)), f_on_x(math.exp(ld))
    for _ in range(max_iter):
        if abs(lb - la) < log_tol:
            break
        if fc > fd:
            lb, ld, fd = ld, lc, fc
            lc = lb - (lb - la) / gr
            fc = f_on_x(math.exp(lc))
        else:
            la, lc, fc = lc, ld, fd
            ld = la + (lb - la) / gr
            fd = f_on_x(math.exp(ld))
    return math.exp(0.5 * (la + lb))


def _estimate_kappa_G_soft(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    kappa_lo: float = _KAPPA_LO,
    kappa_hi: float = _KAPPA_HI,
    min_eff_regions: float = 2.0,
    cache: "_StrandBBCache | None" = None,
    prev_kappa: float | None = None,
) -> float:
    """γ-weighted soft MLE for κ_G on the G-class symmetric Beta-Binomial.

    When the effective γ-weight is below ``min_eff_regions`` (not
    enough data attributed to G to constrain overdispersion), falls
    back to the Binomial limit ``kappa_hi`` — the strongest, least
    conservative default.  This prevents a run-away low-κ fit from
    absorbing extreme k values through the symmetric BB's U-shape.

    If ``prev_kappa`` is supplied, narrows the initial search bracket
    to [prev/4, prev*4] in log-κ before widening to [lo, hi] only
    when the optimum is at the bracket edge.  κ trajectories within
    an EM run change by < 2× per iteration once the mixture firms up,
    so this warm-start typically halves the search cost.
    """
    if cache is None:
        cache = _StrandBBCache(stats)
    if cache.k.size < 2:
        return kappa_hi
    w = np.asarray(gamma[cache.valid], dtype=np.float64)
    if float(w.sum()) < min_eff_regions:
        return kappa_hi

    def obj(kappa: float) -> float:
        return float(np.sum(w * _bb_g_loglik_cached(cache, kappa)))

    return _bracketed_log_search(obj, kappa_lo, kappa_hi, prev_kappa)


def _estimate_kappa_R_soft(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    strand_specificity: float,
    kappa_lo: float = _KAPPA_LO,
    kappa_hi: float = _KAPPA_HI,
    min_eff_regions: float = 2.0,
    cache: "_StrandBBCache | None" = None,
    prev_kappa: float | None = None,
) -> float:
    """(1-γ)-weighted soft MLE for κ_R on the ρ-integrated R-class BB.

    Same low-weight fallback and warm-start strategy as
    :func:`_estimate_kappa_G_soft`.
    """
    if cache is None:
        cache = _StrandBBCache(stats)
        cache.set_ss(strand_specificity)
    if cache.k.size < 2:
        return kappa_hi
    w = np.asarray(1.0 - gamma[cache.valid], dtype=np.float64)
    if float(w.sum()) < min_eff_regions:
        return kappa_hi

    def obj(kappa: float) -> float:
        return float(np.sum(w * _bb_r_loglik_cached(cache, kappa)))

    return _bracketed_log_search(obj, kappa_lo, kappa_hi, prev_kappa)


def _bracketed_log_search(
    obj,
    kappa_lo: float,
    kappa_hi: float,
    prev_kappa: float | None,
    warm_factor: float = 4.0,
) -> float:
    """Golden-section maximisation in log-κ with an optional warm-start bracket.

    When ``prev_kappa`` is supplied, searches [prev/warm_factor, prev*warm_factor]
    first; if the optimum lies at the edge (suggesting the warm window
    is too narrow), widens to [kappa_lo, kappa_hi] and re-searches.
    In the typical case (after iter 2) this requires ~7 evaluations;
    a full widen adds ~9 more only when warranted.
    """
    if prev_kappa is None or not (kappa_lo < prev_kappa < kappa_hi):
        return _golden_section_max_logx(obj, kappa_lo, kappa_hi)
    a = max(prev_kappa / warm_factor, kappa_lo)
    b = min(prev_kappa * warm_factor, kappa_hi)
    if b <= a:
        return _golden_section_max_logx(obj, kappa_lo, kappa_hi)
    best = _golden_section_max_logx(obj, a, b)
    # Edge check: if the optimum sits on the warm-window boundary
    # (within 20% of an endpoint in log space) fall through to the
    # wider search.
    log_span = math.log(b / a)
    if log_span <= 0.0:
        return best
    edge_margin = 0.20 * log_span
    log_best = math.log(best)
    if log_best - math.log(a) < edge_margin and a > kappa_lo * 1.01:
        return _golden_section_max_logx(obj, kappa_lo, kappa_hi)
    if math.log(b) - log_best < edge_margin and b < kappa_hi * 0.99:
        return _golden_section_max_logx(obj, kappa_lo, kappa_hi)
    return best


def _init_kappa_pooled(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
) -> tuple[float, float]:
    """Initialise (κ_G, κ_R) at the Binomial limit ``_KAPPA_HI``.

    The pooled MLE fit at γ = 0.5 is unstable in homogeneous-population
    cases: a pure-RNA library would force κ_G to the lower bound (a
    U-shaped symmetric BB happily absorbs any extreme k), producing a
    degenerate iter-0 E-step where strand LLR collapses to zero.  The
    correct default before any classification is "no region-level
    overdispersion" — i.e., Binomial(0.5) for G and Binomial(p(ρ)) for R.
    The γ-weighted M-step then adapts both κ's as the classification
    firms up; the low-weight fallback in ``_estimate_kappa_*_soft``
    keeps the Binomial default whenever a component has too little
    effective data to constrain overdispersion.
    """
    return _KAPPA_HI, _KAPPA_HI


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
    lam_G: float,
    mu_R: float,
    sigma_R: float,
    pi_soft: float,
    strand_specificity: float,
    kappa_G: float,
    kappa_R: float,
    llr_fl: np.ndarray | None,
    cache: "_StrandBBCache | None" = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute γ_i = P(z=G | data_i) and the per-region strand LLR.

    Hard:
      * k^s > 0 → γ = 0
      * ineligible → γ = π_soft

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

    llr_count = _count_llr_poisson_ln(
        stats["n_unspliced"][soft],
        stats["mappable_bp"][soft],
        lam_G,
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
) -> tuple[float, float, float, float, float]:
    """Update (λ_G, μ_R, σ_R, π, π_soft) from posteriors.

    * λ_G: γ-weighted Poisson MLE = Σ γ_i k^u_i / Σ γ_i E_i over eligible.
    * (μ_R, σ_R²): moments of ``log(max(k^u_i/E_i − λ_G, eps_μ))`` taken
      **only over hard-spliced anchor regions** (``k^s > 0 & eligible``).
      Those regions are definitively RNA — their per-bp rate directly
      reflects the expressed-gene distribution.  Using ``(1−γ)``-weights
      over *all* eligible regions (the naive Gaussian-mixture update)
      is circular: it lets class R absorb empty/low-count regions via
      a broadening σ_R, which then reclassifies more regions as R,
      collapsing λ_G.  Hard anchors break that loop.
    * π: γ-mean over eligible.
    * π_soft: γ-mean over eligible & k^s == 0.
    """
    n_u = stats["n_unspliced"]
    E = stats["mappable_bp"]
    n_s = stats["n_spliced"]

    if not eligible.any():
        return 0.0, 0.0, 1.0, 0.5, 0.5

    g = gamma[eligible]
    ku = n_u[eligible]
    ee = E[eligible]
    ns = n_s[eligible]

    # λ_G update
    num = float(np.sum(g * ku))
    den = float(np.sum(g * ee))
    lam_G = num / den if den > _EPS else 0.0

    # μ_R, σ_R: anchor-only update.  If no anchors are available
    # (e.g. tiny test cases), fall back to a broad uninformative
    # LogNormal — the E-step count LLR will then depend mostly on
    # Poisson(λ_G).
    anchor_mask = (n_s > 0) & eligible
    if anchor_mask.any():
        ee_a = E[anchor_mask]
        ku_a = n_u[anchor_mask]
        rate_a = np.where(ee_a > 0, ku_a / ee_a, 0.0)
        med_E = float(np.median(ee_a[ee_a > 0])) if np.any(ee_a > 0) else 1.0
        eps_mu = 1.0 / max(med_E, 1.0)
        mu_hat = np.maximum(rate_a - lam_G, eps_mu)
        log_mu = np.log(mu_hat)
        if log_mu.size >= 2:
            mu_R = float(np.mean(log_mu))
            var_R = float(np.var(log_mu))
            sigma_R = math.sqrt(max(var_R, 1e-4))
        else:
            # Single anchor: center on it but keep a broad σ.
            mu_R = float(log_mu.item())
            sigma_R = 1.5
    else:
        # No anchors — uninformative broad prior centered well above λ_G.
        mu_R = math.log(max(10.0 * lam_G, 1e-6))
        sigma_R = 2.0

    # Mixing proportions
    pi = float(np.mean(g))
    soft_mask = ns == 0
    if soft_mask.any():
        pi_soft = float(np.mean(g[soft_mask]))
    else:
        pi_soft = pi
    pi = float(np.clip(pi, _EPS, 1.0 - _EPS))
    pi_soft = float(np.clip(pi_soft, _EPS, 1.0 - _EPS))

    return lam_G, mu_R, sigma_R, pi, pi_soft


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
    diagnostics: bool = False,
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
    """
    n_u = stats["n_unspliced"]
    n_s = stats["n_spliced"]
    E = stats["mappable_bp"]
    n_regions = n_u.shape[0]

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

    # Initial λ_G: robust median-split pooled rate over soft eligible regions.
    soft_mask = eligible & (n_s == 0)
    if soft_mask.any():
        lam_G = _robust_initial_lam_G(n_u[soft_mask], E[soft_mask])
    else:
        lam_G = 0.0

    # Initial μ_R, σ_R from spliced-seed regions (definitively RNA).
    expressed_seed_mask = (n_s > 0) & eligible
    if expressed_seed_mask.any():
        ee = E[expressed_seed_mask]
        rates = np.where(ee > 0, n_u[expressed_seed_mask] / ee, 0.0)
        med_E = float(np.median(ee[ee > 0])) if np.any(ee > 0) else 1.0
        eps_mu = 1.0 / max(med_E, 1.0)
        mu_hat = np.maximum(rates - lam_G, eps_mu)
        log_mu = np.log(mu_hat)
        mu_R = float(np.mean(log_mu))
        sigma_R = max(float(np.std(log_mu)), 0.5)
    else:
        mu_R = math.log(max(10.0 * lam_G, 1e-6))
        sigma_R = 2.0

    # --- Initial (κ_G, κ_R): pooled MLE at γ = 0.5.  No magic numbers.
    kappa_G, kappa_R = _init_kappa_pooled(stats, strand_specificity)

    n_soft = int(soft_mask.sum())
    n_elig = int(eligible.sum())
    pi = pi_init
    pi_soft = pi_init

    logger.info(
        "EM calibration init: n_elig=%d, n_soft=%d, n_spliced_hard=%d, "
        "π=%.3f, π_soft=%.3f, λ_G=%.3e, μ_R=%.2f, σ_R=%.2f, "
        "κ_G=%.2f, κ_R=%.2f  (SS=%.3f)",
        n_elig, n_soft, n_spliced_hard,
        pi, pi_soft, lam_G, mu_R, sigma_R,
        kappa_G, kappa_R, float(strand_specificity),
    )

    history: list[dict] = []
    converged = False
    n_iter = 0
    gdna_fl_model = None

    prev_pi_soft = pi_soft
    llr_strand_full = np.zeros(n_regions, dtype=np.float64)

    # One-shot cache for BB strand log-likelihoods (k, n, log_binom are
    # invariant across EM iterations and across both κ_G / κ_R searches).
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
            stats, eligible, lam_G, mu_R, sigma_R,
            pi_soft, strand_specificity, kappa_G, kappa_R, llr_fl,
            cache=bb_cache,
        )

        # --- M-step ---
        lam_G, mu_R, sigma_R, pi, pi_soft = _m_step(stats, eligible, gamma)

        # --- κ updates (γ-weighted soft MLEs, self-consistent within EM) ---
        # Warm-start the log-κ search with the previous iteration's
        # value once κ has moved off the Binomial-limit default.
        prev_kG = kappa_G if kappa_G < _KAPPA_HI else None
        prev_kR = kappa_R if kappa_R < _KAPPA_HI else None
        kappa_G = _estimate_kappa_G_soft(
            stats, gamma, cache=bb_cache, prev_kappa=prev_kG,
        )
        kappa_R = _estimate_kappa_R_soft(
            stats, gamma, strand_specificity, cache=bb_cache,
            prev_kappa=prev_kR,
        )

        # --- FL model update (γ-weighted gDNA FL) ---
        if fl_region_ids.size > 0:
            gdna_fl_model = _build_fl_histogram(
                fl_region_ids, fl_frag_lens, gamma, n_regions,
            )

        delta = abs(pi_soft - prev_pi_soft)
        if diagnostics:
            history.append({
                "iter": n_iter,
                "lam_G": lam_G,
                "mu_R": mu_R,
                "sigma_R": sigma_R,
                "pi": pi,
                "pi_soft": pi_soft,
                "kappa_G": kappa_G,
                "kappa_R": kappa_R,
                "delta_pi_soft": delta,
                "mean_gamma_soft": float(np.mean(gamma[eligible & (n_s == 0)]))
                    if (eligible & (n_s == 0)).any() else float("nan"),
            })
        logger.debug(
            "EM iter %d: λ_G=%.3e μ_R=%.2f σ_R=%.2f π=%.3f π_soft=%.3f "
            "κ_G=%.2f κ_R=%.2f Δ=%.2e",
            n_iter, lam_G, mu_R, sigma_R, pi, pi_soft, kappa_G, kappa_R, delta,
        )

        if delta < convergence_tol:
            converged = True
            break
        prev_pi_soft = pi_soft

    logger.info(
        "EM calibration %s in %d iters: λ_G=%.3e, μ_R=%.2f, σ_R=%.2f, π=%.3f, "
        "κ_G=%.2f, κ_R=%.2f",
        "converged" if converged else "did NOT converge",
        n_iter, lam_G, mu_R, sigma_R, pi, kappa_G, kappa_R,
    )

    # Strand-only diagnostic posterior (what does strand alone say?).
    region_gamma_strand = _strand_posterior(llr_strand_full, pi_0=0.5)

    return EMFit(
        lam_G=float(lam_G),
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
    )
