"""Mappability-aware EM mixture deconvolution for gDNA calibration (v4).

Two-component mixture over exonic regions:

* Class G (gDNA-only):  k^u_i ~ Poisson(λ_G · E_i),  k^s_i = 0 (hard).
* Class R (expressed):   k^u_i ~ Poisson((λ_G + (1-ρ)·μ_i) · E_i),
                         k^s_i ~ Poisson(ρ·μ_i·E_i),
                         log μ_i ~ N(μ_R, σ_R²).

Signals fused as additive LLRs on logit γ:
  1. Count LLR via Poisson-LogNormal marginal (Gauss-Hermite quadrature).
  2. Strand LLR via Binomial(n, 0.5) vs Binomial(n, p_i), gated by a
     data-driven aggregate z-test on observed sense-bias.
  3. Fragment-length LLR via shape-normalised histogram ratios (from iter 1+).

See ``docs/calibration/calibration_v4_em_mixture_plan.md``.
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
    # Whether the strand channel passed the aggregate z-test and was used
    strand_used: bool
    # Aggregate sense-bias z-score log: z = (p̂_sense - 0.5)/sqrt(0.25/N)
    strand_z: float
    # Strand LLR mode actually used ("binomial" or "betabinom")
    strand_llr_mode: str = "binomial"
    # Beta-Binomial shared concentration κ (0.0 when binomial mode)
    kappa: float = 0.0
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
# Strand LLR (Binomial: 0.5 vs p_i)
# ---------------------------------------------------------------------------


def _aggregate_strand_z(
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
) -> tuple[float, float, float]:
    """Aggregate one-sided z-score for "library is stranded".

    Pools sense/antisense counts across all eligible regions with a
    known tx_strand and tests H0: p_sense = 0.5 vs H1: p_sense > 0.5
    under a Binomial null.  Returns ``(z, p_hat, N)``; ``z = 0`` when
    ``N < 1``.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    tx_strand = stats["tx_strand"]
    mask = eligible & (tx_strand != 0) & (n_unspliced > 0)
    if not mask.any():
        return 0.0, 0.5, 0.0
    k_sense_i = np.zeros_like(n_unspliced)
    pos_mask = mask & (tx_strand == 1)
    neg_mask = mask & (tx_strand == -1)
    k_sense_i[pos_mask] = n_unspliced[pos_mask] - n_pos[pos_mask]
    k_sense_i[neg_mask] = n_pos[neg_mask]
    K = float(k_sense_i[mask].sum())
    N = float(n_unspliced[mask].sum())
    if N < 1.0:
        return 0.0, 0.5, 0.0
    p_hat = K / N
    se = math.sqrt(0.25 / N)
    z = (p_hat - 0.5) / se if se > 0 else 0.0
    return z, p_hat, N


def _strand_llr(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    enabled: bool,
    noise_floor: float = 1e-12,
) -> np.ndarray:
    """Binomial strand LLR: log Binom(k_sense | n, 0.5) − log Binom(k_sense | n, p_i).

    Returns zeros unless ``enabled`` is True (the caller decides via an
    aggregate z-test on observed strand bias).  Zero contribution also
    when the region's gene orientation is ambiguous or n_unspliced < 2.

    Parameters
    ----------
    noise_floor
        Effective floor applied to ``1 − ss`` under the RNA model.
        This is ``max(ε_bio, ε_CI)`` computed by the caller and
        captures structural antisense (ε_bio) and strand-trainer
        measurement uncertainty (ε_CI).  Defaults to 1e-12 for
        backwards compatibility in callers that don't yet pass it;
        production callers should pass the config-derived value.
        See ``docs/calibration/strand_llr_noise_floor_design.md``.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    tx_strand = stats["tx_strand"]
    n = n_unspliced.shape[0]
    llr = np.zeros(n, dtype=np.float64)

    if not enabled:
        return llr

    valid = (n_unspliced >= 2) & (tx_strand != 0)
    if not valid.any():
        return llr

    eps = float(max(noise_floor, _EPS))
    ss_c = float(np.clip(strand_specificity, 0.5, 1.0 - eps))
    # Sense-read count per region (R1-antisense convention).
    k_sense = np.zeros(n, dtype=np.float64)
    k_sense[tx_strand == 1] = n_unspliced[tx_strand == 1] - n_pos[tx_strand == 1]
    k_sense[tx_strand == -1] = n_pos[tx_strand == -1]

    k = k_sense[valid]
    tot = n_unspliced[valid]
    log_ratio_sense = math.log(0.5 / ss_c)
    log_ratio_anti = math.log(0.5 / (1.0 - ss_c))
    llr[valid] = k * log_ratio_sense + (tot - k) * log_ratio_anti
    return llr


# ---------------------------------------------------------------------------
# Beta-Binomial strand LLR (shared κ)
# ---------------------------------------------------------------------------


def _betaln_scalar(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _betaln_vec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return _vlgamma(a) + _vlgamma(b) - _vlgamma(a + b)


def _compute_k_sense_valid(
    stats: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (valid_mask, k_sense[valid], n_unspliced[valid]).

    Same gating as the Binomial strand LLR: tx_strand != 0 and
    n_unspliced >= 2.  ``k_sense`` is the R1-antisense-convention
    sense-read count.
    """
    n_unspliced = stats["n_unspliced"]
    n_pos = stats["n_pos"]
    tx_strand = stats["tx_strand"]
    valid = (n_unspliced >= 2) & (tx_strand != 0)
    if not valid.any():
        return valid, np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
    k_sense = np.zeros(n_unspliced.shape[0], dtype=np.float64)
    k_sense[tx_strand == 1] = n_unspliced[tx_strand == 1] - n_pos[tx_strand == 1]
    k_sense[tx_strand == -1] = n_pos[tx_strand == -1]
    return valid, k_sense[valid].astype(np.float64), n_unspliced[valid].astype(np.float64)


def _strand_llr_betabinom(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    kappa: float,
    enabled: bool,
    noise_floor: float = 1e-12,
) -> np.ndarray:
    """Beta-Binomial strand LLR with shared dispersion κ.

    gDNA: BetaBinom(k_sense | n, κ/2, κ/2)  — symmetric around 0.5.
    RNA:  BetaBinom(k_sense | n, κ·ss, κ·(1−ss)).

    LLR = log P_G − log P_R (the C(n,k) combinatorial cancels exactly).
    Degenerates to the Binomial LLR as κ → ∞ and vanishes identically
    when ss = 0.5 (G and R coincide).

    Applies the same gating (``tx_strand != 0``, ``n_unspliced ≥ 2``)
    as the Binomial channel, and the same ``(0.5, 1-eps)`` clip on
    ``ss`` for numerical safety.
    """
    n = stats["n_unspliced"].shape[0]
    llr = np.zeros(n, dtype=np.float64)
    if not enabled or kappa <= 0.0:
        return llr
    valid, k, tot = _compute_k_sense_valid(stats)
    if k.size == 0:
        return llr

    eps = float(max(noise_floor, _EPS))
    ss_c = float(np.clip(strand_specificity, 0.5, 1.0 - eps))

    a_g = kappa / 2.0
    a_e = kappa * ss_c
    b_e = kappa * (1.0 - ss_c)

    ll_g = _betaln_vec(k + a_g, tot - k + a_g) - _betaln_scalar(a_g, a_g)
    ll_r = _betaln_vec(k + a_e, tot - k + b_e) - _betaln_scalar(a_e, b_e)

    llr[valid] = ll_g - ll_r
    return llr


def _golden_section_max(
    f,
    a: float,
    b: float,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> float:
    """Argmax of *f* on [*a*, *b*] via golden-section search (no scipy)."""
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


def _estimate_kappa_marginal(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    strand_specificity: float,
    *,
    noise_floor: float = 1e-12,
    kappa_lo: float = 0.01,
    kappa_hi: float = 500.0,
) -> float | None:
    """Estimate shared κ via exact mixture marginal log-likelihood.

    Maximises

        ℓ(κ) = Σ_r log[ γ_r · P_G(k_r | n_r, κ) + (1 − γ_r) · P_R(k_r | n_r, κ) ]

    where P_G = BetaBin(·, κ/2, κ/2) and P_R = BetaBin(·, κ·ss, κ·(1-ss)).

    When one component carries less than one effective region of weight
    (``sum(γ) < 1`` or ``sum(1−γ) < 1``), the responsibilities are too
    extreme for a stable two-component fit, so γ is reset to 0.5 — the
    objective collapses to a single symmetric Beta-Binomial, which is
    label-free and always stable.

    Returns ``None`` when fewer than 3 eligible (n ≥ 2, tx_strand != 0)
    regions are available.
    """
    valid, k, n = _compute_k_sense_valid(stats)
    if k.size < 3:
        return None

    g = gamma[valid].astype(np.float64)
    if float(g.sum()) < 1.0 or float((1.0 - g).sum()) < 1.0:
        g = np.full_like(g, 0.5)

    eps = float(max(noise_floor, _EPS))
    ss = float(np.clip(strand_specificity, 0.5, 1.0 - eps))

    with np.errstate(divide="ignore"):
        log_g = np.log(np.clip(g, _EPS, 1.0))
        log_1mg = np.log(np.clip(1.0 - g, _EPS, 1.0))

    def _ll(kappa: float) -> float:
        a_g = kappa / 2.0
        ll_g = _betaln_vec(k + a_g, n - k + a_g) - _betaln_scalar(a_g, a_g)
        a_e = kappa * ss
        b_e = kappa * (1.0 - ss)
        ll_r = _betaln_vec(k + a_e, n - k + b_e) - _betaln_scalar(a_e, b_e)
        return float(np.sum(np.logaddexp(log_g + ll_g, log_1mg + ll_r)))

    kappa = _golden_section_max(_ll, kappa_lo, kappa_hi)
    return max(kappa, 0.0)


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
    strand_enabled: bool,
    llr_fl: np.ndarray | None,
    strand_noise_floor: float = 1e-12,
    strand_llr_mode: str = "binomial",
    kappa: float = 0.0,
) -> np.ndarray:
    """Compute γ_i = P(z=G | data_i) for eligible soft regions.

    Hard:
      * k^s > 0 → γ = 0
      * ineligible → γ = π_soft
    """
    n = stats["n_unspliced"].shape[0]
    n_s = stats["n_spliced"]
    gamma = np.full(n, pi_soft, dtype=np.float64)

    hard_expr = (n_s > 0) & eligible
    gamma[hard_expr] = 0.0

    soft = eligible & ~hard_expr
    if not soft.any():
        return gamma

    pi_safe = float(np.clip(pi_soft, _EPS, 1.0 - _EPS))
    log_prior_odds = math.log(pi_safe / (1.0 - pi_safe))

    llr_count = _count_llr_poisson_ln(
        stats["n_unspliced"][soft],
        stats["mappable_bp"][soft],
        lam_G,
        mu_R,
        sigma_R,
    )
    if strand_llr_mode == "betabinom" and kappa > 0.0:
        llr_strand_full = _strand_llr_betabinom(
            stats, strand_specificity, kappa, strand_enabled,
            noise_floor=strand_noise_floor,
        )
    else:
        llr_strand_full = _strand_llr(
            stats, strand_specificity, strand_enabled, noise_floor=strand_noise_floor
        )
    llr_strand = llr_strand_full[soft]

    if llr_fl is not None:
        llr_fl_soft = llr_fl[soft]
    else:
        llr_fl_soft = 0.0

    log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft
    log_odds = np.clip(log_odds, -500.0, 500.0)
    gamma[soft] = 1.0 / (1.0 + np.exp(-log_odds))
    return gamma


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
    strand_noise_floor: float = 1e-12,
    strand_llr_mode: str = "binomial",
    max_iterations: int = 50,
    convergence_tol: float = 1e-4,
    strand_z_threshold: float = 3.0,
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
        Library SS ∈ [0.5, 1.0].  Used to compute the per-region
        strand LLR *when the aggregate sense-bias z-test passes*.
    mean_frag_len
        Library mean fragment length.  Used only as the eligibility
        floor for ``mappable_bp_i`` (a region must carry at least one
        fragment-worth of mappable length).
    strand_z_threshold
        Minimum aggregate z for ``p̂_sense = Σk_sense / Σn_unspliced``
        to enable the strand LLR channel.  At ``z = 3`` this
        corresponds to ~99.9% one-sided confidence that the library
        is stranded given the observed calibration-region counts.
        Replaces the old hard SS threshold.
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
            strand_used=False,
            strand_z=0.0,
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
        # Validate bounds ONCE up-front so the per-iteration hot-path
        # (_build_fl_histogram, _fl_llr — both called on up to ~15M items
        # per EM iter) can skip the mask entirely. This saves ~4.6 s/iter
        # on VCaP-scale inputs.
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

    # --- Data-driven strand channel gating ---
    # One-sided binomial z-test on the *aggregate* sense fraction.
    # This replaces a hard SS threshold: the question is not "is the
    # library stated to be stranded?" but "do the observed counts
    # carry enough evidence to support a sense-biased model over a
    # symmetric one?"  At z = 3 we require ~3σ confidence given the
    # pooled count.  With, e.g., N=1,000 fragments and true SS=0.55
    # this just clears; with N=100 it does not.  No per-library
    # tuning required.
    strand_z, p_sense_hat, n_strand_pool = _aggregate_strand_z(stats, eligible)
    strand_used = strand_z >= float(strand_z_threshold)
    logger.info(
        "EM calibration strand test: N=%.0f  p_hat=%.3f  z=%.2f  threshold=%.2f  → %s",
        n_strand_pool, p_sense_hat, strand_z, strand_z_threshold,
        "ENABLED" if strand_used else "disabled",
    )
    if strand_used:
        ss_clipped = float(np.clip(strand_specificity, 0.5, 1.0 - max(strand_noise_floor, _EPS)))
        logger.info(
            "EM calibration strand LLR: ss=%.6f → ss_clipped=%.6f  floor=%.4g  "
            "log_ratio_sense=%+.3f  log_ratio_anti=%+.3f nats (max_per_anti_read)",
            float(strand_specificity), ss_clipped, strand_noise_floor,
            math.log(0.5 / ss_clipped), math.log(0.5 / (1.0 - ss_clipped)),
        )

    # --- Anchor-only seed ---
    gamma, pi_init, n_spliced_hard = _seed_gamma(stats, eligible)

    # Initial λ_G: robust median-split pooled rate over soft eligible
    # regions.  Completely parameter-free and well-defined for any
    # region count ≥ 1.
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
        # No spliced regions → broad LogNormal prior; EM will adapt.
        mu_R = math.log(max(10.0 * lam_G, 1e-6))
        sigma_R = 2.0

    n_soft = int(soft_mask.sum())
    n_elig = int(eligible.sum())
    pi = pi_init
    pi_soft = pi_init  # uninformative at start; M-step updates.

    logger.info(
        "EM calibration init: n_elig=%d, n_soft=%d, n_spliced_hard=%d, "
        "π=%.3f, π_soft=%.3f, λ_G=%.3e, μ_R=%.2f, σ_R=%.2f",
        n_elig, n_soft, n_spliced_hard,
        pi, pi_soft, lam_G, mu_R, sigma_R,
    )

    history: list[dict] = []
    converged = False
    n_iter = 0
    gdna_fl_model = None  # built after iter 0

    # Shared Beta-Binomial κ (only updated when strand_llr_mode == "betabinom"
    # and the strand channel is enabled).
    kappa: float = 0.0
    if strand_llr_mode == "betabinom" and strand_used:
        k0 = _estimate_kappa_marginal(
            stats, gamma, strand_specificity,
            noise_floor=strand_noise_floor,
        )
        kappa = float(k0) if k0 is not None else 0.0
        logger.info(
            "EM calibration strand channel: mode=betabinom  κ_init=%.2f",
            kappa,
        )
    elif strand_used:
        logger.info("EM calibration strand channel: mode=binomial")

    prev_pi_soft = pi_soft

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
        gamma = _e_step(
            stats, eligible, lam_G, mu_R, sigma_R,
            pi_soft, strand_specificity, strand_used, llr_fl,
            strand_noise_floor=strand_noise_floor,
            strand_llr_mode=strand_llr_mode,
            kappa=kappa,
        )

        # --- M-step ---
        lam_G, mu_R, sigma_R, pi, pi_soft = _m_step(stats, eligible, gamma)

        # --- κ update (self-consistent within EM) ---
        if strand_llr_mode == "betabinom" and strand_used:
            k_new = _estimate_kappa_marginal(
                stats, gamma, strand_specificity,
                noise_floor=strand_noise_floor,
            )
            if k_new is not None:
                kappa = float(k_new)

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
                "kappa": kappa,
                "delta_pi_soft": delta,
                "mean_gamma_soft": float(np.mean(gamma[eligible & (n_s == 0)]))
                    if (eligible & (n_s == 0)).any() else float("nan"),
            })
        logger.debug(
            "EM iter %d: λ_G=%.3e μ_R=%.2f σ_R=%.2f π=%.3f π_soft=%.3f κ=%.2f Δ=%.2e",
            n_iter, lam_G, mu_R, sigma_R, pi, pi_soft, kappa, delta,
        )

        if delta < convergence_tol:
            converged = True
            break
        prev_pi_soft = pi_soft

    logger.info(
        "EM calibration %s in %d iters: λ_G=%.3e, μ_R=%.2f, σ_R=%.2f, π=%.3f, "
        "mode=%s, κ=%.2f",
        "converged" if converged else "did NOT converge",
        n_iter, lam_G, mu_R, sigma_R, pi, strand_llr_mode, kappa,
    )

    return EMFit(
        lam_G=float(lam_G),
        mu_R=float(mu_R),
        sigma_R=float(sigma_R),
        pi=float(pi),
        pi_soft=float(pi_soft),
        gamma=gamma,
        strand_used=bool(strand_used),
        strand_z=float(strand_z),
        strand_llr_mode=str(strand_llr_mode),
        kappa=float(kappa),
        n_iter=n_iter,
        converged=converged,
        n_eligible=n_elig,
        n_soft=n_soft,
        n_spliced_hard=n_spliced_hard,
        history=history,
    )
