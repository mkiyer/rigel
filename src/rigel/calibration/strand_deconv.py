"""Per-region strand-based gDNA/RNA deconvolution.

Phase 3 of the v5 strand-deconvolution plan
(``docs/fineregions/strand_model_impl_plan_v5.md``).

This module is a self-contained math layer:

* ``build_strand_region_counts`` projects the 12-channel fractional counts
  onto per-region ``(k_sense, k_antisense, n_total)`` triples using the
  region's transcript-strand class.
* ``deconvolve_regions_by_strand`` decomposes each region's total count
  into a (mean, upper, RNA-lower) gDNA estimate using the strand
  contrast ``p_r1_sense`` and the gDNA strand-balance overdispersion
  ``kappa_d``.
* ``screen_no_rna_exons`` identifies exon regions whose strand pattern
  is consistent with pure gDNA at confidence ``rna_lower_confidence``.
* ``estimate_kappa_d`` produces a refit of ``kappa_d`` using high-purity
  seed regions (intergenic + intron-only) plus self-trained no-RNA
  exons.

The module never mutates its inputs, never calls any orchestrator, and
never depends on EM. Wiring into the calibration pipeline is Phase 4.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import betaln, gammaln
from scipy.stats import norm

from ._arrays import PayloadArrays, RegionArrays
from .strand_balance import (
    StrandBalanceEstimate,
    estimate_strand_balance,
)
from .fractional_evidence import (
    TS_NEG,
    TS_POS,
    is_exon_any,
    is_intergenic,
    is_intron_only,
)
from .signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
)
from .strand_summary import STRAND_CONTRAST_NUMERICAL_FLOOR, StrandSummary


__all__ = [
    "DNA_STRAND_UNIFORM_PRIOR",
    "MAX_EXACT_POSTERIOR_N",
    "FLAG_INELIGIBLE",
    "FLAG_NEAR_UNSTRANDED",
    "FLAG_KAPPA_FALLBACK",
    "FLAG_APPROX_NORMAL",
    "FLAG_EXON_SELF_TRAIN",
    "StrandRegionCounts",
    "CompartmentStrandCounts",
    "RegionGdnaEstimate",
    "RegionGdnaChannelEstimate",
    "KappaDEstimate",
    "build_strand_region_counts",
    "build_compartment_strand_counts",
    "deconvolve_regions_by_strand",
    "deconvolve_compartments_by_strand",
    "strand_log_likelihood_d_grid",
    "screen_no_rna_exons",
    "estimate_kappa_d",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Fallback ``kappa_d`` when no seed regions are available. With a Beta(1, 1)
# prior on the gDNA strand-success probability the equivalent beta-binomial
# concentration is alpha + beta = 2, i.e. completely uninformative.
DNA_STRAND_UNIFORM_PRIOR: float = 2.0

# Maximum ``n_total`` (rounded up) for which we run the exact discrete
# posterior; larger regions fall back to the normal approximation. The
# exact path is O(n_total**2) so this caps each region at ~40k ops.
MAX_EXACT_POSTERIOR_N: int = 200

# uint8 flag bits.
FLAG_INELIGIBLE: int = 1 << 0
FLAG_NEAR_UNSTRANDED: int = 1 << 1
FLAG_KAPPA_FALLBACK: int = 1 << 2
FLAG_APPROX_NORMAL: int = 1 << 3
FLAG_EXON_SELF_TRAIN: int = 1 << 4


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StrandRegionCounts:
    """Per-region strand-folded counts ready for deconvolution.

    ``k_sense`` is the count on the channel matching the region's
    transcript strand (POS channel for ``TS_POS``, NEG for ``TS_NEG``);
    ``k_antisense`` is the opposite channel. Regions with ts_class in
    ``{TS_NONE, TS_AMBIG}`` are marked ineligible and carry zeros.
    """

    k_sense: np.ndarray  # float32[R]
    k_antisense: np.ndarray  # float32[R]
    n_total: np.ndarray  # float32[R]
    eligible: np.ndarray  # bool[R]
    p_r1_sense: float


@dataclass(frozen=True, slots=True)
class CompartmentStrandCounts:
    """Per-region strand-folded unspliced counts by region compartment."""

    contained_sense: np.ndarray  # float32[R]
    contained_antisense: np.ndarray  # float32[R]
    contained_total: np.ndarray  # float32[R]

    boundary_left_sense: np.ndarray  # float32[R]
    boundary_left_antisense: np.ndarray  # float32[R]
    boundary_left_total: np.ndarray  # float32[R]

    boundary_right_sense: np.ndarray  # float32[R]
    boundary_right_antisense: np.ndarray  # float32[R]
    boundary_right_total: np.ndarray  # float32[R]

    eligible: np.ndarray  # bool[R]
    p_r1_sense: float


@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    """Per-region gDNA fragment-count posterior summary."""

    n_total: np.ndarray  # float32[R]
    mean_count: np.ndarray  # float32[R] — E[D | data]
    upper_count: np.ndarray  # float32[R] — paired upper bound on D
    rna_lower_count: np.ndarray  # float32[R] — lower bound on R at confidence c
    precision: np.ndarray  # float32[R] in [0, 1]
    flags: np.ndarray  # uint8[R]
    kappa_d: float
    kappa_d_n_seed_regions: int
    kappa_d_n_exon_self_training: int
    p_r1_sense: float
    rna_lower_confidence: float


@dataclass(frozen=True, slots=True)
class RegionGdnaChannelEstimate:
    """Compartment-aware per-region gDNA/RNA deconvolution summary."""

    contained_mean: np.ndarray  # float32[R]
    contained_upper: np.ndarray
    contained_rna_lower: np.ndarray
    contained_precision: np.ndarray

    boundary_left_mean: np.ndarray  # float32[R]
    boundary_left_upper: np.ndarray
    boundary_left_rna_lower: np.ndarray
    boundary_left_precision: np.ndarray

    boundary_right_mean: np.ndarray  # float32[R]
    boundary_right_upper: np.ndarray
    boundary_right_rna_lower: np.ndarray
    boundary_right_precision: np.ndarray

    flags: np.ndarray  # uint16[R], union of compartment flags
    kappa_d: float
    p_r1_sense: float
    rna_lower_confidence: float


@dataclass(frozen=True, slots=True)
class KappaDEstimate:
    """Output of :func:`estimate_kappa_d`.

    Wraps the refit :class:`StrandBalanceEstimate` together with the seed
    and exon self-training region counts so callers can record both
    populations in diagnostics.
    """

    balance: StrandBalanceEstimate
    n_seed_regions: int
    n_exon_self_training: int

    @property
    def kappa(self) -> float:
        return float(self.balance.kappa)

    @property
    def fallback_used(self) -> bool:
        return _is_hard_fallback(self.balance)


# ---------------------------------------------------------------------------
# Build per-region strand counts
# ---------------------------------------------------------------------------


def _aggregate_pos_neg_all_compartments(
    payload_arrays: PayloadArrays,
) -> tuple[np.ndarray, np.ndarray]:
    """Sum POS / NEG channels across the three compartments × both splice states."""
    rc = payload_arrays.region_counts_sorted
    pos = np.zeros(rc.shape[0], dtype=np.float64)
    neg = np.zeros(rc.shape[0], dtype=np.float64)
    for compartment in (
        COMPARTMENT_CONTAINED,
        COMPARTMENT_BOUNDARY_LEFT,
        COMPARTMENT_BOUNDARY_RIGHT,
    ):
        for splice in (0, SPLICE_SPLICED):
            pos += rc[:, channel_index(compartment, splice, CHANNEL_STRAND_POS)].astype(
                np.float64, copy=False
            )
            neg += rc[:, channel_index(compartment, splice, CHANNEL_STRAND_NEG)].astype(
                np.float64, copy=False
            )
    return pos, neg


def _fold_pos_neg_by_transcript_strand(
    pos: np.ndarray,
    neg: np.ndarray,
    transcript_strand_class: np.ndarray,
    *,
    preserve_ineligible_total: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fold absolute POS/NEG channels into transcript-relative counts."""
    pos_arr = np.asarray(pos, dtype=np.float64)
    neg_arr = np.asarray(neg, dtype=np.float64)
    ts = np.asarray(transcript_strand_class)
    if pos_arr.shape != neg_arr.shape or pos_arr.shape != ts.shape:
        raise ValueError(
            "_fold_pos_neg_by_transcript_strand: pos, neg, and ts_class shapes must match; "
            f"got {pos_arr.shape}, {neg_arr.shape}, and {ts.shape}."
        )

    k_sense = np.zeros(pos_arr.shape, dtype=np.float32)
    k_antisense = np.zeros(pos_arr.shape, dtype=np.float32)
    is_pos = ts == TS_POS
    is_neg = ts == TS_NEG
    k_sense[is_pos] = pos_arr[is_pos].astype(np.float32, copy=False)
    k_antisense[is_pos] = neg_arr[is_pos].astype(np.float32, copy=False)
    k_sense[is_neg] = neg_arr[is_neg].astype(np.float32, copy=False)
    k_antisense[is_neg] = pos_arr[is_neg].astype(np.float32, copy=False)
    if preserve_ineligible_total:
        n_total = (pos_arr + neg_arr).astype(np.float32, copy=False)
    else:
        n_total = (k_sense + k_antisense).astype(np.float32, copy=False)
    eligible = (is_pos | is_neg).astype(bool, copy=False)
    return k_sense, k_antisense, n_total, eligible


def _compartment_pos_neg(
    payload_arrays: PayloadArrays,
    compartment: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return POS/NEG unspliced mass for one compartment."""
    rc = payload_arrays.region_counts_sorted
    pos = rc[:, channel_index(compartment, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)].astype(
        np.float64, copy=False
    )
    neg = rc[:, channel_index(compartment, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)].astype(
        np.float64, copy=False
    )
    return pos, neg


def build_strand_region_counts(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    p_r1_sense: float,
) -> StrandRegionCounts:
    """Fold POS/NEG channel counts onto sense/antisense per the region's ts_class."""
    if not 0.0 <= float(p_r1_sense) <= 1.0:
        raise ValueError(
            f"build_strand_region_counts: p_r1_sense must be in [0, 1]; got {p_r1_sense!r}"
        )

    pos, neg = _aggregate_pos_neg_all_compartments(payload_arrays)
    k_sense, k_antisense, n_total, strand_eligible = _fold_pos_neg_by_transcript_strand(
        pos, neg, region_arrays.ts_class
    )
    eligible = strand_eligible & (n_total > 0.0)

    return StrandRegionCounts(
        k_sense=k_sense,
        k_antisense=k_antisense,
        n_total=n_total,
        eligible=eligible.astype(bool, copy=False),
        p_r1_sense=float(p_r1_sense),
    )


def build_compartment_strand_counts(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    p_r1_sense: float,
) -> CompartmentStrandCounts:
    """Fold unspliced POS/NEG counts by contained/left-boundary/right-boundary slots."""
    if not 0.0 <= float(p_r1_sense) <= 1.0:
        raise ValueError(
            "build_compartment_strand_counts: p_r1_sense must be in [0, 1]; "
            f"got {p_r1_sense!r}"
        )

    contained_pos, contained_neg = _compartment_pos_neg(payload_arrays, COMPARTMENT_CONTAINED)
    left_pos, left_neg = _compartment_pos_neg(payload_arrays, COMPARTMENT_BOUNDARY_LEFT)
    right_pos, right_neg = _compartment_pos_neg(payload_arrays, COMPARTMENT_BOUNDARY_RIGHT)

    contained_sense, contained_anti, contained_total, eligible = _fold_pos_neg_by_transcript_strand(
        contained_pos,
        contained_neg,
        region_arrays.ts_class,
        preserve_ineligible_total=True,
    )
    left_sense, left_anti, left_total, _ = _fold_pos_neg_by_transcript_strand(
        left_pos,
        left_neg,
        region_arrays.ts_class,
        preserve_ineligible_total=True,
    )
    right_sense, right_anti, right_total, _ = _fold_pos_neg_by_transcript_strand(
        right_pos,
        right_neg,
        region_arrays.ts_class,
        preserve_ineligible_total=True,
    )

    return CompartmentStrandCounts(
        contained_sense=contained_sense,
        contained_antisense=contained_anti,
        contained_total=contained_total,
        boundary_left_sense=left_sense,
        boundary_left_antisense=left_anti,
        boundary_left_total=left_total,
        boundary_right_sense=right_sense,
        boundary_right_antisense=right_anti,
        boundary_right_total=right_total,
        eligible=eligible.astype(bool, copy=False),
        p_r1_sense=float(p_r1_sense),
    )


# ---------------------------------------------------------------------------
# Posterior helpers
# ---------------------------------------------------------------------------


def _round_n(n_total: float) -> int:
    """Round float counts to the nearest non-negative integer for the exact path."""
    return max(0, int(round(float(n_total))))


def _log_binom_pmf(j: np.ndarray, n: int, p: float) -> np.ndarray:
    """Log PMF of Binomial(n, p) at integer points ``j`` (vectorised)."""
    if not np.isfinite(p) or not 0.0 <= float(p) <= 1.0:
        raise ValueError(f"_log_binom_pmf: p must be finite and in [0, 1]; got {p!r}")
    j_arr = np.asarray(j)
    if n <= 0:
        return np.where(j_arr == 0, 0.0, -np.inf)
    if p == 0.0:
        return np.where(j_arr == 0, 0.0, -np.inf)
    if p == 1.0:
        return np.where(j_arr == n, 0.0, -np.inf)
    log_c = gammaln(n + 1) - gammaln(j_arr + 1) - gammaln(n - j_arr + 1)
    return log_c + j_arr * np.log(p) + (n - j_arr) * np.log1p(-p)


def _log_beta_binom_pmf(k: np.ndarray, n: int, alpha: float, beta: float) -> np.ndarray:
    """Log PMF of BetaBinomial(n, alpha, beta) at integer points ``k``."""
    if not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError(f"_log_beta_binom_pmf: alpha must be finite and positive; got {alpha!r}")
    if not np.isfinite(beta) or beta <= 0.0:
        raise ValueError(f"_log_beta_binom_pmf: beta must be finite and positive; got {beta!r}")
    if n <= 0:
        return np.where(k == 0, 0.0, -np.inf)
    log_c = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return log_c + betaln(k + alpha, n - k + beta) - betaln(alpha, beta)


def strand_log_likelihood_d_grid(
    k_sense_obs: int,
    n_total: int,
    d_grid: np.ndarray,
    *,
    kappa_d: float,
    p_r1_sense: float,
) -> np.ndarray:
    """Return ``log L_strand(K | D=d, N)`` for integer ``d_grid``.

    ``D`` is the gDNA count and ``R = N - D`` is the RNA count. RNA sense
    observations follow ``Binom(R, p_r1_sense)`` and gDNA sense observations
    follow the symmetric beta-binomial strand-balance model parameterized by
    ``kappa_d``.
    """
    n = int(n_total)
    if n < 0:
        raise ValueError(f"strand_log_likelihood_d_grid: n_total must be >= 0; got {n}")
    if not np.isfinite(kappa_d) or kappa_d <= 0.0:
        raise ValueError(
            f"strand_log_likelihood_d_grid: kappa_d must be finite and positive; got {kappa_d!r}"
        )
    if not np.isfinite(p_r1_sense) or not 0.0 <= float(p_r1_sense) <= 1.0:
        raise ValueError(
            "strand_log_likelihood_d_grid: p_r1_sense must be finite and in [0, 1]; "
            f"got {p_r1_sense!r}"
        )

    k_obs = int(max(0, min(n, int(k_sense_obs))))
    d_arr = np.rint(np.asarray(d_grid, dtype=np.float64)).astype(np.int64)
    out = np.full(d_arr.shape, -np.inf, dtype=np.float64)
    a = float(kappa_d) / 2.0

    for pos, d_raw in np.ndenumerate(d_arr):
        d = int(d_raw)
        if d < 0 or d > n:
            continue
        r = n - d
        j_lo = max(0, k_obs - d)
        j_hi = min(r, k_obs)
        if j_lo > j_hi:
            continue
        js = np.arange(j_lo, j_hi + 1)
        terms = _log_binom_pmf(js, r, float(p_r1_sense)) + _log_beta_binom_pmf(
            k_obs - js,
            d,
            a,
            a,
        )
        finite_terms = terms[np.isfinite(terms)]
        if finite_terms.size == 0:
            continue
        m = float(np.max(finite_terms))
        out[pos] = m + np.log(np.sum(np.exp(finite_terms - m)))

    return out


def _exact_posterior_R(
    k_sense_obs: int,
    n: int,
    *,
    kappa_d: float,
    p_r1_sense: float,
) -> np.ndarray:
    """Exact posterior PMF over R ∈ {0,...,n} with uniform prior on R.

    Likelihood:
        P(k_sense | R=r, D=n-r)
          = sum_{j=max(0, k_sense-(n-r))..min(r, k_sense)}
              Binom(j | r, p_r1_sense) * BB(k_sense-j | n-r, kappa_d/2, kappa_d/2)
    """
    if n < 0:
        raise ValueError(f"_exact_posterior_R: n must be >= 0; got {n}")
    if k_sense_obs < 0 or k_sense_obs > n:
        # Snap out-of-range observed counts (can occur from float rounding).
        k_sense_obs = int(max(0, min(n, k_sense_obs)))

    r_grid = np.arange(n + 1, dtype=np.int64)
    d_for_r = n - r_grid
    log_post = strand_log_likelihood_d_grid(
        k_sense_obs,
        n,
        d_for_r,
        kappa_d=kappa_d,
        p_r1_sense=p_r1_sense,
    )

    finite = np.isfinite(log_post)
    if not np.any(finite):
        # Degenerate (should not happen). Return a flat prior.
        out = np.ones(n + 1, dtype=np.float64) / float(n + 1)
        return out
    m = float(np.max(log_post[finite]))
    post = np.exp(log_post - m)
    post[~finite] = 0.0
    z = float(post.sum())
    if z <= 0.0:
        out = np.ones(n + 1, dtype=np.float64) / float(n + 1)
        return out
    return post / z


def _summarize_exact(
    posterior: np.ndarray,
    *,
    rna_lower_confidence: float,
    n_int: int,
) -> tuple[float, float, float, float]:
    """Return (mean_R, R_lower_count, R_upper_count, sd_R) from a discrete PMF.

    ``R_lower_count`` = largest integer L with P(R >= L) >= c.
    ``R_upper_count`` = smallest integer U with P(R <= U) >= c.
    """
    r_grid = np.arange(n_int + 1, dtype=np.float64)
    mean_r = float(np.sum(r_grid * posterior))
    var_r = float(np.sum((r_grid - mean_r) ** 2 * posterior))
    sd_r = float(np.sqrt(max(var_r, 0.0)))

    # Lower bound on R at confidence c: largest L s.t. P(R >= L) >= c.
    surv = np.cumsum(posterior[::-1])[::-1]  # surv[L] = P(R >= L)
    L_candidates = np.where(surv >= rna_lower_confidence)[0]
    r_lower = float(L_candidates.max()) if L_candidates.size > 0 else 0.0

    # Upper bound on R at confidence c: smallest U s.t. P(R <= U) >= c.
    cdf = np.cumsum(posterior)
    U_candidates = np.where(cdf >= rna_lower_confidence)[0]
    r_upper = float(U_candidates.min()) if U_candidates.size > 0 else float(n_int)

    return mean_r, r_lower, r_upper, sd_r


def _summarize_normal(
    k_sense: float,
    n: float,
    *,
    kappa_d: float,
    p_r1_sense: float,
    rna_lower_confidence: float,
) -> tuple[float, float, float, float]:
    """Normal-approximation summary of the posterior over R given (k_sense, n).

    Treats R as continuous on [0, n] and linearises k_sense ≈ mu(R) with
    mu(R) = 0.5 n + R (p − 0.5). Variance is evaluated at R̂ to avoid
    needing to invert the heteroscedastic relationship.
    """
    delta = float(p_r1_sense) - 0.5
    if abs(delta) < STRAND_CONTRAST_NUMERICAL_FLOOR:
        # Identifiability guard — caller has already set FLAG_NEAR_UNSTRANDED.
        return 0.0, 0.0, float(n), 0.0

    r_hat = (k_sense - 0.5 * n) / delta
    r_hat = float(np.clip(r_hat, 0.0, n))
    d_hat = float(max(0.0, n - r_hat))

    var_rna = r_hat * p_r1_sense * (1.0 - p_r1_sense)
    if 1.0 + kappa_d > 0.0:
        var_dna = 0.25 * d_hat * (d_hat + kappa_d) / (1.0 + kappa_d)
    else:
        var_dna = 0.0
    var_k = max(var_rna + var_dna, 0.0)
    var_r = var_k / (delta * delta)
    sd_r = float(np.sqrt(var_r))

    z = float(norm.ppf(float(rna_lower_confidence)))
    r_lower = float(np.clip(r_hat - z * sd_r, 0.0, n))
    r_upper = float(np.clip(r_hat + z * sd_r, 0.0, n))
    return r_hat, r_lower, r_upper, sd_r


# ---------------------------------------------------------------------------
# Deconvolution
# ---------------------------------------------------------------------------


def deconvolve_regions_by_strand(
    counts: StrandRegionCounts,
    *,
    kappa_d: float,
    rna_lower_confidence: float,
    kappa_d_n_seed_regions: int = 0,
    kappa_d_n_exon_self_training: int = 0,
    kappa_d_fallback_used: bool = False,
) -> RegionGdnaEstimate:
    """Per-region gDNA fragment-count posterior given ``kappa_d`` and ``p_r1_sense``."""
    if not 0.5 <= float(rna_lower_confidence) < 1.0:
        raise ValueError(
            "deconvolve_regions_by_strand: rna_lower_confidence must be in [0.5, 1.0); "
            f"got {rna_lower_confidence!r}"
        )
    if not np.isfinite(kappa_d) or kappa_d <= 0.0:
        raise ValueError(
            f"deconvolve_regions_by_strand: kappa_d must be finite and positive; got {kappa_d!r}"
        )

    R = int(counts.k_sense.shape[0])
    n_total = counts.n_total.astype(np.float32, copy=True)
    mean_count = np.zeros(R, dtype=np.float32)
    upper_count = np.zeros(R, dtype=np.float32)
    rna_lower_count = np.zeros(R, dtype=np.float32)
    precision = np.zeros(R, dtype=np.float32)
    flags = np.zeros(R, dtype=np.uint8)

    p_r1_sense = float(counts.p_r1_sense)
    if not np.isfinite(p_r1_sense) or not 0.0 <= p_r1_sense <= 1.0:
        raise ValueError(
            "deconvolve_regions_by_strand: counts.p_r1_sense must be finite and in [0, 1]; "
            f"got {counts.p_r1_sense!r}"
        )
    near_unstranded = abs(p_r1_sense - 0.5) < STRAND_CONTRAST_NUMERICAL_FLOOR

    eligible = np.asarray(counts.eligible, dtype=bool)

    for i in range(R):
        if not eligible[i] or n_total[i] <= 0.0:
            flags[i] |= FLAG_INELIGIBLE
            mean_count[i] = n_total[i]
            upper_count[i] = n_total[i]
            rna_lower_count[i] = 0.0
            precision[i] = 0.0
            if kappa_d_fallback_used:
                flags[i] |= FLAG_KAPPA_FALLBACK
            continue

        if near_unstranded:
            flags[i] |= FLAG_NEAR_UNSTRANDED
            mean_count[i] = n_total[i]
            upper_count[i] = n_total[i]
            rna_lower_count[i] = 0.0
            precision[i] = 0.0
            if kappa_d_fallback_used:
                flags[i] |= FLAG_KAPPA_FALLBACK
            continue

        n_int = _round_n(float(n_total[i]))
        k_int = int(max(0, min(n_int, round(float(counts.k_sense[i])))))

        if n_int <= MAX_EXACT_POSTERIOR_N:
            posterior = _exact_posterior_R(k_int, n_int, kappa_d=kappa_d, p_r1_sense=p_r1_sense)
            r_hat, r_lower, r_upper, sd_r = _summarize_exact(
                posterior, rna_lower_confidence=rna_lower_confidence, n_int=n_int
            )
        else:
            flags[i] |= FLAG_APPROX_NORMAL
            r_hat, r_lower, r_upper, sd_r = _summarize_normal(
                float(counts.k_sense[i]),
                float(n_total[i]),
                kappa_d=kappa_d,
                p_r1_sense=p_r1_sense,
                rna_lower_confidence=rna_lower_confidence,
            )

        n_val = float(n_total[i])
        r_hat_out = float(np.clip(r_hat, 0.0, n_val))
        r_lower_out = float(np.clip(r_lower, 0.0, n_val))
        mean_count[i] = n_val - r_hat_out
        upper_count[i] = n_val - r_lower_out
        rna_lower_count[i] = r_lower_out
        # Precision in [0, 1]: 1 ⇔ posterior SD is zero, 0 ⇔ SD ≥ n/2.
        if n_val > 0.0:
            precision[i] = float(max(0.0, min(1.0, 1.0 - 2.0 * sd_r / n_val)))
        else:
            precision[i] = 0.0

        if kappa_d_fallback_used:
            flags[i] |= FLAG_KAPPA_FALLBACK

    return RegionGdnaEstimate(
        n_total=n_total,
        mean_count=mean_count,
        upper_count=upper_count,
        rna_lower_count=rna_lower_count,
        precision=precision,
        flags=flags,
        kappa_d=float(kappa_d),
        kappa_d_n_seed_regions=int(kappa_d_n_seed_regions),
        kappa_d_n_exon_self_training=int(kappa_d_n_exon_self_training),
        p_r1_sense=p_r1_sense,
        rna_lower_confidence=float(rna_lower_confidence),
    )


def _single_compartment_estimate(
    counts: CompartmentStrandCounts,
    *,
    k_sense: np.ndarray,
    k_antisense: np.ndarray,
    n_total: np.ndarray,
    kappa_d: float,
    rna_lower_confidence: float,
) -> RegionGdnaEstimate:
    return deconvolve_regions_by_strand(
        StrandRegionCounts(
            k_sense=k_sense.astype(np.float32, copy=False),
            k_antisense=k_antisense.astype(np.float32, copy=False),
            n_total=n_total.astype(np.float32, copy=False),
            eligible=(counts.eligible & (n_total > 0.0)).astype(bool, copy=False),
            p_r1_sense=counts.p_r1_sense,
        ),
        kappa_d=kappa_d,
        rna_lower_confidence=rna_lower_confidence,
    )


def deconvolve_compartments_by_strand(
    counts: CompartmentStrandCounts,
    *,
    kappa_d: float,
    rna_lower_confidence: float,
) -> RegionGdnaChannelEstimate:
    """Run the strand deconvolution independently for each unspliced compartment."""
    contained = _single_compartment_estimate(
        counts,
        k_sense=counts.contained_sense,
        k_antisense=counts.contained_antisense,
        n_total=counts.contained_total,
        kappa_d=kappa_d,
        rna_lower_confidence=rna_lower_confidence,
    )
    boundary_left = _single_compartment_estimate(
        counts,
        k_sense=counts.boundary_left_sense,
        k_antisense=counts.boundary_left_antisense,
        n_total=counts.boundary_left_total,
        kappa_d=kappa_d,
        rna_lower_confidence=rna_lower_confidence,
    )
    boundary_right = _single_compartment_estimate(
        counts,
        k_sense=counts.boundary_right_sense,
        k_antisense=counts.boundary_right_antisense,
        n_total=counts.boundary_right_total,
        kappa_d=kappa_d,
        rna_lower_confidence=rna_lower_confidence,
    )
    flags = (
        contained.flags.astype(np.uint16)
        | boundary_left.flags.astype(np.uint16)
        | boundary_right.flags.astype(np.uint16)
    )
    return RegionGdnaChannelEstimate(
        contained_mean=contained.mean_count,
        contained_upper=contained.upper_count,
        contained_rna_lower=contained.rna_lower_count,
        contained_precision=contained.precision,
        boundary_left_mean=boundary_left.mean_count,
        boundary_left_upper=boundary_left.upper_count,
        boundary_left_rna_lower=boundary_left.rna_lower_count,
        boundary_left_precision=boundary_left.precision,
        boundary_right_mean=boundary_right.mean_count,
        boundary_right_upper=boundary_right.upper_count,
        boundary_right_rna_lower=boundary_right.rna_lower_count,
        boundary_right_precision=boundary_right.precision,
        flags=flags,
        kappa_d=float(kappa_d),
        p_r1_sense=float(counts.p_r1_sense),
        rna_lower_confidence=float(rna_lower_confidence),
    )


# ---------------------------------------------------------------------------
# Exon self-training screen
# ---------------------------------------------------------------------------


def _has_spliced_mass(payload_arrays: PayloadArrays) -> np.ndarray:
    """Per-region boolean mask: any SPLICE_SPLICED mass across compartments/strands."""
    rc = payload_arrays.region_counts_sorted
    spliced = np.zeros(rc.shape[0], dtype=np.float64)
    for compartment in (
        COMPARTMENT_CONTAINED,
        COMPARTMENT_BOUNDARY_LEFT,
        COMPARTMENT_BOUNDARY_RIGHT,
    ):
        for strand in (CHANNEL_STRAND_POS, CHANNEL_STRAND_NEG):
            spliced += rc[:, channel_index(compartment, SPLICE_SPLICED, strand)].astype(
                np.float64, copy=False
            )
    return spliced > 0.0


def screen_no_rna_exons(
    counts: StrandRegionCounts,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    kappa_d_seed: float,
    *,
    rna_lower_confidence: float,
) -> np.ndarray:
    """Identify exon regions consistent with pure gDNA at the seed ``kappa_d``.

    A region is accepted iff:
      * its signature contains any exon bit;
      * ``ts_class`` is informative (TS_POS or TS_NEG);
      * both strand channels have observed mass (``k_sense > 0`` and
        ``k_antisense > 0``);
      * the region has no SPLICE_SPLICED mass;
      * the two-sided BetaBinomial p-value of ``k_sense`` under the pure-gDNA
        null ``k_sense ~ BB(n, kappa_d/2, kappa_d/2)`` exceeds the threshold
        ``1 - rna_lower_confidence`` (i.e., the observation is not rejected
        as inconsistent with H₀: R = 0 at the requested confidence level).

    The original Phase 3 plan phrased the test as
    ``P(R > 0 | data, kappa_d_seed) <= 1 - c``. Implemented literally with a
    flat prior on ``R ∈ {0,..,n}`` that singleton posterior is dilute and the
    criterion is essentially unreachable for finite ``n``. We instead use the
    frequentist tail test under H₀, which is the standard reading of "data
    consistent with pure gDNA" and matches the intent of the spec. See the
    Phase 3 progress log for details.

    Returns a boolean mask over the sorted region table.
    """
    if not 0.5 <= float(rna_lower_confidence) < 1.0:
        raise ValueError(
            "screen_no_rna_exons: rna_lower_confidence must be in [0.5, 1.0); "
            f"got {rna_lower_confidence!r}"
        )

    p_r1_sense = float(counts.p_r1_sense)
    if not np.isfinite(p_r1_sense) or not 0.0 <= p_r1_sense <= 1.0:
        raise ValueError(
            "screen_no_rna_exons: counts.p_r1_sense must be finite and in [0, 1]; "
            f"got {counts.p_r1_sense!r}"
        )
    if abs(p_r1_sense - 0.5) < STRAND_CONTRAST_NUMERICAL_FLOOR:
        return np.zeros(region_arrays.ts_class.shape, dtype=bool)

    is_exon = is_exon_any(region_arrays.signature)
    has_strand = (region_arrays.ts_class == TS_POS) | (region_arrays.ts_class == TS_NEG)
    has_both = (counts.k_sense > 0.0) & (counts.k_antisense > 0.0)
    no_spliced = ~_has_spliced_mass(payload_arrays)

    candidate = is_exon & has_strand & has_both & no_spliced
    accepted = np.zeros(candidate.shape, dtype=bool)
    if not np.any(candidate):
        return accepted

    if not np.isfinite(kappa_d_seed) or kappa_d_seed <= 0.0:
        raise ValueError(
            f"screen_no_rna_exons: kappa_d_seed must be finite and positive; got {kappa_d_seed!r}"
        )
    alpha = float(kappa_d_seed) / 2.0
    threshold = 1.0 - float(rna_lower_confidence)

    idxs = np.where(candidate)[0]
    for i in idxs:
        n_int = _round_n(float(counts.n_total[i]))
        if n_int <= 0:
            continue
        k_int = int(max(0, min(n_int, round(float(counts.k_sense[i])))))
        ks = np.arange(n_int + 1)
        log_pmf = _log_beta_binom_pmf(ks, n_int, alpha, alpha)
        log_p_obs = float(log_pmf[k_int])
        tail = log_pmf <= log_p_obs + 1e-12
        if not np.any(tail):
            continue
        # Normalise (log_pmf sums to 0 already, but be safe against rounding).
        m_all = float(np.max(log_pmf))
        z = m_all + np.log(np.sum(np.exp(log_pmf - m_all)))
        m_tail = float(np.max(log_pmf[tail]))
        p_value = float(np.exp(m_tail + np.log(np.sum(np.exp(log_pmf[tail] - m_tail))) - z))
        if p_value >= threshold:
            accepted[i] = True

    return accepted


# ---------------------------------------------------------------------------
# kappa_d refit
# ---------------------------------------------------------------------------


_HARD_FALLBACK_REASONS: frozenset[str] = frozenset(
    {
        "n_regions < MIN_REGIONS_FOR_STRAND_MOM",
        "no positive strand exposure",
    }
)


def _is_hard_fallback(estimate: StrandBalanceEstimate) -> bool:
    """True when the MoM estimate is uninformative (no data, not just tight data)."""
    reason = str(estimate.fallback_reason).split(";", maxsplit=1)[0].strip()
    return bool(estimate.fallback_used) and reason in _HARD_FALLBACK_REASONS


def _refit_with_fallback(
    pos: np.ndarray,
    neg: np.ndarray,
    mask: np.ndarray,
) -> StrandBalanceEstimate:
    """Run MoM strand-balance estimate; substitute uniform prior on *hard* fallback only.

    Soft fallbacks (e.g. "residual variance <= binomial expectation") keep the
    MoM kappa estimate — they signal data that is *tighter* than Binomial(0.5),
    which is information, not uncertainty.
    """
    estimate = estimate_strand_balance(pos, neg, mask)
    if not _is_hard_fallback(estimate) and np.isfinite(estimate.kappa) and estimate.kappa > 0.0:
        return estimate
    # Replace the kappa with the uninformative Beta(1, 1) prior. The other
    # fields (n_regions, residual_sum, …) are preserved so diagnostics still
    # show what the orchestrator saw.
    return StrandBalanceEstimate(
        kappa=DNA_STRAND_UNIFORM_PRIOR,
        n_regions=estimate.n_regions,
        n_fragments=estimate.n_fragments,
        n_pos=estimate.n_pos,
        n_neg=estimate.n_neg,
        observed_pos_fraction=estimate.observed_pos_fraction,
        residual_sum=estimate.residual_sum,
        binomial_variance_sum=estimate.binomial_variance_sum,
        max_overdispersed_variance_sum=estimate.max_overdispersed_variance_sum,
        overdispersion_factor=estimate.overdispersion_factor,
        fallback_used=True,
        fallback_reason=(estimate.fallback_reason or "kappa unidentifiable")
        + "; using DNA_STRAND_UNIFORM_PRIOR",
    )


def estimate_kappa_d(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    counts: StrandRegionCounts,
    strand_summary: StrandSummary,
    *,
    rna_lower_confidence: float,
) -> KappaDEstimate:
    """Estimate the gDNA strand-balance overdispersion ``kappa_d``.

    Procedure:
      1. Seed from high-purity regions (intergenic + intron-only) using
         the existing contained-unspliced POS/NEG counts.
      2. Run :func:`screen_no_rna_exons` with the seed ``kappa_d`` to
         self-identify exon regions that look like pure gDNA.
      3. Refit MoM on the union (seed + accepted exons), using the
         folded ``k_sense / k_antisense`` from ``counts`` for the exon
         additions.

    ``strand_summary.p_r1_sense`` must match ``counts.p_r1_sense`` so Phase 4
    wiring cannot silently mix summaries from different strand-model fits.
    """
    if abs(float(strand_summary.p_r1_sense) - float(counts.p_r1_sense)) > 1.0e-9:
        raise ValueError(
            "estimate_kappa_d: strand_summary.p_r1_sense does not match "
            f"counts.p_r1_sense ({strand_summary.p_r1_sense!r} != {counts.p_r1_sense!r})."
        )

    seed_mask = is_intergenic(region_arrays.signature) | is_intron_only(region_arrays.signature)
    seed_pos = payload_arrays.contained_unspliced_pos
    seed_neg = payload_arrays.contained_unspliced_neg

    seed_estimate = _refit_with_fallback(seed_pos, seed_neg, seed_mask)
    n_seed = int(((seed_pos + seed_neg) > 0.0)[seed_mask].sum())

    # If the seed itself is a *hard* fallback (no data) the protocol is
    # unidentifiable; skip self-training. Soft fallbacks (data tighter than
    # Binomial) are kept and used to seed the screen.
    if _is_hard_fallback(seed_estimate):
        return KappaDEstimate(
            balance=seed_estimate,
            n_seed_regions=n_seed,
            n_exon_self_training=0,
        )

    # Self-training: accept exons whose strand pattern is consistent with pure gDNA.
    accepted_exons = screen_no_rna_exons(
        counts,
        region_arrays,
        payload_arrays,
        seed_estimate.kappa,
        rna_lower_confidence=rna_lower_confidence,
    )
    n_exons = int(accepted_exons.sum())

    if n_exons == 0:
        # No exons to add — keep the seed estimate.
        return KappaDEstimate(
            balance=seed_estimate,
            n_seed_regions=n_seed,
            n_exon_self_training=0,
        )

    # Build a per-region POS/NEG vector that mixes the seed channels with
    # the folded counts on accepted exons. For TS_POS exons sense == POS,
    # for TS_NEG exons sense == NEG; we recover POS/NEG by inverting the fold.
    pos_combined = np.array(seed_pos, dtype=np.float32, copy=True)
    neg_combined = np.array(seed_neg, dtype=np.float32, copy=True)

    is_pos = region_arrays.ts_class == TS_POS
    is_neg = region_arrays.ts_class == TS_NEG
    add_pos_exons = accepted_exons & is_pos
    add_neg_exons = accepted_exons & is_neg
    pos_combined[add_pos_exons] = counts.k_sense[add_pos_exons]
    neg_combined[add_pos_exons] = counts.k_antisense[add_pos_exons]
    pos_combined[add_neg_exons] = counts.k_antisense[add_neg_exons]
    neg_combined[add_neg_exons] = counts.k_sense[add_neg_exons]

    refit_mask = seed_mask | accepted_exons
    refit_estimate = _refit_with_fallback(pos_combined, neg_combined, refit_mask)

    return KappaDEstimate(
        balance=refit_estimate,
        n_seed_regions=n_seed,
        n_exon_self_training=n_exons,
    )
