"""Fusion of regional density evidence with local strand evidence."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import truncnorm

from .density_model import density_logpmf_grid
from .strand_deconv import (
    MAX_EXACT_POSTERIOR_N,
    _approx_log_likelihood_minor_beta_binom,
    _minor_count_from_folded,
    _posterior_weights_from_log,
    _uniform_grid_log_weights,
    strand_log_likelihood_d_grid_minor_beta_binom,
)
from .strand_summary import STRAND_CONTRAST_NUMERICAL_FLOOR

if TYPE_CHECKING:  # pragma: no cover - annotations only.
    from ._arrays import RegionArrays
    from .density_model import DensityEvidence
    from .density_observation import DensityObservation
    from .region_count_ledger import RegionCountLedger
    from .strand_deconv import StrandRegionCounts
    from .strand_summary import StrandSummary


__all__ = [
    "ADAPTIVE_EXACT_MAX",
    "FUSED_APPROX",
    "FUSED_BOUNDARY_FALLBACK",
    "FUSED_DEGENERATE",
    "FUSED_DENSITY_ONLY",
    "FUSED_DENSITY_TAIL",
    "FUSED_EXACT",
    "FUSED_FALLBACK",
    "FUSED_STRAND_USED",
    "FusedRegionGdnaEvidence",
    "fuse_density_and_strand",
]


FUSED_DENSITY_ONLY: int = 1 << 0
FUSED_STRAND_USED: int = 1 << 1
FUSED_EXACT: int = 1 << 2
FUSED_APPROX: int = 1 << 3
FUSED_DENSITY_TAIL: int = 1 << 4
FUSED_DEGENERATE: int = 1 << 5
FUSED_FALLBACK: int = 1 << 6
FUSED_BOUNDARY_FALLBACK: int = 1 << 7

ADAPTIVE_EXACT_MAX: int = 1000
_BOUNDARY_WINDOW: int = 200
_MIN_VARIANCE: float = 1.0e-6


@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    """Per-region bounded posterior over local gDNA count."""

    mean_count: np.ndarray
    upper_count: np.ndarray
    variance_count: np.ndarray
    rna_lower_count: np.ndarray
    observed_compatible_count: np.ndarray
    density_weight: np.ndarray
    strand_weight: np.ndarray
    density_applicable: np.ndarray
    strand_applicable: np.ndarray
    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    flags: np.ndarray

    def to_summary_dict(self) -> dict[str, object]:
        """Return compact JSON-safe fusion diagnostics."""
        flags = np.asarray(self.flags, dtype=np.uint8)
        return {
            "n_regions": int(flags.size),
            "mean_count": _summary_stats(self.mean_count),
            "upper_count": _summary_stats(self.upper_count),
            "rna_lower_count": _summary_stats(self.rna_lower_count),
            "density_weight": _summary_stats(self.density_weight),
            "strand_weight": _summary_stats(self.strand_weight),
            "tail_probability": _summary_stats(self.tail_probability),
            "expected_tail_count": _summary_stats(self.expected_tail_count),
            "flags": {
                "n_density_only": _flag_count(flags, FUSED_DENSITY_ONLY),
                "n_strand_used": _flag_count(flags, FUSED_STRAND_USED),
                "n_exact": _flag_count(flags, FUSED_EXACT),
                "n_approx": _flag_count(flags, FUSED_APPROX),
                "n_density_tail": _flag_count(flags, FUSED_DENSITY_TAIL),
                "n_degenerate": _flag_count(flags, FUSED_DEGENERATE),
                "n_fallback": _flag_count(flags, FUSED_FALLBACK),
                "n_boundary_fallback": _flag_count(flags, FUSED_BOUNDARY_FALLBACK),
            },
        }


def fuse_density_and_strand(
    *,
    region_arrays: "RegionArrays",
    ledger: "RegionCountLedger",
    density_observation: "DensityObservation",
    density_evidence: "DensityEvidence",
    strand_counts: "StrandRegionCounts",
    strand_summary: "StrandSummary",
    kappa_d: float,
    confidence: float,
) -> FusedRegionGdnaEvidence:
    """Fuse density predictive counts with strand likelihoods region-by-region."""
    conf = float(confidence)
    if not 0.5 <= conf < 1.0:
        raise ValueError(f"fuse_density_and_strand: confidence must be in [0.5, 1.0); got {conf}.")
    if not np.isfinite(kappa_d) or kappa_d <= 0.0:
        raise ValueError(f"fuse_density_and_strand: kappa_d must be positive; got {kappa_d!r}.")

    _ = ledger
    observed = np.asarray(density_observation.observed_compatible_count, dtype=np.float64)
    n_regions = int(observed.size)
    if np.asarray(region_arrays.start).size != n_regions:
        raise ValueError("fuse_density_and_strand: region array and observation sizes differ.")

    mean_count = np.zeros(n_regions, dtype=np.float32)
    upper_count = np.zeros(n_regions, dtype=np.float32)
    variance_count = np.zeros(n_regions, dtype=np.float32)
    rna_lower_count = np.zeros(n_regions, dtype=np.float32)
    density_weight = np.zeros(n_regions, dtype=np.float32)
    strand_weight = np.zeros(n_regions, dtype=np.float32)
    density_applicable = np.ones(n_regions, dtype=bool)
    strand_applicable = np.zeros(n_regions, dtype=bool)
    flags = np.zeros(n_regions, dtype=np.uint8)

    tail_probability = _density_array_or_zeros(density_evidence.tail_probability, n_regions)
    expected_tail_count = _density_array_or_zeros(density_evidence.expected_tail_count, n_regions)
    flags[tail_probability > 0.0] |= np.uint8(FUSED_DENSITY_TAIL)

    p_r1_sense = float(strand_counts.p_r1_sense)
    if not np.isclose(float(strand_summary.p_r1_sense), p_r1_sense, rtol=0.0, atol=1e-12):
        raise ValueError(
            "fuse_density_and_strand: strand_summary.p_r1_sense must match strand_counts; "
            f"got {strand_summary.p_r1_sense!r} vs {p_r1_sense!r}."
        )
    near_unstranded = abs(p_r1_sense - 0.5) < STRAND_CONTRAST_NUMERICAL_FLOOR
    strand_eligible = np.asarray(strand_counts.eligible, dtype=bool) & (not near_unstranded)

    for region_idx in range(n_regions):
        n_observed = int(max(0, round(float(observed[region_idx]))))
        if n_observed == 0:
            flags[region_idx] |= np.uint8(FUSED_DEGENERATE | FUSED_DENSITY_ONLY)
            density_weight[region_idx] = 1.0
            continue

        n_strand = int(max(n_observed, round(float(strand_counts.n_total[region_idx]))))
        k_sense = int(max(0, min(n_strand, round(float(strand_counts.k_sense[region_idx])))))
        use_strand = bool(strand_eligible[region_idx] and n_strand > 0)
        strand_applicable[region_idx] = use_strand

        if n_observed <= MAX_EXACT_POSTERIOR_N:
            result = _fuse_exact_region(
                density_evidence=density_evidence,
                region_idx=region_idx,
                n_observed=n_observed,
                n_strand=n_strand,
                k_sense=k_sense,
                use_strand=use_strand,
                kappa_d=float(kappa_d),
                p_r1_sense=p_r1_sense,
                strand_summary=strand_summary,
                confidence=conf,
            )
        else:
            result = _fuse_approx_region(
                density_evidence=density_evidence,
                region_idx=region_idx,
                n_observed=n_observed,
                n_strand=n_strand,
                k_sense=k_sense,
                use_strand=use_strand,
                kappa_d=float(kappa_d),
                p_r1_sense=p_r1_sense,
                strand_summary=strand_summary,
                confidence=conf,
            )

        mean_count[region_idx] = np.float32(result[0])
        upper_count[region_idx] = np.float32(result[1])
        variance_count[region_idx] = np.float32(result[2])
        rna_lower_count[region_idx] = np.float32(max(float(n_observed) - result[1], 0.0))
        density_weight[region_idx] = np.float32(result[3])
        strand_weight[region_idx] = np.float32(result[4])
        flags[region_idx] |= np.uint8(result[5])

    return FusedRegionGdnaEvidence(
        mean_count=mean_count,
        upper_count=upper_count,
        variance_count=variance_count,
        rna_lower_count=rna_lower_count,
        observed_compatible_count=observed.astype(np.float32, copy=False),
        density_weight=density_weight,
        strand_weight=strand_weight,
        density_applicable=density_applicable,
        strand_applicable=strand_applicable,
        tail_probability=tail_probability.astype(np.float32, copy=False),
        expected_tail_count=expected_tail_count.astype(np.float32, copy=False),
        flags=flags,
    )


def _fuse_exact_region(
    *,
    density_evidence: "DensityEvidence",
    region_idx: int,
    n_observed: int,
    n_strand: int,
    k_sense: int,
    use_strand: bool,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
    confidence: float,
) -> tuple[float, float, float, float, float, int]:
    d_grid = np.arange(n_observed + 1, dtype=np.int64)
    log_density = density_logpmf_grid(density_evidence, region_idx, d_grid)
    if use_strand:
        log_strand = _strand_beta_binom_log_likelihood_d_grid(
            k_sense,
            n_strand,
            d_grid,
            kappa_d=kappa_d,
            p_r1_sense=p_r1_sense,
            strand_summary=strand_summary,
        )
        base_flag = FUSED_EXACT | FUSED_STRAND_USED
        density_w = strand_w = 0.5
    else:
        log_strand = np.zeros(n_observed + 1, dtype=log_density.dtype)
        base_flag = FUSED_EXACT | FUSED_DENSITY_ONLY
        density_w = 1.0
        strand_w = 0.0

    return (
        *_summarize_grid_posterior(d_grid, log_density + log_strand, confidence),
        density_w,
        strand_w,
        base_flag,
    )


def _fuse_approx_region(
    *,
    density_evidence: "DensityEvidence",
    region_idx: int,
    n_observed: int,
    n_strand: int,
    k_sense: int,
    use_strand: bool,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
    confidence: float,
) -> tuple[float, float, float, float, float, int]:
    density_mean, density_var = _density_normal_params(density_evidence, region_idx, n_observed)
    if density_var <= _MIN_VARIANCE:
        mean = float(np.clip(density_mean, 0.0, n_observed))
        return mean, mean, 0.0, 1.0, 0.0, FUSED_DEGENERATE | FUSED_DENSITY_ONLY

    density_tau = 1.0 / max(density_var, _MIN_VARIANCE)
    strand_tau = 0.0
    strand_mean = 0.0

    if use_strand:
        strand_mean, strand_var = _strand_normal_params(
            k_sense,
            n_strand,
            n_observed,
            kappa_d=kappa_d,
            p_r1_sense=p_r1_sense,
            strand_summary=strand_summary,
        )
        strand_tau = 1.0 / max(strand_var, _MIN_VARIANCE)

    tau = density_tau + strand_tau
    if tau <= 0.0 or not np.isfinite(tau):
        mean = float(np.clip(density_mean, 0.0, n_observed))
        return mean, mean, 0.0, 1.0, 0.0, FUSED_FALLBACK | FUSED_DENSITY_ONLY

    loc = (density_tau * density_mean + strand_tau * strand_mean) / tau
    variance = 1.0 / tau
    sigma = float(np.sqrt(max(variance, _MIN_VARIANCE)))
    d_hat = _bounded_mode(
        loc,
        sigma,
        n_observed,
        density_evidence,
        region_idx,
        use_strand,
        k_sense,
        n_strand,
        kappa_d,
        p_r1_sense,
        strand_summary,
    )

    if d_hat < 5.0 or d_hat > float(n_observed - 5):
        if n_observed <= ADAPTIVE_EXACT_MAX:
            exact = _fuse_exact_region(
                density_evidence=density_evidence,
                region_idx=region_idx,
                n_observed=n_observed,
                n_strand=n_strand,
                k_sense=k_sense,
                use_strand=use_strand,
                kappa_d=kappa_d,
                p_r1_sense=p_r1_sense,
                strand_summary=strand_summary,
                confidence=confidence,
            )
            return (*exact[:5], exact[5] | FUSED_BOUNDARY_FALLBACK)
        boundary = _boundary_window_posterior(
            density_evidence=density_evidence,
            region_idx=region_idx,
            n_observed=n_observed,
            n_strand=n_strand,
            k_sense=k_sense,
            use_strand=use_strand,
            kappa_d=kappa_d,
            p_r1_sense=p_r1_sense,
            strand_summary=strand_summary,
            confidence=confidence,
            upper=d_hat > float(n_observed / 2),
        )
        return (
            *boundary[:3],
            _density_weight(density_tau, strand_tau),
            _strand_weight(density_tau, strand_tau),
            boundary[3],
        )

    if not np.isfinite(loc) or not np.isfinite(sigma) or sigma <= 0.0:
        mean = float(np.clip(loc, 0.0, n_observed))
        flags = FUSED_FALLBACK | (FUSED_STRAND_USED if use_strand else FUSED_DENSITY_ONLY)
        return (
            mean,
            mean,
            0.0,
            _density_weight(density_tau, strand_tau),
            _strand_weight(density_tau, strand_tau),
            flags,
        )

    a = (0.0 - loc) / sigma
    b = (float(n_observed) - loc) / sigma
    try:
        dist = truncnorm(a, b, loc=loc, scale=sigma)
        mean = float(dist.mean())
        var = float(dist.var())
        upper = float(dist.ppf(confidence))
    except (FloatingPointError, ValueError):
        mean = float(np.clip(loc, 0.0, n_observed))
        var = 0.0
        upper = mean
        flag = FUSED_FALLBACK
    else:
        flag = FUSED_APPROX

    if not np.isfinite(mean) or not np.isfinite(var) or not np.isfinite(upper):
        mean = float(np.clip(loc, 0.0, n_observed))
        var = 0.0
        upper = mean
        flag = FUSED_FALLBACK

    flag |= FUSED_STRAND_USED if use_strand else FUSED_DENSITY_ONLY
    return (
        float(np.clip(mean, 0.0, n_observed)),
        float(np.clip(upper, 0.0, n_observed)),
        float(max(var, 0.0)),
        _density_weight(density_tau, strand_tau),
        _strand_weight(density_tau, strand_tau),
        flag,
    )


def _bounded_mode(
    loc: float,
    sigma: float,
    n_observed: int,
    density_evidence: "DensityEvidence",
    region_idx: int,
    use_strand: bool,
    k_sense: int,
    n_strand: int,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
) -> float:
    if not np.isfinite(loc) or not np.isfinite(sigma) or sigma <= 0.0:
        return float(np.clip(loc, 0.0, n_observed))

    def objective(d_value: float) -> float:
        d_int = np.array([int(round(float(np.clip(d_value, 0.0, n_observed))))])
        log_density = density_logpmf_grid(density_evidence, region_idx, d_int)[0]
        if use_strand:
            log_strand = _strand_beta_binom_log_likelihood_d_grid(
                k_sense,
                n_strand,
                d_int,
                kappa_d=kappa_d,
                p_r1_sense=p_r1_sense,
                strand_summary=strand_summary,
            )[0]
        else:
            log_strand = 0.0
        value = log_density + log_strand
        return -float(value) if np.isfinite(value) else np.inf

    try:
        opt = minimize_scalar(objective, bounds=(0.0, float(n_observed)), method="bounded")
    except ValueError:
        return float(np.clip(loc, 0.0, n_observed))
    if not opt.success or not np.isfinite(opt.x):
        return float(np.clip(loc, 0.0, n_observed))
    return float(np.clip(opt.x, 0.0, n_observed))


def _boundary_window_posterior(
    *,
    density_evidence: "DensityEvidence",
    region_idx: int,
    n_observed: int,
    n_strand: int,
    k_sense: int,
    use_strand: bool,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
    confidence: float,
    upper: bool,
) -> tuple[float, float, float, int]:
    if upper:
        d_grid = np.arange(max(0, n_observed - _BOUNDARY_WINDOW), n_observed + 1, dtype=np.int64)
    else:
        d_grid = np.arange(0, min(n_observed, _BOUNDARY_WINDOW) + 1, dtype=np.int64)
    log_density = density_logpmf_grid(density_evidence, region_idx, d_grid)
    if use_strand:
        log_strand = _strand_beta_binom_log_likelihood_d_grid(
            k_sense,
            n_strand,
            d_grid,
            kappa_d=kappa_d,
            p_r1_sense=p_r1_sense,
            strand_summary=strand_summary,
        )
        flag = FUSED_BOUNDARY_FALLBACK | FUSED_STRAND_USED
    else:
        log_strand = np.zeros(d_grid.shape, dtype=log_density.dtype)
        flag = FUSED_BOUNDARY_FALLBACK | FUSED_DENSITY_ONLY
    mean, upper_count, var = _summarize_grid_posterior(d_grid, log_density + log_strand, confidence)
    return mean, upper_count, var, flag


def _summarize_grid_posterior(
    d_grid: np.ndarray,
    log_post: np.ndarray,
    confidence: float,
) -> tuple[float, float, float]:
    finite = np.isfinite(log_post)
    if not np.any(finite):
        mean = float(np.mean(d_grid)) if d_grid.size else 0.0
        return mean, mean, 0.0
    log_norm = float(logsumexp(log_post[finite]))
    probs = np.zeros(d_grid.shape, dtype=np.float64)
    probs[finite] = np.exp(log_post[finite] - log_norm)
    total = float(probs.sum())
    if total <= 0.0 or not np.isfinite(total):
        mean = float(np.mean(d_grid)) if d_grid.size else 0.0
        return mean, mean, 0.0
    probs /= total
    d_float = d_grid.astype(np.float64, copy=False)
    mean = float(np.sum(d_float * probs))
    var = float(np.sum((d_float - mean) ** 2 * probs))
    cdf = np.cumsum(probs)
    q_idx = int(np.searchsorted(cdf, float(confidence), side="left"))
    q_idx = min(max(q_idx, 0), d_grid.size - 1)
    return mean, float(d_grid[q_idx]), max(var, 0.0)


def _strand_beta_binom_log_likelihood_d_grid(
    k_sense: int,
    n_strand: int,
    d_grid: np.ndarray,
    *,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
) -> np.ndarray:
    k_minor = _minor_count_from_folded(k_sense, n_strand, p_r1_sense)
    return strand_log_likelihood_d_grid_minor_beta_binom(
        k_minor,
        n_strand,
        d_grid,
        kappa_d=kappa_d,
        minor_rate_alpha=float(strand_summary.minor_rate_alpha),
        minor_rate_beta=float(strand_summary.minor_rate_beta),
    )


def _density_normal_params(
    evidence: "DensityEvidence",
    region_idx: int,
    n_observed: int,
) -> tuple[float, float]:
    mean_arr = np.asarray(evidence.mean_unbounded, dtype=np.float64)
    mean = float(mean_arr[region_idx]) if region_idx < mean_arr.size else 0.0
    if evidence.variance_unbounded is None:
        variance = max(mean, 1.0)
    else:
        var_arr = np.asarray(evidence.variance_unbounded, dtype=np.float64)
        variance = float(var_arr[region_idx]) if region_idx < var_arr.size else max(mean, 1.0)
    if not np.isfinite(mean):
        mean = 0.0
    if not np.isfinite(variance) or variance < 0.0:
        variance = max(abs(mean), 1.0)
    return float(np.clip(mean, 0.0, n_observed)), variance


def _strand_normal_params(
    k_sense: int,
    n_strand: int,
    n_observed: int,
    *,
    kappa_d: float,
    p_r1_sense: float,
    strand_summary: "StrandSummary",
) -> tuple[float, float]:
    if abs(float(p_r1_sense) - 0.5) < STRAND_CONTRAST_NUMERICAL_FLOOR:
        return float(n_observed), float(max(n_observed, 1))

    d_grid = np.unique(np.rint(np.linspace(0, n_observed, 257)).astype(np.int64))
    k_minor = _minor_count_from_folded(k_sense, n_strand, p_r1_sense)
    log_likelihood = _approx_log_likelihood_minor_beta_binom(
        k_minor,
        n_strand,
        d_grid,
        kappa_d=kappa_d,
        minor_rate_alpha=float(strand_summary.minor_rate_alpha),
        minor_rate_beta=float(strand_summary.minor_rate_beta),
    )
    log_prior = _uniform_grid_log_weights(d_grid, n_observed)
    posterior = _posterior_weights_from_log(log_likelihood + log_prior)
    d_float = d_grid.astype(np.float64, copy=False)
    mean = float(np.sum(d_float * posterior, dtype=np.float64))
    variance = float(np.sum(((d_float - mean) ** 2) * posterior, dtype=np.float64))
    return float(np.clip(mean, 0.0, n_observed)), float(max(variance, _MIN_VARIANCE))


def _density_weight(density_tau: float, strand_tau: float) -> float:
    total = density_tau + strand_tau
    if total <= 0.0 or not np.isfinite(total):
        return 1.0
    return float(np.clip(density_tau / total, 0.0, 1.0))


def _strand_weight(density_tau: float, strand_tau: float) -> float:
    total = density_tau + strand_tau
    if total <= 0.0 or not np.isfinite(total):
        return 0.0
    return float(np.clip(strand_tau / total, 0.0, 1.0))


def _density_array_or_zeros(values: np.ndarray | None, size: int) -> np.ndarray:
    if values is None:
        return np.zeros(size, dtype=np.float64)
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != (size,):
        return np.zeros(size, dtype=np.float64)
    return np.nan_to_num(arr, nan=0.0, posinf=np.inf, neginf=0.0)


def _flag_count(flags: np.ndarray, bit: int) -> int:
    return int(np.sum((flags & np.uint8(bit)) != 0))


def _summary_stats(values: np.ndarray) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return {"min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0, "mean": 0.0}
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float64)
    return {
        "min": float(np.min(finite)),
        "p50": float(np.quantile(finite, 0.50)),
        "p95": float(np.quantile(finite, 0.95)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
    }
