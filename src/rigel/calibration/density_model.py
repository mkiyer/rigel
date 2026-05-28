"""gDNA density evidence model (v4 Phase 3).

This module turns a transient :class:`DensityObservation` into compact
per-region density evidence. It is intentionally independent of the
calibration orchestrator so Phase 3 can land before Phase 4 rewires the
pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

import math
from typing import Mapping

import numpy as np
from scipy.stats import nbinom

from .density_observation import DensityObservation


__all__ = [
    "DENSITY_FALLBACK_PRIOR_ALPHA",
    "DENSITY_MIN_BOUNDARY_INFO",
    "DENSITY_MIN_BOUNDARY_OPPORTUNITY",
    "DENSITY_MIN_EFF_LENGTH",
    "DENSITY_MIN_PRIOR_REGIONS",
    "DENSITY_PHI_FLOOR",
    "DENSITY_PHI_TRIM_UPPER",
    "DENSITY_PRIOR_MAX_PRECISION",
    "DENSITY_PRIOR_MIN_CV",
    "DENSITY_TAIL_PROBABILITY_WARN",
    "FLAG_FALLBACK_USED",
    "FLAG_HIGH_TAIL_TENSION",
    "FLAG_LOW_BOUNDARY_OPPORTUNITY",
    "FLAG_LOW_CONTAINED_OPPORTUNITY",
    "FLAG_NON_ANCHOR",
    "FLAG_PRIOR_DOMINATED",
    "PRIOR_FAMILY_ALL",
    "PRIOR_FAMILY_DETERMINISTIC_ZERO",
    "PRIOR_FAMILY_FALLBACK_BROAD",
    "PRIOR_FAMILY_INTERGENIC",
    "PRIOR_FAMILY_INTRON",
    "DensityEvidence",
    "GammaRatePrior",
    "InsufficientAnchors",
    "compute_beta_cap",
    "density_logpmf_grid",
    "fit_density_evidence",
    "fit_gamma_prior",
    "select_rho_ref",
]


# Prior precision floor: CV of the Gamma prior cannot fall below 5%.
DENSITY_PRIOR_MIN_CV: float = 0.05
DENSITY_PRIOR_MAX_PRECISION: float = 1.0 / (DENSITY_PRIOR_MIN_CV**2)  # 400.0

# Robust MoM controls.
DENSITY_PHI_FLOOR: float = 1.0e-6
DENSITY_PHI_TRIM_UPPER: float = 0.05
DENSITY_FALLBACK_PRIOR_ALPHA: float = 1.0

# Opportunity / information thresholds.
DENSITY_MIN_EFF_LENGTH: float = 1.0
DENSITY_INFO_VAR_FLOOR: float = 1.0e-12
DENSITY_MIN_BOUNDARY_OPPORTUNITY: float = 1.0
DENSITY_MIN_BOUNDARY_INFO: float = 0.05
DENSITY_MIN_PRIOR_REGIONS: int = 20

# Predictive tension threshold.
DENSITY_TAIL_PROBABILITY_WARN: float = 0.95

_SCORE_FLOOR: float = 1.0e-6
_PHI_EPS: float = 1.0e-12
_RHO_EPS: float = 1.0e-12
_STAT_EPS: float = 1.0e-12


# uint8 prior family codes.
PRIOR_FAMILY_INTERGENIC: int = 0
PRIOR_FAMILY_INTRON: int = 1
PRIOR_FAMILY_ALL: int = 2
PRIOR_FAMILY_FALLBACK_BROAD: int = 3
PRIOR_FAMILY_DETERMINISTIC_ZERO: int = 4

_PRIOR_FAMILY_LABELS = {
    PRIOR_FAMILY_INTERGENIC: "INTERGENIC",
    PRIOR_FAMILY_INTRON: "INTRON",
    PRIOR_FAMILY_ALL: "ALL",
    PRIOR_FAMILY_FALLBACK_BROAD: "FALLBACK_BROAD",
    PRIOR_FAMILY_DETERMINISTIC_ZERO: "DETERMINISTIC_ZERO",
}


# uint8 flag bits.
FLAG_LOW_CONTAINED_OPPORTUNITY: int = 0x01
FLAG_LOW_BOUNDARY_OPPORTUNITY: int = 0x02
FLAG_PRIOR_DOMINATED: int = 0x04
FLAG_HIGH_TAIL_TENSION: int = 0x08
FLAG_NON_ANCHOR: int = 0x10
FLAG_FALLBACK_USED: int = 0x20


class InsufficientAnchors(ValueError):
    """Raised when an anchor family cannot fit a positive Gamma-rate prior."""


@dataclass(frozen=True, slots=True)
class GammaRatePrior:
    """Gamma prior over a gDNA fragment density rate."""

    family: str
    alpha: float
    beta: float
    mean_density: float
    phi: float
    beta_raw: float
    beta_cap: float
    cap_applied: bool
    n_regions: int
    n_fragments: float
    eff_length: float
    residual_sum: float
    poisson_variance_sum: float
    extra_variance_basis_sum: float
    trim_upper: float
    n_trimmed: int
    trimmed_mu_fraction: float
    fallback_depth: int
    fit_status: str
    pearson_chi2_trimmed: float

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "family": self.family,
            "alpha": float(self.alpha),
            "beta": float(self.beta),
            "mean_density": float(self.mean_density),
            "phi": float(self.phi),
            "beta_raw": float(self.beta_raw),
            "beta_cap": float(self.beta_cap),
            "cap_applied": bool(self.cap_applied),
            "n_regions": int(self.n_regions),
            "n_fragments": float(self.n_fragments),
            "eff_length": float(self.eff_length),
            "residual_sum": float(self.residual_sum),
            "poisson_variance_sum": float(self.poisson_variance_sum),
            "extra_variance_basis_sum": float(self.extra_variance_basis_sum),
            "trim_upper": float(self.trim_upper),
            "n_trimmed": int(self.n_trimmed),
            "trimmed_mu_fraction": float(self.trimmed_mu_fraction),
            "fallback_depth": int(self.fallback_depth),
            "fit_status": self.fit_status,
            "pearson_chi2_trimmed": float(self.pearson_chi2_trimmed),
        }


@dataclass(frozen=True, slots=True)
class DensityEvidence:
    """Compact per-region gDNA density evidence.

    v4 intentionally stores only the arrays downstream consumers need.
    Full predictive diagnostics are computed transiently during fitting and
    reduced to flags or summary statistics.
    """

    rho_post: np.ndarray
    relative_exposure: np.ndarray
    mean_unbounded: np.ndarray
    upper_unbounded: np.ndarray
    prior_family: np.ndarray
    fallback_depth: np.ndarray
    flags: np.ndarray
    confidence: float
    priors: dict[str, GammaRatePrior]
    rho_ref: float
    rho_ref_source: str
    alpha_post: np.ndarray | None = None
    beta_post: np.ndarray | None = None
    contained_leff: np.ndarray | None = None
    boundary_count: np.ndarray | None = None
    variance_unbounded: np.ndarray | None = None
    tail_probability: np.ndarray | None = None
    expected_tail_count: np.ndarray | None = None
    information: np.ndarray | None = None
    applicable: np.ndarray | None = None

    def to_summary_dict(self) -> dict[str, object]:
        """Return a compact JSON-safe summary of stored evidence arrays."""
        flags = np.asarray(self.flags, dtype=np.uint8)
        return {
            "confidence": float(self.confidence),
            "n_regions": int(flags.size),
            "rho_ref": float(self.rho_ref),
            "rho_ref_source": self.rho_ref_source,
            "priors": {
                name: self.priors[name].to_summary_dict() if name in self.priors else None
                for name in ("INTERGENIC", "INTRON", "ALL")
            },
            "fallback_depth_histogram": _uint_histogram(self.fallback_depth, range(4)),
            "prior_family_histogram": {
                label: int(np.sum(self.prior_family == code))
                for code, label in _PRIOR_FAMILY_LABELS.items()
            },
            "flags": {
                "n_low_contained_opportunity": _flag_count(flags, FLAG_LOW_CONTAINED_OPPORTUNITY),
                "n_low_boundary_opportunity": _flag_count(flags, FLAG_LOW_BOUNDARY_OPPORTUNITY),
                "n_prior_dominated": _flag_count(flags, FLAG_PRIOR_DOMINATED),
                "n_high_tail_tension": _flag_count(flags, FLAG_HIGH_TAIL_TENSION),
                "n_fallback_used": _flag_count(flags, FLAG_FALLBACK_USED),
                "n_non_anchor": _flag_count(flags, FLAG_NON_ANCHOR),
            },
            "relative_exposure": _summary_stats(self.relative_exposure),
            "rho_post": _summary_stats(self.rho_post, include_min=False, include_mean=False),
            "mean_unbounded": _summary_stats(
                self.mean_unbounded, include_min=False, include_mean=False
            ),
            "upper_unbounded": _summary_stats(
                self.upper_unbounded, include_min=False, include_mean=False
            ),
            "tail_probability": _summary_stats(
                self.tail_probability if self.tail_probability is not None else np.zeros(0),
                include_min=False,
            ),
            "expected_tail_count": _summary_stats(
                self.expected_tail_count if self.expected_tail_count is not None else np.zeros(0),
                include_min=False,
            ),
        }


@dataclass(frozen=True, slots=True)
class _PredictiveDiagnostics:
    alpha_post: np.ndarray
    beta_post: np.ndarray
    rho_post: np.ndarray
    mean_c: np.ndarray
    p_nb: np.ndarray
    mean_unbounded: np.ndarray
    upper_unbounded: np.ndarray
    tail_probability: np.ndarray
    expected_tail_count: np.ndarray
    density_over_observed_ratio: np.ndarray
    w_boundary_opportunity: np.ndarray
    w_boundary_count: np.ndarray
    variance_unbounded: np.ndarray
    precision_unbounded: np.ndarray


def compute_beta_cap(rho: float) -> float:
    """Cap fitted prior precision so prior CV cannot fall below 5%."""
    rho_f = float(rho)
    if rho_f <= 0.0 or not np.isfinite(rho_f):
        return float("inf")
    return float(DENSITY_PRIOR_MAX_PRECISION / rho_f)


def density_logpmf_grid(
    evidence: DensityEvidence,
    region_idx: int,
    d_grid: np.ndarray,
) -> np.ndarray:
    """Return ``log P_density(D=d)`` for one region and integer ``d_grid``.

    The density model is a Gamma-rate posterior predictive distribution over
    contained gDNA counts, shifted by the observed boundary-compatible count.
    Boundary counts can be fractional because calibration accumulates fractional
    evidence; for the integer fusion grid we evaluate the nearest contained-count
    support point after subtracting the boundary shift.
    """
    idx = int(region_idx)
    grid = np.asarray(d_grid, dtype=np.float64)

    if idx < 0 or idx >= np.asarray(evidence.mean_unbounded).size:
        raise IndexError(
            f"density_logpmf_grid: region_idx {idx} out of bounds for "
            f"{np.asarray(evidence.mean_unbounded).size} regions."
        )

    alpha_arr = evidence.alpha_post
    beta_arr = evidence.beta_post
    leff_arr = evidence.contained_leff
    boundary_arr = evidence.boundary_count
    if alpha_arr is None or beta_arr is None or leff_arr is None or boundary_arr is None:
        center = int(round(float(np.asarray(evidence.mean_unbounded, dtype=np.float64)[idx])))
        return np.where(np.rint(grid).astype(np.int64) == center, 0.0, -np.inf)

    alpha = float(np.asarray(alpha_arr, dtype=np.float64)[idx])
    beta = float(np.asarray(beta_arr, dtype=np.float64)[idx])
    leff = float(np.asarray(leff_arr, dtype=np.float64)[idx])
    boundary = float(np.asarray(boundary_arr, dtype=np.float64)[idx])

    if leff <= 0.0 or alpha <= 0.0 or beta <= 0.0 or not np.isfinite(alpha + beta + leff):
        center = 0 if alpha <= 0.0 else int(round(boundary))
        return np.where(np.rint(grid).astype(np.int64) == center, 0.0, -np.inf)

    contained_grid = np.rint(grid - boundary).astype(np.int64)
    out = np.full(grid.shape, -np.inf, dtype=np.float64)
    valid = contained_grid >= 0
    if not np.any(valid):
        return out

    p_nb = beta / (beta + leff)
    p_nb = float(np.clip(p_nb, np.finfo(np.float64).tiny, 1.0))
    out[valid] = nbinom.logpmf(contained_grid[valid], alpha, p_nb)
    return np.nan_to_num(out, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)


def fit_gamma_prior(
    counts: np.ndarray,
    opportunities: np.ndarray,
    *,
    family: str,
    beta_cap: float | None = None,
    phi_floor: float = DENSITY_PHI_FLOOR,
    trim_upper: float = DENSITY_PHI_TRIM_UPPER,
    min_regions: int = DENSITY_MIN_PRIOR_REGIONS,
    min_eff_length: float = DENSITY_MIN_EFF_LENGTH,
    fallback_depth: int = 0,
) -> GammaRatePrior:
    """Fit a robust Pearson-trimmed Gamma-rate prior for one anchor family."""
    C_all = np.asarray(counts, dtype=np.float64)
    L_all = np.asarray(opportunities, dtype=np.float64)
    if C_all.shape != L_all.shape:
        raise ValueError(
            f"fit_gamma_prior: counts/opportunities shape mismatch: {C_all.shape} vs {L_all.shape}."
        )
    eligible = (L_all >= float(min_eff_length)) & np.isfinite(C_all) & np.isfinite(L_all)
    n_eligible = int(np.sum(eligible))
    if n_eligible < int(min_regions):
        raise InsufficientAnchors(
            f"fit_gamma_prior({family}): {n_eligible} eligible anchors < min_regions={min_regions}."
        )

    C = C_all[eligible]
    L = L_all[eligible]
    n_fragments = float(np.sum(C))
    eff_length = float(np.sum(L))
    if eff_length <= 0.0 or n_fragments <= 0.0:
        raise InsufficientAnchors(
            f"fit_gamma_prior({family}): non-positive density evidence "
            f"(n_fragments={n_fragments}, eff_length={eff_length})."
        )

    rho = float(n_fragments / eff_length)
    mu = rho * L
    raw_residual = (C - mu) ** 2
    score = raw_residual / np.maximum(mu, _SCORE_FLOOR)

    trim = min(max(float(trim_upper), 0.0), 1.0)
    k_keep = max(int(min_regions), int(math.ceil((1.0 - trim) * score.size)))
    k_keep = min(k_keep, score.size)
    if k_keep == score.size:
        keep = np.arange(score.size)
    else:
        keep = np.argpartition(score, k_keep - 1)[:k_keep]

    S_trim = float(np.sum(raw_residual[keep]))
    B_trim = float(np.sum(mu[keep]))
    A_trim = float(np.sum(mu[keep] ** 2))
    pearson_chi2_trimmed = float(np.sum(score[keep]))

    phi_hat = max(0.0, (S_trim - B_trim) / max(A_trim, _PHI_EPS))
    phi = float(max(phi_hat, float(phi_floor)))
    alpha_raw = float(1.0 / phi)
    beta_raw = float(alpha_raw / max(rho, _RHO_EPS))

    cap = compute_beta_cap(rho) if beta_cap is None else float(beta_cap)
    beta = float(min(beta_raw, cap))
    alpha = float(rho * beta)
    cap_applied = bool(np.isfinite(cap) and beta_raw > cap)

    kept = np.zeros(score.size, dtype=bool)
    kept[keep] = True
    mu_total = float(np.sum(mu))
    trimmed_mu_fraction = (
        float(np.sum(mu[~kept]) / mu_total) if mu_total > 0.0 and np.any(~kept) else 0.0
    )

    return GammaRatePrior(
        family=str(family),
        alpha=alpha,
        beta=beta,
        mean_density=float(alpha / beta) if beta > 0.0 else 0.0,
        phi=float(1.0 / alpha) if alpha > 0.0 else float("inf"),
        beta_raw=beta_raw,
        beta_cap=cap,
        cap_applied=cap_applied,
        n_regions=n_eligible,
        n_fragments=n_fragments,
        eff_length=eff_length,
        residual_sum=S_trim,
        poisson_variance_sum=B_trim,
        extra_variance_basis_sum=A_trim,
        trim_upper=trim,
        n_trimmed=int(score.size - k_keep),
        trimmed_mu_fraction=trimmed_mu_fraction,
        fallback_depth=int(fallback_depth),
        fit_status="ok",
        pearson_chi2_trimmed=pearson_chi2_trimmed,
    )


def select_rho_ref(priors: Mapping[str, GammaRatePrior]) -> tuple[float, str]:
    """Return the v4 reference density and source label."""
    all_prior = priors.get("ALL")
    if _prior_ok(all_prior):
        return float(all_prior.mean_density), "ALL"

    weighted_num = 0.0
    weighted_den = 0.0
    for name in ("INTERGENIC", "INTRON"):
        prior = priors.get(name)
        if not _prior_ok(prior):
            continue
        weighted_num += float(prior.alpha) * float(prior.mean_density)
        weighted_den += float(prior.alpha)
    if weighted_den > 0.0:
        return float(weighted_num / weighted_den), "WEIGHTED_FAMILIES"
    return 0.0, "ZERO"


def fit_density_evidence(
    observation: DensityObservation,
    *,
    confidence: float = 0.95,
    min_eff_length: float = DENSITY_MIN_EFF_LENGTH,
) -> DensityEvidence:
    """Fit compact density evidence for all regions in a `DensityObservation`."""
    conf = float(confidence)
    if not 0.5 <= conf < 1.0:
        raise ValueError(f"fit_density_evidence: confidence must be in [0.5, 1.0); got {conf}.")

    counts = np.asarray(observation.contained_count, dtype=np.float64)
    contained_leff = np.asarray(observation.contained_leff, dtype=np.float64)
    boundary_count = np.asarray(observation.boundary_count, dtype=np.float64)
    boundary_leff = np.asarray(observation.boundary_leff, dtype=np.float64)
    R = counts.size

    priors = _fit_anchor_priors(observation, min_eff_length=float(min_eff_length))
    rho_ref, rho_ref_source = select_rho_ref(priors)
    if rho_ref <= 0.0 or not np.isfinite(rho_ref):
        sparse_intergenic = _sparse_intergenic_prior(
            observation,
            min_eff_length=float(min_eff_length),
        )
        if sparse_intergenic is None:
            return _deterministic_zero_evidence(observation, confidence=conf)
        priors = {"INTERGENIC": sparse_intergenic}
        rho_ref, rho_ref_source = select_rho_ref(priors)

    alpha_prior, beta_prior, prior_family, fallback_depth = _assign_region_priors(
        observation,
        priors,
        rho_ref=rho_ref,
    )

    diag = _compute_predictive_diagnostics(
        contained_count=counts,
        boundary_count=boundary_count,
        contained_leff=contained_leff,
        boundary_leff=boundary_leff,
        alpha_prior=alpha_prior,
        beta_prior=beta_prior,
        confidence=conf,
    )
    relative_exposure = np.divide(
        diag.rho_post,
        rho_ref,
        out=np.ones(R, dtype=np.float64),
        where=rho_ref > 0.0,
    )

    flags = np.zeros(R, dtype=np.uint8)
    flags[contained_leff < float(min_eff_length)] |= FLAG_LOW_CONTAINED_OPPORTUNITY
    flags[boundary_leff < DENSITY_MIN_BOUNDARY_OPPORTUNITY] |= FLAG_LOW_BOUNDARY_OPPORTUNITY
    flags[diag.w_boundary_opportunity < DENSITY_MIN_BOUNDARY_INFO] |= FLAG_PRIOR_DOMINATED
    flags[diag.tail_probability > DENSITY_TAIL_PROBABILITY_WARN] |= FLAG_HIGH_TAIL_TENSION
    flags[~np.asarray(observation.is_anchor, dtype=bool)] |= FLAG_NON_ANCHOR
    flags[fallback_depth > 0] |= FLAG_FALLBACK_USED

    return DensityEvidence(
        rho_post=np.asarray(diag.rho_post, dtype=np.float64),
        relative_exposure=np.asarray(relative_exposure, dtype=np.float64),
        mean_unbounded=np.asarray(diag.mean_unbounded, dtype=np.float64),
        upper_unbounded=np.asarray(diag.upper_unbounded, dtype=np.float64),
        prior_family=np.asarray(prior_family, dtype=np.uint8),
        fallback_depth=np.asarray(fallback_depth, dtype=np.uint8),
        flags=flags,
        confidence=conf,
        priors=dict(priors),
        rho_ref=float(rho_ref),
        rho_ref_source=rho_ref_source,
        alpha_post=np.asarray(diag.alpha_post, dtype=np.float64),
        beta_post=np.asarray(diag.beta_post, dtype=np.float64),
        contained_leff=np.asarray(contained_leff, dtype=np.float64),
        boundary_count=np.asarray(boundary_count, dtype=np.float64),
        variance_unbounded=np.asarray(diag.variance_unbounded, dtype=np.float64),
        tail_probability=np.asarray(diag.tail_probability, dtype=np.float64),
        expected_tail_count=np.asarray(diag.expected_tail_count, dtype=np.float64),
        information=_density_information(diag.variance_unbounded),
        applicable=np.asarray(contained_leff >= float(min_eff_length), dtype=bool),
    )


def _density_information(variance_unbounded: np.ndarray) -> np.ndarray:
    """Per-region density precision: 1 / max(variance, var_floor).

    Zero variance (deterministic prior or empty opportunity) yields zero
    information, signalling to the fusion engine that this region's density
    channel carries no usable evidence.
    """
    var = np.asarray(variance_unbounded, dtype=np.float64)
    info = np.zeros_like(var)
    np.divide(1.0, var, out=info, where=var > DENSITY_INFO_VAR_FLOOR)
    return info


def _fit_anchor_priors(
    observation: DensityObservation,
    *,
    min_eff_length: float,
) -> dict[str, GammaRatePrior]:
    counts = np.asarray(observation.contained_count, dtype=np.float64)
    leff = np.asarray(observation.contained_leff, dtype=np.float64)
    masks = {
        "INTERGENIC": np.asarray(observation.anchor_intergenic, dtype=bool),
        "INTRON": np.asarray(observation.anchor_intron, dtype=bool),
        "ALL": np.asarray(observation.is_anchor, dtype=bool),
    }
    priors: dict[str, GammaRatePrior] = {}
    for family, mask in masks.items():
        try:
            priors[family] = fit_gamma_prior(
                counts[mask],
                leff[mask],
                family=family,
                min_eff_length=min_eff_length,
                fallback_depth=0 if family in {"INTERGENIC", "INTRON"} else 1,
            )
        except InsufficientAnchors:
            continue
    return priors


def _sparse_intergenic_prior(
    observation: DensityObservation,
    *,
    min_eff_length: float,
) -> GammaRatePrior | None:
    counts = np.asarray(observation.contained_count, dtype=np.float64)
    leff = np.asarray(observation.contained_leff, dtype=np.float64)
    mask = (
        np.asarray(observation.anchor_intergenic, dtype=bool)
        & np.isfinite(counts)
        & np.isfinite(leff)
        & (leff >= float(min_eff_length))
    )
    if not np.any(mask):
        return None
    n_fragments = float(np.sum(counts[mask]))
    eff_length = float(np.sum(leff[mask]))
    if n_fragments <= 0.0 or eff_length <= 0.0:
        return None

    rho = float(n_fragments / eff_length)
    alpha = float(DENSITY_FALLBACK_PRIOR_ALPHA)
    beta = float(alpha / max(rho, _RHO_EPS))
    return GammaRatePrior(
        family="INTERGENIC_SPARSE",
        alpha=alpha,
        beta=beta,
        mean_density=rho,
        phi=float(1.0 / alpha),
        beta_raw=beta,
        beta_cap=float("inf"),
        cap_applied=False,
        n_regions=int(np.sum(mask)),
        n_fragments=n_fragments,
        eff_length=eff_length,
        residual_sum=0.0,
        poisson_variance_sum=n_fragments,
        extra_variance_basis_sum=0.0,
        trim_upper=0.0,
        n_trimmed=0,
        trimmed_mu_fraction=0.0,
        fallback_depth=1,
        fit_status="sparse_intergenic_fallback",
        pearson_chi2_trimmed=0.0,
    )


def _assign_region_priors(
    observation: DensityObservation,
    priors: Mapping[str, GammaRatePrior],
    *,
    rho_ref: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    R = int(np.asarray(observation.contained_count).size)
    alpha = np.zeros(R, dtype=np.float64)
    beta = np.zeros(R, dtype=np.float64)
    family = np.full(R, PRIOR_FAMILY_FALLBACK_BROAD, dtype=np.uint8)
    depth = np.full(R, 2, dtype=np.uint8)
    unassigned = np.ones(R, dtype=bool)

    def assign(mask: np.ndarray, prior: GammaRatePrior, depth_value: int, family_code: int) -> None:
        selected = unassigned & np.asarray(mask, dtype=bool)
        if not np.any(selected):
            return
        alpha[selected] = float(prior.alpha)
        beta[selected] = float(prior.beta)
        family[selected] = np.uint8(family_code)
        depth[selected] = np.uint8(depth_value)
        unassigned[selected] = False

    intergenic = np.asarray(observation.anchor_intergenic, dtype=bool)
    intron = np.asarray(observation.anchor_intron, dtype=bool)
    intergenic_prior = priors.get("INTERGENIC")
    intron_prior = priors.get("INTRON")
    all_prior = priors.get("ALL")

    if _prior_ok(intergenic_prior):
        assign(intergenic, intergenic_prior, 0, PRIOR_FAMILY_INTERGENIC)
    if _prior_ok(intron_prior):
        assign(intron, intron_prior, 0, PRIOR_FAMILY_INTRON)
    if _prior_ok(all_prior):
        assign(unassigned, all_prior, 1, PRIOR_FAMILY_ALL)
    if np.any(unassigned):
        fallback = _broad_fallback_prior(float(rho_ref))
        assign(unassigned, fallback, 2, PRIOR_FAMILY_FALLBACK_BROAD)

    return alpha, beta, family, depth


def _broad_fallback_prior(rho_ref: float) -> GammaRatePrior:
    if rho_ref <= 0.0 or not np.isfinite(rho_ref):
        raise ValueError(f"_broad_fallback_prior: rho_ref must be positive; got {rho_ref!r}")
    alpha = float(DENSITY_FALLBACK_PRIOR_ALPHA)
    beta = float(alpha / rho_ref)
    return GammaRatePrior(
        family="FALLBACK_BROAD",
        alpha=alpha,
        beta=beta,
        mean_density=float(alpha / beta),
        phi=float(1.0 / alpha),
        beta_raw=beta,
        beta_cap=float("inf"),
        cap_applied=False,
        n_regions=0,
        n_fragments=0.0,
        eff_length=0.0,
        residual_sum=0.0,
        poisson_variance_sum=0.0,
        extra_variance_basis_sum=0.0,
        trim_upper=0.0,
        n_trimmed=0,
        trimmed_mu_fraction=0.0,
        fallback_depth=2,
        fit_status="fallback_broad",
        pearson_chi2_trimmed=0.0,
    )


def _deterministic_zero_evidence(
    observation: DensityObservation,
    *,
    confidence: float,
) -> DensityEvidence:
    counts = np.asarray(observation.contained_count, dtype=np.float64)
    boundary = np.asarray(observation.boundary_count, dtype=np.float64)
    R = int(counts.size)
    flags = np.full(R, FLAG_FALLBACK_USED | FLAG_PRIOR_DOMINATED, dtype=np.uint8)
    return DensityEvidence(
        rho_post=np.zeros(R, dtype=np.float64),
        relative_exposure=np.ones(R, dtype=np.float64),
        mean_unbounded=np.zeros(R, dtype=np.float64),
        upper_unbounded=np.zeros(R, dtype=np.float64),
        prior_family=np.full(R, PRIOR_FAMILY_DETERMINISTIC_ZERO, dtype=np.uint8),
        fallback_depth=np.full(R, 3, dtype=np.uint8),
        flags=flags,
        confidence=float(confidence),
        priors={},
        rho_ref=0.0,
        rho_ref_source="ZERO",
        alpha_post=np.zeros(R, dtype=np.float64),
        beta_post=np.ones(R, dtype=np.float64),
        contained_leff=np.asarray(observation.contained_leff, dtype=np.float64),
        boundary_count=np.asarray(boundary, dtype=np.float64),
        variance_unbounded=np.zeros(R, dtype=np.float64),
        tail_probability=np.zeros(R, dtype=np.float64),
        expected_tail_count=np.zeros(R, dtype=np.float64),
        information=np.zeros(R, dtype=np.float64),
        applicable=np.zeros(R, dtype=bool),
    )


def _compute_predictive_diagnostics(
    *,
    contained_count: np.ndarray,
    boundary_count: np.ndarray,
    contained_leff: np.ndarray,
    boundary_leff: np.ndarray,
    alpha_prior: np.ndarray,
    beta_prior: np.ndarray,
    confidence: float,
) -> _PredictiveDiagnostics:
    C = np.asarray(contained_count, dtype=np.float64)
    B = np.asarray(boundary_count, dtype=np.float64)
    L_c = np.asarray(contained_leff, dtype=np.float64)
    L_b = np.asarray(boundary_leff, dtype=np.float64)
    alpha0 = np.asarray(alpha_prior, dtype=np.float64)
    beta0 = np.asarray(beta_prior, dtype=np.float64)

    alpha_post = alpha0 + B
    beta_post = beta0 + L_b
    rho_post = np.divide(
        alpha_post,
        beta_post,
        out=np.zeros_like(alpha_post, dtype=np.float64),
        where=beta_post > 0.0,
    )
    mean_c = rho_post * L_c
    p_nb = np.divide(
        beta_post,
        beta_post + L_c,
        out=np.ones_like(beta_post, dtype=np.float64),
        where=(beta_post + L_c) > 0.0,
    )
    p_nb = np.clip(p_nb, np.finfo(np.float64).tiny, 1.0)

    upper_c = nbinom.ppf(float(confidence), alpha_post, p_nb)
    upper_c = np.where(L_c <= 0.0, 0.0, upper_c)
    upper_c = np.nan_to_num(upper_c, nan=0.0, posinf=np.inf, neginf=0.0)
    mean_unbounded = B + mean_c
    upper_unbounded = B + upper_c

    c_int = np.floor(np.maximum(C, 0.0)).astype(np.int64)
    sf_c = nbinom.sf(c_int, alpha_post, p_nb)
    sf_y = nbinom.sf(c_int - 1, alpha_post + 1.0, p_nb)
    tail_probability = np.nan_to_num(sf_c, nan=0.0, posinf=1.0, neginf=0.0)
    expected_tail_count = mean_c * np.nan_to_num(sf_y, nan=0.0) - c_int.astype(np.float64) * sf_c
    expected_tail_count = np.maximum(
        np.nan_to_num(expected_tail_count, nan=0.0, posinf=np.inf, neginf=0.0),
        0.0,
    )

    observed_total = B + C
    density_over_observed_ratio = np.divide(
        mean_unbounded,
        np.maximum(observed_total, _STAT_EPS),
        out=np.zeros_like(mean_unbounded),
    )
    variance = mean_c * (
        1.0 + np.divide(L_c, beta_post, out=np.zeros_like(L_c), where=beta_post > 0.0)
    )
    cv = np.divide(
        np.sqrt(np.maximum(variance, 0.0)),
        np.maximum(mean_c, _STAT_EPS),
        out=np.full_like(mean_c, np.inf),
    )
    precision = 1.0 / (1.0 + cv)
    precision = np.where(np.isfinite(precision), precision, 0.0)
    w_boundary_opportunity = np.divide(
        L_b,
        beta0 + L_b,
        out=np.zeros_like(L_b),
        where=(beta0 + L_b) > 0.0,
    )
    w_boundary_count = np.divide(
        B,
        alpha0 + B,
        out=np.zeros_like(B),
        where=B > 0.0,
    )

    return _PredictiveDiagnostics(
        alpha_post=np.asarray(alpha_post, dtype=np.float64),
        beta_post=np.asarray(beta_post, dtype=np.float64),
        rho_post=np.asarray(rho_post, dtype=np.float64),
        mean_c=np.asarray(mean_c, dtype=np.float64),
        p_nb=np.asarray(p_nb, dtype=np.float64),
        mean_unbounded=np.asarray(mean_unbounded, dtype=np.float64),
        upper_unbounded=np.asarray(upper_unbounded, dtype=np.float64),
        tail_probability=np.asarray(tail_probability, dtype=np.float64),
        expected_tail_count=np.asarray(expected_tail_count, dtype=np.float64),
        density_over_observed_ratio=np.asarray(density_over_observed_ratio, dtype=np.float64),
        w_boundary_opportunity=np.asarray(w_boundary_opportunity, dtype=np.float64),
        w_boundary_count=np.asarray(w_boundary_count, dtype=np.float64),
        variance_unbounded=np.asarray(variance, dtype=np.float64),
        precision_unbounded=np.asarray(precision, dtype=np.float64),
    )


def _prior_ok(prior: GammaRatePrior | None) -> bool:
    return (
        prior is not None
        and prior.fit_status in {"ok", "sparse_intergenic_fallback"}
        and prior.mean_density > 0.0
        and prior.alpha > 0.0
        and prior.beta > 0.0
        and np.isfinite(prior.mean_density)
        and np.isfinite(prior.alpha)
        and np.isfinite(prior.beta)
    )


def _flag_count(flags: np.ndarray, bit: int) -> int:
    return int(np.sum((flags & np.uint8(bit)) != 0))


def _uint_histogram(values: np.ndarray, bins) -> dict[int, int]:
    arr = np.asarray(values, dtype=np.uint8)
    return {int(v): int(np.sum(arr == int(v))) for v in bins}


def _summary_stats(
    values: np.ndarray,
    *,
    include_min: bool = True,
    include_mean: bool = True,
) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        out: dict[str, float] = {"p50": 0.0, "p95": 0.0, "max": 0.0}
        if include_min:
            out = {"min": 0.0, **out}
        if include_mean:
            out["mean"] = 0.0
        return out
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float64)
    out = {
        "p50": float(np.quantile(finite, 0.50)),
        "p95": float(np.quantile(finite, 0.95)),
        "max": float(np.max(finite)),
    }
    if include_min:
        out = {"min": float(np.min(finite)), **out}
    if include_mean:
        out["mean"] = float(np.mean(finite))
    return out
