"""Boundary-to-contained gDNA imputation for the v6 calibration path."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import nbinom

from .background_model import BackgroundModel
from .boundaries import BoundaryTable
from .density_observation import DensityObservation
from .strand_deconv import RegionGdnaChannelEstimate


__all__ = [
    "FLAG_BOUNDARY_LOW_OPPORTUNITY",
    "FLAG_BOUNDARY_NO_EVIDENCE",
    "FLAG_BOUNDARY_STRAND_RNA",
    "BoundaryLocalPosterior",
    "build_boundary_local_posterior",
    "predict_contained_gdna_from_excess",
    "region_boundary_indices",
]


FLAG_BOUNDARY_LOW_OPPORTUNITY: int = 1 << 0
FLAG_BOUNDARY_NO_EVIDENCE: int = 1 << 1
FLAG_BOUNDARY_STRAND_RNA: int = 1 << 2

_MIN_BETA: float = 1.0e-12


@dataclass(frozen=True, slots=True)
class BoundaryLocalPosterior:
    """Per-region boundary-derived Gamma excess and contained prediction."""

    alpha_excess: np.ndarray  # float64[R]
    beta_excess: np.ndarray  # float64[R]
    mu_local: np.ndarray  # float32[R]
    upper_local: np.ndarray  # float32[R]
    flags: np.ndarray  # uint16[R]


def _validate_confidence(confidence: float) -> float:
    out = float(confidence)
    if not np.isfinite(out) or not 0.5 <= out < 1.0:
        raise ValueError(f"confidence must be finite and in [0.5, 1.0); got {confidence!r}.")
    return out


def region_boundary_indices(boundaries: BoundaryTable) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(left_boundary, right_boundary)`` arrays in sorted region order."""
    ref_region_offsets = np.asarray(boundaries.ref_region_offsets, dtype=np.int64)
    ref_boundary_offsets = np.asarray(boundaries.ref_boundary_offsets, dtype=np.int64)
    n_refs = int(ref_region_offsets.shape[0] - 1)
    region_count = int(ref_region_offsets[-1])
    left_boundary = np.empty(region_count, dtype=np.int64)
    right_boundary = np.empty(region_count, dtype=np.int64)
    for ref_idx in range(n_refs):
        region_start = int(ref_region_offsets[ref_idx])
        region_end = int(ref_region_offsets[ref_idx + 1])
        boundary_start = int(ref_boundary_offsets[ref_idx])
        count = region_end - region_start
        if count <= 0:
            continue
        local = np.arange(count, dtype=np.int64)
        left_boundary[region_start:region_end] = boundary_start + local
        right_boundary[region_start:region_end] = boundary_start + local + 1
    return left_boundary, right_boundary


def _validate_region_shapes(
    observation: DensityObservation,
    boundaries: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> int:
    region_count = int(np.asarray(observation.contained_leff).shape[0])
    if int(boundaries.ref_region_offsets[-1]) != region_count:
        raise ValueError(
            "build_boundary_local_posterior: boundary/observation region counts differ; "
            f"got {int(boundaries.ref_region_offsets[-1])} and {region_count}."
        )
    if strand_channels is not None:
        for name in (
            "boundary_left_mean",
            "boundary_right_mean",
            "boundary_left_rna_lower",
            "boundary_right_rna_lower",
        ):
            arr = np.asarray(getattr(strand_channels, name))
            if arr.shape != (region_count,):
                raise ValueError(
                    f"build_boundary_local_posterior: strand_channels.{name} has shape "
                    f"{arr.shape}; expected {(region_count,)}."
                )
    return region_count


def predict_contained_gdna_from_excess(
    background: BackgroundModel,
    contained_leff: np.ndarray,
    alpha_excess: np.ndarray,
    beta_excess: np.ndarray,
    *,
    confidence: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict contained gDNA from background prior plus excess Gamma evidence."""
    conf = _validate_confidence(confidence)
    leff = np.asarray(contained_leff, dtype=np.float64)
    alpha_x = np.asarray(alpha_excess, dtype=np.float64)
    beta_x = np.asarray(beta_excess, dtype=np.float64)
    if leff.shape != alpha_x.shape or leff.shape != beta_x.shape:
        raise ValueError(
            "predict_contained_gdna_from_excess: leff/alpha/beta shapes must match; "
            f"got {leff.shape}, {alpha_x.shape}, and {beta_x.shape}."
        )

    alpha_post = float(background.rho_off_alpha) + np.maximum(alpha_x, 0.0)
    beta_post = float(background.rho_off_beta) + np.maximum(beta_x, 0.0)
    beta_post = np.maximum(beta_post, _MIN_BETA)
    leff_pos = np.maximum(leff, 0.0)

    mean = alpha_post * leff_pos / beta_post
    p_nb = beta_post / (beta_post + leff_pos)
    p_nb = np.clip(p_nb, np.finfo(np.float64).tiny, 1.0)
    upper = nbinom.ppf(conf, alpha_post, p_nb)
    upper = np.where(leff_pos <= 0.0, 0.0, upper)
    upper = np.maximum(np.nan_to_num(upper, nan=0.0, posinf=np.inf, neginf=0.0), mean)

    return mean.astype(np.float32), upper.astype(np.float32)


def _raw_boundary_counts(boundaries: BoundaryTable) -> tuple[np.ndarray, np.ndarray]:
    left_count = (
        np.asarray(boundaries.right_region_unspliced_pos, dtype=np.float64)
        + np.asarray(boundaries.right_region_unspliced_neg, dtype=np.float64)
    )
    right_count = (
        np.asarray(boundaries.left_region_unspliced_pos, dtype=np.float64)
        + np.asarray(boundaries.left_region_unspliced_neg, dtype=np.float64)
    )
    return left_count, right_count


def build_boundary_local_posterior(
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    confidence: float = 0.95,
    min_boundary_opportunity: float = 1.0,
) -> BoundaryLocalPosterior:
    """Build per-region boundary-to-contained predictions.

    ``alpha_excess`` and ``beta_excess`` are excess evidence only.  The weak
    background posterior is used when converting that evidence into local
    contained-count predictions, but it is not stored in the transferable arrays.
    """
    _validate_confidence(confidence)
    if not np.isfinite(min_boundary_opportunity) or float(min_boundary_opportunity) < 0.0:
        raise ValueError(
            "min_boundary_opportunity must be finite and >= 0; "
            f"got {min_boundary_opportunity!r}."
        )
    region_count = _validate_region_shapes(observation, boundaries, strand_channels)
    left_boundary, right_boundary = region_boundary_indices(boundaries)

    if strand_channels is None:
        boundary_left_count, boundary_right_count = _raw_boundary_counts(boundaries)
        left_alpha = boundary_left_count[left_boundary]
        right_alpha = boundary_right_count[right_boundary]
        strand_rna = np.zeros(region_count, dtype=bool)
    else:
        left_alpha = np.asarray(strand_channels.boundary_left_mean, dtype=np.float64)
        right_alpha = np.asarray(strand_channels.boundary_right_mean, dtype=np.float64)
        strand_rna = (
            (np.asarray(strand_channels.boundary_left_rna_lower) > 0.0)
            | (np.asarray(strand_channels.boundary_right_rna_lower) > 0.0)
        )

    left_beta = np.asarray(boundaries.right_region_boundary_leff, dtype=np.float64)[left_boundary]
    right_beta = np.asarray(boundaries.left_region_boundary_leff, dtype=np.float64)[right_boundary]
    alpha_excess = np.maximum(left_alpha, 0.0) + np.maximum(right_alpha, 0.0)
    beta_excess = np.maximum(left_beta, 0.0) + np.maximum(right_beta, 0.0)

    mu_local, upper_local = predict_contained_gdna_from_excess(
        background,
        observation.contained_leff,
        alpha_excess,
        beta_excess,
        confidence=confidence,
    )

    flags = np.zeros(region_count, dtype=np.uint16)
    flags[beta_excess <= float(min_boundary_opportunity)] |= FLAG_BOUNDARY_LOW_OPPORTUNITY
    flags[alpha_excess <= 0.0] |= FLAG_BOUNDARY_NO_EVIDENCE
    flags[strand_rna] |= FLAG_BOUNDARY_STRAND_RNA

    return BoundaryLocalPosterior(
        alpha_excess=alpha_excess.astype(np.float64, copy=False),
        beta_excess=beta_excess.astype(np.float64, copy=False),
        mu_local=mu_local,
        upper_local=upper_local,
        flags=flags,
    )