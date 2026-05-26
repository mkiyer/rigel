"""Sequential boundary evidence sweeps for the v6 calibration path."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .background_model import BackgroundModel
from .boundaries import BoundaryTable
from .boundary_model import (
    BoundaryLocalPosterior,
    predict_contained_gdna_from_excess,
)
from .density_observation import DensityObservation
from .strand_deconv import RegionGdnaChannelEstimate


__all__ = [
    "FLAG_SWEEP_FROM_LEFT",
    "FLAG_SWEEP_FROM_RIGHT",
    "BoundarySweepResult",
    "compute_boundary_transfer_weight",
    "run_boundary_sweep",
]


FLAG_SWEEP_FROM_LEFT: int = 1 << 0
FLAG_SWEEP_FROM_RIGHT: int = 1 << 1

_TRANSFER_MIN_OPPORTUNITY: float = 1.0
_TRANSFER_MASS_PSEUDOCOUNT: float = 1.0


@dataclass(frozen=True, slots=True)
class BoundarySweepResult:
    """Per-region boundary evidence after left/right sequential sweeps."""

    alpha_excess: np.ndarray  # float64[R]
    beta_excess: np.ndarray  # float64[R]
    forward_alpha_excess: np.ndarray  # float64[R]
    forward_beta_excess: np.ndarray
    reverse_alpha_excess: np.ndarray
    reverse_beta_excess: np.ndarray
    mu_sweep: np.ndarray  # float32[R]
    upper_sweep: np.ndarray  # float32[R]
    transfer_weight: np.ndarray  # float64[B]
    flags: np.ndarray  # uint16[R]


def _validate_region_count(
    boundaries: BoundaryTable,
    observation: DensityObservation,
    local: BoundaryLocalPosterior,
) -> int:
    region_count = int(np.asarray(observation.contained_leff).shape[0])
    if int(boundaries.ref_region_offsets[-1]) != region_count:
        raise ValueError(
            "run_boundary_sweep: boundary/observation region counts differ; "
            f"got {int(boundaries.ref_region_offsets[-1])} and {region_count}."
        )
    for name in ("alpha_excess", "beta_excess"):
        arr = np.asarray(getattr(local, name))
        if arr.shape != (region_count,):
            raise ValueError(
                f"run_boundary_sweep: local.{name} has shape {arr.shape}; "
                f"expected {(region_count,)}."
            )
    return region_count


def _boundary_raw_counts(boundaries: BoundaryTable) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    unspliced = (
        np.asarray(boundaries.left_region_unspliced_pos, dtype=np.float64)
        + np.asarray(boundaries.left_region_unspliced_neg, dtype=np.float64)
        + np.asarray(boundaries.right_region_unspliced_pos, dtype=np.float64)
        + np.asarray(boundaries.right_region_unspliced_neg, dtype=np.float64)
    )
    spliced = (
        np.asarray(boundaries.left_region_spliced_pos, dtype=np.float64)
        + np.asarray(boundaries.left_region_spliced_neg, dtype=np.float64)
        + np.asarray(boundaries.right_region_spliced_pos, dtype=np.float64)
        + np.asarray(boundaries.right_region_spliced_neg, dtype=np.float64)
    )
    opportunity = (
        np.asarray(boundaries.left_region_boundary_leff, dtype=np.float64)
        + np.asarray(boundaries.right_region_boundary_leff, dtype=np.float64)
    )
    return unspliced, spliced, opportunity


def _strand_boundary_counts(
    boundaries: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate,
) -> tuple[np.ndarray, np.ndarray]:
    boundary_count = boundaries.n_boundaries
    gdna = np.zeros(boundary_count, dtype=np.float64)
    rna = np.zeros(boundary_count, dtype=np.float64)
    left_region = boundaries.left_region_index()
    right_region = boundaries.right_region_index()

    has_left = left_region >= 0
    has_right = right_region >= 0
    gdna[has_left] += np.asarray(strand_channels.boundary_right_mean, dtype=np.float64)[
        left_region[has_left]
    ]
    rna[has_left] += np.asarray(strand_channels.boundary_right_rna_lower, dtype=np.float64)[
        left_region[has_left]
    ]
    gdna[has_right] += np.asarray(strand_channels.boundary_left_mean, dtype=np.float64)[
        right_region[has_right]
    ]
    rna[has_right] += np.asarray(strand_channels.boundary_left_rna_lower, dtype=np.float64)[
        right_region[has_right]
    ]
    return gdna, rna


def compute_boundary_transfer_weight(
    boundaries: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate | None,
    background: BackgroundModel,
) -> np.ndarray:
    """Return boundary evidence-reliability weights in ``[0, 1]``."""
    raw_unspliced, raw_spliced, opportunity = _boundary_raw_counts(boundaries)
    if strand_channels is None:
        evidence_mass = raw_unspliced
        rna_penalty_mass = raw_spliced
    else:
        evidence_mass, rna_penalty_mass = _strand_boundary_counts(boundaries, strand_channels)

    opportunity_factor = opportunity / (opportunity + _TRANSFER_MIN_OPPORTUNITY)
    expected_background = max(float(background.rho_off_mean), 0.0) * opportunity
    information_factor = evidence_mass / (
        evidence_mass + expected_background + _TRANSFER_MASS_PSEUDOCOUNT
    )
    rna_factor = 1.0 / (1.0 + np.maximum(rna_penalty_mass, 0.0))
    weight = opportunity_factor * information_factor * rna_factor
    weight = np.clip(np.nan_to_num(weight, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    weight[np.asarray(boundaries.is_terminal, dtype=bool)] = 0.0
    return weight.astype(np.float64, copy=False)


def _validate_transfer_weight(boundaries: BoundaryTable, weight: np.ndarray) -> np.ndarray:
    out = np.asarray(weight, dtype=np.float64)
    if out.shape != (boundaries.n_boundaries,):
        raise ValueError(
            "run_boundary_sweep: transfer_weight must have shape (n_boundaries,); "
            f"got {out.shape}, expected {(boundaries.n_boundaries,)}."
        )
    out = np.clip(np.nan_to_num(out, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    out[np.asarray(boundaries.is_terminal, dtype=bool)] = 0.0
    return out


def _forward_scan(
    local_alpha: np.ndarray,
    local_beta: np.ndarray,
    boundaries: BoundaryTable,
    weight: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    region_count = int(local_alpha.shape[0])
    from_left_alpha = np.zeros(region_count, dtype=np.float64)
    from_left_beta = np.zeros(region_count, dtype=np.float64)
    ref_region_offsets = np.asarray(boundaries.ref_region_offsets, dtype=np.int64)
    ref_boundary_offsets = np.asarray(boundaries.ref_boundary_offsets, dtype=np.int64)

    for ref_idx in range(int(ref_region_offsets.shape[0] - 1)):
        region_start = int(ref_region_offsets[ref_idx])
        region_end = int(ref_region_offsets[ref_idx + 1])
        boundary_start = int(ref_boundary_offsets[ref_idx])
        inclusive_alpha = 0.0
        inclusive_beta = 0.0
        for region_idx in range(region_start, region_end):
            local_i = region_idx - region_start
            if local_i > 0:
                w = float(weight[boundary_start + local_i])
                from_left_alpha[region_idx] = w * inclusive_alpha
                from_left_beta[region_idx] = w * inclusive_beta
            inclusive_alpha = float(local_alpha[region_idx]) + from_left_alpha[region_idx]
            inclusive_beta = float(local_beta[region_idx]) + from_left_beta[region_idx]

    return from_left_alpha, from_left_beta


def _reverse_scan(
    local_alpha: np.ndarray,
    local_beta: np.ndarray,
    boundaries: BoundaryTable,
    weight: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    region_count = int(local_alpha.shape[0])
    from_right_alpha = np.zeros(region_count, dtype=np.float64)
    from_right_beta = np.zeros(region_count, dtype=np.float64)
    ref_region_offsets = np.asarray(boundaries.ref_region_offsets, dtype=np.int64)
    ref_boundary_offsets = np.asarray(boundaries.ref_boundary_offsets, dtype=np.int64)

    for ref_idx in range(int(ref_region_offsets.shape[0] - 1)):
        region_start = int(ref_region_offsets[ref_idx])
        region_end = int(ref_region_offsets[ref_idx + 1])
        boundary_start = int(ref_boundary_offsets[ref_idx])
        inclusive_alpha = 0.0
        inclusive_beta = 0.0
        for region_idx in range(region_end - 1, region_start - 1, -1):
            local_i = region_idx - region_start
            if local_i < region_end - region_start - 1:
                w = float(weight[boundary_start + local_i + 1])
                from_right_alpha[region_idx] = w * inclusive_alpha
                from_right_beta[region_idx] = w * inclusive_beta
            inclusive_alpha = float(local_alpha[region_idx]) + from_right_alpha[region_idx]
            inclusive_beta = float(local_beta[region_idx]) + from_right_beta[region_idx]

    return from_right_alpha, from_right_beta


def run_boundary_sweep(
    local: BoundaryLocalPosterior,
    boundaries: BoundaryTable,
    observation: DensityObservation,
    background: BackgroundModel,
    *,
    transfer_weight: np.ndarray | None = None,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    confidence: float = 0.95,
) -> BoundarySweepResult:
    """Propagate boundary excess evidence left-to-right and right-to-left."""
    region_count = _validate_region_count(boundaries, observation, local)
    if transfer_weight is None:
        weight = compute_boundary_transfer_weight(boundaries, strand_channels, background)
    else:
        weight = _validate_transfer_weight(boundaries, transfer_weight)

    local_alpha = np.maximum(np.asarray(local.alpha_excess, dtype=np.float64), 0.0)
    local_beta = np.maximum(np.asarray(local.beta_excess, dtype=np.float64), 0.0)
    forward_alpha, forward_beta = _forward_scan(local_alpha, local_beta, boundaries, weight)
    reverse_alpha, reverse_beta = _reverse_scan(local_alpha, local_beta, boundaries, weight)
    alpha_total = local_alpha + forward_alpha + reverse_alpha
    beta_total = local_beta + forward_beta + reverse_beta
    mu_sweep, upper_sweep = predict_contained_gdna_from_excess(
        background,
        observation.contained_leff,
        alpha_total,
        beta_total,
        confidence=confidence,
    )

    flags = np.zeros(region_count, dtype=np.uint16)
    flags[forward_alpha > 0.0] |= FLAG_SWEEP_FROM_LEFT
    flags[reverse_alpha > 0.0] |= FLAG_SWEEP_FROM_RIGHT

    return BoundarySweepResult(
        alpha_excess=alpha_total.astype(np.float64, copy=False),
        beta_excess=beta_total.astype(np.float64, copy=False),
        forward_alpha_excess=forward_alpha,
        forward_beta_excess=forward_beta,
        reverse_alpha_excess=reverse_alpha,
        reverse_beta_excess=reverse_beta,
        mu_sweep=mu_sweep,
        upper_sweep=upper_sweep,
        transfer_weight=weight,
        flags=flags,
    )