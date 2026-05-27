"""Immutable calibration result for the v6 RegionCalibration cutover."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..frag_length_model import FragmentLengthModels
from ._diagnostics import Diagnostics
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from .fl import POOL_EB_PRIOR_ESS, FLModels, build_fl_models
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:  # pragma: no cover - imported for annotations only.
    from .background_model import BackgroundModel
    from .boundary_model import BoundaryLocalPosterior
    from .boundary_sweep import BoundarySweepResult
    from .calibration_iteration import RegionCalibration
    from .prior import PriorTable
    from .strand_deconv import RegionGdnaChannelEstimate


__all__ = ["CalibrationResult", "build_calibration_result"]


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Calibration handoff produced before locus-level EM.

    The production calibration contract is now :class:`RegionCalibration`.
    Legacy density/fusion/exposure fields are intentionally absent; downstream
    prior assembly consumes ``region_calibration.mu_gdna``, ``upper_gdna``, and
    ``A_r`` directly.
    """

    fl_models: FLModels
    diagnostics: Diagnostics
    region_calibration: "RegionCalibration"
    strand_channels: "RegionGdnaChannelEstimate | None"
    background_model: "BackgroundModel"
    boundary_local: "BoundaryLocalPosterior"
    boundary_sweep: "BoundarySweepResult"
    prior_table: "PriorTable | None" = None
    n_multi_loci: int = 0

    @property
    def global_fl_mean(self) -> float:
        return float(self.fl_models.global_.mean)

    @property
    def rna_fl_mean(self) -> float:
        return float(self.fl_models.rna.mean)

    @property
    def gdna_fl_mean(self) -> float:
        return float(self.fl_models.gdna.mean)

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "fl_models": self.fl_models.to_summary_dict(),
            "diagnostics": self.diagnostics.to_summary_dict(),
            "n_multi_loci": int(self.n_multi_loci),
            "region_calibration": _region_calibration_summary(self.region_calibration),
            "background_model": _background_model_summary(self.background_model),
            "boundary_local": _boundary_local_summary(self.boundary_local),
            "boundary_sweep": _boundary_sweep_summary(self.boundary_sweep),
            "strand_channels": _strand_channels_summary(self.strand_channels),
            "prior_table": (
                self.prior_table.to_summary_dict() if self.prior_table is not None else None
            ),
        }


def build_calibration_result(
    *,
    payload: CalibrationScanPayload,
    scan_trained: FragmentLengthModels,
    region_calibration: "RegionCalibration | None" = None,
    strand_channels: "RegionGdnaChannelEstimate | None" = None,
    background_model: "BackgroundModel | None" = None,
    boundary_local: "BoundaryLocalPosterior | None" = None,
    boundary_sweep: "BoundarySweepResult | None" = None,
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
    fl_models: FLModels | None = None,
    region_signature=None,
    prior_table: "PriorTable | None" = None,
) -> CalibrationResult:
    """Assemble the v6 calibration result around ``RegionCalibration``."""
    if region_calibration is None:
        raise ValueError("build_calibration_result: region_calibration is required.")
    if background_model is None:
        background_model = getattr(region_calibration, "background_model", None)
    if boundary_local is None:
        boundary_local = getattr(region_calibration, "boundary_local", None)
    if boundary_sweep is None:
        boundary_sweep = getattr(region_calibration, "boundary_sweep", None)
    if background_model is None:
        raise ValueError("build_calibration_result: background_model is required.")
    if boundary_local is None:
        raise ValueError("build_calibration_result: boundary_local is required.")
    if boundary_sweep is None:
        raise ValueError("build_calibration_result: boundary_sweep is required.")
    if fl_models is None:
        fl_models = build_fl_models(
            global_counts=extract_global_counts(scan_trained),
            rna_counts=extract_rna_counts(scan_trained),
            gdna_counts=extract_gdna_counts(payload),
            max_size=scan_trained.max_size,
            prior_ess=fl_prior_ess,
        )
    diagnostics = Diagnostics.from_payload(payload, signature=region_signature)
    return CalibrationResult(
        fl_models=fl_models,
        diagnostics=diagnostics,
        region_calibration=region_calibration,
        strand_channels=strand_channels,
        background_model=background_model,
        boundary_local=boundary_local,
        boundary_sweep=boundary_sweep,
        prior_table=prior_table,
        n_multi_loci=0,
    )


def _region_calibration_summary(region_calibration: "RegionCalibration") -> dict[str, object]:
    from .latent_states import STATE_NAMES

    states = np.asarray(region_calibration.p_states, dtype=np.float64)
    state_mass = {
        name: {
            "sum": float(np.sum(states[:, idx], dtype=np.float64)),
            "mean": float(np.mean(states[:, idx])) if states.shape[0] else 0.0,
        }
        for idx, name in enumerate(STATE_NAMES)
    }
    return {
        "n_regions": int(states.shape[0]),
        "rho_off": float(region_calibration.rho_off),
        "kappa_d": (
            None if region_calibration.kappa_d is None else float(region_calibration.kappa_d)
        ),
        "capture_enrichment_target": float(region_calibration.capture_enrichment_target),
        "n_passes": int(region_calibration.n_passes),
        "converged": bool(region_calibration.converged),
        "state_mass": state_mass,
        "p_expressed": _summary_stats(region_calibration.p_expressed),
        "p_captured": _summary_stats(region_calibration.p_captured),
        "mu_gdna": _summary_stats(region_calibration.mu_gdna),
        "upper_gdna": _summary_stats(region_calibration.upper_gdna),
        "rna_lower": _summary_stats(region_calibration.rna_lower),
        "prior_mass": _prior_mass_summary(region_calibration.prior_mass),
        "A_r": _summary_stats(region_calibration.A_r),
        "gamma_r": _summary_stats(region_calibration.gamma_r),
        "flag_histogram": _uint_histogram(region_calibration.flags),
        "pass_diagnostics": _json_safe(region_calibration.pass_diagnostics),
    }


def _prior_mass_summary(prior_mass) -> dict[str, object]:
    return {
        "unspliced_total": _summary_stats(prior_mass.unspliced_total),
        "gdna_unspliced_mean": _summary_stats(prior_mass.gdna_unspliced_mean),
        "rna_unspliced_mean": _summary_stats(prior_mass.rna_unspliced_mean),
        "max_abs_mass_conservation_error": float(
            np.max(prior_mass.mass_conservation_error())
        )
        if prior_mass.unspliced_total.size
        else 0.0,
        "method_histogram": _uint_histogram(prior_mass.method),
        "precision": _summary_stats(prior_mass.precision),
        "flag_histogram": _uint_histogram(prior_mass.flags),
    }


def _background_model_summary(background: "BackgroundModel") -> dict[str, object]:
    flags = np.asarray(background.flags, dtype=np.uint16)
    seed_mask = np.asarray(background.seed_mask, dtype=bool)
    top_t_mask = np.asarray(background.top_t_exclusion_mask, dtype=bool)
    return {
        "rho_off_alpha": float(background.rho_off_alpha),
        "rho_off_beta": float(background.rho_off_beta),
        "rho_off_mean": float(background.rho_off_mean),
        "n_seed_regions": int(background.n_seed_regions),
        "n_fragments": float(background.n_fragments),
        "eff_length": float(background.eff_length),
        "fit_status": str(background.fit_status),
        "n_regions": int(seed_mask.size),
        "n_seed_mask_true": int(np.count_nonzero(seed_mask)),
        "n_top_t_excluded": int(np.count_nonzero(top_t_mask)),
        "flag_histogram": _uint_histogram(flags),
    }


def _boundary_local_summary(local: "BoundaryLocalPosterior") -> dict[str, object]:
    alpha = np.asarray(local.alpha_excess)
    return {
        "n_regions": int(alpha.size),
        "alpha_excess": _summary_stats(local.alpha_excess),
        "beta_excess": _summary_stats(local.beta_excess),
        "mu_local": _summary_stats(local.mu_local),
        "upper_local": _summary_stats(local.upper_local),
        "n_regions_with_evidence": int(np.count_nonzero(alpha > 0.0)),
        "flag_histogram": _uint_histogram(local.flags),
    }


def _boundary_sweep_summary(sweep: "BoundarySweepResult") -> dict[str, object]:
    alpha = np.asarray(sweep.alpha_excess)
    return {
        "n_regions": int(alpha.size),
        "alpha_excess": _summary_stats(sweep.alpha_excess),
        "beta_excess": _summary_stats(sweep.beta_excess),
        "forward_alpha_excess": _summary_stats(sweep.forward_alpha_excess),
        "reverse_alpha_excess": _summary_stats(sweep.reverse_alpha_excess),
        "mu_sweep": _summary_stats(sweep.mu_sweep),
        "upper_sweep": _summary_stats(sweep.upper_sweep),
        "transfer_weight": _summary_stats(sweep.transfer_weight),
        "n_regions_with_swept_evidence": int(np.count_nonzero(alpha > 0.0)),
        "flag_histogram": _uint_histogram(sweep.flags),
    }


def _strand_channels_summary(
    strand_channels: "RegionGdnaChannelEstimate | None",
) -> dict[str, object] | None:
    if strand_channels is None:
        return None
    contained = np.asarray(strand_channels.contained_mean)
    return {
        "n_regions": int(contained.size),
        "internal_rna_lower_ci": float(strand_channels.internal_rna_lower_ci),
        "p_r1_sense": float(strand_channels.p_r1_sense),
        "kappa_d": float(strand_channels.kappa_d),
        "contained_mean": _summary_stats(strand_channels.contained_mean),
        "contained_upper": _summary_stats(strand_channels.contained_upper),
        "contained_rna_lower": _summary_stats(strand_channels.contained_rna_lower),
        "contained_precision": _summary_stats(strand_channels.contained_precision),
        "boundary_left_mean": _summary_stats(strand_channels.boundary_left_mean),
        "boundary_left_rna_lower": _summary_stats(strand_channels.boundary_left_rna_lower),
        "boundary_right_mean": _summary_stats(strand_channels.boundary_right_mean),
        "boundary_right_rna_lower": _summary_stats(strand_channels.boundary_right_rna_lower),
        "flag_histogram": _uint_histogram(strand_channels.flags),
    }


def _summary_stats(values: np.ndarray) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return {"min": 0.0, "p50": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0, "mean": 0.0}
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float64)
    return {
        "min": float(np.min(finite)),
        "p50": float(np.quantile(finite, 0.50)),
        "p95": float(np.quantile(finite, 0.95)),
        "p99": float(np.quantile(finite, 0.99)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
    }


def _uint_histogram(values: np.ndarray) -> dict[int, int]:
    arr = np.asarray(values)
    if arr.size == 0:
        return {}
    return {int(value): int(np.sum(arr == value)) for value in np.unique(arr)}


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return [_json_safe(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value
