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
from .fl import POOL_EB_PRIOR_ESS, POOL_SCORING_PRIOR_ESS, FLModels, build_fl_models
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
    Legacy density/fusion fields are intentionally absent. PR 04 exposes
    ``region_calibration.region_exposure`` and PR 05 projects it onto
    source-agnostic component EM effective lengths.
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
    fl_scoring_prior_ess: float = POOL_SCORING_PRIOR_ESS,
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
            scoring_prior_ess=fl_scoring_prior_ess,
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
        "kappa_d": (
            None if region_calibration.kappa_d is None else float(region_calibration.kappa_d)
        ),
        "n_passes": int(region_calibration.n_passes),
        "converged": bool(region_calibration.converged),
        "rho_off": float(region_calibration.rho_off),
        "state_mass": state_mass,
        "p_unexpressed": _summary_stats(region_calibration.p_unexpressed),
        "p_expressed": _summary_stats(region_calibration.p_expressed),
        "mu_gdna": _summary_stats(region_calibration.mu_gdna),
        "upper_gdna": _summary_stats(region_calibration.upper_gdna),
        "rna_lower": _summary_stats(region_calibration.rna_lower),
        "region_exposure": _region_exposure_summary(region_calibration.region_exposure),
        "flag_histogram": _uint_histogram(region_calibration.flags),
        "pass_diagnostics": _json_safe(region_calibration.pass_diagnostics),
        "region_unspliced_mass": _region_unspliced_mass_summary(
            region_calibration.region_unspliced_mass
        ),
        "background_density": _background_density_summary(
            region_calibration.background_density
        ),
    }


def _region_unspliced_mass_summary(rum) -> dict[str, object]:
    total = np.asarray(rum.total_mass, dtype=np.float64)
    gdna = np.asarray(rum.gdna_mass, dtype=np.float64)
    rna = np.asarray(rum.rna_mass, dtype=np.float64)
    conservation_err = (
        float(np.max(np.abs(gdna + rna - total))) if total.size else 0.0
    )
    return {
        "n_regions": int(total.size),
        "total_mass": _summary_stats(total),
        "gdna_mass": _summary_stats(gdna),
        "rna_mass": _summary_stats(rna),
        "region_size_bp": _summary_stats(np.asarray(rum.region_size_bp, dtype=np.float64)),
        "precision": _summary_stats(np.asarray(rum.precision, dtype=np.float64)),
        "unspliced_counts": _summary_stats(
            np.asarray(rum.unspliced_counts, dtype=np.float64)
        ),
        "method_histogram": _uint_histogram(np.asarray(rum.method, dtype=np.uint8)),
        "flag_histogram": _uint_histogram(np.asarray(rum.flags, dtype=np.uint16)),
        "max_abs_mass_conservation_error": conservation_err,
    }


def _background_density_summary(bg_density) -> dict[str, object]:
    return {
        "rho0_mean": float(bg_density.rho0_mean),
        "alpha0": float(bg_density.alpha0),
        "beta0": float(bg_density.beta0),
        "log_dispersion": float(bg_density.log_dispersion),
        "n_effective_regions": float(bg_density.n_effective_regions),
        "n_regions_in_pool": int(bg_density.n_regions_in_pool),
        "info_histogram": list(bg_density.info_histogram),
        "fit_status": str(bg_density.fit_status),
    }


def _region_exposure_summary(exposure) -> dict[str, object]:
    return {
        "n_regions": int(np.asarray(exposure.omega).size),
        "rho0": float(exposure.rho0),
        "tau2": float(exposure.tau2),
        "tau2_hat": float(exposure.tau2_hat),
        "tau2_method": str(exposure.tau2_method),
        "tau2_pool_size": int(exposure.tau2_pool_size),
        "omega": _summary_stats(exposure.omega),
        "log_omega": _summary_stats(exposure.log_omega),
        "raw_ratio": _summary_stats(exposure.raw_ratio),
        "log_raw_ratio": _summary_stats(exposure.log_raw_ratio),
        "shrink_weight": _summary_stats(exposure.shrink_weight),
        "v_obs": _summary_stats(exposure.v_obs),
        "lambda_global": _summary_stats(exposure.lambda_global),
        "support_count": _summary_stats(np.asarray(exposure.support_count, dtype=np.float64)),
        "flag_histogram": _uint_histogram(exposure.flags),
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
    from .strand_deconv import (
        FLAG_LOW_STRAND_RELIABILITY,
        FLAG_NEAR_UNSTRANDED,
        FLAG_RELIABILITY_APPROX,
    )

    contained = np.asarray(strand_channels.contained_mean)
    flags = np.asarray(strand_channels.flags, dtype=np.uint16)
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
        "contained_reliability": _summary_stats(strand_channels.contained_reliability),
        "contained_log_bayes_factor": _summary_stats(strand_channels.contained_log_bayes_factor),
        "boundary_left_reliability": _summary_stats(strand_channels.boundary_left_reliability),
        "boundary_left_log_bayes_factor": _summary_stats(
            strand_channels.boundary_left_log_bayes_factor
        ),
        "boundary_right_reliability": _summary_stats(strand_channels.boundary_right_reliability),
        "boundary_right_log_bayes_factor": _summary_stats(
            strand_channels.boundary_right_log_bayes_factor
        ),
        "n_regions_low_reliability": int(
            np.count_nonzero((flags & FLAG_LOW_STRAND_RELIABILITY) != 0)
        ),
        "n_regions_approx_reliability": int(
            np.count_nonzero((flags & FLAG_RELIABILITY_APPROX) != 0)
        ),
        "n_regions_near_unstranded": int(np.count_nonzero((flags & FLAG_NEAR_UNSTRANDED) != 0)),
        "flag_histogram": _uint_histogram(flags),
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
