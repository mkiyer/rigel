"""Calibration E/M iteration for the v6 four-state region model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from ._arrays import RegionArrays
from .background_model import BackgroundModel
from .boundaries import BoundaryTable
from .boundary_model import (
    BoundaryLocalPosterior,
    build_boundary_local_posterior,
    predict_contained_gdna_from_excess,
)
from .boundary_sweep import BoundarySweepResult, run_boundary_sweep
from .density_observation import DensityObservation
from .latent_states import (
    N_STATES,
    STATE_UNEXPRESSED_CAPTURE,
    STATE_UNEXPRESSED_OFFTARGET,
    STATE_IS_CAPTURED,
    STATE_IS_EXPRESSED,
    build_logbf_capture,
    build_logbf_expression,
    build_logbf_gdna_density,
    build_logbf_strand,
    build_state_log_prior,
    build_state_log_tensor,
    normalize_state_log_tensor,
)
from .strand_deconv import RegionGdnaChannelEstimate

__all__ = [
    "FLAG_STATE_AMBIGUOUS",
    "FLAG_STRAND_UNINFORMATIVE",
    "FLAG_BOUNDARY_SPARSE",
    "FLAG_EXPRESSED_UNCERTAIN",
    "FLAG_EXACT_STRAND_POSTERIOR",
    "PRIOR_MASS_METHOD_DENSITY",
    "PRIOR_MASS_METHOD_STRAND",
    "PriorMassDeconvolution",
    "RegionCalibration",
    "CalibrationStepResult",
    "build_prior_mass_deconvolution",
    "calibration_e_step",
    "calibration_m_step",
    "run_calibration_iteration",
]

FLAG_STATE_AMBIGUOUS: int = 1 << 0
FLAG_STRAND_UNINFORMATIVE: int = 1 << 1
FLAG_BOUNDARY_SPARSE: int = 1 << 2
FLAG_EXPRESSED_UNCERTAIN: int = 1 << 3
FLAG_EXACT_STRAND_POSTERIOR: int = 1 << 4

PRIOR_MASS_METHOD_DENSITY: int = 1
PRIOR_MASS_METHOD_STRAND: int = 2

_EPS: float = 1.0e-12
_MAX_EXPOSURE_VALUE: float = 1.0e30


@dataclass(frozen=True, slots=True)
class PriorMassDeconvolution:
    """Mass-conserving unspliced RNA/gDNA prior evidence per region."""

    unspliced_total: np.ndarray
    gdna_unspliced_mean: np.ndarray
    rna_unspliced_mean: np.ndarray
    method: np.ndarray
    precision: np.ndarray
    flags: np.ndarray

    def __post_init__(self) -> None:
        total = np.asarray(self.unspliced_total, dtype=np.float32)
        region_count = int(total.shape[0])
        if total.ndim != 1:
            raise ValueError(f"unspliced_total must be 1D; got shape {total.shape}.")
        object.__setattr__(self, "unspliced_total", total)

        for field_name in ("gdna_unspliced_mean", "rna_unspliced_mean", "precision"):
            values = _as_float32_vector(field_name, getattr(self, field_name), region_count)
            object.__setattr__(self, field_name, values)

        method = np.asarray(self.method, dtype=np.uint8)
        if method.shape != (region_count,):
            raise ValueError(f"method must have shape ({region_count},); got {method.shape}.")
        object.__setattr__(self, "method", method)

        flags = np.asarray(self.flags, dtype=np.uint16)
        if flags.shape != (region_count,):
            raise ValueError(f"flags must have shape ({region_count},); got {flags.shape}.")
        object.__setattr__(self, "flags", flags)

        if np.any(total < 0.0):
            raise ValueError("unspliced_total must be non-negative.")
        if np.any(self.gdna_unspliced_mean < 0.0) or np.any(self.rna_unspliced_mean < 0.0):
            raise ValueError("prior mass means must be non-negative.")
        if not np.all(np.isfinite(self.precision)) or np.any(self.precision < 0.0):
            raise ValueError("precision must be finite and non-negative.")
        if not np.allclose(
            self.gdna_unspliced_mean + self.rna_unspliced_mean,
            total,
            rtol=1.0e-5,
            atol=1.0e-5,
        ):
            raise ValueError("prior mass must conserve unspliced_total per region.")

    def mass_conservation_error(self) -> np.ndarray:
        """Return per-region absolute conservation residual."""
        return np.abs(
            self.gdna_unspliced_mean + self.rna_unspliced_mean - self.unspliced_total
        ).astype(np.float32, copy=False)


@dataclass(frozen=True, slots=True)
class RegionCalibration:
    """Region-level four-state calibration output for downstream v6 wiring."""

    p_states: np.ndarray
    mu_gdna: np.ndarray
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    prior_mass: PriorMassDeconvolution
    A_r: np.ndarray
    gamma_r: np.ndarray
    rho_off: float
    kappa_d: float | None
    capture_enrichment_target: float
    n_passes: int
    converged: bool
    flags: np.ndarray
    pass_diagnostics: tuple[dict[str, object], ...]
    background_model: BackgroundModel | None = None
    boundary_local: BoundaryLocalPosterior | None = None
    boundary_sweep: BoundarySweepResult | None = None

    def __post_init__(self) -> None:
        p_states = np.asarray(self.p_states, dtype=np.float32)
        if p_states.ndim != 2 or p_states.shape[1] != N_STATES:
            raise ValueError(f"p_states must have shape (R, {N_STATES}); got {p_states.shape}.")
        region_count = int(p_states.shape[0])
        row_sums = p_states.sum(axis=1)
        if not np.allclose(row_sums, 1.0, rtol=1.0e-5, atol=1.0e-5):
            raise ValueError("p_states rows must sum to 1 within tolerance.")
        object.__setattr__(self, "p_states", p_states)

        for field_name in ("mu_gdna", "upper_gdna", "rna_lower", "A_r", "gamma_r"):
            values = _as_float32_vector(field_name, getattr(self, field_name), region_count)
            object.__setattr__(self, field_name, values)

        if np.any(self.mu_gdna < 0.0):
            raise ValueError("mu_gdna must be non-negative.")
        if np.any(self.upper_gdna + 1.0e-5 < self.mu_gdna):
            raise ValueError("upper_gdna must be >= mu_gdna within tolerance.")
        if np.any(self.rna_lower < 0.0):
            raise ValueError("rna_lower must be non-negative.")
        prior_mass = self.prior_mass
        if not isinstance(prior_mass, PriorMassDeconvolution):
            prior_mass = PriorMassDeconvolution(
                unspliced_total=getattr(prior_mass, "unspliced_total"),
                gdna_unspliced_mean=getattr(prior_mass, "gdna_unspliced_mean"),
                rna_unspliced_mean=getattr(prior_mass, "rna_unspliced_mean"),
                method=getattr(prior_mass, "method"),
                precision=getattr(prior_mass, "precision"),
                flags=getattr(prior_mass, "flags"),
            )
        if prior_mass.unspliced_total.shape != (region_count,):
            raise ValueError(
                "prior_mass arrays must match p_states region count; "
                f"got {prior_mass.unspliced_total.shape}, expected {(region_count,)}."
            )
        object.__setattr__(self, "prior_mass", prior_mass)
        for field_name in ("A_r", "gamma_r"):
            values = np.asarray(getattr(self, field_name), dtype=np.float32)
            if not np.all(np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError(f"{field_name} must be finite and non-negative.")

        flags = np.asarray(self.flags, dtype=np.uint16)
        if flags.shape != (region_count,):
            raise ValueError(f"flags must have shape ({region_count},); got {flags.shape}.")
        object.__setattr__(self, "flags", flags)

        if not np.isfinite(self.rho_off) or float(self.rho_off) < 0.0:
            raise ValueError(f"rho_off must be finite and non-negative; got {self.rho_off!r}.")
        if self.kappa_d is not None and (
            not np.isfinite(self.kappa_d) or float(self.kappa_d) <= 0.0
        ):
            raise ValueError(f"kappa_d must be None or finite and positive; got {self.kappa_d!r}.")
        if (
            not np.isfinite(self.capture_enrichment_target)
            or float(self.capture_enrichment_target) < 1.0
        ):
            raise ValueError(
                "capture_enrichment_target must be finite and >= 1; "
                f"got {self.capture_enrichment_target!r}."
            )
        if int(self.n_passes) < 1:
            raise ValueError(f"n_passes must be >= 1; got {self.n_passes!r}.")
        diagnostics = tuple(dict(item) for item in self.pass_diagnostics)
        object.__setattr__(self, "pass_diagnostics", diagnostics)

    @property
    def p_unexpressed_offtarget(self) -> np.ndarray:
        return self.p_states[:, STATE_UNEXPRESSED_OFFTARGET]

    @property
    def p_unexpressed_capture(self) -> np.ndarray:
        return self.p_states[:, STATE_UNEXPRESSED_CAPTURE]

    @property
    def p_expressed_capture(self) -> np.ndarray:
        return self.p_states[:, 2]

    @property
    def p_expressed_offtarget(self) -> np.ndarray:
        return self.p_states[:, 3]

    @property
    def p_expressed(self) -> np.ndarray:
        return self.p_states[:, STATE_IS_EXPRESSED].sum(axis=1)

    @property
    def p_captured(self) -> np.ndarray:
        return self.p_states[:, STATE_IS_CAPTURED].sum(axis=1)


@dataclass(frozen=True, slots=True)
class CalibrationStepResult:
    """Single E-step output before scalar M-step refits."""

    p_states: np.ndarray
    mu_gdna: np.ndarray
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    prior_mass: PriorMassDeconvolution
    A_r: np.ndarray
    gamma_r: np.ndarray
    flags: np.ndarray
    local_posterior: BoundaryLocalPosterior
    sweep: BoundarySweepResult
    log_tensor: np.ndarray
    sum_log_evidence: float

    @property
    def p_expressed(self) -> np.ndarray:
        return self.p_states[:, STATE_IS_EXPRESSED].sum(axis=1)

    @property
    def p_captured(self) -> np.ndarray:
        return self.p_states[:, STATE_IS_CAPTURED].sum(axis=1)



def _as_float32_vector(name: str, values: np.ndarray, region_count: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float32)
    if array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _as_float64_vector(name: str, values: np.ndarray, region_count: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _validate_region_inputs(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> int:
    region_count = int(np.asarray(observation.contained_leff).shape[0])
    if np.asarray(region_arrays.signature).shape != (region_count,):
        raise ValueError(
            "region_arrays.signature must match observation region count; "
            f"got {np.asarray(region_arrays.signature).shape}, expected {(region_count,)}."
        )
    if int(boundaries.ref_region_offsets[-1]) != region_count:
        raise ValueError(
            "boundaries and observation disagree on region count; "
            f"got {int(boundaries.ref_region_offsets[-1])} and {region_count}."
        )
    if np.asarray(background.seed_mask).shape != (region_count,):
        raise ValueError(
            "background.seed_mask must match observation region count; "
            f"got {np.asarray(background.seed_mask).shape}, expected {(region_count,)}."
        )
    if strand_channels is not None:
        _as_float64_vector("strand_channels.contained_mean", strand_channels.contained_mean, region_count)
    return region_count


def _validate_local_posterior(local_posterior: BoundaryLocalPosterior, region_count: int) -> None:
    _as_float64_vector("local_posterior.alpha_excess", local_posterior.alpha_excess, region_count)
    _as_float64_vector("local_posterior.beta_excess", local_posterior.beta_excess, region_count)


def _safe_exposure(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    ratio = np.divide(
        np.asarray(numerator, dtype=np.float64),
        np.maximum(np.asarray(denominator, dtype=np.float64), _EPS),
        out=np.zeros_like(np.asarray(numerator, dtype=np.float64)),
        where=np.asarray(denominator, dtype=np.float64) >= 0.0,
    )
    ratio = np.nan_to_num(ratio, nan=0.0, posinf=_MAX_EXPOSURE_VALUE, neginf=0.0)
    ratio = np.clip(ratio, 0.0, _MAX_EXPOSURE_VALUE)
    return ratio.astype(np.float32)


def build_prior_mass_deconvolution(
    observation: DensityObservation,
    *,
    mu_gdna: np.ndarray,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> PriorMassDeconvolution:
    """Build mass-conserving unspliced RNA/gDNA prior evidence per region."""
    unspliced_total = _as_float64_vector(
        "observation.observed_compatible_count",
        observation.observed_compatible_count,
        int(np.asarray(observation.observed_compatible_count).shape[0]),
    )
    region_count = int(unspliced_total.shape[0])
    density_mu = _as_float64_vector("mu_gdna", mu_gdna, region_count)

    method = np.full(region_count, PRIOR_MASS_METHOD_DENSITY, dtype=np.uint8)
    precision = np.zeros(region_count, dtype=np.float32)
    flags = np.zeros(region_count, dtype=np.uint16)

    if strand_channels is None:
        gdna = density_mu.copy()
    else:
        contained = _as_float64_vector(
            "strand_channels.contained_mean", strand_channels.contained_mean, region_count
        )
        left = _as_float64_vector(
            "strand_channels.boundary_left_mean", strand_channels.boundary_left_mean, region_count
        )
        right = _as_float64_vector(
            "strand_channels.boundary_right_mean", strand_channels.boundary_right_mean, region_count
        )
        gdna = contained + left + right
        method.fill(PRIOR_MASS_METHOD_STRAND)
        precision = np.maximum.reduce(
            [
                np.asarray(strand_channels.contained_precision, dtype=np.float32),
                np.asarray(strand_channels.boundary_left_precision, dtype=np.float32),
                np.asarray(strand_channels.boundary_right_precision, dtype=np.float32),
            ]
        ).astype(np.float32, copy=False)
        flags = np.asarray(strand_channels.flags, dtype=np.uint16).copy()

    gdna = np.clip(np.nan_to_num(gdna, nan=0.0, posinf=0.0, neginf=0.0), 0.0, unspliced_total)
    # Preserve exact conservation after float32 conversion by deriving RNA from
    # the float32 total and gDNA arrays.
    total32 = unspliced_total.astype(np.float32)
    gdna32 = gdna.astype(np.float32)
    rna32 = np.maximum(total32 - gdna32, 0.0).astype(np.float32)

    return PriorMassDeconvolution(
        unspliced_total=total32,
        gdna_unspliced_mean=gdna32,
        rna_unspliced_mean=rna32,
        method=method,
        precision=precision,
        flags=flags,
    )


def _derive_region_flags(
    p_states: np.ndarray,
    local_posterior: BoundaryLocalPosterior,
    strand_channels: RegionGdnaChannelEstimate | None,
) -> np.ndarray:
    region_count = int(p_states.shape[0])
    flags = np.zeros(region_count, dtype=np.uint16)
    flags[np.max(p_states, axis=1) < 0.6] |= FLAG_STATE_AMBIGUOUS

    p_expressed = p_states[:, STATE_IS_EXPRESSED].sum(axis=1)
    uncertain_expression = (p_expressed > 0.3) & (p_expressed < 0.7)
    flags[uncertain_expression] |= FLAG_EXPRESSED_UNCERTAIN

    sparse_boundary = (np.asarray(local_posterior.alpha_excess) <= 0.0) | (
        np.asarray(local_posterior.beta_excess) <= 0.0
    )
    flags[sparse_boundary] |= FLAG_BOUNDARY_SPARSE

    if strand_channels is None or not np.isfinite(strand_channels.kappa_d):
        flags |= FLAG_STRAND_UNINFORMATIVE
    return flags


def calibration_e_step(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    local_posterior: BoundaryLocalPosterior | None = None,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    pass_index: int = 0,
    confidence: float = 0.95,
    background_boost: float = 1.0,
    transfer_weight: np.ndarray | None = None,
) -> CalibrationStepResult:
    """Run one four-state calibration E-step."""
    region_count = _validate_region_inputs(
        region_arrays,
        observation,
        boundaries,
        background,
        strand_channels,
    )
    if local_posterior is None:
        local_posterior = build_boundary_local_posterior(
            observation,
            boundaries,
            background,
            strand_channels=strand_channels,
            confidence=confidence,
        )
    else:
        _validate_local_posterior(local_posterior, region_count)

    sweep = run_boundary_sweep(
        local_posterior,
        boundaries,
        observation,
        background,
        transfer_weight=transfer_weight,
        strand_channels=strand_channels,
        confidence=confidence,
    )
    state_log_prior = build_state_log_prior(
        region_arrays,
        background,
        pass_index=pass_index,
        background_boost=background_boost,
    )
    logbf_expression = build_logbf_expression(observation, strand_channels)
    logbf_capture = build_logbf_capture(sweep, observation, background)
    logbf_density = build_logbf_gdna_density(observation, background, sweep, strand_channels)
    logbf_strand = build_logbf_strand(strand_channels, region_count)
    log_tensor = build_state_log_tensor(
        state_log_prior,
        logbf_expression,
        logbf_capture,
        logbf_density,
        logbf_strand,
    )
    p_states = normalize_state_log_tensor(log_tensor)

    contained_leff = _as_float64_vector(
        "observation.contained_leff", observation.contained_leff, region_count
    )
    p_captured = p_states[:, STATE_IS_CAPTURED].sum(axis=1).astype(np.float64)
    off_target_mu = np.maximum(float(background.rho_off_mean), 0.0) * np.maximum(
        contained_leff,
        0.0,
    )
    upper_off_target = predict_contained_gdna_from_excess(
        background,
        contained_leff,
        np.zeros(region_count, dtype=np.float64),
        np.zeros(region_count, dtype=np.float64),
        confidence=confidence,
    )[1].astype(np.float64, copy=False)
    captured_mu = _as_float64_vector("sweep.mu_sweep", sweep.mu_sweep, region_count)
    captured_upper = _as_float64_vector("sweep.upper_sweep", sweep.upper_sweep, region_count)

    mu_gdna = (p_captured * captured_mu + (1.0 - p_captured) * off_target_mu).astype(np.float32)
    upper_gdna = (
        p_captured * captured_upper + (1.0 - p_captured) * upper_off_target
    ).astype(np.float32)
    upper_gdna = np.maximum(upper_gdna, mu_gdna)

    if strand_channels is None:
        rna_lower = np.zeros(region_count, dtype=np.float32)
    else:
        rna_lower = np.maximum(
            _as_float64_vector(
                "strand_channels.contained_rna_lower",
                strand_channels.contained_rna_lower,
                region_count,
            ),
            0.0,
        ).astype(np.float32)

    prior_mass = build_prior_mass_deconvolution(
        observation,
        mu_gdna=mu_gdna,
        strand_channels=strand_channels,
    )

    denominator = np.maximum(float(background.rho_off_mean), 0.0) * np.maximum(contained_leff, 0.0)
    exposure = _safe_exposure(mu_gdna, denominator)
    captured_exposure = _safe_exposure(captured_mu, denominator)
    flags = _derive_region_flags(p_states, local_posterior, strand_channels)
    sum_log_evidence = float(np.sum(logsumexp(log_tensor, axis=1), dtype=np.float64))

    return CalibrationStepResult(
        p_states=p_states,
        mu_gdna=mu_gdna,
        upper_gdna=upper_gdna,
        rna_lower=rna_lower,
        prior_mass=prior_mass,
        A_r=exposure,
        gamma_r=captured_exposure,
        flags=flags,
        local_posterior=local_posterior,
        sweep=sweep,
        log_tensor=log_tensor,
        sum_log_evidence=sum_log_evidence,
    )


def calibration_m_step(
    observation: DensityObservation,
    background: BackgroundModel,
    p_states: np.ndarray,
    sweep: BoundarySweepResult,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    damping: float = 0.5,
    capture_enrichment_target: float = 1.0,
    alpha_floor: float = 1.0,
    beta_floor: float = 1.0,
) -> tuple[BackgroundModel, float | None, float]:
    """Refit scalar calibration parameters from state posteriors."""
    if not np.isfinite(damping) or not 0.0 <= float(damping) <= 1.0:
        raise ValueError(f"damping must be finite and in [0, 1]; got {damping!r}.")
    if not np.isfinite(capture_enrichment_target) or float(capture_enrichment_target) < 1.0:
        raise ValueError(
            "capture_enrichment_target must be finite and >= 1; "
            f"got {capture_enrichment_target!r}."
        )
    if not np.isfinite(alpha_floor) or float(alpha_floor) <= 0.0:
        raise ValueError(f"alpha_floor must be finite and > 0; got {alpha_floor!r}.")
    if not np.isfinite(beta_floor) or float(beta_floor) <= 0.0:
        raise ValueError(f"beta_floor must be finite and > 0; got {beta_floor!r}.")

    states = np.asarray(p_states, dtype=np.float64)
    if states.ndim != 2 or states.shape[1] != N_STATES:
        raise ValueError(f"p_states must have shape (R, {N_STATES}); got {states.shape}.")
    region_count = int(states.shape[0])
    contained_leff = _as_float64_vector(
        "observation.contained_leff", observation.contained_leff, region_count
    )
    if strand_channels is None:
        gdna_count = _as_float64_vector(
            "observation.contained_count", observation.contained_count, region_count
        )
        kappa_d = None
    else:
        gdna_count = _as_float64_vector(
            "strand_channels.contained_mean", strand_channels.contained_mean, region_count
        )
        kappa_d = float(strand_channels.kappa_d)

    seed_mask = np.asarray(background.seed_mask, dtype=bool)
    if seed_mask.shape != (region_count,):
        raise ValueError(
            "background.seed_mask must match region count; "
            f"got {seed_mask.shape}, expected {(region_count,)}."
        )
    unexpressed_offtarget_weight = (
        states[:, STATE_UNEXPRESSED_OFFTARGET] * seed_mask.astype(np.float64)
    )
    weighted_fragments = float(
        np.sum(unexpressed_offtarget_weight * gdna_count, dtype=np.float64)
    )
    weighted_eff_length = float(
        np.sum(unexpressed_offtarget_weight * contained_leff, dtype=np.float64)
    )
    alpha_hat = float(alpha_floor) + weighted_fragments
    beta_hat = float(beta_floor) + weighted_eff_length
    eta = float(damping)
    alpha_next = (1.0 - eta) * float(background.rho_off_alpha) + eta * alpha_hat
    beta_next = (1.0 - eta) * float(background.rho_off_beta) + eta * beta_hat
    beta_next = max(beta_next, np.finfo(np.float64).tiny)
    rho_next = alpha_next / beta_next

    n_seed_regions = int(np.count_nonzero((unexpressed_offtarget_weight > 0.0) & seed_mask))
    if n_seed_regions == 0 or weighted_eff_length <= 0.0:
        fit_status = "prior_only"
    elif background.fit_status == "ok":
        fit_status = "ok"
    else:
        fit_status = "sparse"

    next_background = BackgroundModel(
        rho_off_alpha=float(alpha_next),
        rho_off_beta=float(beta_next),
        rho_off_mean=float(rho_next),
        seed_mask=seed_mask.copy(),
        top_t_exclusion_mask=np.asarray(background.top_t_exclusion_mask, dtype=bool).copy(),
        n_seed_regions=n_seed_regions,
        n_fragments=weighted_fragments,
        eff_length=weighted_eff_length,
        fit_status=fit_status,
        flags=np.asarray(background.flags, dtype=np.uint16).copy(),
    )

    captured_weight = states[:, STATE_IS_CAPTURED].sum(axis=1)
    off_target_mu = np.maximum(float(background.rho_off_mean), 0.0) * np.maximum(
        contained_leff,
        0.0,
    )
    captured_mu = _as_float64_vector("sweep.mu_sweep", sweep.mu_sweep, region_count)
    denominator = float(np.sum(captured_weight * off_target_mu, dtype=np.float64))
    if denominator > 0.0:
        ratio_hat = float(np.sum(captured_weight * captured_mu, dtype=np.float64) / denominator)
        ratio_hat = max(1.0, ratio_hat)
        next_capture_target = (1.0 - eta) * float(capture_enrichment_target) + eta * ratio_hat
    else:
        next_capture_target = float(capture_enrichment_target)
    next_capture_target = max(1.0, float(next_capture_target))

    return next_background, kappa_d, next_capture_target


def _relative_change(current_value: float, previous_value: float | None) -> float:
    if previous_value is None:
        return float("inf")
    denominator = max(abs(float(previous_value)), _EPS)
    return abs(float(current_value) - float(previous_value)) / denominator


def run_calibration_iteration(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    local_posterior: BoundaryLocalPosterior | None = None,
    transfer_weight: np.ndarray | None = None,
    max_calibration_passes: int = 5,
    p_tol: float = 0.01,
    rho_tol: float = 0.02,
    damping: float = 0.5,
    confidence: float = 0.95,
    background_boost: float = 1.0,
    capture_enrichment_target: float = 1.0,
) -> RegionCalibration:
    """Run the additive v6 four-state calibration loop."""
    if int(max_calibration_passes) < 1:
        raise ValueError(
            "max_calibration_passes must be >= 1; "
            f"got {max_calibration_passes!r}."
        )
    if not np.isfinite(p_tol) or float(p_tol) < 0.0:
        raise ValueError(f"p_tol must be finite and >= 0; got {p_tol!r}.")
    if not np.isfinite(rho_tol) or float(rho_tol) < 0.0:
        raise ValueError(f"rho_tol must be finite and >= 0; got {rho_tol!r}.")

    current_background = background
    current_capture_target = float(capture_enrichment_target)
    current_kappa = float(strand_channels.kappa_d) if strand_channels is not None else None
    previous_p_states: np.ndarray | None = None
    previous_rho: float | None = None
    diagnostics: list[dict[str, object]] = []
    last_step: CalibrationStepResult | None = None
    converged = False

    for pass_index in range(int(max_calibration_passes)):
        step = calibration_e_step(
            region_arrays,
            observation,
            boundaries,
            current_background,
            local_posterior=local_posterior,
            strand_channels=strand_channels,
            pass_index=pass_index,
            confidence=confidence,
            background_boost=background_boost,
            transfer_weight=transfer_weight,
        )
        if previous_p_states is None:
            max_state_shift = float("inf")
        else:
            max_state_shift = float(np.max(np.abs(step.p_states - previous_p_states)))
        rho_shift = _relative_change(current_background.rho_off_mean, previous_rho)
        converged = bool(
            previous_p_states is not None
            and max_state_shift < float(p_tol)
            and rho_shift < float(rho_tol)
        )
        diagnostics.append(
            {
                "pass_index": int(pass_index),
                "max_state_shift": float(max_state_shift),
                "rho_off": float(current_background.rho_off_mean),
                "relative_rho_shift": float(rho_shift),
                "kappa_d": current_kappa,
                "capture_enrichment_target": float(current_capture_target),
                "sum_log_evidence": float(step.sum_log_evidence),
                "converged": bool(converged),
                "n_regions_expressed": int(np.count_nonzero(step.p_expressed > 0.5)),
                "n_regions_captured": int(np.count_nonzero(step.p_captured > 0.5)),
                "n_regions_unexpressed_offtarget": int(
                    np.count_nonzero(step.p_states[:, STATE_UNEXPRESSED_OFFTARGET] > 0.5)
                ),
            }
        )
        last_step = step
        if converged or pass_index == int(max_calibration_passes) - 1:
            break

        previous_p_states = step.p_states.copy()
        previous_rho = float(current_background.rho_off_mean)
        current_background, current_kappa, current_capture_target = calibration_m_step(
            observation,
            current_background,
            step.p_states,
            step.sweep,
            strand_channels,
            damping=damping,
            capture_enrichment_target=current_capture_target,
        )

    if last_step is None:  # pragma: no cover - max_calibration_passes guard prevents this.
        raise RuntimeError("run_calibration_iteration: no calibration pass was executed.")

    return RegionCalibration(
        p_states=last_step.p_states,
        mu_gdna=last_step.mu_gdna,
        upper_gdna=last_step.upper_gdna,
        rna_lower=last_step.rna_lower,
        prior_mass=last_step.prior_mass,
        A_r=last_step.A_r,
        gamma_r=last_step.gamma_r,
        rho_off=float(current_background.rho_off_mean),
        kappa_d=current_kappa,
        capture_enrichment_target=float(current_capture_target),
        n_passes=len(diagnostics),
        converged=bool(converged),
        flags=last_step.flags,
        pass_diagnostics=tuple(diagnostics),
        background_model=current_background,
        boundary_local=last_step.local_posterior,
        boundary_sweep=last_step.sweep,
    )
