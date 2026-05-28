"""Per-region empirical-Bayes gDNA exposure factors."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - imported for annotations only.
    from .calibration_iteration import BackgroundDensity, RegionUnsplicedMass

__all__ = [
    "FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL",
    "FLAG_EXPOSURE_IMPUTED_TIER3",
    "FLAG_EXPOSURE_NOT_TAU2_POOL",
    "FLAG_EXPOSURE_NO_SUPPORT",
    "FLAG_EXPOSURE_NUMERIC_CEILING",
    "FLAG_EXPOSURE_NUMERIC_FLOOR",
    "RegionExposure",
    "estimate_region_exposure",
]

FLAG_EXPOSURE_NO_SUPPORT: int = 1 << 0
FLAG_EXPOSURE_NOT_TAU2_POOL: int = 1 << 1
FLAG_EXPOSURE_IMPUTED_TIER3: int = 1 << 2
FLAG_EXPOSURE_NUMERIC_FLOOR: int = 1 << 3
FLAG_EXPOSURE_NUMERIC_CEILING: int = 1 << 4
FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL: int = 1 << 5

_METHOD_STRAND: int = 1
_METHOD_BOUNDARY: int = 2
_METHOD_BACKGROUND_FALLBACK: int = 3
_REAL_TAU2_METHODS: frozenset[str] = frozenset({"moment", "moment_damped"})
_TINY: float = np.finfo(np.float64).tiny
_MIN_EXPOSURE_POOL_P_UNEXPRESSED: float = 0.80


@dataclass(frozen=True, slots=True)
class RegionExposure:
    """Per-region multiplicative gDNA sampling exposure.

    PR 04 produces this object for diagnostics and future PR 05 denominator
    normalization. PR 04 itself keeps locus EM exposure-neutral.
    """

    omega: np.ndarray
    log_omega: np.ndarray
    raw_ratio: np.ndarray
    log_raw_ratio: np.ndarray
    shrink_weight: np.ndarray
    v_obs: np.ndarray
    lambda_global: np.ndarray
    rho0: float
    tau2: float
    tau2_hat: float
    support_count: np.ndarray
    tau2_pool_size: int
    tau2_method: str
    flags: np.ndarray

    _VALID_METHODS: ClassVar[frozenset[str]] = frozenset(
        {"bootstrap_neutral", "no_pool_neutral", "moment", "moment_damped"}
    )

    def __post_init__(self) -> None:
        omega = _as_float64_array("omega", self.omega)
        region_count = int(omega.shape[0])
        object.__setattr__(self, "omega", omega)

        for name in (
            "log_omega",
            "raw_ratio",
            "log_raw_ratio",
            "shrink_weight",
            "v_obs",
            "lambda_global",
        ):
            values = _as_float64_array(name, getattr(self, name), region_count)
            object.__setattr__(self, name, values)

        support = np.ascontiguousarray(np.asarray(self.support_count, dtype=np.uint64))
        if support.shape != (region_count,):
            raise ValueError(
                f"support_count must have shape ({region_count},); got {support.shape}."
            )
        object.__setattr__(self, "support_count", support)

        flags = np.ascontiguousarray(np.asarray(self.flags, dtype=np.uint16))
        if flags.shape != (region_count,):
            raise ValueError(f"flags must have shape ({region_count},); got {flags.shape}.")
        object.__setattr__(self, "flags", flags)

        rho0 = float(self.rho0)
        tau2 = float(self.tau2)
        tau2_hat = float(self.tau2_hat)
        pool_size = int(self.tau2_pool_size)
        method = str(self.tau2_method)
        if not np.isfinite(rho0) or rho0 <= 0.0:
            raise ValueError(f"rho0 must be finite and > 0; got {self.rho0!r}.")
        if not np.isfinite(tau2) or tau2 < 0.0:
            raise ValueError(f"tau2 must be finite and >= 0; got {self.tau2!r}.")
        if not np.isfinite(tau2_hat) or tau2_hat < 0.0:
            raise ValueError(f"tau2_hat must be finite and >= 0; got {self.tau2_hat!r}.")
        if pool_size < 0:
            raise ValueError(f"tau2_pool_size must be >= 0; got {self.tau2_pool_size!r}.")
        if method not in self._VALID_METHODS:
            raise ValueError(
                f"tau2_method must be one of {sorted(self._VALID_METHODS)}; got {method!r}."
            )
        object.__setattr__(self, "rho0", rho0)
        object.__setattr__(self, "tau2", tau2)
        object.__setattr__(self, "tau2_hat", tau2_hat)
        object.__setattr__(self, "tau2_pool_size", pool_size)
        object.__setattr__(self, "tau2_method", method)

        if not np.all(np.isfinite(self.omega)) or np.any(self.omega <= 0.0):
            raise ValueError("omega must be finite and strictly positive.")
        if not np.all(np.isfinite(self.log_omega)):
            raise ValueError("log_omega must be finite.")
        if not np.allclose(self.log_omega, np.log(self.omega), rtol=1.0e-10, atol=1.0e-12):
            raise ValueError("log_omega must be consistent with log(omega).")
        if not np.all(np.isfinite(self.raw_ratio)) or np.any(self.raw_ratio <= 0.0):
            raise ValueError("raw_ratio must be finite and strictly positive.")
        if not np.all(np.isfinite(self.log_raw_ratio)):
            raise ValueError("log_raw_ratio must be finite.")
        if not np.allclose(
            self.log_raw_ratio,
            np.log(self.raw_ratio),
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ValueError("log_raw_ratio must be consistent with log(raw_ratio).")
        if (
            not np.all(np.isfinite(self.shrink_weight))
            or np.any(self.shrink_weight < 0.0)
            or np.any(self.shrink_weight > 1.0)
        ):
            raise ValueError("shrink_weight must be finite and in [0, 1].")
        if not np.all(np.isfinite(self.v_obs)) or np.any(self.v_obs < 0.0):
            raise ValueError("v_obs must be finite and non-negative.")
        if not np.all(np.isfinite(self.lambda_global)) or np.any(self.lambda_global < 0.0):
            raise ValueError("lambda_global must be finite and non-negative.")


def estimate_region_exposure(
    region_unspliced_mass: "RegionUnsplicedMass",
    background_density: "BackgroundDensity",
    p_unexpressed: np.ndarray,
    *,
    previous: RegionExposure | None = None,
    alpha_floor: float = 1.0,
    tau2_damping: float = 0.5,
    winsorize_k: float = 4.0,
    no_support_v_obs: float = 1.0e6,
    omega_floor: float = 1.0e-6,
    omega_ceiling: float = 1.0e6,
) -> RegionExposure:
    """Estimate per-region exposure by shrinking raw density ratios toward one."""
    _validate_scalar_params(
        alpha_floor=alpha_floor,
        tau2_damping=tau2_damping,
        winsorize_k=winsorize_k,
        no_support_v_obs=no_support_v_obs,
        omega_floor=omega_floor,
        omega_ceiling=omega_ceiling,
    )

    gdna = _as_float64_array("region_unspliced_mass.gdna_mass", region_unspliced_mass.gdna_mass)
    region_count = int(gdna.shape[0])
    region_bp = _as_float64_array(
        "region_unspliced_mass.region_size_bp",
        region_unspliced_mass.region_size_bp,
        region_count,
    )
    precision = _as_float64_array(
        "region_unspliced_mass.precision",
        region_unspliced_mass.precision,
        region_count,
    )
    support = np.ascontiguousarray(
        np.asarray(region_unspliced_mass.unspliced_counts, dtype=np.uint64)
    )
    if support.shape != (region_count,):
        raise ValueError(f"unspliced_counts must have shape ({region_count},); got {support.shape}.")
    method = np.ascontiguousarray(np.asarray(region_unspliced_mass.method, dtype=np.uint8))
    if method.shape != (region_count,):
        raise ValueError(f"method must have shape ({region_count},); got {method.shape}.")
    p_unx = _as_float64_array("p_unexpressed", p_unexpressed, region_count)
    if not np.all(np.isfinite(p_unx)):
        raise ValueError("p_unexpressed must be finite.")
    p_unx = np.clip(p_unx, 0.0, 1.0)

    rho0 = float(background_density.rho0_mean)
    if not np.isfinite(rho0) or rho0 <= 0.0:
        raise ValueError(f"background_density.rho0_mean must be finite and > 0; got {rho0!r}.")
    if np.any(region_bp <= 0.0):
        raise ValueError("region_size_bp must be strictly positive.")

    alpha = float(alpha_floor)
    pseudo_size = alpha / rho0
    rho_hat = (gdna + alpha) / (region_bp + pseudo_size)
    raw_ratio = np.maximum(rho_hat / rho0, _TINY)
    log_raw_ratio = np.log(raw_ratio)
    lambda_global = rho0 * region_bp

    support_float = support.astype(np.float64, copy=False)
    v_obs = np.full(region_count, float(no_support_v_obs), dtype=np.float64)
    positive_support = support > np.uint64(0)
    v_obs[positive_support] = 1.0 / support_float[positive_support]

    flags = np.zeros(region_count, dtype=np.uint16)
    flags[~positive_support] |= np.uint16(FLAG_EXPOSURE_NO_SUPPORT)
    tier3_mask = method == np.uint8(_METHOD_BACKGROUND_FALLBACK)
    flags[tier3_mask] |= np.uint16(FLAG_EXPOSURE_IMPUTED_TIER3)

    if _is_bootstrap_density(background_density):
        flags |= np.uint16(FLAG_EXPOSURE_BOOTSTRAP_NEUTRAL | FLAG_EXPOSURE_NOT_TAU2_POOL)
        return _neutral_exposure(
            raw_ratio=raw_ratio,
            log_raw_ratio=log_raw_ratio,
            v_obs=v_obs,
            lambda_global=lambda_global,
            rho0=rho0,
            support=support,
            flags=flags,
            tau2_method="bootstrap_neutral",
        )

    pool_mask = (
        ((method == np.uint8(_METHOD_STRAND)) | (method == np.uint8(_METHOD_BOUNDARY)))
        & (region_bp >= 1.0)
        & positive_support
        & (p_unx >= _MIN_EXPOSURE_POOL_P_UNEXPRESSED)
    )
    weights = precision * support_float * p_unx
    weights = np.where(pool_mask, weights, 0.0)
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    weights = np.maximum(weights, 0.0)
    active_pool = weights > 0.0
    flags[~active_pool] |= np.uint16(FLAG_EXPOSURE_NOT_TAU2_POOL)
    pool_size = int(np.count_nonzero(active_pool))

    if pool_size == 0:
        return _neutral_exposure(
            raw_ratio=raw_ratio,
            log_raw_ratio=log_raw_ratio,
            v_obs=v_obs,
            lambda_global=lambda_global,
            rho0=rho0,
            support=support,
            flags=flags,
            tau2_method="no_pool_neutral",
        )

    dispersion = max(float(background_density.log_dispersion), _TINY)
    clip_radius = float(winsorize_k) * dispersion
    y_fit = np.clip(log_raw_ratio, -clip_radius, clip_radius)
    empirical_var = _weighted_mean(y_fit * y_fit, weights)
    mean_v_obs = _weighted_mean(v_obs, weights)
    tau2_hat = max(empirical_var - mean_v_obs, 0.0)

    if previous is not None and previous.tau2_method in _REAL_TAU2_METHODS:
        tau2 = (1.0 - float(tau2_damping)) * float(previous.tau2) + float(
            tau2_damping
        ) * tau2_hat
        tau2_method = "moment_damped"
    else:
        tau2 = tau2_hat
        tau2_method = "moment"

    if tau2 <= 0.0:
        shrink_weight = np.zeros(region_count, dtype=np.float64)
    else:
        shrink_weight = tau2 / (tau2 + v_obs)
    shrink_weight = np.asarray(shrink_weight, dtype=np.float64)
    shrink_weight[~active_pool] = 0.0

    log_omega_unclipped = shrink_weight * log_raw_ratio
    log_omega_unclipped[~active_pool] = 0.0
    log_floor = float(np.log(omega_floor))
    log_ceiling = float(np.log(omega_ceiling))
    clipped_low = log_omega_unclipped < log_floor
    clipped_high = log_omega_unclipped > log_ceiling
    log_omega = np.clip(log_omega_unclipped, log_floor, log_ceiling)
    omega = np.exp(log_omega)
    omega[clipped_low] = float(omega_floor)
    omega[clipped_high] = float(omega_ceiling)
    log_omega = np.log(omega)
    flags[clipped_low] |= np.uint16(FLAG_EXPOSURE_NUMERIC_FLOOR)
    flags[clipped_high] |= np.uint16(FLAG_EXPOSURE_NUMERIC_CEILING)

    return RegionExposure(
        omega=omega,
        log_omega=log_omega,
        raw_ratio=raw_ratio,
        log_raw_ratio=log_raw_ratio,
        shrink_weight=shrink_weight,
        v_obs=v_obs,
        lambda_global=lambda_global,
        rho0=rho0,
        tau2=float(tau2),
        tau2_hat=float(tau2_hat),
        support_count=support,
        tau2_pool_size=pool_size,
        tau2_method=tau2_method,
        flags=flags,
    )


def _neutral_exposure(
    *,
    raw_ratio: np.ndarray,
    log_raw_ratio: np.ndarray,
    v_obs: np.ndarray,
    lambda_global: np.ndarray,
    rho0: float,
    support: np.ndarray,
    flags: np.ndarray,
    tau2_method: str,
) -> RegionExposure:
    region_count = int(raw_ratio.shape[0])
    return RegionExposure(
        omega=np.ones(region_count, dtype=np.float64),
        log_omega=np.zeros(region_count, dtype=np.float64),
        raw_ratio=raw_ratio,
        log_raw_ratio=log_raw_ratio,
        shrink_weight=np.zeros(region_count, dtype=np.float64),
        v_obs=v_obs,
        lambda_global=lambda_global,
        rho0=float(rho0),
        tau2=0.0,
        tau2_hat=0.0,
        support_count=support,
        tau2_pool_size=0,
        tau2_method=tau2_method,
        flags=flags,
    )


def _is_bootstrap_density(background_density: "BackgroundDensity") -> bool:
    method_histogram = tuple(int(x) for x in getattr(background_density, "method_histogram", ()))
    return (
        str(background_density.fit_status) == "fallback_bootstrap"
        and int(background_density.n_regions_in_pool) == 0
        and sum(method_histogram) == 0
    )


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    w_sum = float(np.sum(weights, dtype=np.float64))
    if w_sum <= 0.0:
        return 0.0
    return float(np.sum(np.asarray(values, dtype=np.float64) * weights, dtype=np.float64) / w_sum)


def _as_float64_array(
    name: str,
    values: np.ndarray,
    region_count: int | None = None,
) -> np.ndarray:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    if array.ndim != 1:
        raise ValueError(f"{name} must be 1D; got shape {array.shape}.")
    if region_count is not None and array.shape != (region_count,):
        raise ValueError(f"{name} must have shape ({region_count},); got {array.shape}.")
    return array


def _validate_scalar_params(
    *,
    alpha_floor: float,
    tau2_damping: float,
    winsorize_k: float,
    no_support_v_obs: float,
    omega_floor: float,
    omega_ceiling: float,
) -> None:
    if not np.isfinite(alpha_floor) or float(alpha_floor) <= 0.0:
        raise ValueError(f"alpha_floor must be finite and > 0; got {alpha_floor!r}.")
    if not np.isfinite(tau2_damping) or not 0.0 <= float(tau2_damping) <= 1.0:
        raise ValueError(f"tau2_damping must be finite and in [0, 1]; got {tau2_damping!r}.")
    if not np.isfinite(winsorize_k) or float(winsorize_k) <= 0.0:
        raise ValueError(f"winsorize_k must be finite and > 0; got {winsorize_k!r}.")
    if not np.isfinite(no_support_v_obs) or float(no_support_v_obs) <= 0.0:
        raise ValueError(
            f"no_support_v_obs must be finite and > 0; got {no_support_v_obs!r}."
        )
    if (
        not np.isfinite(omega_floor)
        or not np.isfinite(omega_ceiling)
        or not 0.0 < float(omega_floor) < 1.0 < float(omega_ceiling)
    ):
        raise ValueError(
            "omega_floor/omega_ceiling must satisfy 0 < floor < 1 < ceiling; "
            f"got {omega_floor!r}, {omega_ceiling!r}."
        )
