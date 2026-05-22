"""rigel.calibration._regional_exposure \u2014 fail-fast stub (fractional cutover).

The legacy regional gDNA exposure estimator consumed the integer
per-region count matrix and the orientation channel split that the
fractional accumulator no longer emits. Until the fractional
replacement lands :meth:`RegionalGdnaExposure.build` raises
:class:`FractionalCutoverPending`. ``RegionalGdnaExposure.uniform`` is
preserved because the orchestrator still uses it as a degenerate
default while the prior pipeline is stubbed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from ._arrays import RegionArrays
from .errors import FractionalCutoverPending


__all__ = [
    "ExposureMode",
    "RegionalWeightApplicationStats",
    "RegionalGdnaExposure",
    "REFERENCE_QUANTILE",
]


ExposureMode = Literal["uniform", "kappa_eb"]

#: Default reference quantile for the ``A_r = exp(...)`` normalisation
#: anchor. Kept here so legacy summary writers compile.
REFERENCE_QUANTILE: float = 0.95


@dataclass(frozen=True, slots=True)
class RegionalWeightApplicationStats:
    """Per-call counters for the regional weight application step."""

    n_units_seen: int = 0
    n_units_weighted: int = 0
    n_units_no_gdna: int = 0
    n_units_missing_midpoint: int = 0
    n_units_cross_ref_skipped: int = 0

    def to_dict(self) -> dict[str, int]:
        return {
            "n_units_seen": int(self.n_units_seen),
            "n_units_weighted": int(self.n_units_weighted),
            "n_units_no_gdna": int(self.n_units_no_gdna),
            "n_units_missing_midpoint": int(self.n_units_missing_midpoint),
            "n_units_cross_ref_skipped": int(self.n_units_cross_ref_skipped),
        }


@dataclass(frozen=True, slots=True)
class RegionalGdnaExposure:
    """Per-region gDNA exposure weight ``A_r`` and lookup table.

    Dataclass shape preserved across the fractional cutover so that
    ``CalibrationResult`` and ``summary.json`` writers keep compiling.
    Only :meth:`uniform` is implemented; :meth:`build` raises.
    """

    rho_hat: np.ndarray
    log_weight: np.ndarray
    weight: np.ndarray
    mode: ExposureMode
    rho_ref: float
    n_at_floor: int
    reference_quantile: float = REFERENCE_QUANTILE
    per_class: dict[str, dict[str, float]] = field(default_factory=dict)
    rho_global: float = 0.0
    kappa_alpha_global: float = 0.0
    kappa_opportunity_bp: float = 0.0
    kappa_fallback_used: bool = False
    kappa_fallback_reason: str = ""
    observed_log_spread: float = 0.0
    null_log_spread: float = 0.0
    n_negative_rho_floored: int = 0

    ref_offsets: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    ref_id: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    start: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    end: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))

    @classmethod
    def uniform(cls, region_arrays: RegionArrays) -> "RegionalGdnaExposure":
        """Identity exposure: ``A_r == 1`` and ``log A_r == 0`` for all rows."""
        R = int(region_arrays.start.size)
        zeros = np.zeros(R, dtype=np.float64)
        ones = np.ones(R, dtype=np.float64)
        return cls(
            rho_hat=zeros,
            log_weight=zeros.copy(),
            weight=ones,
            mode="uniform",
            rho_ref=0.0,
            reference_quantile=REFERENCE_QUANTILE,
            n_at_floor=0,
            per_class={},
            rho_global=0.0,
            kappa_alpha_global=0.0,
            kappa_opportunity_bp=0.0,
            kappa_fallback_used=False,
            kappa_fallback_reason="",
            observed_log_spread=0.0,
            null_log_spread=0.0,
            n_negative_rho_floored=0,
            ref_offsets=np.asarray(region_arrays.ref_offsets, dtype=np.int32).copy(),
            ref_id=np.asarray(region_arrays.ref_id, dtype=np.int32).copy(),
            start=np.asarray(region_arrays.start, dtype=np.int64).copy(),
            end=np.asarray(region_arrays.end, dtype=np.int64).copy(),
        )

    @classmethod
    def build(cls, *_args: object, **_kwargs: object) -> "RegionalGdnaExposure":
        """Removed under the fractional cutover. Raises :class:`FractionalCutoverPending`."""
        raise FractionalCutoverPending(
            "RegionalGdnaExposure.build: regional gDNA exposure consumes the "
            "integer per-region count matrix and orientation channels that "
            "the fractional accumulator no longer emits. The replacement "
            "fractional exposure pipeline has not yet landed."
        )

    def to_summary_dict(self) -> dict[str, object]:
        """Minimal summary block consumed by ``CalibrationResult.to_summary_dict``.

        Under the cutover the orchestrator only emits a uniform exposure
        field, so the summary records the mode + bookkeeping only.
        """
        return {
            "mode": str(self.mode),
            "n_regions": int(self.weight.size),
            "rho_ref": float(self.rho_ref),
            "rho_global": float(self.rho_global),
            "reference_quantile": float(self.reference_quantile),
            "n_at_floor": int(self.n_at_floor),
            "n_negative_rho_floored": int(self.n_negative_rho_floored),
            "kappa_opportunity_bp": float(self.kappa_opportunity_bp),
            "kappa_alpha_global": float(self.kappa_alpha_global),
            "kappa_fallback_used": bool(self.kappa_fallback_used),
            "kappa_fallback_reason": str(self.kappa_fallback_reason),
            "observed_log_spread": float(self.observed_log_spread),
            "null_log_spread": float(self.null_log_spread),
            "per_class": dict(self.per_class),
        }
