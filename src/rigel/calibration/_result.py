"""Immutable calibration result for the Phase 1 cleanup path."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..frag_length_model import FragmentLengthModels
from ._diagnostics import Diagnostics
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from .density_global import GlobalDensityTable
from .fl import POOL_EB_PRIOR_ESS, FLModels, build_fl_models
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:  # pragma: no cover - Phase 3/7 types land later.
    from .exposure import RegionExposure
    from .strand_deconv import RegionGdnaEstimate


__all__ = ["CalibrationResult", "build_calibration_result"]


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Calibration hand-off produced before locus-level EM.

    Phase 1 deliberately removes the legacy prior and regional-exposure
    scaffolding. The two region fields are placeholders populated by later
    phases of the strand-deconvolution implementation.
    """

    global_densities: GlobalDensityTable | None
    fl_models: FLModels
    diagnostics: Diagnostics
    n_multi_loci: int = 0
    region_gdna: RegionGdnaEstimate | None = None
    region_exposure: RegionExposure | None = None

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
            "global_densities": (
                self.global_densities.to_summary_dict()
                if self.global_densities is not None
                else None
            ),
            "fl_models": self.fl_models.to_summary_dict(),
            "diagnostics": self.diagnostics.to_summary_dict(),
            "n_multi_loci": int(self.n_multi_loci),
            "region_gdna": None,
            "region_exposure": None,
        }


def build_calibration_result(
    *,
    payload: CalibrationScanPayload,
    scan_trained: FragmentLengthModels,
    global_densities: GlobalDensityTable | None = None,
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
    fl_models: FLModels | None = None,
    region_signature=None,
    region_gdna: "RegionGdnaEstimate | None" = None,
    region_exposure: "RegionExposure | None" = None,
) -> CalibrationResult:
    """Assemble the calibration result without legacy prior scaffolding."""
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
        global_densities=global_densities,
        fl_models=fl_models,
        diagnostics=diagnostics,
        n_multi_loci=0,
        region_gdna=region_gdna,
        region_exposure=region_exposure,
    )
