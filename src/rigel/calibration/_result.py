"""Immutable calibration result for the v5 calibration path."""

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
from .density_global import GlobalDensityTable
from .fl import POOL_EB_PRIOR_ESS, FLModels, build_fl_models
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:  # pragma: no cover - imported for annotations only.
    from .exposure import RegionExposure
    from .strand_deconv import RegionGdnaEstimate


__all__ = ["CalibrationResult", "build_calibration_result"]


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Calibration hand-off produced before locus-level EM.

    Phase 4 populates strand-deconvolved per-region gDNA counts and a uniform
    region-exposure surface. Locus prior assembly consumes these in Phase 5.
    """

    global_densities: GlobalDensityTable | None
    fl_models: FLModels
    diagnostics: Diagnostics
    region_gdna: "RegionGdnaEstimate"
    region_exposure: "RegionExposure"
    n_multi_loci: int = 0
    rna_lower_confidence: float = 0.95

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
            "calibration_config": {
                "rna_lower_confidence": float(self.rna_lower_confidence),
            },
            "global_densities": (
                self.global_densities.to_summary_dict()
                if self.global_densities is not None
                else None
            ),
            "fl_models": self.fl_models.to_summary_dict(),
            "diagnostics": self.diagnostics.to_summary_dict(),
            "n_multi_loci": int(self.n_multi_loci),
            "strand_deconv": _strand_deconv_summary(self.region_gdna),
            "region_exposure": self.region_exposure.to_summary_dict(),
        }


def _strand_deconv_summary(region_gdna: "RegionGdnaEstimate") -> dict[str, object]:
    """Return the compact Phase 4 strand-deconvolution summary block."""
    from .strand_deconv import (
        FLAG_APPROX_NORMAL,
        FLAG_INELIGIBLE,
        FLAG_KAPPA_FALLBACK,
        FLAG_NEAR_UNSTRANDED,
    )

    flags = np.asarray(region_gdna.flags, dtype=np.uint8)
    return {
        "rna_lower_confidence": float(region_gdna.rna_lower_confidence),
        "p_r1_sense": float(region_gdna.p_r1_sense),
        "kappa_d": float(region_gdna.kappa_d),
        "kappa_d_n_seed_regions": int(region_gdna.kappa_d_n_seed_regions),
        "kappa_d_n_exon_self_training": int(region_gdna.kappa_d_n_exon_self_training),
        "kappa_d_fallback_used": bool(np.any((flags & FLAG_KAPPA_FALLBACK) != 0)),
        "n_regions": int(flags.size),
        "n_regions_eligible": int(np.sum((flags & FLAG_INELIGIBLE) == 0)),
        "n_regions_approx_normal": int(np.sum((flags & FLAG_APPROX_NORMAL) != 0)),
        "n_regions_near_unstranded": int(np.sum((flags & FLAG_NEAR_UNSTRANDED) != 0)),
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
    rna_lower_confidence: float = 0.95,
) -> CalibrationResult:
    """Assemble the calibration result without legacy prior scaffolding."""
    if region_gdna is None:
        raise ValueError("build_calibration_result: region_gdna is required in Phase 4.")
    if region_exposure is None:
        raise ValueError("build_calibration_result: region_exposure is required in Phase 4.")
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
        rna_lower_confidence=float(rna_lower_confidence),
    )
