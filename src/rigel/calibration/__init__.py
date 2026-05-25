"""Rigel gDNA calibration public surface."""

from ._diagnostics import Diagnostics
from ._exposure import l_eff_contained
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._orchestrator import calibrate
from ._result import CalibrationResult, build_calibration_result
from .density_model import DensityEvidence, GammaRatePrior, fit_density_evidence
from .exposure import RegionExposure
from .strand_balance import StrandBalanceEstimate, estimate_strand_balance
from .fl import (
    POOL_EB_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    FLModels,
    Quality,
    build_fl_models,
)
from .strand_summary import StrandSummary

__all__ = [
    "calibrate",
    "CalibrationResult",
    "build_calibration_result",
    "RegionExposure",
    "DensityEvidence",
    "GammaRatePrior",
    "fit_density_evidence",
    "StrandBalanceEstimate",
    "estimate_strand_balance",
    "l_eff_contained",
    "StrandSummary",
    "FLModels",
    "Quality",
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "POOL_EB_PRIOR_ESS",
    "build_fl_models",
    "Diagnostics",
    "extract_global_counts",
    "extract_rna_counts",
    "extract_gdna_counts",
]
