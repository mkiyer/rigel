"""Rigel gDNA calibration (v4: mappability-aware EM mixture deconvolution)."""

from ._calibrate import calibrate_gdna
from ._em import EMFit, run_em
from ._fl_model import build_gdna_fl_model
from ._result import CalibrationResult
from ._stats import compute_region_stats, compute_sense_fraction


__all__ = [
    "CalibrationResult",
    "EMFit",
    "build_gdna_fl_model",
    "calibrate_gdna",
    "compute_region_stats",
    "compute_sense_fraction",
    "run_em",
]
