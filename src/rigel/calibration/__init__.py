"""Rigel gDNA calibration (v4: mappability-aware EM mixture deconvolution)."""

from ._annotate import annotate_capture_class
from ._calibrate import calibrate_gdna
from ._em import EMFit, run_em
from ._fl_model import build_gdna_fl_model
from ._result import CalibrationResult
from ._stats import compute_region_stats, compute_sense_fraction


__all__ = [
    "CalibrationResult",
    "EMFit",
    "annotate_capture_class",
    "build_gdna_fl_model",
    "calibrate_gdna",
    "compute_region_stats",
    "compute_sense_fraction",
    "run_em",
]
