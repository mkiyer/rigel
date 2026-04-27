"""Rigel gDNA calibration (SRD v1: Simple Regional Deconvolution).

The public surface is intentionally tiny:

* :class:`CalibrationResult` — dataclass returned by ``calibrate_gdna``.
* :func:`calibrate_gdna` — single buffer-walk orchestrator implementing
  Pass 0 (geometric categorization) → Pass 1 (pool assembly) → Pass 2
  (1-D RNA/gDNA mixture) → Pass 3 (Empirical-Bayes fragment-length
  models).  See ``docs/calibration/srd_v1_implementation.md``.
"""

from ._result import CalibrationResult, GdnaFlQuality
from ._simple import calibrate_gdna


__all__ = [
    "CalibrationResult",
    "GdnaFlQuality",
    "calibrate_gdna",
]
