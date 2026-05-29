"""Rigel calibration — joint fractional-accumulator + calibration-v6.

This package recovers the library-wide gDNA contamination model and the
per-region deconvoluted mass from the accumulator payload. ``rigel quant``
calls :func:`calibrate`; the result feeds :func:`assemble_priors`, which
produces the per-locus EM priors.

The package is being rebuilt PR-by-PR (see
``docs/acc_caljointmodel/00_implementation_plan.md``); the inference bodies
(:func:`calibrate`, :class:`CalibrationSubstrate`, :func:`assemble_priors`)
are skeletons that raise :class:`NotImplementedError` until their PR lands.

This module is pure re-exports; implementations live in the submodules.
"""

from __future__ import annotations

from ..config import CalibrationConfig
from .calibrate import calibrate
from .errors import CalibrationConvergenceError, CalibrationSubstrateError
from .priors import assemble_priors
from .result import CalibrationResult
from .substrate import CalibrationSubstrate

__all__ = [
    "CalibrationConfig",
    "CalibrationResult",
    "CalibrationSubstrate",
    "CalibrationSubstrateError",
    "CalibrationConvergenceError",
    "calibrate",
    "assemble_priors",
]
