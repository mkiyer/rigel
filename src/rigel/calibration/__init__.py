"""Rigel calibration — joint fractional-accumulator + calibration-v6.

This package recovers the library-wide gDNA contamination model and the
per-region deconvoluted mass from the accumulator payload. ``rigel quant``
calls :func:`calibrate`; the result feeds :func:`assemble_priors`, which
produces the per-locus EM priors.

This module is pure re-exports; implementations live in the submodules.
"""

from __future__ import annotations

from ..config import CalibrationConfig
from .calibrate import calibrate
from .errors import CalibrationSubstrateError
from .priors import LocusPriors, assemble_priors
from .result import CalibrationResult
from .strand_balance import StrandBalance
from .substrate import CalibrationSubstrate, SubstrateView

__all__ = [
    "CalibrationConfig",
    "CalibrationResult",
    "CalibrationSubstrate",
    "SubstrateView",
    "StrandBalance",
    "CalibrationSubstrateError",
    "calibrate",
    "assemble_priors",
    "LocusPriors",
]
