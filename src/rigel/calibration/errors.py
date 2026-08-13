"""Calibration error types."""

from __future__ import annotations


class CalibrationSubstrateError(ValueError):
    """Raised when the calibration substrate is malformed or misaligned.

    Examples: the region geometry does not boundary up 1:1 with the accumulator
    payload, or a payload array has an unexpected shape.
    """


class CalibrationStrandError(ValueError):
    """Raised when the library carries no spliced reads to train the strand model.

    The RNA strand orientation (``rna_sense_frac``) is identifiable only from spliced unique
    mappers — gDNA cannot splice, so spliced reads are the unambiguous-RNA anchor.
    With zero spliced observations the strand clue degenerates to 0.5 (symmetric),
    indistinguishable from gDNA, and the deconvolution cannot separate sense RNA from
    contamination. A real RNA-seq library always has spliced reads; a library without
    any is out of scope, so we fail loudly rather than silently mis-deconvolve.
    """


__all__ = ["CalibrationSubstrateError", "CalibrationStrandError"]
