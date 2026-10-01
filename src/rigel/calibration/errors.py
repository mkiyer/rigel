"""Calibration error types."""

from __future__ import annotations


class CalibrationSubstrateError(ValueError):
    """Raised when the calibration substrate is malformed or misaligned.

    Examples: the region geometry does not line up 1:1 with the accumulator
    payload, or a payload array has an unexpected shape.
    """


__all__ = ["CalibrationSubstrateError"]
