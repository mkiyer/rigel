"""Calibration error types (doc ``04_interface_contract.md`` §9)."""

from __future__ import annotations


class CalibrationSubstrateError(ValueError):
    """Raised when the calibration substrate is malformed or misaligned.

    Examples: the region geometry does not line up 1:1 with the accumulator
    payload, or a payload array has an unexpected shape.
    """


class CalibrationConvergenceError(RuntimeError):
    """Raised when the calibration outer loop fails to converge."""


__all__ = ["CalibrationSubstrateError", "CalibrationConvergenceError"]
