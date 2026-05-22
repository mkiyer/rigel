"""Typed errors for the Phase 3 fractional cutover.

The single ``FractionalCutoverPending`` exception is raised by every
code path that depends on the Phase 4 gDNA density / per-locus prior
estimator. Phase 4 will delete this module in a single grep-and-replace.
"""

from __future__ import annotations

from collections.abc import Mapping

#: Single shared message. All fail-fast sites raise with exactly this
#: text so users see one consistent error.
FRACTIONAL_CUTOVER_PENDING_MESSAGE = (
    "gDNA density and per-locus priors are pending the Phase 4 estimator "
    "(see docs/fineregions/fractional_accumulator_python_cutover.md \u00a75). "
    "`rigel quant` runs scan and calibration to completion and writes "
    "summary.json, then exits here. Use the calibration payload and "
    "FractionalEvidenceView for ad-hoc analysis in the meantime."
)


class FractionalCutoverPending(RuntimeError):
    """Phase 4 estimator dependency not yet implemented.

    Raised by ``compute_global_densities``, ``RegionalGdnaExposure.build``,
    ``assemble_priors``, ``quant_from_buffer`` (after calibration), and
    ``CalibrationResult.with_priors``. The CLI catches this and exits
    with code 64 (EX_USAGE) without a stack trace.
    """

    def __init__(
        self,
        message: str | None = None,
        *,
        calibration_summary: Mapping[str, object] | None = None,
    ) -> None:
        super().__init__(message or FRACTIONAL_CUTOVER_PENDING_MESSAGE)
        self.calibration_summary = (
            dict(calibration_summary) if calibration_summary is not None else None
        )


__all__ = ["FractionalCutoverPending", "FRACTIONAL_CUTOVER_PENDING_MESSAGE"]
