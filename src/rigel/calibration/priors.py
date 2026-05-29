"""assemble_priors() — bridge from CalibrationResult to per-locus EM priors.

The D3 prior bridge turns the calibration output (exposure posterior + library
hyperparameters) into the per-locus Dirichlet prior scalars the locus EM
consumes, via exposure-weighted effective lengths. The exact inputs/outputs
are firmed up in PR 6; this is a PR 1 skeleton so the public surface exists.
See ``docs/acc_caljointmodel/00_implementation_plan.md``.
"""

from __future__ import annotations


def assemble_priors(*args, **kwargs):
    raise NotImplementedError(
        "rigel.calibration.assemble_priors is a PR 1 skeleton; the CalibrationResult → "
        "per-locus prior bridge lands in PR 6 — see "
        "docs/acc_caljointmodel/00_implementation_plan.md."
    )


__all__ = ["assemble_priors"]
