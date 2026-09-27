"""CalibrationDiagnostics — the exported view of the fitted total-density landscape.

:class:`CalibrationResult` is the frozen prior the EM consumes; this is what ``rigel quant`` writes
beside it: the total-density curve ``P(log ρ)`` on its grid (``gdna_density_kde.feather``) and the
training population it was fitted on (``gdna_density_regions.feather``).

The calibrator builds one only when the landscape was fit (it needs wall inputs); otherwise it is
``None`` and neither file is written.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CalibrationDiagnostics:
    """The total-density curve and its training population."""

    kde_x: np.ndarray  # log ρ grid
    kde_logp: np.ndarray  # log P̂ on the grid (the plottable curve)
    #: per-region training log-densities — EVERY training region, not a sample
    rug_log_rho: np.ndarray
    rug_kind: np.ndarray  # int region-kind codes (0=intergenic,1=intron,2=exon)

    @classmethod
    def from_abundance_landscape(cls, al) -> "CalibrationDiagnostics":
        """Build from a fitted :class:`~rigel.calibration.abundance_landscape.AbundanceLandscape`."""
        ls = al.landscape
        return cls(
            kde_x=np.asarray(ls.log_rho, dtype=np.float64),
            kde_logp=np.asarray(ls.logP, dtype=np.float64),
            rug_log_rho=np.asarray(al.train_log_rho, dtype=np.float64),
            rug_kind=np.asarray(al.train_class, dtype=np.int64),
        )


__all__ = ["CalibrationDiagnostics"]
