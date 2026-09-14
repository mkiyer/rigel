"""CalibrationDiagnostics — the report-facing view of the fitted density landscape.

:class:`CalibrationResult` is the frozen prior the EM consumes; this is what the QC report may
*show*. It carries the total-density curve ``P(log ρ)`` on its grid, the two labelled dominant
modes, and their separation in nats. On a hybrid-capture library the density is bimodal — a low
"depleted" (off-target) mode and a high "enriched" (on-target) mode — and the separation is the log
enrichment factor.

Two deliberate limits. No categorical "capture worked" verdict is assigned; the threshold is the
analyst's. And the field described is a TOTAL density, so the panel reads ENRICHMENT and never the
gDNA/RNA split.

The calibrator builds one only when the landscape was fit (it needs wall inputs); otherwise it is
``None`` and the report omits the panel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CalibrationDiagnostics:
    """gDNA-density KDE curve + labeled modes for the report (no verdict)."""

    kde_x: np.ndarray  # log ρ_g grid
    kde_logp: np.ndarray  # log P̂ on the grid (the plottable curve)
    bandwidth: float
    n_eff: float
    n_modes: int  # total local maxima (descriptive; bimodal ⇒ likely capture)
    depleted_mode: float | None  # log ρ_g of the lower-density dominant mode
    enriched_mode: float | None  # log ρ_g of the higher-density dominant mode
    separation_nats: float | None  # enriched − depleted (0 if unimodal)
    enrichment_factor: float | None  # exp(separation_nats)
    rug_log_rho: (
        np.ndarray
    )  # per-region training log-densities — EVERY training region, not a sample
    rug_kind: np.ndarray  # int region-kind codes (0=intergenic,1=intron,2=exon,3=boundary)
    # ``rug_log_rho`` is the FULL training population (tens of thousands of rows at panel scale), not
    # a sample: `rigel report` writes it to `gdna_density_regions.feather` as a data export and no
    # chart spec reads it, so there is nothing to downsample for. Sample at the consumer if one ever
    # plots it directly.

    @classmethod
    def from_abundance_landscape(cls, al) -> "CalibrationDiagnostics":
        """Build from a fitted :class:`~rigel.calibration.abundance_landscape.AbundanceLandscape`.

        Every number here is READ FROM THE CENSUS rather than re-derived from the curve, so the mode
        labels mean something a curve alone cannot say: *depleted* is the basin containing the pooled
        intergenic anchor rate — an independent measurement of the same level — and *enriched* is the
        largest-mass basin strictly above it.

        ``bandwidth`` is the smoothing ACTUALLY IN FORCE — the grid step in decades — not a fitted
        kernel width. Nearly every per-region kernel is clamped to one grid step, so the knn width is
        not the resolution and reporting it would mislead.
        ``n_eff`` is the training-region count.

        ``separation_nats`` is the census's mode ratio and is the one field here that is
        RESOLUTION-SENSITIVE: it is displayed
        rather than consumed, while the depleted level beside it is grid-robust. A reader must not
        treat it as a calibrated enrichment factor."""
        ls = al.landscape
        x = np.asarray(ls.log_rho, dtype=np.float64)
        logp = np.asarray(ls.logP, dtype=np.float64)
        depleted = float(al.depleted.log_rho)
        enriched = depleted if al.enriched is None else float(al.enriched.log_rho)
        separation = enriched - depleted
        # the grid step in DECADES — the resolution the curve is rendered at, which is the smoothing
        # this estimator actually applies (see the docstring).
        step_dec = float(x[1] - x[0]) / float(np.log(10.0)) if x.size > 1 else 0.0
        return cls(
            kde_x=x,
            kde_logp=logp,
            bandwidth=step_dec,
            n_eff=float(al.n_train),
            n_modes=int(len(al.modes)),
            depleted_mode=depleted,
            enriched_mode=enriched,
            separation_nats=separation,
            enrichment_factor=float(np.exp(separation)),
            rug_log_rho=np.asarray(al.train_log_rho, dtype=np.float64),
            rug_kind=np.asarray(al.train_class, dtype=np.int64),
        )


__all__ = ["CalibrationDiagnostics"]
