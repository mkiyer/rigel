"""CalibrationDiagnostics — report-facing diagnostics from the calibrator.

Distinct from :class:`CalibrationResult` (the frozen prior the EM consumes):
this carries what the QC report wants to *show*, chiefly the fitted gDNA-density
KDE ``P(log ρ_g)``. On a hybrid-capture library that density is bimodal — a low
"depleted" (off-target) mode and a high "enriched" (on-target) mode; their
separation in nats is the log enrichment factor. We surface the curve, the two
dominant modes, and the separation — but deliberately assign **no** categorical
"capture worked" verdict (that threshold is left to the analyst).

The calibrator builds one only when it actually fits the Phase-2 KDE (enough
training nodes); otherwise it is ``None`` and the report omits the panel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Cap on the number of training-node "rug" points carried for the report. Real
#: genomes have millions of teacher nodes; the rug only needs enough to show the
#: shape, so we stride-downsample to this many.
_RUG_CAP = 3000


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
    rug_log_rho: np.ndarray  # downsampled per-node training log-densities
    rug_kind: np.ndarray  # int node-kind codes (0=intergenic,1=intron,2=exon,3=boundary)

    @classmethod
    def from_prior(cls, prior) -> "CalibrationDiagnostics":
        """Build from a fitted :class:`GdnaDensityPrior`."""
        # prior.modes is all local maxima sorted by height (desc). Use the two
        # tallest as the dominant depleted/enriched pair, then order them by x.
        modes = list(prior.modes)
        top_x = sorted(float(m[0]) for m in modes[:2])
        if len(top_x) >= 2:
            depleted, enriched = top_x[0], top_x[1]
            separation = enriched - depleted
        elif len(top_x) == 1:
            depleted = enriched = top_x[0]
            separation = 0.0
        else:
            depleted = enriched = separation = None
        enrichment = float(np.exp(separation)) if separation is not None else None

        tx = np.asarray(prior.train_x, dtype=np.float64)
        tk = np.asarray(prior.train_kind)
        if tx.size > _RUG_CAP:
            idx = np.linspace(0, tx.size - 1, _RUG_CAP).astype(int)
            tx, tk = tx[idx], tk[idx]

        return cls(
            kde_x=np.asarray(prior.x_grid, dtype=np.float64),
            kde_logp=np.asarray(prior.logP_grid, dtype=np.float64),
            bandwidth=float(prior.bandwidth),
            n_eff=float(prior.n_eff),
            n_modes=len(modes),
            depleted_mode=depleted,
            enriched_mode=enriched,
            separation_nats=separation,
            enrichment_factor=enrichment,
            rug_log_rho=tx,
            rug_kind=tk.astype(np.int64),
        )


__all__ = ["CalibrationDiagnostics"]
