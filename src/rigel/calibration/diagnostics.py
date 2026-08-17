"""CalibrationDiagnostics — report-facing diagnostics from the calibrator.

Distinct from :class:`CalibrationResult` (the frozen prior the EM consumes):
this carries what the QC report wants to *show*, chiefly the fitted gDNA-density
KDE ``P(log ρ_g)``. On a hybrid-capture library that density is bimodal — a low
"depleted" (off-target) mode and a high "enriched" (on-target) mode; their
separation in nats is the log enrichment factor. We surface the curve, the two
dominant modes, and the separation — but deliberately assign **no** categorical
"capture worked" verdict (that threshold is left to the analyst).

The calibrator builds one only when it actually fits the Phase-2 KDE (enough
training regions); otherwise it is ``None`` and the report omits the panel.
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
    rug_log_rho: np.ndarray  # downsampled per-region training log-densities
    rug_kind: np.ndarray  # int region-kind codes (0=intergenic,1=intron,2=exon,3=boundary)

    @classmethod
    def from_prior(cls, prior) -> "CalibrationDiagnostics":
        """Build from a fitted :class:`~rigel.calibration.npmle.DensityNPMLE` enrichment landscape (Role A).
        Modes are the local maxima of the fitted log-density curve (it carries no per-region training points,
        so the rug is empty).

        ⛔ **It does NOT accept a**
        :class:`~rigel.calibration.landscape.DensityLandscape` — that landscape has no ``bandwidth`` and no
        ``n_cells``, so this raises ``AttributeError`` on one. This docstring claimed it did; the claim was
        never exercised, because `calibrate` only ever calls this with the enrichment prior.
        ⚠ **So the QC report's "bimodal ⇒ capture enrichment" caption is computed from the TOTAL-density
        landscape, which is composition-vacuous** — it reads enrichment, never the gDNA split. Pointing it at
        the gDNA hyperprior is a real change with a real audience and is left as an owner call, not smuggled
        in behind a docstring."""
        x = np.asarray(prior.log_rho, dtype=np.float64)
        logp = np.asarray(prior.logP, dtype=np.float64)
        # local maxima of the log-density curve, tallest first → the dominant depleted/enriched pair.
        interior = np.where((logp[1:-1] > logp[:-2]) & (logp[1:-1] >= logp[2:]))[0] + 1
        order = interior[np.argsort(logp[interior])[::-1]]
        top_x = sorted(float(x[i]) for i in order[:2])
        if len(top_x) >= 2:
            depleted, enriched = top_x[0], top_x[1]
            separation = enriched - depleted
        elif len(top_x) == 1:
            depleted = enriched = top_x[0]
            separation = 0.0
        else:
            depleted = enriched = separation = None
        enrichment = float(np.exp(separation)) if separation is not None else None
        empty = np.zeros(0)
        return cls(
            kde_x=x,
            kde_logp=logp,
            bandwidth=float(prior.bandwidth),
            n_eff=float(prior.n_cells),  # collapsed-cell count (the NPMLE has no Kish n_eff)
            n_modes=int(interior.size),
            depleted_mode=depleted,
            enriched_mode=enriched,
            separation_nats=separation,
            enrichment_factor=enrichment,
            rug_log_rho=empty,
            rug_kind=empty.astype(np.int64),
        )


__all__ = ["CalibrationDiagnostics"]
