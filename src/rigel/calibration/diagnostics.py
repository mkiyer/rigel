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
    rug_log_rho: np.ndarray  # per-region training log-densities — EVERY training region, not a sample
    rug_kind: np.ndarray  # int region-kind codes (0=intergenic,1=intron,2=exon,3=boundary)
    # ⚠ ``rug_log_rho`` said "downsampled" until 2026-08-21, when it was always EMPTY (the npmle
    # carried no training points). It is now the full training population — ~30.7 k rows at panel
    # scale, which `rigel report` writes to `gdna_density_regions.feather` as a data export and no
    # chart spec reads, so there is nothing to downsample FOR. Sample at the consumer if one ever
    # plots it directly.

    @classmethod
    def from_abundance_landscape(cls, al) -> "CalibrationDiagnostics":
        """Build from a fitted :class:`~rigel.calibration.abundance_landscape.AbundanceLandscape`.

        ⭐⭐ **Every number here is READ FROM THE CENSUS rather than re-derived from the curve, and that
        is the improvement over the ``DensityNPMLE`` version this replaces.** That one took the two
        TALLEST local maxima and called them depleted and enriched; this takes the census's basins, where
        *depleted* is the basin containing the pooled intergenic ANCHOR rate — an independent
        measurement of the same level — and *enriched* is the largest-mass basin strictly above it. So
        the labels mean something a curve alone cannot say.

        ⭐ **The rug is REAL.** The npmle carried no per-region training points, so the report's rug was
        always empty and the CLI wrote a zero-row feather. The landscape publishes its training centres
        and their classes, so the panel can show the population under the fit.

        ⚠ **``bandwidth`` is the smoothing ACTUALLY IN FORCE — the grid step in decades — not a fitted
        kernel width.** On the benchmark ~99 % of this estimator's kernels are clamped to one grid step,
        so the per-region knn width is not the resolution and reporting it would mislead
        (`TRAPS: a-floored-knob-is-not-the-bandwidth`). ``n_eff`` is the training-region count, a real
        ``n`` where the npmle could only offer a collapsed-cell count.

        ⚠ **``separation_nats`` is the census's mode ratio, and it is the one field here that is
        RESOLUTION-SENSITIVE** (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`; the ruling is
        `DESIGN.md` §3.1a-iii). It is displayed rather than consumed, and the depleted level beside it
        is grid-robust — but a reader must not treat it as a calibrated enrichment factor.
        ⚠ And the field this describes is a TOTAL density, so the panel reads ENRICHMENT and never the
        gDNA split — composition-vacuous, exactly as the npmle version was."""
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
