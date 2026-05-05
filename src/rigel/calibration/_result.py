"""rigel.calibration._result — immutable calibration result.

The frozen :class:`CalibrationResult` is the public hand-off between
:func:`rigel.calibration.calibrate` (the v6 orchestrator) and
:func:`rigel.pipeline.quant_from_buffer` (which consumes
``fl_models``/``global_densities`` and back-fills the per-locus
:class:`PriorTable` via :meth:`CalibrationResult.with_priors`).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..frag_length_model import FragmentLengthModels
from ._diagnostics import Diagnostics
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from .density_global import GlobalDensityTable
from .fl import POOL_EB_PRIOR_ESS, FLModels, build_fl_models
from .locus_prior import LocusGdnaEstimate, MultiLocusPrior, PriorTable
from .scan_payload import CalibrationScanPayload


__all__ = [
    "CalibrationResult",
    "build_calibration_result",
    "build_multi_locus_prior_df",
    "build_per_locus_gdna_df",
]


# ---------------------------------------------------------------------------
# Diagnostic-dataframe builders (locked column order; pinned by tests)
# ---------------------------------------------------------------------------

_MULTI_LOCUS_COLUMNS: tuple[str, ...] = (
    "multi_locus_id", "n_obs", "n_gdna", "n_rna", "pi_gdna", "n_loci",
)

_PER_LOCUS_COLUMNS: tuple[str, ...] = (
    "multi_locus_id", "ref", "start", "end", "span",
    "n_obs", "n_gdna",
    "n_gdna_intergenic", "n_gdna_intron", "n_gdna_exon_intron",
    "pi_gdna", "n_eligible_boundaries", "fallback_flags",
)


def build_multi_locus_prior_df(
    mlps: tuple[MultiLocusPrior, ...],
) -> pd.DataFrame:
    if not mlps:
        return pd.DataFrame({c: [] for c in _MULTI_LOCUS_COLUMNS})
    return pd.DataFrame(
        {
            "multi_locus_id": [m.multi_locus_id for m in mlps],
            "n_obs":          [m.n_obs          for m in mlps],
            "n_gdna":         [m.n_gdna         for m in mlps],
            "n_rna":          [m.n_rna          for m in mlps],
            "pi_gdna":        [m.pi_gdna        for m in mlps],
            "n_loci":         [len(m.per_locus) for m in mlps],
        },
        columns=list(_MULTI_LOCUS_COLUMNS),
    )


def build_per_locus_gdna_df(
    mlps: tuple[MultiLocusPrior, ...],
) -> pd.DataFrame:
    if not mlps:
        return pd.DataFrame({c: [] for c in _PER_LOCUS_COLUMNS})
    rows: list[dict[str, object]] = []
    for ml in mlps:
        e: LocusGdnaEstimate
        for e in ml.per_locus:
            rows.append(
                {
                    "multi_locus_id":        ml.multi_locus_id,
                    "ref":                   e.locus.ref,
                    "start":                 e.locus.start,
                    "end":                   e.locus.end,
                    "span":                  e.locus.span,
                    "n_obs":                 e.n_obs,
                    "n_gdna":                e.n_gdna,
                    "n_gdna_intergenic":     e.n_gdna_intergenic,
                    "n_gdna_intron":         e.n_gdna_intron,
                    "n_gdna_exon_intron":    e.n_gdna_exon_intron,
                    "pi_gdna":               e.pi_gdna,
                    "n_eligible_boundaries": e.n_eligible_boundaries,
                    "fallback_flags":        e.fallback_flags,
                }
            )
    return pd.DataFrame(rows, columns=list(_PER_LOCUS_COLUMNS))


# ---------------------------------------------------------------------------
# CalibrationResult schema
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Immutable v6 calibration result."""

    global_densities: GlobalDensityTable    # M4
    fl_models:        FLModels              # M7 — sole finalized FL surface
    prior_table:      PriorTable            # M6

    diagnostics:  Diagnostics               # named breakdown of n_observed
    n_multi_loci: int

    # Eager diagnostic dataframes (locked schemas)
    multi_locus_prior_df: pd.DataFrame
    per_locus_gdna_df:    pd.DataFrame

    # ---- Convenience zero-copy aliases ----
    @property
    def alpha_gdna(self) -> np.ndarray:
        return self.prior_table.alpha_gdna

    @property
    def alpha_rna(self) -> np.ndarray:
        return self.prior_table.alpha_rna

    @property
    def prior_weight_rna(self) -> list[np.ndarray]:
        return self.prior_table.prior_weight_rna

    @property
    def global_fl_mean(self) -> float:
        return float(self.fl_models.global_.mean)

    @property
    def rna_fl_mean(self) -> float:
        return float(self.fl_models.rna.mean)

    @property
    def gdna_fl_mean(self) -> float:
        return float(self.fl_models.gdna.mean)

    # ---- Mutator-style helper (frozen-safe) ----
    def with_priors(self, prior_table: PriorTable) -> "CalibrationResult":
        return dataclasses.replace(
            self,
            prior_table=prior_table,
            n_multi_loci=len(prior_table.multi_locus_priors),
            multi_locus_prior_df=build_multi_locus_prior_df(
                prior_table.multi_locus_priors
            ),
            per_locus_gdna_df=build_per_locus_gdna_df(
                prior_table.multi_locus_priors
            ),
        )

    def to_summary_dict(self) -> dict[str, object]:
        mean_pi = (
            float(np.mean([m.pi_gdna for m in self.prior_table.multi_locus_priors]))
            if self.prior_table.multi_locus_priors
            else 0.0
        )
        return {
            "global_densities": self.global_densities.to_summary_dict(),
            "fl_models":        self.fl_models.to_summary_dict(),
            "diagnostics":      self.diagnostics.to_summary_dict(),
            "n_multi_loci":     self.n_multi_loci,
            "c_base":           float(self.prior_table.c_base_value),
            "mean_pi_gdna":     mean_pi,
        }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def build_calibration_result(
    *,
    payload:          CalibrationScanPayload,
    scan_trained:     FragmentLengthModels,
    global_densities: GlobalDensityTable,
    prior_table:      PriorTable | None = None,
    fl_prior_ess:     float = POOL_EB_PRIOR_ESS,
    fl_models:        FLModels | None = None,
) -> CalibrationResult:
    """Assemble the immutable v6 calibration result.

    Six explicit lines, zero ambiguity about ownership: the FL submodule
    owns the FL pipeline, the diagnostics submodule owns the named
    breakdown, the locus_prior submodule owns the priors, and this
    function just wires them together.

    ``prior_table`` defaults to :meth:`PriorTable.empty` so callers
    that have not yet built the locus graph (e.g. the calibration
    orchestrator) can produce a result and later swap the real table
    in via :meth:`CalibrationResult.with_priors`.

    ``fl_models`` may be passed by callers that already built the
    FL models (e.g. to seed ``compute_global_densities`` with the
    gDNA-FL mean).  When ``None`` (default), this function builds them
    via :func:`build_fl_models`.
    """
    if prior_table is None:
        prior_table = PriorTable.empty()
    if fl_models is None:
        fl_models = build_fl_models(
            global_counts = extract_global_counts(scan_trained),
            rna_counts    = extract_rna_counts(scan_trained),
            gdna_counts   = extract_gdna_counts(payload),
            max_size      = scan_trained.max_size,
            prior_ess     = fl_prior_ess,
        )
    diagnostics = Diagnostics.from_payload(payload)
    return CalibrationResult(
        global_densities=global_densities,
        fl_models=fl_models,
        prior_table=prior_table,
        diagnostics=diagnostics,
        n_multi_loci=len(prior_table.multi_locus_priors),
        multi_locus_prior_df=build_multi_locus_prior_df(prior_table.multi_locus_priors),
        per_locus_gdna_df=build_per_locus_gdna_df(prior_table.multi_locus_priors),
    )
