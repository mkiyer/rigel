"""rigel.calibration._result_v6 — Immutable v6 calibration result.

Per ``docs/calibration/calibration_v6_plan.md`` §2.10 and
``docs/calibration/m7_implementation_plan.md`` §4.

Coexists with the legacy SRD-v1 :class:`rigel.calibration._result.CalibrationResult`
until M8c, when the legacy schema is deleted and this module is renamed
to ``result.py``.

Field discipline:

* No ``version`` field (the class is the version).
* No ``active`` / ``enabled`` flag (existence implies use).
* No mutable cache slots (``with_priors`` returns a new instance via
  ``dataclasses.replace``).
* No reference to the raw ``CalibrationScanPayload`` (downstream code
  reads diagnostic dataframes, not ``fl_hist``).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .density_global import GlobalDensityTable
from ._fl_pool import PoolFLModels
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

#: Column order for ``CalibrationResult.multi_locus_prior_df``.
_MULTI_LOCUS_COLUMNS: tuple[str, ...] = (
    "multi_locus_id", "n_obs", "n_gdna", "n_rna", "pi_gdna", "n_loci",
)

#: Column order for ``CalibrationResult.per_locus_gdna_df``.
_PER_LOCUS_COLUMNS: tuple[str, ...] = (
    "multi_locus_id", "ref", "start", "end", "span",
    "n_obs", "n_gdna",
    "n_gdna_intergenic", "n_gdna_intron", "n_gdna_exon_intron",
    "pi_gdna", "n_eligible_boundaries", "fallback_flags",
)


def build_multi_locus_prior_df(
    mlps: tuple[MultiLocusPrior, ...],
) -> pd.DataFrame:
    """One row per ``MultiLocusPrior`` with the locked schema."""
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
    """One row per ``LocusGdnaEstimate`` (across all MultiLoci)."""
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

    # ---- Block 1: global gDNA densities (M4) ----
    global_densities: GlobalDensityTable

    # ---- Block 2: pool FL models (M7) ----
    pool: PoolFLModels

    # ---- Block 3: per-MultiLocus priors (M6) ----
    prior_table: PriorTable

    # ---- Block 4: provenance ----
    payload_summary: dict[str, int]
    n_multi_loci: int

    # ---- Block 5: derived diagnostic dataframes (eager) ----
    multi_locus_prior_df: pd.DataFrame
    per_locus_gdna_df: pd.DataFrame

    # ---- Convenience accessors (zero-copy) ----
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
    def gdna_fl_mean(self) -> float:
        return self.pool.gdna_fl_mean

    # ---- Mutator-style helpers (return a new instance; frozen-safe) ----
    def with_priors(self, prior_table: PriorTable) -> "CalibrationResult":
        """Return a copy with ``prior_table`` replaced.

        Both diagnostic dataframes are rebuilt from the new
        ``prior_table.multi_locus_priors``.
        """
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
        """JSON-serializable summary for ``summary.json``."""
        mean_pi = (
            float(np.mean([m.pi_gdna for m in self.prior_table.multi_locus_priors]))
            if self.prior_table.multi_locus_priors
            else 0.0
        )
        return {
            "global_densities": self.global_densities.to_summary_dict(),
            "pool":             self.pool.to_summary_dict(),
            "n_multi_loci":     self.n_multi_loci,
            "payload_summary":  dict(self.payload_summary),
            "c_base":           float(self.prior_table.c_base_value),
            "mean_pi_gdna":     mean_pi,
        }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def build_calibration_result(
    *,
    payload: CalibrationScanPayload,
    global_densities: GlobalDensityTable,
    pool: PoolFLModels,
    prior_table: PriorTable,
) -> CalibrationResult:
    """Assemble the immutable v6 calibration result."""
    payload_summary = {
        "n_observed":           int(payload.n_observed),
        "n_excluded_multimap":  int(payload.n_excluded_multimap),
        "n_excluded_chimera":   int(payload.n_excluded_chimera),
        "n_excluded_artifact":  int(payload.n_excluded_artifact),
        "n_unobserved":         int(payload.n_unobserved),
        "n_unannotated_ref":    int(payload.n_unannotated_ref),
    }

    return CalibrationResult(
        global_densities=global_densities,
        pool=pool,
        prior_table=prior_table,
        payload_summary=payload_summary,
        n_multi_loci=len(prior_table.multi_locus_priors),
        multi_locus_prior_df=build_multi_locus_prior_df(prior_table.multi_locus_priors),
        per_locus_gdna_df=build_per_locus_gdna_df(prior_table.multi_locus_priors),
    )
