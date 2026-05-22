"""rigel.calibration.locus_prior \u2014 fail-fast stubs (fractional cutover).

The integer 8-mask per-locus prior assembler is removed during the
fractional cutover. Dataclasses (``LocusGdnaEstimate``,
``MultiLocusPrior``, ``PriorTable``, ``ExpectedGdnaPriorParts``,
``ExonGdnaDiagnostics``) and their fallback bitflags are preserved so
result-serialization code and tests compile. Entry-point functions
raise :class:`FractionalCutoverPending`.

``PriorTable.empty()`` is retained because the orchestrator uses it as
the lazy-backfill seed in :class:`CalibrationResult`.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

from .errors import FractionalCutoverPending


__all__ = [
    "INTERGENIC_FLANK_BP_DEFAULT",
    "FLAG_INTERGENIC_ZERO_LEFF",
    "FLAG_INTRON_ZERO_LEFF",
    "FLAG_EXON_INTRON_NO_ELIGIBLE",
    "FLAG_PI_CLIPPED",
    "ExpectedGdnaPriorParts",
    "ExonGdnaDiagnostics",
    "LocusGdnaEstimate",
    "MultiLocusPrior",
    "PriorTable",
    "estimate_locus_gdna",
    "expected_gdna_count_global",
    "enable_gdna_for_multilocus",
    "assemble_multilocus_prior",
    "assemble_priors",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INTERGENIC_FLANK_BP_DEFAULT = 5_000

FLAG_INTERGENIC_ZERO_LEFF = 1 << 0
FLAG_INTRON_ZERO_LEFF = 1 << 1
FLAG_EXON_INTRON_NO_ELIGIBLE = 1 << 2
FLAG_PI_CLIPPED = 1 << 3


# ---------------------------------------------------------------------------
# Result schemas
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ExonGdnaDiagnostics:
    """EXON branch diagnostics for locoregional gDNA estimates."""

    rho_exon_boundary: float = 0.0
    rho_exon_contained: float = 0.0
    rho_exon_composite: float = 0.0
    exon_boundary_precision: float = 0.0
    exon_contained_precision: float = 0.0
    n_exon_contained_observed: float = 0.0
    n_exon_contained_estimated: float = 0.0


@dataclass(frozen=True, slots=True)
class LocusGdnaEstimate:
    """Per-``Locus`` gDNA mass estimate (shape preserved across cutover)."""

    # ``locus`` is typed as object to avoid a load-time import of
    # ``rigel.locus`` from the calibration package.
    locus: object
    n_obs: int
    n_gdna_intergenic: float
    n_gdna_intron: float
    n_gdna_boundary_observed: float
    n_gdna_exon_only: float
    n_gdna: float
    pi_gdna: float
    rho_loco: tuple[float, float, float]
    leff_loco: tuple[float, float, float]
    n_eligible_boundaries: int
    n_boundary_events: float
    nrna_active: bool
    fallback_flags: int
    exon_gdna: ExonGdnaDiagnostics = dataclasses.field(default_factory=ExonGdnaDiagnostics)


@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    """Per-``MultiLocus`` prior aggregation."""

    multi_locus_id: int
    n_obs: int
    n_gdna: float
    n_rna: float
    pi_gdna: float
    gdna_prior_count: float
    per_locus: tuple[LocusGdnaEstimate, ...]


@dataclass(frozen=True, slots=True)
class PriorTable:
    """The full prior table consumed by the batch EM."""

    multi_locus_priors: tuple[MultiLocusPrior, ...]
    gdna_prior_count: np.ndarray
    gdna_prior_count_em: np.ndarray
    gdna_eff_len: np.ndarray
    enable_gdna: np.ndarray
    gdna_eff_len_unweighted: np.ndarray = dataclasses.field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    gdna_em_exposure_weight: np.ndarray = dataclasses.field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )

    @classmethod
    def empty(cls) -> "PriorTable":
        return cls(
            multi_locus_priors=(),
            gdna_prior_count=np.empty(0, dtype=np.float64),
            gdna_prior_count_em=np.empty(0, dtype=np.float64),
            gdna_eff_len=np.empty(0, dtype=np.float64),
            enable_gdna=np.empty(0, dtype=np.uint8),
            gdna_eff_len_unweighted=np.empty(0, dtype=np.float64),
            gdna_em_exposure_weight=np.empty(0, dtype=np.float64),
        )


@dataclass(frozen=True, slots=True)
class ExpectedGdnaPriorParts:
    """Decomposition of the global-only expected gDNA pseudocount."""

    total: float
    intergenic_contained: float
    intron_contained: float
    boundary_crossing_expected: float
    exon_contained_expected: float
    density_exposure_intergenic: float
    density_exposure_intron: float
    density_exposure_boundary: float


# ---------------------------------------------------------------------------
# Fail-fast entry points
# ---------------------------------------------------------------------------


def estimate_locus_gdna(*_args: object, **_kwargs: object) -> LocusGdnaEstimate:
    raise FractionalCutoverPending(
        "estimate_locus_gdna: per-locus gDNA estimation consumes the "
        "integer 8-mask payload and is removed during the fractional cutover."
    )


def expected_gdna_count_global(*_args: object, **_kwargs: object) -> ExpectedGdnaPriorParts:
    raise FractionalCutoverPending(
        "expected_gdna_count_global: global-only gDNA pseudocount consumes "
        "the legacy GlobalDensityTable and is removed during the fractional "
        "cutover."
    )


def enable_gdna_for_multilocus(*_args: object, **_kwargs: object) -> bool:
    raise FractionalCutoverPending(
        "enable_gdna_for_multilocus: gDNA eligibility flag construction is "
        "gated on the fractional prior pipeline."
    )


def assemble_multilocus_prior(*_args: object, **_kwargs: object) -> MultiLocusPrior:
    raise FractionalCutoverPending(
        "assemble_multilocus_prior: integer prior aggregation is removed "
        "during the fractional cutover."
    )


def assemble_priors(*_args: object, **_kwargs: object) -> PriorTable:
    raise FractionalCutoverPending(
        "assemble_priors: the full integer prior pipeline is removed during "
        "the fractional cutover. The replacement fractional prior pipeline "
        "has not yet landed."
    )
