"""Per-locus and per-MultiLocus priors for the EM solver.

This module is the assembly point of the v6 calibration prior pipeline.
It consumes:

* :class:`GlobalDensityTable` from :mod:`rigel.calibration.density_global`
  (M4),
* :class:`CalibrationScanPayload` from :mod:`rigel.calibration.scan_payload`
  (M3),
* the :class:`MultiLocus` / :class:`Locus` graph from :mod:`rigel.locus`
  (M5),

and produces:

* :class:`LocusGdnaEstimate` — per-``Locus`` gDNA mass estimate.
* :class:`MultiLocusPrior` — aggregation across the constituent ``Locus``
  intervals of a ``MultiLocus``.
* :class:`PriorTable` — the full table consumed by the batch EM
  (``alpha_gdna``, ``alpha_rna``, ``prior_weight_rna``).

The M5 helper :func:`build_prior_weight_rna` is preserved unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..locus import Locus, MultiLocus
from ..scored_fragments import ScoredFragments
from ._arrays import PayloadArrays, RegionArrays
from ._locus_n_obs import build_t_to_local_locus, partition_units_to_loci
from ._region_index_py import RegionIndexPy
from .density_global import GlobalDensityTable
from .density_loco import shrink_to_loco
from .regions import RegionType
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:
    from ..index import TranscriptIndex


__all__ = [
    "C_BASE_DEFAULT",
    "FLAG_INTERGENIC_ZERO_LEFF",
    "FLAG_INTRON_ZERO_LEFF",
    "FLAG_EXON_INTRON_NO_ELIGIBLE",
    "FLAG_PI_CLIPPED",
    "LocusGdnaEstimate",
    "MultiLocusPrior",
    "PriorTable",
    "estimate_locus_gdna",
    "assemble_multilocus_prior",
    "assemble_priors",
    "build_prior_weight_rna",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default Dirichlet evidence strength for the per-MultiLocus
#: ``(α_gdna, α_rna)`` prior.  Matches the SRD-v1 default for
#: continuity; tuned via the synthetic siphon harness.
C_BASE_DEFAULT = 10.0

#: Bitmask flags recorded on :class:`LocusGdnaEstimate.fallback_flags`.
FLAG_INTERGENIC_ZERO_LEFF = 1 << 0
FLAG_INTRON_ZERO_LEFF = 1 << 1
FLAG_EXON_INTRON_NO_ELIGIBLE = 1 << 2
FLAG_PI_CLIPPED = 1 << 3


# ---------------------------------------------------------------------------
# Result schemas
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class LocusGdnaEstimate:
    """Per-``Locus`` gDNA mass estimate."""

    locus: Locus
    n_obs: int
    n_gdna_intergenic: float
    n_gdna_intron: float
    n_gdna_exon_intron: float
    n_gdna: float
    pi_gdna: float                                 # clipped to [0, 1]
    rho_loco: tuple[float, float, float]           # (intergenic, intron, exon-intron)
    leff_loco: tuple[float, float, float]
    n_eligible_boundaries: int
    fallback_flags: int                             # bitmask of FLAG_*


@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    """Per-``MultiLocus`` prior (aggregated across its ``Locus`` intervals)."""

    multi_locus_id: int
    n_obs: int
    n_gdna: float
    n_rna: float
    pi_gdna: float
    per_locus: tuple[LocusGdnaEstimate, ...]


@dataclass(frozen=True, slots=True)
class PriorTable:
    """The full prior table consumed by the batch EM.

    ``alpha_gdna[i]`` and ``alpha_rna[i]`` correspond to
    ``multi_locus_priors[i].multi_locus_id == i`` (priors are stored
    in the same order as the input ``multi_loci`` list, indexed by
    ``multi_locus_id``).
    """

    multi_locus_priors: tuple[MultiLocusPrior, ...]
    alpha_gdna: np.ndarray                          # float64, (n_loci,)
    alpha_rna: np.ndarray                           # float64, (n_loci,)
    prior_weight_rna: list[np.ndarray]              # float32, [n_loci][n_components_i]
    c_base_value: float

    @classmethod
    def empty(cls, *, c_base_value: float = C_BASE_DEFAULT) -> "PriorTable":
        """Return a degenerate prior table with no loci.

        Used as the lazy-backfill seed in the calibration orchestrator:
        :func:`rigel.calibration.calibrate` returns a
        :class:`CalibrationResult` with this empty table, and
        :meth:`CalibrationResult.with_priors` swaps in the real table
        once the locus graph has been built.
        """
        return cls(
            multi_locus_priors=(),
            alpha_gdna=np.empty(0, dtype=np.float64),
            alpha_rna=np.empty(0, dtype=np.float64),
            prior_weight_rna=[],
            c_base_value=float(c_base_value),
        )


# ---------------------------------------------------------------------------
# M5 helper: per-component nRNA-suppression weights (preserved)
# ---------------------------------------------------------------------------

def build_prior_weight_rna(
    multi_locus: MultiLocus,
    em_data: ScoredFragments | None = None,  # noqa: ARG001 (reserved for future use)
    *,
    is_synthetic: np.ndarray | None = None,
    nrna_weight: float = 0.0,
) -> np.ndarray:
    """Construct the per-component nRNA-suppression weight vector.

    Components are laid out by the EM as ``[t_0, t_1, ..., t_{n_t-1},
    gDNA]``.  Synthetic nRNA spans live in the same ``index.t_df``
    transcript table as real transcripts (one component each), flagged
    by the ``is_synthetic`` boolean column.  This helper returns a
    ``float32`` array of length ``n_t + 1`` with real-mRNA entries set
    to ``1.0``, synthetic-nRNA entries set to ``nrna_weight``, and the
    trailing gDNA entry set to ``1.0`` (the gDNA component is not
    affected by ``prior_weight_rna``; the solver routes it through
    ``alpha_gdna`` / ``gdna_idx``).

    Parameters
    ----------
    multi_locus : MultiLocus
    em_data : ScoredFragments, optional
        Reserved for future per-fragment weighting.  Unused.
    is_synthetic : np.ndarray, optional
        Global ``index.t_df["is_synthetic"]`` boolean column.  When
        ``None``, every component is treated as mRNA (the legacy M5
        all-ones behavior — preserves bit-identical output for
        callers that have not yet plumbed ``nrna_weight``).
    nrna_weight : float
        Per-component weight applied to synthetic-nRNA entries.
        ``0.0`` (default) zeros out the nRNA prior contribution;
        ``1.0`` puts nRNA on equal footing with mRNA.  Ignored when
        ``is_synthetic`` is ``None``.
    """
    n_t = int(multi_locus.transcript_indices.shape[0])
    weights = np.ones(n_t + 1, dtype=np.float32)
    if is_synthetic is None:
        return weights
    t_idx = multi_locus.transcript_indices
    synth_mask = np.asarray(is_synthetic, dtype=bool)[t_idx]
    weights[:n_t] = np.where(synth_mask, np.float32(nrna_weight), np.float32(1.0))
    return weights


# ---------------------------------------------------------------------------
# Per-Locus core
# ---------------------------------------------------------------------------

def _shrink_one_type(
    type_mask: np.ndarray,
    region_ids: np.ndarray,
    count_col: np.ndarray,
    leff: np.ndarray,
    rho_global: float,
    kappa: float,
) -> tuple[float, float, float]:
    """INTERGENIC / INTRON branch: one count column, one denominator.

    Returns ``(n_loco, leff_loco, rho_loco)``.
    """
    if not type_mask.any():
        return 0.0, 0.0, rho_global
    rids = region_ids[type_mask]
    n_loco = float(count_col[rids].sum())
    leff_loco = float(leff[type_mask].sum())
    rho_loco = shrink_to_loco(n_loco, leff_loco, rho_global, kappa)
    return n_loco, leff_loco, rho_loco


def _shrink_exon_intron(
    types: np.ndarray,
    region_ids: np.ndarray,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    leff: np.ndarray,
    rho_global: float,
    kappa: float,
    gdna_fl_mean: float,
) -> tuple[float, float, int, float]:
    """EXON-INTRON branch: capture-aware boundary-flux density × full
    exonic L_eff in this locus.

    Returns ``(n_loco_eligible, leff_full_or_zero, n_eligible_boundaries, rho_loco)``.
    The L_eff of the predicted-count step is the full exonic L_eff inside
    the locus, NOT the eligible-boundary slice (per plan §2.7.1 step 3).
    """
    m = types == int(RegionType.EXON)
    if not m.any():
        return 0.0, 0.0, 0, rho_global
    rids = region_ids[m]
    bf_l = region_arrays.bf_left[rids]
    bf_r = region_arrays.bf_right[rids]
    n_eligible = int(bf_l.sum() + bf_r.sum())
    leff_full = float(leff[m].sum())

    if n_eligible == 0:
        return 0.0, leff_full, 0, rho_global

    n_loco_eligible = float(
        (payload_arrays.u_left[rids] * bf_l.astype(np.int64)).sum()
        + (payload_arrays.u_right[rids] * bf_r.astype(np.int64)).sum()
    )
    leff_loco_eligible = float(n_eligible) * gdna_fl_mean
    rho_loco = shrink_to_loco(n_loco_eligible, leff_loco_eligible, rho_global, kappa)
    return n_loco_eligible, leff_full, n_eligible, rho_loco


def estimate_locus_gdna(
    locus: Locus,
    n_obs: int,
    region_index: RegionIndexPy,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    global_densities: GlobalDensityTable,
    gdna_fl_mean: float,
) -> LocusGdnaEstimate:
    """Estimate the gDNA mass and π̂_gdna for one ``Locus``.

    Per ``docs/calibration/calibration_v6_plan.md`` §2.7.1.
    """
    region_ids = region_index.overlap(locus.ref_id, locus.start, locus.end)
    if region_ids.size == 0:
        raise RuntimeError(
            f"estimate_locus_gdna: Locus(ref={locus.ref!r}, start={locus.start}, "
            f"end={locus.end}) overlaps no regions. BAM reference does not "
            f"match index — rebuild the index or check the BAM header."
        )

    types = region_arrays.type[region_ids]
    starts = region_arrays.start[region_ids]
    ends = region_arrays.end[region_ids]
    cl_lo = np.maximum(starts, locus.start)
    cl_hi = np.minimum(ends, locus.end)
    cl_len = (cl_hi - cl_lo).astype(np.float64)
    leff = cl_len + (gdna_fl_mean - 1.0)            # l_eff_overlap, vectorized

    fallback_flags = 0

    # INTERGENIC
    n_ig, leff_ig, rho_loco_ig = _shrink_one_type(
        types == int(RegionType.INTERGENIC),
        region_ids,
        payload_arrays.intergenic_per_region,
        leff,
        global_densities.intergenic.rho,
        global_densities.intergenic.kappa.value,
    )
    n_gdna_ig = rho_loco_ig * leff_ig
    if leff_ig <= 0.0:
        fallback_flags |= FLAG_INTERGENIC_ZERO_LEFF

    # INTRON
    n_in, leff_in, rho_loco_in = _shrink_one_type(
        types == int(RegionType.INTRON),
        region_ids,
        payload_arrays.intron_per_region,
        leff,
        global_densities.intron.rho,
        global_densities.intron.kappa.value,
    )
    n_gdna_in = rho_loco_in * leff_in
    if leff_in <= 0.0:
        fallback_flags |= FLAG_INTRON_ZERO_LEFF

    # EXON-INTRON
    n_ei, leff_ei, n_eligible_boundaries, rho_loco_ei = _shrink_exon_intron(
        types,
        region_ids,
        region_arrays,
        payload_arrays,
        leff,
        global_densities.exon_intron.rho,
        global_densities.exon_intron.kappa.value,
        gdna_fl_mean,
    )
    n_gdna_ei = rho_loco_ei * leff_ei
    if n_eligible_boundaries == 0:
        fallback_flags |= FLAG_EXON_INTRON_NO_ELIGIBLE

    n_gdna_total = n_gdna_ig + n_gdna_in + n_gdna_ei
    pi_unclipped = (n_gdna_total / n_obs) if n_obs > 0 else 0.0
    pi_gdna = float(min(1.0, max(0.0, pi_unclipped)))
    if pi_unclipped > 1.0:
        fallback_flags |= FLAG_PI_CLIPPED

    return LocusGdnaEstimate(
        locus=locus,
        n_obs=int(n_obs),
        n_gdna_intergenic=float(n_gdna_ig),
        n_gdna_intron=float(n_gdna_in),
        n_gdna_exon_intron=float(n_gdna_ei),
        n_gdna=float(n_gdna_total),
        pi_gdna=pi_gdna,
        rho_loco=(float(rho_loco_ig), float(rho_loco_in), float(rho_loco_ei)),
        leff_loco=(float(leff_ig), float(leff_in), float(leff_ei)),
        n_eligible_boundaries=int(n_eligible_boundaries),
        fallback_flags=int(fallback_flags),
    )


# ---------------------------------------------------------------------------
# MultiLocus aggregation + orchestrator
# ---------------------------------------------------------------------------

def assemble_multilocus_prior(
    ml: MultiLocus,
    per_locus_estimates: tuple[LocusGdnaEstimate, ...],
) -> MultiLocusPrior:
    """Aggregate per-``Locus`` estimates into one ``MultiLocusPrior``."""
    n_obs = sum(e.n_obs for e in per_locus_estimates)
    n_gdna = sum(e.n_gdna for e in per_locus_estimates)
    n_rna = max(0.0, n_obs - n_gdna)
    pi_raw = (n_gdna / n_obs) if n_obs > 0 else 0.0
    return MultiLocusPrior(
        multi_locus_id=int(ml.multi_locus_id),
        n_obs=int(n_obs),
        n_gdna=float(n_gdna),
        n_rna=float(n_rna),
        pi_gdna=float(min(1.0, max(0.0, pi_raw))),
        per_locus=per_locus_estimates,
    )


def assemble_priors(
    multi_loci: list[MultiLocus],
    em_data: ScoredFragments,
    index: "TranscriptIndex",
    payload: CalibrationScanPayload,
    global_densities: GlobalDensityTable,
    *,
    gdna_fl_mean: float | None = None,
    c_base: float = C_BASE_DEFAULT,
    nrna_weight: float = 0.0,
) -> PriorTable:
    """Build the full :class:`PriorTable` for the batch EM.

    Pure-Python pass; no I/O, no mutation of inputs.

    Parameters
    ----------
    nrna_weight : float
        Per-component weight applied to synthetic-nRNA transcripts in
        each ``MultiLocus``.  See :func:`build_prior_weight_rna`.
    """
    if gdna_fl_mean is None:
        gdna_fl_mean = global_densities.gdna_fl_mean
    if not (gdna_fl_mean > 0.0):
        raise ValueError(
            f"assemble_priors: gdna_fl_mean must be > 0; got {gdna_fl_mean!r}."
        )

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    region_index = RegionIndexPy(arrays=region_arrays)

    # Hoist per-transcript columns once (slow path only consults these).
    t_starts = index.t_df["start"].values
    t_ref_codes = index.t_df["ref"].cat.codes.values
    cat_to_ref_id = np.array(
        [index.ref_name_to_id[str(name)] for name in index.t_df["ref"].cat.categories],
        dtype=np.int32,
    )
    t_ref = cat_to_ref_id[t_ref_codes]

    # ``is_synthetic`` flags rigel-generated nRNA spans (separate
    # transcripts in ``index.t_df``, but each its own EM component).
    # Older indices predate this column; treat absence as "no synthetic
    # transcripts" so legacy callers still get all-ones weights.
    if "is_synthetic" in index.t_df.columns:
        is_synthetic_arr = index.t_df["is_synthetic"].to_numpy(dtype=bool)
    else:
        is_synthetic_arr = None
    n_ml = len(multi_loci)
    alpha_gdna = np.zeros(n_ml, dtype=np.float64)
    alpha_rna = np.zeros(n_ml, dtype=np.float64)
    multi_locus_priors: list[MultiLocusPrior | None] = [None] * n_ml
    prior_weight_rna: list[np.ndarray | None] = [None] * n_ml

    for ml in multi_loci:
        t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
        units_per_locus = partition_units_to_loci(ml, em_data, t_to_local)

        per_locus_est = tuple(
            estimate_locus_gdna(
                locus=loc,
                n_obs=int(units_per_locus[j].size),
                region_index=region_index,
                region_arrays=region_arrays,
                payload_arrays=payload_arrays,
                global_densities=global_densities,
                gdna_fl_mean=gdna_fl_mean,
            )
            for j, loc in enumerate(ml.loci)
        )

        ml_prior = assemble_multilocus_prior(ml, per_locus_est)
        idx = int(ml.multi_locus_id)
        if multi_locus_priors[idx] is not None:
            raise RuntimeError(
                f"assemble_priors: duplicate multi_locus_id={idx} in multi_loci."
            )
        multi_locus_priors[idx] = ml_prior
        alpha_gdna[idx] = c_base * ml_prior.pi_gdna
        alpha_rna[idx] = c_base * (1.0 - ml_prior.pi_gdna)
        prior_weight_rna[idx] = build_prior_weight_rna(
            ml,
            em_data,
            is_synthetic=is_synthetic_arr,
            nrna_weight=nrna_weight,
        )

    if any(p is None for p in multi_locus_priors):
        missing = [i for i, p in enumerate(multi_locus_priors) if p is None]
        raise RuntimeError(
            f"assemble_priors: multi_loci list does not cover multi_locus_ids "
            f"{missing}. Caller must pass exactly one MultiLocus per id 0..n-1."
        )

    return PriorTable(
        multi_locus_priors=tuple(multi_locus_priors),  # type: ignore[arg-type]
        alpha_gdna=alpha_gdna,
        alpha_rna=alpha_rna,
        prior_weight_rna=prior_weight_rna,  # type: ignore[arg-type]
        c_base_value=float(c_base),
    )
