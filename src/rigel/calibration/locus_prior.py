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

from ..frag_length_model import FragmentLengthModel
from ..locus import Locus, MultiLocus
from ..scored_fragments import ScoredFragments
from ._arrays import PayloadArrays, RegionArrays
from ._exposure import (
    boundary_crossing_exposure,
    boundary_side_in_window,
    contained_exposure_clipped,
)
from ._locus_n_obs import build_t_to_local_locus, partition_units_to_loci
from ._region_index_py import RegionIndexPy
from .density_global import GlobalDensityTable
from .density_loco import shrink_to_loco
from .regions import RegionType
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:
    from ..index import TranscriptIndex


__all__ = [
    "INTERGENIC_FLANK_BP_DEFAULT",
    "FLAG_INTERGENIC_ZERO_LEFF",
    "FLAG_INTRON_ZERO_LEFF",
    "FLAG_EXON_INTRON_NO_ELIGIBLE",
    "FLAG_PI_CLIPPED",
    "ExpectedGdnaPriorParts",
    "LocusGdnaEstimate",
    "MultiLocusPrior",
    "PriorTable",
    "estimate_locus_gdna",
    "expected_gdna_count_global",
    "enable_gdna_for_multilocus",
    "assemble_multilocus_prior",
    "assemble_priors",
    "build_prior_weight_rna",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Phase 4 default: bp added to both sides of the intergenic-density
#: query window. Mass terms still use the unflanked locus interval; the
#: flank only widens the *evidence* window for the intergenic branch.
#: Phase 4 ships this at 5_000 (Phase 3 shipped 0).
INTERGENIC_FLANK_BP_DEFAULT = 5_000

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
    """Per-``Locus`` gDNA mass estimate.

    The ``n_gdna_*`` mass fields exactly partition ``n_gdna``:

    * ``n_gdna_intergenic`` — local intergenic density × core-locus
      intergenic exposure (typically 0 for transcript-defined loci).
    * ``n_gdna_intron`` — local intron density × core-locus intron exposure.
    * ``n_gdna_boundary_observed`` — *observed* boundary-crossing fragment
      count restricted to eligible sides inside the locus.
    * ``n_gdna_exon_only`` — local boundary density × core-locus
      exon-contained exposure, the imputed exon-contained gDNA mass that
      does not enter the C++ scanner's boundary-flux numerator.

    All four describe the same EM-unit population as ``n_obs``, so
    ``pi_gdna = n_gdna / n_obs`` is arithmetically consistent with the
    fragments that downstream EM will route.
    """

    locus: Locus
    n_obs: int
    n_gdna_intergenic: float
    n_gdna_intron: float
    n_gdna_boundary_observed: float
    n_gdna_exon_only: float
    n_gdna: float
    pi_gdna: float                                 # clipped to [0, 1]
    rho_loco: tuple[float, float, float]           # (intergenic, intron, boundary)
    leff_loco: tuple[float, float, float]          # (L_core_ig, L_core_in, L_core_exon)
    n_eligible_boundaries: int                      # eligible sides inside locus
    n_boundary_events: float                        # raw observed boundary count in locus
    nrna_active: bool                               # set in Phase 6; False here
    fallback_flags: int                             # bitmask of FLAG_*

    @property
    def n_gdna_exon_intron(self) -> float:
        """Compatibility aggregate kept for one release.

        Equals ``n_gdna_boundary_observed + n_gdna_exon_only``. The two
        split fields have replaced this aggregate as the canonical
        decomposition; see
        ``docs/calibration/locoregional_gdna_redesign_2026-05-07.md``.
        """
        return self.n_gdna_boundary_observed + self.n_gdna_exon_only


@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    """Per-``MultiLocus`` prior (aggregated across its ``Locus`` intervals).

    Bayesian-prior redesign Phase 1 schema:

    * ``gdna_prior_count`` (= :math:`\\eta_g`) is the canonical
      asymmetric pseudocount the EM consumes through ``alpha_gdna``.
    * ``rna_prior_count`` (= :math:`\\eta_r`) is the symmetric companion
      and is wired to ``0`` once Phase 2 lands.
    * ``n_obs``/``n_gdna``/``n_rna``/``pi_gdna`` remain as diagnostic
      summaries of the locoregional pipeline; they do not affect the EM
      prior beyond this transitional phase.
    """

    multi_locus_id: int
    n_obs: int
    n_gdna: float
    n_rna: float
    pi_gdna: float
    gdna_prior_count: float
    rna_prior_count: float
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
    #: Bayesian-prior redesign Phase 1 canonical fields.
    gdna_prior_count: np.ndarray                    # float64, (n_loci,)
    rna_prior_count: np.ndarray                     # float64, (n_loci,)
    #: ``uint8`` flag forwarded to the native EM as ``locus_enable_gdna``.
    enable_gdna: np.ndarray                         # uint8, (n_loci,)

    @classmethod
    def empty(cls) -> "PriorTable":
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
            gdna_prior_count=np.empty(0, dtype=np.float64),
            rna_prior_count=np.empty(0, dtype=np.float64),
            enable_gdna=np.empty(0, dtype=np.uint8),
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
    """Construct the per-component RNA-prior allocation weight vector.

    .. note:: **Policy lock (Bayesian-prior redesign Phase 0.5).**

       ``prior_weight_rna`` is an *RNA-prior allocation weight*, not a
       likelihood term. The native EM uses it only inside
       ``compute_ovr_prior_and_warm_start`` to distribute a positive
       ``alpha_rna`` budget across RNA components proportional to
       coverage × weight. When ``alpha_rna == 0`` (the asymmetric-prior
       default in the redesign), this allocation path collapses and
       ``prior_weight_rna`` becomes objective-neutral by construction.

       Do not migrate this weight into the E-step likelihood as a
       backdoor way to keep nRNA suppression alive under
       ``alpha_rna == 0`` — that would replace one heuristic prior
       force with a heuristic likelihood force. nRNA-versus-gDNA
       arbitration is a separate modeling problem and will be
       addressed once the gDNA-versus-total-RNA prior foundation is
       stable. See ``docs/bayesian_prior/bayesian_prior_plan_v3.md``
       §5 Phase 0.5 and the corresponding test in
       ``tests/test_prior_weight_rna_policy.py``.

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

def _density_term_prorated(
    *,
    type_in_locus: np.ndarray,        # (R',) bool: regions in locus of this type
    region_ids_locus: np.ndarray,     # (R',) sorted-position region IDs in locus
    full_count_col: np.ndarray,       # (R,) per-region count (whole-region, payload)
    eff_full: np.ndarray,             # (R',) full-region FL-weighted exposure
    eff_clip_evidence: np.ndarray,    # (R',) evidence-window clipped exposure
    eff_clip_core: np.ndarray,        # (R',) core-locus clipped exposure
    rho_global: float,
    kappa: float,
) -> tuple[float, float, float, float, float]:
    """Prorated INTERGENIC / INTRON locus contribution.

    Returns ``(rho_loco, l_core, l_evidence, n_evidence, n_mass)`` where:

    * ``n_evidence`` and ``l_evidence`` parameterise the locoregional
      shrinkage — both use the (possibly flank-padded) evidence window.
    * ``n_mass = rho_loco * l_core`` is the predicted mass that lives
      *inside* the EM-unit window (``[locus.start, locus.end]``), keeping
      ``π_gdna`` arithmetic consistent with ``n_obs``.

    Per-region count proration: the C++ scanner emits one whole-region
    contained count per region; multiplying by
    ``ratio_r = eff_clip_r / eff_full_r`` gives the FL-aware fraction of
    that count attributable to the clipped window without double-counting
    fragments whose start position lies outside the window.
    """
    if not type_in_locus.any():
        return rho_global, 0.0, 0.0, 0.0, 0.0
    rids = region_ids_locus[type_in_locus]
    full_counts = full_count_col[rids].astype(np.float64)
    eff_full_t = eff_full[type_in_locus]
    eff_evidence_t = eff_clip_evidence[type_in_locus]
    eff_core_t = eff_clip_core[type_in_locus]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_evidence = np.where(
            eff_full_t > 0.0, eff_evidence_t / eff_full_t, 0.0
        )
    n_evidence = float((full_counts * ratio_evidence).sum())
    l_evidence = float(eff_evidence_t.sum())
    l_core = float(eff_core_t.sum())
    rho_loco = shrink_to_loco(n_evidence, l_evidence, rho_global, kappa)
    n_mass = rho_loco * l_core
    return rho_loco, l_core, l_evidence, n_evidence, n_mass


def _boundary_term_prorated(
    *,
    exon_in_locus: np.ndarray,       # (R',) bool: type == EXON for in-locus regions
    region_ids_locus: np.ndarray,    # (R',) sorted-position region IDs in locus
    starts_locus: np.ndarray,        # (R',) region starts (sorted-position)
    ends_locus: np.ndarray,          # (R',) region ends
    bf_left_locus: np.ndarray,       # (R',) eligibility flag (any region type)
    bf_right_locus: np.ndarray,      # (R',) eligibility flag
    u_left: np.ndarray,              # (R,) per-region left-edge flux (payload)
    u_right: np.ndarray,             # (R,) per-region right-edge flux
    eff_clip_core: np.ndarray,       # (R',) core-locus clipped exposure
    locus_start: int,
    locus_end: int,
    b_cross: float,
    rho_global: float,
    kappa: float,
) -> tuple[float, int, float, float, float]:
    """Boundary branch (Issues 1+2+3).

    Eligible boundary sides are those that (a) are flagged
    ``boundary_flux_*`` on an EXON region (capture-tile eligibility from
    the index) and (b) lie inside the locus interval. The boundary
    density is

        rho_boundary = N_crossing_in_locus
                       / (n_eligible_sides_in_locus * B_cross)

    Returns ``(rho_loco, n_eligible_in_locus, n_boundary_observed,
    l_core_exon_contained, n_exon_only_mass)``. The two predicted-mass
    fields are kept disjoint: ``n_boundary_observed`` is a raw fragment
    count and ``n_exon_only_mass = rho * l_core_exon`` predicts the
    *contained* exonic gDNA mass that the boundary-flux numerator does
    not see.
    """
    if not exon_in_locus.any():
        return rho_global, 0, 0.0, 0.0, 0.0

    rids_exon = region_ids_locus[exon_in_locus]
    starts_exon = starts_locus[exon_in_locus]
    ends_exon = ends_locus[exon_in_locus]
    bf_l_exon = bf_left_locus[exon_in_locus]
    bf_r_exon = bf_right_locus[exon_in_locus]

    # Restrict to boundary sides that are eligible *and* inside the locus.
    in_left, in_right = boundary_side_in_window(
        starts_exon, ends_exon, locus_start, locus_end
    )
    elig_left = bf_l_exon & in_left
    elig_right = bf_r_exon & in_right

    n_eligible = int(elig_left.sum() + elig_right.sum())
    l_core_exon = float(eff_clip_core[exon_in_locus].sum())

    if n_eligible == 0:
        # No local boundary evidence: shrink to global density and
        # impute exon-contained mass from it.
        n_exon_only = rho_global * l_core_exon
        return rho_global, 0, 0.0, l_core_exon, n_exon_only

    n_boundary_observed = float(
        (u_left[rids_exon] * elig_left.astype(np.int64)).sum()
        + (u_right[rids_exon] * elig_right.astype(np.int64)).sum()
    )
    denom = float(n_eligible) * b_cross
    rho_loco = shrink_to_loco(n_boundary_observed, denom, rho_global, kappa)
    n_exon_only = rho_loco * l_core_exon
    return rho_loco, n_eligible, n_boundary_observed, l_core_exon, n_exon_only


@dataclass(frozen=True, slots=True)
class _LocusScratch:
    """Per-Locus shared work between :func:`estimate_locus_gdna` and
    :func:`expected_gdna_count_global`.

    Computed once by :func:`_compute_locus_scratch` and threaded into
    both helpers so the per-Locus ``region_index.overlap`` query and
    FL-aware ``contained_exposure_clipped`` evaluations are not
    duplicated.  Internal optimisation used by :func:`assemble_priors`;
    the public helpers retain their original signatures and recompute
    the scratch on demand when ``_scratch`` is ``None``.
    """

    pad_lo: int
    pad_hi: int
    region_ids: np.ndarray
    types: np.ndarray
    starts: np.ndarray
    ends: np.ndarray
    bf_l: np.ndarray
    bf_r: np.ndarray
    eff_full: np.ndarray
    eff_clip_core: np.ndarray
    eff_clip_evidence: np.ndarray
    is_ig: np.ndarray
    is_in: np.ndarray
    is_ex: np.ndarray


def _compute_locus_scratch(
    locus: Locus,
    region_index: RegionIndexPy,
    region_arrays: RegionArrays,
    gdna_fl: FragmentLengthModel,
    *,
    intergenic_flank_bp: int,
    ref_length: int | None,
) -> _LocusScratch | None:
    """Build a :class:`_LocusScratch` for this locus.

    Returns ``None`` when the (padded) region overlap is empty so the
    caller can preserve the legacy "no overlap" behaviour
    (``estimate_locus_gdna`` raises; ``expected_gdna_count_global``
    returns a zeroed :class:`ExpectedGdnaPriorParts`).
    """
    locus_start = int(locus.start)
    locus_end = int(locus.end)
    pad = max(0, int(intergenic_flank_bp))
    pad_lo = max(0, locus_start - pad)
    pad_hi = locus_end + pad if ref_length is None else min(int(ref_length), locus_end + pad)

    region_ids = region_index.overlap(locus.ref_id, pad_lo, pad_hi)
    if region_ids.size == 0:
        return None

    types = region_arrays.type[region_ids]
    starts = region_arrays.start[region_ids]
    ends = region_arrays.end[region_ids]
    bf_l = region_arrays.bf_left[region_ids]
    bf_r = region_arrays.bf_right[region_ids]

    eff_full, eff_clip_core = contained_exposure_clipped(
        starts, ends, locus_start, locus_end, gdna_fl
    )
    if pad == 0:
        eff_clip_evidence = eff_clip_core
    else:
        _, eff_clip_evidence = contained_exposure_clipped(
            starts, ends, pad_lo, pad_hi, gdna_fl
        )

    is_ig = types == int(RegionType.INTERGENIC)
    is_in = types == int(RegionType.INTRON)
    is_ex = types == int(RegionType.EXON)

    return _LocusScratch(
        pad_lo=pad_lo,
        pad_hi=pad_hi,
        region_ids=region_ids,
        types=types,
        starts=starts,
        ends=ends,
        bf_l=bf_l,
        bf_r=bf_r,
        eff_full=eff_full,
        eff_clip_core=eff_clip_core,
        eff_clip_evidence=eff_clip_evidence,
        is_ig=is_ig,
        is_in=is_in,
        is_ex=is_ex,
    )


def estimate_locus_gdna(
    locus: Locus,
    n_obs: int,
    region_index: RegionIndexPy,
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    global_densities: GlobalDensityTable,
    gdna_fl: FragmentLengthModel,
    *,
    intergenic_flank_bp: int = INTERGENIC_FLANK_BP_DEFAULT,
    ref_length: int | None = None,
    splicing_anchor_tolerance: int = 0,
    b_cross: float | None = None,
    _scratch: _LocusScratch | None = None,
) -> LocusGdnaEstimate:
    """Estimate the gDNA mass and π̂_gdna for one ``Locus``.

    Implements §3.5 of
    ``docs/calibration/locoregional_gdna_redesign_2026-05-07.md``:

    1. Single padded query against ``region_index``.
    2. Two FL-aware exposure vectors per region: core (mass) and
       evidence (intergenic flank).
    3. Four-term decomposition (intergenic, intron,
       boundary-observed, exon-only-imputed) that exactly sums to
       ``n_gdna``.

    Parameters
    ----------
    intergenic_flank_bp : int
        Phase 4 padding added on both sides of the locus for the
        *intergenic* density's evidence window. Mass terms always use
        the unflanked locus interval. Defaults to 0 (no flank); raise
        in Phase 4.
    ref_length : int, optional
        Reference contig length used to clamp the padded window. When
        ``None`` (or when ``intergenic_flank_bp == 0``) no clamp is
        applied beyond the natural ``max(0, ...)`` lower bound.
    """
    locus_start = int(locus.start)
    locus_end = int(locus.end)

    # 1. Padded region query (§3.5 step 1). Mass terms always use the
    #    unflanked window; the pad widens only the *intergenic* evidence
    #    branch via ratio_evidence below.
    if _scratch is None:
        _scratch = _compute_locus_scratch(
            locus,
            region_index,
            region_arrays,
            gdna_fl,
            intergenic_flank_bp=intergenic_flank_bp,
            ref_length=ref_length,
        )
        if _scratch is None:
            raise RuntimeError(
                f"estimate_locus_gdna: Locus(ref={locus.ref!r}, start={locus_start}, "
                f"end={locus_end}) overlaps no regions. BAM reference does not "
                f"match index — rebuild the index or check the BAM header."
            )

    region_ids = _scratch.region_ids
    starts = _scratch.starts
    ends = _scratch.ends
    bf_l = _scratch.bf_l
    bf_r = _scratch.bf_r
    eff_full = _scratch.eff_full
    eff_clip_core = _scratch.eff_clip_core
    eff_clip_evidence = _scratch.eff_clip_evidence
    is_ig = _scratch.is_ig
    is_in = _scratch.is_in
    is_ex = _scratch.is_ex

    fallback_flags = 0

    # 3. INTERGENIC (§3.5 step 3).
    rho_ig, l_core_ig, _, n_ev_ig, n_mass_ig = _density_term_prorated(
        type_in_locus=is_ig,
        region_ids_locus=region_ids,
        full_count_col=payload_arrays.intergenic_per_region,
        eff_full=eff_full,
        eff_clip_evidence=eff_clip_evidence,
        eff_clip_core=eff_clip_core,
        rho_global=global_densities.intergenic.rho,
        kappa=global_densities.intergenic.kappa.value,
    )
    if l_core_ig <= 0.0:
        fallback_flags |= FLAG_INTERGENIC_ZERO_LEFF

    # 4. INTRON (§3.5 step 4); evidence == core (no flank for intron).
    rho_in, l_core_in, _, n_ev_in, n_mass_in = _density_term_prorated(
        type_in_locus=is_in,
        region_ids_locus=region_ids,
        full_count_col=payload_arrays.intron_per_region,
        eff_full=eff_full,
        eff_clip_evidence=eff_clip_core,    # no flank for intron
        eff_clip_core=eff_clip_core,
        rho_global=global_densities.intron.rho,
        kappa=global_densities.intron.kappa.value,
    )
    if l_core_in <= 0.0:
        fallback_flags |= FLAG_INTRON_ZERO_LEFF

    # 5. BOUNDARY (§3.5 step 5).
    if b_cross is None:
        b_cross = boundary_crossing_exposure(
            gdna_fl, splicing_anchor_tolerance=splicing_anchor_tolerance
        )
    rho_b, n_eligible, n_b_observed, l_core_exon, n_exon_only = _boundary_term_prorated(
        exon_in_locus=is_ex,
        region_ids_locus=region_ids,
        starts_locus=starts,
        ends_locus=ends,
        bf_left_locus=bf_l,
        bf_right_locus=bf_r,
        u_left=payload_arrays.u_left,
        u_right=payload_arrays.u_right,
        eff_clip_core=eff_clip_core,
        locus_start=locus_start,
        locus_end=locus_end,
        b_cross=b_cross,
        rho_global=global_densities.exon_intron.rho,
        kappa=global_densities.exon_intron.kappa.value,
    )
    if n_eligible == 0:
        fallback_flags |= FLAG_EXON_INTRON_NO_ELIGIBLE

    # 6. Total + decomposition invariant.
    n_gdna_total = n_mass_ig + n_mass_in + n_b_observed + n_exon_only

    # 7. π̂ with zero-n_obs safeguard (§3.5 step 7).
    if n_obs <= 0:
        pi_gdna = 0.0
    else:
        pi_unclipped = n_gdna_total / n_obs
        if pi_unclipped > 1.0:
            fallback_flags |= FLAG_PI_CLIPPED
        pi_gdna = float(min(1.0, max(0.0, pi_unclipped)))

    return LocusGdnaEstimate(
        locus=locus,
        n_obs=int(n_obs),
        n_gdna_intergenic=float(n_mass_ig),
        n_gdna_intron=float(n_mass_in),
        n_gdna_boundary_observed=float(n_b_observed),
        n_gdna_exon_only=float(n_exon_only),
        n_gdna=float(n_gdna_total),
        pi_gdna=pi_gdna,
        rho_loco=(float(rho_ig), float(rho_in), float(rho_b)),
        leff_loco=(float(l_core_ig), float(l_core_in), float(l_core_exon)),
        n_eligible_boundaries=int(n_eligible),
        n_boundary_events=float(n_b_observed),
        nrna_active=False,
        fallback_flags=int(fallback_flags),
    )


# ---------------------------------------------------------------------------
# Bayesian-prior redesign Phase 2: global-only η_g(ℓ) helper
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class ExpectedGdnaPriorParts:
    """Decomposition of the global-only expected gDNA pseudocount η_g(ℓ).

    Implements §2 of ``docs/bayesian_prior/bayesian_prior_plan_v3.md``:

    .. math::

        \\eta_g(\\ell) = \\rho^{ig} L^{ig}_\\ell
                       + \\rho^{in} L^{in}_\\ell
                       + \\rho^{b} ( s_\\ell B_{\\text{cross}} + L^{ex}_\\ell ).

    All four mass terms are *expected* counts under the relevant global
    density. None of them depend on locus-local fragment observations,
    making the resulting ``total`` a function of ``C_ℓ`` only.
    """

    total: float
    intergenic_contained: float
    intron_contained: float
    boundary_crossing_expected: float
    exon_contained_expected: float
    #: Density-estimation exposures consumed by the LOO trigger denominator.
    density_exposure_intergenic: float
    density_exposure_intron: float
    density_exposure_boundary: float


def expected_gdna_count_global(
    locus: Locus,
    region_index: RegionIndexPy,
    region_arrays: RegionArrays,
    global_densities: GlobalDensityTable,
    gdna_fl: FragmentLengthModel,
    *,
    ref_length: int | None = None,  # noqa: ARG001 (reserved; global-only path is unflanked)
    splicing_anchor_tolerance: int = 0,
    b_cross: float | None = None,
    _scratch: _LocusScratch | None = None,
) -> ExpectedGdnaPriorParts:
    """Pure global-only expected gDNA pseudocount for one ``Locus``.

    The function projects the *global* gDNA densities onto the locus's
    own FL-aware exposure. It is the v3 plan §2 prior formula and has
    **no dependency on per-locus fragment observations** (no
    ``payload_arrays`` argument). This is enforced at the type level so
    that local counts cannot leak into the ordinary global-only prior.

    Returns
    -------
    ExpectedGdnaPriorParts
        ``total`` is the canonical η_g(ℓ); the four mass-component
        fields and three density-exposure fields are exposed for
        diagnostics and the leave-one-locus-out wrapper.
    """
    if not (gdna_fl.mean > 0.0):
        raise ValueError(
            f"expected_gdna_count_global: gdna_fl.mean must be > 0; got {gdna_fl.mean!r}."
        )

    locus_start = int(locus.start)
    locus_end = int(locus.end)

    if _scratch is not None:
        # Reuse caller's scratch.  The flanked overlap is a superset of
        # the unflanked overlap; mass terms use ``eff_clip_core`` which
        # already clips to ``[locus_start, locus_end]`` and is therefore
        # 0 for any region that does not overlap the unflanked locus.
        # ``boundary_side_in_window`` likewise tests ``[locus_start,
        # locus_end]`` and naturally rejects flanked-only regions.
        starts = _scratch.starts
        ends = _scratch.ends
        bf_l = _scratch.bf_l
        bf_r = _scratch.bf_r
        eff_clip_core = _scratch.eff_clip_core
        is_ig = _scratch.is_ig
        is_in = _scratch.is_in
        is_ex = _scratch.is_ex
    else:
        region_ids = region_index.overlap(locus.ref_id, locus_start, locus_end)
        if region_ids.size == 0:
            # No overlapping regions ⇒ no exposure ⇒ η_g = 0. The caller
            # decides whether this is anomalous (it usually means a BAM
            # ref/index mismatch, which ``estimate_locus_gdna`` flags
            # loudly).
            return ExpectedGdnaPriorParts(
                total=0.0,
                intergenic_contained=0.0,
                intron_contained=0.0,
                boundary_crossing_expected=0.0,
                exon_contained_expected=0.0,
                density_exposure_intergenic=0.0,
                density_exposure_intron=0.0,
                density_exposure_boundary=0.0,
            )

        types = region_arrays.type[region_ids]
        starts = region_arrays.start[region_ids]
        ends = region_arrays.end[region_ids]
        bf_l = region_arrays.bf_left[region_ids]
        bf_r = region_arrays.bf_right[region_ids]

        # FL-PMF-weighted contained exposure: full + clipped to the locus.
        # The pure global-only path uses no flank — mass terms always live
        # inside the unflanked locus interval.
        _, eff_clip_core = contained_exposure_clipped(
            starts, ends, locus_start, locus_end, gdna_fl
        )

        is_ig = types == int(RegionType.INTERGENIC)
        is_in = types == int(RegionType.INTRON)
        is_ex = types == int(RegionType.EXON)

    l_core_ig = float(eff_clip_core[is_ig].sum()) if is_ig.any() else 0.0
    l_core_in = float(eff_clip_core[is_in].sum()) if is_in.any() else 0.0
    l_core_ex = float(eff_clip_core[is_ex].sum()) if is_ex.any() else 0.0

    # Boundary expected count: ρ_b · (s_ℓ · B_cross + L^ex_ℓ).
    # s_ℓ counts eligible boundary sides whose crossing point lies
    # inside the locus interval. Eligibility flags come from the
    # capture-tile annotation in the index (``bf_left``/``bf_right``).
    if is_ex.any():
        starts_ex = starts[is_ex]
        ends_ex = ends[is_ex]
        bf_l_ex = bf_l[is_ex]
        bf_r_ex = bf_r[is_ex]
        in_left, in_right = boundary_side_in_window(
            starts_ex, ends_ex, locus_start, locus_end
        )
        s_ell = int((bf_l_ex & in_left).sum() + (bf_r_ex & in_right).sum())
    else:
        s_ell = 0

    b_cross = boundary_crossing_exposure(
        gdna_fl, splicing_anchor_tolerance=splicing_anchor_tolerance
    ) if b_cross is None else float(b_cross)

    rho_ig = float(global_densities.intergenic.rho)
    rho_in = float(global_densities.intron.rho)
    rho_b = float(global_densities.exon_intron.rho)

    intergenic_contained = rho_ig * l_core_ig
    intron_contained = rho_in * l_core_in
    boundary_crossing_expected = rho_b * float(s_ell) * b_cross
    exon_contained_expected = rho_b * l_core_ex

    total = (
        intergenic_contained
        + intron_contained
        + boundary_crossing_expected
        + exon_contained_expected
    )

    return ExpectedGdnaPriorParts(
        total=float(total),
        intergenic_contained=float(intergenic_contained),
        intron_contained=float(intron_contained),
        boundary_crossing_expected=float(boundary_crossing_expected),
        exon_contained_expected=float(exon_contained_expected),
        density_exposure_intergenic=float(l_core_ig),
        density_exposure_intron=float(l_core_in),
        density_exposure_boundary=float(s_ell) * b_cross,
    )


def enable_gdna_for_multilocus(
    ml: MultiLocus,
    em_data: ScoredFragments,
) -> bool:
    """Decouples gDNA component eligibility from prior strength.

    A multi-locus is gDNA-eligible iff at least one of its EM units is
    *unspliced* and has a finite gDNA log-likelihood. This is the same
    condition the C++ extractor previously inferred from
    ``alpha_gdna > 0``; under the Bayesian-prior redesign it must be
    computed explicitly because ``alpha_gdna`` is now a pure
    pseudocount that may legitimately be zero.

    Must be called **before** ``partition_and_free`` because it reads
    the global per-unit ``is_spliced`` and ``gdna_log_liks`` arrays.
    """
    units = ml.unit_indices
    if units.size == 0:
        return False
    spliced = np.asarray(em_data.is_spliced[units], dtype=bool)
    return bool(np.any((~spliced) & np.isfinite(em_data.gdna_log_liks[units])))


# ---------------------------------------------------------------------------
# MultiLocus aggregation + orchestrator
# ---------------------------------------------------------------------------

def assemble_multilocus_prior(
    ml: MultiLocus,
    per_locus_estimates: tuple[LocusGdnaEstimate, ...],
    *,
    gdna_prior_count: float,
    rna_prior_count: float,
) -> MultiLocusPrior:
    """Aggregate per-``Locus`` estimates into one ``MultiLocusPrior``.

    ``gdna_prior_count`` / ``rna_prior_count`` are the canonical
    asymmetric Dirichlet pseudocounts the EM consumes (Phase 1
    schema).
    """
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
        gdna_prior_count=float(gdna_prior_count),
        rna_prior_count=float(rna_prior_count),
        per_locus=per_locus_estimates,
    )


def assemble_priors(
    multi_loci: list[MultiLocus],
    em_data: ScoredFragments,
    index: "TranscriptIndex",
    payload: CalibrationScanPayload,
    global_densities: GlobalDensityTable,
    *,
    gdna_fl: FragmentLengthModel | None = None,
    nrna_weight: float = 0.0,
    intergenic_flank_bp: int = INTERGENIC_FLANK_BP_DEFAULT,
    splicing_anchor_tolerance: int = 0,
) -> PriorTable:
    """Build the full :class:`PriorTable` for the batch EM.

    Pure-Python pass; no I/O, no mutation of inputs.

    Bayesian-prior redesign Phase\u00a02 semantics:

    * The canonical prior is the *global-only* expected gDNA pseudocount
      :math:`\\eta_g(\\ell)` from
      :func:`expected_gdna_count_global`. ``alpha_gdna`` is set to that
      count summed across the multi-locus's constituent loci;
      ``alpha_rna`` is fixed at ``0`` (uninformative RNA prior).
    * gDNA component *eligibility* is a separate ``uint8`` flag,
      :func:`enable_gdna_for_multilocus`, computed from per-unit
      ``is_spliced`` and ``gdna_log_liks`` arrays before
      ``partition_and_free`` consumes them. This decouples eligibility
      from prior strength (Phase\u00a00 native ABI change).
    * The legacy locoregional ``estimate_locus_gdna`` is still called
      to populate the per-locus diagnostic dataframe (``pi_gdna``,
      mass decomposition, fallback flags), but its output is **not**
      consumed by ``alpha_gdna``/``alpha_rna``.

    Parameters
    ----------
    gdna_fl : FragmentLengthModel, optional
        gDNA fragment-length distribution used for the FL-PMF-weighted
        containment effective length. Defaults to
        ``global_densities.gdna_fl``.
    nrna_weight : float
        Per-component weight applied to synthetic-nRNA transcripts in
        each ``MultiLocus``.  See :func:`build_prior_weight_rna`.
    intergenic_flank_bp : int
        Padding (bp) added on both sides of each locus when querying
        regions for the *intergenic* density's evidence window in the
        diagnostic ``estimate_locus_gdna`` pass. The Phase\u00a02 prior
        itself is unflanked; flank-based priors arrive in Phase\u00a04
        as ``gdna_prior_source = "independent_flank"``.
    """
    if gdna_fl is None:
        gdna_fl = global_densities.gdna_fl
    if not (gdna_fl.mean > 0.0):
        raise ValueError(
            f"assemble_priors: gdna_fl.mean must be > 0; got {gdna_fl.mean!r}."
        )

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    region_index = RegionIndexPy(arrays=region_arrays)

    # Hoist boundary-crossing exposure out of the per-locus loop:
    # depends only on (gdna_fl, splicing_anchor_tolerance), both
    # constant across the multi-loci. Eliminates ~one O(max_size)
    # numpy reduction per Locus.
    b_cross = boundary_crossing_exposure(
        gdna_fl, splicing_anchor_tolerance=splicing_anchor_tolerance
    )

    # Per-ref length lookup for clamping the intergenic flank window
    # (diagnostic path only). Phase\u00a02's canonical prior is unflanked.
    if intergenic_flank_bp > 0:
        ref_lengths_arr = np.asarray(
            list(index.ref_lengths.values()), dtype=np.int64
        )
    else:
        ref_lengths_arr = None

    # Hoist per-transcript columns once (slow path only consults these).
    t_starts = index.t_df["start"].values
    t_ref_codes = index.t_df["ref"].cat.codes.values
    cat_to_ref_id = np.array(
        [index.ref_name_to_id[str(name)] for name in index.t_df["ref"].cat.categories],
        dtype=np.int32,
    )
    t_ref = cat_to_ref_id[t_ref_codes]

    if "is_synthetic" in index.t_df.columns:
        is_synthetic_arr = index.t_df["is_synthetic"].to_numpy(dtype=bool)
    else:
        is_synthetic_arr = None

    n_ml = len(multi_loci)
    alpha_gdna = np.zeros(n_ml, dtype=np.float64)
    alpha_rna = np.zeros(n_ml, dtype=np.float64)
    gdna_prior_count_arr = np.zeros(n_ml, dtype=np.float64)
    rna_prior_count_arr = np.zeros(n_ml, dtype=np.float64)
    enable_gdna_arr = np.zeros(n_ml, dtype=np.uint8)
    multi_locus_priors: list[MultiLocusPrior | None] = [None] * n_ml
    prior_weight_rna: list[np.ndarray | None] = [None] * n_ml

    for ml in multi_loci:
        idx = int(ml.multi_locus_id)
        if multi_locus_priors[idx] is not None:
            raise RuntimeError(
                f"assemble_priors: duplicate multi_locus_id={idx} in multi_loci."
            )

        # 1. Per-locus diagnostic (legacy locoregional pi_gdna estimate).
        t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
        units_per_locus = partition_units_to_loci(ml, em_data, t_to_local)
        # Pre-compute per-Locus shared scratch (region overlap + FL-aware
        # contained exposures) once per locus.  Both
        # ``estimate_locus_gdna`` and ``expected_gdna_count_global``
        # consume it via the ``_scratch`` kwarg, halving the per-locus
        # ``region_index.overlap`` and ``contained_exposure_clipped``
        # calls.  When the locus has no overlapping regions the scratch
        # is ``None`` and each helper falls back to its legacy path
        # (raise / zeroed result).
        loci_scratch = tuple(
            _compute_locus_scratch(
                loc,
                region_index,
                region_arrays,
                gdna_fl,
                intergenic_flank_bp=intergenic_flank_bp,
                ref_length=int(ref_lengths_arr[loc.ref_id])
                if ref_lengths_arr is not None
                and 0 <= loc.ref_id < ref_lengths_arr.size
                else None,
            )
            for loc in ml.loci
        )
        per_locus_est = tuple(
            estimate_locus_gdna(
                locus=loc,
                n_obs=int(units_per_locus[j].size),
                region_index=region_index,
                region_arrays=region_arrays,
                payload_arrays=payload_arrays,
                global_densities=global_densities,
                gdna_fl=gdna_fl,
                intergenic_flank_bp=intergenic_flank_bp,
                ref_length=int(ref_lengths_arr[loc.ref_id])
                if ref_lengths_arr is not None
                and 0 <= loc.ref_id < ref_lengths_arr.size
                else None,
                splicing_anchor_tolerance=splicing_anchor_tolerance,
                b_cross=b_cross,
                _scratch=scratch,
            )
            for j, (loc, scratch) in enumerate(zip(ml.loci, loci_scratch, strict=True))
        )

        # 2. Canonical Phase\u00a02 prior: η_g summed across constituent
        #    loci; α_rna fixed at 0; eligibility from the unit-level
        #    helper (independent of prior strength).
        eta_g = float(
            sum(
                expected_gdna_count_global(
                    locus=loc,
                    region_index=region_index,
                    region_arrays=region_arrays,
                    global_densities=global_densities,
                    gdna_fl=gdna_fl,
                    splicing_anchor_tolerance=splicing_anchor_tolerance,
                    b_cross=b_cross,
                    _scratch=scratch,
                ).total
                for loc, scratch in zip(ml.loci, loci_scratch, strict=True)
            )
        )
        eta_r = 0.0

        ml_prior = assemble_multilocus_prior(
            ml, per_locus_est,
            gdna_prior_count=eta_g,
            rna_prior_count=eta_r,
        )
        multi_locus_priors[idx] = ml_prior
        gdna_prior_count_arr[idx] = eta_g
        rna_prior_count_arr[idx] = eta_r
        alpha_gdna[idx] = eta_g
        alpha_rna[idx] = eta_r
        enable_gdna_arr[idx] = np.uint8(
            1 if enable_gdna_for_multilocus(ml, em_data) else 0
        )
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
        gdna_prior_count=gdna_prior_count_arr,
        rna_prior_count=rna_prior_count_arr,
        enable_gdna=enable_gdna_arr,
    )
