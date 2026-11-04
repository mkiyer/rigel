"""What a chain slot is before any message passing: its geometry, its static masks and its initial belief.

Layer 3, beneath the sweep. Everything here is a pure function of the accumulator substrate, the region geometry, the
sj axis and the fragment-length pmfs, with no sweep state and no global prior. There is no per-face quantity: a
contiguous boundary is 0 bp wide, so it has one set of numbers seen identically from both sides and one crossing
divisor, and the index states each sj's ``(src, dst, strand)``, so the chain's adjacency says which flank a spliced
deposit belongs to. The accumulator deposits ``+1`` on every object a fragment touched, so a slot's count is both
its density numerator and its Poisson ``n``, and ``Var(log rho) = 1/n`` holds against it. The boundary flags pass
through :class:`RegionStatics` uninterpreted; the transfer builders read them (``native/transfer_kernel.h``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import Strand
from .effective_length import (
    crossing_eff_length,
    observable_contained_eff_length,
    observable_crossing_eff_length,
)
from .region_chain import BOUNDARY, REGION, RegionChain
from .signature import (
    TS_AMBIG,
    TS_NEG,
    TS_POS,
    rna_active_strands,
)
from .simplex_logodds import _solve_regions_logodds_all

__all__ = [
    "RegionGeometry",
    "build_region_geometry",
    "region_gdna_geometry",
    "RegionBelief",
    "RegionStatics",
    "build_region_statics",
    "init_beliefs",
    "g1_locked",
]


@dataclass(frozen=True, slots=True)
class RegionGeometry:
    """Per-slot static geometry, all arrays of length ``n_slots``: one set of numbers per slot.

    The three populations carry the names of the accumulator's three banks (``boundary_unspliced``,
    ``boundary_spliced``, ``sj_count`` in ``tests/native/_accumulator_reference.py``). Of the molecules that touched
    one boundary::

        unspliced_count   crossed it contiguously, spliced nowhere      a gDNA/RNA mixture
        spliced_count     crossed it contiguously, spliced elsewhere    certified RNA
        sj_count          never crossed it: it jumped from here         certified RNA

    The last two are different molecules: at a donor the sj flux is all the RNA that splices there, while the
    spliced crossing is the few molecules that read through; inside an exon there is no sj, so ``sj_count`` is 0.

    The strand axes differ. ``unspliced_count`` and ``spliced_count`` are keyed by genome strand (where the read
    aligned); ``sj_count`` is keyed by transcript strand, derived from each sj's annotated strand.
    """

    n_slots: int

    #: float64[n_slots, 2] — ``boundary_unspliced`` at a BOUNDARY, ``region_contained`` at a REGION, by genome
    #: strand: the gDNA/RNA mixture being deconvolved, the only population that is one.
    unspliced_count: np.ndarray
    #: float64[n_slots] — the gDNA divisor: contained placements at a REGION, crossing placements at a BOUNDARY.
    eff_gdna: np.ndarray
    #: float64[n_slots] — the RNA divisor, the same frames on the RNA length pmf. Both divisors are 0 where there is
    #: no opportunity, never floored: a consumer divides only where the divisor is positive.
    eff_rna: np.ndarray
    #: float64[n_slots, 2] — ``boundary_spliced``, by genome strand: molecules that crossed this boundary
    #: contiguously and spliced elsewhere. Certified RNA, since gDNA cannot splice, so it floors the RNA in this
    #: boundary's population. 0 on every REGION slot, because ``region_contained`` is credited only to a fragment
    #: that used no sj. Its divisor is :attr:`eff_rna`: a contiguous crossing is one whatever the molecule did
    #: elsewhere.
    spliced_count: np.ndarray
    #: float64[n_slots, 2] — ``sj_count`` gathered onto the two boundaries each sj leaves and enters, by transcript
    #: strand: the flux that changes template here.
    sj_count: np.ndarray
    #: float64[n_slots, 2] — the same flux split by which genomic end of its sj this boundary is: ``_lo`` for sj
    #: whose low end is here, ``_hi`` for the high end. The halves belong to different flanks: a molecule that
    #: splices here has its body on one side only, so a face carries ``_lo`` toward the low neighbour and ``_hi``
    #: toward the high one; their sum, ``sj_count``, is the whole flux leaving the boundary. The names must stay
    #: genomic, never donor/acceptor: the index's ``FLAG_DONOR_s`` marks the genomic-low end of an intron on both
    #: strands, which on ``−`` is the transcript's acceptor.
    sj_count_lo: np.ndarray
    sj_count_hi: np.ndarray
    #: float64[n_slots, 2] — the certified rate per face, ``Σ_J flux_J / A_J`` over the face's junctions, by
    #: transcript strand. The junctions at one face are disjoint routes (each molecule crosses one), so their rates
    #: sum. No pseudocount: a zero-flux face reads exactly 0.
    route_rate_lo: np.ndarray
    route_rate_hi: np.ndarray


def build_region_geometry(
    chain: RegionChain,
    substrate,
    region_arrays,
    sj,
    gdna_fl_pmf: np.ndarray,
    rna_fl_pmf: np.ndarray,
) -> RegionGeometry:
    """Assemble the per-slot geometry from the substrate's populations onto the chain.

    The divisors:

    ================  ==========================================================================
    REGION, both      the contained placements over the region's observable pieces
                      (``observable_contained_eff_length``); a region masked whole has none
    BOUNDARY, gDNA    the crossing placements within the bases to the nearest masked base either
                      side, ``mu_g − 1`` where nothing is masked (``observable_crossing_eff_length``);
                      gDNA's template is the chromosome, so it is untapered
    BOUNDARY, RNA     the same observable reaches, ``mu_r − 1`` where nothing is masked. The RNA
                      half of a crossing is bounded by its transcript, so this untapered form
                      overstates RNA's opportunity and leans the split toward gDNA
    sj, RNA           the exonic per-strand reach (``crossing_eff_length``): only a molecule that
                      spliced across a sj uses it
    ================  ==========================================================================

    A sj's donor is the boundary right of ``src_region`` and its acceptor the boundary left of ``dst_region``,
    read off the chain's adjacency. A reference terminal links to ``-1``, so a sj can never attach across
    references: this function raises instead.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    n = chain.n_slots

    # ── the two contiguous populations: the mixture, and the certified-RNA floor beside it ───────
    unspliced_count = np.zeros((n, 2), dtype=np.float64)
    unspliced_count[is_region] = np.asarray(substrate.region_contained.count, np.float64)[
        obj[is_region]
    ]
    unspliced_count[is_boundary] = np.asarray(substrate.boundary_unspliced.count, np.float64)[
        obj[is_boundary]
    ]
    # REGION slots stay 0: a contained fragment used no sj, so it cannot be a spliced molecule
    spliced_count = np.zeros((n, 2), dtype=np.float64)
    spliced_count[is_boundary] = np.asarray(substrate.boundary_spliced.count, np.float64)[
        obj[is_boundary]
    ]

    # ── the two per-component divisors ───────────────────────────────────────────────────────────
    n_regions = int(region_arrays.n_regions)

    def divisor(pmf: np.ndarray) -> np.ndarray:
        """Per-slot effective length over the observable genome: contained placements at a REGION, crossing
        placements at a BOUNDARY."""
        out = np.zeros(n, dtype=np.float64)
        if n_regions:
            out[is_region] = observable_contained_eff_length(region_arrays, pmf)[obj[is_region]]
        if is_boundary.any():
            out[is_boundary] = observable_crossing_eff_length(region_arrays, pmf)[obj[is_boundary]]
        return out

    eff_gdna = divisor(gdna_fl_pmf)
    eff_rna = divisor(rna_fl_pmf)

    # ── the jumping population: a sj is a factor on the boundaries it leaves and enters ───
    sj_count = np.zeros((n, 2), dtype=np.float64)
    jc_lo = np.zeros((n, 2), dtype=np.float64)
    jc_hi = np.zeros((n, 2), dtype=np.float64)
    route_rate_lo = np.zeros((n, 2), dtype=np.float64)
    route_rate_hi = np.zeros((n, 2), dtype=np.float64)
    if sj.n_sj:
        slot_of_region = np.zeros(int(chain.n_regions_total), dtype=np.int64)
        slot_of_region[obj[is_region]] = np.flatnonzero(is_region)
        donor = np.asarray(chain.right)[slot_of_region[np.asarray(sj.src_region, np.int64)]]
        acceptor = np.asarray(chain.left)[slot_of_region[np.asarray(sj.dst_region, np.int64)]]
        if np.any(donor < 0) or np.any(acceptor < 0):
            raise ValueError(
                "a sj boundary attaches to a reference terminal, which has no boundary beside it. "
                "Both of an intron's endpoints are interior interfaces of the same reference by "
                "construction (splice_graph I5), so this is a sj axis addressing a different "
                "partition than the payload was scanned on."
            )
        # the sj's count summed over the genome-strand columns, filed under the sj's transcript strand
        flux = np.asarray(substrate.sj.count, np.float64).sum(axis=1)
        eff = crossing_eff_length(rna_fl_pmf, sj.reach_lo, sj.reach_hi)
        column = np.where(np.asarray(sj.strand) == np.int8(Strand.POS), 0, 1)
        # ``donor`` is the boundary at the sj's genomic-low end and ``acceptor`` the high one on both strands:
        # the chain's links are genomic and a sj runs ``src < dst``
        live_route = eff > 0.0
        for boundary, rr in ((donor, route_rate_lo), (acceptor, route_rate_hi)):
            np.add.at(
                rr,
                (boundary[live_route], column[live_route]),
                flux[live_route] / eff[live_route],
            )
        for boundary, jc in ((donor, jc_lo), (acceptor, jc_hi)):
            np.add.at(sj_count, (boundary, column), flux)
            np.add.at(jc, (boundary, column), flux)

    return RegionGeometry(
        n_slots=int(n),
        unspliced_count=unspliced_count,
        eff_gdna=eff_gdna,
        eff_rna=eff_rna,
        spliced_count=spliced_count,
        sj_count=sj_count,
        sj_count_lo=jc_lo,
        sj_count_hi=jc_hi,
        route_rate_lo=route_rate_lo,
        route_rate_hi=route_rate_hi,
    )


def region_gdna_geometry(geometry: RegionGeometry):
    """Per-slot gDNA support ``(unspliced count, eff_gdna)``.

    ``calibrate`` fits the landscape prior on it and derives the λ bracket from it, and the sweep's diagnostics
    capture publishes it as ``mass_global``. The kernel (``native/solve_kernel.cpp``) forms the same pair from the
    arrays :func:`sweep.solve_chain` hands it, so the two definitions must stay equal. The divisor is ``eff_gdna``,
    so ``rho_g = f_g·M/E_g`` and ``f_g`` is a count share, not a density share. A total density must never be formed
    as ``mass / eff_gdna``: that divides a total by one component's opportunity.
    """
    return np.asarray(geometry.unspliced_count, np.float64).sum(axis=1), np.asarray(
        geometry.eff_gdna, np.float64
    )


@dataclass(frozen=True, slots=True)
class RegionBelief:
    """Per-slot solved state on the chain, all arrays of length ``n_slots`` (the region+boundary chain, not the
    region axis).

    The composition ``(f_pos, f_neg, f_g)`` of the slot's unspliced mass, stored as fractions because a boundary
    has two faces but one composition, and ``var_gdna = Var(log f_g)``, a grid moment over the λ lattice. The
    variance is in log-fraction space and not bounded by ¼; the delta method gives ``Var(f_g) ≈ f_g²·Var(log f_g)``.
    It is the precision state: 0 is locked and certain, ``inf`` is no information. The landscape prior reads it as
    its training weight and location floor (``calibrate._fit_gdna_hyperprior``); no message reads it.
    """

    f_pos: np.ndarray
    f_neg: np.ndarray
    f_g: np.ndarray
    var_gdna: np.ndarray
    #: Whether the slot holds a composition: its own composition channel, structural certainty, or a composition
    #: row from a neighbour. A level, a ceiling or a cube row is only a bound and does not count, because a
    #: bound-only slot settles where the prior puts it, and training the prior on it would train it on its own
    #: echo. Published by ``sweep.solve_chain``; ``None`` before any solve, which the landscape's selector reads as
    #: the annotation alone.
    #: Gate: ``tests/calibration/test_landscape_training_population.py``.
    has_composition: np.ndarray | None = None


# ---------------------------------------------------------------------------
# Initialization — the signature-binary G1/G2/G3 belief on the chain.
# ---------------------------------------------------------------------------


def g1_locked(free_pos, free_neg) -> np.ndarray:
    """The G1 class: a structurally pure-gDNA slot, where neither RNA strand is admissible.

    :func:`_type_belief` pins these slots at ``{0,0,1}`` with ``Var(log f_g) = 0``. It applies to regions and
    boundaries alike: an intergenic region and an intergenic|exon boundary are both G1. This is the one definition,
    and every instrument imports it. The predicate says the composition is certain, not that the count is clean: a
    G1 boundary can still sit beside RNA-carrying exons.
    """
    return ~np.asarray(free_pos, bool) & ~np.asarray(free_neg, bool)


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """The per-slot composition ``(f_pos, f_neg, f_g, var_gdna)`` from the signature class and the strand-only solve.

    ``free_pos``/``free_neg`` say whether each strand's RNA is admissible; ``deconv`` is the strand-only
    :class:`RegionDeconv`, with no global prior. The default is all-gDNA ``{0,0,1}``:

    * G1 (neither strand free): a locked gDNA sink, ``{0,0,1}`` at variance 0.
    * G2 (exactly one strand free, with data): the strand solve alone resolves it, with ``f_active = 1 − f_g``.
    * G3 (both strands free): unresolvable by strand, so ``{0,0,1}`` at variance ``inf``; the sweep resolves it.
    """
    n = free_pos.shape[0]
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g = np.ones(n)  # the signature-binary all-gDNA default; the count plays no role in it
    var_g = np.full(n, np.inf)

    g1 = g1_locked(free_pos, free_neg)
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2 with data takes the strand-only solve; G3 keeps the default at variance inf
    fgv = np.asarray(deconv.gdna_frac_var, dtype=np.float64)
    f_g[g2_active] = np.asarray(deconv.gdna_frac, dtype=np.float64)[g2_active]
    f_pos[g2_active] = np.asarray(deconv.rna_pos_frac, dtype=np.float64)[g2_active]
    f_neg[g2_active] = np.asarray(deconv.rna_neg_frac, dtype=np.float64)[g2_active]
    var_g[g2_active] = fgv[g2_active]

    var_g[g1] = 0.0
    return f_pos, f_neg, f_g, var_g


@dataclass(frozen=True, slots=True)
class RegionStatics:
    """Per-slot structural masks (length ``n_slots``) that never change during the sweep; no count lives here.

    ``free_pos``/``free_neg`` are the strands on which RNA may be present at all: a region's own transcript bits,
    a boundary's strand continuity. ``boundary_flags`` carries the splice graph's 8 structural bits
    (``TSS_s``/``TES_s``/``DONOR_s``/``ACCEPTOR_s``) at each BOUNDARY slot and 0 on REGION slots. They stay raw bits,
    because each consumer reads a different combination.
    """

    n_slots: int
    free_pos: np.ndarray  # bool — RNA-active (transcript continuity); the RNA-crossing gate
    free_neg: np.ndarray  # bool
    boundary_flags: np.ndarray  # uint16 — graph structural bits; 0 on REGION slots


def build_region_statics(
    chain: RegionChain,
    region_arrays,
    boundary_flags: np.ndarray | None = None,
) -> RegionStatics:
    """Gather the structural masks onto the chain in one slot-keyed pass over the signatures and adjacency.

    A boundary always has a region on both sides, so ``chain.left``/``chain.right`` give its flanks. The mask is a
    continuity gate: a strand-``s`` unspliced crossing can be RNA only where strand ``s`` is present on both flanks,
    which makes a TSS/TES boundary (intergenic|exon) a gDNA sink.

    ``boundary_flags`` is the per-contiguous-boundary ``uint16[E]`` from
    :func:`~rigel.calibration.splice_graph.build_boundary_flags_array`; ``None`` leaves the field zero, read as no
    structural information.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    n = int(chain.n_slots)
    flags = _check_boundary_flags(boundary_flags, int(chain.n_boundaries_total))

    sig = np.asarray(region_arrays.signature).astype(np.int64)
    n_regions = sig.shape[0]
    region_idx = np.clip(obj, 0, max(n_regions - 1, 0))

    # the signature at each REGION slot, then read through the chain's adjacency at each BOUNDARY slot
    slot_sig = np.where(is_region, sig[region_idx] if n_regions else 0, 0)
    left = np.clip(np.asarray(chain.left), 0, max(n - 1, 0))
    right = np.clip(np.asarray(chain.right), 0, max(n - 1, 0))
    sig_l = np.where(is_boundary, slot_sig[left], 0)
    sig_r = np.where(is_boundary, slot_sig[right], 0)

    ts = np.where(
        is_region, np.asarray(region_arrays.strand_class)[region_idx] if n_regions else 0, -1
    )
    pos_l, neg_l = rna_active_strands(sig_l)
    pos_r, neg_r = rna_active_strands(sig_r)

    free_pos = np.where(is_region, (ts == TS_POS) | (ts == TS_AMBIG), pos_l & pos_r)
    free_neg = np.where(is_region, (ts == TS_NEG) | (ts == TS_AMBIG), neg_l & neg_r)

    return RegionStatics(
        n_slots=n,
        # No per-slot spliced floor: the message layer owns spliced RNA, and a local floor would count it twice,
        # putting RNA that spliced elsewhere into a boundary's unspliced f_pos.
        free_pos=free_pos,
        free_neg=free_neg,
        boundary_flags=np.where(
            is_boundary, flags[np.clip(obj, 0, max(flags.shape[0] - 1, 0))], 0
        ).astype(np.uint16),
    )


def _check_boundary_flags(boundary_flags, n_boundaries: int) -> np.ndarray:
    """Validate the per-contiguous-boundary flags against the chain; a mis-sized array would shift every flag."""
    if boundary_flags is None:
        return np.zeros(max(n_boundaries, 1), dtype=np.uint16)
    flags = np.asarray(boundary_flags, dtype=np.uint16)
    if flags.shape != (n_boundaries,):
        raise ValueError(
            f"boundary_flags has shape {flags.shape}; expected ({n_boundaries},), one per contiguous boundary. "
            f"Build it with splice_graph.build_boundary_flags_array(index) against the SAME index the "
            f"payload was scanned on. There are no terminal slots: a reference with k regions owns "
            f"k-1 boundaries, not k+1."
        )
    return flags if flags.size else np.zeros(1, dtype=np.uint16)


def init_beliefs(
    geometry: RegionGeometry,
    statics: RegionStatics,
    *,
    rna_sense_frac: float,
    n_grid: int,
    logodds_window: float,
    n_threads: int = 1,
) -> RegionBelief:
    """The initial G1/G2/G3 :class:`RegionBelief` on the chain.

    Only single-strand slots are solved (:mod:`simplex_logodds`, one tilt column): the binomial strand likelihood
    plus the Jeffreys reference, with no global prior. :func:`_type_belief` then sets the class beliefs: a
    single-strand slot with no gDNA resolves to ``f_g≈0`` from the strand alone, G1 sinks lock at ``{0,0,1}``, and
    both-strand slots hold ``{0,0,1}`` at variance ``inf`` for the sweep, so they are masked out of the solve.
    """
    st = statics
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    single = fp ^ fn
    count = np.asarray(geometry.unspliced_count, np.float64)
    deconv = _solve_regions_logodds_all(
        count[:, 0],
        count[:, 1],
        fp & single,
        fn & single,
        kappa=float(rna_sense_frac),
        n_grid=n_grid,
        L=logodds_window,
        n_threads=n_threads,
    )
    f_pos, f_neg, f_g, var_g = _type_belief(fp, fn, deconv, count.sum(axis=1))
    return RegionBelief(f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_gdna=var_g)


# gDNA density carries no per-capture-class scaling: a gDNA level crosses a composition-unlicensed hop unscaled (the
# transfer policy's level lane), because a pure-gDNA object's own observed total already is its gDNA density at its
# own capture.
