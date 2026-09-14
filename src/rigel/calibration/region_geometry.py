"""rigel.calibration.region_geometry — per-slot chain geometry, beliefs, densities, statics, and init.

The lower layer beneath the belief-propagation backbone (`sweep`): everything that describes *what a
chain slot is* before any message passing. Pure functions of the accumulator substrate + the region
geometry + the sj axis + the FL pmfs — no sweep state, no global prior.

There is no per-FACE concept anywhere here, and that is structural rather than a simplification. A
contiguous boundary is a 0-bp boundary: one set of numbers, seen identically from both sides, so its two
sides cannot carry different divisors. What that dissolves:

* sj-strand routing and exon-bit flank gating, which would exist only to *guess* which flank a spliced
  deposit belonged to. The index states ``(src, dst, strand)`` explicitly, so the chain's own adjacency
  answers it;
* any choice between two per-face spliced divisors — there is one crossing formula;
* ``mass`` versus ``n``. The accumulator deposits ``+1`` on every object the fragment touched, so
  ``count`` is both the density numerator and the Poisson ``n``, and ``Var(log rho) = 1/n`` is honest
  against it. A fractional mass would need a separate integer flux carried beside it.

Contents:

* `RegionGeometry` / `build_region_geometry` — the per-slot static geometry: the unspliced count being
  deconvolved, its two per-component divisors, and the mature (sj) flux with its own.
* `RegionBelief` — the per-slot pie ``(f_pos, f_neg, f_g)``; the per-component message densities
  ``rho = f·M/E`` are computed inline in ``sweep.solve_chain``.
* `RegionStatics` / `build_region_statics` — the static per-slot solver inputs (per-strand counts, masks).
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.

Layering: LAYER 3. It imports DOWN to `region_chain` and `signature` (0) and `effective_length` (2),
and SIDEWAYS to `simplex_logodds` (3) — never `sweep` (layer 6) or `landscape` (layer 5), so it sits
cleanly below both. The flags array travels through `RegionStatics` without this module knowing what
any bit means: the terminus and junction bits are read by `messages.transfer_rows`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import Strand
from .effective_length import UNBOUNDED_REACH, contained_eff_length, crossing_eff_length
from .region_chain import BOUNDARY, REGION, RegionChain
from .signature import (
    TS_AMBIG,
    TS_NEG,
    TS_POS,
    mrna_active_strands,
    nrna_active_strands,
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

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class RegionGeometry:
    """Per-chain-slot static geometry (length ``n_slots``). One set of numbers per slot.

    Three populations, named for the accumulator's own three banks: ``unspliced`` / ``spliced`` / ``sj``
    is what ``boundary_unspliced``, ``boundary_spliced`` and ``sj_count`` are called in the executable
    specification (`tests/native/_accumulator_reference.py`), and a consumer that renames a bank on the
    way in is how one quantity comes to have two names. At one boundary, of the molecules that touched
    it::

        unspliced_count   crossed it CONTIGUOUSLY, spliced nowhere       a gDNA/RNA MIXTURE
        spliced_count     crossed it CONTIGUOUSLY, spliced elsewhere     certified RNA
        sj_count    never crossed it -- it JUMPED from here        certified RNA

    The last two are different molecules and routinely differ by two orders of magnitude at the same
    boundary. At a donor boundary the sj flux is the gene's whole mature output while the spliced
    crossing is the handful of molecules that read through without splicing; at a region bound inside an
    exon there is no sj at all, so ``sj_count`` is 0 while ``spliced_count`` carries everything. The
    word *mature* fits both and therefore distinguishes neither, which is why it is not used.

    Two strand axes, and they are not the same axis. ``unspliced_count`` and ``spliced_count`` are keyed
    by GENOME strand — where the read aligned, the accumulator's one storage convention.
    ``sj_count``/``eff_sj`` are keyed by TRANSCRIPT strand, derived from each sj's own annotated strand.
    That derivation is why the accumulator stores no sense/antisense column, and putting the two axes in
    one array under one name would be two conventions in one schema.
    """

    n_slots: int

    #: float64[n_slots, 2] — ``boundary_unspliced`` at a BOUNDARY, ``region_contained`` at a REGION. By GENOME
    #: strand. The mixture being deconvolved, and the only population that is one.
    #: It is both the density numerator and the Poisson ``n``: the accumulator deposits ``+1`` on every
    #: object the fragment touched, so there is no fractional mass to carry separately.
    unspliced_count: np.ndarray
    #: float64[n_slots] — the gDNA divisor. Contained placements at a REGION, crossing placements at a BOUNDARY.
    eff_gdna: np.ndarray
    #: float64[n_slots] — the RNA divisor, the same frames on the RNA length pmf.
    #: Both are 0 where there is no opportunity, never floored: a consumer divides only where
    #: the divisor is positive, and an object with no opportunity for a component emits nothing.
    eff_rna: np.ndarray
    #: float64[n_slots, 2] — ``boundary_spliced``: molecules that crossed this boundary CONTIGUOUSLY
    #: having spliced somewhere else in the same molecule. By GENOME strand. Certified RNA, since gDNA
    #: cannot splice, which is what makes it a floor on the RNA inside this boundary's own population.
    #: 0 on every REGION slot, structurally: ``region_contained`` is credited only when the fragment used no
    #: sj at all (`_accumulator_reference.deposit`).
    #: Its divisor is :attr:`eff_rna` — a contiguous crossing is a contiguous crossing whatever the
    #: molecule did elsewhere — so it needs no effective length of its own.
    spliced_count: np.ndarray
    #: float64[n_slots, 2] — ``sj_count``, gathered onto the two boundaries each sj leaves and enters,
    #: by TRANSCRIPT strand. The flux that changes template here: what the SPLICE IN hands an exon, and
    #: what the SPLICE OUT measures the continuing share against.
    sj_count: np.ndarray
    #: float64[n_slots, 2] — the SUMMED sj divisor, same keying. Several sj on one boundary are
    #: several estimates of one rate, so the pooled statement is ``Σcount / ΣE`` — the ratio of sums, never
    # the mean of ratios (``ρ_bg = Σg/ΣE``).
    eff_sj: np.ndarray
    #: float64[n_slots, 2] — the SAME flux, split by which genomic END of its sj this BOUNDARY is.
    #: ``_lo`` is the flux of sj whose genomic-LOW end is here; ``_hi`` the genomic-HIGH end's.
    #: The two halves belong to DIFFERENT FLANKS of the same BOUNDARY: a molecule that splices at
    #: this position has its body in the exon on ONE side of it — the low side if this is the sj's low
    #: end, the high side if this is its high end — and it never enters the other flank at all. So the
    #: flux a face carries is ``_lo`` toward the low neighbour and ``_hi`` toward the high one, which is
    #: how the transfer policy reads them; summing them, as ``sj_count`` above does, is right for a
    #: claim about the whole flux leaving this boundary.
    #: Written in GENOMIC terms, never donor/acceptor: the index's ``FLAG_DONOR_s`` bit marks the
    #: genomic-LOW end of an ``s``-strand intron on BOTH strands, so on ``−`` it sits at the transcript's
    #: biological ACCEPTOR. Naming these ``_lo``/``_hi`` is what stops that being a sign error.
    #: The reciprocal-opportunity total, counts/bp — model-free at a BOUNDARY, truncated at a REGION.
    #: The accumulator deposits ``1/A(w)`` per fragment, so ``E[sum] = rho * P(A > 0)``
    #: (`tests/native/_accumulator_reference.py`): at a BOUNDARY ``A = w−1`` and ``P(w ≥ 2) = 1`` on any
    #: real library, so the boundary slots ARE the density, exactly, for any length distribution; at a
    #: REGION ``A = (ell−w+1)₊`` and the slot reads ``rho * P(w ≤ ell)``, a per-component pmf functional
    #: (an order of magnitude at a short exon, exactly 0 below ``frag_min``) — a density SHAPE, not a
    #: level.
    #: ⛔ A total abundance must still never be formed as ``mass / effective_length``: that divisor
    #: depends on which component the fragments came from — 100 counts in a 500 bp region reads 0.25 as
    #: pure gDNA and 0.33 as pure RNA at mean lengths 100/200 — which is circular. The boundary-form bank
    #: reads the same number either way; the REGION form moves the same circularity into its SUPPORT
    #: rather than removing it. A truly untruncated REGION total is
    #: `total_abundance.region_counts_and_exposure`'s START/END pair, not this bank.
    inv_abundance: np.ndarray
    #: The sj flux's own model-free abundance ``[n, 2]`` BY TRANSCRIPT STRAND, split by which
    #: genomic END of its sj this BOUNDARY is — the same reciprocal-opportunity deposit
    #: (`sj_inv_length_sum`), so it is in the SAME units as :attr:`inv_abundance`: sum the strands and
    #: add for a FACE's total, or read one column for that strand's CERTIFIED-RNA measurement (a
    #: spliced fragment cannot be gDNA, so this arm needs no deconvolution).
    #: Without it a face's total is not a total, and the consequence is not subtle: mature RNA cannot
    #: cross an exon|intron boundary contiguously, so an exon and its boundary hold genuinely different
    #: fragment populations and their unspliced totals differ by the whole sj flux, with no enrichment
    #: involved. A transport that reads that difference as enrichment scales gDNA by it, which costs more
    #: than an order of magnitude even on a condition with no probes at all.
    inv_sj_lo: np.ndarray
    inv_sj_hi: np.ndarray
    sj_count_lo: np.ndarray
    sj_count_hi: np.ndarray
    #: float64[n_slots, 2] — the matching divisors, split the same way.
    eff_sj_lo: np.ndarray
    eff_sj_hi: np.ndarray
    #: float64[n_slots, 2] — the ROUTE-SUMMED certified rate per face, ``Σ_J flux_J / A_J`` over the
    #: face's junctions, by transcript strand. The junctions at one face are DISJOINT ROUTES (each
    #: molecule crosses exactly one), so their rates SUM; the pooled ``sj_count/eff_sj`` ratio is the
    #: opportunity-weighted MEAN instead and under-reads a k-route face about k-fold. Summing here, at
    #: the source, is what makes every consumer inherit the right form. A raw observation with no
    #: pseudocount: a zero-flux face reads exactly 0.
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
    """Assemble the per-slot geometry from the substrate's five populations onto the chain.

    The divisors:

    ================  ==========================================================================
    REGION, both      ``contained_eff_length(region_len, pmf)`` — no reach argument exists
    BOUNDARY, gDNA    ``UNBOUNDED_REACH`` both sides ⇒ ``mu_g − 1``. gDNA's template is the
                      chromosome, so ``taper_g = 1``: physics, not a choice
    BOUNDARY, RNA     ``UNBOUNDED_REACH`` by ruling ⇒ ``mu_r − 1``. An unspliced crossing is a
                      MIXTURE whose RNA half alone is bounded, so reach there is per COMPONENT;
                      the untapered form carries a known genome-wide gDNA over-call of about a
                      tenth. A per-boundary taper would be an A/B, never run; its switch was
                      removed (2026-09-13) rather than carried unfed
    sj, RNA           the real exonic per-strand reach. A sj is used only by a molecule that
                      spliced across it, so its divisor is the spliced one
    ================  ==========================================================================

    Where a sj attaches: its donor is the boundary to the RIGHT of ``src_region`` and its acceptor the
    boundary to the LEFT of ``dst_region``, since molecules leave the template at the first and arrive
    at the second, so both boundaries genuinely saw the flux. Both are read off the chain's own
    adjacency (``right``/``left``) rather than re-derived from per-reference offsets, which is also what
    makes cross-reference leakage impossible: a reference terminal links to ``-1``, and this function
    raises rather than attaching a sj there.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    n = chain.n_slots

    # ── the two CONTIGUOUS populations: the mixture, and the certified-RNA floor beside it ───────
    unspliced_count = np.zeros((n, 2), dtype=np.float64)
    unspliced_count[is_region] = np.asarray(substrate.region_contained.count, np.float64)[
        obj[is_region]
    ]
    unspliced_count[is_boundary] = np.asarray(substrate.boundary_unspliced.count, np.float64)[
        obj[is_boundary]
    ]
    # REGION slots stay 0, and that is structural rather than a shortcut: a contained fragment used no
    # sj, so a region's contained population cannot hold a spliced molecule.
    spliced_count = np.zeros((n, 2), dtype=np.float64)
    spliced_count[is_boundary] = np.asarray(substrate.boundary_spliced.count, np.float64)[
        obj[is_boundary]
    ]

    # ── the two per-component divisors ───────────────────────────────────────────────────────────
    region_len = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_regions = region_len.shape[0]
    unbounded = np.full(1, UNBOUNDED_REACH)

    def divisor(pmf: np.ndarray) -> np.ndarray:
        """Per-slot effective length: contained at a REGION, crossing at a BOUNDARY, the crossing at
        :data:`UNBOUNDED_REACH` on both sides for either component."""
        contained = contained_eff_length(region_len, pmf) if n_regions else np.zeros(0)
        n_boundaries = max(int(chain.n_boundaries_total), 1)
        crossing = np.full(n_boundaries, float(crossing_eff_length(pmf, unbounded, unbounded)[0]))
        out = np.zeros(n, dtype=np.float64)
        if n_regions:
            out[is_region] = contained[obj[is_region]]
        if is_boundary.any():
            out[is_boundary] = crossing[np.clip(obj[is_boundary], 0, crossing.shape[0] - 1)]
        return out

    # gDNA takes NO reach argument, ever: its template is the chromosome, so ``taper_g = 1``. That is
    # physics; the reach ruling is only about the RNA component.
    eff_gdna = divisor(gdna_fl_pmf)
    eff_rna = divisor(rna_fl_pmf)

    # ── the reciprocal-opportunity totals, straight off the banks ─────────────────────────────────
    # No divisor is applied here. A BOUNDARY slot carries ``boundary_unspliced``'s inv-length sum,
    # whose expectation IS the density (``rho * P(w >= 2) = rho``). A REGION slot carries
    # ``region_contained``'s inv-opportunity sum, whose expectation is ``rho * P(w <= ell)`` — a
    # per-component truncation, not a level, so
    # every REGION<->BOUNDARY ratio downstream carries that factor. Record the bias; do NOT divide it
    # out here — the pooled ``P_hat(w <= ell)`` resurrects none of the zero banks and re-imports the
    # per-component pmf the channel exists to avoid. Swapping in a truncation-free start bank
    # (`region_start_count / ell`) here was priced and refused: it helps pass-0 exons and regresses the
    # deliverable, the zero controls worst of all.
    inv_abundance = np.zeros(n, dtype=np.float64)
    _r_inv = getattr(substrate.region_contained, "inv_opportunity_sum", None)
    _b_inv = getattr(substrate.boundary_unspliced, "inv_length_sum", None)
    if _r_inv is not None and n_regions:
        inv_abundance[is_region] = np.asarray(_r_inv, np.float64)[obj[is_region]]
    if _b_inv is not None and is_boundary.any():
        _bi = np.asarray(_b_inv, np.float64)
        inv_abundance[is_boundary] = _bi[np.clip(obj[is_boundary], 0, _bi.shape[0] - 1)]

    # ── the JUMPING population: a sj boundary is a FACTOR on the boundaries it leaves and enters ───
    sj_count = np.zeros((n, 2), dtype=np.float64)
    eff_sj = np.zeros((n, 2), dtype=np.float64)
    #: the same flux kept apart by which genomic END of its sj the boundary is — see the dataclass.
    jc_lo = np.zeros((n, 2), dtype=np.float64)
    jc_hi = np.zeros((n, 2), dtype=np.float64)
    #: the flux's MODEL-FREE abundance per face — summed over the sj attached at that end, no divisor
    inv_sj_lo = np.zeros((n, 2), dtype=np.float64)
    inv_sj_hi = np.zeros((n, 2), dtype=np.float64)
    route_rate_lo = np.zeros((n, 2), dtype=np.float64)
    route_rate_hi = np.zeros((n, 2), dtype=np.float64)
    ej_lo = np.zeros((n, 2), dtype=np.float64)
    ej_hi = np.zeros((n, 2), dtype=np.float64)
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
        # The flux is the sj's own count summed over the GENOME-strand columns, filed under the sj's
        # TRANSCRIPT strand. That join is the whole of "sense is derived, never stored".
        flux = np.asarray(substrate.sj.count, np.float64).sum(axis=1)
        eff = crossing_eff_length(rna_fl_pmf, sj.reach_lo, sj.reach_hi)
        column = np.where(np.asarray(sj.strand) == np.int8(Strand.POS), 0, 1)
        # ``donor`` is the boundary at the sj's genomic-LOW end and ``acceptor`` the genomic-HIGH one,
        # for BOTH strands, because ``chain.right``/``chain.left`` are genomic and a sj runs
        # ``src < dst`` (`splice_graph`). The names are the index's; the meaning is genomic.
        _sj_inv = getattr(substrate.sj, "inv_length_sum", None)
        _sj_inv = None if _sj_inv is None else np.asarray(_sj_inv, np.float64)
        # One column per sj — the executable specification's shape. Asserted rather than reshaped:
        # a 2-column array here would broadcast into the face totals and silently double them.
        if _sj_inv is not None and _sj_inv.shape != (int(sj.n_sj),):
            raise ValueError(
                f"sj inv_length_sum has shape {_sj_inv.shape}, expected ({int(sj.n_sj)},) — one "
                "reciprocal-opportunity column per sj (tests/native/_accumulator_reference.py)"
            )
        live_route = eff > 0.0
        for boundary, rr in ((donor, route_rate_lo), (acceptor, route_rate_hi)):
            np.add.at(
                rr,
                (boundary[live_route], column[live_route]),
                flux[live_route] / eff[live_route],
            )
        for boundary, jc, ej, inv_face in (
            (donor, jc_lo, ej_lo, inv_sj_lo),
            (acceptor, jc_hi, ej_hi, inv_sj_hi),
        ):
            np.add.at(sj_count, (boundary, column), flux)
            np.add.at(eff_sj, (boundary, column), eff)
            np.add.at(jc, (boundary, column), flux)
            np.add.at(ej, (boundary, column), eff)
            if _sj_inv is not None:
                # filed under the sj's TRANSCRIPT strand, exactly as its count is — the flux is
                # certified RNA of that strand, so a policy can read it as a per-component measurement.
                np.add.at(inv_face, (boundary, column), _sj_inv)

    return RegionGeometry(
        n_slots=int(n),
        unspliced_count=unspliced_count,
        inv_abundance=inv_abundance,
        inv_sj_lo=inv_sj_lo,
        inv_sj_hi=inv_sj_hi,
        eff_gdna=eff_gdna,
        eff_rna=eff_rna,
        spliced_count=spliced_count,
        sj_count=sj_count,
        eff_sj=eff_sj,
        sj_count_lo=jc_lo,
        sj_count_hi=jc_hi,
        eff_sj_lo=ej_lo,
        eff_sj_hi=ej_hi,
        route_rate_lo=route_rate_lo,
        route_rate_hi=route_rate_hi,
    )


def region_gdna_geometry(geometry: RegionGeometry):
    """Per-slot gDNA support ``(unspliced count, eff_gdna)``, shared by :func:`sweep.solve_chain` and
    ``calibrate`` so every consumer divides by ONE definition.

    It returns ``eff_gdna`` and nothing else, so ``rho_g = f_g·M/E_g`` — which is exactly what makes
    ``sum_c rho_c·E_c = M`` hold and therefore what makes ``f_g`` a COUNT share rather than a density
    share. The name must keep saying "gDNA": a reader auditing whether the EM prior mixes the two should
    not have to open a second file to find out which component this divisor belongs to.

    No total density may be formed from this pair as ``mass / eff_gdna``: that is a total over one
    component's opportunity model, and :mod:`.abundance_landscape`'s measured totals are what a density
    field is fitted on instead.

    There is nothing to sum at a boundary. One set of numbers per slot means no ``mass_l + mass_r`` over
    ``E_l + E_r``, and therefore no half cancelling against a half missing from a per-face length.
    """
    return np.asarray(geometry.unspliced_count, np.float64).sum(axis=1), np.asarray(
        geometry.eff_gdna, np.float64
    )


@dataclass(frozen=True, slots=True)
class RegionBelief:
    """Per-region solved state on the chain: the composition pie `(f_pos, f_neg, f_g)` over the region's UNSPLICED
    mass + the gDNA share's posterior variance in LOG-FRACTION space, `var_gdna` =
    `Var(log f_c)`, never `Var(f_c)`. All length ``n_slots``.

    The first axis is the unified region+boundary CHAIN, not the region axis: :func:`init_beliefs`
    builds every array from ``geometry.unspliced_count``, which is ``float64[n_slots, 2]``. Sizing a new
    array off "regions" builds the wrong shape.

    The variances are log-space — grid moments of `log f_c` over the λ lattice
    (`simplex_logodds._solve_logodds`), matching the log-density message currency. They are
    therefore not bounded by ¼ and routinely exceed it; a consumer needing the linear `Var(f_c)` must
    convert (delta method `Var(f_c) ≈ f_c²·Var(log f_c)`, as `sweep.solve_chain` does for
    `composition_logvar`).

    The variance is the precision state: `Var(log f_c)=0` is locked and certain (a forbidden strand, say)
    and `=∞` is no information (unsolved). It feeds the honest message send — a source's outgoing
    precision is degraded from its own `Var_own` by the communication noise, so an unsure region speaks
    quietly. The composition is stored as a FRACTION, the face-invariant quantity, because a boundary has
    two faces but one composition; the density `ρ=f·M_face/E_face` is the message currency (computed
    inline in `sweep.solve_chain`) and the mass `m=f·M_face` (`RegionDeconv`) is the output."""

    f_pos: np.ndarray
    f_neg: np.ndarray
    f_g: np.ndarray
    var_gdna: np.ndarray
    #: Does this slot hold a COMPOSITION? — an own composition channel, structural certainty, or a
    #: composition row received from a neighbour. A level, a ceiling or a cube row is a BOUND and does
    #: not count: a node whose only evidence is a bound does not train the landscape prior, because a
    #: bound-only slot settles where the prior puts it inside the admitted half-line, so training the
    #: prior on it is training it on its own echo. Published by `sweep.solve_chain` from the held
    #: messages; ``None`` on a belief no solve has produced (`init_beliefs`), which the landscape's
    #: selector reads as the annotation alone.
    #: Gate: `tests/calibration/test_landscape_training_population.py`.
    has_composition: np.ndarray | None = None


# ---------------------------------------------------------------------------
# Initialization — the signature-binary G1/G2/G3 belief on the chain.
# ---------------------------------------------------------------------------
#
# A strand axis is hard-LOCKED (a forbidden strand, an intergenic sink) by the per-region ``allow_pos`` /
# ``allow_neg`` forbid mask in the solve. The init ALSO sets the per-component precision state ``var(f_c)``
# ``0`` = locked/certain (a forbidden strand, an intergenic gDNA sink), ``inf`` =
# no information (an admissible-but-unsolved axis — it will listen to messages, and emits none until solved).
# A solved single-strand (G2) region takes the strand-solve posterior variance.


def g1_locked(free_pos, free_neg) -> np.ndarray:
    """The G1 class — a structurally pure-gDNA slot: neither RNA strand admissible, so the composition
    is structurally certain.

    This is the predicate :func:`_type_belief` pins ``{0,0,1}`` at ``Var(log f_g) = 0`` on, and it
    applies to both axes: an intergenic region and an intergenic|exon boundary are both G1, because RNA
    cannot cross a gene boundary any more than it can occupy intergenic space.

    One definition, living here beside the code that applies it, so that every instrument that
    classifies objects by it reads one predicate instead of re-deriving its own. Two homes for one
    predicate is how a region-only variant survives in two scripts and a test at once.

    A boundary that is structurally gDNA still sits between RNA-carrying exons, so its crossing mass is
    RNA-contaminated: the predicate says the COMPOSITION is certain, not that the count is clean.
    """
    return ~np.asarray(free_pos, bool) & ~np.asarray(free_neg, bool)


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """Build the per-region composition ``(f_pos, f_neg, f_g)`` for ONE region type (regions OR boundaries) from its
    signature-binary classification + its strand-only solve.

    ``free_pos``/``free_neg`` are the per-region booleans for whether each strand's RNA axis is
    admissible (a region's own ±transcript bits; a boundary's ±strand CONTINUITY across the boundary).
    ``deconv`` is the strand-only :class:`RegionDeconv`: no global prior, no imputation. The
    signature-binary default is all-gDNA ``{0,0,1}``, and the three classes override it:

    * G1 (neither strand free — an intergenic region or a no-RNA-crossing boundary): a locked gDNA sink,
      keeping ``{0,0,1}``.
    * G2 (exactly one strand free, with data): the strand deconvolution alone resolves the pie, a
      single-strand region being 1-D with ``f_active = 1 − f_g``.
    * G3 (both strands free — AMBIG): unresolvable by strand, so the ``{0,0,1}`` default is kept at
      maximum (``inf``) variance and the sweep resolves it from neighbour messages and the prior.

    Returns the four per-region arrays ``(f_pos, f_neg, f_g, var_gdna)`` — the
    composition plus the precision state: ``var=0`` locked, ``inf`` no information, else the strand-solve
    posterior variance. The variances are ``Var(log f_c)``, log-space and unbounded above, never
    ``Var(f_c)``.
    """
    n = free_pos.shape[0]
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g = np.ones(n)  # the signature-binary all-gDNA default; the count plays no role in it
    # precision state: gDNA unsolved (inf); a strand axis is locked (0) iff forbidden, else unsolved (inf).
    var_g = np.full(n, np.inf)

    g1 = g1_locked(free_pos, free_neg)
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2-active: take the strand-only solve (median f_g, mean f±, and the posterior variances). G1 sinks
    # and G3 AMBIG slots keep the {0,0,1} default at maximum variance.
    fgv = np.asarray(deconv.gdna_frac_var, dtype=np.float64)
    f_g[g2_active] = np.asarray(deconv.gdna_frac, dtype=np.float64)[g2_active]
    f_pos[g2_active] = np.asarray(deconv.rna_pos_frac, dtype=np.float64)[g2_active]
    f_neg[g2_active] = np.asarray(deconv.rna_neg_frac, dtype=np.float64)[g2_active]
    var_g[g2_active] = fgv[g2_active]

    # G1 sink: lock the gDNA axis (the fractions are already the {0,0,1} default).
    var_g[g1] = 0.0
    return f_pos, f_neg, f_g, var_g


@dataclass(frozen=True, slots=True)
class RegionStatics:
    """Per-slot STATIC structural masks (length ``n_slots``). The sweep mutates only the dynamic
    :class:`RegionBelief`; these never change.

    No count lives here — this class is structure only. All three fragment populations sit together
    on :class:`RegionGeometry`, where their difference is visible, and a consumer slices them there. One
    quantity, one place.

    ``free_pos``/``free_neg`` are the axes on which RNA may be present at all (a region's own
    ±transcript bits; a boundary's ±continuity — the RNA-crossing gate);
    ``mrna_active_pos``/``mrna_active_neg`` are the tighter mature-RNA axes (a region's ±exon bits; a
    boundary's ±contiguous exon) that select the per-region solver prior.

    ``boundary_flags`` carries the splice graph's 8 structural bits (``TSS_s``/``TES_s``/``DONOR_s``/
    ``ACCEPTOR_s``) at each BOUNDARY slot and ``0`` on REGION slots, including when no graph was
    supplied.

    Raw bits, not pre-derived predicates. Every consumer wants a different combination of them, and a
    single pre-derived "is this a splice site" predicate was measured to be nearly the complement of
    what it was meant to replace. Compose with
    :func:`~rigel.calibration.splice_graph.is_terminus` /
    :func:`~rigel.calibration.splice_graph.is_splice_site`.
    """

    n_slots: int
    free_pos: np.ndarray  # bool — nascent-RNA-active (transcript continuity); the RNA-crossing gate
    free_neg: np.ndarray  # bool
    mrna_active_pos: (
        np.ndarray
    )  # bool — mature-RNA-active (contiguous exon); selects the region prior
    mrna_active_neg: np.ndarray  # bool
    boundary_flags: np.ndarray  # uint16 — graph structural bits; 0 on REGION slots


def build_region_statics(
    chain: RegionChain,
    region_arrays,
    boundary_flags: np.ndarray | None = None,
) -> RegionStatics:
    """Gather the structural masks onto the chain, in ONE slot-keyed pass.

    It takes no substrate: with the counts on :class:`RegionGeometry` this reads nothing but the
    signatures and the chain's own adjacency, and an unused parameter would only invite a caller to
    think the two are coupled.

    Regions and boundaries share one keyed pass rather than a twin pair of helpers, and the flank lookup
    is the chain's own adjacency. Boundary endpoints are implicit — boundary ``i`` lies between region
    ``i`` and region ``i+1``, and a boundary always has a region on both sides — so
    ``chain.left``/``chain.right`` answer it and no ``-1`` terminal branch has any cases.

    The allow mask is the transcript-structure CONTINUITY gate: a strand-``s`` unspliced crossing can be
    RNA only where strand ``s`` is present on BOTH flanks. That blocks RNA at a TSS/TES (intergenic|exon
    leaves neither strand continuous, so the slot is a gDNA sink) and at a mixed exon|AMBIG boundary.
    ``mrna_active_s`` is the tighter mature-crossing gate: contiguous exon on both flanks.

    ``boundary_flags`` is the per-contiguous-boundary ``uint16[E]`` from
    :func:`~rigel.calibration.splice_graph.build_boundary_flags_array`; ``None`` leaves the field zero,
    which every current consumer treats as "no structural information".
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
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    mrp_l, mrn_l = mrna_active_strands(sig_l)
    mrp_r, mrn_r = mrna_active_strands(sig_r)
    mr_self_p, mr_self_n = mrna_active_strands(slot_sig)

    free_pos = np.where(is_region, (ts == TS_POS) | (ts == TS_AMBIG), nrp_l & nrp_r)
    free_neg = np.where(is_region, (ts == TS_NEG) | (ts == TS_AMBIG), nrn_l & nrn_r)

    return RegionStatics(
        n_slots=n,
        # No per-region spliced floor. Spliced mature handling is owned by the message layer, so a
        # region-local floor would double-count it and would inflate a boundary's unspliced f_pos with
        # mature RNA, pushing phantom RNA into introns. Measured at least as good as keeping it in every
        # κ × capture × ±gDNA regime.
        free_pos=free_pos,
        free_neg=free_neg,
        mrna_active_pos=np.where(is_region, mr_self_p, mrp_l & mrp_r),
        mrna_active_neg=np.where(is_region, mr_self_n, mrn_l & mrn_r),
        boundary_flags=np.where(
            is_boundary, flags[np.clip(obj, 0, max(flags.shape[0] - 1, 0))], 0
        ).astype(np.uint16),
    )


def _check_boundary_flags(boundary_flags, n_boundaries: int) -> np.ndarray:
    """Validate the per-contiguous-boundary flags against the chain, BEFORE anything else is computed.

    A mis-sized array would shift every flag by one boundary — a defect invisible in aggregate and
    undetectable by a bit-identity gate while nothing reads the flags. Refuse it at the door.
    """
    if boundary_flags is None:
        return np.zeros(max(n_boundaries, 1), dtype=np.uint16)
    flags = np.asarray(boundary_flags, dtype=np.uint16)
    if flags.shape != (n_boundaries,):
        raise ValueError(
            f"boundary_flags has shape {flags.shape}; expected ({n_boundaries},), one per contiguous boundary. "
            f"Build it with splice_graph.build_boundary_flags_array(index) against the SAME index the "
            f"payload was scanned on. ⚠ There are no terminal slots: a reference with k regions owns "
            f"k-1 boundaries, not k+1."
        )
    return flags if flags.size else np.zeros(1, dtype=np.uint16)


def init_beliefs(
    chain: RegionChain,
    geometry: RegionGeometry,
    statics: RegionStatics,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
    logodds_window: float = 10.0,
) -> RegionBelief:
    """The signature-binary G1/G2/G3 initial :class:`RegionBelief` on the unified chain.

    The single-strand slots are strand-solved by the log-density log-odds solver (:mod:`simplex_logodds`,
    a one-cell tilt): the bare strand likelihood plus the Jeffreys reference, with no global prior and no
    imputation, both of which enter later in the sweep. The signature-binary class overrides
    (:func:`_type_belief`) then set the G1/G2/G3 belief. Single-strand introns resolve to ``f_g≈0`` from
    the Beta-Binomial tilt alone, which is the zero-gDNA gate; intergenic and TSS sinks lock at
    ``{0,0,1}``; AMBIG regions hold ``{0,0,1}`` at maximum variance for the sweep — so an AMBIG slot is
    not solved here at all (its cube would be discarded), which is why each slot's admissible strands are
    masked to the single-strand case before the solve."""
    st = statics
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    single = fp ^ fn
    count = np.asarray(geometry.unspliced_count, np.float64)
    deconv = _solve_regions_logodds_all(
        count[:, 0],
        count[:, 1],
        fp & single,
        fn & single,
        count.sum(axis=1),
        # the strand solve's certified-RNA floor is strand-agnostic, so the two GENOME-strand columns
        # are summed here rather than stored pre-summed — the geometry keeps the axis it was deposited on.
        np.asarray(geometry.spliced_count, np.float64).sum(axis=1),
        kappa=float(rna_sense_frac),
        od_g=gdna_strand_overdispersion,
        od_r=rna_strand_overdispersion,
        n_grid=n_grid,
        L=logodds_window,
    )
    f_pos, f_neg, f_g, var_g = _type_belief(fp, fn, deconv, count.sum(axis=1))
    return RegionBelief(f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_gdna=var_g)


# ---------------------------------------------------------------------------
# ⛔ Do not add a per-capture-class gDNA scaling here. A gDNA LEVEL crosses a composition-unlicensed hop
# UNSCALED, and the transfer policy's level lane is where that rule lives. Capture-OFF is not a separate
# case — it is the same expression.
#
# The tempting derivation is that gDNA is uniform BEFORE capture, so a gDNA density claim should
# transport between adjacent objects by the ratio of their capture efficiencies; that efficiency is
# fixed by probe geometry, probes are designed from the annotation, so it looks like a property of the
# object's structural class, directly observable on any class with structurally pure-gDNA members. It was
# built that way (REGION/BOUNDARY x off-probe / half-covered / fully-covered, a pooled ``Σcount/ΣE`` per
# class) and priced inert: a pure-gDNA object's own observed total already IS its gDNA density at its own
# capture stratum, so a pooled class ratio only re-derives locally-available information, worse.
