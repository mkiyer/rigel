"""rigel.second_pass — score the side buffer's held fragments, draw one hypothesis each, and drain them.

Pass 1 holds two kinds of fragment whole: one whose unsequenced gap has more than one surviving explanation, and a
multimapper more than one of whose placements admits a surviving (placement, hypothesis) pair. Three functions resolve
them::

    score_held_fragments  ->  choose_hypotheses  ->  drain
      the factors, apart       one draw per fragment   re-deposit the winner alone

The prior is pass 1's own tally: a hypothesis is a path, a path is a set of accumulator objects, and pass 1 counted the
molecules at each. So ``score(h) = rho(h) x f(L_h) x s(h)``, and on a multi-placement record also ``A_p``
(:func:`placement_abundance`) and ``e(p)``, the hit's sequence at the EM's prices. Every input comes from pass 1, so the
second pass runs before calibration, and calibration runs once, on the complete tally. ``L`` comes from
``Accumulator.length_under``, the C++ the deposit itself uses, so the tool has one definition of fragment length.
``drain`` is pure: it rebuilds the accumulators from the payload, so a cached scan can be drained at any seed without
re-reading the BAM.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .calibration.effective_length import observable_contained_eff_length
from .calibration.signature import RegionType, coarse_type_array
from .native import PLACE_BOUNDARY, PLACE_REGION, PLACE_SJ
from .native import Accumulator as NativeAccumulator
from .native import HeldBank as NativeHeldBank
from .scoring import sequence_log_penalties
from .scan_payload import (
    ADDITIVE_AXES,
    AccumulatorPayload,
    DrainQC,
)
from .types import Strand

if TYPE_CHECKING:
    from .calibration.region_arrays import RegionArrays


__all__ = [
    "HeldScores",
    "HypothesisTerms",
    "PlacementPath",
    "choose_hypotheses",
    "combine_factors",
    "drain",
    "score_held_fragments",
    "strand_terms",
]


@dataclass(frozen=True, slots=True)
class HypothesisTerms:
    """The factors of each hypothesis's score, kept apart so a change in the score can be traced to one of them.

    Every array is flat and aligned with ``payload.deferred.hypothesis_offsets``: one entry per hypothesis, in the
    bank's canonical order.
    """

    length: np.ndarray  # int64 — L under each hypothesis, from the C++
    density: np.ndarray  # float64 — rho, in fragments per base
    length_likelihood: np.ndarray  # float64 — f(L_h) under the pmf appropriate to the hypothesis
    strand: np.ndarray  # float64 — P(align_strand | hypothesis)
    #: float64 — A_p (:func:`placement_abundance`), carried per hypothesis; 1 on a one-placement record, where it
    #: would be a common factor and is not applied.
    placement_abundance: np.ndarray
    #: float64 — log e(p), the placement's hit sequence at the EM's prices, ``nm · m + clip · c``; 0 on a
    #: one-placement record.
    sequence: np.ndarray


@dataclass(frozen=True, slots=True)
class HeldScores:
    """Normalised scores over each held fragment's hypotheses, plus the terms that produced them."""

    score: np.ndarray  # float64, flat, sums to 1 within each record's run
    terms: HypothesisTerms
    #: Records with more than one candidate tied for the lead (:func:`combine_factors`).
    n_undecided: int


def _distinguishing_boundaries(
    region_bounds: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The local boundary range that separates ∅ from a path splicing ``[start, end)``, endpoints included.

    By the deposit rule (``tests/native/_accumulator_reference.py``) a boundary is crossed iff it lies strictly inside
    a contiguous segment. Over a fragment ``[s, e)``, ∅ crosses every boundary in ``(s, e)`` and the spliced path those
    in ``(s, start)`` and ``(end, e)``, so the two differ at the region bounds ``start <= c <= end``. Both endpoints of
    an annotated intron are region bounds, so the range is never empty.

    The endpoints must stay included: ``start < c < end`` drops the two guaranteed discriminators and returns an empty
    range whenever the intron spans one region, giving the unspliced path a ``rho`` of 0.
    """
    first = np.searchsorted(region_bounds, starts, side="left")
    last = np.searchsorted(region_bounds, ends, side="right")
    return first, last


#: P(alignment orientation | gDNA): double-stranded DNA has no sense direction, so each orientation has probability
#: one half. A biological constant, not a tunable.
P_ORIENTATION_GIVEN_GDNA = 0.5


def _strand_terms(
    align: np.ndarray, implied_strand: np.ndarray, rna_sense_frac: float
) -> tuple[np.ndarray, np.ndarray]:
    """:func:`strand_terms` over arrays, element by element."""
    align = np.asarray(align, dtype=np.int64)
    implied = np.asarray(implied_strand, dtype=np.int64)
    named = (implied == int(Strand.POS)) | (implied == int(Strand.NEG))
    spliced = np.where(
        named,
        np.where(implied == align, rna_sense_frac, 1.0 - rna_sense_frac),
        P_ORIENTATION_GIVEN_GDNA,
    )
    return spliced, np.full(spliced.shape, P_ORIENTATION_GIVEN_GDNA)


def strand_terms(*, align: int, implied_strand: int, rna_sense_frac: float) -> tuple[float, float]:
    """``(spliced, genomic)``: the likelihood of the observed alignment orientation under each hypothesis.

    With ``t`` the strand a candidate implies, ``a`` the aligned strand and ``p = P(a == t | RNA) = rna_sense_frac``
    (≈ 0.01 on a dUTP library, ≈ 0.99 on an R1-sense one)::

        H_spliced   uses a junction, so the molecule is RNA on strand t:   L = p if a == t, else 1 - p
        H_genomic   crosses the gap contiguously; its gDNA component:     L = 0.5 either way

    The genomic hypothesis's unspliced-RNA part would give the same ``p`` / ``1 - p`` as the spliced candidate and
    cancel, so only its gDNA component discriminates. That term must stay orientation-independent: a fitted
    library-wide strand marginal would move it with the spliced term and cancel most of the discrimination. An
    ambiguous or absent implied strand names no orientation to compare against.
    """
    spliced, genomic = _strand_terms(np.array([align]), np.array([implied_strand]), rna_sense_frac)
    return float(spliced[0]), float(genomic[0])


def _expand(starts: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """The flat indices ``starts[i] + k`` for ``0 ≤ k < counts[i]``, run after run."""
    counts = np.asarray(counts, dtype=np.int64)
    total = int(counts.sum())
    return np.repeat(
        np.asarray(starts, dtype=np.int64) - (np.cumsum(counts) - counts), counts
    ) + np.arange(total, dtype=np.int64)


def _segment_min(values: np.ndarray, offsets: np.ndarray, empty: float) -> np.ndarray:
    """The least of ``values[offsets[i]:offsets[i + 1]]`` per segment, *empty* where the segment is empty: a path's
    bottleneck."""
    out = np.full(offsets.shape[0] - 1, empty, dtype=np.float64)
    filled = offsets[1:] > offsets[:-1]
    if filled.any():
        out[filled] = np.minimum.reduceat(values, offsets[:-1][filled])
    return out


def _combine_runs(
    runs: np.ndarray, factors: list[tuple[np.ndarray, np.ndarray | None]]
) -> tuple[np.ndarray, np.ndarray]:
    """:func:`combine_factors` over every run of ``runs`` at once: ``factors`` is the list, in the order they are
    applied, of ``(log factor, in_run)``, where ``in_run`` marks per run whether the factor belongs to its list (``None``
    where it belongs to every run's). Returns ``(scores, undecided)``, scores flat over the runs and ``undecided`` per
    run."""
    lengths, starts = np.diff(runs), runs[:-1]
    n = int(runs[-1])
    alive = np.ones(n, dtype=bool)
    total = np.zeros(n, dtype=np.float64)
    for factor, in_run in factors:
        live = np.isfinite(factor)
        applies = np.logical_or.reduceat(
            live & alive, starts
        )  # zero for every survivor: uninformative
        if in_run is not None:
            applies &= in_run
        applied = np.repeat(applies, lengths)
        alive &= live | ~applied
        total += np.where(applied & live, factor, 0.0)
    lead = np.maximum.reduceat(np.where(alive, total, -np.inf), starts)
    weight = np.where(alive, np.exp(total - np.repeat(lead, lengths)), 0.0)
    scores = weight / np.repeat(np.add.reduceat(weight, starts), lengths)
    top = scores == np.repeat(np.maximum.reduceat(scores, starts), lengths)
    return scores, np.add.reduceat(top.astype(np.int64), starts) > 1


def combine_factors(
    density: np.ndarray,
    length_likelihood: np.ndarray,
    strand: np.ndarray,
    placement_abundance: np.ndarray | None = None,
    sequence: np.ndarray | None = None,
) -> tuple[np.ndarray, bool]:
    """Combine one fragment's factors into a normalised posterior over its candidates.

    Returns ``(scores, undecided)``: the scores sum to 1, and ``undecided`` is true when more than one candidate ties
    for the lead::

        score(h)  =  rho(h)  x  f(L_h)  x  s(h),   normalised over the candidate set

    A fragment held on its placements adds ``e(p)`` (*sequence*, in logs) and ``A_p`` (*placement_abundance*); both
    are ``None`` on a one-placement record, where they would be common factors. The product is compared only within
    one fragment, so it needs no common scale, and it is formed in logs with the maximum subtracted, so it cannot
    underflow.

    A factor that is zero for every surviving candidate is uninformative and skipped, since zero would destroy the
    product rather than cancel; a factor zero for only some candidates stays decisive, and nothing is floored. Factors
    are applied from the most evidence behind them to the least (length pmf, strand fraction, sequence price, then
    ``A_p`` and ``rho``), so a record whose traffic is zero for every candidate falls to its length term, not to a
    coin toss. One record of :func:`score_held_fragments`' combination, which runs every record at once.
    """
    n = density.shape[0]
    if n == 0:
        return np.zeros(0), False
    with np.errstate(divide="ignore"):
        factors = [(np.log(length_likelihood), None), (np.log(strand), None)]
        if sequence is not None:
            factors.append((np.asarray(sequence, dtype=np.float64), None))
        if placement_abundance is not None:
            factors.append((np.log(placement_abundance), None))
        factors.append((np.log(density), None))
    scores, undecided = _combine_runs(np.array([0, n], dtype=np.int64), factors)
    return scores, bool(undecided[0])


def resolve_intron_lookups(payload, accumulators) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Answer, one reference at a time, the two questions the scoring loop asks of each hypothesis intron.

    Which annotated junction it is, if any (the kernel's own ``Accumulator.sj_edge_ids``), and the local boundary
    range that separates the genomic path from a path splicing it (:func:`_distinguishing_boundaries`). Batched per
    reference rather than per intron, which on a deep library is millions of numpy calls.

    Returns three arrays over the flat hypothesis-intron axis (``hypothesis_intron_offsets`` indexes it): the
    payload's sj id or -1, and the range's two ends, local to the reference.
    """
    d = payload.deferred
    n_intron = int(d.hypothesis_introns.shape[0]) // 2
    sj_of = np.full(n_intron, -1, dtype=np.int64)
    first_of = np.zeros(n_intron, dtype=np.int64)
    last_of = np.zeros(n_intron, dtype=np.int64)
    if n_intron == 0:
        return sj_of, first_of, last_of

    starts = np.ascontiguousarray(d.hypothesis_introns[0::2], dtype=np.int64)
    ends = np.ascontiguousarray(d.hypothesis_introns[1::2], dtype=np.int64)
    # the hypothesis each intron belongs to, and the placement each hypothesis belongs to
    hyp_of = np.repeat(
        np.arange(d.n_hypotheses, dtype=np.int64), np.diff(d.hypothesis_intron_offsets)
    )
    placement_of = np.repeat(
        np.arange(d.n_placements, dtype=np.int64), np.diff(d.hypothesis_offsets)
    )[hyp_of]
    # the motif the lookup filters on: the observed one when the aligner wrote it, else the implied strand
    observed = np.asarray(d.sj_strand, dtype=np.int64)[placement_of]
    implied = np.asarray(d.hypothesis_sj_strand, dtype=np.int64)[hyp_of]
    motif = np.where(observed != int(Strand.NONE), observed, implied).astype(np.int32)

    refs = np.asarray(d.ref, dtype=np.int64)[placement_of]
    bounds = np.asarray(payload.region_bounds, dtype=np.int64)
    for ref in np.unique(refs):
        here = refs == ref
        ref = int(ref)
        lo = int(payload.ref_region_bound_offsets[ref])
        hi = int(payload.ref_region_bound_offsets[ref + 1])
        s_here, e_here = starts[here], ends[here]
        first_of[here], last_of[here] = _distinguishing_boundaries(bounds[lo:hi], s_here, e_here)
        # the accumulator answers in ref-local slots, whose base `_Accumulators` checked against `ref_sj_offsets`
        local = np.asarray(
            accumulators[ref].sj_edge_ids(
                starts=s_here, ends=e_here, sj_strand=np.ascontiguousarray(motif[here])
            ),
            dtype=np.int64,
        )
        base = int(payload.ref_sj_offsets[ref])
        sj_of[here] = np.where(local < 0, -1, local + base)
    return sj_of, first_of, last_of


class _Accumulators:
    """A lazy ``ref -> Accumulator`` map built from the payload's own region-bound axis, junctions installed.

    This lets the second pass score and drain a cached scan without re-reading the BAM. ``region_types`` comes from
    the index: the payload carries the region bounds but not their typing, and the typing assigns a deposit to a
    length pool.

    ``set_sj`` must always be called. Without a junction table every observed intron reads as unannotated, the
    junction banks stay zero, and the tally still looks well-formed.

    :attr:`built` holds the references actually built, which is the drain delta's support, so :func:`_gather_delta`
    iterates it instead of every reference.
    """

    def __init__(self, payload: AccumulatorPayload, region_types: np.ndarray, sj) -> None:
        self._payload = payload
        self._region_types = region_types
        self._sj = sj
        self.built: dict[int, NativeAccumulator] = {}

    def __getitem__(self, ref: int) -> NativeAccumulator:
        if ref in self.built:
            return self.built[ref]
        payload, sj = self._payload, self._sj
        region_bound_lo = int(payload.ref_region_bound_offsets[ref])
        region_bound_hi = int(payload.ref_region_bound_offsets[ref + 1])
        region_lo = int(payload.ref_region_offsets[ref])
        n_regions = max(region_bound_hi - region_bound_lo - 1, 0)
        accumulator = NativeAccumulator(
            region_bounds=np.ascontiguousarray(
                payload.region_bounds[region_bound_lo:region_bound_hi], dtype=np.int64
            ),
            region_types=np.ascontiguousarray(
                self._region_types[region_lo : region_lo + n_regions], dtype=np.uint8
            ),
            max_length=payload.max_length,
            ref=ref,
        )
        # The junction CSR sliced to this reference. A junction id is its rank, so the slot base must equal the
        # payload's `ref_sj_offsets[ref]`; a mismatch would write one reference's traffic onto another's slots.
        empty = region_bound_hi <= region_bound_lo
        slot_lo = 0 if empty else int(sj.offsets[region_bound_lo])
        slot_hi = 0 if empty else int(sj.offsets[region_bound_hi])
        if slot_lo != int(payload.ref_sj_offsets[ref]):
            raise ValueError(
                f"reference {ref}'s sj start at slot {slot_lo} in the index but the payload's "
                f"ref_sj_offsets says {int(payload.ref_sj_offsets[ref])}; a sj id IS its rank, so "
                f"the two axes must agree."
            )
        accumulator.set_sj(
            offsets=np.zeros(1, np.int32)
            if empty
            else np.ascontiguousarray(
                sj.offsets[region_bound_lo : region_bound_hi + 1] - slot_lo, dtype=np.int32
            ),
            boundary_right=np.ascontiguousarray(
                sj.boundary_right[slot_lo:slot_hi] - region_bound_lo, dtype=np.int32
            ),
            sj_strand=np.ascontiguousarray(sj.strand[slot_lo:slot_hi], dtype=np.int8),
        )
        self.built[ref] = accumulator
        return accumulator


@dataclass(frozen=True, slots=True)
class PlacementPath:
    """The accumulator objects a placement's observed path deposits on, as flat ids on the payload's axes, read by
    the deposit's own rule (``Accumulator.placement_masses``): the annotated junctions its observed introns resolve
    to (``sj``), the contiguous boundaries it crosses (``boundaries``), and the one region that contains it
    (``regions``, empty unless the path neither crosses nor splices)."""

    sj: np.ndarray
    boundaries: np.ndarray
    regions: np.ndarray

    @property
    def spliced(self) -> bool:
        """The deposit books the placement spliced: certified RNA under every one of its paths."""
        return bool(self.sj.size)


@dataclass(frozen=True, slots=True)
class _Paths:
    """Every held placement's :class:`PlacementPath` at once: per placement, its junctions' count, their
    ``sj_inv_length_sum`` bottleneck and the strands of the transcripts holding them (OR-ed); its boundaries' count
    and their ``boundary_unspliced_inv_length_sum`` bottleneck; and its containing region (-1 for none)."""

    n_sj: np.ndarray
    sj_traffic: np.ndarray
    sj_strand: np.ndarray
    n_boundaries: np.ndarray
    boundary_traffic: np.ndarray
    region: np.ndarray

    @classmethod
    def from_objects(
        cls,
        payload: AccumulatorPayload,
        offsets: np.ndarray,
        kinds: np.ndarray,
        ids: np.ndarray,
        junction_strands: np.ndarray | None = None,
    ) -> "_Paths":
        """From the placements' objects, CSR over placements by ``offsets``: each object's kind
        (``rigel.native.PLACE_*``) and flat id, the objects of one kind of a placement contiguous. Without
        ``junction_strands`` no strand is read (every ``sj_strand`` 0)."""
        n = offsets.shape[0] - 1
        owner = np.repeat(np.arange(n, dtype=np.int64), np.diff(offsets))

        def of_kind(kind):
            mine = kinds == kind
            count = np.bincount(owner[mine], minlength=n).astype(np.int64)
            return ids[mine], np.r_[0, np.cumsum(count)], count

        sj, sj_offsets, n_sj = of_kind(PLACE_SJ)
        boundaries, boundary_offsets, n_boundaries = of_kind(PLACE_BOUNDARY)
        regions, _region_offsets, n_regions = of_kind(PLACE_REGION)
        if int(n_regions.max(initial=0)) > 1:
            raise ValueError("a placement's path is contained in at most one region")
        sj_strand = np.zeros(n, dtype=np.int64)
        has_sj = n_sj > 0
        if junction_strands is not None and has_sj.any():
            sj_strand[has_sj] = np.bitwise_or.reduceat(
                junction_strands[sj], sj_offsets[:-1][has_sj]
            )
        region = np.full(n, -1, dtype=np.int64)
        region[n_regions > 0] = regions
        return cls(
            n_sj=n_sj,
            sj_traffic=_segment_min(
                np.asarray(payload.sj_inv_length_sum, dtype=np.float64)[sj], sj_offsets, 0.0
            ),
            sj_strand=sj_strand,
            n_boundaries=n_boundaries,
            boundary_traffic=_segment_min(
                np.asarray(payload.boundary_unspliced_inv_length_sum, dtype=np.float64)[boundaries],
                boundary_offsets,
                0.0,
            ),
            region=region,
        )

    @property
    def spliced(self) -> np.ndarray:
        return self.n_sj > 0


def pass_one_bulk(
    payload: AccumulatorPayload, region_types: np.ndarray, region_support: np.ndarray
) -> float:
    """gDNA's bulk density on pass one: unique contained fragments per observable start over the intergenic regions,
    the unit :func:`placement_abundance` reads a contained placement in. 0 when there is no intergenic support."""
    intergenic = region_types == RegionType.INTERGENIC
    support = float(region_support[intergenic].sum())
    if support <= 0.0:
        return 0.0
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)[intergenic]
    return float(counts.sum()) / support


def _abundances(
    payload: AccumulatorPayload, paths: _Paths, region_support: np.ndarray, bulk: float
) -> np.ndarray:
    """:func:`placement_abundance` for every placement of *paths*."""
    region = paths.region
    contained = region >= 0
    traffic = np.zeros(region.shape[0], dtype=np.float64)
    if contained.any():
        count = np.asarray(payload.region_contained_count, dtype=np.float64)[region[contained]].sum(
            axis=1
        )
        opportunity = np.asarray(region_support, dtype=np.float64)[region[contained]]
        with np.errstate(divide="ignore", invalid="ignore"):
            traffic[contained] = np.where(opportunity > 0.0, count / opportunity, 0.0)
    crossing = paths.n_boundaries > 0
    traffic[crossing] = paths.boundary_traffic[crossing]
    # gDNA's floor, where gDNA can be: the larger of the two, the traffic on a tie
    floored = np.where(bulk > traffic, bulk, traffic)
    return np.where(paths.spliced, paths.sj_traffic, floored)


def placement_abundance(
    payload: AccumulatorPayload,
    path: PlacementPath,
    region_support: np.ndarray,
    bulk: float,
) -> float:
    """``A_p``: the pass-one unique traffic at the objects the placement's own path deposits on, held at least at
    gDNA's *bulk* wherever gDNA can be. By the deposit rule the path (:class:`PlacementPath`) deposits on:

    * the annotated junctions it uses: the bottleneck of their ``sj_inv_length_sum``;
    * else the contiguous boundaries it crosses: the bottleneck of their ``boundary_unspliced_inv_length_sum``;
    * else the one region that contains it: its ``region_contained_count`` per observable start. *region_support* is
      each region's contained support over its observable pieces under the unconditional anchor
      (``effective_length.observable_contained_eff_length``), since no fragment is observed on a masked base;
    * else nothing (an unannotated intron swallowed every boundary between two regions), which reads 0.

    Every form is a start density per base, so one bound serves all three.

    gDNA's floor: pass one counts only what it deposited, so a copy whose every fragment multimaps reads no traffic
    however much gDNA it holds. gDNA is one template at a uniform rate and capture only enriches, so a placement not
    booked spliced reads ``max(A_p, bulk)``: a max, because ``A_p`` already holds gDNA's share where pass one saw it.
    """
    kinds = np.concatenate(
        [
            np.full(path.sj.shape[0], PLACE_SJ, dtype=np.int8),
            np.full(path.boundaries.shape[0], PLACE_BOUNDARY, dtype=np.int8),
            np.full(path.regions.shape[0], PLACE_REGION, dtype=np.int8),
        ]
    )
    ids = np.concatenate([path.sj, path.boundaries, path.regions]).astype(np.int64)
    paths = _Paths.from_objects(payload, np.array([0, ids.shape[0]]), kinds, ids)
    return float(_abundances(payload, paths, region_support, bulk)[0])


def _held_bank(deferred) -> NativeHeldBank:
    """The held bank as the accumulator's batch reads take it, over the bank's own arrays."""
    return NativeHeldBank(
        start=deferred.start,
        end=deferred.end,
        align_strand=deferred.align_strand,
        sj_strand=deferred.sj_strand,
        observed_intron_offsets=deferred.observed_intron_offsets,
        observed_introns=deferred.observed_introns,
        hypothesis_offsets=deferred.hypothesis_offsets,
        hypothesis_sj_strand=deferred.hypothesis_sj_strand,
        hypothesis_intron_offsets=deferred.hypothesis_intron_offsets,
        hypothesis_introns=deferred.hypothesis_introns,
    )


@dataclass(frozen=True, slots=True)
class _HeldScorer:
    """What scoring the held fragments reads beyond the fragments themselves: the length laws, the strand fraction,
    the pass-one tallies' derived terms and the per-reference accumulators (:func:`score_held_fragments`)."""

    payload: AccumulatorPayload
    rna_pmf: np.ndarray
    global_pmf: np.ndarray
    max_size: int
    rna_sense_frac: float
    #: A_p's divisor per region: the unconditional anchor over the region's observable pieces
    region_support: np.ndarray
    bulk: float
    #: the EM's sequence prices at the library's mismatch rate
    mismatch_price: float
    clip_price: float
    junction_strands: np.ndarray
    accumulators: _Accumulators
    #: over the flat hypothesis-intron axis: each intron's sj id or -1, and its distinguishing boundary range
    sj_of: np.ndarray
    first_of: np.ndarray
    last_of: np.ndarray
    intron_offsets: np.ndarray

    @classmethod
    def build(
        cls, payload: AccumulatorPayload, fl_models, rna_sense_frac: float, regions, sj
    ) -> "_HeldScorer":
        n_regions = int(np.asarray(payload.region_contained_count).shape[0])
        if regions.n_regions != n_regions:
            raise ValueError(
                f"score_held_fragments: the frame has {regions.n_regions:,} regions and the payload {n_regions:,} — "
                "they address different partitions. Pass RegionArrays.from_index of the index the scan used."
            )
        region_types = coarse_type_array(regions.signature)
        region_support = observable_contained_eff_length(regions, fl_models.global_pmf)
        mismatch_price, clip_price = sequence_log_penalties(payload.mismatch_rate)
        accumulators = _Accumulators(payload, region_types, sj)
        sj_of, first_of, last_of = resolve_intron_lookups(payload, accumulators)
        return cls(
            payload=payload,
            rna_pmf=fl_models.rna_pmf,
            global_pmf=fl_models.global_pmf,
            max_size=int(fl_models.max_size),
            rna_sense_frac=rna_sense_frac,
            region_support=region_support,
            bulk=pass_one_bulk(payload, region_types, region_support),
            mismatch_price=mismatch_price,
            clip_price=clip_price,
            junction_strands=np.asarray(sj.strand, dtype=np.int64),
            accumulators=accumulators,
            sj_of=sj_of,
            first_of=first_of,
            last_of=last_of,
            intron_offsets=np.asarray(payload.deferred.hypothesis_intron_offsets, dtype=np.int64),
        )

    def reads(self) -> tuple[np.ndarray, _Paths]:
        """Every hypothesis's ``L`` (``Accumulator.length_under``) and every placement's :class:`PlacementPath`
        (``Accumulator.placement_masses``), the deposit's own reads, one batch per reference."""
        payload = self.payload
        d = payload.deferred
        bank = _held_bank(d)
        length = np.zeros(d.n_hypotheses, dtype=np.int64)
        n_objects = np.zeros(d.n_placements, dtype=np.int64)
        batches = []
        for ref in np.unique(d.ref):
            ref = int(ref)
            placements = np.flatnonzero(d.ref == ref)
            lengths, offsets, kinds, ids = self.accumulators[ref].held_reads(placements, bank)
            length[
                _expand(d.hypothesis_offsets[placements], np.diff(d.hypothesis_offsets)[placements])
            ] = lengths
            n_objects[placements] = np.diff(offsets)
            kinds = np.asarray(kinds)
            base = np.select(
                [kinds == PLACE_SJ, kinds == PLACE_BOUNDARY],
                [int(payload.ref_sj_offsets[ref]), int(payload.ref_boundary_offsets[ref])],
                int(payload.ref_region_offsets[ref]),
            )
            batches.append((placements, kinds, np.asarray(ids, dtype=np.int64) + base))
        offsets = np.r_[0, np.cumsum(n_objects)]
        kinds = np.empty(int(offsets[-1]), dtype=np.int8)
        ids = np.empty(int(offsets[-1]), dtype=np.int64)
        for placements, batch_kinds, batch_ids in batches:
            at = _expand(offsets[placements], n_objects[placements])
            kinds[at], ids[at] = batch_kinds, batch_ids
        return length, _Paths.from_objects(payload, offsets, kinds, ids, self.junction_strands)

    def terms(self) -> HypothesisTerms:
        """Every hypothesis's factors; the placement's own terms (``A_p`` and the sequence price) only on a
        multi-placement record."""
        payload = self.payload
        d = payload.deferred
        length, paths = self.reads()
        per_placement = np.diff(d.hypothesis_offsets)
        placement = np.repeat(np.arange(d.n_placements, dtype=np.int64), per_placement)
        implies = np.diff(self.intron_offsets) > 0
        spliced = paths.spliced[placement]
        # the length pmf: RNA's where the path is certified spliced, the unconditional anchor's otherwise
        inside = (length >= 0) & (length <= self.max_size)
        at = np.where(inside, length, 0)
        length_likelihood = np.where(
            inside,
            np.where(implies | spliced, self.rna_pmf[at], self.global_pmf[at]),
            0.0,
        )
        # the strand: the observed motif's, else the junctions' the path uses, else the hypothesis's own
        motif = d.sj_strand[placement]
        observed = motif != int(Strand.NONE)
        reading = np.where(
            observed, motif, np.where(spliced, paths.sj_strand[placement], d.hypothesis_sj_strand)
        )
        spliced_term, genomic_term = _strand_terms(
            d.align_strand[placement], reading, self.rna_sense_frac
        )
        strand = np.where(observed | spliced | implies, spliced_term, genomic_term)
        # the placement's own terms, on a multi-placement record
        per_record = np.diff(d.placement_offsets)
        multi = np.repeat(per_record > 1, per_record)
        nm, clip = d.nm.astype(np.float64), d.clip.astype(np.float64)
        with np.errstate(invalid="ignore"):  # a price of -inf at no edit is never read
            sequence = np.where(d.nm != 0, nm * self.mismatch_price, 0.0) + np.where(
                d.clip != 0, clip * self.clip_price, 0.0
            )
        abundance = _abundances(payload, paths, self.region_support, self.bulk)
        return HypothesisTerms(
            length=length,
            density=self._density(implies, placement),
            length_likelihood=length_likelihood,
            strand=strand,
            placement_abundance=np.where(multi, abundance, 1.0)[placement],
            sequence=np.where(multi, sequence, 0.0)[placement],
        )

    def _density(self, implies: np.ndarray, placement: np.ndarray) -> np.ndarray:
        """``rho``: a spliced path's evidence is its junctions' ``sj_inv_length_sum``, deposited by the same rule as a
        contiguous boundary's, so it is on ∅'s scale; the genomic path's is the unspliced crossing density where the
        others jump. ∅ is scored only over the union of the competing implied introns, the stretch it claims
        contiguous and every spliced path jumps, so it is not penalised for touching more objects than they do; that
        union is the placement's contiguous run of the flat hypothesis-intron axis. Each is the bottleneck of its
        objects, 0 where it has none."""
        payload, d, offsets = self.payload, self.payload.deferred, self.intron_offsets
        sj_traffic = np.asarray(payload.sj_inv_length_sum, dtype=np.float64)
        annotated = self.sj_of >= 0
        junction = np.zeros(self.sj_of.shape[0], dtype=np.float64)
        junction[annotated] = sj_traffic[self.sj_of[annotated]]
        spliced = _segment_min(junction, offsets, 0.0)
        # ∅: every boundary in each intron's distinguishing range, over the placement's run of introns
        run = offsets[d.hypothesis_offsets]
        n_boundaries = np.maximum(self.last_of - self.first_of, 0)
        boundary = _expand(self.first_of, n_boundaries)
        intron = np.repeat(np.arange(n_boundaries.shape[0], dtype=np.int64), n_boundaries)
        intron_placement = np.repeat(np.arange(d.n_placements, dtype=np.int64), np.diff(run))
        base = np.asarray(payload.ref_boundary_offsets, dtype=np.int64)[d.ref][intron_placement]
        crossing = np.asarray(payload.boundary_unspliced_inv_length_sum, dtype=np.float64)[
            base[intron] + boundary - 1
        ]
        by_placement = np.r_[
            0, np.cumsum(np.bincount(intron_placement[intron], minlength=d.n_placements))
        ]
        genomic = _segment_min(crossing, by_placement, 0.0)
        return np.where(implies, spliced, genomic[placement])


def score_held_fragments(
    payload: AccumulatorPayload,
    *,
    fl_models,
    rna_sense_frac: float,
    regions: "RegionArrays",
    sj,
) -> HeldScores:
    """Score every held fragment's hypotheses. Pure: reads pass-1 state, writes nothing.

    Parameters
    ----------
    payload
        Pass 1's tally and its side buffer.
    fl_models
        ``pipeline.library_fl_models(payload, index)``. ``rna_pmf`` scores a spliced hypothesis (certified RNA);
        ``global_pmf``, the unconditional anchor, scores the genomic one, whose component is unknown.
    rna_sense_frac
        ``P(align_strand == transcript strand | RNA)``, from ``StrandModels.exonic_spliced``. Near 0 on a dUTP
        library; used as the probability of agreement, never inverted.
    regions, sj
        From the index, since neither is on the payload: ``RegionArrays.from_index`` (the scanner's partition with
        the unobservable mask, giving region types and observable support) and ``build_sj_arrays``.
    """
    deferred = payload.deferred
    terms = _HeldScorer.build(payload, fl_models, rna_sense_frac, regions, sj).terms()
    # Combined over each record's whole run, every placement's hypotheses at once.
    per_record = np.diff(deferred.placement_offsets)
    multi = per_record > 1
    with np.errstate(divide="ignore"):
        factors = [
            (np.log(terms.length_likelihood), None),
            (np.log(terms.strand), None),
            (terms.sequence, multi),
            (np.log(terms.placement_abundance), multi),
            (np.log(terms.density), None),
        ]
    score, undecided = _combine_runs(deferred.record_hypothesis_offsets, factors)
    return HeldScores(score=score, terms=terms, n_undecided=int(undecided.sum()))


def choose_hypotheses(scores: HeldScores, payload: AccumulatorPayload, *, seed: int) -> np.ndarray:
    """One multinomial draw per held fragment; returns a local hypothesis index per record.

    One hypothesis wins the whole fragment, so counts stay integers. The bank's canonical order is identical at any
    worker count, and one ``rng.random(n)`` call draws once per record in that order, so the draw is reproducible.

    The draw must never be keyed on the fragment's content: identical fragments would then draw identically, and a
    60/40 posterior over a hundred duplicates would collapse to 100/0.
    """
    offsets = payload.deferred.record_hypothesis_offsets
    starts, ends = offsets[:-1], offsets[1:]
    n = payload.deferred.n_fragments
    if n == 0:
        return np.zeros(0, np.int64)

    # Scores are normalised within each run, so a run's cumulative sum spans [base, base + 1].
    cumulative = np.cumsum(scores.score)
    base = cumulative[starts] - scores.score[starts]
    draws = np.random.default_rng(seed).random(n)
    picked = np.searchsorted(cumulative, base + draws, side="right")
    # Clamp onto the run's last positive slot: a uniform near 1 can overshoot the run's final cumulative value by one
    # rounding and land on the next record's hypothesis, or on a trailing hypothesis the score excluded.
    record = np.repeat(np.arange(n, dtype=np.int64), ends - starts)
    slot = np.arange(scores.score.shape[0], dtype=np.int64)
    positive = scores.score > 0.0
    last = np.full(n, -1, dtype=np.int64)
    np.maximum.at(last, record[positive], slot[positive])
    return np.clip(picked, starts, last) - starts


def drain(
    payload: AccumulatorPayload,
    choices: np.ndarray,
    *,
    region_types: np.ndarray,
    sj,
) -> AccumulatorPayload:
    """Replay every held fragment with its chosen hypothesis; return the drained payload.

    ``choices[i]`` is a local hypothesis index for the ``i``-th record in the bank's canonical order, as
    :func:`choose_hypotheses` produces. Each fragment re-enters ``Accumulator.deposit`` with that hypothesis alone, so
    arbitration is degenerate and the ordinary deposit rules apply: there is one tally path.

    Pure: the delta accumulates into fresh per-reference accumulators and is added to pass one's arrays, so both
    tallies exist afterwards. The bank is already merged and ordered, so the drain is single-threaded and reproducible:
    each reference's records are replayed in the bank's order.
    """
    deferred = payload.deferred
    n = deferred.n_fragments
    choices = np.ascontiguousarray(choices, dtype=np.int64)
    if choices.shape != (n,):
        raise ValueError(
            f"drain got {choices.shape[0]} choices for {n} held fragments. One choice per held record, "
            f"in the bank's canonical order — a mismatch means the scores and the draw walked different "
            f"queues."
        )
    # The choice indexes the record's run of hypotheses over all its placements, so it names a placement and a
    # hypothesis at once; the replay goes onto that placement's reference.
    record_runs = deferred.record_hypothesis_offsets
    chosen = record_runs[:-1] + choices
    outside = np.flatnonzero((choices < 0) | (chosen >= record_runs[1:]))
    if outside.size:
        i = int(outside[0])
        raise ValueError(f"choice {choices[i]} is outside record {i}'s hypothesis set")
    placement = np.searchsorted(deferred.hypothesis_offsets, chosen, side="right") - 1
    accumulators = _Accumulators(payload, region_types, sj)
    bank = _held_bank(deferred)
    counts = {
        "deposited": 0,
        "dropped_too_long": 0,
        "dropped_empty": 0,
        "dropped_strand_undefined": 0,
    }
    ref = deferred.ref[placement]
    for r in np.unique(ref):
        mine = ref == r
        for outcome, k in (
            accumulators[int(r)].deposit_held(placement[mine], chosen[mine], bank).items()
        ):
            if outcome == "deferred_undetermined_gap":
                raise AssertionError(
                    "a hypothesis set of size one deferred, which is unreachable: arbitration defers only "
                    "when two or more hypotheses survive."
                )
            counts[outcome] += k
    chose_no_implied_intron = int((np.diff(deferred.hypothesis_intron_offsets)[chosen] == 0).sum())
    delta = _gather_delta(payload, accumulators)
    return payload.with_drain(
        delta,
        DrainQC(
            offered=n,
            offered_multimapper=int((np.diff(deferred.placement_offsets) > 1).sum()),
            **counts,
            chose_no_implied_intron=chose_no_implied_intron,
            chose_implied_intron=n - chose_no_implied_intron,
            census_before=payload.gap_resolution,
        ),
    )


def _gather_delta(
    payload: AccumulatorPayload, accumulators: _Accumulators
) -> dict[str, np.ndarray]:
    """The drained accumulators' channels, placed into globally shaped arrays. Only the references the drain touched
    were built; the rest stay zero."""
    axis_offsets = {
        "region": payload.ref_region_offsets,
        "boundary": payload.ref_boundary_offsets,
        "sj": payload.ref_sj_offsets,
    }
    delta = {name: np.zeros_like(getattr(payload, name)) for name, _axis in ADDITIVE_AXES}
    for ref, accumulator in accumulators.built.items():
        for name, axis in ADDITIVE_AXES:
            block = np.asarray(getattr(accumulator, name))
            if axis == "library":
                # The pool histograms and the unconditional anchor are library-wide, so every
                # reference's contribution lands on the same rows.
                delta[name] += block.reshape(delta[name].shape)
            else:
                lo = int(axis_offsets[axis][ref])
                delta[name][lo : lo + block.shape[0]] += block
    return delta
