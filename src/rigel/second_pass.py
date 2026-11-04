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
from .calibration.signature import coarse_type_array
from .native import PLACE_BOUNDARY, PLACE_REGION, PLACE_SJ
from .native import Accumulator as NativeAccumulator
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


def strand_terms(*, align: int, implied_strand: int, rna_sense_frac: float) -> tuple[float, float]:
    """``(spliced, genomic)``: the likelihood of the observed alignment orientation under each hypothesis.

    With ``t`` the strand a candidate implies, ``a`` the aligned strand and ``p = P(a == t | RNA) = rna_sense_frac``
    (≈ 0.01 on a dUTP library, ≈ 0.99 on an R1-sense one)::

        H_spliced   uses a junction, so the molecule is RNA on strand t:   L = p if a == t, else 1 - p
        H_genomic   crosses the gap contiguously; its gDNA component:     L = 0.5 either way

    The genomic hypothesis's unspliced-RNA part would give the same ``p`` / ``1 - p`` as the spliced candidate and
    cancel, so only its gDNA component discriminates. That term must stay orientation-independent: a fitted
    library-wide strand marginal would move it with the spliced term and cancel most of the discrimination.
    """
    if implied_strand in (int(Strand.POS), int(Strand.NEG)):
        spliced = rna_sense_frac if implied_strand == align else 1.0 - rna_sense_frac
    else:
        # An ambiguous or absent implied strand names no orientation to compare against.
        spliced = P_ORIENTATION_GIVEN_GDNA
    return spliced, P_ORIENTATION_GIVEN_GDNA


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
    coin toss.
    """
    n = density.shape[0]
    if n == 0:
        return np.zeros(0), False
    with np.errstate(divide="ignore"):
        logs = [np.log(length_likelihood), np.log(strand)]
        if sequence is not None:
            logs.append(np.asarray(sequence, dtype=np.float64))
        if placement_abundance is not None:
            logs.append(np.log(placement_abundance))
        logs.append(np.log(density))
    alive = np.ones(n, dtype=bool)
    total = np.zeros(n, dtype=np.float64)
    for factor in logs:
        live = np.isfinite(factor)
        if not live[alive].any():
            continue  # zero for every survivor: uninformative
        alive &= live
        total += np.where(live, factor, 0.0)
    weight = np.where(alive, np.exp(total - total[alive].max()), 0.0)
    scores = weight / weight.sum()
    # Exact equality: identical inputs give bit-identical scores, and a 0.50001 / 0.49999 split is decided.
    return scores, bool(int((scores == scores.max()).sum()) > 1)


def _bottleneck(values: list[float]) -> float:
    """The density of a path that needs several objects: the scarcest one's, since a molecule on the path was present
    at every object on it."""
    return float(min(values)) if values else 0.0


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


def _placement_path(
    payload: AccumulatorPayload, p: int, accumulator: NativeAccumulator
) -> PlacementPath:
    """Placement ``p``'s :class:`PlacementPath`, from the accumulator of its reference."""
    d = payload.deferred
    ref = int(d.ref[p])
    kinds, ids, _masses = accumulator.placement_masses(
        start=int(d.start[p]),
        end=int(d.end[p]),
        observed_introns=[tuple(q) for q in d.observed_introns_of(p).tolist()],
        align_strand=int(d.align_strand[p]),
        sj_strand=int(d.sj_strand[p]),
    )
    kinds, ids = np.asarray(kinds), np.asarray(ids, dtype=np.int64)
    return PlacementPath(
        sj=ids[kinds == PLACE_SJ] + int(payload.ref_sj_offsets[ref]),
        boundaries=ids[kinds == PLACE_BOUNDARY] + int(payload.ref_boundary_offsets[ref]),
        regions=ids[kinds == PLACE_REGION] + int(payload.ref_region_offsets[ref]),
    )


def pass_one_bulk(
    payload: AccumulatorPayload, region_types: np.ndarray, region_support: np.ndarray
) -> float:
    """gDNA's bulk density on pass one: unique contained fragments per observable start over the intergenic regions,
    the unit :func:`placement_abundance` reads a contained placement in. 0 when there is no intergenic support."""
    intergenic = region_types == 0
    support = float(region_support[intergenic].sum())
    if support <= 0.0:
        return 0.0
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)[intergenic]
    return float(counts.sum()) / support


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
    if path.spliced:
        return _bottleneck(payload.sj_inv_length_sum[path.sj].tolist())
    if path.boundaries.size:
        traffic = _bottleneck(payload.boundary_unspliced_inv_length_sum[path.boundaries].tolist())
    elif path.regions.size:
        (region,) = path.regions.tolist()
        count = float(np.asarray(payload.region_contained_count[region], dtype=np.float64).sum())
        opportunity = float(region_support[region])
        traffic = count / opportunity if opportunity > 0.0 else 0.0
    else:
        traffic = 0.0
    return max(traffic, bulk)


@dataclass(frozen=True, slots=True)
class _HeldScorer:
    """What scoring a held fragment reads beyond the fragment itself: the length laws, the strand fraction, the
    pass-one tallies' derived terms and the per-reference accumulators (:func:`score_held_fragments`)."""

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

    def score_placement(self, p: int, multi: bool, terms: HypothesisTerms) -> None:
        """Write placement ``p``'s factors into *terms*, one slot per hypothesis; the placement's own terms (``A_p``
        and the sequence price) only on a multi-placement record."""
        d = self.payload.deferred
        accumulator = self.accumulators[int(d.ref[p])]
        h0, h1 = int(d.hypothesis_offsets[p]), int(d.hypothesis_offsets[p + 1])
        path = _placement_path(self.payload, p, accumulator)
        if multi:
            terms.placement_abundance[h0:h1] = placement_abundance(
                self.payload, path, self.region_support, self.bulk
            )
            nm, clip = int(d.nm[p]), int(d.clip[p])
            terms.sequence[h0:h1] = (nm * self.mismatch_price if nm else 0.0) + (
                clip * self.clip_price if clip else 0.0
            )
        hypotheses = [
            (
                [tuple(q) for q in d.hypothesis_introns_of(h).tolist()],
                int(d.hypothesis_sj_strand[h]),
            )
            for h in range(h0, h1)
        ]
        # `length_under` needs the hypothesis objects back in the shape the binding reads.
        spans = [_Span(introns, implied) for introns, implied in hypotheses]
        observed = [tuple(q) for q in d.observed_introns_of(p).tolist()]
        align, observed_motif = int(d.align_strand[p]), int(d.sj_strand[p])
        for local, (introns, implied_strand) in enumerate(hypotheses):
            slot = h0 + local
            L = int(
                accumulator.length_under(
                    start=int(d.start[p]),
                    end=int(d.end[p]),
                    observed_introns=observed,
                    align_strand=align,
                    sj_strand=observed_motif,
                    hypotheses=spans,
                    hypothesis_index=local,
                )
            )
            terms.length[slot] = L
            terms.density[slot] = self._density(p, slot, implies_introns=bool(introns))
            pmf = self.rna_pmf if (introns or path.spliced) else self.global_pmf
            terms.length_likelihood[slot] = float(pmf[L]) if 0 <= L <= self.max_size else 0.0
            terms.strand[slot] = self._strand(
                align, observed_motif, path, implied_strand, implies_introns=bool(introns)
            )

    def _density(self, p: int, slot: int, *, implies_introns: bool) -> float:
        """``rho``: a spliced path's evidence is its junctions' ``sj_inv_length_sum``, deposited by the same rule as a
        contiguous boundary's, so it is on ∅'s scale; the genomic path's is the unspliced crossing density where the
        others jump. ∅ is scored only over the union of the competing implied introns, the stretch it claims
        contiguous and every spliced path jumps, so it is not penalised for touching more objects than they do; that
        union is the placement's contiguous run of the flat hypothesis-intron axis."""
        payload, offsets = self.payload, self.intron_offsets
        if implies_introns:
            observed_densities = []
            for t in range(int(offsets[slot]), int(offsets[slot + 1])):
                jid = int(self.sj_of[t])
                observed_densities.append(0.0 if jid < 0 else float(payload.sj_inv_length_sum[jid]))
            return _bottleneck(observed_densities)
        d = payload.deferred
        h0, h1 = int(d.hypothesis_offsets[p]), int(d.hypothesis_offsets[p + 1])
        boundary_base = int(payload.ref_boundary_offsets[int(d.ref[p])])
        boundary_densities = []
        for t in range(int(offsets[h0]), int(offsets[h1])):
            for boundary in range(int(self.first_of[t]), int(self.last_of[t])):
                boundary_densities.append(
                    float(payload.boundary_unspliced_inv_length_sum[boundary_base + boundary - 1])
                )
        return _bottleneck(boundary_densities)

    def _strand(
        self,
        align: int,
        observed_motif: int,
        path: PlacementPath,
        implied_strand: int,
        *,
        implies_introns: bool,
    ) -> float:
        """The probability of the observed orientation under the placement's own reading. An observed motif makes
        every hypothesis RNA on the motif's strand, since its intron is cut under all of them: constant within the
        placement, not across a multimapper's, where an unspliced hit carries gDNA's ½. An annotated splice the
        aligner wrote no motif for is RNA on the strand of the transcripts holding its junctions, as the deposit books
        it; both strands say nothing."""
        if observed_motif != int(Strand.NONE):
            reading = observed_motif
        elif path.spliced:
            reading = int(np.bitwise_or.reduce(self.junction_strands[path.sj]))
        else:
            spliced_term, genomic_term = strand_terms(
                align=align, implied_strand=implied_strand, rna_sense_frac=self.rna_sense_frac
            )
            return spliced_term if implies_introns else genomic_term
        term, _ = strand_terms(
            align=align, implied_strand=reading, rna_sense_frac=self.rna_sense_frac
        )
        return term


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
    n_hyp = deferred.n_hypotheses
    terms = HypothesisTerms(
        length=np.zeros(n_hyp, np.int64),
        density=np.zeros(n_hyp, np.float64),
        length_likelihood=np.zeros(n_hyp, np.float64),
        strand=np.ones(n_hyp, np.float64),
        placement_abundance=np.ones(n_hyp, np.float64),
        sequence=np.zeros(n_hyp, np.float64),
    )
    scorer = _HeldScorer.build(payload, fl_models, rna_sense_frac, regions, sj)
    score = np.zeros(n_hyp, np.float64)
    record_runs = deferred.record_hypothesis_offsets
    n_undecided = 0
    for i in range(deferred.n_fragments):
        placements = deferred.placements_of(i)
        multi = len(placements) > 1
        for p in placements:
            scorer.score_placement(p, multi, terms)
        # Combined over the record's whole run, every placement's hypotheses at once.
        r0, r1 = int(record_runs[i]), int(record_runs[i + 1])
        score[r0:r1], undecided = combine_factors(
            terms.density[r0:r1],
            terms.length_likelihood[r0:r1],
            terms.strand[r0:r1],
            terms.placement_abundance[r0:r1] if multi else None,
            terms.sequence[r0:r1] if multi else None,
        )
        n_undecided += undecided
    return HeldScores(score=score, terms=terms, n_undecided=n_undecided)


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
    tallies exist afterwards. The bank is already merged and ordered, so the drain is single-threaded and reproducible.
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
    accumulators = _Accumulators(payload, region_types, sj)
    counts = {
        "deposited": 0,
        "dropped_too_long": 0,
        "dropped_empty": 0,
        "dropped_strand_undefined": 0,
    }
    chose_no_implied_intron = 0
    offered_multimapper = 0
    record_runs = deferred.record_hypothesis_offsets

    for i in range(n):
        # The choice indexes the record's run of hypotheses over all its placements, so it names a placement and a
        # hypothesis at once; the replay goes onto that placement's reference.
        h = int(record_runs[i]) + int(choices[i])
        if not (int(record_runs[i]) <= h < int(record_runs[i + 1])):
            raise ValueError(f"choice {choices[i]} is outside record {i}'s hypothesis set")
        p = deferred.placement_of_hypothesis(h)
        offered_multimapper += int(len(deferred.placements_of(i)) > 1)
        introns = [tuple(pair) for pair in deferred.hypothesis_introns_of(h).tolist()]
        chose_no_implied_intron += not introns
        outcome = accumulators[int(deferred.ref[p])].deposit(
            start=int(deferred.start[p]),
            end=int(deferred.end[p]),
            observed_introns=[tuple(q) for q in deferred.observed_introns_of(p).tolist()],
            align_strand=int(deferred.align_strand[p]),
            sj_strand=int(deferred.sj_strand[p]),
            hypotheses=[_Span(introns, int(deferred.hypothesis_sj_strand[h]))],
        )
        if outcome == "deferred_undetermined_gap":
            raise AssertionError(
                "a hypothesis set of size one deferred, which is unreachable: arbitration defers only "
                "when two or more hypotheses survive."
            )
        counts[outcome] += 1

    delta = _gather_delta(payload, accumulators)
    return payload.with_drain(
        delta,
        DrainQC(
            offered=n,
            offered_multimapper=offered_multimapper,
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


@dataclass(frozen=True, slots=True)
class _Span:
    """A hypothesis in the shape ``Accumulator.length_under`` and ``deposit`` read.

    The binding reads ``introns`` and ``sj_strand`` off whatever it is handed (the specification's ``GapHypothesis``
    names), so a test can pass the reference's own objects.
    """

    introns: list
    sj_strand: int
