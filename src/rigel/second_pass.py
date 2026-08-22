"""rigel.second_pass — scoring the side buffer's held fragments, and draining it.

     §3 the score (P2), §5 the draw and §6 the drain (P3)

Pass 1 holds every fragment whose unsequenced gap has more than one surviving explanation — 2–3.5 % of a
library, and systematically the **long** ones, because a longer gap admits more hypotheses. This module
resolves them, in three steps that are three functions::

    score_held_fragments  ->  choose_hypotheses  ->  drain
        the three factors      one multinomial       re-enter deposit with the
        of §3, kept apart      draw per fragment     winner ALONE, and consume the bank

⭐ **``drain`` is PURE: payload in, payload out.** It rebuilds one accumulator per reference from the
payload's own region_bound axis, so the whole second pass runs off a *cached* scan — which is what lets one scan be
drained repeatedly at different seeds without re-reading the BAM, and what and P6
both need. It also means the drain never sees a thread, so reproducibility is structural.

⛔ **There is exactly ONE tally path.** The drain re-enters ``Accumulator.deposit`` with a hypothesis set
of size one, so arbitration is degenerate and every crossing rule, quantum, pool and ``L`` is the code that
ran in pass one. `tests/native/test_accumulator_drain.py` checks the consequence directly: a drained
fragment gives byte-identically the tally that offering only that hypothesis would have given.

⭐ **THE PRIOR IS THE ACCUMULATOR ITSELF.** A hypothesis is a path, and a path is a set of accumulator
objects; pass 1 already counted how many molecules used each of them. So the score needs no transcript
abundance, no calibration output, and no second pass over the BAM::

    score(h)  =  rho(h)  x  f(L_h)  x  s(h)

⚠ **Every input comes from pass 1**, which is what lets this run BEFORE calibration and lets calibration
then run exactly once, on the complete tally.

⛔ **``L`` IS NOT COMPUTED HERE.** It comes from ``Accumulator.length_under`` — the same C++ that will
compute it again at drain time. The tool has ONE definition of
fragment length, and a scorer with its own would be a second definition of exactly the quantity that audit
unified.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._bam_impl import Accumulator as NativeAccumulator
from .scan_payload import (
    _DEFERRED_RECORD_FIELDS,
    ADDITIVE_AXES,
    AccumulatorPayload,
    DrainQC,
)
from .types import Strand


__all__ = [
    "HeldScores",
    "HypothesisTerms",
    "choose_hypotheses",
    "combine_factors",
    "drain",
    "lift_choices",
    "score_held_fragments",
    "strand_terms",
]


@dataclass(frozen=True, slots=True)
class HypothesisTerms:
    """The three factors, kept apart. ⭐ Diagnostics, not bookkeeping.

    the umbrella census partitions the held set by *which question is open*,
    which is the same thing as *which factor discriminates*. Reporting the terms separately is what lets a
    regression be attributed to one of them instead of to "the score moved".

    Every array is flat and aligned with ``payload.deferred.hypothesis_offsets`` — one entry per
    hypothesis, in the bank's own canonical order.
    """

    length: np.ndarray  # int64 — L under each hypothesis, from the C++
    density: np.ndarray  # float64 — rho, in fragments per base
    length_likelihood: np.ndarray  # float64 — f(L_h) under the pmf appropriate to the hypothesis
    strand: np.ndarray  # float64 — P(align_strand | hypothesis)


@dataclass(frozen=True, slots=True)
class HeldScores:
    """Normalised scores over each held fragment's hypotheses, plus the terms that produced them."""

    score: np.ndarray  # float64, flat, sums to 1 within each record's run
    terms: HypothesisTerms
    #: ⚠ Records whose hypotheses all scored zero. Their score run is left UNIFORM rather than NaN —
    #: a fragment the evidence cannot separate is a fragment to pick from at random, not one to drop.
    #: ⭐ The count is the honest denominator for "how much did the evidence actually decide?" and P2's
    #: measurement of D-3 reads it directly.
    n_undecided: int

    @property
    def n_hypotheses(self) -> int:
        return int(self.score.shape[0])


def _exact_region_bound(region_bounds: np.ndarray, lo: int, hi: int, position: int) -> int:
    """Flat region_bound index of ``position`` within ``region_bounds[lo:hi]``, or -1 if it is not a region_bound there."""
    k = lo + int(np.searchsorted(region_bounds[lo:hi], position))
    return k if k < hi and int(region_bounds[k]) == position else -1


def _sj_id(sj, region_bounds, lo, hi, start, end, sj_strand) -> int:
    """The annotated sj slot for one intron, or -1.

    ⚠ Mirrors ``Accumulator::sj_edge_id`` — same CSR, same strand rule (the filter applies only when the
    motif strand is DEFINITE; a non-definite one matches on coordinates alone). ⭐ Not a duplicate of the
    length definition: this is a lookup into the index's own table, and the slot IS the payload's sj
    axis index.
    """
    donor = _exact_region_bound(region_bounds, lo, hi, start)
    if donor < 0:
        return -1
    acceptor = _exact_region_bound(region_bounds, lo, hi, end)
    if acceptor < 0:
        return -1
    definite = sj_strand in (int(Strand.POS), int(Strand.NEG))
    for k in range(int(sj.offsets[donor]), int(sj.offsets[donor + 1])):
        if int(sj.boundary_right[k]) != acceptor:
            continue
        if definite and int(sj.strand[k]) != sj_strand:
            continue
        return k
    return -1


def _distinguishing_boundaries(
    region_bounds: np.ndarray, lo: int, hi: int, start: int, end: int
) -> tuple[int, int]:
    """The LOCAL boundary range that separates ∅ from a path splicing ``[start, end)``. Endpoints INCLUDED.

    ⭐ **Read off the deposit rule, not chosen** (`_accumulator_reference.py`, the per-segment crossing
    loop): a boundary is crossed iff it lies **strictly inside a contiguous segment**. So over a fragment
    ``[s, e)``:

    * ∅ is one segment and crosses every boundary in ``(s, e)``;
    * the spliced path is ``[s, start)`` and ``[end, e)`` and crosses those in ``(s, start)`` and
      ``(end, e)``.

    The difference is the boundaries at region_bounds ``start <= c <= end`` — and since both endpoints of an annotated
    intron are cut by construction, those two boundaries always exist and always discriminate.

    ⛔ **This asked for ``start < c < end`` until 2026-08-02 (D-6).** That drops exactly the two
    guaranteed discriminators, and returns an EMPTY range whenever the intron spans one region — handing ∅
    a structural ``rho`` of 0. Measured over the pilot: on a gdna100 library the shipped rule read
    **43.4 %** of ∅ hypotheses as zero-density where the derived rule reads ⭐ **0.0000**.
    """
    first = int(np.searchsorted(region_bounds[lo:hi], start, side="left"))
    last = int(np.searchsorted(region_bounds[lo:hi], end, side="right"))
    return first, last


#: ⭐ **P(alignment orientation | gDNA) — biologically fixed, not estimated.** Double-stranded DNA has no
#: sense direction, so either orientation is equally likely. ⚠ Not a tunable and not a fitted marginal: it
#: is the same value `StrandModels` already records for gDNA everywhere else — *"gDNA is scored with a fixed
#: strand probability of 0.5 (no strand bias), not learned from intergenic data"*.
P_ORIENTATION_GIVEN_GDNA = 0.5


def strand_terms(*, align: int, implied_strand: int, rna_sense_frac: float) -> tuple[float, float]:
    """``(spliced, genomic)`` — the likelihood of the OBSERVED ALIGNMENT ORIENTATION under each hypothesis.

    Let ``t`` be the strand a candidate implies, ``a`` the strand the fragment aligned to, and
    ``p = P(a == t | RNA) = rna_sense_frac`` the library's directional sense fraction (≈ 0.01 on a dUTP
    library, ≈ 0.99 on an R1-sense one; the direction-agnostic specificity is ``max(p, 1 - p)``)::

        H_spliced   used a sj, and gDNA cannot splice, so the molecule is RNA on strand t
                    L = p        if a == t
                      = 1 - p    otherwise
        H_genomic   crossed the gap contiguously; the component that DISCRIMINATES is gDNA
                    L = 0.5      either orientation, because DNA is double-stranded

    ⭐ **Why ``0.5`` and not a fitted mixture.** ∅ also covers *unspliced RNA* — but unspliced RNA would
    give the same ``p`` / ``1 - p`` as the spliced candidate and therefore cancel, contributing nothing. The
    only part of ∅ that can separate it from a spliced path is its gDNA component, and that component's
    orientation likelihood is a biological constant.

    ⛔ **A GLOBAL MIXTURE MARGINAL IS WORSE THAN USELESS HERE, and the algebra says so.** With any
    *constant* ``c`` for ∅ the orientation discrimination is::

        [c / p] / [c / (1 - p)]  =  (1 - p) / p  =  98.0        <- c cancels

    so ``0.5`` and ``1.0`` discriminate identically; ``0.5`` is simply the value that makes it a
    probability, since the two orientations must sum to 1 under gDNA while ``1.0`` twice sums to 2.
    ⚠ But an orientation-*dependent* ∅ term ``q / (1 - q)`` taken from the library-wide genic marginal
    (`StrandModels.exonic`, measured at ``q = 0.1825`` on a gdna100 pilot) gives
    ``q(1-p) / p(1-q) = 21.9`` — it moves ∅ in the SAME direction as the spliced term and destroys **78 %**
    of the signal. Measured, it also cost accuracy. ⛔ Do not reintroduce a fitted marginal here: a global
    value says nothing about an individual fragment, whose gene may be silent (pure gDNA) or highly
    expressed.
    """
    if implied_strand in (int(Strand.POS), int(Strand.NEG)):
        spliced = rna_sense_frac if implied_strand == align else 1.0 - rna_sense_frac
    else:
        # ⚠ An AMBIGUOUS or absent implied strand (D-5: one path claimed by both strands) names no
        # orientation to compare against, so it says nothing either way.
        spliced = P_ORIENTATION_GIVEN_GDNA
    return spliced, P_ORIENTATION_GIVEN_GDNA


def combine_factors(
    density: np.ndarray, length_likelihood: np.ndarray, strand: np.ndarray
) -> tuple[np.ndarray, bool]:
    """Combine one fragment's three factors into a normalised posterior over its candidates.

    Returns ``(scores, undecided)``. The scores sum to 1 across the candidate set; ``undecided`` says the
    evidence separated nothing, so the run was left uniform.

        score(h)  =  rho(h)  x  f(L_h)  x  s(h),   normalised over the candidate set

    ⭐ **THE FACTORS ARE ONLY EVER COMPARED WITHIN ONE FRAGMENT.** The product is not a rate and not a
    calibrated likelihood — ``rho`` carries units of fragments per base while ``f`` and ``s`` are
    probabilities — and it does not need to be, because the normalisation makes any common scale cancel.

    ⛔ **AN ALL-ZERO FACTOR IS UNINFORMATIVE, NOT DECISIVE** (fixed 2026-08-03). That normalisation is
    exactly why: a factor taking the *same* value for every candidate cancels and cannot affect the answer.
    Zero is the one value where the arithmetic loses that property — instead of cancelling it destroys the
    product, the normalisation cannot run, and the record collapses to a coin toss that discards the other
    two factors. So a factor that is zero for every candidate is dropped for that fragment.

    ⚠ Measured before the fix, on 171,534 pilot fragments scored against the simulator's per-fragment
    truth: of the 10 records that fell back to uniform, the **length term alone would have picked the
    correct candidate 8 times out of the 8 it could decide** — 100 % — and the coin got them right by
    chance instead. Both of P4's impossible-length fragments came from those records, where the length term
    had already scored the impossible answer at zero.

    ⛔ **THE PARTIAL-ZERO CASE IS UNTOUCHED, AND THAT IS THE POINT.** A factor that is zero for *some*
    candidates and positive for others is highly informative — the zero says "no evidence for this path" —
    and it stays decisive. That is the owner's D-3 ruling ("no fallback, a hard zero stays hard") left
    exactly as it was. ⭐ No constant is introduced here; nothing is floored, softened or smoothed.

    ⭐ **AND "UNINFORMATIVE" IS JUDGED AMONG THE SURVIVORS, NOT GLOBALLY.** Factors are applied in order of
    how much evidence stands behind them — the length pmf is fitted from millions of fragments, the strand
    fraction from the library's whole spliced population, a single object's traffic from however many
    fragments happened to touch it — and each is skipped when it is flat-zero across the candidates *still
    standing*. ⚠ Judging it globally is not enough, and the pilot showed why: on record 155262 the length
    term left two candidates possible at ``9.9e-06`` and ``1.4e-03``, traffic was zero for both, and a rule
    that only checked the whole set fell through to a **fair** coin — discarding a 143-fold difference in
    likelihood. Weighting by the surviving factor is the same statement as dropping an uninformative one,
    applied one level down.
    """
    n = density.shape[0]
    if n == 0:
        return np.zeros(0), False
    # ⚠ Ordered by weight of evidence, strongest first: an ``f`` of 0 rests on the whole library's length
    # distribution, an ``s`` of 0 on its whole spliced population, a ``rho`` of 0 on one object's traffic.
    alive = np.ones(n, dtype=bool)
    scores = np.ones(n, dtype=np.float64)
    for factor in (length_likelihood, strand, density):
        if not (factor[alive] > 0.0).any():
            continue  # flat zero among the survivors: uninformative here, so it says nothing
        alive &= factor > 0.0
        scores = scores * factor
    scores = np.where(alive, scores, 0.0)
    total = float(scores.sum())
    scores = scores / total if total > 0.0 else np.full(n, 1.0 / n)
    # ⭐ `undecided` says what it means: more than one candidate is tied for the lead, so the draw between
    # them is a coin toss. ⚠ Stated as a property of the ANSWER rather than of which factors contributed,
    # because that is what the diagnostics want and it cannot drift from the arithmetic above. Exact
    # equality is right: identical inputs give bit-identical outputs, and a 0.50001/0.49999 split is
    # decided, barely.
    return scores, bool(int((scores == scores.max()).sum()) > 1)


def _bottleneck(values: list[float]) -> float:
    """Combine the densities of the several objects one path needs. ⭐ The scarcest one bounds the path.

    A molecule that took this path was present at **every** object on it, so the path's rate is not the
    sum — it is bounded by the object that saw the least traffic.

    ✅ **D-1 IS CLOSED (P2.2, 2026-08-02) and the parameter is GONE.** ``geometric`` was implemented
    alongside this so the choice could be measured rather than taken on taste. Measured over the 8 pilot
    conditions, the two pick a **different winner on 0.47–0.59 % of held records** and leave the zero mask
    **bit-identical** (both return 0 whenever any input is 0, so D-1 never interacted with D-3). No
    accuracy argument separates them at that scale and P2 has no truth criterion that could, so the
    derivation decides: this one has an argument behind it and ``geometric`` had only a robustness
    intuition. ⚠ Shares do move — over 0.10 on 1.9–3.4 % of records — so if P4's drained-tail gate ever
    misses against truth, this is a documented place to look.
    """
    return float(min(values)) if values else 0.0


class _Accumulators:
    """A lazy ``ref -> Accumulator`` map built from the PAYLOAD's own region_bound axis, sj installed.

    ⭐ **This is what lets the second pass run off a payload rather than off a live scanner.** The
    accumulator is described by its region_bound positions, and the payload carries them — so a cached scan can be
    scored and drained without re-reading the BAM, which is what keeps `build_scan_cache.py` useful for
    every re-measurement and P6 ask for.

    ⛔ **`set_sj` is not optional here and the failure is invisible.** With no sj table every
    observed intron reads as unannotated, so the sj banks stay at zero and the tally still looks
    well-formed — the C++ refuses the same omission on the scan path for exactly this reason. The scorer
    only needs ``length_under`` and would survive without it; the drain would silently credit no sj
    boundary at all.

    ⚠ ``region_types`` comes from the INDEX. The payload echoes the region_bound axis but not the typing, and the
    typing is what assigns a deposit to a length pool.

    ⚠ Lazy, and :attr:`built` is the set that actually exists — which is exactly the delta's support, so
    :func:`_gather_delta` iterates it instead of every reference in the genome.
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
        # ⭐ The sj CSR sliced to THIS reference. A sj id IS its rank, so the slot base must
        # agree with the payload's own `ref_sj_offsets[ref]` — checked rather than assumed, because a
        # silent disagreement would write one reference's sj traffic onto another's slots.
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


def score_held_fragments(
    payload: AccumulatorPayload,
    *,
    fl_models,
    rna_sense_frac: float,
    region_types: np.ndarray,
    sj,
) -> HeldScores:
    """Score every held fragment's hypotheses. Pure: reads pass-1 state, writes nothing.

    Parameters
    ----------
    payload
        Pass 1's tally **and** its side buffer.
    fl_models
        ``build_fl_models(payload)``. ⭐ ``rna_pmf`` scores a spliced hypothesis, which is certified RNA;
        ``global_pmf`` — the unconditional anchor — scores the genomic one, whose component is unknown and
        must therefore be marginalised over the library's own composition.
    rna_sense_frac
        ``P(align_strand == the transcript's strand | RNA)`` — ``StrandModels.exonic_spliced``, which is
        certified RNA because an annotated splice proves it. ⚠ On an R1-antisense (dUTP) library this is
        ≈ 0.01, so **disagreement is the likely case**. It is used as the probability
        of agreement it is, never inverted.
    region_types, sj
        From the index: ``build_region_partition_arrays`` and ``build_sj_arrays``. ⚠ The sj
        CSR is genuinely not on the payload, so this dependency is real rather than an oversight.
    """
    deferred = payload.deferred
    n_hyp = deferred.n_hypotheses
    length = np.zeros(n_hyp, np.int64)
    density = np.zeros(n_hyp, np.float64)
    length_likelihood = np.zeros(n_hyp, np.float64)
    strand = np.ones(n_hyp, np.float64)
    score = np.zeros(n_hyp, np.float64)

    region_bounds = payload.region_bounds
    rna_pmf, global_pmf = fl_models.rna_pmf, fl_models.global_pmf
    max_size = int(fl_models.max_size)

    accumulators = _Accumulators(payload, region_types, sj)

    n_undecided = 0
    for i in range(deferred.n_fragments):
        ref = int(deferred.ref[i])
        start, end = int(deferred.start[i]), int(deferred.end[i])
        align = int(deferred.align_strand[i])
        observed_motif = int(deferred.sj_strand[i])
        observed = [tuple(p) for p in deferred.observed_introns_of(i).tolist()]
        h0, h1 = int(deferred.hypothesis_offsets[i]), int(deferred.hypothesis_offsets[i + 1])
        hypotheses = [
            (
                [tuple(p) for p in deferred.hypothesis_introns_of(h).tolist()],
                int(deferred.hypothesis_sj_strand[h]),
            )
            for h in range(h0, h1)
        ]
        acc = accumulators[ref]
        region_bound_lo, region_bound_hi = (
            int(payload.ref_region_bound_offsets[ref]),
            int(payload.ref_region_bound_offsets[ref + 1]),
        )
        boundary_base = int(payload.ref_boundary_offsets[ref])
        # ⭐ The region the GENOMIC hypothesis claims is contiguous and every spliced one jumps: the union
        # of the competing implied introns. Scoring `∅` over exactly this — rather than over its whole
        # path — is what keeps the comparison symmetric, since otherwise `∅` is penalised simply for
        # touching more objects than a path that jumps them.
        contested = [intron for introns, _ in hypotheses for intron in introns]

        # `length_under` needs the hypothesis objects back in the shape the binding reads.
        spans = [_Span(introns, sj) for introns, sj in hypotheses]

        for local, (introns, implied_strand) in enumerate(hypotheses):
            slot = h0 + local
            L = int(
                acc.length_under(
                    start=start,
                    end=end,
                    observed_introns=observed,
                    align_strand=align,
                    sj_strand=observed_motif,
                    hypotheses=spans,
                    hypothesis_index=local,
                )
            )
            length[slot] = L

            # -- rho ------------------------------------------------------------------------------
            if introns:
                # A spliced path's evidence is the sj it uses. ⚠ `sj_inv_length_sum` is deposited
                # by the SAME rule as a contiguous boundary, so the two are the same quantity on the same
                # scale — that is what makes this comparable to `∅`'s number at all.
                motif = observed_motif if observed_motif != int(Strand.NONE) else implied_strand
                observed_densities = []
                for a, b in introns:
                    jid = _sj_id(sj, region_bounds, region_bound_lo, region_bound_hi, a, b, motif)
                    observed_densities.append(
                        0.0 if jid < 0 else float(payload.sj_inv_length_sum[jid])
                    )
                density[slot] = _bottleneck(observed_densities)
            else:
                # The genomic path's evidence is the unspliced crossing density where the others jump.
                boundary_densities = []
                for a, b in contested:
                    first, last = _distinguishing_boundaries(
                        region_bounds, region_bound_lo, region_bound_hi, a, b
                    )
                    for boundary in range(first, last):
                        boundary_densities.append(
                            float(
                                payload.boundary_unspliced_inv_length_sum[
                                    boundary_base + boundary - 1
                                ]
                            )
                        )
                density[slot] = _bottleneck(boundary_densities)

            # -- f(L) -----------------------------------------------------------------------------
            pmf = rna_pmf if introns else global_pmf
            length_likelihood[slot] = float(pmf[L]) if 0 <= L <= max_size else 0.0

            # -- strand ---------------------------------------------------------------------------
            # ⭐ An OBSERVED motif pins the fragment's strand, and pass 1 has already constrained the
            # hypotheses to it, so the term is constant across the set and cancels.
            if observed_motif == int(Strand.NONE):
                spliced_term, genomic_term = strand_terms(
                    align=align,
                    implied_strand=implied_strand,
                    rna_sense_frac=rna_sense_frac,
                )
                strand[slot] = spliced_term if introns else genomic_term

        # ⭐ The three factors are combined in ONE place, so the rule that an all-zero factor is
        # uninformative cannot be stated differently anywhere else. See :func:`combine_factors`.
        score[h0:h1], undecided = combine_factors(
            density[h0:h1], length_likelihood[h0:h1], strand[h0:h1]
        )
        n_undecided += undecided

    return HeldScores(
        score=score,
        terms=HypothesisTerms(
            length=length,
            density=density,
            length_likelihood=length_likelihood,
            strand=strand,
        ),
        n_undecided=n_undecided,
    )


def lift_choices(whole: AccumulatorPayload, parts, choices: np.ndarray):
    """⭐⭐ Carry hypothesis choices made on the WHOLE library onto its origin PARTITIONS.

    ⛔ **WHY THIS EXISTS, AND IT IS AN IDENTITY NOT A CONVENIENCE.** Splitting a BAM by fragment origin
    and re-scanning reconstructs the pass-one payload exactly, because pass one deposits each fragment
    independently — that identity is what makes an origin-split oracle a valid truth source. The second
    pass breaks it: its multinomial is scored against the *whole* payload's densities, so draining three
    partitions separately is a different operation from draining the whole and
    ``Sum(partitions) != whole``.

    ⭐ The constraint is on **WHERE the choice is made, not on whether it can be made at all.** Score and
    draw ONCE on the whole library, then replay each fragment's already-chosen hypothesis inside whichever
    partition holds it. Every fragment then deposits independently *given its choice*, the identity is
    restored, and each partition's drained tally is comparable with the whole's::

        scores  = score_held_fragments(whole, ...)
        choices = choose_hypotheses(scores, whole, seed=...)      # <- ONCE, on the whole
        drained = drain(whole, choices, ...)                      # the tally the solver reads
        lifted, ambiguous = lift_choices(whole, partitions, choices)     # the truth, per origin
        for p, ch in zip(partitions, lifted):
            drain(p, ch, ...)

    ⛔⛔ **EVERY PARTITION IS LIFTED IN ONE CALL, AND THAT IS THE WHOLE POINT OF THE SIGNATURE.** A record's
    choice is consumed from a per-key queue, so the queue's state has to be shared across the partitions or
    two of them can take the SAME entry and leave another unused — ``Sum(partitions) != whole`` again, by a
    different route. Taking a sequence makes the identity a property of this function rather than of a
    caller's discipline. ⚠ A one-partition caller passes ``[p]``.
    ⭐ *This was found by perturbation, not by design*: the first version took a single partition and a
    per-call queue, and the gate set could not see the defect (TRAPS: perturb-every-gate).

    ⭐ **The key is the bank's own canonical sort key** — ``_DEFERRED_RECORD_FIELDS``, the tuple the C++
    sorts on before the bank crosses the ABI, imported rather than restated so there is one definition of
    record identity (TRAPS: a-test-that-redefines). :class:`DeferredFragments` guarantees the property this rests on:
    *"two records that tie on that key are identical records, so no tie-break is needed or possible."*
    Identical records have identical hypothesis SETS — enumeration reads the span and the annotation, never
    the origin — so a LOCAL hypothesis index transfers between them unchanged.

    ⚠ **The one ambiguity, counted rather than hidden.** When several whole-library records share a key and
    are split across origins, no partition can know which of them it holds; the assignment is greedy in
    canonical order. The deposits are interchangeable (identical records, identical hypothesis sets) but the
    ORIGIN attribution is not — swapping a choice between a gDNA and an RNA fragment of the same span would
    credit one origin's mass to the other. ⛔ So the count is RETURNED, never swallowed: a caller must
    report it, and it bounds the truth error exactly. On distinct-span substrates it is 0 and the lift is
    exact.

    Returns ``(lifted, n_ambiguous)`` — ``lifted`` is one ``int64`` array per partition, each aligned with
    that partition's ``deferred.hypothesis_offsets``.

    :raises ValueError: if a partition holds a record ``whole`` does not, or holds more copies of one key
        than ``whole`` does. Either means the payloads are not a partition of one scan, and nothing
        downstream is meaningful.
    """
    w = whole.deferred
    parts = list(parts)
    choices = np.asarray(choices, np.int64)
    if choices.shape[0] != w.n_fragments:
        raise ValueError(
            f"lift_choices got {choices.shape[0]} choices for {w.n_fragments} whole-library held "
            "records. One choice per held record, in the bank's canonical order."
        )

    def _keys(d):
        return list(
            zip(*(np.asarray(getattr(d, f), np.int64).tolist() for f in _DEFERRED_RECORD_FIELDS))
        )

    # key -> the whole's choices for that key, in canonical order. Consumed left to right ACROSS all
    # partitions, so each entry is handed out exactly once.
    pool: dict[tuple, list[int]] = {}
    for i, k in enumerate(_keys(w)):
        pool.setdefault(k, []).append(int(choices[i]))
    taken: dict[tuple, int] = {}

    out, ambiguous = [], 0
    for p in parts:
        d = p.deferred
        lifted = np.zeros(d.n_fragments, np.int64)
        for j, k in enumerate(_keys(d)):
            run = pool.get(k)
            if run is None:
                raise ValueError(
                    f"lift_choices: a partition holds record {j} with key {k}, which is not in the whole "
                    "library's held bank. The payloads are not a partition of one scan — check they came "
                    "from the same BAM."
                )
            at = taken.get(k, 0)
            if at >= len(run):
                raise ValueError(
                    f"lift_choices: the partitions hold more records with key {k} ({at + 1}) than the "
                    f"whole library does ({len(run)}). They are not a partition of one scan."
                )
            if len(run) > 1:
                ambiguous += 1
            lifted[j] = run[at]
            taken[k] = at + 1
        out.append(lifted)
    return out, ambiguous


def choose_hypotheses(scores: HeldScores, payload: AccumulatorPayload, *, seed: int) -> np.ndarray:
    """⭐ One multinomial draw per held fragment. Returns a LOCAL hypothesis index per record.

    one hypothesis wins the whole fragment. No fractional deposit — integers
    stay integers and no transcript is fractionated.

    ⭐ **Seeding: order by identity, then ONE stream** — S2.1's rule, applied verbatim. The deferred bank
    carries a canonical sort (S1) that is bit-identical at any worker count, so "the i-th record" is
    well defined and a single stream consumed in that order is reproducible by construction.

    ⛔ **Never key the draw on the fragment's CONTENT.** A content hash ties on exactly the duplicates it
    would harm: 100 identical fragments would draw identically, and a 60/40 posterior would collapse to
    100/0 instead of splitting. S2.1 established the rule and the reason.

    ⚠ Vectorised rather than a loop, and that is not only for speed: one ``rng.random(n)`` call is one
    draw per record in queue order, which makes the correspondence between stream position and queue index
    a property of the code rather than of a loop body that could accidentally consume twice.
    """
    offsets = payload.deferred.hypothesis_offsets
    starts, ends = offsets[:-1], offsets[1:]
    n = payload.deferred.n_fragments
    if n == 0:
        return np.zeros(0, np.int64)

    # Scores are normalised within each run, so a run's cumulative sum spans [base, base + 1].
    cumulative = np.cumsum(scores.score)
    base = cumulative[starts] - scores.score[starts]
    draws = np.random.default_rng(seed).random(n)
    picked = np.searchsorted(cumulative, base + draws, side="right")
    # ⚠ Clamped into the run. A uniform arbitrarily close to 1 can exceed a run's final cumulative value
    # by one float rounding and land on the next run's first slot — which would be another fragment's
    # hypothesis, silently.
    return np.clip(picked, starts, ends - 1) - starts


def drain(
    payload: AccumulatorPayload,
    choices: np.ndarray,
    *,
    region_types: np.ndarray,
    sj,
) -> AccumulatorPayload:
    """⭐ **THE DRAIN.** Replay every held fragment with its chosen hypothesis; return the drained payload.

     ``choices[i]`` is a local hypothesis index for the ``i``-th record of the
    bank's canonical order — :func:`choose_hypotheses` produces exactly that.

    ⭐ **ONE TALLY PATH.** Each fragment re-enters ``Accumulator.deposit`` with its chosen hypothesis
    **alone**: a set of size one, so arbitration is degenerate and the ordinary rules decide whether it
    deposits or is rejected. There is no second deposit implementation and no duplicated crossing logic, so
    byte-identity with `tests/native/_accumulator_reference.py` holds for free rather than by argument.

    ⭐ **Pure: payload in, payload out.** The delta is accumulated into fresh per-reference accumulators
    and added to pass one's arrays, so both tallies exist afterwards and the drain's contribution to every
    channel is a subtraction. It also means the drain never sees a thread — the bank is already merged and
    canonically ordered — so reproducibility is structural rather than something a gate has to establish.
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
    chose_genomic = 0

    for i in range(n):
        h = int(deferred.hypothesis_offsets[i]) + int(choices[i])
        if not (int(deferred.hypothesis_offsets[i]) <= h < int(deferred.hypothesis_offsets[i + 1])):
            raise ValueError(f"choice {choices[i]} is outside record {i}'s hypothesis set")
        introns = [tuple(pair) for pair in deferred.hypothesis_introns_of(h).tolist()]
        chose_genomic += not introns
        outcome = accumulators[int(deferred.ref[i])].deposit(
            start=int(deferred.start[i]),
            end=int(deferred.end[i]),
            observed_introns=[tuple(p) for p in deferred.observed_introns_of(i).tolist()],
            align_strand=int(deferred.align_strand[i]),
            sj_strand=int(deferred.sj_strand[i]),
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
            **counts,
            chose_genomic=chose_genomic,
            chose_spliced=n - chose_genomic,
            census_before=payload.gap_resolution,
        ),
    )


def _gather_delta(
    payload: AccumulatorPayload, accumulators: _Accumulators
) -> dict[str, np.ndarray]:
    """The drained accumulators' channels, placed into globally-shaped arrays.

    ⚠ Only references the drain actually touched were built, so everything else stays zero — which is the
    correct delta and also the cheap path when a library's held fragments sit on a few contigs.
    """
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
                # ⚠ The pool histograms and the unconditional anchor are library-wide, so every
                # reference's contribution lands on the same rows.
                delta[name] += block.reshape(delta[name].shape)
            else:
                lo = int(axis_offsets[axis][ref])
                delta[name][lo : lo + block.shape[0]] += block
    return delta


@dataclass(frozen=True, slots=True)
class _Span:
    """The shape ``Accumulator.length_under`` reads a hypothesis in.

    ⚠ Duck-typed on purpose, matching the specification's ``GapHypothesis`` attribute names: the binding
    reads ``introns`` / ``sj_strand`` / ``supporting_t_inds`` off whatever it is handed, so the parity gate
    can pass the reference's own objects. A second tuple convention here would be a second representation
    to keep in step.
    """

    introns: list
    sj_strand: int

    @property
    def supporting_t_inds(self):  # the first pass never reads these, and neither does the scorer
        return ()
