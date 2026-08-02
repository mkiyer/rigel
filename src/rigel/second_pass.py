"""rigel.second_pass — scoring the side buffer's held fragments.

    Spec: ``docs/SPEC_SECOND_PASS.md`` (§3 the score, §4 the classes, P2 the phase this implements)

Pass 1 holds every fragment whose unsequenced gap has more than one surviving explanation — 2–3.5 % of a
library, and systematically the **long** ones, because a longer gap admits more hypotheses. This module
scores those explanations. It does **not** drain: the selection and the re-deposit are P3.

⭐ **THE PRIOR IS THE ACCUMULATOR ITSELF.** A hypothesis is a path, and a path is a set of accumulator
objects; pass 1 already counted how many molecules used each of them. So the score needs no transcript
abundance, no calibration output, and no second pass over the BAM::

    score(h)  =  rho(h)  x  f(L_h)  x  s(h)

⚠ **Every input comes from pass 1**, which is what lets this run BEFORE calibration and lets calibration
then run exactly once, on the complete tally. `SPEC_SECOND_PASS.md` §2.

⛔ **``L`` IS NOT COMPUTED HERE.** It comes from ``Accumulator.length_under`` — the same C++ that will
compute it again at drain time. `docs/FRAGMENT_LENGTH_AUDIT.md` C0/C2 left the tool with ONE definition of
fragment length, and a scorer with its own would be a second definition of exactly the quantity that audit
unified.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._bam_impl import Accumulator as NativeAccumulator
from .scan_payload import AccumulatorPayload
from .types import Strand


__all__ = ["HeldScores", "HypothesisTerms", "score_held_fragments"]


@dataclass(frozen=True, slots=True)
class HypothesisTerms:
    """The three factors, kept apart. ⭐ Diagnostics, not bookkeeping.

    ``SPEC_SECOND_PASS.md`` §4: the umbrella census partitions the held set by *which question is open*,
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


def _exact_cut(cuts: np.ndarray, lo: int, hi: int, position: int) -> int:
    """Flat cut index of ``position`` within ``cuts[lo:hi]``, or -1 if it is not a cut there."""
    k = lo + int(np.searchsorted(cuts[lo:hi], position))
    return k if k < hi and int(cuts[k]) == position else -1


def _junction_id(junctions, cuts, lo, hi, start, end, sj_strand) -> int:
    """The annotated junction slot for one intron, or -1.

    ⚠ Mirrors ``Accumulator::sj_edge_id`` — same CSR, same strand rule (the filter applies only when the
    motif strand is DEFINITE; a non-definite one matches on coordinates alone). ⭐ Not a duplicate of the
    length definition: this is a lookup into the index's own table, and the slot IS the payload's junction
    axis index.
    """
    donor = _exact_cut(cuts, lo, hi, start)
    if donor < 0:
        return -1
    acceptor = _exact_cut(cuts, lo, hi, end)
    if acceptor < 0:
        return -1
    definite = sj_strand in (int(Strand.POS), int(Strand.NEG))
    for k in range(int(junctions.offsets[donor]), int(junctions.offsets[donor + 1])):
        if int(junctions.acceptor_cut[k]) != acceptor:
            continue
        if definite and int(junctions.strand[k]) != sj_strand:
            continue
        return k
    return -1


def _lines_inside(cuts: np.ndarray, lo: int, hi: int, start: int, end: int) -> tuple[int, int]:
    """The half-open range of LOCAL line indices strictly inside ``[start, end)``.

    A line at cut ``c`` lies inside the interval when ``start < c < end``; the deposit's own crossing rule
    (`accumulator.cpp`) resolves it with the same two searches.
    """
    first = int(np.searchsorted(cuts[lo:hi], start, side="right"))
    last = int(np.searchsorted(cuts[lo:hi], end, side="left"))
    return first, last


def _aggregate(values: list[float], how: str) -> float:
    """Combine the densities of the several objects one path needs.

    ⛔ **OPEN DECISION D-1** (`SPEC_SECOND_PASS.md` §8). The molecule must be present at *every* object on
    its path, so the path's rate is not the sum. ``min`` is the bottleneck reading — the scarcest object
    bounds the path — and is the one with an argument behind it; ``geometric`` is less dominated by a
    single near-zero. ⚠ **Both are implemented so P2 can MEASURE which to keep.** Picked by taste this is
    a tunable; picked by measurement it is a derivation, and the measurement is the whole point of the
    parameter's existence. It is expected to be removed once the answer is in.
    """
    if not values:
        return 0.0
    if how == "min":
        return float(min(values))
    if how == "geometric":
        # ⚠ In log space: these are densities on the order of 1e-3 and a product over a long path
        # underflows float64 well before the path gets biologically implausible.
        positive = [v for v in values if v > 0.0]
        if len(positive) != len(values):
            return 0.0
        return float(np.exp(np.mean(np.log(positive))))
    raise ValueError(f"unknown aggregate {how!r}; expected 'min' or 'geometric'")


def score_held_fragments(
    payload: AccumulatorPayload,
    *,
    fl_models,
    rna_sense_frac: float,
    node_types: np.ndarray,
    junctions,
    aggregate: str = "min",
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
        ``P(align_strand == the transcript's strand | RNA)``. ⚠ On an R1-antisense (dUTP) library this is
        ≈ 0.01, so **disagreement is the likely case** — see `LEDGER.md` P0. It is used as the probability
        of agreement it is, never inverted.
    node_types, junctions
        From the index: ``build_node_partition_arrays`` and ``build_junction_edge_arrays``. ⚠ The junction
        CSR is genuinely not on the payload, so this dependency is real rather than an oversight.
    aggregate
        D-1, and temporary. See :func:`_aggregate`.
    """
    deferred = payload.deferred
    n_hyp = deferred.n_hypotheses
    length = np.zeros(n_hyp, np.int64)
    density = np.zeros(n_hyp, np.float64)
    length_likelihood = np.zeros(n_hyp, np.float64)
    strand = np.ones(n_hyp, np.float64)
    score = np.zeros(n_hyp, np.float64)

    cuts = payload.cut_positions
    scale = float(1 << 32)  # the fixed-point INV_LENGTH_SCALE; decoded once, here
    rna_pmf, global_pmf = fl_models.rna_pmf, fl_models.global_pmf
    max_size = int(fl_models.max_size)

    # One Accumulator per reference, for `length_under` alone. ⚠ Built from the INDEX's partition, which
    # is where node_types lives; the payload echoes the cut axis but not the typing.
    accumulators: dict[int, NativeAccumulator] = {}

    def accumulator_for(ref: int) -> NativeAccumulator:
        if ref not in accumulators:
            lo, hi = int(payload.ref_cut_offsets[ref]), int(payload.ref_cut_offsets[ref + 1])
            node_lo = int(payload.ref_node_offsets[ref])
            n_nodes = max(hi - lo - 1, 0)
            accumulators[ref] = NativeAccumulator(
                cuts=np.ascontiguousarray(cuts[lo:hi], dtype=np.int64),
                node_types=np.ascontiguousarray(
                    node_types[node_lo : node_lo + n_nodes], dtype=np.uint8
                ),
                max_length=payload.max_length,
                ref=ref,
            )
        return accumulators[ref]

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
        acc = accumulator_for(ref)
        cut_lo, cut_hi = int(payload.ref_cut_offsets[ref]), int(payload.ref_cut_offsets[ref + 1])
        edge_base = int(payload.ref_edge_offsets[ref])
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
                # A spliced path's evidence is the junctions it uses. ⚠ `sj_inv_length_sum` is deposited
                # with the SAME quantum as a contiguous edge, so the two are the same quantity on the
                # same scale — that is what makes this comparable to `∅`'s number at all.
                motif = observed_motif if observed_motif != int(Strand.NONE) else implied_strand
                observed_densities = []
                for a, b in introns:
                    jid = _junction_id(junctions, cuts, cut_lo, cut_hi, a, b, motif)
                    observed_densities.append(
                        0.0 if jid < 0 else float(payload.sj_inv_length_sum[jid].sum()) / scale
                    )
                density[slot] = _aggregate(observed_densities, aggregate)
            else:
                # The genomic path's evidence is the unspliced crossing density where the others jump.
                edge_densities = []
                for a, b in contested:
                    first, last = _lines_inside(cuts, cut_lo, cut_hi, a, b)
                    for line in range(first, last):
                        edge_densities.append(
                            float(payload.edge_unspliced_inv_length_sum[edge_base + line - 1].sum())
                            / scale
                        )
                density[slot] = _aggregate(edge_densities, aggregate)

            # -- f(L) -----------------------------------------------------------------------------
            pmf = rna_pmf if introns else global_pmf
            length_likelihood[slot] = float(pmf[L]) if 0 <= L <= max_size else 0.0

            # -- strand ---------------------------------------------------------------------------
            # ⭐ An OBSERVED motif pins the fragment's strand, and pass 1 has already constrained the
            # hypotheses to it (`LEDGER.md` P1), so the term is constant across the set and cancels.
            if observed_motif == int(Strand.NONE) and introns:
                if implied_strand in (int(Strand.POS), int(Strand.NEG)):
                    agrees = implied_strand == align
                    strand[slot] = rna_sense_frac if agrees else 1.0 - rna_sense_frac
                # ⚠ An AMBIGUOUS implied strand (D-5: one path, two strands) says nothing either way.
            # `∅` may be gDNA, which is strand-symmetric, so it keeps 1.0 and lets rho and f decide.

            score[slot] = density[slot] * length_likelihood[slot] * strand[slot]

        total = float(score[h0:h1].sum())
        if total > 0.0:
            score[h0:h1] /= total
        else:
            # ⛔ Not a drop and not a NaN. The evidence did not separate these hypotheses, so the honest
            # posterior is uniform and the draw is a coin toss — which is exactly what the multinomial in
            # P3 will do with it.
            score[h0:h1] = 1.0 / max(h1 - h0, 1)
            n_undecided += 1

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
