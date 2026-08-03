"""⭐ The accumulator arbitrates: one surviving hypothesis deposits, two or more are deferred.

       Owner ruling, 2026-08-01

A fragment's unsequenced mate gap may hold no intron, one, or several, and which it is **cannot be
observed** — the bases are not there. It is a likelihood question, and the likelihood needs a
fragment-length distribution that does not exist until the first pass is over. So the first pass does not
guess: it enumerates, and it either finds one answer or holds the fragment for the second pass.

⛔ **The decision moved.** `bam_scanner.cpp` used to compute a bool `path_ambiguous` and the accumulator
obeyed it, and `accumulator.h` said why: *"the accumulator cannot decide it — only the caller has the
candidate list."* It could not decide because the caller **collapsed the answer before handing it over**.
Given the hypothesis set, the decision belongs where the outcome was already reported.
"""

from __future__ import annotations

import numpy as np

from rigel.types import Strand

from ._accumulator_reference import (
    Accumulator,
    DepositOutcome,
    GapHypothesis,
    Partition,
)


#: One reference, cut at every coordinate the owner's §1 example names, so both paths land on real nodes.
CUTS = [0, 1000, 2000, 3000, 3050, 4000, 5000, 6000]
TYPES = [0, 2, 1, 2, 1, 2, 0]


def _acc(**kw):
    return Accumulator(Partition.from_cuts([CUTS], node_types=[TYPES]), **kw)


# ── the owner's §1 example, to the base pair ───────────────────────────────────────────────────────
#
#     TA  exons (1000,2000)  (3000,3050)  (4000,5000)   introns (2000,3000) (3050,4000)
#     TB  exons (1000,2000)               (4000,5000)   intron  (2000,4000)
#
#     fragment  blocks [1800,1950)                [4050,4200)      span 2400
#                      ==========|~~ unsequenced ~~|=========
#
#: ⭐ TA's path crosses an exon — (3000,3050) — that no read ever touched. A gap hypothesis is a PATH
#: through the annotation, not "an intron", and that is the whole reason the first pass cannot resolve it.
TA = GapHypothesis(((2000, 3000), (3050, 4000)), sj_strand=Strand.POS, supporting_t_inds=(11,))
TB = GapHypothesis(((2000, 4000),), sj_strand=Strand.POS, supporting_t_inds=(22,))
L_TA = 2400 - (1000 + 950)  # 450
L_TB = 2400 - 2000  # 400


def test_the_two_compatible_paths_have_the_lengths_the_owner_computed():
    """⚠ Pinned separately from the arbitration, because a deferred queue holding the right fragment with the
    wrong lengths would still deferred queue it — the second pass would then choose between two wrong answers."""
    assert (L_TA, L_TB) == (450, 400)
    acc = _acc()
    for path, expected in ((TA, L_TA), (TB, L_TB)):
        length, _introns, _absorbed = acc._hypothesis_length(1800, 4200, (), path)
        assert length == expected


def test_TWO_COMPATIBLE_PATHS_are_BUFFERED_and_neither_is_chosen():
    """⭐ The headline. 450 bp and 400 bp are both real molecules and nothing sequenced says which."""
    acc = _acc()
    outcome = acc.deposit(0, 1800, 4200, hypotheses=(TA, TB))
    assert outcome is DepositOutcome.DEFERRED
    t = acc.tally
    assert int(t.node_start_count.sum()) == 0, "an undetermined fragment locates nowhere"
    assert int(t.deposited_lengths.sum()) == 0, "and has no length to bin"
    assert int(t.pool_lengths.sum()) == 0
    assert len(t.deferred) == 1
    held = t.deferred[0]
    # ⭐ Stored WHOLE, with both paths and the transcripts supporting each: the second pass cannot
    # choose between answers it was not given, and it weights them by those transcripts' abundance.
    assert held.hypotheses == (TA, TB)
    assert [path.supporting_t_inds for path in held.hypotheses] == [(11,), (22,)]
    # ⛔ Both are spliced and neither is the unspliced hypothesis — the molecule is certified RNA either way and
    # the open question is purely WHICH STRUCTURE.
    assert t.gap_resolution["gap_deferred_which_introns"] == 1


def test_the_SAME_fragment_deposits_when_only_ONE_transcript_is_compatible():
    """The discriminating arm: the ambiguity is a property of the annotation, not of the fragment.

    Drop TB from the candidate set and nothing about the molecule changes — but there is now one answer,
    so it deposits at TA's 450 bp.
    """
    acc = _acc()
    assert acc.deposit(0, 1800, 4200, hypotheses=(TA,)) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(np.nonzero(t.deposited_lengths)[0][0]) == L_TA
    assert not t.deferred
    assert t.gap_resolution["gap_resolved_spliced"] == 1


def test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST():
    """⭐ The property the census's shape rests on, pinned as the REASON and not as a consequence.

    A spliced hypothesis cuts bases the genomic one keeps, so ``L_spliced <= L_genomic`` — always, for any
    extent, any observed introns and any implied path. The single arbitration filter is
    ``L <= max_fragment_length``, so it follows that **if the genomic path survives, every spliced path
    survives too**, and the survivor set can never be exactly ``{genomic}`` while a spliced path was
    offered.

    ⛔ That is why :class:`GapResolution` has no ``RESOLVED_UNSPLICED``: the class existed, was documented
    as "every spliced path was ruled out by length", and **no fragment could enter it**. Deleting the class
    and pinning the ordering here is what stops it coming back — a future filter that made the ordering
    false would fail this test rather than quietly need a class that is gone.

    ⚠ Randomised over the whole coordinate space rather than over hand-picked cases, because the claim is
    universal. Introns are drawn unsorted and may be reversed, zero-length, overlapping, nested, or
    entirely outside the fragment — every configuration ``_normalise_introns`` exists for.
    """
    import random

    acc = _acc()
    rng = random.Random(20260801)
    checked = 0
    for _ in range(20_000):
        start, end = sorted(rng.sample(range(0, 6001), 2))
        if end == start:
            continue
        observed = tuple(tuple(rng.sample(range(0, 6001), 2)) for _ in range(rng.randint(0, 2)))
        genomic, _introns, _absorbed = acc._hypothesis_length(start, end, observed, GapHypothesis())
        for _ in range(rng.randint(1, 3)):
            spliced_path = GapHypothesis(
                tuple(tuple(rng.sample(range(0, 6001), 2)) for _ in range(rng.randint(1, 3)))
            )
            spliced, _introns, _absorbed = acc._hypothesis_length(
                start, end, observed, spliced_path
            )
            checked += 1
            assert spliced <= genomic, (
                f"[{start},{end}) observed={observed} path={spliced_path.introns}: the spliced path "
                f"L={spliced} exceeds the genomic L={genomic}. Cutting bases cannot lengthen a molecule, "
                f"so this breaks the ordering GapResolution's shape depends on."
            )
    assert checked > 15_000, f"only {checked} pairs compared; the sweep is not doing its job"


def test_the_genomic_hypothesis_is_the_EMPTY_path_and_needs_no_flag():
    """⭐ Cutting nothing means the gap is real template — the molecule is gDNA, or nascent, which is the
    same unspliced span. So "could this be genomic?" needs no separate flag, and the nascent shadow
    transcript is not a candidate: it IS this hypothesis.

    ⚠ The limit is raised above the 2400 bp span on purpose. At the default 1000 the unspliced hypothesis is
    ruled out by length and the fragment deposits — which is the NEXT test, and keeping the two apart is
    what stops this one passing for the wrong reason.
    """
    acc = _acc(max_fragment_length=3000)
    outcome = acc.deposit(0, 1800, 4200, hypotheses=(TA, GapHypothesis()))
    assert outcome is DepositOutcome.DEFERRED
    # ⛔ The genomic path against exactly one spliced path: the open question is RNA or gDNA — one bit,
    # and it is the composition question calibration exists to answer.
    assert acc.tally.gap_resolution["gap_deferred_rna_or_gdna"] == 1
    assert acc.tally.gap_resolution["gap_deferred_which_introns"] == 0


def test_a_fragment_with_NO_gap_hypothesis_is_not_in_the_umbrella_at_all():
    """⚠ Non-vacuity for the census: a fragment that never had a question to answer must not be counted
    as having answered one, or the umbrella's denominator silently becomes the whole library."""
    acc = _acc()
    assert acc.deposit(0, 1800, 2400) is DepositOutcome.DEPOSITED
    assert sum(acc.tally.gap_resolution.values()) == 0


# ── the span rule is the max_fragment_length filter applied to the unspliced hypothesis ─────────────────────


def test_a_span_over_the_limit_RULES_OUT_the_genomic_hypothesis():
    """⭐ The owner's rule "if the genomic span exceeds max_fragment_length, assume it is RNA" is not a
    separate rule: the unspliced hypothesis's ``L`` **is** that span, so the ordinary hypothesis filter deletes
    it. Span 2400 against a limit of 1000 leaves TA alone, and it deposits at 450.
    """
    acc = _acc(max_fragment_length=1000)
    assert acc.deposit(0, 1800, 4200, hypotheses=(TA, GapHypothesis())) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(np.nonzero(t.deposited_lengths)[0][0]) == L_TA
    assert t.gap_resolution["gap_resolved_spliced"] == 1
    assert not t.deferred


def test_under_the_limit_the_unspliced_hypothesis_SURVIVES_and_the_fragment_is_DEFERRED():
    """⛔ The other side of the same rule, and the population the owner says is always deferred: a short
    unspliced fragment with an intron in its gap is compatible with DNA **and** with RNA."""
    acc = _acc(max_fragment_length=1000)
    # span 900, and a 400 bp intron inside the gap -> unspliced L = 900, spliced L = 500. Both possible.
    outcome = acc.deposit(
        0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())
    )
    assert outcome is DepositOutcome.DEFERRED
    assert acc.tally.gap_resolution["gap_deferred_rna_or_gdna"] == 1


# ── conservation — G1's gate, and it is the only thing G1 may be judged by ──────────────────────────


def test_EVERY_OFFERED_FRAGMENT_IS_ACCOUNTED_FOR_EXACTLY_ONCE():
    """⭐ **THE G1 GATE.** ``deposited + deferred + dropped_* == offered``, exactly.

    ⛔ G1 must NOT be judged by a calibration A/B. Between G1 and G2 the tally is deliberately thinner —
    the ambiguous mass is retained in the deferred queue but not yet deposited — so accuracy would read worse for
    a reason that is the design working. Conservation is what says nothing was lost.

    ⚠ Every branch below is exercised on purpose. A conservation identity over a population that only
    ever takes one route is satisfied by any bookkeeping at all.
    """
    acc = _acc(max_fragment_length=1000)
    offered = 0

    def offer(*args, **kw):
        nonlocal offered
        offered += 1
        return acc.deposit(*args, **kw)

    assert offer(0, 1800, 2400) is DepositOutcome.DEPOSITED
    assert offer(0, 1800, 4200, hypotheses=(TA,)) is DepositOutcome.DEPOSITED
    assert offer(0, 1800, 4200, hypotheses=(TA, TB)) is DepositOutcome.DEFERRED
    assert offer(0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())) is (
        DepositOutcome.DEFERRED
    )
    assert offer(0, 1800, 4200) is DepositOutcome.TOO_LONG
    assert offer(0, 1800, 2400, align_strand=Strand.AMBIGUOUS) is DepositOutcome.STRAND_UNDEFINED
    assert offer(0, 9000, 9500) is DepositOutcome.EMPTY

    qc = acc.tally.qc
    assert sum(qc[outcome.value] for outcome in DepositOutcome) == offered
    # ⭐ And the deferred queue is not merely counted — it holds the fragments the counter claims.
    assert len(acc.tally.deferred) == qc["deferred_undetermined_gap"] == 2
    assert int(acc.tally.node_start_count.sum()) == qc["deposited"] == 2


def test_the_gap_resolution_SUBCLASSES_CLOSE_against_the_umbrella_and_the_deferred_queue():
    """⭐ Owner ruling: one umbrella, carefully partitioned, so all of it is counted.

    Two identities, and they are different statements. The first says the subclasses partition the
    umbrella; the second says the three DEFERRED_* subclasses are exactly what went into the deferred queue.
    """
    acc = _acc(max_fragment_length=1000)
    acc.deposit(0, 1800, 4200, hypotheses=(TA,))  # determined, spliced
    acc.deposit(0, 1800, 4200, hypotheses=(TA, TB))  # deferred, path only
    acc.deposit(
        0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())
    )  # deferred, component
    acc.deposit(
        0,
        1800,
        2700,
        hypotheses=(
            GapHypothesis(((2100, 2500),)),
            GapHypothesis(((2200, 2600),)),
            GapHypothesis(),
        ),
    )
    acc.deposit(0, 1800, 2400)  # no hypothesis — not in the umbrella

    census = acc.tally.gap_resolution
    assert census == {
        "gap_resolved_spliced": 1,
        "gap_deferred_rna_or_gdna": 1,
        "gap_deferred_which_introns": 1,
        "gap_deferred_both": 1,
    }
    deferred = sum(v for k, v in census.items() if k.startswith("gap_deferred"))
    assert deferred == acc.tally.qc["deferred_undetermined_gap"] == len(acc.tally.deferred)


# ── the flattened deferred queue is what the payload carries, and it is specified in the reference ──────────


#: The two deferred fragments the flattening tests use. ⚠ ``SHORT`` sorts BEFORE ``LONG`` — same start,
#: smaller end — so the expected arrays below are in that order and not in deposit order.
LONG = (0, 1800, 4200, (TA, TB))
SHORT = (0, 1800, 2700, (GapHypothesis(((2100, 2500),)), GapHypothesis()))


def test_the_deferred_queue_FLATTENS_to_a_CSR_that_round_trips():
    """⚠ The reference stores records for readability; the payload carries arrays. Both must be ONE
    representation, so the flattening is specified in the reference and gated here rather than being
    argued equal across two languages."""
    acc = _acc(max_fragment_length=1000)
    for ref, start, end, hypotheses in (LONG, SHORT):
        acc.deposit(ref, start, end, hypotheses=hypotheses)

    arrays = acc.tally.deferred_arrays()
    assert arrays["start"].tolist() == [1800, 1800]
    assert arrays["end"].tolist() == [2700, 4200]
    # Offsets are cumulative and start at 0, so an empty deferred queue is [0] and never [].
    assert arrays["hypothesis_offsets"].tolist() == [0, 2, 4]
    # SHORT's spliced path has one intron and its genomic path none; then TA has two and TB one.
    assert arrays["hypothesis_intron_offsets"].tolist() == [0, 1, 1, 3, 4]
    assert arrays["hypothesis_introns"].tolist() == [
        2100, 2500,               # SHORT, spliced
                                  # SHORT, genomic — no introns, which is what makes it genomic
        2000, 3000, 3050, 4000,   # TA
        2000, 4000,               # TB
    ]  # fmt: skip
    assert arrays["hypothesis_t"].tolist() == [11, 22]
    assert arrays["hypothesis_t_offsets"].tolist() == [0, 0, 0, 1, 2]
    assert Accumulator(Partition.from_cuts([CUTS], node_types=[TYPES])).tally.deferred_arrays()[
        "hypothesis_offsets"
    ].tolist() == [0]


def test_the_FLATTENED_queue_DOES_NOT_DEPEND_ON_DEPOSIT_ORDER():
    """⭐ What the sort is FOR, and the reason it is in the reference rather than in the exporter.

    Every other bank is a sum of integers and integer addition is associative, so a per-worker merge is
    exact whatever order the chunks arrived in. ⛔ The deferred queue is a LIST — concatenating per-worker
    deferred queues would give a different byte sequence at 1, 2, 4 and 8 workers with identical contents, and
    `test_accumulator_worker_determinism.py` would fail on a difference that means nothing. Sorting on
    the record's own content is the canonical form.
    """
    forward, backward = _acc(max_fragment_length=1000), _acc(max_fragment_length=1000)
    for ref, start, end, hypotheses in (LONG, SHORT):
        forward.deposit(ref, start, end, hypotheses=hypotheses)
    for ref, start, end, hypotheses in (SHORT, LONG):
        backward.deposit(ref, start, end, hypotheses=hypotheses)

    assert [f.end for f in forward.tally.deferred] != [f.end for f in backward.tally.deferred], (
        "the two accumulators must actually differ in deposit order, or this passes vacuously"
    )
    a, b = forward.tally.deferred_arrays(), backward.tally.deferred_arrays()
    assert a.keys() == b.keys()
    for name in a:
        assert np.array_equal(a[name], b[name]), f"{name} depends on the order fragments arrived in"
