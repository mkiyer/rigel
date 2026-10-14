"""``Accumulator.offer`` — one fragment, several placements — in the SPECIFICATION, case by case.

A multimapping fragment is one molecule with several alignments, each on its own reference with its own
gap hypotheses. The rule (``_accumulator_reference.Accumulator.offer``): a placement excluded by its
strand or its clip excludes itself alone; identical placements are one; one placement remaining goes
through ``deposit`` bit for bit; several remaining are arbitrated as the UNION of their (placement,
hypothesis) pairs under the one length filter, one survivor deposits at its placement, more are held
WHOLE as ``deferred_multiple_placements``; the drain re-enters ``deposit`` at the chosen placement with
the chosen hypothesis alone. Every case here is a brute-force statement of that rule on a two-reference
partition, so the native ``AccumulatorSet::offer`` can be held to it array by array
(``test_accumulator_native_parity.py``).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.types import Strand

from ._accumulator_reference import (
    UNSPLICED_ONLY,
    Accumulator,
    DepositOutcome,
    GapHypothesis,
    Partition,
    Placement,
    Tally,
)

#: Two references with the same geometry, so a placement on reference 1 exercises every per-reference
#: offset and a constant stamp cannot pass.
BOUNDS = [0, 1000, 1600, 1700, 1900, 2000, 2200, 4000]
TYPES = [0, 2, 1, 1, 1, 2, 0]  # intergenic, exon, intron ×3, exon, intergenic
SJ = [
    (0, 1600, 2000, int(Strand.POS)),
    (0, 1700, 1900, int(Strand.NEG)),
    (1, 1600, 2000, int(Strand.POS)),
    (1, 1700, 1900, int(Strand.NEG)),
]
MAX_LENGTH = 1000
WIDE = GapHypothesis(introns=((1600, 2000),), sj_strand=int(Strand.POS), supporting_t_inds=(1,))
NARROW = GapHypothesis(introns=((1700, 1900),), sj_strand=int(Strand.NEG), supporting_t_inds=(2,))
GENOMIC = GapHypothesis()


def _fresh() -> Accumulator:
    partition = Partition.from_region_bounds([BOUNDS, BOUNDS], region_types=[TYPES, TYPES], sj=SJ)
    return Accumulator(partition, max_fragment_length=MAX_LENGTH)


def _pl(
    ref,
    start,
    end,
    align=int(Strand.POS),
    sj=int(Strand.NONE),
    observed=(),
    hypotheses=UNSPLICED_ONLY,
):
    return Placement(ref, start, end, align, sj, tuple(observed), tuple(hypotheses))


def _state(acc: Accumulator) -> dict:
    """Every field of the tally in a comparable form (the bank through its own canonical flattening)."""
    out = {}
    for f in dataclasses.fields(Tally):
        v = getattr(acc.tally, f.name)
        out[f.name] = acc.tally.deferred_arrays() if f.name == "deferred" else v
    return out


def _assert_same(a: Accumulator, b: Accumulator, label: str) -> None:
    sa, sb = _state(a), _state(b)
    for name in sa:
        x, y = sa[name], sb[name]
        if isinstance(x, dict):
            assert x.keys() == y.keys(), f"{label}: {name} keys"
            for k in x:
                if isinstance(x[k], np.ndarray):
                    assert x[k].dtype == y[k].dtype and np.array_equal(x[k], y[k]), (
                        f"{label}: {name}[{k}]"
                    )
                else:
                    assert x[k] == y[k], f"{label}: {name}[{k}] {x[k]} != {y[k]}"
        else:
            assert x.dtype == y.dtype and np.array_equal(x, y), f"{label}: {name}"


def _conserved(acc: Accumulator, offered: int) -> None:
    qc = acc.tally.qc
    total = (
        qc["deposited"]
        + qc["deferred_undetermined_gap"]
        + qc["deferred_multiple_placements"]
        + qc["dropped_too_long"]
        + qc["dropped_empty"]
        + qc["dropped_strand_undefined"]
    )
    assert total == offered, f"the identity does not close: {qc} against {offered} offered"


# ── ONE PLACEMENT IS `deposit`, BIT FOR BIT ────────────────────────────────────────────────────────

ONE_PLACEMENT_CASES = [
    ("contained in an exon", _pl(1, 1100, 1400)),
    ("crossing one boundary", _pl(1, 900, 1300)),
    (
        "observed spliced on the annotated sj",
        _pl(0, 1400, 2150, sj=int(Strand.POS), observed=[(1600, 2000)]),
    ),
    ("the three-hypothesis gap, held", _pl(1, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC))),
    ("strand undefined, dropped", _pl(0, 1100, 1400, align=int(Strand.NONE))),
    ("over the limit, dropped", _pl(0, 100, 1500)),
    ("off the reference, dropped empty", _pl(1, 5000, 5200)),
]


@pytest.mark.parametrize(
    "label,placement", ONE_PLACEMENT_CASES, ids=[c[0] for c in ONE_PLACEMENT_CASES]
)
def test_one_placement_is_deposit_bit_for_bit(label, placement):
    direct, offered = _fresh(), _fresh()
    want = direct.deposit(
        placement.ref,
        placement.start,
        placement.end,
        observed_introns=placement.observed_introns,
        align_strand=placement.align_strand,
        sj_strand=placement.sj_strand,
        hypotheses=placement.hypotheses,
    )
    got = offered.offer((placement,))
    assert got is want, f"{label}: {got} != {want}"
    _assert_same(direct, offered, label)
    _conserved(offered, 1)


def test_the_one_placement_battery_reaches_every_outcome_it_claims():
    outcomes = {_fresh().offer((pl,)) for _label, pl in ONE_PLACEMENT_CASES}
    assert outcomes == {
        DepositOutcome.DEPOSITED,
        DepositOutcome.DEFERRED,
        DepositOutcome.STRAND_UNDEFINED,
        DepositOutcome.TOO_LONG,
        DepositOutcome.EMPTY,
    }


# ── THE EXCLUSIONS ─────────────────────────────────────────────────────────────────────────────────


def test_identical_placements_are_one_placement():
    p = _pl(1, 1100, 1400)
    twice, once = _fresh(), _fresh()
    assert twice.offer((p, p)) is DepositOutcome.DEPOSITED
    once.deposit(1, 1100, 1400)
    _assert_same(twice, once, "a duplicated placement")
    assert twice.tally.qc["deferred_multiple_placements"] == 0


def test_a_strand_undefined_placement_excludes_only_itself():
    ok = _pl(1, 1100, 1400)
    bad = _pl(0, 1100, 1400, align=int(Strand.AMBIGUOUS))
    offered, direct = _fresh(), _fresh()
    assert offered.offer((bad, ok)) is DepositOutcome.DEPOSITED
    direct.deposit(1, 1100, 1400)
    _assert_same(offered, direct, "the defined placement decides alone")
    assert offered.tally.qc["dropped_strand_undefined"] == 0, (
        "an excluded placement is not a dropped fragment"
    )


def test_every_placement_excluded_rejects_the_fragment_once_and_names_the_strand_only_when_it_was_the_strand():
    both_strand = _fresh()
    assert (
        both_strand.offer(
            (_pl(0, 1100, 1400, align=int(Strand.NONE)), _pl(1, 1100, 1400, align=int(Strand.NONE)))
        )
        is DepositOutcome.STRAND_UNDEFINED
    )
    assert both_strand.tally.qc["dropped_strand_undefined"] == 1
    _conserved(both_strand, 1)
    mixed = _fresh()
    assert (
        mixed.offer((_pl(0, 1100, 1400, align=int(Strand.NONE)), _pl(1, 5000, 5100)))
        is DepositOutcome.EMPTY
    )
    assert mixed.tally.qc["dropped_empty"] == 1 and mixed.tally.qc["dropped_strand_undefined"] == 0
    _conserved(mixed, 1)


# ── THE UNION ARBITRATION ──────────────────────────────────────────────────────────────────────────


def test_every_pair_over_the_limit_means_all_stand_and_the_fragment_is_held():
    """Two contiguous placements of 1,500 bp, both over the 1,000 bp limit: the filter would empty the
    union, so every pair stands and the fragment is held with both — exactly the gap rule's clause."""
    acc = _fresh()
    assert acc.offer((_pl(0, 100, 1600), _pl(1, 100, 1600))) is DepositOutcome.DEFERRED_PLACEMENTS
    assert len(acc.tally.deferred) == 1 and acc.tally.deferred[0].n_placements == 2
    assert acc.tally.qc["deferred_multiple_placements"] == 1
    assert acc.tally.qc["dropped_too_long"] == 0, "nothing is rejected until the drain chooses"


def test_one_survivor_across_placements_deposits_at_its_placement_and_leaves_the_gap_census_alone():
    """Placement A is over the limit, placement B is a spliced path within it: B alone survives, the
    fragment deposits at B with that hypothesis, and the umbrella census is not touched — a direct
    deposit of B would have recorded ``gap_resolved_spliced``; a multi-placement fragment never does."""
    too_long = _pl(0, 100, 1600)
    spliced = _pl(1, 1400, 2150, hypotheses=(WIDE,))
    offered, direct = _fresh(), _fresh()
    assert offered.offer((too_long, spliced)) is DepositOutcome.DEPOSITED
    direct.deposit(1, 1400, 2150, hypotheses=(WIDE,))
    assert direct.tally.gap_resolution["gap_resolved_spliced"] == 1, (
        "the control: a direct deposit censuses"
    )
    assert all(v == 0 for v in offered.tally.gap_resolution.values())
    for name in (
        "region_contained_count",
        "boundary_unspliced_count",
        "boundary_spliced_count",
        "sj_count",
        "sj_inv_length_sum",
        "pool_lengths",
        "deposited_lengths",
    ):
        assert np.array_equal(getattr(offered.tally, name), getattr(direct.tally, name)), name
    assert offered.tally.qc == direct.tally.qc
    _conserved(offered, 1)


def test_two_contained_placements_are_held_whole_in_canonical_order_whatever_the_hit_order():
    """The paralog case: a fragment contained in an exon on each reference, nothing to tell them apart.
    Held whole, both placements, each with its unspliced hypothesis, the placements in their own order."""
    a, b = _pl(1, 1100, 1400), _pl(0, 1150, 1450)
    forward, backward = _fresh(), _fresh()
    assert forward.offer((a, b)) is DepositOutcome.DEFERRED_PLACEMENTS
    assert backward.offer((b, a)) is DepositOutcome.DEFERRED_PLACEMENTS
    held = forward.tally.deferred[0]
    assert [p.ref for p in held.placements] == [0, 1], (
        "reference 0's placement first, whatever the hit order"
    )
    assert held.placements[1] == a and held.placements[0] == b
    _assert_same(forward, backward, "hit order")
    assert forward.tally.qc["deferred_multiple_placements"] == 1
    assert forward.tally.qc["deferred_undetermined_gap"] == 0
    assert all(v == 0 for v in forward.tally.gap_resolution.values())
    for name in (
        "region_contained_count",
        "region_start_count",
        "deposited_lengths",
        "pool_lengths",
    ):
        assert int(getattr(forward.tally, name).sum()) == 0, (
            f"{name}: a held fragment locates nowhere"
        )
    _conserved(forward, 1)


def test_the_held_placement_carries_the_clipped_extent():
    acc = _fresh()
    assert acc.offer((_pl(0, -50, 400), _pl(1, 3800, 4300))) is DepositOutcome.DEFERRED_PLACEMENTS
    held = acc.tally.deferred[0]
    assert [(p.ref, p.start, p.end) for p in held.placements] == [(0, 0, 400), (1, 3800, 4000)]


# ── THE FLATTENING ─────────────────────────────────────────────────────────────────────────────────


def test_the_bank_flattens_with_a_placement_level_and_round_trips():
    """A one-placement gap record beside a two-placement record: three offset levels, every per-placement
    array holding exactly what the per-fragment array held before placements existed."""
    acc = _fresh()
    acc.deposit(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC))  # held on its gap, one placement
    acc.offer((_pl(1, 1100, 1400), _pl(0, 1150, 1450)))  # held on its placements
    arrays = acc.tally.deferred_arrays()
    assert arrays["placement_offsets"].tolist() == [0, 2, 3], (
        "records: (ref 0 @1150, ref 1 @1100) then (ref 0 @1500)"
    )
    assert arrays["ref"].tolist() == [0, 1, 0]
    assert arrays["start"].tolist() == [1150, 1100, 1500]
    assert arrays["hypothesis_offsets"].tolist() == [0, 1, 2, 5], (
        "one hypothesis per contained placement, three on the gap"
    )
    assert arrays["hypothesis_intron_offsets"].tolist() == [0, 0, 0, 1, 2, 2]
    assert arrays["hypothesis_t_offsets"].tolist() == [0, 0, 0, 1, 2, 2]
    assert all(v.dtype == np.int64 for v in arrays.values())
    assert set(arrays) == {
        "placement_offsets",
        "ref",
        "start",
        "end",
        "align_strand",
        "sj_strand",
        "observed_intron_offsets",
        "observed_introns",
        "hypothesis_offsets",
        "hypothesis_sj_strand",
        "hypothesis_intron_offsets",
        "hypothesis_introns",
        "hypothesis_t_offsets",
        "hypothesis_t",
    }


def test_the_canonical_order_does_not_depend_on_deposit_order():
    forward, backward = _fresh(), _fresh()
    forward.offer((_pl(1, 1100, 1400), _pl(0, 1150, 1450)))
    forward.deposit(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC))
    backward.deposit(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC))
    backward.offer((_pl(0, 1150, 1450), _pl(1, 1100, 1400)))
    assert [f.n_placements for f in forward.tally.deferred] != [
        f.n_placements for f in backward.tally.deferred
    ]
    _assert_same(forward, backward, "deposit order")


# ── THE DRAIN ──────────────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("choice,ref", [(0, 0), (1, 1)])
def test_draining_a_placement_choice_equals_depositing_that_placement_directly(choice, ref):
    a, b = _pl(1, 1100, 1400), _pl(0, 1150, 1450)
    held, direct = _fresh(), _fresh()
    held.offer((a, b))
    counters = held.drain([choice])
    chosen = {0: b, 1: a}[choice]
    assert chosen.ref == ref
    direct.deposit(chosen.ref, chosen.start, chosen.end)
    _assert_same(held, direct, f"choice {choice}")
    assert counters["offered"] == 1 and counters["offered_multimapper"] == 1
    assert counters["deposited"] == 1 and counters["chose_genomic"] == 1
    assert held.tally.deferred == [] and held.tally.qc["deferred_multiple_placements"] == 0
    _conserved(held, 1)


def test_a_record_local_choice_names_a_placement_and_a_hypothesis():
    """Placement A carries the three gap hypotheses and placement B its unspliced one: the record's run
    is A's three then B's one, so choice 3 is B and choice 1 is A's second."""
    acc = _fresh()
    acc.offer((_pl(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)), _pl(1, 1100, 1400)))
    held = acc.tally.deferred[0]
    assert [len(p.hypotheses) for p in held.placements] == [3, 1]
    placement, path = Accumulator._locate(held, 3)
    assert placement.ref == 1 and path is UNSPLICED_ONLY[0]
    placement, path = Accumulator._locate(held, 1)
    assert placement.ref == 0 and path is NARROW
    with pytest.raises(IndexError):
        Accumulator._locate(held, 4)
    drained, direct = _fresh(), _fresh()
    drained.offer((_pl(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)), _pl(1, 1100, 1400)))
    drained.drain([1])
    direct.deposit(0, 1500, 2100, hypotheses=(NARROW,))
    for name in (
        "region_contained_count",
        "boundary_unspliced_count",
        "sj_count",
        "sj_inv_length_sum",
        "deposited_lengths",
    ):
        assert np.array_equal(getattr(drained.tally, name), getattr(direct.tally, name)), name


def test_a_drained_choice_over_the_limit_is_rejected_as_too_long_and_still_conserves():
    acc = _fresh()
    acc.offer((_pl(0, 100, 1600), _pl(1, 100, 1600)))
    counters = acc.drain([0])
    assert counters["dropped_too_long"] == 1 and counters["deposited"] == 0
    assert counters["deposited"] + counters["dropped_too_long"] == counters["offered"] == 1
    _conserved(acc, 1)
