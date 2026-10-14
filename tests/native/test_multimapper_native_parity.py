"""``AccumulatorSet::offer`` against the specification, case by case and on a random battery.

The native set is built from the specification's own two-reference ``Partition`` (the flat partition the
scanner installs), both sides are offered the SAME ``Placement`` objects, and after every fragment every
reference's banks, the counters, the umbrella census and the ONE canonical bank (every reference's
gap-held records and the set's multi-placement records in one order) must be byte-identical. The
one-placement path is the parity battery's ``deposit`` reached through ``offer``; the multi-placement
cases are ``test_multimapper_offer.py``'s, stated once there and held here against the C++.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel.native import AccumulatorSet as NativeSet
from rigel.types import Strand

from ._accumulator_reference import (
    UNSPLICED_ONLY,
    Accumulator,
    GapHypothesis,
    Partition,
    Tally,
)
from .test_multimapper_offer import (
    BOUNDS,
    GENOMIC,
    MAX_LENGTH,
    NARROW,
    ONE_PLACEMENT_CASES,
    SJ,
    TYPES,
    WIDE,
    _pl,
)

#: The Tally's library-wide arrays: summed over the native set's references for the comparison.
_LIBRARY_WIDE = {"pool_lengths", "deposited_lengths"}


def _pair():
    partition = Partition.from_region_bounds([BOUNDS, BOUNDS], region_types=[TYPES, TYPES], sj=SJ)
    reference = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    native = NativeSet(
        region_bounds=np.ascontiguousarray(partition.region_bounds, dtype=np.int64),
        ref_region_bound_offsets=np.ascontiguousarray(
            partition.ref_region_bound_offsets, dtype=np.int64
        ),
        region_types=np.ascontiguousarray(partition.region_types, dtype=np.uint8),
        max_length=MAX_LENGTH,
    )
    native.set_sj(
        np.ascontiguousarray(partition.sj_offsets, dtype=np.int64),
        np.ascontiguousarray(partition.sj_boundary_right, dtype=np.int64),
        np.ascontiguousarray(partition.sj_strand, dtype=np.int8),
        np.ascontiguousarray(partition.ref_region_bound_offsets, dtype=np.int64),
    )
    return reference, native


def _slices(partition: Partition, f: int) -> dict[str, slice]:
    c0, c1 = (
        int(partition.ref_region_bound_offsets[f]),
        int(partition.ref_region_bound_offsets[f + 1]),
    )
    return {
        "region": slice(
            int(partition.ref_region_offsets[f]), int(partition.ref_region_offsets[f + 1])
        ),
        "boundary": slice(
            int(partition.ref_boundary_offsets[f]), int(partition.ref_boundary_offsets[f + 1])
        ),
        "sj": slice(int(partition.sj_offsets[c0]), int(partition.sj_offsets[c1])),
    }


def _assert_parity(reference: Accumulator, native: NativeSet, label: str) -> None:
    tally, partition = reference.tally, reference.partition
    for field in dataclasses.fields(Tally):
        want = getattr(tally, field.name)
        if field.name == "deferred":
            got = dict(native.deferred)
            want = tally.deferred_arrays()
            assert got.keys() == want.keys(), (
                f"{label}: deferred keys {sorted(got)} != {sorted(want)}"
            )
            for k in want:
                assert got[k].dtype == want[k].dtype, f"{label}: deferred[{k}] dtype"
                assert np.array_equal(got[k], want[k]), (
                    f"{label}: deferred[{k}] {got[k].tolist()} != {want[k].tolist()}"
                )
            continue
        if field.name in ("qc", "gap_resolution"):
            got = dict(getattr(native, field.name))
            assert got == dict(want), f"{label}: {field.name} {got} != {dict(want)}"
            continue
        if field.name in _LIBRARY_WIDE:
            got = sum(
                np.asarray(getattr(native.at(f), field.name), dtype=np.int64)
                for f in range(native.n_refs)
            )
            assert np.array_equal(got, np.asarray(want, dtype=np.int64)), f"{label}: {field.name}"
            continue
        axis = (
            "region"
            if field.name.startswith("region_")
            else "boundary"
            if field.name.startswith("boundary_")
            else "sj"
        )
        for f in range(native.n_refs):
            got = np.asarray(getattr(native.at(f), field.name))
            expected = np.asarray(want)[_slices(partition, f)[axis]]
            assert got.dtype == expected.dtype, f"{label}: {field.name} dtype on ref {f}"
            assert np.array_equal(got, expected), (
                f"{label}: {field.name} on ref {f}: {got.tolist()} != {expected.tolist()}"
            )


def _offer_both(reference, native, label, placements) -> None:
    want = reference.offer(placements)
    got = native.offer(placements)
    assert got == want.value, f"{label}: outcome {got!r} != {want.value!r}"
    _assert_parity(reference, native, label)


MULTI_CASES = [
    ("identical placements", (_pl(1, 1100, 1400), _pl(1, 1100, 1400))),
    (
        "a strand-undefined placement beside a defined one",
        (_pl(0, 1100, 1400, align=int(Strand.AMBIGUOUS)), _pl(1, 1100, 1400)),
    ),
    (
        "every placement strand-undefined",
        (_pl(0, 1100, 1400, align=int(Strand.NONE)), _pl(1, 1100, 1400, align=int(Strand.NONE))),
    ),
    (
        "an empty clip beside a strand-undefined one",
        (_pl(0, 1100, 1400, align=int(Strand.NONE)), _pl(1, 5000, 5100)),
    ),
    ("every pair over the limit", (_pl(0, 100, 1600), _pl(1, 100, 1600))),
    (
        "one survivor across placements, a spliced path",
        (_pl(0, 100, 1600), _pl(1, 1400, 2150, hypotheses=(WIDE,))),
    ),
    ("the paralog pair", (_pl(1, 1100, 1400), _pl(0, 1150, 1450))),
    ("clipped extents", (_pl(0, -50, 400), _pl(1, 3800, 4300))),
    (
        "a gap-ambiguous placement beside a contained one",
        (_pl(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)), _pl(1, 1100, 1400)),
    ),
    (
        "three placements, two identical",
        (_pl(1, 1100, 1400), _pl(0, 2050, 2150), _pl(1, 1100, 1400)),
    ),
    (
        "observed spliced placements on both references",
        (
            _pl(0, 1400, 2150, sj=int(Strand.POS), observed=[(1600, 2000)]),
            _pl(1, 1400, 2150, sj=int(Strand.POS), observed=[(1600, 2000)]),
        ),
    ),
]


def test_one_placement_through_offer_is_byte_identical():
    reference, native = _pair()
    _assert_parity(reference, native, "empty")
    for label, placement in ONE_PLACEMENT_CASES:
        _offer_both(reference, native, label, (placement,))


def test_every_multi_placement_case_is_byte_identical():
    reference, native = _pair()
    _offer_both(
        reference,
        native,
        "a gap-held unique mapper first",
        (_pl(1, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)),),
    )
    for label, placements in MULTI_CASES:
        _offer_both(reference, native, label, placements)
    assert reference.tally.qc["deferred_multiple_placements"] > 0
    assert reference.tally.qc["deferred_undetermined_gap"] > 0, (
        "the gap-held path must be in the mix"
    )
    bank = dict(native.deferred)
    assert (np.diff(bank["placement_offsets"]) > 1).any(), (
        "the union bank must hold a multi-placement record"
    )


def test_the_union_bank_interleaves_gap_held_and_multimapping_records_in_one_order():
    """Gap-held records on both references and multi-placement records between them: the native union
    (every reference's bank and the set's own, sorted once) is the specification's one list."""
    reference, native = _pair()
    _offer_both(
        reference, native, "gap on ref 1", (_pl(1, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)),)
    )
    _offer_both(reference, native, "paralog pair", (_pl(1, 1100, 1400), _pl(0, 1150, 1450)))
    _offer_both(
        reference, native, "gap on ref 0", (_pl(0, 1500, 2100, hypotheses=(WIDE, NARROW, GENOMIC)),)
    )
    _offer_both(
        reference, native, "pair starting on ref 0 later", (_pl(0, 1550, 1580), _pl(1, 2050, 2150))
    )
    bank = dict(native.deferred)
    assert bank["placement_offsets"].tolist() == [0, 2, 3, 5, 6]
    assert bank["ref"].tolist() == [0, 1, 0, 0, 1, 1]


def test_a_random_battery_is_byte_identical():
    """Two thousand random fragments with one to three placements on two references, random strands,
    observed introns and hypothesis sets; parity asserted every fifty fragments and at the end."""
    rng = np.random.default_rng(2026_10_11)
    reference, native = _pair()
    introns = [(1600, 2000), (1700, 1900), (1650, 1950)]
    hyps = (
        WIDE,
        NARROW,
        GENOMIC,
        GapHypothesis(introns=((1650, 1950),), sj_strand=int(Strand.POS), supporting_t_inds=(5,)),
    )
    for k in range(2000):
        n = int(rng.integers(1, 4))
        placements = []
        for _ in range(n):
            ref = int(rng.integers(0, 2))
            start = int(rng.integers(-100, 4000))
            end = start + int(rng.integers(1, 1800))
            align = int(
                rng.choice(
                    [
                        int(Strand.POS),
                        int(Strand.NEG),
                        int(Strand.POS),
                        int(Strand.NEG),
                        int(Strand.NONE),
                    ]
                )
            )
            observed = []
            sj = int(Strand.NONE)
            if rng.uniform() < 0.25:
                observed = [introns[int(rng.integers(0, 3))]]
                sj = int(rng.choice([int(Strand.POS), int(Strand.NEG), int(Strand.AMBIGUOUS)]))
            size = int(rng.integers(0, 4))
            chosen = (
                tuple(hyps[i] for i in sorted(rng.choice(4, size=size, replace=False).tolist()))
                or UNSPLICED_ONLY
            )
            if rng.uniform() < 0.3 and size:
                chosen = tuple(chosen) + (GENOMIC,)
            placements.append(
                _pl(ref, start, end, align=align, sj=sj, observed=observed, hypotheses=chosen)
            )
        want = reference.offer(tuple(placements))
        got = native.offer(tuple(placements))
        assert got == want.value, f"fragment {k}: {got!r} != {want.value!r} for {placements}"
        if k % 50 == 49:
            _assert_parity(reference, native, f"after fragment {k}")
    _assert_parity(reference, native, "the end")
    qc = reference.tally.qc
    assert qc["deferred_multiple_placements"] > 50 and qc["deferred_undetermined_gap"] > 50, qc
    assert (
        qc["deposited"] > 200
        and qc["dropped_too_long"] > 50
        and qc["dropped_strand_undefined"] > 10
    ), qc
