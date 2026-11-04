"""``lift_choices``: carrying a second-pass hypothesis choice from the whole library onto a subset.

Splitting a BAM by origin and re-scanning reconstructs the pass-one payload exactly, because pass one deposits each
fragment independently; that makes an origin-split oracle valid. The second pass's multinomial is scored against the
whole payload's densities, so draining each partition on its own would give ``sum(partitions) != whole``. Instead the
choices are scored and drawn once on the whole library and replayed inside whichever partition holds each fragment.
These tests pin what that rests on: a partition that is the whole gets its own choices back exactly; a proper subset
gets the choices its records drew on the whole, checked per record because an aggregate cannot see two records swapped;
a complete split consumes every choice exactly once. A record the whole does not hold and a wrong-length choice array
both raise, and the ambiguity count is pinned, since a duplicate key split across origins is the one case the lift
cannot resolve and the count bounds the error.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.scan_payload import DeferredFragments
from _oracle import lift_choices


def _bank(records, n_hyp_each=2):
    """A ``DeferredFragments`` holding ``records`` = list of (ref, start, end, astrand, sjstrand, nm, clip), built in
    the bank's own canonical order (sorted on the record key), the order ``lift_choices`` reads. Every record gets
    ``n_hyp_each`` hypotheses, which makes a local hypothesis index meaningful and transferable."""
    recs = sorted(records)
    n = len(recs)
    cols = list(zip(*recs)) if n else [()] * 7
    return DeferredFragments(
        placement_offsets=np.arange(0, n + 1, dtype=np.int64),  # one placement per record
        ref=np.asarray(cols[0], np.int64),
        start=np.asarray(cols[1], np.int64),
        end=np.asarray(cols[2], np.int64),
        align_strand=np.asarray(cols[3], np.int64),
        sj_strand=np.asarray(cols[4], np.int64),
        observed_intron_offsets=np.zeros(n + 1, np.int64),
        observed_introns=np.zeros(0, np.int64),
        hypothesis_offsets=np.arange(0, n_hyp_each * n + 1, n_hyp_each, dtype=np.int64),
        hypothesis_sj_strand=np.zeros(n_hyp_each * n, np.int64),
        hypothesis_intron_offsets=np.zeros(n_hyp_each * n + 1, np.int64),
        hypothesis_introns=np.zeros(0, np.int64),
        nm=np.asarray(cols[5], np.int64),
        clip=np.asarray(cols[6], np.int64),
    )


class _P:
    """The one attribute ``lift_choices`` reads off a payload."""

    def __init__(self, bank):
        self.deferred = bank


#: Distinct records in canonical order. Near-duplicate pairs differ only in `align_strand`, only in `sj_strand`, only in
#: `nm` and only in `clip`; if every record differed in (ref, start, end), dropping a trailing key field would change
#: nothing.
_RECS = [
    (0, 100, 300, 0, 0, 0, 0),
    (0, 400, 700, 0, 1, 0, 0),
    (0, 900, 1200, 1, 0, 0, 0),
    (1, 50, 260, 0, 0, 0, 0),
    (1, 50, 260, 1, 0, 0, 0),  # <- differs from the previous only in align_strand
    (1, 800, 1000, 0, 0, 0, 0),
    (1, 800, 1000, 0, 1, 0, 0),  # <- differs from the previous only in sj_strand
    (2, 10, 210, 0, 0, 0, 0),
    (2, 10, 210, 0, 0, 1, 0),  # <- differs from the previous only in nm
    (2, 10, 210, 0, 0, 1, 3),  # <- differs from the previous only in clip
]
_CHOICES = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 0], np.int64)


def test_a_partition_that_is_the_whole_gets_its_own_choices_back():
    """The degenerate case, and it must be exact: if the identity fails here it fails everywhere."""
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [whole], _CHOICES)
    assert ambiguous == 0
    assert np.array_equal(lifted[0], _CHOICES)


def test_a_proper_subset_gets_exactly_the_choices_its_records_drew():
    """Checked per record against the whole's own choice for the same key. An aggregate check (the same multiset of
    choices) would pass with two records swapped, which would misattribute one origin's mass to another."""
    keep = [_RECS[0], _RECS[2], _RECS[4], _RECS[6]]
    whole = _P(_bank(_RECS))
    part = _P(_bank(keep))
    lifted, ambiguous = lift_choices(whole, [part], _CHOICES)
    assert ambiguous == 0
    want = {r: int(_CHOICES[i]) for i, r in enumerate(sorted(_RECS))}
    for j, r in enumerate(sorted(keep)):
        assert int(lifted[0][j]) == want[r], (j, r, lifted[0][j], want[r])


def test_a_complete_split_consumes_every_choice_exactly_once():
    """Two disjoint partitions covering the whole must between them use each whole-library choice once, which gives
    ``Sum(partitions) == whole`` after the drain."""
    a = [_RECS[0], _RECS[3], _RECS[5], _RECS[8]]
    b = [_RECS[1], _RECS[2], _RECS[4], _RECS[6], _RECS[7], _RECS[9]]
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [_P(_bank(a)), _P(_bank(b))], _CHOICES)
    assert ambiguous == 0
    got = [int(x) for arr in lifted for x in arr]
    assert sorted(got) == sorted(int(x) for x in _CHOICES)
    assert len(got) == len(_RECS)


def test_a_record_the_whole_does_not_hold_raises():
    """A partition scanned from a different BAM, or a drifted payload pair, would otherwise produce a plausible tally,
    so it must be an error."""
    whole = _P(_bank(_RECS[:3]))
    part = _P(_bank([_RECS[5]]))
    with pytest.raises(ValueError, match="not in the whole library's held bank"):
        lift_choices(whole, [part], _CHOICES[:3])


def test_a_wrong_length_choice_array_raises():
    """One choice per held record is the contract; a mismatch means the choices came from a different
    scan and every index below would be off by an unknown amount."""
    whole = _P(_bank(_RECS))
    with pytest.raises(ValueError, match="One choice per held record"):
        lift_choices(whole, [whole], _CHOICES[:2])


def test_duplicate_keys_split_across_partitions_are_counted():
    """Two records with the same key are identical (``DeferredFragments``' guarantee), so no partition can know which of
    them it holds: the deposits are interchangeable, the origin attribution is not. The count is returned, and it bounds
    the truth error exactly. It is per record of the partition, so a partition holding one of a duplicate pair reports
    1; rounding up is the safe direction for a bound.
    """
    dup = (0, 100, 300, 0, 0, 0, 0)
    recs = [dup, dup, (0, 500, 700, 0, 0, 0, 0)]
    whole = _P(_bank(recs))
    ch = np.array([0, 1, 1], np.int64)
    # each partition takes one of the identical pair; neither can know which
    lifted, ambiguous = lift_choices(whole, [_P(_bank([dup])), _P(_bank([dup]))], ch)
    assert ambiguous == 2, ambiguous
    # The queue is consumed across partitions, which is the identity itself: the duplicate group's two choices are {0,
    # 1}, and per-call state would let both partitions take run[0] while one entry went unused.
    assert sorted([int(lifted[0][0]), int(lifted[1][0])]) == [0, 1]
    # and a partition of distinct records reports none
    l3, a3 = lift_choices(whole, [_P(_bank([recs[2]]))], ch)
    assert a3 == 0
    assert int(l3[0][0]) == 1


def test_the_key_uses_every_field_of_the_record_identity():
    """``_RECS`` contains pairs differing only in ``align_strand``, only in ``sj_strand``, only in ``nm`` and only in
    ``clip``, so a lift ignoring a trailing field would collapse each pair into one queue and hand the same choice to
    both, a mismatch against the whole's own per-record choice.
    """
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [whole], _CHOICES)
    assert ambiguous == 0, (
        "the fixture must have no duplicate FULL keys, or this gate proves nothing"
    )
    assert np.array_equal(lifted[0], _CHOICES)
    # the two near-duplicate pairs must have drawn different choices, or the test is vacuous
    order = sorted(_RECS)
    for a, b in (
        ((1, 50, 260, 0, 0, 0, 0), (1, 50, 260, 1, 0, 0, 0)),
        ((1, 800, 1000, 0, 0, 0, 0), (1, 800, 1000, 0, 1, 0, 0)),
        ((2, 10, 210, 0, 0, 0, 0), (2, 10, 210, 0, 0, 1, 0)),
        ((2, 10, 210, 0, 0, 1, 0), (2, 10, 210, 0, 0, 1, 3)),
    ):
        ca, cb = _CHOICES[order.index(a)], _CHOICES[order.index(b)]
        assert ca != cb, (a, b, ca, cb)


def test_partitions_holding_more_copies_of_a_key_than_the_whole_raises():
    """Two partitions each holding a record whose key is unique in the whole are not a partition of one scan (two scans,
    or a double-counted fragment). Without the guard the second partition would re-use the first's choice, giving
    ``sum(partitions) > whole``.
    """
    only = (0, 100, 300, 0, 0, 0, 0)
    whole = _P(_bank([only, (0, 500, 700, 0, 0, 0, 0)]))
    ch = np.array([0, 1], np.int64)
    with pytest.raises(ValueError, match="more records with key"):
        lift_choices(whole, [_P(_bank([only])), _P(_bank([only]))], ch)
