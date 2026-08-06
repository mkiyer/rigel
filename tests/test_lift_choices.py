"""``lift_choices`` — carrying a second-pass hypothesis choice from the WHOLE library onto a SUBSET.

⛔ **THE IDENTITY THIS PROTECTS.** Splitting a BAM by fragment origin and re-scanning reconstructs the
pass-one payload exactly, because pass one deposits each fragment independently — which is what makes an
origin-split oracle a valid truth source. The second pass breaks it: its multinomial is scored against the
*whole* payload's densities, so ``Sum(partitions) != whole`` if each partition is drained on its own.

⭐ The repair is to score and draw ONCE on the whole library and replay each fragment's chosen hypothesis
inside whichever partition holds it. These gates pin the three properties that makes rest on:

===  ===============================================================================================
L1   a partition that IS the whole gets its own choices back, exactly
L2   a proper subset gets exactly the choices its records drew on the whole — checked per record,
     not in aggregate, because an aggregate cannot see two records swapping
L3   the partitions PARTITION: every whole-library choice is consumed exactly once across a
     complete split.  ⭐ This is the identity itself
===  ===============================================================================================

Plus the two failure modes, both of which must be loud rather than silent: a record the whole does not
hold (L4), and a wrong-length choice array (L5). And L6 pins the ambiguity COUNT, because a duplicate key
split across origins is the one case the lift cannot resolve and the count is what bounds it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.scan_payload import DeferredFragments
from rigel.second_pass import lift_choices


def _bank(records, n_hyp_each=2):
    """A ``DeferredFragments`` holding ``records`` = list of (ref, start, end, astrand, sjstrand).

    ⚠ Built in the bank's own canonical order — sorted on the record key — because that order is the
    contract ``lift_choices`` reads. Every record gets ``n_hyp_each`` hypotheses, which is what makes a
    LOCAL hypothesis index meaningful and transferable."""
    recs = sorted(records)
    n = len(recs)
    cols = list(zip(*recs)) if n else [()] * 5
    return DeferredFragments(
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
        hypothesis_t_offsets=np.zeros(n_hyp_each * n + 1, np.int64),
        hypothesis_t=np.zeros(0, np.int64),
    )


class _P:
    """The one attribute ``lift_choices`` reads off a payload."""

    def __init__(self, bank):
        self.deferred = bank


#: Distinct records in canonical order. ⭐⭐ **THE LAST TWO PAIRS DIFFER ONLY IN `align_strand` AND ONLY
#: IN `sj_strand`.** That is deliberate and it is what gives the key-completeness gate teeth: with every
#: record differing in (ref, start, end) already, dropping a trailing key field changes nothing and the
#: fixture is invariant under a weaker identity — which is exactly the hole perturbation P1 found in the
#: first version of this file (`TRAPS.md` A2's "one length bin" shape).
_RECS = [
    (0, 100, 300, 0, 0),
    (0, 400, 700, 0, 1),
    (0, 900, 1200, 1, 0),
    (1, 50, 260, 0, 0),
    (1, 50, 260, 1, 0),  # <- differs from the previous ONLY in align_strand
    (1, 800, 1000, 0, 0),
    (1, 800, 1000, 0, 1),  # <- differs from the previous ONLY in sj_strand
]
_CHOICES = np.array([1, 0, 1, 1, 0, 0, 1], np.int64)


def test_L1_a_partition_that_is_the_whole_gets_its_own_choices_back():
    """The degenerate case, and it must be EXACT — not merely close. If the identity fails here it fails
    everywhere, and nothing below is worth reading."""
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [whole], _CHOICES)
    assert ambiguous == 0
    assert np.array_equal(lifted[0], _CHOICES)


def test_L2_a_proper_subset_gets_exactly_the_choices_ITS_records_drew():
    """⭐ Checked PER RECORD against the whole's own choice for the same key. An aggregate check (same
    multiset of choices) would pass with two records swapped, which is precisely the bug that would
    misattribute one origin's mass to another."""
    keep = [_RECS[0], _RECS[2], _RECS[4], _RECS[6]]
    whole = _P(_bank(_RECS))
    part = _P(_bank(keep))
    lifted, ambiguous = lift_choices(whole, [part], _CHOICES)
    assert ambiguous == 0
    want = {r: int(_CHOICES[i]) for i, r in enumerate(sorted(_RECS))}
    for j, r in enumerate(sorted(keep)):
        assert int(lifted[0][j]) == want[r], (j, r, lifted[0][j], want[r])


def test_L3_a_complete_split_consumes_every_choice_EXACTLY_ONCE():
    """⭐⭐ THE IDENTITY. Two disjoint partitions covering the whole must, between them, use each
    whole-library choice once — which is what makes ``Sum(partitions) == whole`` after the drain."""
    a = [_RECS[0], _RECS[3], _RECS[5]]
    b = [_RECS[1], _RECS[2], _RECS[4], _RECS[6]]
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [_P(_bank(a)), _P(_bank(b))], _CHOICES)
    assert ambiguous == 0
    got = [int(x) for arr in lifted for x in arr]
    assert sorted(got) == sorted(int(x) for x in _CHOICES)
    assert len(got) == len(_RECS)


def test_L4_a_record_the_whole_does_NOT_hold_RAISES():
    """⛔ Silence here would be the worst outcome: a partition scanned from a different BAM, or a payload
    pair that drifted, would quietly produce a plausible tally. It must be an error."""
    whole = _P(_bank(_RECS[:3]))
    part = _P(_bank([_RECS[5]]))
    with pytest.raises(ValueError, match="not in the whole library's held bank"):
        lift_choices(whole, [part], _CHOICES[:3])


def test_L5_a_wrong_length_choice_array_RAISES():
    """One choice per held record is the contract; a mismatch means the choices came from a different
    scan and every index below would be off by an unknown amount."""
    whole = _P(_bank(_RECS))
    with pytest.raises(ValueError, match="One choice per held record"):
        lift_choices(whole, [whole], _CHOICES[:2])


def test_L6_duplicate_keys_split_across_partitions_are_COUNTED():
    """⭐ THE ONE AMBIGUITY, and the gate is that it is REPORTED rather than hidden.

    Two records with the same key are identical records (``DeferredFragments``' own guarantee), so no
    partition can know which of them it holds. The deposits are interchangeable; the ORIGIN attribution is
    not. So the count is returned and it bounds the truth error exactly — a caller that drops it is the
    defect this gate exists to make visible.

    ⚠ Note the count is per RECORD of the partition, so a partition holding one of a duplicate pair reports
    1 — it cannot report 'half a fragment', and rounding up is the safe direction for a bound."""
    dup = (0, 100, 300, 0, 0)
    recs = [dup, dup, (0, 500, 700, 0, 0)]
    whole = _P(_bank(recs))
    ch = np.array([0, 1, 1], np.int64)
    # each partition takes ONE of the identical pair — neither can know which
    lifted, ambiguous = lift_choices(whole, [_P(_bank([dup])), _P(_bank([dup]))], ch)
    assert ambiguous == 2, ambiguous
    # ⭐⭐ AND THE QUEUE IS CONSUMED ACROSS PARTITIONS, WHICH IS THE IDENTITY ITSELF. The duplicate group's
    #    two choices are {0, 1}; if the state were per-call both partitions would take run[0] and one entry
    #    would go unused. Perturbation P2 found exactly that in the first version of this file.
    assert sorted([int(lifted[0][0]), int(lifted[1][0])]) == [0, 1]
    # and a partition of DISTINCT records reports none
    l3, a3 = lift_choices(whole, [_P(_bank([recs[2]]))], ch)
    assert a3 == 0
    assert int(l3[0][0]) == 1


def test_L7_the_key_uses_EVERY_field_of_the_record_identity():
    """⭐⭐ KEY COMPLETENESS. ``_RECS`` contains two pairs that differ ONLY in ``align_strand`` and ONLY in
    ``sj_strand``, so a lift that ignored a trailing field would collapse each pair into one queue and hand
    the same choice to both — detected here as a mismatch against the whole's own per-record choice.

    ⛔ Written because perturbation P1 (truncate the key to four fields) passed every other gate: with all
    records differing in (ref, start, end), the fixture was invariant under a weaker identity."""
    whole = _P(_bank(_RECS))
    lifted, ambiguous = lift_choices(whole, [whole], _CHOICES)
    assert ambiguous == 0, "the fixture must have no duplicate FULL keys, or this gate proves nothing"
    assert np.array_equal(lifted[0], _CHOICES)
    # the two near-duplicate pairs must have drawn DIFFERENT choices, or the gate is vacuous (A14)
    order = sorted(_RECS)
    for a, b in (((1, 50, 260, 0, 0), (1, 50, 260, 1, 0)),
                 ((1, 800, 1000, 0, 0), (1, 800, 1000, 0, 1))):
        ca, cb = _CHOICES[order.index(a)], _CHOICES[order.index(b)]
        assert ca != cb, (a, b, ca, cb)


def test_L8_partitions_holding_MORE_copies_of_a_key_than_the_whole_RAISES():
    """⛔ THE OVER-CONSUMPTION GUARD. Two partitions each holding a record whose key is UNIQUE in the whole
    is not a partition of one scan — it is two scans, or a double-counted fragment. Without the guard the
    second partition would silently re-use the first's choice, and ``Sum(partitions) > whole``.

    ⭐ Written because perturbation P5 (delete the guard) passed every other gate: no fixture over-consumed,
    so the guard was never reached (`TRAPS.md` A14 — count the opportunities the change had to fire)."""
    only = (0, 100, 300, 0, 0)
    whole = _P(_bank([only, (0, 500, 700, 0, 0)]))
    ch = np.array([0, 1], np.int64)
    with pytest.raises(ValueError, match="more records with key"):
        lift_choices(whole, [_P(_bank([only])), _P(_bank([only]))], ch)
