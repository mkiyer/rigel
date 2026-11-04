"""The accumulator's fragment length ``L`` and its deposit geometry, proved against an independent oracle.

``L`` is the tool's one definition of fragment length, so it is checked against a different algorithm, not a different
spelling of the same one: the molecule's covered bases are re-derived by integer set arithmetic, ``set(range(start,
end))`` minus each intron's range, with no sorting, merging, cursor walk or searchsorted. Malformed intron lists
(overlapping, nested, abutting, duplicated, zero-length, out of range) that production handles by explicit logic, the
oracle handles by construction, because set subtraction is idempotent and order-free. ``L``, the crossed boundaries and
the containing region all follow from that one set, so the tests prove the geometry as well as the length: whatever
counts toward ``L`` must also count as coverage for crossing, or the density estimator is biased.

The space is enumerated exhaustively up to two introns over a 12 bp reference and sampled at BAM-scale coordinates with
up to four. The fixture carries no annotated sj; spliced routing is ``test_accumulator_spec.py``'s.
``test_accumulator_native_parity.py`` holds the C++ byte-identical to this reference, so proving the reference proves
the C++.
"""

from __future__ import annotations

import itertools

import numpy as np

import pytest

from ._accumulator_reference import Accumulator, DepositOutcome, GapHypothesis, Partition


# --- the oracle -------------------------------------------------------------------------------------


def covered_bases(ref_len: int, start: int, end: int, introns) -> set[int]:
    """The molecule's covered genomic bases, by set arithmetic. A fragment overhanging a reference end is clipped, and
    ``L`` is the clipped length."""
    s, e = max(int(start), 0), min(int(end), ref_len)
    if e <= s:
        return set()
    out = set(range(s, e))
    for a, b in introns:
        out -= set(range(max(int(a), s), min(int(b), e)))
    return out


def oracle_deposits(region_bounds, ref_len: int, start: int, end: int, introns, spliced: bool):
    """Every population the deposit rule credits, derived from the covered set alone.

    Returns ``(L, crossed_boundaries, contained_region)`` with local indices, or ``None`` for a fragment that deposits
    nothing. ``region_bounds`` are this reference's region-bound positions. The oracle derives only what an assertion
    reads.
    """
    covered = covered_bases(ref_len, start, end, introns)
    if not covered:
        return None
    length = len(covered)
    # A boundary at region-bound position p is crossed iff the molecule holds the bases on both sides of it. Reported as
    # boundary indices: region bound ``i`` is boundary ``i − 1`` (the deposit writes ``boundary_base + boundary - 1``),
    # because bound 0 is the reference start and owns no interior boundary. The oracle re-derives that offset rather
    # than importing it, so an offset error shows.
    crossed = {i - 1 for i, p in enumerate(region_bounds) if (p - 1) in covered and p in covered}

    def region_of(pos):
        return min(
            max(int(np.searchsorted(region_bounds, pos, side="right")) - 1, 0),
            len(region_bounds) - 2,
        )

    lo, hi = region_of(min(covered)), region_of(max(covered))
    contained = lo if (not spliced and lo == hi) else None
    return length, crossed, contained


# --- the fixture ------------------------------------------------------------------------------------

_REGION_BOUNDS = np.array(
    [0, 3, 4, 9, 12], dtype=np.int64
)  # 4 regions: widths 3, 1, 5, 3 — includes a 1 bp region
_REF_LEN = 12


def _acc(max_fragment_length: int = 10_000) -> Accumulator:
    n_regions = _REGION_BOUNDS.size - 1
    part = Partition(
        region_bounds=_REGION_BOUNDS.copy(),
        ref_region_bound_offsets=np.array([0, _REGION_BOUNDS.size], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_regions - 1], dtype=np.int64),
        sj_offsets=np.zeros(_REGION_BOUNDS.size + 1, dtype=np.int64),
        sj_boundary_right=np.zeros(0, dtype=np.int64),
        sj_strand=np.zeros(0, dtype=np.int8),
        region_types=np.zeros(n_regions, dtype=np.int8),
    )
    return Accumulator(part, max_fragment_length=max_fragment_length)


def _observed(acc: Accumulator):
    t = acc.tally
    return (
        {
            i
            for i in range(t.boundary_unspliced_count.shape[0])
            if t.boundary_unspliced_count[i].sum()
        },
        next(
            (
                i
                for i in range(t.region_contained_count.shape[0])
                if t.region_contained_count[i].sum()
            ),
            None,
        ),
    )


def _check_one(start: int, end: int, introns) -> None:
    """Deposit one fragment and assert every population against the set oracle."""
    acc = _acc()
    outcome = acc.deposit(0, start, end, observed_introns=introns)
    want = oracle_deposits(_REGION_BOUNDS, _REF_LEN, start, end, introns, spliced=False)

    if want is None:
        assert outcome is DepositOutcome.EMPTY, (
            f"[{start},{end}) introns={introns}: oracle says the molecule is empty, "
            f"the accumulator returned {outcome}"
        )
        return
    length, crossed, contained = want
    assert outcome is DepositOutcome.DEPOSITED, f"[{start},{end}) introns={introns} -> {outcome}"

    got_crossed, got_contained = _observed(acc)
    ctx = f"[{start},{end}) introns={introns}"

    # L is read back from the deposit itself, the number that reached the tally.
    t = acc.tally
    n_cross = int(t.boundary_unspliced_count.sum())
    if n_cross:
        assert _close(
            float(t.boundary_unspliced_inv_length_sum.sum()), n_cross / (length - 1), n_cross
        ), f"{ctx}: the deposited L disagrees with the oracle's {length}"
    if contained is not None:
        assert int(t.deposited_lengths[length]) == 1, (
            f"{ctx}: the deposited L disagrees with the oracle's {length}"
        )

    assert got_crossed == crossed, f"{ctx}: crossed boundaries {got_crossed} != oracle {crossed}"
    assert got_contained == contained, f"{ctx}: contained {got_contained} != oracle {contained}"
    # the start-count invariant: exactly one, at the region holding the first covered base
    assert int(t.region_start_count.sum()) == 1, ctx


# --- Exhaustive: every configuration with up to two introns ------------------------------------------


def _all_intervals(n: int):
    return [(a, b) for a in range(n + 1) for b in range(a, n + 1)]


@pytest.mark.parametrize("n_introns", [0, 1, 2])
def test_exhaustive_against_the_set_oracle(n_introns):
    """Every fragment and every intron list with ``n_introns`` introns over a 12 bp reference: the complete space, so
    zero-length, outside, nested, overlapping, abutting and duplicated introns all occur by enumeration.
    """
    spans = [(s, e) for s, e in _all_intervals(_REF_LEN) if e > s]
    intron_lists = list(
        itertools.combinations_with_replacement(_all_intervals(_REF_LEN), n_introns)
    )
    checked = 0
    for start, end in spans:
        for introns in intron_lists:
            _check_one(start, end, list(introns))
            checked += 1
    assert checked > 0
    print(f"\n  n_introns={n_introns}: {checked:,} configurations checked against the set oracle")


# --- Randomised: realistic coordinates, up to four introns -------------------------------------------


def test_randomised_at_realistic_scale():
    """The same oracle at coordinates a real BAM produces, with up to 4 introns. Only this scale can produce a 300 bp
    molecule spanning a 10 kb intron, the shape that shows the length limit applies to ``L`` and never to the span. The
    seed is fixed so the proof is reproducible.
    """
    rng = np.random.default_rng(20260731)
    region_bounds = np.array([0, 137, 138, 400, 401, 1200, 5000, 5001, 9000], dtype=np.int64)
    ref_len = 9000
    n_regions = region_bounds.size - 1
    part = Partition(
        region_bounds=region_bounds.copy(),
        ref_region_bound_offsets=np.array([0, region_bounds.size], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_regions - 1], dtype=np.int64),
        sj_offsets=np.zeros(region_bounds.size + 1, dtype=np.int64),
        sj_boundary_right=np.zeros(0, dtype=np.int64),
        sj_strand=np.zeros(0, dtype=np.int8),
        region_types=np.zeros(n_regions, dtype=np.int8),
    )
    n_checked = 0
    for _ in range(4000):
        start = int(rng.integers(-50, ref_len + 50))
        end = start + int(rng.integers(1, 3000))
        introns = []
        for _ in range(int(rng.integers(0, 5))):
            a = int(rng.integers(start - 20, end + 20))
            introns.append((a, a + int(rng.integers(0, 900))))
        acc = Accumulator(part, max_fragment_length=10**9)
        outcome = acc.deposit(0, start, end, observed_introns=introns)
        want = oracle_deposits(region_bounds, ref_len, start, end, introns, spliced=False)
        ctx = f"[{start},{end}) introns={introns}"
        if want is None:
            assert outcome is DepositOutcome.EMPTY, ctx
            continue
        length, crossed, contained = want
        assert outcome is DepositOutcome.DEPOSITED, ctx
        got_crossed, got_contained = _observed(acc)
        assert got_crossed == crossed, f"{ctx}: {got_crossed} != {crossed}"
        assert got_contained == contained, ctx
        t = acc.tally
        n_cross = int(t.boundary_unspliced_count.sum())
        if n_cross:
            assert _close(
                float(t.boundary_unspliced_inv_length_sum.sum()), n_cross / (length - 1), n_cross
            ), f"{ctx}: deposited L != oracle {length}"
        n_checked += 1
    assert n_checked > 3000, (
        f"only {n_checked} fragments actually deposited — the sweep is degenerate"
    )


# --- Coverage and length share one set ----------------------------------------------------------------


def test_L_equals_the_covered_base_count_and_crossings_use_THAT_SAME_set():
    """Whatever counts toward ``L`` must also count as coverage for crossing, or the density estimator is biased. A mate
    gap counts toward ``L`` and crosses boundaries; an intron does neither. Both are asserted from one set.
    """
    # a paired-end molecule with an unsequenced mate gap: the gap is part of the molecule
    acc = _acc()
    acc.deposit(0, 1, 11, observed_introns=[])
    assert int(acc.tally.region_contained_count.sum()) == 0
    got_crossed, _ = _observed(acc)
    assert got_crossed == {0, 1, 2}, (
        "the mate gap must carry the molecule across every interior boundary"
    )

    # the same span with the middle excised as an intron: L shrinks and the boundaries under it go uncrossed
    acc2 = _acc()
    acc2.deposit(0, 1, 11, observed_introns=[(3, 9)])
    got2, _ = _observed(acc2)
    assert got2 == set(), (
        "an intron must not carry the molecule across the boundaries it splices over"
    )
    assert len(covered_bases(_REF_LEN, 1, 11, [(3, 9)])) == 4


# --- The unconditional length histogram and its invariant -------------------------------------------


def test_deposited_lengths_bins_every_accepted_fragment_exactly_once():
    """``Σ deposited_lengths == Σ region_start_count == qc.deposited``.

    Three counters over one population. The start-count invariant says every fragment was located in space; this one
    says every fragment was binned by length. The histogram anchors every fragment-length model in the tool, so it must
    not run one fragment short.
    """
    acc = _acc(max_fragment_length=64)
    lengths = []
    for start, end, introns in (
        (0, 5, []),  # contained
        (1, 11, []),  # crossing, mate gap
        (1, 11, [(3, 9)]),  # crossing nothing: an intron swallows the interior
        (2, 12, [(4, 5)]),
        (-3, 15, []),  # clipped to the reference on both sides
        (0, 12, [(0, 3)]),  # a leading intron: the path starts at 3
    ):
        assert acc.deposit(0, start, end, observed_introns=introns) is DepositOutcome.DEPOSITED
        lengths.append(len(covered_bases(_REF_LEN, start, end, introns)))
    t = acc.tally
    assert int(t.deposited_lengths.sum()) == len(lengths)
    assert int(t.deposited_lengths.sum()) == int(t.region_start_count.sum())
    assert int(t.deposited_lengths.sum()) == t.qc[DepositOutcome.DEPOSITED.value]
    # and it is binned at L, the oracle's L rather than the accumulator's own
    expected = np.zeros_like(t.deposited_lengths)
    for length in lengths:
        expected[length] += 1
    np.testing.assert_array_equal(t.deposited_lengths, expected)


def test_a_REJECTED_fragment_is_not_binned():
    """The histogram is unconditional given deposit. A fragment over the length limit, or one whose path is
    undetermined, deposits nothing and bins nothing; each is counted in ``qc`` instead. So the histogram describes
    exactly the population the pure pools are drawn from; an anchor over a wider population would mismatch the pools'
    frame.
    """
    acc = _acc(max_fragment_length=4)
    assert acc.deposit(0, 0, 12, observed_introns=[]) is DepositOutcome.TOO_LONG
    assert int(acc.tally.deposited_lengths.sum()) == 0
    assert acc.tally.qc[DepositOutcome.TOO_LONG.value] == 1

    # A deferred fragment: two hypotheses about its unsequenced gap survive, so its L is not one number and there is
    # nothing to bin. The deferral is the accumulator's own verdict, so the test hands it a hypothesis set with two
    # answers.
    acc2 = _acc()
    outcome = acc2.deposit(0, 0, 5, hypotheses=(GapHypothesis(((1, 3),)), GapHypothesis()))
    assert outcome is DepositOutcome.DEFERRED
    assert int(acc2.tally.deposited_lengths.sum()) == 0
    assert len(acc2.tally.deferred) == 1, "and it is HELD, not dropped"

    acc3 = _acc()
    assert acc3.deposit(0, 4, 4, observed_introns=[]) is DepositOutcome.EMPTY
    assert int(acc3.tally.deposited_lengths.sum()) == 0


def test_the_unconditional_histogram_is_a_SUPERSET_of_the_pure_pools():
    """Every pooled fragment is also binned here, at the same length: the pools are conditioned subsets of this
    population. The superset is strict (an exonic contained fragment or a multi-boundary crossing enters no pool, since
    an impure pool is worse than a missing one), which is why the pools cannot anchor themselves.
    """
    acc = _acc(max_fragment_length=64)
    for start, end in ((0, 2), (0, 5), (1, 11), (3, 4), (5, 9), (2, 12)):
        acc.deposit(0, start, end, observed_introns=[])
    t = acc.tally
    pooled = t.pool_lengths.sum(axis=0)
    assert pooled.sum() > 0, "the fixture must actually populate some pools or this proves nothing"
    assert np.all(t.deposited_lengths >= pooled), (
        "a fragment in a pure pool must also be binned in the unconditional histogram, at the same L"
    )
    assert int(t.deposited_lengths.sum()) == t.qc[DepositOutcome.DEPOSITED.value]


#: Fractions are float64, so `n` round-to-nearest additions differ from the exact answer by at most `n` ulp.
EPS = float(np.finfo(np.float64).eps)


def _close(got: float, want: float, deposits: int) -> bool:
    return abs(got - want) <= max(abs(want), 1.0) * deposits * EPS
