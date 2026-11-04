"""The conserved mass: the deposit rule's own tests, on the executable specification.

``boundary_unspliced_count`` is ``+1`` on every boundary a fragment crosses, so a fragment books ``max(K, 1)`` of them
and a sum over objects counts object incidences, not fragments. The prior the EM reads is a fragment count, and deriving
one from a density over-calls it under capture; the mass bank carries the count directly, one fragment depositing
exactly one. Conservation alone cannot pin the rule (a ``1/K`` share conserves too), so this module makes four separate
claims:

===  ==========================================================================================
 1   a fragment's deposits sum to one, within the float64 representation
 2   the mass equals the per-base attribution, which ``1/K`` cannot reproduce because a base does
     not know how many boundaries its fragment crossed; this claim pins the rule
 3   a spliced fragment's mass reaches the spliced bank and nothing else
 4   every crossed object shares the one fragment, sj boundaries included
===  ==========================================================================================

Claim 2 restates the rule independently: the specification uses slices and two ``searchsorted`` calls, while the
per-base walk asks of every base which region holds it and which of that region's bounds lie strictly inside the
fragment. Spec: ``_accumulator_reference.py``; C++ parity: ``test_accumulator_native_parity.py``.
"""

from __future__ import annotations

import bisect
from fractions import Fraction

import numpy as np

from rigel.types import Strand

from ._accumulator_reference import Accumulator, DepositOutcome, Partition

#: The mass banks are float64 fractions. A fragment makes a bounded number of deposits, each rounded to nearest, so its
#: total differs from 1 by at most `deposits/2` ulp, with `np.finfo(np.float64).eps` the machine's own.
EPS = float(np.finfo(np.float64).eps)


def budget(deposits: int) -> float:
    """The representation error of `deposits` round-to-nearest additions summing to about 1."""
    return deposits * EPS


#: Boundaries 1 and 5 have both flanks wider than every fragment below: the control, where a crossing fragment can cross
#: only that one boundary. Boundaries 2, 3 and 4 bound a cluster of 50 bp regions, so one 300 bp fragment crosses three
#: at once; under capture such clusters hold most gDNA crossings.
REGION_BOUNDS = [0, 400, 1000, 1050, 1100, 1500, 1900]
LENGTHS = (60, 150, 300)
MAX_LENGTH = 1000


def _partition(region_bounds=REGION_BOUNDS) -> Partition:
    return Partition.from_region_bounds(
        [region_bounds], region_types=[[0] * (len(region_bounds) - 1)]
    )


def _mass_by_base(region_bounds, start: int, end: int) -> dict[int, Fraction]:
    """The per-boundary mass of one unspliced fragment, attributed base by base.

    Each base is split equally between the bounds of its region that lie strictly inside the fragment; divided by the
    fragment's length, that is the mass, reached without forming a slice. ``1/K`` is not expressible this way, which is
    why this attribution pins the rule.
    """
    length = end - start
    out: dict[int, Fraction] = {}
    for position in range(start, end):
        region = bisect.bisect_right(region_bounds, position) - 1
        boundaries = [i for i in (region, region + 1) if start < region_bounds[i] < end]
        for boundary in boundaries:
            out[boundary] = out.get(boundary, Fraction(0)) + Fraction(1, length * len(boundaries))
    return out


def _enumerate(partition, region_bounds=REGION_BOUNDS, lengths=LENGTHS, perturb=None):
    """Every fragment of every length at every start position. Returns ``(accumulator, n)``."""
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    n = 0
    for length in lengths:
        for start in range(region_bounds[0], region_bounds[-1] - length + 1):
            outcome = accumulator.deposit(0, start, start + length)
            assert outcome is DepositOutcome.DEPOSITED, f"[{start},{start + length}) was {outcome}"
            n += 1
    return accumulator, n


# ---------------------------------------------------------------------------
# claim 1: conservation, exactly
# ---------------------------------------------------------------------------


def _exact_fragment_mass(partition, region_bounds, start, end) -> tuple[Fraction, int]:
    """``(Σ mass, contained)`` for one fragment, as the float64 bank holds it. The accumulator is fresh, so the sum over
    boundaries after one deposit is that fragment's deposit.
    """
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(0, start, end)
    assert outcome is DepositOutcome.DEPOSITED
    total = float(np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64).sum())
    contained = int(np.asarray(accumulator.tally.region_contained_count, np.int64).sum())
    return total, contained


def test_the_RULE_conserves_exactly_before_any_quantisation():
    """The law in exact rational arithmetic: ``Σ share == 1`` per crossing fragment.

    Every unspliced path is one segment. If it stays inside a region it is ``contained``; otherwise it has ``n >= 2``
    slices, the first crossing right only, the last left only, every interior one both. Every slice has at least one
    bounding boundary, so ``Σ slice_len / L = 1`` exactly and the "no boundary here" guard drops nothing. This is stated
    on the rule, through the per-base attribution, not on the bank, which cannot represent a third exactly; the next
    test checks the bank, so a rounding residue and a modelling error stay distinguishable.
    """
    for length in LENGTHS:
        for start in range(REGION_BOUNDS[0], REGION_BOUNDS[-1] - length + 1):
            by_base = _mass_by_base(REGION_BOUNDS, start, start + length)
            if not by_base:
                continue  # contained: the region bank holds the whole fragment
            assert sum(by_base.values()) == 1, (
                f"fragment [{start},{start + length}) shares sum to {sum(by_base.values())}"
            )


def test_the_BANK_conserves_to_within_the_REPRESENTATION_and_no_more():
    """The same law as stored, within the price of float64. A fragment crossing ``K`` boundaries makes ``2K`` deposits
    (``K + 1`` slices: the ends deposit once, each interior one twice), each rounded to nearest, so the bank can differ
    from 1 by at most ``2K`` half-ulp. The worst drift is also asserted as a range: zero would mean the representation
    changed, more than the budget that the deposits are not round-to-nearest.
    """
    partition = _partition()
    worst = 0.0
    for length in LENGTHS:
        for start in range(REGION_BOUNDS[0], REGION_BOUNDS[-1] - length + 1):
            total, contained = _exact_fragment_mass(partition, REGION_BOUNDS, start, start + length)
            crossings = len(_mass_by_base(REGION_BOUNDS, start, start + length))
            allowed = budget(2 * max(crossings, 1))
            deviation = abs(total + contained - 1)
            worst = max(worst, deviation)
            assert deviation <= allowed, (
                f"fragment [{start},{start + length}) deposited {total} + {contained} "
                f"contained, off by {deviation:.3e} against a budget of {allowed:.3e}"
            )
    # A few ulp at worst, asserted as a range for the reason the docstring gives.
    assert 0 < worst <= budget(6), (
        f"the worst per-fragment representation drift is {worst:.3e}; if it is now zero the "
        f"representation changed, and if it is larger the deposits are not rounding as assumed"
    )


def test_a_CONTAINED_fragment_deposits_NO_mass_at_all():
    """Its whole fragment is already ``region_contained_count``; mass at a boundary too would count it twice."""
    partition = _partition()
    total, contained = _exact_fragment_mass(
        partition, REGION_BOUNDS, 100, 160
    )  # well inside region 0
    assert contained == 1
    assert total == 0


def test_the_two_populations_together_are_EXHAUSTIVE_for_an_unspliced_fragment():
    """A contiguous interval either stays inside one region or crosses at least one boundary, so over the enumeration
    ``Σ contained + Σ mass`` is the fragment count. The budget is one rounding per crossing, as in the per-fragment
    test.
    """
    accumulator, n = _enumerate(_partition())
    mass = float(np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64).sum())
    contained = int(np.asarray(accumulator.tally.region_contained_count, np.int64).sum())
    crossings = int(np.asarray(accumulator.tally.boundary_unspliced_count, np.int64).sum())
    assert abs(mass + contained - n) <= budget(crossings), (
        f"{n} fragments deposited {float(mass)} of mass plus {contained} contained"
    )
    # It is a fragment count to well under one fragment, which is what the bank is for.
    assert abs(float(mass) + contained - n) < 1e-6


# ---------------------------------------------------------------------------
# claim 2: the value, which separates this rule from 1/K
# ---------------------------------------------------------------------------


def test_the_mass_is_the_PER_BASE_attribution():
    """The mass is coverage-weighted, not ``1/K``. Both rules conserve, so no conservation test can tell them apart, and
    at ``K == 1`` they coincide; the per-base walk splits each base between its own region's bounds, which no ``1/K``
    rule can match because a base does not know how many boundaries the fragment crossed. The tolerance is the
    representation: each of a boundary's ``<= 2 * count`` deposits rounds to nearest, so the bank can differ by at most
    that many half-ulp.
    """
    partition = _partition()
    accumulator, _n = _enumerate(partition)
    expected: dict[int, Fraction] = {}
    for length in LENGTHS:
        for start in range(REGION_BOUNDS[0], REGION_BOUNDS[-1] - length + 1):
            for boundary, mass in _mass_by_base(REGION_BOUNDS, start, start + length).items():
                expected[boundary] = expected.get(boundary, Fraction(0)) + mass

    bank = np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64)
    count = np.asarray(accumulator.tally.boundary_unspliced_count, np.int64).sum(axis=1)
    for boundary in range(1, len(REGION_BOUNDS) - 1):
        boundary = boundary - 1
        got = float(bank[boundary])
        want = expected.get(boundary, Fraction(0))
        budget = float(count[boundary])
        assert abs(got - want) <= budget, (
            f"boundary {boundary} @ {REGION_BOUNDS[boundary]}: bank {float(got)} vs per-base {float(want)}, "
            f"budget {float(budget):.3e}"
        )


def test_the_mass_EQUALS_the_count_where_both_flanks_exceed_every_fragment():
    """Where both flanking regions are longer than every fragment, a crossing fragment crosses only that boundary, so
    its whole 1.0 lands there and ``mass == count`` exactly. Where regions are shorter, one fragment spreads over
    several boundaries and ``mass`` falls below ``count``; that gap is the K-inflation that keeps a per-object sum from
    being a fragment count under capture.
    """
    partition = _partition()
    accumulator, _n = _enumerate(partition)
    bank = np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64)
    count = np.asarray(accumulator.tally.boundary_unspliced_count, np.int64).sum(axis=1)

    wide = [
        boundary
        for boundary in range(1, len(REGION_BOUNDS) - 1)
        if min(
            REGION_BOUNDS[boundary] - REGION_BOUNDS[boundary - 1],
            REGION_BOUNDS[boundary + 1] - REGION_BOUNDS[boundary],
        )
        > max(LENGTHS)
    ]
    narrow = [boundary for boundary in range(1, len(REGION_BOUNDS) - 1) if boundary not in wide]
    assert wide and narrow, "the fixture must contain both cases or this gate tests one of them"

    for boundary in wide:
        assert float(bank[boundary - 1]) == float(count[boundary - 1]), (
            f"boundary {boundary} has flanks wider than every fragment, so mass must equal count exactly"
        )
    for boundary in narrow:
        assert float(bank[boundary - 1]) < float(count[boundary - 1]), (
            f"boundary {boundary} sits in the short cluster, so its count must over-state the fragments"
        )


# ---------------------------------------------------------------------------
# claim 3 — routing
# ---------------------------------------------------------------------------


def test_a_SPLICED_fragment_s_mass_reaches_the_SPLICED_bank_ALONE():
    """A spliced fragment is certified RNA at that boundary (gDNA cannot splice), so its mass must never enter the
    unspliced competition the prior arbitrates. The geometry is the staggered case the second bank exists for: the
    fragment splices over one exon boundary, then runs contiguously across another isoform's boundary inside its second
    block.
    """
    region_bounds = [0, 1000, 2000, 9000, 9050, 10000]
    partition = Partition.from_region_bounds(
        [region_bounds],
        region_types=[[0, 2, 1, 2, 2]],
        sj=[(0, 2000, 9000, Strand.POS)],
    )
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(
        0, 1900, 9200, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
    )
    assert outcome is DepositOutcome.DEPOSITED

    tally = accumulator.tally
    assert int(np.asarray(tally.sj_count, np.int64).sum()) == 1, "the sj must be credited"
    assert float(np.asarray(tally.boundary_unspliced_mass, np.float64).sum()) == 0.0, (
        "certified RNA leaked into the bank being deconvolved"
    )
    spliced = float(np.asarray(tally.boundary_spliced_mass, np.float64).sum())
    sj = float(np.asarray(tally.sj_mass, np.float64).sum())
    assert spliced > 0, "the contiguous crossing at 9050 deposited nothing"
    # A partial, exactly as the rule specifies: the boundary takes only the bases adjacent to it and the sj
    # takes the rest. Derived slice by slice from the rule::
    #
    #   L = (2000-1900) + (9200-9000) = 300
    #   block 0 [1900,2000)  crosses nothing; its only boundary is the sj   -> 100/300 sj
    #   block 1 [9000,9200)  crosses the boundary at 9050, giving two slices:
    #       [9000,9050)  bounded by the sj (block start) and the boundary       ->  25/300 each
    #       [9050,9200)  bounded by the boundary alone (block end is the frag end)    -> 150/300 boundary
    #
    # 200/300 for the boundary would mean a boundary-crossing block gives all its bases to boundaries,
    # leaving the sj that bounds it with nothing.
    assert spliced < 1
    want_boundary = Fraction(25, 300) + Fraction(150, 300)
    want_sj = Fraction(100, 300) + Fraction(25, 300)
    assert abs(spliced - want_boundary) <= budget(2)
    assert abs(sj - want_sj) <= budget(2)
    # Together they are the whole fragment.
    assert abs(spliced + sj - 1) <= budget(4)


def test_the_spliced_mass_is_a_PARTIAL_and_never_a_conservation_ledger():
    """Read as "the number of spliced fragments here", this bank would be wrong, and a consumer could read it that way,
    so the partial is pinned."""
    region_bounds = [0, 1000, 2000, 9000, 9050, 10000]
    partition = Partition.from_region_bounds(
        [region_bounds], region_types=[[0, 2, 1, 2, 2]], sj=[(0, 2000, 9000, Strand.POS)]
    )
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    n = 0
    for start in range(1900, 1990):
        outcome = accumulator.deposit(
            0, start, 9200, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
        )
        assert outcome is DepositOutcome.DEPOSITED
        n += 1
    total = float(np.asarray(accumulator.tally.boundary_spliced_mass, np.float64).sum())
    assert 0 < total < n


# ---------------------------------------------------------------------------
# claim 4: the law over all three axes, which makes a library count computable
# ---------------------------------------------------------------------------

#: The staggered geometry claim 3 uses: one annotated sj, and a short region (9000–9050) whose far boundary a fragment
#: may or may not reach, so one partition carries both a spliced fragment that crosses a boundary and one that crosses
#: none.
SPLICED_REGION_BOUNDS = [0, 1000, 2000, 9000, 9050, 10000]


def _spliced_partition() -> Partition:
    return Partition.from_region_bounds(
        [SPLICED_REGION_BOUNDS], region_types=[[0, 2, 1, 2, 2]], sj=[(0, 2000, 9000, Strand.POS)]
    )


def _total_deposited_mass(tally) -> Fraction:
    """Every conserved-mass bank summed: unspliced boundaries, spliced boundaries, and sj. ``sj_mass`` is read with
    :func:`getattr` so a missing bank reports the law failing by a whole fragment rather than an ``AttributeError``.
    """
    total = Fraction(0)
    for name in ("boundary_unspliced_mass", "boundary_spliced_mass", "sj_mass"):
        bank = getattr(tally, name, None)
        if bank is not None:
            total += float(np.asarray(bank, np.float64).sum())
    return total


def test_a_spliced_fragment_crossing_NO_BOUNDARY_still_deposits_a_WHOLE_FRAGMENT_of_mass():
    """Both blocks lie inside one region each: the fragment is not ``contained`` (its path spans a sj) and crosses no
    boundary. Without a sj mass it would deposit nothing anywhere while ``sj_count`` credits it. It is all RNA, since
    gDNA cannot splice, which is why the gDNA side's conserved count is complete without this bank and the RNA side's is
    not.
    """
    accumulator = Accumulator(_spliced_partition(), max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(
        0, 1900, 9040, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
    )
    assert outcome is DepositOutcome.DEPOSITED

    tally = accumulator.tally
    assert int(np.asarray(tally.sj_count, np.int64).sum()) == 1, "the sj must be credited"
    assert int(np.asarray(tally.region_contained_count, np.int64).sum()) == 0, (
        "a path spanning a sj is not contained in ONE region"
    )
    assert _total_deposited_mass(tally) == 1, (
        "a deposited fragment must place exactly one fragment of mass somewhere; this one placed "
        "none, so no sum over conserved banks can count it"
    )


#: Two sj, and the middle exon is a whole region. A one-sj fixture cannot exercise the equal-share rule: every block it
#: produces is bounded on one side only, so a rule that forgot to share would pass.
TWO_SJ_REGION_BOUNDS = [0, 1000, 2000, 5000, 5100, 8000, 9000, 10000]


def _two_sj_partition() -> Partition:
    return Partition.from_region_bounds(
        [TWO_SJ_REGION_BOUNDS],
        region_types=[[0, 2, 1, 2, 1, 2, 0]],
        sj=[(0, 2000, 5000, Strand.POS), (0, 5100, 8000, Strand.POS)],
    )


#: The three blocks of ``[1950, 8050)`` spliced over both introns: 50 bp, the whole 100 bp middle exon, then 50 bp. None
#: contains an interior boundary, so all three are sj-bounded and the middle one is bounded twice.
TWO_SJ_FRAGMENT = (1950, 8050, [(2000, 5000), (5100, 8000)])


def test_a_block_bounded_by_TWO_sj_SPLITS_its_mass_between_them():
    """The test that pins the sj rule, as claim 2 pins the boundary rule. Conservation cannot tell sharing a
    doubly-bounded block from giving it whole to each bound, so this states where the mass went, per sj, in exact
    rational arithmetic::

        block 0  [1950,2000)   50/200 -> j0            (fragment start on the left: one bound)
        block 1  [5000,5100)  100/200 -> j0 and j1     (a sj at both ends: shared)
        block 2  [8000,8050)   50/200 -> j1            (fragment end on the right: one bound)

    so ``j0 == j1 == 1/2`` and the two sum to exactly one fragment.
    """
    start, end, introns = TWO_SJ_FRAGMENT
    accumulator = Accumulator(_two_sj_partition(), max_fragment_length=MAX_LENGTH)
    assert (
        accumulator.deposit(0, start, end, observed_introns=introns, sj_strand=Strand.POS)
        is DepositOutcome.DEPOSITED
    )
    # Summed over the strand columns. `sj_mass` is per strand for artifact detection, and conservation is
    # strand-agnostic (a fragment's 1.0 is shared across the objects it touched whatever its strand), so every test here
    # reads the total.
    bank = np.asarray(accumulator.tally.sj_mass, np.float64).sum(axis=1)
    assert bank.size == 2, "the fixture must carry both sj or this gate tests one of them"
    got = [float(v) for v in bank]

    length = (2000 - 1950) + (5100 - 5000) + (8050 - 8000)
    want = [
        Fraction(2000 - 1950, length) + Fraction(5100 - 5000, length) / 2,
        Fraction(5100 - 5000, length) / 2 + Fraction(8050 - 8000, length),
    ]
    for jid, (g, w) in enumerate(zip(got, want)):
        assert abs(g - w) <= budget(1), (
            f"sj {jid} holds {float(g)} but its bases give it {float(w)}"
        )
    assert sum(got) == 1, f"the three blocks sum to {float(sum(got))}, not one fragment"


def test_a_fragment_crossing_BOTH_sj_AND_boundaries_gives_EVERY_crossed_object_a_share():
    """A spliced fragment may cross several sj and several boundaries, and all of them share the one fragment.

    Conservation alone does not imply sharing: a rule giving a boundary-crossing block's bases entirely to boundaries
    leaves a sj flanked by two such blocks with zero while ``sj_count`` credits it, and the fragment still sums to 1.0.
    The fixture is that case::

        block 0  [1900,2500)  crosses boundary @2000     and is bounded on the right by sj 0
        block 1  [5000,5200)  crosses boundary @5100     bounded by sj 0 left, sj 1 right
        block 2  [8000,8100)  crosses nothing            bounded by sj 1 on the left

    so sj 0 is flanked by two boundary-crossing blocks.
    """
    region_bounds = [0, 1000, 2000, 2500, 5000, 5100, 5200, 8000, 9000, 10000]
    partition = Partition.from_region_bounds(
        [region_bounds],
        region_types=[[0, 2, 2, 1, 2, 2, 1, 2, 0]],
        sj=[(0, 2500, 5000, Strand.POS), (0, 5200, 8000, Strand.POS)],
    )
    accumulator = Accumulator(partition, max_fragment_length=2000)
    assert (
        accumulator.deposit(
            0, 1900, 8100, observed_introns=[(2500, 5000), (5200, 8000)], sj_strand=Strand.POS
        )
        is DepositOutcome.DEPOSITED
    )
    tally = accumulator.tally

    crossed_boundaries = np.flatnonzero(
        np.asarray(tally.boundary_spliced_count, np.int64).sum(axis=1)
    )
    crossed_sj = np.flatnonzero(np.asarray(tally.sj_count, np.int64).sum(axis=1))
    assert crossed_boundaries.size == 2, (
        "the fixture must cross two boundaries or it tests one case"
    )
    assert crossed_sj.size == 2, "the fixture must use both sj"

    boundary_mass = np.asarray(tally.boundary_spliced_mass, np.float64)
    sj_mass = np.asarray(tally.sj_mass, np.float64).sum(axis=1)  # strand-agnostic, see above
    for boundary in crossed_boundaries:
        assert boundary_mass[boundary] > 0, f"boundary {boundary} was crossed but holds no mass"
    for jid in crossed_sj:
        assert sj_mass[jid] > 0, (
            f"sj {jid} was crossed but holds NO mass — the fragment's 1.0 is conserved yet this "
            f"object got none of it, which is not a sharing"
        )

    total = _total_deposited_mass(tally)
    deposits = int(len(crossed_boundaries) * 2 + len(crossed_sj) * 2)
    assert abs(total - 1) <= budget(deposits), f"the fragment deposited {float(total)}, not one"


def test_a_sj_claims_at_BOTH_its_positions_and_never_ALSO_as_a_contiguous_crossing():
    """Two isoforms sharing a donor and differing by 50 bp at the acceptor::

        TA+  exons (1000,2000) (5000,10000)
        TB+  exons (1000,2000) (5050,10000)

    A fragment with blocks ``(1800,2000)`` and ``(5000,5200)`` (``L = 400``) splices TA's ``(2000,5000)`` and then runs
    contiguously across TB's acceptor boundary at 5050, touching three boundary positions, two of them the sj's:

        ===============  ======  ==================================================  ==========
        bases            n       bounded by                                          share
        ===============  ======  ==================================================  ==========
        ``[1800,2000)``  200     the sj's donor side (fragment start on left)        200 sj
        ``[5000,5050)``   50     the sj's acceptor side, and the boundary @5050      25 / 25
        ``[5050,5200)``  150     the boundary @5050 (fragment end on right)          150 boundary
        ===============  ======  ==================================================  ==========

    so the sj holds ``225/400`` and the boundary ``175/400``. A rule giving a boundary-crossing block's bases entirely
    to boundaries scores 200/200, because the sj never collects on its acceptor side; both sum to 1.0, so only this test
    tells them apart.

    A sj must not also register as a contiguous crossing, which would count one traversal twice. Structurally it cannot
    (a boundary is crossed iff it lies strictly inside a block, and a sj's donor and acceptor are block endpoints), and
    the test asserts it.
    """
    region_bounds = [0, 1000, 2000, 5000, 5050, 10000, 12000]
    partition = Partition.from_region_bounds(
        [region_bounds],
        region_types=[[0, 2, 1, 2, 2, 0]],
        sj=[(0, 2000, 5000, Strand.POS), (0, 2000, 5050, Strand.POS)],
    )
    accumulator = Accumulator(partition, max_fragment_length=2000)
    assert (
        accumulator.deposit(0, 1800, 5200, observed_introns=[(2000, 5000)], sj_strand=Strand.POS)
        is DepositOutcome.DEPOSITED
    )
    tally = accumulator.tally
    length = (2000 - 1800) + (5200 - 5000)

    spliced_by_boundary = np.asarray(tally.boundary_spliced_count, np.int64).sum(axis=1)
    assert spliced_by_boundary.tolist() == [0, 0, 0, 1, 0], (
        "exactly one CONTIGUOUS crossing — at the boundary @5050 — and none at the sj's own two "
        f"positions, which would be the splice counted twice; got {spliced_by_boundary.tolist()}"
    )
    assert np.asarray(tally.sj_count, np.int64).sum(axis=1).tolist() == [1, 0], (
        "the fragment splices TA's sj, not TB's"
    )

    boundary_at_5050 = float(np.asarray(tally.boundary_spliced_mass, np.float64)[3])
    sj = float(np.asarray(tally.sj_mass, np.float64)[0].sum())
    want_sj = Fraction(200 + 25, length)
    want_boundary = Fraction(25 + 150, length)
    assert abs(sj - want_sj) <= budget(2), (
        f"sj holds {float(sj):.6f}, want {float(want_sj):.6f} — it must claim on "
        f"BOTH its donor and acceptor sides"
    )
    assert abs(boundary_at_5050 - want_boundary) <= budget(2)
    assert abs(sj + boundary_at_5050 - 1) <= budget(4)


def test_EVERY_deposited_fragment_places_exactly_ONE_fragment_of_mass_spliced_or_not():
    """The law over a mixed population: ``Σ mass + Σ contained == n``.

    Claim 1 is scoped to unspliced fragments because the boundary rule alone does not hold across a sj; this is the same
    law once the sj axis carries its share, and a library fragment count rests on it. The enumeration sweeps the second
    block's end across the short region's far boundary, so it holds spliced fragments that cross a boundary and ones
    that cross none, plus the unspliced control. The budget is claim 1's: ``n`` fragments making ``d`` deposits drift by
    at most ``d/2`` ulp.
    """
    partition = _spliced_partition()
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    n = 0
    for start in range(1900, 1990):
        for end in (9040, 9200):  # crosses no boundary / crosses the boundary at 9050
            outcome = accumulator.deposit(
                0, start, end, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
            )
            assert outcome is DepositOutcome.DEPOSITED
            n += 1
    for start in range(100, 190):  # the unspliced control, contained in region 0
        assert accumulator.deposit(0, start, start + 60) is DepositOutcome.DEPOSITED
        n += 1

    tally = accumulator.tally
    contained = int(np.asarray(tally.region_contained_count, np.int64).sum())
    deposits = int(
        np.asarray(tally.boundary_unspliced_count, np.int64).sum()
        + np.asarray(tally.boundary_spliced_count, np.int64).sum()
        + np.asarray(tally.sj_count, np.int64).sum()
    )
    total = _total_deposited_mass(tally) + contained
    assert abs(total - n) <= budget(deposits), (
        f"{n} fragments deposited {float(total)} fragments of mass — {float(n - total)} unaccounted"
    )


# ---------------------------------------------------------------------------
# the schema
# ---------------------------------------------------------------------------


def test_the_mass_has_ONE_column_while_every_other_boundary_bank_has_TWO():
    """The mass has no strand axis: it converts an object-incidence total into a fragment count, a question with no
    strand in it, and a strand column would claim the two genome strands are separable for this channel."""
    tally = Accumulator(_partition()).tally
    n_boundaries = _partition().n_boundaries
    assert tally.boundary_unspliced_mass.shape == (n_boundaries,)
    assert tally.boundary_spliced_mass.shape == (n_boundaries,)
    assert tally.boundary_unspliced_count.shape == (n_boundaries, 2)
    for name in ("boundary_unspliced_mass", "boundary_spliced_mass", "sj_mass"):
        assert getattr(tally, name).dtype == np.float64, (
            f"{name} must be float64: counts are integers and fractions are float64."
        )
