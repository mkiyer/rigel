"""⭐⭐⭐ THE CONSERVED MASS — the deposit rule's own gates, on the executable specification.

    Spec: ``_accumulator_reference.py``   ·   C++ parity: ``test_accumulator_native_parity.py``

**Why this bank exists.** ``boundary_unspliced_count`` is ``+1`` on every boundary a fragment crosses, so a
fragment books ``max(K, 1)`` of them and a sum over objects is an object-INCIDENCE count, not a fragment
count. The prior the EM reads is a fragment count, and manufacturing one from a density over-calls by
+15.1 % under capture. This bank carries the count directly: **one fragment deposits exactly one.**

⭐⭐ **FOUR CLAIMS, AND THEY ARE NOT THE SAME CLAIM.** Conservation alone is far too weak to pin the
rule — the ``1/K`` rule an earlier design proposed conserves exactly too. So:

⚠ Claim 1 was the fixed-point rounding contract and is GONE with the fixed point (owner, 2026-08-10:
one numeric convention). The banks are float64 fractions; there is no quantum and no scale.

===  ==========================================================================================
 2   a fragment's deposits sum to ONE, to within the float64 representation
 3   the mass equals the PER-BASE attribution — which ``1/K`` cannot reproduce, because a base
     does not know how many boundaries its fragment crossed. **This is the claim that pins the rule.**
 4   a spliced fragment's mass reaches the SPLICED bank and nothing else
 5   every crossed object SHARES the one fragment, sj boundaries included
===  ==========================================================================================

⚠ Claim 3 states the rule a second time, deliberately and in one place only. The specification reaches
the answer through slices and two ``searchsorted`` calls; this reaches it by asking of every base "which
region holds you, and which of that region's bounding region_bounds lie strictly inside this fragment?". Agreement is
not automatic — the same shape ``reference_on_real_data.py``'s ``bisect`` walk has against the
specification's index ranges (``TRAPS: a-test-that-redefines``).
"""

from __future__ import annotations

import bisect
from fractions import Fraction

import numpy as np

from rigel.types import Strand

from ._accumulator_reference import Accumulator, DepositOutcome, Partition

#: ⭐ ONE CONVENTION: the mass banks are float64 fractions. A fragment makes a bounded number of deposits,
#: each rounded to nearest, so its total differs from 1 by at most `deposits/2` ulp. ⛔ DERIVED, not
#: fitted — `np.finfo(np.float64).eps` is the machine's, not a tolerance anyone chose.
EPS = float(np.finfo(np.float64).eps)


def budget(deposits: int) -> float:
    """The representation error of `deposits` round-to-nearest additions summing to about 1."""
    return deposits * EPS


#: ⭐ The fixture, and every part of it is load-bearing. Boundaries 1 and 5 have BOTH flanks wider than every
#: fragment below — the control, where a crossing fragment can only cross that one boundary. Boundaries 2, 3 and 4
#: bound a cluster of 50 bp regions, so one 300 bp fragment crosses three at once. Under capture that
#: cluster IS the geometry: 40.2 % of gDNA fragments are crossings there against 4.8 % off it.
REGION_BOUNDS = [0, 400, 1000, 1050, 1100, 1500, 1900]
LENGTHS = (60, 150, 300)
MAX_LENGTH = 1000


def _partition(region_bounds=REGION_BOUNDS) -> Partition:
    return Partition.from_region_bounds(
        [region_bounds], region_types=[[0] * (len(region_bounds) - 1)]
    )


def _mass_by_base(region_bounds, start: int, end: int) -> dict[int, Fraction]:
    """The per-boundary mass of ONE unspliced fragment, attributed BASE BY BASE.

    Every base is asked which region holds it; the region's two bounding region_bounds are boundaries iff they lie
    strictly inside the fragment; the base is split equally between the boundaries that do bound it. Divided
    by the fragment's length, that is the mass — reached without ever forming a slice.

    ⭐ **``1/K`` is not expressible this way**, which is exactly why this is the gate that pins the rule
    rather than merely re-checking conservation.
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
# claim 2 — conservation, exactly
# ---------------------------------------------------------------------------


def _exact_fragment_mass(partition, region_bounds, start, end) -> tuple[Fraction, int]:
    """``(Σ mass, contained)`` for ONE fragment, in EXACT rational arithmetic.

    ⚠ Read out of a FRESH accumulator by differencing nothing — the bank starts at zero, so the sum over
    boundaries after one deposit IS that fragment's deposit. The fixed point is converted back to a rational
    exactly (the quantum is an integer over a power of two), so this loses nothing.
    """
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(0, start, end)
    assert outcome is DepositOutcome.DEPOSITED
    total = float(np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64).sum())
    contained = int(np.asarray(accumulator.tally.region_contained_count, np.int64).sum())
    return total, contained


def test_the_RULE_conserves_exactly_before_any_quantisation():
    """⭐⭐⭐ **THE LAW ITSELF, in exact rational arithmetic.** ``Σ share == 1`` per crossing fragment.

    Every unspliced path is ONE segment. If it stays inside a region it is ``contained``. Otherwise it has
    ``n >= 2`` slices: the first crosses right only, the last crosses left only, every interior one
    crosses both. **Every slice has at least one bounding boundary**, so ``Σ slice_len / L = 1`` exactly and
    nothing is dropped by the "no boundary here" guard.

    ⛔ Stated on the RULE, through the per-base attribution, and NOT on the bank — the bank is a fixed
    point and cannot represent a third exactly. Those are two different claims and the next test makes
    the other one. Conflating them is how a quantisation residue gets mistaken for a modelling error,
    or worse, how a real error gets absorbed into a tolerance that was sized for rounding.
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
    """⭐⭐ The same law as stored, and the price of storing it in float64.

    ⛔ The budget is derived, not fitted (``TRAPS: no-magic-numbers``). A fragment crossing ``K`` boundaries
    makes ``2K`` deposits — ``K + 1`` slices, of which the two ends deposit once and each interior one
    twice — each rounded to nearest, so the bank can differ from 1 by at most ``2K`` half-ulp.

    ⭐ The measured worst is reported as an assertion of its own, because a budget that is never
    approached says nothing about whether it is the right budget. ⚠ Under the fixed point this read
    **4.657e-10**; float64 brings it to a few ulp, which is the whole reason the convention changed.
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
    # ⭐ A few ulp at worst. ⛔ Asserted as a RANGE: zero would mean the representation changed again,
    # and more than the budget would mean the deposits are not round-to-nearest.
    assert 0 < worst <= budget(6), (
        f"the worst per-fragment representation drift is {worst:.3e}; if it is now zero the "
        f"representation changed, and if it is larger the deposits are not rounding as assumed"
    )


def test_a_CONTAINED_fragment_deposits_NO_mass_at_all():
    """⛔ Its whole fragment is already ``region_contained_count``; mass at a boundary as well would be it
    counted twice. This is the branch the law's ``contained`` term exists for."""
    partition = _partition()
    total, contained = _exact_fragment_mass(
        partition, REGION_BOUNDS, 100, 160
    )  # well inside region 0
    assert contained == 1
    assert total == 0


def test_the_two_populations_together_are_EXHAUSTIVE_for_an_unspliced_fragment():
    """⭐ A contiguous interval either stays inside one region or crosses at least one boundary — so over the
    enumeration, ``Σ contained + Σ mass`` is the fragment COUNT, with no third case to account for.

    ⚠ At the bank's precision: the budget is the total number of crossings, one quantum each, by the
    same derivation as the per-fragment test above.
    """
    accumulator, n = _enumerate(_partition())
    mass = float(np.asarray(accumulator.tally.boundary_unspliced_mass, np.float64).sum())
    contained = int(np.asarray(accumulator.tally.region_contained_count, np.int64).sum())
    crossings = int(np.asarray(accumulator.tally.boundary_unspliced_count, np.int64).sum())
    assert abs(mass + contained - n) <= budget(crossings), (
        f"{n} fragments deposited {float(mass)} of mass plus {contained} contained"
    )
    # ⭐ And it is a FRAGMENT COUNT to well under one fragment, which is the whole point of the bank.
    assert abs(float(mass) + contained - n) < 1e-6


# ---------------------------------------------------------------------------
# claim 3 — the VALUE, which is what separates this rule from 1/K
# ---------------------------------------------------------------------------


def test_the_mass_is_the_PER_BASE_attribution():
    """⭐⭐⭐ **THE CLAIM THAT PINS THE RULE.** Coverage-weighted, not ``1/K``.

    Both rules conserve, so no conservation gate can tell them apart. This one can: the per-base walk
    splits each base between the boundaries bounding its own region, and no per-base rule can produce ``1/K``
    because a base carries no knowledge of how many boundaries the whole fragment crossed.

    ⭐ **Measured, not argued** (2026-08-08). Injecting an exact ``1/K`` deposit into the specification
    leaves every other gate in this module GREEN — both conservation gates, the exhaustiveness sum, and
    even the wide-flank closed form, where ``K == 1`` makes the two rules identical — and fails **only
    this one**. It is the sole geometric gate separating the shipped rule from the design draft it
    replaced.

    The tolerance is the quantisation and nothing else — each of a boundary's ``<= 2 * count`` deposits is
    rounded onto the 2^-32 grid, so the bank can differ by at most ``count / 2^32`` fragments. ⛔ Not a
    fitted tolerance (``TRAPS: no-magic-numbers``).
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
            f"quantisation budget {float(budget):.3e}"
        )


def test_the_mass_EQUALS_the_count_where_both_flanks_exceed_every_fragment():
    """⭐⭐ The closed form, and the reason the shipped assembler is right OFF capture and wrong ON it.

    Where both flanking regions are longer than every fragment, a crossing fragment can cross only that
    one boundary, so its whole 1.0 lands there and ``mass == count`` exactly. Where the regions are shorter,
    one fragment is spread over several boundaries and ``mass`` falls far below ``count`` — that gap is the
    K-inflation, and it is why a per-object sum is not a fragment count.
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
# claim 4 — routing
# ---------------------------------------------------------------------------


def test_a_SPLICED_fragment_s_mass_reaches_the_SPLICED_bank_ALONE():
    """⛔ A spliced fragment is certified RNA at that boundary — gDNA cannot splice — so its mass must never
    enter the unspliced competition the prior arbitrates.

    ⭐ The geometry is the staggered case the second bank exists for: the fragment splices over one
    exon boundary and then runs CONTIGUOUSLY across another isoform's boundary inside its second block.
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
    # ⭐ A PARTIAL, and the exact one the design specifies: the boundary takes only the bases adjacent to it,
    # and the sj takes the rest. Derived from the rule, slice by slice, NOT read off the code::
    #
    #   L = (2000-1900) + (9200-9000) = 300
    #   block 0 [1900,2000)  crosses nothing; its only boundary is the sj   -> 100/300 sj
    #   block 1 [9000,9200)  crosses the boundary at 9050, giving two slices:
    #       [9000,9050)  bounded by the SJ (block start) and the BOUNDARY       ->  25/300 each
    #       [9050,9200)  bounded by the BOUNDARY alone (block end is the frag end)    -> 150/300 boundary
    #
    # ⚠ The predecessor of this assertion wanted 200/300 for the boundary — it encoded the rule in which a
    # boundary-crossing block gave its bases ENTIRELY to boundaries, so the sj bounding it got nothing.
    assert spliced < 1
    want_boundary = Fraction(25, 300) + Fraction(150, 300)
    want_sj = Fraction(100, 300) + Fraction(25, 300)
    assert abs(spliced - want_boundary) <= budget(2)
    assert abs(sj - want_sj) <= budget(2)
    # ⭐ And together they are the whole fragment, which the predecessor could not state at all.
    assert abs(spliced + sj - 1) <= budget(4)


def test_the_spliced_mass_is_a_PARTIAL_and_never_a_conservation_ledger():
    """⚠ Read as "the number of spliced fragments here", this bank is simply wrong — and that is a
    reading a consumer can make, so it is pinned rather than left to the docstring."""
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
# claim 5 — the law over ALL THREE axes, which is what makes a LIBRARY count computable
# ---------------------------------------------------------------------------

#: The staggered geometry claim 4 uses: one annotated sj, and a short region (9000–9050) whose far
#: boundary a fragment may or may not reach. That choice is what lets the SAME partition carry both a
#: spliced fragment that crosses a boundary and one that crosses none.
SPLICED_REGION_BOUNDS = [0, 1000, 2000, 9000, 9050, 10000]


def _spliced_partition() -> Partition:
    return Partition.from_region_bounds(
        [SPLICED_REGION_BOUNDS], region_types=[[0, 2, 1, 2, 2]], sj=[(0, 2000, 9000, Strand.POS)]
    )


def _total_deposited_mass(tally) -> Fraction:
    """Every conserved-mass bank summed — unspliced boundaries, spliced boundaries, and sj.

    ⚠ ``sj_mass`` is read with :func:`getattr` so that when the bank does not exist this reports the
    LAW failing by a whole fragment, rather than an ``AttributeError`` that says nothing about how much
    mass went missing. A falsification has to fail for the reason it is about.
    """
    total = Fraction(0)
    for name in ("boundary_unspliced_mass", "boundary_spliced_mass", "sj_mass"):
        bank = getattr(tally, name, None)
        if bank is not None:
            total += float(np.asarray(bank, np.float64).sum())
    return total


def test_a_spliced_fragment_crossing_NO_BOUNDARY_still_deposits_a_WHOLE_FRAGMENT_of_mass():
    """⛔⛔ **THE POPULATION A CONSERVED COUNT CANNOT SEE.** Both blocks lie inside one region each.

    The fragment is not ``contained`` — its path spans a sj, so it is not inside ONE region — and it
    crosses no boundary, so the slice loop never runs. Under the rule as it stands it deposits **nothing
    anywhere** while ``sj_count`` credits it, which is a fragment that exists on the incidence axis and
    on no conserved one.

    ⭐ Measured on the origin-split oracle at ladder g50 capture_off: **1,222,375 of 4,830,713 RNA
    fragments (25.3 %)** are in this population, against **0 of 4,997,761** gDNA fragments — gDNA cannot
    splice, so its conserved count is already complete. That asymmetry is why the library ``f_gdna``
    is computable from the gDNA side alone and the RNA side is not.
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


#: ⛔⛔ **TWO sj, and the middle exon is a WHOLE region.** A one-sj fixture cannot exercise
#: the equal-share rule at all — every block it produces is bounded on one side only, so ``len(bounds)``
#: is always 1 and a rule that forgot to share would pass every gate. Measured: deleting the share made
#: 18/18 gates pass until this fixture existed (`TRAPS: could-the-arm-have-fired`).
TWO_SJ_REGION_BOUNDS = [0, 1000, 2000, 5000, 5100, 8000, 9000, 10000]


def _two_sj_partition() -> Partition:
    return Partition.from_region_bounds(
        [TWO_SJ_REGION_BOUNDS],
        region_types=[[0, 2, 1, 2, 1, 2, 0]],
        sj=[(0, 2000, 5000, Strand.POS), (0, 5100, 8000, Strand.POS)],
    )


#: The three blocks of ``[1950, 8050)`` spliced over both introns: 50 bp, the whole 100 bp middle exon,
#: then 50 bp. None contains an interior boundary, so all three are sj-bounded and the middle one is
#: bounded TWICE.
TWO_SJ_FRAGMENT = (1950, 8050, [(2000, 5000), (5100, 8000)])


def test_a_block_bounded_by_TWO_sj_SPLITS_its_mass_between_them():
    """⭐⭐⭐ **THE CLAIM THAT PINS THE SJ RULE**, the way claim 3 pins the boundary rule.

    Conservation alone cannot see the difference between sharing a doubly-bounded block and giving it
    whole to each bound — one of those sums to 1 only because the fixture never produced such a block.
    This states WHERE the mass went, per sj, in exact rational arithmetic::

        block 0  [1950,2000)   50/200 -> j0            (fragment start on the left: one bound)
        block 1  [5000,5100)  100/200 -> j0 and j1     (a sj at BOTH ends: shared)
        block 2  [8000,8050)   50/200 -> j1            (fragment end on the right: one bound)

    so ``j0 == j1 == 1/2`` and the two sum to exactly one fragment.
    """
    start, end, introns = TWO_SJ_FRAGMENT
    accumulator = Accumulator(_two_sj_partition(), max_fragment_length=MAX_LENGTH)
    assert (
        accumulator.deposit(0, start, end, observed_introns=introns, sj_strand=Strand.POS)
        is DepositOutcome.DEPOSITED
    )
    # ⭐ SUMMED over the strand columns. `sj_mass` went per strand on 2026-08-13 for artifact
    # detection, and CONSERVATION is a strand-agnostic property — a fragment's 1.0 is shared across the
    # objects it touched whatever strand it aligned to — so every gate in this file reads the total.
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
    """⭐⭐⭐ **THE OWNER'S REQUIREMENT, STATED AS A GATE** (2026-08-10): a spliced fragment may cross
    several sj *and* several boundaries, and **all of them share the one fragment**.

    ⛔ **Conservation alone does NOT imply sharing, and the first implementation proved it.** That rule
    gave a boundary-crossing block's bases entirely to boundaries, so a sj flanked by two boundary-crossing
    blocks received exactly **0.000000** while ``sj_count`` credited it — and the fragment still summed
    to 1.0, so every conservation gate passed. Measured on real geometry: 35 of 8,436 crossed sj
    held zero mass. This gate asserts the property conservation cannot see.

    The fixture is the awkward case on purpose::

        block 0  [1900,2500)  crosses boundary @2000     — and is bounded on the right by sj 0
        block 1  [5000,5200)  crosses boundary @5100     — bounded by sj 0 left, sj 1 right
        block 2  [8000,8100)  crosses NOTHING        — bounded by sj 1 on the left

    so sj 0 is flanked by two boundary-crossing blocks: under the old rule it got nothing.
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
    """⭐⭐⭐ **THE OWNER'S CASE (2026-08-10), and it discriminates where conservation cannot.**

    Two isoforms sharing a donor and differing by 50 bp at the acceptor::

        TA+  exons (1000,2000) (5000,10000)
        TB+  exons (1000,2000) (5050,10000)

    A fragment with blocks ``(1800,2000)`` and ``(5000,5200)`` — ``L = 400`` — splices TA's
    ``(2000,5000)`` and then runs CONTIGUOUSLY across TB's acceptor boundary at 5050. It therefore touches
    THREE boundary positions, of which the sj accounts for two:

    ===============  ======  ==================================================  ==========
    bases            n       bounded by                                          share
    ===============  ======  ==================================================  ==========
    ``[1800,2000)``  200     the sj's DONOR side (fragment start on left)   200 sj
    ``[5000,5050)``   50     the sj's ACCEPTOR side, and the boundary @5050      25 / 25
    ``[5050,5200)``  150     the boundary @5050 (fragment end on right)               150 boundary
    ===============  ======  ==================================================  ==========

    so the sj holds ``225/400`` and the boundary ``175/400``.

    ⛔ **The predecessor rule scored 200/200 here** — it gave a boundary-crossing block's bases entirely to
    boundaries, so the sj never collected on its ACCEPTOR side. Both rules sum to 1.0, so no
    conservation gate could tell them apart; this one can, which is why the case is pinned.

    ⛔ **And a sj must not ALSO register as a contiguous crossing** — that would be the same
    traversal counted twice. It cannot, structurally: a boundary is crossed iff it lies STRICTLY inside a
    block, and a sj's donor and acceptor are block endpoints. Asserted rather than assumed.
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
    """⭐⭐⭐ **THE LAW, over a population that mixes both.** ``Σ mass + Σ contained == n``.

    Claim 2 states this for unspliced fragments only, and that scoping was correct rather than an
    oversight — the rule genuinely did not hold across a sj. This is the same law once the
    sj axis carries its share, and it is the identity a library fragment count rests on.

    The enumeration sweeps the second block's end across the short region's far boundary, so it contains
    spliced fragments that cross a boundary and spliced fragments that cross none, plus the unspliced
    population for the control.

    ⚠ The budget is the quantisation and nothing else, by claim 2's derivation: every deposit is
    rounded onto the 2^-32 grid, so ``n`` fragments making ``d`` deposits drift by at most ``d/2``
    quanta. ⛔ Not a fitted tolerance (``TRAPS: no-magic-numbers``).
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
    """⛔ A strand axis nothing reads is half the bank wasted, and worse, it is a claim: it says the
    two genome strands are separable for this channel. They are not — the mass exists to convert an
    object-incidence total into a fragment count, and that question has no strand in it."""
    tally = Accumulator(_partition()).tally
    n_boundaries = _partition().n_boundaries
    assert tally.boundary_unspliced_mass.shape == (n_boundaries,)
    assert tally.boundary_spliced_mass.shape == (n_boundaries,)
    assert tally.boundary_unspliced_count.shape == (n_boundaries, 2)
    for name in ("boundary_unspliced_mass", "boundary_spliced_mass", "sj_mass"):
        assert getattr(tally, name).dtype == np.float64, (
            f"{name} must be float64. ⭐ ONE CONVENTION: a COUNT is an integer, a FRACTION is float64. "
            f"A fixed point here would be a second convention in the same schema, and it was measured "
            f"LESS accurate than float64 on the reciprocal-opportunity theorem (7.0e-10 against "
            f"5.8e-15 at region_len 151)."
        )
