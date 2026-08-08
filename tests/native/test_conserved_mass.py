"""⭐⭐⭐ THE CONSERVED MASS — the deposit rule's own gates, on the executable specification.

    Spec: ``_accumulator_reference.py``   ·   C++ parity: ``test_accumulator_native_parity.py``

**Why this bank exists.** ``edge_unspliced_count`` is ``+1`` on every line a fragment crosses, so a
fragment books ``max(K, 1)`` of them and a sum over objects is an object-INCIDENCE count, not a fragment
count. The prior the EM reads is a fragment count, and manufacturing one from a density over-calls by
+15.1 % under capture. This bank carries the count directly: **one fragment deposits exactly one.**

⭐⭐ **FOUR CLAIMS, AND THEY ARE NOT THE SAME CLAIM.** Conservation alone is far too weak to pin the
rule — the ``1/K`` rule an earlier design proposed conserves exactly too. So:

===  ==========================================================================================
 1   ``mass_quantum`` rounds halves AWAY FROM ZERO — the byte-identity contract with the C++
 2   an unspliced fragment's deposits sum to EXACTLY ONE, in rational arithmetic
 3   the mass equals the PER-BASE attribution — which ``1/K`` cannot reproduce, because a base
     does not know how many lines its fragment crossed. **This is the claim that pins the rule.**
 4   a spliced fragment's mass reaches the SPLICED bank and nothing else
===  ==========================================================================================

⚠ Claim 3 states the rule a second time, deliberately and in one place only. The specification reaches
the answer through slices and two ``searchsorted`` calls; this reaches it by asking of every base "which
node holds you, and which of that node's bounding cuts lie strictly inside this fragment?". Agreement is
not automatic — the same shape ``reference_on_real_data.py``'s ``bisect`` walk has against the
specification's index ranges (``TRAPS: a-test-that-redefines``).
"""

from __future__ import annotations

import bisect
from fractions import Fraction

import numpy as np
import pytest

from rigel.types import Strand

from ._accumulator_reference import (
    INV_LENGTH_SCALE,
    Accumulator,
    DepositOutcome,
    Partition,
    inv_length_quantum,
    mass_quantum,
)


#: ⭐ The fixture, and every part of it is load-bearing. Lines 1 and 5 have BOTH flanks wider than every
#: fragment below — the control, where a crossing fragment can only cross that one line. Lines 2, 3 and 4
#: bound a cluster of 50 bp nodes, so one 300 bp fragment crosses three at once. Under capture that
#: cluster IS the geometry: 40.2 % of gDNA fragments are crossings there against 4.8 % off it.
CUTS = [0, 400, 1000, 1050, 1100, 1500, 1900]
LENGTHS = (60, 150, 300)
MAX_LENGTH = 1000


def _partition(cuts=CUTS) -> Partition:
    return Partition.from_cuts([cuts], node_types=[[0] * (len(cuts) - 1)])


def _mass_by_base(cuts, start: int, end: int) -> dict[int, Fraction]:
    """The per-line mass of ONE unspliced fragment, attributed BASE BY BASE.

    Every base is asked which node holds it; the node's two bounding cuts are lines iff they lie
    strictly inside the fragment; the base is split equally between the lines that do bound it. Divided
    by the fragment's length, that is the mass — reached without ever forming a slice.

    ⭐ **``1/K`` is not expressible this way**, which is exactly why this is the gate that pins the rule
    rather than merely re-checking conservation.
    """
    length = end - start
    out: dict[int, Fraction] = {}
    for position in range(start, end):
        node = bisect.bisect_right(cuts, position) - 1
        lines = [i for i in (node, node + 1) if start < cuts[i] < end]
        for line in lines:
            out[line] = out.get(line, Fraction(0)) + Fraction(1, length * len(lines))
    return out


def _enumerate(partition, cuts=CUTS, lengths=LENGTHS, perturb=None):
    """Every fragment of every length at every start position. Returns ``(accumulator, n)``."""
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    n = 0
    for length in lengths:
        for start in range(cuts[0], cuts[-1] - length + 1):
            outcome = accumulator.deposit(0, start, start + length)
            assert outcome is DepositOutcome.DEPOSITED, f"[{start},{start + length}) was {outcome}"
            n += 1
    return accumulator, n


# ---------------------------------------------------------------------------
# claim 1 — the fixed-point contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "slice_len, length, n_cross, expected",
    [
        (1, 2, 1, INV_LENGTH_SCALE // 2),  # exactly a half — representable, no rounding
        (1, 1, 1, INV_LENGTH_SCALE),  # a whole fragment
        (1, 2, 2, INV_LENGTH_SCALE // 4),  # a quarter
        (1, 3, 1, (2 * INV_LENGTH_SCALE + 3) // 6),  # a third — rounds UP, away from zero
        (2, 3, 1, (4 * INV_LENGTH_SCALE + 3) // 6),  # two thirds — also away from zero
    ],
)
def test_the_mass_quantum_rounds_halves_AWAY_FROM_ZERO(slice_len, length, n_cross, expected):
    """⛔ The rounding mode is part of the contract, not an implementation detail.

    Byte-identity between the C++ and this specification is undefined without it, and Python's built-in
    ``round`` is banker's rounding, which differs at ties. Stated the same way
    :func:`inv_length_quantum`'s is, because they must agree.
    """
    assert mass_quantum(slice_len, length, n_cross) == expected


def test_the_mass_quantum_agrees_with_the_inv_length_quantum_where_they_overlap():
    """⭐ At ``slice_len == 1`` and ``n_cross == 1`` the two quanta are the SAME number by definition.

    A drift between the two rounding implementations would otherwise be invisible: each is only ever
    compared against itself.
    """
    for placements in range(1, 400):
        assert mass_quantum(1, placements, 1) == inv_length_quantum(placements)


def test_a_zero_or_negative_divisor_RAISES_rather_than_dividing():
    for bad in ((1, 0, 1), (1, 5, 0), (1, -3, 1)):
        with pytest.raises(ValueError):
            mass_quantum(*bad)


# ---------------------------------------------------------------------------
# claim 2 — conservation, exactly
# ---------------------------------------------------------------------------


def _exact_fragment_mass(partition, cuts, start, end) -> tuple[Fraction, int]:
    """``(Σ mass, contained)`` for ONE fragment, in EXACT rational arithmetic.

    ⚠ Read out of a FRESH accumulator by differencing nothing — the bank starts at zero, so the sum over
    lines after one deposit IS that fragment's deposit. The fixed point is converted back to a rational
    exactly (the quantum is an integer over a power of two), so this loses nothing.
    """
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(0, start, end)
    assert outcome is DepositOutcome.DEPOSITED
    total = Fraction(int(np.asarray(accumulator.tally.edge_unspliced_mass, np.uint64).sum()),
                     INV_LENGTH_SCALE)
    contained = int(np.asarray(accumulator.tally.node_contained_count, np.int64).sum())
    return total, contained


def test_the_RULE_conserves_exactly_before_any_quantisation():
    """⭐⭐⭐ **THE LAW ITSELF, in exact rational arithmetic.** ``Σ share == 1`` per crossing fragment.

    Every unspliced path is ONE segment. If it stays inside a node it is ``contained``. Otherwise it has
    ``n >= 2`` slices: the first crosses right only, the last crosses left only, every interior one
    crosses both. **Every slice has at least one bounding line**, so ``Σ slice_len / L = 1`` exactly and
    nothing is dropped by the "no line here" guard.

    ⛔ Stated on the RULE, through the per-base attribution, and NOT on the bank — the bank is a fixed
    point and cannot represent a third exactly. Those are two different claims and the next test makes
    the other one. Conflating them is how a quantisation residue gets mistaken for a modelling error,
    or worse, how a real error gets absorbed into a tolerance that was sized for rounding.
    """
    for length in LENGTHS:
        for start in range(CUTS[0], CUTS[-1] - length + 1):
            by_base = _mass_by_base(CUTS, start, start + length)
            if not by_base:
                continue  # contained: the node bank holds the whole fragment
            assert sum(by_base.values()) == 1, (
                f"fragment [{start},{start + length}) shares sum to {sum(by_base.values())}"
            )


def test_the_BANK_conserves_to_within_the_QUANTISATION_and_no_more():
    """⭐⭐ The same law as stored, and the price of storing it as an integer.

    ⛔ The budget is derived, not fitted (``TRAPS: no-magic-numbers``). A fragment crossing ``K`` lines
    makes ``2K`` deposits — ``K + 1`` slices, of which the two ends deposit once and each interior one
    twice — and each is rounded onto the 2^-32 grid, so the bank can differ from 1 by at most
    ``2K * (1/2) = K`` quanta.

    ⭐ The measured worst is reported as an assertion of its own, because a budget that is never
    approached says nothing about whether it is the right budget.
    """
    partition = _partition()
    worst = Fraction(0)
    for length in LENGTHS:
        for start in range(CUTS[0], CUTS[-1] - length + 1):
            total, contained = _exact_fragment_mass(partition, CUTS, start, start + length)
            crossings = len(_mass_by_base(CUTS, start, start + length))
            budget = Fraction(crossings, INV_LENGTH_SCALE)
            deviation = abs(total + contained - 1)
            worst = max(worst, deviation)
            assert deviation <= budget, (
                f"fragment [{start},{start + length}) deposited {float(total)} + {contained} "
                f"contained, off by {float(deviation):.3e} against a budget of {float(budget):.3e}"
            )
    # ⭐ Two quanta out of 2^32 at worst — 4.7e-10 of a fragment. That is the whole cost of
    # TRAPS: integer-channels-reproduce here, and it is four orders below one fragment at any depth.
    assert 0 < worst <= Fraction(2, INV_LENGTH_SCALE), (
        f"the worst per-fragment quantisation drift is {float(worst):.3e}; if it is now zero the "
        f"representation changed, and if it is larger the grid is no longer fine enough"
    )


def test_a_CONTAINED_fragment_deposits_NO_mass_at_all():
    """⛔ Its whole fragment is already ``node_contained_count``; mass at a line as well would be it
    counted twice. This is the branch the law's ``contained`` term exists for."""
    partition = _partition()
    total, contained = _exact_fragment_mass(partition, CUTS, 100, 160)  # well inside node 0
    assert contained == 1
    assert total == 0


def test_the_two_populations_together_are_EXHAUSTIVE_for_an_unspliced_fragment():
    """⭐ A contiguous interval either stays inside one node or crosses at least one line — so over the
    enumeration, ``Σ contained + Σ mass`` is the fragment COUNT, with no third case to account for.

    ⚠ At the bank's precision: the budget is the total number of crossings, one quantum each, by the
    same derivation as the per-fragment test above.
    """
    accumulator, n = _enumerate(_partition())
    mass = Fraction(int(np.asarray(accumulator.tally.edge_unspliced_mass, np.uint64).sum()),
                    INV_LENGTH_SCALE)
    contained = int(np.asarray(accumulator.tally.node_contained_count, np.int64).sum())
    crossings = int(np.asarray(accumulator.tally.edge_unspliced_count, np.int64).sum())
    assert abs(mass + contained - n) <= Fraction(crossings, INV_LENGTH_SCALE), (
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
    splits each base between the lines bounding its own node, and no per-base rule can produce ``1/K``
    because a base carries no knowledge of how many lines the whole fragment crossed.

    ⭐ **Measured, not argued** (2026-08-08). Injecting an exact ``1/K`` deposit into the specification
    leaves every other gate in this module GREEN — both conservation gates, the exhaustiveness sum, and
    even the wide-flank closed form, where ``K == 1`` makes the two rules identical — and fails **only
    this one**. It is the sole geometric gate separating the shipped rule from the design draft it
    replaced.

    The tolerance is the quantisation and nothing else — each of a line's ``<= 2 * count`` deposits is
    rounded onto the 2^-32 grid, so the bank can differ by at most ``count / 2^32`` fragments. ⛔ Not a
    fitted tolerance (``TRAPS: no-magic-numbers``).
    """
    partition = _partition()
    accumulator, _n = _enumerate(partition)
    expected: dict[int, Fraction] = {}
    for length in LENGTHS:
        for start in range(CUTS[0], CUTS[-1] - length + 1):
            for line, mass in _mass_by_base(CUTS, start, start + length).items():
                expected[line] = expected.get(line, Fraction(0)) + mass

    bank = np.asarray(accumulator.tally.edge_unspliced_mass, np.uint64)
    count = np.asarray(accumulator.tally.edge_unspliced_count, np.int64).sum(axis=1)
    for line in range(1, len(CUTS) - 1):
        edge = line - 1
        got = Fraction(int(bank[edge]), INV_LENGTH_SCALE)
        want = expected.get(line, Fraction(0))
        budget = Fraction(int(count[edge]), INV_LENGTH_SCALE)
        assert abs(got - want) <= budget, (
            f"line {line} @ {CUTS[line]}: bank {float(got)} vs per-base {float(want)}, "
            f"quantisation budget {float(budget):.3e}"
        )


def test_the_mass_EQUALS_the_count_where_both_flanks_exceed_every_fragment():
    """⭐⭐ The closed form, and the reason the shipped assembler is right OFF capture and wrong ON it.

    Where both flanking nodes are longer than every fragment, a crossing fragment can cross only that
    one line, so its whole 1.0 lands there and ``mass == count`` exactly. Where the nodes are shorter,
    one fragment is spread over several lines and ``mass`` falls far below ``count`` — that gap is the
    K-inflation, and it is why a per-object sum is not a fragment count.
    """
    partition = _partition()
    accumulator, _n = _enumerate(partition)
    bank = np.asarray(accumulator.tally.edge_unspliced_mass, np.uint64)
    count = np.asarray(accumulator.tally.edge_unspliced_count, np.int64).sum(axis=1)

    wide = [line for line in range(1, len(CUTS) - 1)
            if min(CUTS[line] - CUTS[line - 1], CUTS[line + 1] - CUTS[line]) > max(LENGTHS)]
    narrow = [line for line in range(1, len(CUTS) - 1) if line not in wide]
    assert wide and narrow, "the fixture must contain both cases or this gate tests one of them"

    for line in wide:
        assert Fraction(int(bank[line - 1]), INV_LENGTH_SCALE) == int(count[line - 1]), (
            f"line {line} has flanks wider than every fragment, so mass must equal count exactly"
        )
    for line in narrow:
        assert Fraction(int(bank[line - 1]), INV_LENGTH_SCALE) < int(count[line - 1]), (
            f"line {line} sits in the short cluster, so its count must over-state the fragments"
        )


# ---------------------------------------------------------------------------
# claim 4 — routing
# ---------------------------------------------------------------------------


def test_a_SPLICED_fragment_s_mass_reaches_the_SPLICED_bank_ALONE():
    """⛔ A spliced fragment is certified RNA at that line — gDNA cannot splice — so its mass must never
    enter the unspliced competition the prior arbitrates.

    ⭐ The geometry is the staggered case the second bank exists for: the fragment splices over one
    exon boundary and then runs CONTIGUOUSLY across another isoform's boundary inside its second block.
    """
    cuts = [0, 1000, 2000, 9000, 9050, 10000]
    partition = Partition.from_cuts(
        [cuts],
        node_types=[[0, 2, 1, 2, 2]],
        junctions=[(0, 2000, 9000, Strand.POS)],
    )
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    outcome = accumulator.deposit(
        0, 1900, 9200, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
    )
    assert outcome is DepositOutcome.DEPOSITED

    tally = accumulator.tally
    assert int(np.asarray(tally.sj_count, np.int64).sum()) == 1, "the junction must be credited"
    assert int(np.asarray(tally.edge_unspliced_mass, np.uint64).sum()) == 0, (
        "certified RNA leaked into the bank being deconvolved"
    )
    spliced = Fraction(int(np.asarray(tally.edge_spliced_mass, np.uint64).sum()), INV_LENGTH_SCALE)
    assert spliced > 0, "the contiguous crossing at 9050 deposited nothing"
    # ⭐ A PARTIAL, and the exact one the design specifies: only the second block contains an interior
    # line, so the mass is that block's share of the fragment's bases and never the whole fragment.
    # ⚠ To within one quantum — the bank is a fixed point and 2/3 is not representable on it.
    assert spliced < 1
    want = Fraction(9200 - 9000, (2000 - 1900) + (9200 - 9000))
    assert abs(spliced - want) <= Fraction(1, INV_LENGTH_SCALE)


def test_the_spliced_mass_is_a_PARTIAL_and_never_a_conservation_ledger():
    """⚠ Read as "the number of spliced fragments here", this bank is simply wrong — and that is a
    reading a consumer can make, so it is pinned rather than left to the docstring."""
    cuts = [0, 1000, 2000, 9000, 9050, 10000]
    partition = Partition.from_cuts(
        [cuts], node_types=[[0, 2, 1, 2, 2]], junctions=[(0, 2000, 9000, Strand.POS)]
    )
    accumulator = Accumulator(partition, max_fragment_length=MAX_LENGTH)
    n = 0
    for start in range(1900, 1990):
        outcome = accumulator.deposit(
            0, start, 9200, observed_introns=[(2000, 9000)], sj_strand=Strand.POS
        )
        assert outcome is DepositOutcome.DEPOSITED
        n += 1
    total = Fraction(
        int(np.asarray(accumulator.tally.edge_spliced_mass, np.uint64).sum()), INV_LENGTH_SCALE
    )
    assert 0 < total < n


# ---------------------------------------------------------------------------
# the schema
# ---------------------------------------------------------------------------


def test_the_mass_has_ONE_column_while_every_other_edge_bank_has_TWO():
    """⛔ A strand axis nothing reads is half the bank wasted, and worse, it is a claim: it says the
    two genome strands are separable for this channel. They are not — the mass exists to convert an
    object-incidence total into a fragment count, and that question has no strand in it."""
    tally = Accumulator(_partition()).tally
    n_edges = _partition().n_edges
    assert tally.edge_unspliced_mass.shape == (n_edges,)
    assert tally.edge_spliced_mass.shape == (n_edges,)
    assert tally.edge_unspliced_count.shape == (n_edges, 2)
    for name in ("edge_unspliced_mass", "edge_spliced_mass"):
        assert getattr(tally, name).dtype == np.uint64, (
            f"{name} must be uint64: it is a fixed point at 2^32 and a narrower type wraps at "
            f"realistic depth"
        )
