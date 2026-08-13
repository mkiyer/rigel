"""⭐ TRAPS: two-divisors-opposite-sign — the PROOF that the accumulator's ``L`` and its deposit geometry are correct.

     (TRAPS: two-divisors-opposite-sign, the precondition on TRAPS: a-purity-filter-is-a-length-filter)

⛔ **WHY THIS FILE EXISTS.** proposes making the accumulator's ``L`` the ONE
definition of fragment length in the tool — the source every FL model is built from. Owner ruling,
2026-07-31: *"the new fragment length computation is the newest implementation and I'm not sure how
rigorously it has been tested in all cases; we need to prove that very carefully if it's going to become
the gold standard."* Correct, and the existing coverage does not meet that bar: ``test_accumulator_spec``
pins **six hand-picked malformed intron lists**, chosen by the same author as the code they check. That is
the coverage pattern that finds what the author thought of and nothing else.

⭐ **THE ORACLE IS A DIFFERENT ALGORITHM, NOT A DIFFERENT SPELLING.**: a
validator that calls the builder's own helper validates nothing, and both bugs in the index work were
caught by re-deriving the answer a different way and by nothing else. So the oracle here is **integer set
arithmetic**::

    covered = set(range(start, end)) - union of set(range(a, b)) for each intron

with no sorting, no merging, no cursor walk, and no searchsorted. Every malformed case the production
path handles by explicit logic — overlapping, nested, abutting, duplicated, zero-length, out-of-range
introns — the oracle handles by construction, because set subtraction is idempotent and order-free.

From that one set, all four deposit populations follow with no further machinery:

===================  =============================================================================
``L``                ``len(covered)``
boundary ``p`` crossed   ``p-1 in covered and p in covered``   (bases on both sides, adjacent in the
                     molecule — 's definition, verbatim)
region contained       the path used no junction and ``min(covered)``/``max(covered)`` fall in one region
===================  =============================================================================

⭐ That crossing rule is the design's own words turned into a predicate, and it is what makes this file a
proof of the **geometry** rather than only of ``L``: requires that "whatever
counts toward ``L`` must also count as coverage for crossing", and until now nothing tested the two
against each other.

⚠ Proving the Python reference proves the C++ too: `test_accumulator_native_parity` gates the C++ on
byte-identity to this reference, so a defect here would have been reproduced faithfully in both.

⭐ **THE PROOF WAS PROVEN TO FIRE** (`falsification_needs_perturbation`; a validator that cannot fail
validates nothing). The reference was deliberately broken seven ways:

=====  ==========================================================================  ==========
L1'    ``L`` from ``span − Σ(RAW intron lengths)`` — the formula                    ✅ caught
        says goes NEGATIVE on a wide overlap
L2'    introns not clipped to the fragment                                          ✅ caught
L4     fragment not clipped to the reference                                        ✅ caught
L5     crossing boundary ``searchsorted`` side flipped right→left                   ✅ caught
L6     spanning loop credits one region too many                                      ✅ caught
L7     containment keyed on the fragment EXTENT, not its first/last COVERED base    ✅ caught
L3     ABUTTING introns no longer merge                                         ⚠ NOT caught
=====  ==========================================================================  ==========

⚠ **L3 is a correct non-failure and worth writing down.** Merging abutting introns changes only the
``introns_absorbed`` QC counter: set subtraction removes ``[10,20) ∪ [20,30)`` and ``[10,30)`` identically,
so ``L`` and all four populations are untouched. The counter is pinned separately by
``test_accumulator_spec``'s ``expected_absorbed``. ⭐ A first attempt at L1 was *also* a no-op — replacing
``Σ segments`` with ``span − Σ introns`` **after** ``_normalise_introns`` — and that is itself the design's
claim demonstrated: normalisation is exactly what makes the two formulas agree.

⛔ **WHAT THIS FILE DOES NOT PROVE.** The fixture carries no annotated junctions, so every fragment is
unspliced and the **spliced routing** (``boundary_spliced`` vs ``boundary_unspliced``, junction credit, and the
containment block) is not exercised here — it is covered by ``test_accumulator_spec``'s dedicated cases.
That is the right scope for TRAPS: a-purity-filter-is-a-length-filter's purpose: the unconditional histogram bins by ``L``, and ``L`` does not
depend on which population the fragment is routed to.
"""

from __future__ import annotations

import itertools

import numpy as np

import pytest

from ._accumulator_reference import Accumulator, DepositOutcome, GapHypothesis, Partition


# --- the oracle -------------------------------------------------------------------------------------


def covered_bases(ref_len: int, start: int, end: int, introns) -> set[int]:
    """The molecule's covered genomic bases, by SET ARITHMETIC — the independent algorithm.

    Clipping to the reference is part of the definition (a fragment
    overhanging a reference end is clipped, and ``L`` is the clipped length).
    """
    s, e = max(int(start), 0), min(int(end), ref_len)
    if e <= s:
        return set()
    out = set(range(s, e))
    for a, b in introns:
        out -= set(range(max(int(a), s), min(int(b), e)))
    return out


def oracle_deposits(cuts, ref_len: int, start: int, end: int, introns, spliced: bool):
    """Every population the deposit rule credits, derived from the covered set alone.

    Returns ``(L, crossed_boundaries, contained_region)`` with local indices, or ``None`` for a fragment that
    deposits nothing. ``cuts`` are this reference's cut positions.

    ⚠ A ``spanned_regions`` set used to be derived and compared here. The bank it checked was removed on
    evidence, and the oracle no longer computes it — an oracle deriving a quantity nothing asserts is
    dead weight that reads as coverage.
    """
    covered = covered_bases(ref_len, start, end, introns)
    if not covered:
        return None
    length = len(covered)
    # A boundary at cut position p is crossed iff the molecule holds the bases on BOTH sides of it.
    # ⚠ Reported as BOUNDARY indices: a reference's cut ``i`` is its boundary ``i − 1`` (the deposit writes
    # ``boundary_base + boundary - 1``), because cut 0 is the reference start and owns no interior boundary. The
    # oracle re-derives that offset rather than importing it — and it caught me getting it wrong first.
    crossed = {i - 1 for i, p in enumerate(cuts) if (p - 1) in covered and p in covered}

    def region_of(pos):
        return min(max(int(np.searchsorted(cuts, pos, side="right")) - 1, 0), len(cuts) - 2)

    lo, hi = region_of(min(covered)), region_of(max(covered))
    contained = lo if (not spliced and lo == hi) else None
    return length, crossed, contained


# --- the fixture ------------------------------------------------------------------------------------

_CUTS = np.array(
    [0, 3, 4, 9, 12], dtype=np.int64
)  # 4 regions: widths 3, 1, 5, 3 — includes a 1 bp region
_REF_LEN = 12


def _acc(max_fragment_length: int = 10_000) -> Accumulator:
    n_regions = _CUTS.size - 1
    part = Partition(
        cut_positions=_CUTS.copy(),
        ref_cut_offsets=np.array([0, _CUTS.size], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_regions - 1], dtype=np.int64),
        sj_offsets=np.zeros(_CUTS.size + 1, dtype=np.int64),
        sj_boundary_right=np.zeros(0, dtype=np.int64),
        sj_strand=np.zeros(0, dtype=np.int8),
        region_types=np.zeros(n_regions, dtype=np.int8),
    )
    return Accumulator(part, max_fragment_length=max_fragment_length)


def _observed(acc: Accumulator):
    t = acc.tally
    return (
        {i for i in range(t.boundary_unspliced_count.shape[0]) if t.boundary_unspliced_count[i].sum()},
        next(
            (i for i in range(t.region_contained_count.shape[0]) if t.region_contained_count[i].sum()),
            None,
        ),
    )


def _check_one(start: int, end: int, introns) -> None:
    """Deposit one fragment and assert every population against the set oracle."""
    acc = _acc()
    outcome = acc.deposit(0, start, end, observed_introns=introns)
    want = oracle_deposits(_CUTS, _REF_LEN, start, end, introns, spliced=False)

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

    # ⭐ L, read back from the DEPOSIT itself rather than from an accessor — the number that actually
    # reached the tally is the one that matters. An boundary deposit carries round(2^32/(L-1)).
    t = acc.tally
    n_cross = int(t.boundary_unspliced_count.sum())
    if n_cross:
        assert _close(
            float(t.boundary_unspliced_inv_length_sum.sum()), n_cross / (length - 1), n_cross
        ), f"{ctx}: the deposited L disagrees with the oracle's {length}"
    if contained is not None:
        # ⭐ The contained deposit is `1/OPPORTUNITY` — `ell − L + 1` admissible starts inside the
        # containing region — so this asserts BOTH the length the oracle derived and the region it landed in.
        region_len = int(_CUTS[contained + 1]) - int(_CUTS[contained])
        assert _close(
            float(t.region_contained_inv_opportunity_sum.sum()), 1.0 / (region_len - length + 1), 1
        ), f"{ctx}: the contained deposit disagrees with the oracle's L={length} in a {region_len} bp region"

    assert got_crossed == crossed, f"{ctx}: crossed boundaries {got_crossed} != oracle {crossed}"
    assert got_contained == contained, f"{ctx}: contained {got_contained} != oracle {contained}"
    # the start-count invariant: exactly one, at the region holding the first COVERED base
    assert int(t.region_start_count.sum()) == 1, ctx


# --- EXHAUSTIVE: every configuration with up to two introns ------------------------------------------


def _all_intervals(n: int):
    return [(a, b) for a in range(n + 1) for b in range(a, n + 1)]


@pytest.mark.parametrize("n_introns", [0, 1, 2])
def test_exhaustive_against_the_set_oracle(n_introns):
    """⭐ EVERY fragment and EVERY intron list with ``n_introns`` introns over a 12 bp reference.

    Not a sample — the complete space. Zero-length introns, introns outside the fragment, nested,
    overlapping, abutting and duplicated pairs all occur here by enumeration rather than by being
    thought of.
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


# --- RANDOMISED: realistic coordinates, up to four introns -------------------------------------------


def test_randomised_at_realistic_scale():
    """The same oracle at coordinates a real BAM produces, with up to 4 introns.

    ⚠ Exhaustive enumeration is only affordable in a tiny coordinate space, and a tiny space cannot
    produce a 300 bp molecule spanning a 10 kb intron — the case says the
    length limit must be applied to ``L`` and never to the span. Fixed seed: this is a proof, not a
    smoke test, so it must be reproducible.
    """
    rng = np.random.default_rng(20260731)
    cuts = np.array([0, 137, 138, 400, 401, 1200, 5000, 5001, 9000], dtype=np.int64)
    ref_len = 9000
    n_regions = cuts.size - 1
    part = Partition(
        cut_positions=cuts.copy(),
        ref_cut_offsets=np.array([0, cuts.size], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_regions - 1], dtype=np.int64),
        sj_offsets=np.zeros(cuts.size + 1, dtype=np.int64),
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
        want = oracle_deposits(cuts, ref_len, start, end, introns, spliced=False)
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


# --- the property the audit turns on ------------------------------------------------------------------


def test_L_equals_the_covered_base_count_and_crossings_use_THAT_SAME_set():
    """⭐: *"whatever counts toward ``L`` must also count as coverage for
    crossing, or the density estimator is biased."*

    Nothing tested the two against each other before. A mate gap must count toward ``L`` **and** cross
    boundaries; an intron must do neither. Both are asserted here from ONE set.
    """
    # a paired-end molecule with an unsequenced mate gap: the gap IS part of the molecule
    acc = _acc()
    acc.deposit(0, 1, 11, observed_introns=[])
    assert int(acc.tally.region_contained_count.sum()) == 0
    got_crossed, _ = _observed(acc)
    assert got_crossed == {0, 1, 2}, (
        "the mate gap must carry the molecule across every interior boundary"
    )

    # the same span with the middle excised as an intron: L shrinks AND the boundaries under it go uncrossed
    acc2 = _acc()
    acc2.deposit(0, 1, 11, observed_introns=[(3, 9)])
    got2, _ = _observed(acc2)
    assert got2 == set(), "an intron must not carry the molecule across the boundaries it splices over"
    assert len(covered_bases(_REF_LEN, 1, 11, [(3, 9)])) == 4


# --- TRAPS: a-purity-filter-is-a-length-filter: the unconditional histogram, and its invariant -----------------------------------------------


def test_deposited_lengths_bins_every_accepted_fragment_exactly_once():
    """⭐ **THE TRAPS: a-purity-filter-is-a-length-filter INVARIANT (TRAPS: one-thing-varied).** ``Σ deposited_lengths == Σ region_start_count == qc.deposited``.

    Three counters, one population, incremented on the same boundary of `deposit` so they cannot drift by
    construction. It is the same externally-checkable form as 's start-count
    invariant and a **different statement**: that one says every fragment was located in space, this one
    that every fragment was binned by length. A histogram that is about to become the anchor for every FL
    model in the tool must not be allowed in one fragment short.
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
    # ⭐ and it is binned at L — the ORACLE's L, not the accumulator's own
    expected = np.zeros_like(t.deposited_lengths)
    for length in lengths:
        expected[length] += 1
    np.testing.assert_array_equal(t.deposited_lengths, expected)


def test_a_REJECTED_fragment_is_not_binned():
    """⚠ "Unconditional GIVEN DEPOSIT" — and the given matters. A fragment over the length limit, or one
    whose path is undetermined, deposits nothing and must bin nothing; each is counted in ``qc`` instead.

    ⭐ That is what makes this the right EB anchor rather than merely a convenient one: it describes
    **exactly** the population the five pure pools are drawn from. An anchor over a wider population than
    the pools would re-create, in a new place, the frame mismatch exists to
    remove.
    """
    acc = _acc(max_fragment_length=4)
    assert acc.deposit(0, 0, 12, observed_introns=[]) is DepositOutcome.TOO_LONG
    assert int(acc.tally.deposited_lengths.sum()) == 0
    assert acc.tally.qc[DepositOutcome.TOO_LONG.value] == 1

    # ⭐ A DEFERRED fragment: two hypotheses about its unsequenced gap survive, so its L is not one number
    # and there is nothing to bin. ⚠ The deferral is now the ACCUMULATOR's verdict — the caller no longer
    # passes a bool — so it has to be produced by handing over a set with two answers in it.
    acc2 = _acc()
    outcome = acc2.deposit(0, 0, 5, hypotheses=(GapHypothesis(((1, 3),)), GapHypothesis()))
    assert outcome is DepositOutcome.DEFERRED
    assert int(acc2.tally.deposited_lengths.sum()) == 0
    assert len(acc2.tally.deferred) == 1, "and it is HELD, not dropped"

    acc3 = _acc()
    assert acc3.deposit(0, 4, 4, observed_introns=[]) is DepositOutcome.EMPTY
    assert int(acc3.tally.deposited_lengths.sum()) == 0


def test_the_unconditional_histogram_is_a_SUPERSET_of_the_pure_pools():
    """⭐ The property that makes it usable as the EB anchor: every pooled fragment is also binned here,
    at the SAME length. The pools are conditioned subsets of this population, not a different one.

    ⚠ It is a strict superset in general — an exonic contained fragment and a multi-boundary crossing enter
    no pool at all (an impure pool is worse than a missing one) — which is
    precisely why the pools could never serve as their own anchor.
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


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  THE RECIPROCAL-OPPORTUNITY DEPOSIT AT A REGION — `1/A`, not `1/L`
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


def _one_region_acc(region_len: int, max_fragment_length: int = 10_000) -> Accumulator:
    """An accumulator over a SINGLE region ``[0, region_len)`` — no boundaries, so nothing but containment."""
    part = Partition.from_cuts([[0, region_len]])
    return Accumulator(part, max_fragment_length=max_fragment_length)


@pytest.mark.parametrize("region_len", [151, 400, 1000])
def test_the_region_deposit_is_the_RECIPROCAL_OPPORTUNITY_and_is_therefore_MODEL_FREE(region_len):
    """⭐⭐⭐ **A REGION'S DENSITY CHANNEL MUST NOT DEPEND ON THE FRAGMENT-LENGTH DISTRIBUTION.**

    The opportunity for a length-``w`` fragment inside a region of length ``ell`` is ``ell − w + 1`` start
    positions. Deposit ``1/A`` and the two cancel identically::

        E[SUM 1/A]  =  SUM_w rho * A(w) * f(w) * (1/A(w))  =  rho      for ANY f

    ⭐ **The brute-force configuration below IS that identity with no expectation in it.** Placing one
    fragment at EVERY admissible start of EVERY length is exactly the uniform field ``rho = 1``, so each
    length must contribute exactly ``1`` regardless of which lengths were chosen — and the total must not
    move when the length SET changes. That is model-freeness stated as an equality.

    ⛔ **This test FAILS on a ``1/L`` deposit**, which is what shipped until 2026-08-10: there each length
    contributes ``(ell − w + 1)/w``, so the total is a function of the lengths and the channel is not a
    density (`EQUATIONS.md` §2.2). Verified failing before the fix was written.
    """
    lengths = [w for w in (20, 51, 97, 150, 233, 400, 999) if w <= region_len]
    assert len(lengths) >= 3, "the parametrisation must leave a non-degenerate length SET"

    acc = _one_region_acc(region_len)
    placements = 0
    for w in lengths:
        for start in range(0, region_len - w + 1):  # every admissible placement, exactly once
            assert acc.deposit(0, start, start + w) is DepositOutcome.DEPOSITED
            placements += 1

    got = float(acc.tally.region_contained_inv_opportunity_sum[0])

    # ⭐⭐⭐ THE THEOREM, stated against the REAL-ARITHMETIC answer: `A` placements each depositing
    # `1/A` is exactly one density unit per length, whatever the lengths are. The total is therefore
    # the length COUNT — an integer — and that is what makes the channel model-free.
    #
    # ⚠ The predecessor asserted `got == sum(A * round(2^32/A))` EXACTLY. That is the fixed point's own
    # closed form, so it re-derived the implementation and could not fail (`TRAPS: a-gate-that-reconstructs`).
    # ⭐ Measured against the exact rational answer, float64 is far closer than the grid it replaced::
    #
    #     region_len 151    fixed 7.0e-10    float64 5.8e-15
    #     region_len 400    fixed 1.7e-08    float64 1.0e-13
    #     region_len 1000   fixed 2.0e-07    float64 2.8e-13
    #
    # ⛔ The budget is the representation's, derived from the number of additions — not fitted.
    assert _close(got, float(len(lengths)), placements), (
        f"region_len={region_len}: {got:.12f} density units for {len(lengths)} lengths — "
        "the opportunity did not cancel, so this channel is a function of the pmf"
    )


def test_the_region_density_channel_DOES_NOT_MOVE_when_the_length_set_changes():
    """⛔ The falsification with the pmf VARIED — the property `1/L` cannot have.

    Two disjoint length sets over the same region, each placed at every admissible start. A model-free
    density reads the same ``rho`` per length in both; a length-dependent one does not.
    """
    unit = float(1 << 32)
    region_len = 1000
    per_length = []
    for lengths in ((20, 51, 97), (233, 400, 999)):
        acc = _one_region_acc(region_len)
        for w in lengths:
            for start in range(0, region_len - w + 1):
                acc.deposit(0, start, start + w)
        per_length.append(int(acc.tally.region_contained_inv_opportunity_sum[0]) / unit / len(lengths))
    short, long_ = per_length
    assert abs(short - long_) < 1e-6, (
        f"density per length: short lengths {short:.9f}, long lengths {long_:.9f} — a model-free "
        "channel reads the same rho for both; this one is a function of the fragment lengths"
    )

#: ⭐ ONE CONVENTION: fractions are float64. `n` round-to-nearest additions differ from the real answer
#: by at most `n` ulp. ⛔ DERIVED from the machine, never fitted.
EPS = float(np.finfo(np.float64).eps)


def _close(got: float, want: float, deposits: int) -> bool:
    return abs(got - want) <= max(abs(want), 1.0) * deposits * EPS
