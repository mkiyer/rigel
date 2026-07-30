"""THE S3 GATE — the native accumulator against the executable specification, byte for byte.

    Spec: ``_accumulator_reference.py``   ·   Matrix: ``test_accumulator_spec.py``
    Plan: ``docs/IMPLEMENTATION_PLAN.md`` §0 step 6

``test_accumulator_spec.py`` says what the deposit rule *is*. This module says the C++ implements that
exact rule and no neighbouring one: the same fragments go into both accumulators and **every array, every
dtype and every QC counter must agree**. It is the only gate S3 has.

WHY IT IS DRIVEN THROUGH THE BINDING AND NOT THROUGH ``AccumulatorPayload``
    ``scan_payload.from_scan_result`` raises while ``n_channels != 4``, so until S4 lands every scan dies
    at payload construction. Comparing through the payload would therefore test nothing at all, and would
    do it while looking green. The comparison is against the ``rigel._bam_impl.Accumulator`` class.

WHY THE FIELD LIST IS NOT WRITTEN OUT HERE
    It is read off ``dataclasses.fields(Tally)``. Add a field to the specification and it joins this gate
    automatically; a binding that has not grown the matching property fails loudly. A hand-written list
    would let the two drift, which is the failure mode where a gate reads as coverage it does not have.

⚠ **This module must never be skipped.** The import is plain, so a missing or stale extension is a hard
error rather than a silent pass — a gate that can quietly not run is worse than no gate (see
``CARRY_FORWARD.md`` §3 trap 1).
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel._bam_impl import Accumulator as NativeAccumulator
from rigel.types import Strand

from ._accumulator_reference import (
    Accumulator as ReferenceAccumulator,
    Partition,
    Tally,
)


# ---------------------------------------------------------------------------
# one reference, chosen so every branch of the deposit is reachable
# ---------------------------------------------------------------------------
#
# cuts    0     100    200   201    400        900       1000
# nodes   [ n0 ][ n1  ][n2*][  n3  ][    n4   ][   n5   ]        (* n2 is 1 bp)
# lines         1      2    3      4          5
# types   intergenic, exon, exon, intron, exon, intergenic
#
# So: line 1 has {intergenic, exon} flanks (a splash pool), line 3 has {intron, exon} (the other), n2 is
# a 1 bp node that a fragment can span, and the annotated junction [201, 900) SWALLOWS line 4 — which is
# the case the whole redesign exists for.
#
# ⛔ THREE junctions, not one, and the count is load-bearing. With a single annotated junction no fragment
# can use two, so "credit only the leftmost junction" — the rule the design deliberately REVERSED, and
# which `CARRY_FORWARD.md` §3 trap 21 still recommends — was invisible to this gate: a perturbation
# implementing it passed 5/5. [100,200) and [201,900) are separated by the 1 bp exon n2, so one fragment
# can legitimately use both; [400,900) shares an acceptor with [201,900) but sits on the other strand, so
# the strand filter has something to discriminate that coordinates alone cannot.

CUTS = [0, 100, 200, 201, 400, 900, 1000]
TYPES = [0, 2, 2, 1, 2, 0]
JUNCTIONS = [
    (0, 201, 900, Strand.POS),
    (0, 100, 200, Strand.POS),
    (0, 400, 900, Strand.NEG),
]

MAX_LENGTH = 1000


def _pair(max_length: int = MAX_LENGTH, junctions=JUNCTIONS):
    """A reference accumulator and a native one over the same single-reference partition.

    ⚠ The native junction CSR is taken **from the reference's own ``Partition``** rather than rebuilt
    here. That is deliberate: the agreement between ``Partition.from_cuts`` and the index builder
    ``build_junction_edge_arrays`` is a *different* contract, already pinned by
    ``test_the_csr_slot_order_matches_the_reference_accumulator``. Feeding both sides one CSR isolates
    the thing this module is for — the deposit rule.
    """
    partition = Partition.from_cuts([CUTS], node_types=[TYPES], junctions=junctions)
    reference = ReferenceAccumulator(partition, max_fragment_length=max_length)
    native = NativeAccumulator(
        cuts=np.asarray(CUTS, dtype=np.int64),
        node_types=np.asarray(TYPES, dtype=np.uint8),
        max_length=max_length,
    )
    native.set_junctions(
        np.ascontiguousarray(partition.sj_offsets, dtype=np.int32),
        np.ascontiguousarray(partition.sj_acceptor_cut, dtype=np.int32),
        np.ascontiguousarray(partition.sj_strand, dtype=np.int8),
    )
    return reference, native


def _deposit_both(reference, native, label: str, **kw) -> None:
    """Deposit one fragment into both, then assert full parity while the fragment is still named.

    Comparing after **every** fragment rather than at the end is what makes a failure legible: the first
    disagreement names the case that caused it instead of a summed array that no longer says which
    deposit went wrong.
    """
    want = reference.deposit(0, **kw)
    got = native.deposit(**kw)
    assert got == want.value, f"{label}: outcome {got!r} != {want.value!r}"
    _assert_parity(reference, native, label)


def _assert_parity(reference, native, label: str) -> None:
    for field in dataclasses.fields(Tally):
        expected = getattr(reference.tally, field.name)
        actual = getattr(native, field.name, None)
        assert actual is not None, (
            f"{label}: the native binding has no {field.name!r}. Every field of the specification's "
            f"Tally is part of this gate; a binding that omits one is not comparable."
        )
        if field.name == "qc":
            assert dict(actual) == dict(expected), f"{label}: qc mismatch"
            continue
        assert actual.dtype == expected.dtype, (
            f"{label}: {field.name} dtype {actual.dtype} != {expected.dtype}. The dtype is part of "
            f"byte-identity — a value comparison alone would pass on a widened counter."
        )
        assert actual.shape == expected.shape, (
            f"{label}: {field.name} shape {actual.shape} != {expected.shape}"
        )
        assert np.array_equal(actual, expected), (
            f"{label}: {field.name} differs at "
            f"{np.argwhere(np.asarray(actual) != np.asarray(expected))[:8].tolist()}"
        )


# ---------------------------------------------------------------------------
# the named battery — one entry per branch of the deposit, and per bug it has had
# ---------------------------------------------------------------------------

#: ``(label, deposit kwargs)``. Ordered so that a fragment which changes state (the QC counters, the
#: junction bank) is followed by one that reads it, and every case names what it is FOR.
CASES: list[tuple[str, dict]] = [
    ("contained in an exonic node", dict(start=150, end=190)),
    ("contained, intergenic node (a pure gDNA pool)", dict(start=10, end=90)),
    ("contained, intronic node (the other pure gDNA pool)", dict(start=210, end=390)),
    ("one line crossed, {intergenic, exon} splash", dict(start=50, end=150)),
    ("one line crossed, {intron, exon} splash", dict(start=200, end=210)),
    ("four lines crossed -> no pool, it is a mixture", dict(start=50, end=500)),
    ("spanning the 1 bp node", dict(start=150, end=250)),
    ("ends exactly ON a line: contained, does NOT cross", dict(start=50, end=100)),
    ("starts exactly ON a line", dict(start=100, end=150)),
    ("minus column", dict(start=150, end=190, align_strand=Strand.NEG)),
    (
        "annotated junction, definite motif strand",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.POS),
    ),
    (
        "annotated junction, MISSING motif strand -> coordinates alone",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.NONE),
    ),
    (
        "motif strand DISAGREES with the annotation -> unannotated",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.NEG),
    ),
    (
        "CONTRADICTORY motif strand -> no splice trusted",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.AMBIGUOUS),
    ),
    (
        "annotated junction on the minus column",
        dict(
            start=150, end=950, introns=[(201, 900)], sj_strand=Strand.POS, align_strand=Strand.NEG
        ),
    ),
    (
        "implicit splice -> deposits, barred from the RNA pool",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.POS, sj_implicit=True),
    ),
    (
        "AMBIGUOUS PATH -> deposits on nothing, counted",
        dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.POS, path_ambiguous=True),
    ),
    ("strand NONE -> rejected, counted", dict(start=150, end=300, align_strand=Strand.NONE)),
    (
        "strand AMBIGUOUS -> rejected, counted",
        dict(start=150, end=300, align_strand=Strand.AMBIGUOUS),
    ),
    (
        "unannotated intron -> unspliced bank, nothing across the gap",
        dict(start=150, end=500, introns=[(210, 260)]),
    ),
    (
        "OVERLAPPING introns merge (the MO_3021 chr8 case)",
        dict(start=150, end=500, introns=[(210, 260), (240, 300)]),
    ),
    ("NESTED introns merge", dict(start=150, end=500, introns=[(200, 400), (250, 300)])),
    (
        "wide nested -- the naive L goes NEGATIVE",
        dict(start=150, end=500, introns=[(150, 480), (160, 470)]),
    ),
    (
        "ABUTTING introns merge -- a zero-length exon is impossible",
        dict(start=150, end=500, introns=[(200, 300), (300, 400)]),
    ),
    (
        "leading intron -- the path does not begin at `start`",
        dict(start=150, end=500, introns=[(150, 480)]),
    ),
    ("trailing intron", dict(start=150, end=500, introns=[(200, 500)])),
    ("zero-length intron", dict(start=150, end=500, introns=[(300, 300)])),
    ("intron entirely outside the fragment", dict(start=150, end=300, introns=[(400, 500)])),
    (
        "intron straddling the reference end -> CLIPPED, not dropped",
        dict(start=900, end=1000, introns=[(950, 1200)]),
    ),
    ("clipped at the reference end", dict(start=900, end=1200)),
    ("clipped at the reference start", dict(start=-50, end=150)),
    ("entirely off the reference -> empty, counted", dict(start=2000, end=3000)),
    ("reversed extent -> empty, counted", dict(start=500, end=400)),
    ("L == 1", dict(start=500, end=501)),
    (
        "L == 1 ON AN ANNOTATED JUNCTION -- a count against density 0",
        dict(start=201, end=901, introns=[(201, 900)], sj_strand=Strand.POS),
    ),
    (
        "two junctions credited, one annotated one not",
        dict(start=150, end=950, introns=[(201, 900), (110, 120)], sj_strand=Strand.POS),
    ),
    # ── EVERY annotated junction a path uses is credited, not just the leftmost ───────────────────────
    (
        "TWO annotated junctions on one path -> BOTH credited",
        dict(start=50, end=950, introns=[(100, 200), (201, 900)], sj_strand=Strand.POS),
    ),
    (
        "two annotated junctions, minus column",
        dict(
            start=50,
            end=950,
            introns=[(100, 200), (201, 900)],
            sj_strand=Strand.POS,
            align_strand=Strand.NEG,
        ),
    ),
    (
        "two annotated junctions with a MISSING motif strand",
        dict(start=50, end=950, introns=[(100, 200), (201, 900)], sj_strand=Strand.NONE),
    ),
    (
        "two annotated junctions, motif strand disagrees with BOTH",
        dict(start=50, end=950, introns=[(100, 200), (201, 900)], sj_strand=Strand.NEG),
    ),
    (
        "two annotated junctions, one implicit -> both credited, no RNA pool",
        dict(
            start=50,
            end=950,
            introns=[(100, 200), (201, 900)],
            sj_strand=Strand.POS,
            sj_implicit=True,
        ),
    ),
    # ── a strand-coincident acceptor: only the annotation's own strand may match ──────────────────────
    (
        "the NEG junction, matched on its own strand",
        dict(start=350, end=950, introns=[(400, 900)], sj_strand=Strand.NEG),
    ),
    (
        "the NEG junction, POS motif -> no match, and NOT the POS junction beside it",
        dict(start=350, end=950, introns=[(400, 900)], sj_strand=Strand.POS),
    ),
    (
        "the NEG junction on a coordinates-alone lookup",
        dict(start=350, end=950, introns=[(400, 900)], sj_strand=Strand.NONE),
    ),
    ("the whole reference", dict(start=0, end=1000)),
]


def test_every_named_case_is_byte_identical():
    """The battery, deposited into one accumulator pair in order so the state accumulates."""
    reference, native = _pair()
    _assert_parity(reference, native, "empty")
    for label, kw in CASES:
        _deposit_both(reference, native, label, **kw)


def test_the_fragment_length_limit_agrees_including_the_pool_histogram_width():
    """``max_length`` gates ``L`` *and* sizes the pool histograms, so it must agree on both."""
    reference, native = _pair(max_length=120)
    for label, kw in [
        ("under the limit", dict(start=150, end=260)),
        ("exactly at the limit", dict(start=150, end=270)),
        ("one over the limit -> TOO_LONG, counted", dict(start=150, end=271)),
        ("far over the limit", dict(start=0, end=1000)),
        (
            "long span, short L: the limit is on L, never the span",
            dict(start=150, end=950, introns=[(201, 900)], sj_strand=Strand.POS),
        ),
    ]:
        _deposit_both(reference, native, label, **kw)


def test_node_of_pos_agrees_everywhere_including_outside_the_reference():
    """``node_of_pos`` is public, so its clamp is reachable even though ``deposit`` cannot reach it.

    ⚠ Inside ``deposit`` the clamp is dead by construction — the path is clipped to
    ``[cuts.front(), cuts.back())`` first, so neither end can fall outside — and a perturbation removing the
    upper clamp passed the rest of this module for exactly that reason. But the method is bound, a caller
    may pass anything, and out of range it would index one past the last node. So it is pinned here rather
    than left to the branch that cannot exercise it.
    """
    reference, native = _pair()
    cuts = np.asarray(CUTS, dtype=np.int64)
    for position in [-1000, -1, *CUTS, *[c - 1 for c in CUTS], *[c + 1 for c in CUTS], 5000]:
        want = ReferenceAccumulator._local_node(cuts, position)
        assert native.node_of_pos(position) == want, f"node_of_pos({position})"
    assert reference is not None  # the pair is built for its side effects on the partition


def test_a_reference_with_no_junction_table_agrees():
    """``set_junctions`` is a separate call, so "never called" is a real state and must not differ."""
    partition = Partition.from_cuts([CUTS], node_types=[TYPES])
    reference = ReferenceAccumulator(partition, max_fragment_length=MAX_LENGTH)
    native = NativeAccumulator(
        cuts=np.asarray(CUTS, dtype=np.int64),
        node_types=np.asarray(TYPES, dtype=np.uint8),
        max_length=MAX_LENGTH,
    )
    for label, kw in [
        (
            "no table: an annotated intron is unannotated",
            dict(start=150, end=950, introns=[(201, 900)]),
        ),
        ("no table: a plain crossing", dict(start=50, end=500)),
    ]:
        _deposit_both(reference, native, label, **kw)


# ---------------------------------------------------------------------------
# the randomised arm — where an off-by-one that no named case reaches shows up
# ---------------------------------------------------------------------------


def test_ten_thousand_random_fragments_are_byte_identical():
    """⭐ The arm that actually finds things.

    A named battery tests the cases someone thought of. This one walks the whole coordinate space,
    including positions that are cuts, positions one base either side of a cut, empty and reversed
    extents, and introns that overlap in every configuration. The seed is fixed, so a failure is
    reproducible and a fix is verifiable.

    Parity is asserted on the accumulated tally at the end AND on each fragment's outcome as it goes, so
    a divergence is located to one fragment rather than to a summed array.
    """
    reference, native = _pair()
    rng = np.random.default_rng(0)
    interesting = np.array(CUTS + [c - 1 for c in CUTS] + [c + 1 for c in CUTS] + [-50, 1500])

    for i in range(10_000):
        if rng.random() < 0.5:
            start, end = (int(x) for x in rng.choice(interesting, 2))
        else:
            start, end = (int(x) for x in rng.integers(-20, 1020, 2))
        introns = []
        for _ in range(int(rng.integers(0, 4))):
            a = int(rng.choice(interesting))
            # ⚠ Half the ends are drawn from the interesting set too, so that a random intron can actually
            # LAND on an annotated junction. Drawing the end as `a + U(0, 400)` alone cannot reach the
            # 699 bp junction at all, which left the whole annotated-lookup branch to the named cases.
            # The pair is deliberately left unsorted, so reversed and zero-length introns occur.
            b = (
                int(rng.choice(interesting))
                if rng.random() < 0.5
                else a + int(rng.integers(0, 400))
            )
            introns.append((a, b))
        kw = dict(
            start=start,
            end=end,
            introns=introns,
            align_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE, Strand.AMBIGUOUS])),
            sj_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE, Strand.AMBIGUOUS])),
            sj_implicit=bool(rng.random() < 0.1),
            path_ambiguous=bool(rng.random() < 0.05),
        )
        want = reference.deposit(0, **kw)
        got = native.deposit(**kw)
        assert got == want.value, f"fragment {i} {kw}: outcome {got!r} != {want.value!r}"

    _assert_parity(reference, native, "10,000 random fragments")


def test_the_per_worker_merge_is_bit_identical_at_any_shard_count():
    """⭐ Newly achievable, and the reason every channel is an integer.

    Integer addition is associative, so sharding the same corpus K ways and merging must reproduce the
    single-accumulator answer EXACTLY, on any machine and at any thread count. The float channels this
    replaced differed by ~3.7e-7 per cell across worker counts, which propagated to a ~2.6 % difference
    in the calibration output — the same BAM giving different answers on different machines.
    """
    corpus = [kw for _, kw in CASES] * 3
    _, whole = _pair()
    for kw in corpus:
        whole.deposit(**kw)

    for n_shards in (2, 4, 8):
        shards = []
        for shard in range(n_shards):
            _, native = _pair()
            for kw in corpus[shard::n_shards]:
                native.deposit(**kw)
            shards.append(native)
        merged = shards[0]
        for other in shards[1:]:
            merged.merge_from(other)

        for field in dataclasses.fields(Tally):
            if field.name == "qc":
                assert dict(merged.qc) == dict(whole.qc), f"{n_shards} shards: qc"
                continue
            assert np.array_equal(getattr(merged, field.name), getattr(whole, field.name)), (
                f"{n_shards} shards: {field.name} is not bit-identical to the unsharded run"
            )
