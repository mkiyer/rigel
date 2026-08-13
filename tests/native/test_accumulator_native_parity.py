"""THE S3 GATE — the native accumulator against the executable specification, byte for byte.

    Spec: ``_accumulator_reference.py``   ·   Matrix: ``test_accumulator_spec.py``
    Arbitration: ``test_gap_hypothesis_arbitration.py``

``test_accumulator_spec.py`` says what the deposit rule *is*. This module says the C++ implements that
exact rule and no neighbouring one: the same fragments go into both accumulators and **every array, every
dtype and every QC counter must agree**. It is the only gate S3 has.

⭐ **AND THE DEFERRED QUEUE IS PART OF IT.** A fragment whose gap has more than one surviving explanation
is held WHOLE for the second pass, so the two languages must agree not only on the tally but on *which
fragments were held and with which hypotheses*. That bank is the only one whose ORDER is observable — every
other is a sum of integers, and integer addition is associative — so it is compared through the
canonical flattening the specification itself defines (:meth:`Tally.deferred_arrays`).

WHY IT IS DRIVEN THROUGH THE BINDING AND NOT THROUGH ``AccumulatorPayload``
    The payload is a whole-scan object over a multi-reference partition; this compares one deposit rule on
    one reference, so a divergence is located to one fragment rather than to a summed array. The payload's
    own carriage of the same banks is gated by ``tests/test_accumulator_payload.py``.

WHY THE FIELD LIST IS NOT WRITTEN OUT HERE
    It is read off ``dataclasses.fields(Tally)``. Add a field to the specification and it joins this gate
    automatically; a binding that has not grown the matching property fails loudly. A hand-written list
    would let the two drift, which is the failure mode where a gate reads as coverage it does not have.

⚠ **This module must never be skipped.** The import is plain, so a missing or stale extension is a hard
error rather than a silent pass — a gate that can quietly not run is worse than no gate (see

"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel._bam_impl import Accumulator as NativeAccumulator
from rigel.types import Strand

from ._accumulator_reference import (
    UNSPLICED_ONLY,
    Accumulator as ReferenceAccumulator,
    GapHypothesis,
    Partition,
    Tally,
)


# ---------------------------------------------------------------------------
# one reference, chosen so every branch of the deposit is reachable
# ---------------------------------------------------------------------------
#
# region_bounds    0     100    200   201    400        900       1000
# regions   [ n0 ][ n1  ][n2*][  n3  ][    n4   ][   n5   ]        (* n2 is 1 bp)
# boundaries         1      2    3      4          5
# types   intergenic, exon, exon, intron, exon, intergenic
#
# So: boundary 1 has {intergenic, exon} flanks (a splash pool), boundary 3 has {intron, exon} (the other), n2 is
# a 1 bp region that a fragment can span, and the annotated sj [201, 900) SWALLOWS boundary 4 — which is
# the case the whole redesign exists for.
#
# ⛔ THREE sj, not one, and the count is load-bearing. With a single annotated sj no fragment
# can use two, so "credit only the leftmost sj" — the rule the design deliberately REVERSED, and
# which still recommends — was invisible to this gate: a perturbation
# implementing it passed 5/5. [100,200) and [201,900) are separated by the 1 bp exon n2, so one fragment
# can legitimately use both; [400,900) shares an acceptor with [201,900) but sits on the other strand, so
# the strand filter has something to discriminate that coordinates alone cannot.

REGION_BOUNDS = [0, 100, 200, 201, 400, 900, 1000]
TYPES = [0, 2, 2, 1, 2, 0]

MAX_LENGTH = 1000

#: ⛔ **NOT ZERO, AND THAT IS THE POINT.** Every deferred record is stamped with the reference it came from,
#: because the second pass replays it through ``deposit`` onto that reference's region_bound axis — a wrong stamp
#: drains one chromosome's coordinates onto another's partition. The native accumulator has to be TOLD which
#: reference it is: it is described by its region_bound positions alone and has no other way to know.
#:
#: ⚠ Measured: with ``REF = 0`` a perturbation that hardcodes the stamp to ``0`` instead of reading the
#: accumulator's own id passed the **entire** suite, 1860 tests. A single-reference fixture cannot tell a
#: correct stamp from a constant one. So the fixture's reference is 3, and the leading three references
#: contribute no region_bounds at all — which is legal, and exercises the per-reference offset arithmetic that goes
#: negative when it is written as a plain subtraction.
REF = 3

#: The partition's reference list: three empty references, then the real one at index ``REF``.
_REGION_BOUNDS_PER_REF = [[] for _ in range(REF)] + [REGION_BOUNDS]
_TYPES_PER_REF = [[] for _ in range(REF)] + [TYPES]

SJ = [
    (REF, 201, 900, Strand.POS),
    (REF, 100, 200, Strand.POS),
    (REF, 400, 900, Strand.NEG),
]


#: ⭐ Tally fields that are not ndarrays, and how the SPECIFICATION says to compare each. Every one is
#: read off the reference, so there is no second definition of the comparison — ``deferred`` in particular
#: is a list of records for readability and is compared through the canonical flattening the reference
#: itself specifies. ⚠ A new non-array field that is not listed here fails the assertion in
#: :func:`_assert_parity` rather than quietly dropping out of the gate.
_NON_ARRAY_FIELDS = {
    "qc": lambda tally: dict(tally.qc),
    "gap_resolution": lambda tally: dict(tally.gap_resolution),
    "deferred": lambda tally: tally.deferred_arrays(),
}


def _pair(max_length: int = MAX_LENGTH, sj=SJ):
    """A reference accumulator and a native one over the same single-reference partition.

    ⚠ The native sj CSR is taken **from the reference's own ``Partition``** rather than rebuilt
    here. That is deliberate: the agreement between ``Partition.from_region_bounds`` and the index builder
    ``build_sj_arrays`` is a *different* contract, already pinned by
    ``test_the_csr_slot_order_matches_the_reference_accumulator``. Feeding both sides one CSR isolates
    the thing this module is for — the deposit rule.
    """
    partition = Partition.from_region_bounds(_REGION_BOUNDS_PER_REF, region_types=_TYPES_PER_REF, sj=sj)
    reference = ReferenceAccumulator(partition, max_fragment_length=max_length)
    native = NativeAccumulator(
        region_bounds=np.asarray(REGION_BOUNDS, dtype=np.int64),
        region_types=np.asarray(TYPES, dtype=np.uint8),
        max_length=max_length,
        ref=REF,
    )
    native.set_sj(
        np.ascontiguousarray(partition.sj_offsets, dtype=np.int32),
        np.ascontiguousarray(partition.sj_boundary_right, dtype=np.int32),
        np.ascontiguousarray(partition.sj_strand, dtype=np.int8),
    )
    return reference, native


def _deposit_both(reference, native, label: str, **kw) -> None:
    """Deposit one fragment into both, then assert full parity while the fragment is still named.

    Comparing after **every** fragment rather than at the end is what makes a failure legible: the first
    disagreement names the case that caused it instead of a summed array that no longer says which
    deposit went wrong.

    ⭐ The hypothesis set defaults to the specification's own ``UNSPLICED_ONLY`` rather than to a literal
    written here, so both sides receive the SAME objects and there is no second spelling of "a fragment
    with nothing to arbitrate".
    """
    kw.setdefault("hypotheses", UNSPLICED_ONLY)
    want = reference.deposit(REF, **kw)
    got = native.deposit(**kw)
    assert got == want.value, f"{label}: outcome {got!r} != {want.value!r}"
    _assert_parity(reference, native, label)


def _assert_array(actual, expected, label: str, name: str) -> None:
    assert actual.dtype == expected.dtype, (
        f"{label}: {name} dtype {actual.dtype} != {expected.dtype}. The dtype is part of "
        f"byte-identity — a value comparison alone would pass on a widened counter."
    )
    assert actual.shape == expected.shape, (
        f"{label}: {name} shape {actual.shape} != {expected.shape}"
    )
    assert np.array_equal(actual, expected), (
        f"{label}: {name} differs at "
        f"{np.argwhere(np.asarray(actual) != np.asarray(expected))[:8].tolist()}"
    )


def _assert_parity(reference, native, label: str) -> None:
    for field in dataclasses.fields(Tally):
        expected_raw = getattr(reference.tally, field.name)
        actual = getattr(native, field.name, None)
        assert actual is not None, (
            f"{label}: the native binding has no {field.name!r}. Every field of the specification's "
            f"Tally is part of this gate; a binding that omits one is not comparable."
        )
        if field.name in _NON_ARRAY_FIELDS:
            expected = _NON_ARRAY_FIELDS[field.name](reference.tally)
            got = dict(actual)
            assert set(got) == set(expected), (
                f"{label}: {field.name} keys {sorted(got)} != {sorted(expected)}"
            )
            for key, want in expected.items():
                if isinstance(want, np.ndarray):
                    _assert_array(np.asarray(got[key]), want, label, f"{field.name}[{key!r}]")
                else:
                    assert got[key] == want, (
                        f"{label}: {field.name}[{key!r}] is {got[key]}, expected {want}"
                    )
            continue
        assert isinstance(expected_raw, np.ndarray), (
            f"{label}: {field.name} is neither an ndarray nor listed in _NON_ARRAY_FIELDS, so this "
            f"gate does not know how the specification says to compare it"
        )
        _assert_array(actual, expected_raw, label, field.name)


# ---------------------------------------------------------------------------
# the named battery — one entry per branch of the deposit, and per bug it has had
# ---------------------------------------------------------------------------

#: ⭐ Hypotheses, named, so the cases below read as what they mean. ``()`` — region_bound nothing — is the GENOMIC
#: hypothesis and needs no flag: the gap is real template, so the molecule is gDNA or nascent.
#: ``supporting_t_inds`` are carried but never read by the first pass; they are what the second pass
#: weights a path by, so the deferred queue has to preserve them and this gate has to compare them.
GENOMIC = GapHypothesis()
LONG_SJ = GapHypothesis(((201, 900),), sj_strand=Strand.POS, supporting_t_inds=(7,))
SHORT_SJ = GapHypothesis(((100, 200),), sj_strand=Strand.POS, supporting_t_inds=(9, 11))
NEG_SJ = GapHypothesis(((400, 900),), sj_strand=Strand.NEG, supporting_t_inds=(13,))
INTRONIC_PATH = GapHypothesis(((201, 400),), sj_strand=Strand.POS)
BOTH_SJ = GapHypothesis(
    ((100, 200), (201, 900)), sj_strand=Strand.POS, supporting_t_inds=(7, 9)
)

#: ``(label, deposit kwargs)``. Ordered so that a fragment which changes state (the QC counters, the
#: sj bank, the deferred queue) is followed by one that reads it, and every case names what it is FOR.
CASES: list[tuple[str, dict]] = [
    ("contained in an exonic region", dict(start=150, end=190)),
    ("contained, intergenic region (a pure gDNA pool)", dict(start=10, end=90)),
    ("contained, intronic region (the other pure gDNA pool)", dict(start=210, end=390)),
    ("one boundary crossed, {intergenic, exon} splash", dict(start=50, end=150)),
    ("one boundary crossed, {intron, exon} splash", dict(start=200, end=210)),
    ("four boundaries crossed -> no pool, it is a mixture", dict(start=50, end=500)),
    ("spanning the 1 bp region", dict(start=150, end=250)),
    ("ends exactly ON a boundary: contained, does NOT cross", dict(start=50, end=100)),
    ("starts exactly ON a boundary", dict(start=100, end=150)),
    ("minus column", dict(start=150, end=190, align_strand=Strand.NEG)),
    (
        "annotated sj, definite motif strand",
        dict(start=150, end=950, observed_introns=[(201, 900)], sj_strand=Strand.POS),
    ),
    (
        "annotated sj, MISSING motif strand -> coordinates alone",
        dict(start=150, end=950, observed_introns=[(201, 900)], sj_strand=Strand.NONE),
    ),
    (
        "motif strand DISAGREES with the annotation -> unannotated",
        dict(start=150, end=950, observed_introns=[(201, 900)], sj_strand=Strand.NEG),
    ),
    (
        "CONTRADICTORY motif strand -> no splice trusted",
        dict(start=150, end=950, observed_introns=[(201, 900)], sj_strand=Strand.AMBIGUOUS),
    ),
    (
        "annotated sj on the minus column",
        dict(
            start=150,
            end=950,
            observed_introns=[(201, 900)],
            sj_strand=Strand.POS,
            align_strand=Strand.NEG,
        ),
    ),
    ("strand NONE -> rejected, counted", dict(start=150, end=300, align_strand=Strand.NONE)),
    (
        "strand AMBIGUOUS -> rejected, counted",
        dict(start=150, end=300, align_strand=Strand.AMBIGUOUS),
    ),
    (
        "unannotated intron -> unspliced bank, nothing across the gap",
        dict(start=150, end=500, observed_introns=[(210, 260)]),
    ),
    (
        "OVERLAPPING introns merge (the MO_3021 chr8 case)",
        dict(start=150, end=500, observed_introns=[(210, 260), (240, 300)]),
    ),
    ("NESTED introns merge", dict(start=150, end=500, observed_introns=[(200, 400), (250, 300)])),
    (
        "wide nested -- the naive L goes NEGATIVE",
        dict(start=150, end=500, observed_introns=[(150, 480), (160, 470)]),
    ),
    (
        "ABUTTING introns merge -- a zero-length exon is impossible",
        dict(start=150, end=500, observed_introns=[(200, 300), (300, 400)]),
    ),
    (
        "leading intron -- the path does not begin at `start`",
        dict(start=150, end=500, observed_introns=[(150, 480)]),
    ),
    ("trailing intron", dict(start=150, end=500, observed_introns=[(200, 500)])),
    ("zero-length intron", dict(start=150, end=500, observed_introns=[(300, 300)])),
    (
        "intron entirely outside the fragment",
        dict(start=150, end=300, observed_introns=[(400, 500)]),
    ),
    (
        "intron straddling the reference end -> CLIPPED, not dropped",
        dict(start=900, end=1000, observed_introns=[(950, 1200)]),
    ),
    ("clipped at the reference end", dict(start=900, end=1200)),
    ("clipped at the reference start", dict(start=-50, end=150)),
    ("entirely off the reference -> empty, counted", dict(start=2000, end=3000)),
    ("reversed extent -> empty, counted", dict(start=500, end=400)),
    ("L == 1", dict(start=500, end=501)),
    (
        "L == 1 ON AN ANNOTATED SJ -- a count against density 0",
        dict(start=201, end=901, observed_introns=[(201, 900)], sj_strand=Strand.POS),
    ),
    (
        "two sj credited, one annotated one not",
        dict(start=150, end=950, observed_introns=[(201, 900), (110, 120)], sj_strand=Strand.POS),
    ),
    # ── EVERY annotated sj a path uses is credited, not just the leftmost ───────────────────────
    (
        "TWO annotated sj on one path -> BOTH credited",
        dict(start=50, end=950, observed_introns=[(100, 200), (201, 900)], sj_strand=Strand.POS),
    ),
    (
        "two annotated sj, minus column",
        dict(
            start=50,
            end=950,
            observed_introns=[(100, 200), (201, 900)],
            sj_strand=Strand.POS,
            align_strand=Strand.NEG,
        ),
    ),
    (
        "two annotated sj with a MISSING motif strand",
        dict(start=50, end=950, observed_introns=[(100, 200), (201, 900)], sj_strand=Strand.NONE),
    ),
    (
        "two annotated sj, motif strand disagrees with BOTH",
        dict(start=50, end=950, observed_introns=[(100, 200), (201, 900)], sj_strand=Strand.NEG),
    ),
    # ── a strand-coincident acceptor: only the annotation's own strand may match ──────────────────────
    (
        "the NEG sj, matched on its own strand",
        dict(start=350, end=950, observed_introns=[(400, 900)], sj_strand=Strand.NEG),
    ),
    (
        "the NEG sj, POS motif -> no match, and NOT the POS sj beside it",
        dict(start=350, end=950, observed_introns=[(400, 900)], sj_strand=Strand.POS),
    ),
    (
        "the NEG sj on a coordinates-alone lookup",
        dict(start=350, end=950, observed_introns=[(400, 900)], sj_strand=Strand.NONE),
    ),
    ("the whole reference", dict(start=0, end=1000)),
    # ── the hypothesis set: ONE survivor deposits, TWO OR MORE are held WHOLE ─────────────────────────
    #
    # ⭐ These are what the arbitration is. Everything above carries the unspliced hypothesis alone, so
    # every one of them is the degenerate case — which is the general case, not a branch.
    (
        "ONE implied path -> deposits, and its intron is region_bound from L",
        dict(start=150, end=950, hypotheses=(LONG_SJ,)),
    ),
    (
        "an implied path's strand is used because NOTHING was sequenced",
        dict(start=150, end=950, sj_strand=Strand.NONE, hypotheses=(LONG_SJ,)),
    ),
    (
        "an OBSERVED motif beats the implied strand, even when they disagree",
        dict(start=350, end=950, sj_strand=Strand.NEG, hypotheses=(LONG_SJ,)),
    ),
    (
        "observed AND implied introns are both region_bound, and both sj credited",
        dict(start=50, end=950, observed_introns=[(100, 200)], hypotheses=(LONG_SJ,)),
    ),
    (
        "an implied intron DUPLICATING an observed one is absorbed, not region_bound twice",
        dict(start=50, end=950, observed_introns=[(201, 900)], hypotheses=(LONG_SJ,)),
    ),
    (
        "TWO implied paths -> DEFERRED, held whole, nothing tallied",
        dict(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ)),
    ),
    (
        "the SAME two paths again -> two records that tie on content but not on identity",
        dict(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ)),
    ),
    (
        "the two paths in the OTHER order -> a different record, and the sort must say so",
        dict(start=50, end=950, hypotheses=(SHORT_SJ, LONG_SJ)),
    ),
    (
        "a spliced path against the GENOMIC one -> deferred, RNA or gDNA",
        dict(start=150, end=500, hypotheses=(INTRONIC_PATH, GENOMIC)),
    ),
    (
        "THREE paths, one of them genomic -> deferred, both questions at once",
        dict(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ, GENOMIC)),
    ),
    (
        "a deferred fragment on the minus column, with a motif strand of its own",
        dict(
            start=50,
            end=950,
            sj_strand=Strand.NEG,
            align_strand=Strand.NEG,
            hypotheses=(LONG_SJ, NEG_SJ),
        ),
    ),
    (
        "a path of TWO implied introns, deposited",
        dict(start=50, end=950, hypotheses=(BOTH_SJ,)),
    ),
    (
        "a two-intron path against a one-intron path -> deferred; the PREFIX must not compare equal",
        dict(start=50, end=950, hypotheses=(BOTH_SJ, SHORT_SJ)),
    ),
    # ⛔ The order contract: a fragment can fail several ways and must count exactly ONCE. A fragment with
    # no genome strand is not recoverable by the second pass — that pass resolves which PATH — so the
    # strand rejection wins over the deferral, and the clipped-away fragment never reaches arbitration.
    (
        "STRAND-UNDEFINED beats the deferral: not held, and not counted as held",
        dict(
            start=50, end=950, align_strand=Strand.NONE, hypotheses=(LONG_SJ, SHORT_SJ)
        ),
    ),
    (
        "EMPTY beats the deferral: clipped to nothing before there is anything to arbitrate",
        dict(start=2000, end=3000, hypotheses=(LONG_SJ, SHORT_SJ)),
    ),
]


def test_every_named_case_is_byte_identical():
    """The battery, deposited into one accumulator pair in order so the state accumulates."""
    reference, native = _pair()
    _assert_parity(reference, native, "empty")
    for label, kw in CASES:
        _deposit_both(reference, native, label, **kw)


def test_the_battery_reaches_every_arbitration_OUTCOME():
    """⚠ Non-vacuity for the block above. Byte-identity over a bank nothing ever wrote is free.

    ⛔ Measured: with the hypothesis cases removed the deferred queue stays empty and the gap census stays
    all-zero, and the whole arbitration half of the deposit is compared only against zeros.
    """
    reference, native = _pair()
    for label, kw in CASES:
        _deposit_both(reference, native, label, **kw)

    census = reference.tally.gap_resolution
    for key, count in census.items():
        assert count > 0, f"the battery never produced {key}: {census}"
    assert len(reference.tally.deferred) == reference.tally.qc["deferred_undetermined_gap"] > 0
    assert dict(native.gap_resolution) == census


def test_the_fragment_length_limit_agrees_including_the_pool_histogram_width():
    """``max_length`` gates ``L`` *and* sizes the pool histograms, so it must agree on both.

    ⭐ It is also the ONE hypothesis filter, which is why the arbitration cases belong here rather than in
    the battery above: at the default limit of 1000 the reference's own span cannot exceed it, so the
    filter is unreachable and "the genomic hypothesis was ruled out by length" cannot be exercised at all.
    """
    reference, native = _pair(max_length=120)
    for label, kw in [
        ("under the limit", dict(start=150, end=260)),
        ("exactly at the limit", dict(start=150, end=270)),
        ("one over the limit -> TOO_LONG, counted", dict(start=150, end=271)),
        ("far over the limit", dict(start=0, end=1000)),
        (
            "long span, short L: the limit is on L, never the span",
            dict(start=150, end=950, observed_introns=[(201, 900)], sj_strand=Strand.POS),
        ),
        # ⭐ "if the genomic span exceeds the limit, assume it is RNA" is not a separate rule: the genomic
        # hypothesis's L IS the span, so the ordinary filter deletes it and the spliced path stands alone.
        (
            "the GENOMIC hypothesis is over the limit -> filtered, the spliced one deposits",
            dict(start=150, end=340, hypotheses=(GapHypothesis(((201, 300),)), GENOMIC)),
        ),
        (
            "EVERY hypothesis over the limit -> the survivors stand and TOO_LONG counts them",
            dict(start=150, end=400, hypotheses=(GapHypothesis(((201, 220),)),)),
        ),
        # ⛔ The filter emptying the set must keep ALL of them, not one: two hypotheses that are both too
        # long are still two answers, so the fragment is deferred rather than silently deposited or dropped.
        (
            "every hypothesis over the limit AND there are two -> deferred, not TOO_LONG",
            dict(start=0, end=1000, hypotheses=(GapHypothesis(((201, 220),)), GENOMIC)),
        ),
    ]:
        _deposit_both(reference, native, label, **kw)


def test_region_of_pos_agrees_everywhere_including_outside_the_reference():
    """``region_of_pos`` is public, so its clamp is reachable even though ``deposit`` cannot reach it.

    ⚠ Inside ``deposit`` the clamp is dead by construction — the path is clipped to
    ``[region_bounds.front(), region_bounds.back())`` first, so neither end can fall outside — and a perturbation removing the
    upper clamp passed the rest of this module for exactly that reason. But the method is bound, a caller
    may pass anything, and out of range it would index one past the last region. So it is pinned here rather
    than left to the branch that cannot exercise it.
    """
    reference, native = _pair()
    region_bounds = np.asarray(REGION_BOUNDS, dtype=np.int64)
    for position in [-1000, -1, *REGION_BOUNDS, *[c - 1 for c in REGION_BOUNDS], *[c + 1 for c in REGION_BOUNDS], 5000]:
        want = ReferenceAccumulator._local_region(region_bounds, position)
        assert native.region_of_pos(position) == want, f"region_of_pos({position})"
    assert reference is not None  # the pair is built for its side effects on the partition


def test_a_reference_with_no_sj_table_agrees():
    """``set_sj`` is a separate call, so "never called" is a real state and must not differ."""
    partition = Partition.from_region_bounds(_REGION_BOUNDS_PER_REF, region_types=_TYPES_PER_REF)
    reference = ReferenceAccumulator(partition, max_fragment_length=MAX_LENGTH)
    native = NativeAccumulator(
        region_bounds=np.asarray(REGION_BOUNDS, dtype=np.int64),
        region_types=np.asarray(TYPES, dtype=np.uint8),
        max_length=MAX_LENGTH,
        ref=REF,
    )
    for label, kw in [
        (
            "no table: an annotated intron is unannotated",
            dict(start=150, end=950, observed_introns=[(201, 900)]),
        ),
        ("no table: a plain crossing", dict(start=50, end=500)),
        (
            "no table: an implied path still region_bounds its intron",
            dict(start=150, end=950, hypotheses=(LONG_SJ,)),
        ),
        (
            "no table: two paths still defer",
            dict(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ)),
        ),
    ]:
        _deposit_both(reference, native, label, **kw)


def test_the_deferred_RECORD_carries_the_fragment_WHOLE_and_the_two_agree_on_it():
    """⭐ The bank the second pass reads, field by field, in one place.

    ⚠ Asserted against literals as well as against the reference. The generic comparison above says the
    two languages agree; it cannot say they agree on the *right* thing, and a bank that stored the
    UNCLIPPED extent, or dropped the hypotheses that the length filter removed, would satisfy it.
    """
    reference, native = _pair(max_length=1000)
    # Clipped at the reference start, so `start` must be 0 and not -50; and a third hypothesis is over no
    # limit here, so all three are retained — the record holds what was OFFERED, not what survived.
    _deposit_both(
        reference,
        native,
        "the record is the offered fragment, clipped",
        start=-50,
        end=950,
        observed_introns=[(110, 120)],
        sj_strand=Strand.POS,
        align_strand=Strand.NEG,
        hypotheses=(LONG_SJ, SHORT_SJ, GENOMIC),
    )
    got = dict(native.deferred)
    assert got["ref"].tolist() == [REF], (
        f"the record is stamped {got['ref'].tolist()}, not [{REF}]. The second pass replays it through "
        f"deposit onto THAT reference's region_bound axis, so a constant stamp drains one chromosome's coordinates "
        f"onto another's partition."
    )
    assert got["start"].tolist() == [0], "the CLIPPED extent is what the drain must replay"
    assert got["end"].tolist() == [950]
    assert got["align_strand"].tolist() == [int(Strand.NEG)]
    assert got["sj_strand"].tolist() == [int(Strand.POS)]
    assert got["observed_introns"].tolist() == [110, 120]
    assert got["observed_intron_offsets"].tolist() == [0, 1]
    assert got["hypothesis_offsets"].tolist() == [0, 3], (
        "all three were OFFERED, so all three are held"
    )
    assert got["hypothesis_intron_offsets"].tolist() == [0, 1, 2, 2]
    assert got["hypothesis_introns"].tolist() == [201, 900, 100, 200]
    assert got["hypothesis_t_offsets"].tolist() == [0, 1, 3, 3]
    assert got["hypothesis_t"].tolist() == [7, 9, 11]


def test_TWO_accumulators_STAMP_THEIR_OWN_REFERENCE():
    """⛔ The discriminating arm for the ``ref`` stamp, and it exists because nothing else caught it.

    ⚠ **Measured, not supposed.** A perturbation replacing ``ref_id_`` with a literal ``0`` in the deferred
    append passed **the entire suite** — 1860 tests — because every fixture in it was single-reference or
    happened to defer only on reference 0. A constant is indistinguishable from a correct value until two
    accumulators with different ids are compared side by side, which is what this does.

    ⭐ The scan path builds one Accumulator per reference and hands each its own index, so this is the unit
    of that contract: the same fragment offered to two accumulators must come back stamped differently.
    """
    held = {}
    for ref in (0, REF):
        native = NativeAccumulator(
            region_bounds=np.asarray(REGION_BOUNDS, dtype=np.int64),
            region_types=np.asarray(TYPES, dtype=np.uint8),
            max_length=MAX_LENGTH,
            ref=ref,
        )
        outcome = native.deposit(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ))
        assert outcome == "deferred_undetermined_gap", outcome
        held[ref] = dict(native.deferred)["ref"].tolist()
    assert held == {0: [0], REF: [REF]}, (
        f"the two accumulators stamped {held}; each must record the reference it IS, not a constant"
    )


# ---------------------------------------------------------------------------
# the randomised arm — where an off-by-one that no named case reaches shows up
# ---------------------------------------------------------------------------


def _random_hypotheses(rng, interesting) -> tuple[GapHypothesis, ...]:
    """0–3 hypotheses over the interesting coordinate set, plus a biased shot at the genomic one.

    ⭐ Duplicate paths, prefix paths, empty paths and shared supporting lists all occur, which is what
    exercises the deferred queue's canonical SORT — the one place in the tally where order is observable and
    therefore the one place where a comparator can be subtly wrong.
    """
    out = []
    for _ in range(int(rng.integers(0, 4))):
        if rng.random() < 0.25:
            out.append(GapHypothesis())
            continue
        introns = []
        for _ in range(int(rng.integers(1, 3))):
            a = int(rng.choice(interesting))
            b = (
                int(rng.choice(interesting))
                if rng.random() < 0.5
                else a + int(rng.integers(0, 400))
            )
            introns.append((a, b))
        out.append(
            GapHypothesis(
                tuple(introns),
                sj_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE])),
                supporting_t_inds=tuple(
                    int(t) for t in rng.integers(0, 5, int(rng.integers(0, 3)))
                ),
            )
        )
    return tuple(out) or UNSPLICED_ONLY


def test_ten_thousand_random_fragments_are_byte_identical():
    """⭐ The arm that actually finds things.

    A named battery tests the cases someone thought of. This one walks the whole coordinate space,
    including positions that are region_bounds, positions one base either side of a region_bound, empty and reversed
    extents, and introns that overlap in every configuration. The seed is fixed, so a failure is
    reproducible and a fix is verifiable.

    Parity is asserted on the accumulated tally at the end AND on each fragment's outcome as it goes, so
    a divergence is located to one fragment rather than to a summed array.
    """
    reference, native = _pair()
    rng = np.random.default_rng(0)
    interesting = np.array(REGION_BOUNDS + [c - 1 for c in REGION_BOUNDS] + [c + 1 for c in REGION_BOUNDS] + [-50, 1500])

    for i in range(10_000):
        if rng.random() < 0.5:
            start, end = (int(x) for x in rng.choice(interesting, 2))
        else:
            start, end = (int(x) for x in rng.integers(-20, 1020, 2))
        observed = []
        for _ in range(int(rng.integers(0, 4))):
            a = int(rng.choice(interesting))
            # ⚠ Half the ends are drawn from the interesting set too, so that a random intron can actually
            # LAND on an annotated sj. Drawing the end as `a + U(0, 400)` alone cannot reach the
            # 699 bp sj at all, which left the whole annotated-lookup branch to the named cases.
            # The pair is deliberately left unsorted, so reversed and zero-length introns occur.
            b = (
                int(rng.choice(interesting))
                if rng.random() < 0.5
                else a + int(rng.integers(0, 400))
            )
            observed.append((a, b))
        kw = dict(
            start=start,
            end=end,
            observed_introns=observed,
            align_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE, Strand.AMBIGUOUS])),
            sj_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE, Strand.AMBIGUOUS])),
            hypotheses=_random_hypotheses(rng, interesting),
        )
        want = reference.deposit(REF, **kw)
        got = native.deposit(**kw)
        assert got == want.value, f"fragment {i} {kw}: outcome {got!r} != {want.value!r}"

    assert reference.tally.qc["deferred_undetermined_gap"] > 100, (
        "the random arm barely deferred anything, so it is not exercising the arbitration"
    )
    _assert_parity(reference, native, "10,000 random fragments")


def test_the_per_worker_merge_is_bit_identical_at_any_shard_count():
    """⭐ Newly achievable, and the reason every channel is an integer.

    Integer addition is associative, so sharding the same corpus K ways and merging must reproduce the
    single-accumulator answer EXACTLY, on any machine and at any thread count. The float channels this
    replaced differed by ~3.7e-7 per cell across worker counts, which propagated to a ~2.6 % difference
    in the calibration output — the same BAM giving different answers on different machines.

    ⛔ **The deferred queue is the one bank this is NOT free for.** It is a list, so concatenating per-worker
    queues gives a different byte sequence at 1, 2, 4 and 8 workers with identical contents. The export
    sorts on the record's own content, and that is what this asserts.
    """
    corpus = [{"hypotheses": UNSPLICED_ONLY, **kw} for _, kw in CASES] * 3
    _, whole = _pair()
    for kw in corpus:
        whole.deposit(**kw)
    assert dict(whole.deferred)["start"].size > 0, "no fragment was deferred; the sort is untested"

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
            if field.name in _NON_ARRAY_FIELDS:
                got, want = dict(getattr(merged, field.name)), dict(getattr(whole, field.name))
                assert set(got) == set(want), f"{n_shards} shards: {field.name} keys"
                for key, expected in want.items():
                    if isinstance(expected, np.ndarray):
                        assert np.array_equal(np.asarray(got[key]), expected), (
                            f"{n_shards} shards: {field.name}[{key!r}] is not bit-identical to the "
                            f"unsharded run"
                        )
                    else:
                        assert got[key] == expected, f"{n_shards} shards: {field.name}[{key!r}]"
                continue
            got, want = getattr(merged, field.name), getattr(whole, field.name)
            # ⭐ Integer banks bit-identical (integer addition is associative); float64 fractions to
            # within the representation, because the merge re-associates their sums. See
            # `test_accumulator_worker_determinism.py`.
            if getattr(want, "dtype", None) == np.float64:
                assert np.allclose(got, want, rtol=want.size * float(np.finfo(np.float64).eps),
                                   atol=0.0), (
                    f"{n_shards} shards: {field.name} differs by MORE than the float64 representation"
                )
            else:
                assert np.array_equal(got, want), (
                    f"{n_shards} shards: {field.name} is not bit-identical to the unsharded run"
                )


def test_an_EMPTY_hypothesis_set_is_the_UNSPLICED_ONLY_set_and_does_not_CRASH():
    """⛔⛔ **A REGRESSION GATE FOR A HARD SEGFAULT** — found 2026-08-10, latent since the arbiter landed.

    The executable specification defaults ``hypotheses`` to :data:`UNSPLICED_ONLY` and says why: *"the
    degenerate case is the general case, not a branch"*. The native binding has **no default**, so a
    caller writing the natural ``hypotheses=()`` offered an EMPTY set — and an empty set walked straight
    past the ``survivors.size() > 1`` deferral into ``survivors.front()`` on an EMPTY vector, indexing a
    ``nullptr`` ``hypotheses``. ``EXC_BAD_ACCESS address=0x0``.

    ⭐ **Nothing could reach it from production**, because the scanner always offers the genomic path —
    which is exactly why it survived: it is unreachable from every code path the suite exercises, and it
    crashes the moment anyone constructs an accumulator directly. It was found while trying to build a
    deposit-behaviour digest, and it cost a session.

    ⛔ The assertion is not "does not crash" — a crash fails the test anyway. It is that an empty set
    means the SAME THING as the specification's default, so both must deposit identically.
    """
    reference, native = _pair()
    kw = dict(start=10, end=90, observed_introns=(), align_strand=int(Strand.POS))

    want = reference.deposit(REF, **kw)  # the reference's own default IS UNSPLICED_ONLY
    got = native.deposit(**kw, hypotheses=())
    assert got == want.value, f"empty hypotheses: native {got!r} != reference {want.value!r}"

    reference_default, native_default = _pair()
    reference_default.deposit(REF, **kw, hypotheses=UNSPLICED_ONLY)
    native_default.deposit(**kw, hypotheses=())
    _assert_parity(reference_default, native_default, "empty-hypotheses vs UNSPLICED_ONLY")
