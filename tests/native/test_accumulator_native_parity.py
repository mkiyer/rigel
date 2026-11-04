"""The native accumulator against the executable specification, byte for byte.

``test_accumulator_spec.py`` states what the deposit rule is; this module checks that the C++ implements exactly that
rule. The same fragments go into both accumulators, and every array, dtype and QC counter must agree: first over a named
battery with one entry per branch of the deposit, then over ten thousand random fragments, with parity asserted after
each so a divergence is located to one fragment.

Every bank but one is a per-object sum (uint32 counts, float64 fractions and masses) that both languages add in the same
fragment order, so each must match exactly. The deferred queue (fragments whose gap has more than one surviving
explanation, held whole for the second pass) is a list, so it is compared through the specification's canonical
flattening, :meth:`Tally.deferred_arrays`. The field list is read off ``dataclasses.fields(Tally)``, so a new channel
joins the comparison automatically and a binding missing the matching property fails. The comparison goes through the
binding rather than ``AccumulatorPayload``, which is a whole-scan, multi-reference object, because this compares one
deposit rule on one reference. The import is plain, so a missing or stale extension is an error, not a skip.
Specification: ``_accumulator_reference.py``.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel._bam_impl.testing import Accumulator as NativeAccumulator
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
# So: boundary 1 has {intergenic, exon} flanks (a splash pool), boundary 3 has {intron, exon} (the
# other), n2 is a 1 bp region that a fragment can span, and the annotated sj [201, 900) swallows
# boundary 4.
#
# Three sj, because with a single annotated sj no fragment can use two, and a "credit only the
# leftmost sj" rule would pass. [100,200) and [201,900) are separated by the 1 bp exon n2, so one
# fragment can use both; [400,900) shares an acceptor with [201,900) but sits on the other strand, so
# the strand filter has something to discriminate that coordinates alone cannot.

REGION_BOUNDS = [0, 100, 200, 201, 400, 900, 1000]
TYPES = [0, 2, 2, 1, 2, 0]

MAX_LENGTH = 1000

#: Not zero, on purpose. Every deferred record is stamped with its source reference, because the second pass replays it
#: through ``deposit`` onto that reference's region-bound axis; a wrong stamp would drain one chromosome's coordinates
#: onto another's partition. The native accumulator knows its reference only because it is told. A single-reference
#: fixture at ``REF = 0`` cannot tell a correct stamp from a constant, so the fixture's reference is 3, and the leading
#: three references contribute no bounds, which is legal and exercises per-reference offset arithmetic that would go
#: negative as a plain subtraction.
REF = 3

#: The partition's reference list: three empty references, then the real one at index ``REF``.
_REGION_BOUNDS_PER_REF = [[] for _ in range(REF)] + [REGION_BOUNDS]
_TYPES_PER_REF = [[] for _ in range(REF)] + [TYPES]

SJ = [
    (REF, 201, 900, Strand.POS),
    (REF, 100, 200, Strand.POS),
    (REF, 400, 900, Strand.NEG),
]


#: Tally fields that are not ndarrays, and how the specification says to compare each. Each comparison is read off the
#: reference, so there is no second definition: ``deferred`` in particular is a list of records, compared through the
#: reference's canonical flattening. An unlisted non-array field fails the assertion in :func:`_assert_parity` rather
#: than dropping out of the comparison.
_NON_ARRAY_FIELDS = {
    "qc": lambda tally: dict(tally.qc),
    "gap_resolution": lambda tally: dict(tally.gap_resolution),
    "deferred": lambda tally: tally.deferred_arrays(),
}


def _pair(max_length: int = MAX_LENGTH, sj=SJ):
    """A reference accumulator and a native one over the same single-reference partition.

    The native sj CSR is taken from the reference's own ``Partition`` rather than rebuilt here. Agreement between
    ``Partition.from_region_bounds`` and the index builder ``build_sj_arrays`` is a separate contract tested elsewhere;
    feeding both sides one CSR isolates the deposit rule.
    """
    partition = Partition.from_region_bounds(
        _REGION_BOUNDS_PER_REF, region_types=_TYPES_PER_REF, sj=sj
    )
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
    """Deposit one fragment into both, then assert full parity while the fragment is still named, so the first
    disagreement names the case that caused it. The hypothesis set defaults to the specification's own
    ``UNSPLICED_ONLY``, so both sides receive the same objects.
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
# the named battery: one entry per branch of the deposit
# ---------------------------------------------------------------------------

#: Hypotheses, named, so the cases below read as what they mean. ``()`` (cut nothing) is the genomic hypothesis and
#: needs no flag: the gap is real template, so the molecule is gDNA or unspliced RNA.
GENOMIC = GapHypothesis()
LONG_SJ = GapHypothesis(((201, 900),), sj_strand=Strand.POS)
SHORT_SJ = GapHypothesis(((100, 200),), sj_strand=Strand.POS)
NEG_SJ = GapHypothesis(((400, 900),), sj_strand=Strand.NEG)
INTRONIC_PATH = GapHypothesis(((201, 400),), sj_strand=Strand.POS)
BOTH_SJ = GapHypothesis(((100, 200), (201, 900)), sj_strand=Strand.POS)

#: ``(label, deposit kwargs)``. Ordered so that a fragment which changes state (the QC counters, the sj bank, the
#: deferred queue) is followed by one that reads it, and every case says what it tests.
CASES: list[tuple[str, dict]] = [
    ("contained in an exonic region", dict(start=150, end=190)),
    ("contained, intergenic region (a gDNA pool)", dict(start=10, end=90)),
    ("contained, intronic region (the other contained gDNA pool)", dict(start=210, end=390)),
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
        "OVERLAPPING introns merge (the cfRNA-C chr8 case)",
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
    # ── every annotated sj a path uses is credited, not just the leftmost ───────────────────────
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
    # ── the hypothesis set: one survivor deposits, two or more are held whole ─────────────────────────
    #
    # Everything above carries the unspliced hypothesis alone, the degenerate case of the same rule.
    (
        "ONE implied path -> deposits, and its intron is cut from L",
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
        "observed AND implied introns are both cut, and both sj credited",
        dict(start=50, end=950, observed_introns=[(100, 200)], hypotheses=(LONG_SJ,)),
    ),
    (
        "an implied intron DUPLICATING an observed one is absorbed, not cut twice",
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
    # A fragment can fail several ways and must count exactly once. A fragment with no genome strand cannot be recovered
    # by the second pass, which resolves only the path, so the strand rejection wins over the deferral, and a
    # clipped-away fragment never reaches arbitration.
    (
        "STRAND-UNDEFINED beats the deferral: not held, and not counted as held",
        dict(start=50, end=950, align_strand=Strand.NONE, hypotheses=(LONG_SJ, SHORT_SJ)),
    ),
    (
        "EMPTY beats the deferral: clipped to nothing before there is anything to arbitrate",
        dict(start=2000, end=3000, hypotheses=(LONG_SJ, SHORT_SJ)),
    ),
]


def test_the_battery_is_byte_identical_and_reaches_every_arbitration_outcome():
    """The battery, deposited in order into one accumulator pair so state accumulates, with parity asserted after every
    case. It also checks that the deferred queue and the gap census are non-empty, since byte-identity over banks
    nothing wrote compares only zeros.
    """
    reference, native = _pair()
    _assert_parity(reference, native, "empty")
    for label, kw in CASES:
        _deposit_both(reference, native, label, **kw)

    census = reference.tally.gap_resolution
    for key, count in census.items():
        assert count > 0, f"the battery never produced {key}: {census}"
    assert len(reference.tally.deferred) == reference.tally.qc["deferred_undetermined_gap"] > 0
    assert dict(native.gap_resolution) == census


def test_the_fragment_length_limit_agrees_including_the_pool_histogram_width():
    """``max_length`` gates ``L`` and sizes the pool histograms, so it must agree on both. It is also the one hypothesis
    filter, so the arbitration cases that need it live here: at the default limit of 1000 the reference's span cannot
    exceed it, and "the genomic hypothesis was ruled out by length" would be unreachable.
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
        # "if the genomic span exceeds the limit, it is RNA" is not a separate rule: the genomic hypothesis's L is the
        # span, so the ordinary filter deletes it and the spliced path stands alone.
        (
            "the GENOMIC hypothesis is over the limit -> filtered, the spliced one deposits",
            dict(start=150, end=340, hypotheses=(GapHypothesis(((201, 300),)), GENOMIC)),
        ),
        (
            "EVERY hypothesis over the limit -> the survivors stand and TOO_LONG counts them",
            dict(start=150, end=400, hypotheses=(GapHypothesis(((201, 220),)),)),
        ),
        # When the filter would empty the set it keeps all of them: two hypotheses that are both too long are still two
        # answers, so the fragment is deferred rather than deposited or dropped.
        (
            "every hypothesis over the limit AND there are two -> deferred, not TOO_LONG",
            dict(start=0, end=1000, hypotheses=(GapHypothesis(((201, 220),)), GENOMIC)),
        ),
    ]:
        _deposit_both(reference, native, label, **kw)


def test_region_of_pos_agrees_everywhere_including_outside_the_reference():
    """``region_of_pos`` is public, so its clamp is reachable even though ``deposit`` cannot reach it: inside
    ``deposit`` the path is clipped to ``[region_bounds.front(), region_bounds.back())`` first. A caller may pass
    anything, and out of range the method would index one past the last region.
    """
    reference, native = _pair()
    region_bounds = np.asarray(REGION_BOUNDS, dtype=np.int64)
    for position in [
        -1000,
        -1,
        *REGION_BOUNDS,
        *[c - 1 for c in REGION_BOUNDS],
        *[c + 1 for c in REGION_BOUNDS],
        5000,
    ]:
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
            "no table: an implied path still cuts its intron",
            dict(start=150, end=950, hypotheses=(LONG_SJ,)),
        ),
        (
            "no table: two paths still defer",
            dict(start=50, end=950, hypotheses=(LONG_SJ, SHORT_SJ)),
        ),
    ]:
        _deposit_both(reference, native, label, **kw)


def test_the_deferred_record_carries_the_fragment_whole_and_the_two_agree_on_it():
    """The bank the second pass reads, field by field. Asserted against literals as well as the reference: parity alone
    cannot catch a bank that both sides store wrongly, such as the unclipped extent, or one that drops the hypotheses
    the length filter removed.
    """
    reference, native = _pair(max_length=1000)
    # Clipped at the reference start, so `start` must be 0 and not -50; and the third hypothesis is over no limit here,
    # so all three are retained: the record holds what was offered, not what survived.
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


def test_two_accumulators_stamp_their_own_reference():
    """Two accumulators with different ids must stamp the same fragment differently. Every other fixture is
    single-reference or defers only on reference 0, so a constant stamp would pass them. The scan builds one Accumulator
    per reference and hands each its own index; this is the unit of that contract.
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
# the randomised battery: where an off-by-one no named case reaches shows up
# ---------------------------------------------------------------------------


def _random_hypotheses(rng, interesting) -> tuple[GapHypothesis, ...]:
    """0–3 hypotheses over the interesting coordinate set, plus a biased shot at the genomic one. Duplicate, prefix and
    empty paths all occur, which exercises the deferred queue's canonical sort, the one place in the tally where order
    is observable.
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
                tuple(introns), sj_strand=int(rng.choice([Strand.POS, Strand.NEG, Strand.NONE]))
            )
        )
    return tuple(out) or UNSPLICED_ONLY


def test_ten_thousand_random_fragments_are_byte_identical():
    """The named battery tests the cases someone thought of; this one walks the coordinate space, including region
    bounds, positions one base either side of a bound, empty and reversed extents, and introns overlapping in every
    configuration. The seed is fixed so a failure is reproducible. Parity is asserted on each fragment's outcome as it
    goes and on the accumulated tally at the end.
    """
    reference, native = _pair()
    rng = np.random.default_rng(0)
    interesting = np.array(
        REGION_BOUNDS
        + [c - 1 for c in REGION_BOUNDS]
        + [c + 1 for c in REGION_BOUNDS]
        + [-50, 1500]
    )

    for i in range(10_000):
        if rng.random() < 0.5:
            start, end = (int(x) for x in rng.choice(interesting, 2))
        else:
            start, end = (int(x) for x in rng.integers(-20, 1020, 2))
        observed = []
        for _ in range(int(rng.integers(0, 4))):
            a = int(rng.choice(interesting))
            # Half the ends are drawn from the interesting set too, so a random intron can land on an annotated sj; `a +
            # U(0, 400)` alone cannot reach the 699 bp sj. The pair is left unsorted, so reversed and zero-length
            # introns occur.
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


def test_the_batched_sj_lookup_is_the_scalar_rule_written_out():
    """`sj_edge_ids` answers many introns at once, and the second pass reads it instead of keeping its own copy of the
    rule in Python.

    The oracle here is that rule written out: both endpoints must be region bounds, the CSR is walked from the donor,
    and the strand filter applies only when the motif is definite. The two must agree on every case the fixture can
    pose: the three annotated junctions at their own strand, the opposite strand and NONE; the strand-discriminating
    pair that shares an acceptor; endpoints that are not region bounds; and a reversed pair.
    """
    _, native = _pair()
    partition = Partition.from_region_bounds(
        _REGION_BOUNDS_PER_REF, region_types=_TYPES_PER_REF, sj=SJ
    )
    offsets = np.asarray(partition.sj_offsets, dtype=np.int64)
    right = np.asarray(partition.sj_boundary_right, dtype=np.int64)
    strand = np.asarray(partition.sj_strand, dtype=np.int64)
    bounds = np.asarray(REGION_BOUNDS, dtype=np.int64)

    def oracle(start: int, end: int, motif: int) -> int:
        def bound(position: int) -> int:
            k = int(np.searchsorted(bounds, position))
            return k if k < bounds.size and int(bounds[k]) == position else -1

        donor, acceptor = bound(start), bound(end)
        if donor < 0 or acceptor < 0:
            return -1
        definite = motif in (int(Strand.POS), int(Strand.NEG))
        for k in range(int(offsets[donor]), int(offsets[donor + 1])):
            if int(right[k]) != acceptor:
                continue
            if definite and int(strand[k]) != motif:
                continue
            return k
        return -1

    cases = []
    for a, b in ((201, 900), (100, 200), (400, 900), (900, 201), (0, 1000), (150, 900), (100, 201)):
        for motif in (int(Strand.POS), int(Strand.NEG), int(Strand.NONE)):
            cases.append((a, b, motif))
    starts = np.array([c[0] for c in cases], dtype=np.int64)
    ends = np.array([c[1] for c in cases], dtype=np.int64)
    motifs = np.array([c[2] for c in cases], dtype=np.int32)

    got = np.asarray(native.sj_edge_ids(starts=starts, ends=ends, sj_strand=motifs), dtype=np.int64)
    want = np.array([oracle(*c) for c in cases], dtype=np.int64)
    assert np.array_equal(got, want), dict(zip(cases, zip(got.tolist(), want.tolist())))
    # the grid must find junctions and miss some, or it proves nothing
    assert (want >= 0).any() and (want < 0).any()
    # and the empty call is legal, because a fragment set with no implied intron reaches it
    assert (
        np.asarray(
            native.sj_edge_ids(
                starts=np.zeros(0, np.int64),
                ends=np.zeros(0, np.int64),
                sj_strand=np.zeros(0, np.int32),
            )
        ).size
        == 0
    )


def test_the_per_worker_merge_is_bit_identical_at_any_shard_count():
    """Integer addition is associative, so sharding the same corpus K ways and merging must reproduce the
    single-accumulator answer exactly, on any machine and at any thread count; that is why every count channel is an
    integer. A float count would let the same BAM give different calibration output on different machines. The deferred
    queue is a list, so concatenated per-worker queues differ in order across worker counts; the export sorts on each
    record's own content, and this asserts it.
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
            # Integer banks bit-identical (integer addition is associative); float64 fractions within the
            # representation, because the merge re-associates their sums.
            if getattr(want, "dtype", None) == np.float64:
                assert np.allclose(
                    got, want, rtol=want.size * float(np.finfo(np.float64).eps), atol=0.0
                ), (
                    f"{n_shards} shards: {field.name} differs by MORE than the float64 representation"
                )
            else:
                assert np.array_equal(got, want), (
                    f"{n_shards} shards: {field.name} is not bit-identical to the unsharded run"
                )


def _exported(native) -> dict:
    """Every bank the binding exports, copied. The array properties are views into the C++ storage, so a snapshot that
    kept them would move with the accumulator and compare equal to anything.
    """
    out = {}
    for field in dataclasses.fields(Tally):
        value = getattr(native, field.name)
        if field.name in _NON_ARRAY_FIELDS:
            out[field.name] = {
                key: np.array(v, copy=True) if isinstance(v, np.ndarray) else v
                for key, v in dict(value).items()
            }
        else:
            out[field.name] = np.array(value, copy=True)
    return out


def _assert_unchanged(before: dict, after: dict, label: str) -> None:
    """Two :func:`_exported` snapshots agree on every bank, dtypes included."""
    for name, want in before.items():
        got = after[name]
        if isinstance(want, dict):
            assert got.keys() == want.keys(), f"{label}: {name} keys"
            for key, expected in want.items():
                if isinstance(expected, np.ndarray):
                    _assert_array(got[key], expected, label, f"{name}[{key!r}]")
                else:
                    assert got[key] == expected, f"{label}: {name}[{key!r}] changed"
            continue
        _assert_array(got, want, label, name)


def test_an_accumulator_merging_itself_is_refused_before_any_bank_is_touched():
    """``merge_from`` refuses ``other is this`` before anything else.

    Every other precondition holds trivially between an accumulator and itself, so without the refusal every bank would
    double, and the held-record bank, appended to itself, would loop forever once it holds a record. The first case
    holds no record, so a missing refusal fails cleanly there (banks doubled) before the second case, which would hang,
    is reached.
    """
    single = [(label, kw) for label, kw in CASES if len(kw.get("hypotheses", UNSPLICED_ONLY)) < 2]
    for label, corpus, holds_records in (
        ("no record held", single, False),
        ("records held", CASES, True),
    ):
        _, accumulator = _pair()
        for _case, kw in corpus:
            accumulator.deposit(**{"hypotheses": UNSPLICED_ONLY, **kw})
        held = dict(accumulator.deferred)["start"].size
        assert (held > 0) == holds_records, f"{label}: {held} records held"
        # the census, the first bank merged, has something a doubling would change
        assert dict(accumulator.qc)["deposited"] > 0
        assert any(dict(accumulator.gap_resolution).values()), f"{label}: an all-zero census"

        before = _exported(accumulator)
        with pytest.raises(ValueError, match="cannot merge itself"):
            accumulator.merge_from(accumulator)
        _assert_unchanged(before, _exported(accumulator), f"a refused self-merge, {label}")


def test_a_merge_across_two_widths_is_refused_before_any_bank_is_touched():
    """``merge_from`` checks ``max_length`` with every other precondition, before merging anything.

    ``max_length`` is the width of the pool and deposited-length histograms. A check placed after the bank merges would
    throw on a half-merged accumulator, and an ``other`` narrower than ``this`` would be read past the end of its
    histogram first. So a refusal must leave every exported bank exactly as it was, whichever side is wider. The wider
    ``other`` comes first so a misplaced check fails on the half-merge before the narrower case can read out of bounds.
    """
    for other_length in (MAX_LENGTH + 500, 120):
        _, this = _pair()
        _, other = _pair(max_length=other_length)
        for _label, kw in CASES:
            this.deposit(**{"hypotheses": UNSPLICED_ONLY, **kw})
            other.deposit(**{"hypotheses": UNSPLICED_ONLY, **kw})
        # `other` holds what a half-merge would have summed into `this` first
        assert dict(other.qc)["deposited"] > 0
        assert dict(other.deferred)["start"].size > 0

        before = _exported(this)
        with pytest.raises(ValueError, match="max_length"):
            this.merge_from(other)
        _assert_unchanged(
            before, _exported(this), f"a refused merge with an other of max_length {other_length}"
        )


def test_an_empty_hypothesis_set_is_the_unspliced_only_set_and_does_not_crash():
    """The specification defaults ``hypotheses`` to :data:`UNSPLICED_ONLY`; the native binding has no default, so a
    caller writing ``hypotheses=()`` offers an empty set. An empty set must mean the same thing as the specification's
    default and deposit identically, rather than reach ``survivors.front()`` on an empty vector. Production never offers
    one (the scanner always offers the genomic path), so only a direct caller can reach it.
    """
    reference, native = _pair()
    kw = dict(start=10, end=90, observed_introns=(), align_strand=int(Strand.POS))

    want = reference.deposit(REF, **kw)  # the reference's own default is UNSPLICED_ONLY
    got = native.deposit(**kw, hypotheses=())
    assert got == want.value, f"empty hypotheses: native {got!r} != reference {want.value!r}"

    reference_default, native_default = _pair()
    reference_default.deposit(REF, **kw, hypotheses=UNSPLICED_ONLY)
    native_default.deposit(**kw, hypotheses=())
    _assert_parity(reference_default, native_default, "empty-hypotheses vs UNSPLICED_ONLY")


def test_a_generator_of_hypotheses_is_the_same_fragment_as_the_tuple():
    """The binding reads ``hypotheses`` in one pass, so a one-shot iterable marshals as its tuple does.

    A second pass over a generator would see nothing, and an empty set deposits as the unspliced fragment alone: an
    implied intron would go uncut, a fragment that must be held would be deposited, and ``length_under`` would answer 0.
    So the specification receives each case's tuple and the native accumulator a generator over it; the two must stay
    byte-identical, and ``length_under`` must give the tuple's length under every hypothesis.
    """
    reference, native = _pair()
    lengths = []
    for label, kw in CASES:
        if "hypotheses" not in kw:
            continue
        hypotheses = kw["hypotheses"]
        want = reference.deposit(REF, **kw)
        got = native.deposit(**{**kw, "hypotheses": (h for h in hypotheses)})
        assert got == want.value, f"{label}: a generator gave {got!r}, the tuple {want.value!r}"
        _assert_parity(reference, native, f"{label}, its hypotheses as a generator")
        for index in range(len(hypotheses)):
            listed = native.length_under(**kw, hypothesis_index=index)
            generated = native.length_under(
                **{**kw, "hypotheses": (h for h in hypotheses)}, hypothesis_index=index
            )
            assert generated == listed, (
                f"{label}: L under hypothesis {index} is {generated} from a generator, "
                f"{listed} from the tuple"
            )
            lengths.append(listed)
    # two-path fragments were held, and the lengths compared are not all zero
    assert reference.tally.qc["deferred_undetermined_gap"] > 0
    assert any(length > 0 for length in lengths)
