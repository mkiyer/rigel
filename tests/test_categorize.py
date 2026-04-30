"""Unit tests for SRD v2 Pass 0 categorization (`rigel.calibration._categorize`)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.buffer import _FinalizedChunk
from rigel.calibration._categorize import (
    FragmentCategory,
    FragmentStrand,
    N_CATEGORIES,
    N_FRAGMENT_STRANDS,
    categorize_chunk,
)
from rigel.splice import SpliceType
from rigel.types import Strand


def _make_chunk(fragments: list[dict]) -> _FinalizedChunk:
    """Build a finalized chunk from a list of per-fragment dicts.

    Each dict supplies the SRD v2 strand-bp columns
    (``exon_bp_pos``/``_neg``, ``tx_bp_pos``/``_neg``) plus
    ``splice_type``, ``num_hits``, ``genomic_footprint`` and
    ``read_length``.

    SRD v3 update: the categorizer's EXON_CONTAINED check is now
    per-transcript (uses ``max_T(exon_bp[T])``), not the aggregate
    ``exon_bp_pos/_neg`` scalars.  To keep the legacy scalar-driven
    fixtures working, this helper synthesizes a tiny per-fragment
    candidate set such that ``max_T(exon_bp[T on +]) == exon_bp_pos``
    and ``max_T(exon_bp[T on -]) == exon_bp_neg``.  The synthesized
    candidates use private transcript indices ``0`` (POS strand) and
    ``1`` (NEG strand), and the helper returns a chunk plus the
    matching ``t_strand_arr`` via the module-level
    ``_LAST_T_STRAND_ARR`` global so the test driver
    ``_categorize`` can look it up.

    Tests that drive the per-T check directly (rather than via the
    legacy scalars) may supply explicit ``t_indices`` /
    ``exon_bp`` lists; those override the synthesis path.

    Default behaviour: ``read_length`` falls back to
    ``genomic_footprint`` when not explicitly set, matching the
    implicit assumption of the legacy fixtures (no inner mate gap).
    """
    n = len(fragments)
    splice_type = np.array([f["splice_type"] for f in fragments], dtype=np.uint8)
    exon_strand = np.array([f.get("exon_strand", 0) for f in fragments], dtype=np.uint8)
    sj_strand = np.array([f.get("sj_strand", 0) for f in fragments], dtype=np.uint8)
    num_hits = np.array([f.get("num_hits", 1) for f in fragments], dtype=np.uint16)
    merge_criteria = np.zeros(n, dtype=np.uint8)
    chimera_type = np.zeros(n, dtype=np.uint8)
    ambig_strand = np.array([f.get("ambig_strand", 0) for f in fragments], dtype=np.uint8)

    # Per-T synthesis: for each fragment, build a candidate set whose
    # strand-restricted max(exon_bp) equals the requested scalar.
    auto_t_inds: list[list[int]] = []
    auto_exon_bp: list[list[int]] = []
    for f in fragments:
        if "t_indices" in f or "exon_bp" in f:
            auto_t_inds.append(list(f.get("t_indices", [])))
            auto_exon_bp.append(list(f.get("exon_bp", [])))
            continue
        ti: list[int] = []
        eb: list[int] = []
        if int(f.get("exon_bp_pos", 0)) > 0:
            ti.append(0)
            eb.append(int(f["exon_bp_pos"]))
        if int(f.get("exon_bp_neg", 0)) > 0:
            ti.append(1)
            eb.append(int(f["exon_bp_neg"]))
        auto_t_inds.append(ti)
        auto_exon_bp.append(eb)
    fl_lists = [f.get("frag_lengths", auto_t_inds[i]) for i, f in enumerate(fragments)]
    ib_lists = [f.get("intron_bp", [0] * len(auto_t_inds[i])) for i, f in enumerate(fragments)]

    t_offsets = np.zeros(n + 1, dtype=np.int32)
    for i, lst in enumerate(auto_t_inds):
        t_offsets[i + 1] = t_offsets[i] + len(lst)

    def _flat(lists, default_dtype):
        if any(lists):
            return np.concatenate([np.array(x, dtype=default_dtype) for x in lists])
        return np.zeros(0, dtype=default_dtype)

    flat_t = _flat(auto_t_inds, np.int32)
    flat_fl = _flat(fl_lists, np.int32)
    flat_eb = _flat(auto_exon_bp, np.int32)
    flat_ib = _flat(ib_lists, np.int32)

    chunk = _FinalizedChunk(
        splice_type=splice_type,
        exon_strand=exon_strand,
        sj_strand=sj_strand,
        num_hits=num_hits,
        merge_criteria=merge_criteria,
        chimera_type=chimera_type,
        t_offsets=t_offsets,
        t_indices=flat_t,
        frag_lengths=flat_fl,
        exon_bp=flat_eb,
        intron_bp=flat_ib,
        ambig_strand=ambig_strand,
        frag_id=np.arange(n, dtype=np.int64),
        read_length=np.array(
            [f.get("read_length", f.get("genomic_footprint", 100)) for f in fragments],
            dtype=np.uint32,
        ),
        genomic_footprint=np.array(
            [f.get("genomic_footprint", 100) for f in fragments], dtype=np.int32
        ),
        genomic_start=np.array([f.get("genomic_start", 0) for f in fragments], dtype=np.int32),
        nm=np.zeros(n, dtype=np.uint16),
        exon_bp_pos=np.array([f.get("exon_bp_pos", 0) for f in fragments], dtype=np.int32),
        exon_bp_neg=np.array([f.get("exon_bp_neg", 0) for f in fragments], dtype=np.int32),
        tx_bp_pos=np.array([f.get("tx_bp_pos", 0) for f in fragments], dtype=np.int32),
        tx_bp_neg=np.array([f.get("tx_bp_neg", 0) for f in fragments], dtype=np.int32),
        size=n,
    )
    return chunk


# Synthetic transcript-strand lookup matching the per-T synthesis in
# ``_make_chunk``: t_index 0 → POS strand, t_index 1 → NEG strand.
_T_STRAND_ARR = np.array([int(Strand.POS), int(Strand.NEG)], dtype=np.int32)


def _categorize(chunk, **kwargs):
    """Test driver: forwards to ``categorize_chunk`` with the synthetic
    ``t_strand_arr`` understood by ``_make_chunk``.
    """
    kwargs.setdefault("t_strand_arr", _T_STRAND_ARR)
    return categorize_chunk(chunk, **kwargs)


# ---------------------------------------------------------------------------
# Pre-filter: only unique-mapper UNSPLICED fragments get categorised
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "splice_type",
    [
        SpliceType.SPLICED_ANNOT,
        SpliceType.SPLICED_UNANNOT,
        SpliceType.SPLICED_IMPLICIT,
        SpliceType.SPLICE_ARTIFACT,
    ],
)
def test_non_unspliced_filtered_out(splice_type):
    chunk = _make_chunk([
        {
            "splice_type": int(splice_type),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "tx_bp_pos": 200,
        },
    ])
    out = _categorize(chunk)
    assert not out.keep[0]
    assert out.category[0] == 255  # sentinel
    assert out.counts_by_category_strand.sum() == 0


def test_multimapper_filtered_out():
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "num_hits": 2,
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "tx_bp_pos": 200,
        },
    ])
    out = _categorize(chunk)
    assert not out.keep[0]
    assert out.counts_by_category_strand.sum() == 0


@pytest.mark.parametrize(
    "read_strand",
    [int(Strand.NONE), int(Strand.AMBIGUOUS)],
)
def test_unknown_or_ambig_read_strand_filtered_out(read_strand):
    """SRD v3 Phase 1: a fragment whose read strand is not POS or NEG
    cannot be expressed in transcript frame and is dropped from
    calibration entirely."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": read_strand,
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "tx_bp_pos": 200,
        },
    ])
    out = _categorize(chunk)
    assert not out.keep[0]
    assert out.counts_by_category_strand.sum() == 0


# ---------------------------------------------------------------------------
# Category × strand truth-table cases
# ---------------------------------------------------------------------------


def test_intergenic_zero_strand_bp():
    """INTERGENIC: read strand POS by convention maps to SENSE."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 175,
            # All four strand-bp columns zero → INTERGENIC.
        },
    ])
    out = _categorize(chunk)
    assert out.keep[0]
    assert out.category[0] == int(FragmentCategory.INTERGENIC)
    assert out.strand[0] == int(FragmentStrand.SENSE)
    assert out.frag_length[0] == 175


@pytest.mark.parametrize(
    "read_strand,expected_label",
    [
        (Strand.POS, FragmentStrand.SENSE),
        (Strand.NEG, FragmentStrand.ANTISENSE),
    ],
)
def test_intergenic_strand_naming_convention(read_strand, expected_label):
    """For INTERGENIC fragments the SENSE/ANTISENSE label is a pure
    convention based on the read's genomic strand."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(read_strand),
            "genomic_footprint": 175,
        },
    ])
    out = _categorize(chunk)
    assert out.category[0] == int(FragmentCategory.INTERGENIC)
    assert out.strand[0] == int(expected_label)


# ---------------------------------------------------------------------------
# Transcript-frame strand labeling: the core SRD v3 Phase 1 logic
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "tx_pos,tx_neg,read_strand,expected_label",
    [
        # tx on (+) only:
        (200, 0, Strand.POS, FragmentStrand.SENSE),
        (200, 0, Strand.NEG, FragmentStrand.ANTISENSE),
        # tx on (-) only:
        (0, 200, Strand.POS, FragmentStrand.ANTISENSE),
        (0, 200, Strand.NEG, FragmentStrand.SENSE),
        # tx on both strands → AMBIG regardless of read strand:
        (200, 200, Strand.POS, FragmentStrand.AMBIG),
        (200, 200, Strand.NEG, FragmentStrand.AMBIG),
    ],
)
def test_intronic_strand_labels(tx_pos, tx_neg, read_strand, expected_label):
    """INTRONIC: label by comparing read strand vs derived ibp_*."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(read_strand),
            "genomic_footprint": 200,
            "exon_bp_pos": 0,
            "exon_bp_neg": 0,
            "tx_bp_pos": tx_pos,
            "tx_bp_neg": tx_neg,
        },
    ])
    out = _categorize(chunk)
    assert out.category[0] == int(FragmentCategory.INTRONIC)
    assert out.strand[0] == int(expected_label)


@pytest.mark.parametrize(
    "ebp_pos,ebp_neg,read_strand,expected_label",
    [
        (200, 0, Strand.POS, FragmentStrand.SENSE),
        (200, 0, Strand.NEG, FragmentStrand.ANTISENSE),
        (0, 200, Strand.POS, FragmentStrand.ANTISENSE),
        (0, 200, Strand.NEG, FragmentStrand.SENSE),
        (200, 200, Strand.POS, FragmentStrand.AMBIG),
        (200, 200, Strand.NEG, FragmentStrand.AMBIG),
    ],
)
def test_exon_contained_strand_labels(ebp_pos, ebp_neg, read_strand, expected_label):
    """EXON_CONTAINED: label by comparing read strand vs ebp_*."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(read_strand),
            "genomic_footprint": 200,
            "exon_bp_pos": ebp_pos,
            "exon_bp_neg": ebp_neg,
            "tx_bp_pos": ebp_pos,
            "tx_bp_neg": ebp_neg,
        },
    ])
    out = _categorize(chunk)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)
    assert out.strand[0] == int(expected_label)


@pytest.mark.parametrize(
    "ebp_pos,ebp_neg,read_strand,expected_label",
    [
        (50, 0, Strand.POS, FragmentStrand.SENSE),
        (50, 0, Strand.NEG, FragmentStrand.ANTISENSE),
        (0, 50, Strand.POS, FragmentStrand.ANTISENSE),
        (0, 50, Strand.NEG, FragmentStrand.SENSE),
        (50, 50, Strand.POS, FragmentStrand.AMBIG),
        (50, 50, Strand.NEG, FragmentStrand.AMBIG),
    ],
)
def test_exon_incompatible_strand_labels(ebp_pos, ebp_neg, read_strand, expected_label):
    """EXON_INCOMPATIBLE: label by comparing read strand vs ebp_*."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(read_strand),
            "genomic_footprint": 200,
            "exon_bp_pos": ebp_pos,
            "exon_bp_neg": ebp_neg,
            "tx_bp_pos": max(ebp_pos, 100),
            "tx_bp_neg": max(ebp_neg, 100),
        },
    ])
    out = _categorize(chunk, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(expected_label)


# ---------------------------------------------------------------------------
# EXON_CONTAINED tolerance edges
# ---------------------------------------------------------------------------


def test_exon_contained_at_tolerance_edge():
    """Exactly tolerance overhang ⇒ EXON_CONTAINED."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 200,
            "exon_bp_pos": 197,  # 3 bp overhang; tol=3
            "tx_bp_pos": 200,
        },
    ])
    out = _categorize(chunk, exon_fit_tolerance_bp=3)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)


def test_exon_incompatible_one_bp_past_tolerance():
    """One bp past tolerance ⇒ EXON_INCOMPATIBLE."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 200,
            "exon_bp_pos": 196,  # 4 bp overhang; tol=3
            "tx_bp_pos": 200,
        },
    ])
    out = _categorize(chunk, exon_fit_tolerance_bp=3)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)


# ---------------------------------------------------------------------------
# Strand-overlap zone (a position covered by both + and - transcripts)
# ---------------------------------------------------------------------------


def test_strand_overlap_zone_exon_contained():
    """A 200 bp fragment fully in a region covered by both + and -
    transcripts: ebp_pos == ebp_neg == 200, tbp_* == 200 each.
    Should be EXON_CONTAINED with AMBIG strand label regardless of
    read strand (sense to one tx, antisense to the other)."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "exon_bp_neg": 200,
            "tx_bp_pos": 200,
            "tx_bp_neg": 200,
        },
    ])
    out = _categorize(chunk)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)
    assert out.strand[0] == int(FragmentStrand.AMBIG)


# ---------------------------------------------------------------------------
# Invariants: rows × strand-bp columns
# ---------------------------------------------------------------------------


def test_invariants_exon_le_tx():
    """The C++ scanner must satisfy exon_bp_strand <= tx_bp_strand on
    every fragment.  We don't enforce it in `categorize_chunk` but
    `intron_bp_strand = max(tx - exon, 0)` clamps to non-negative and
    classification still produces a valid (category, strand) output."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 200,
            "exon_bp_pos": 100,
            "tx_bp_pos": 200,
            "exon_bp_neg": 0,
            "tx_bp_neg": 0,
        },
    ])
    out = _categorize(chunk, exon_fit_tolerance_bp=5)
    # 100 bp exonic + 100 bp intronic on (+); read on (+) ⇒
    # EXON_INCOMPATIBLE / SENSE.
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(FragmentStrand.SENSE)


# ---------------------------------------------------------------------------
# Aggregate output shape & flat tally
# ---------------------------------------------------------------------------


def test_counts_by_category_strand_shape():
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 100,
            "exon_bp_pos": 100,
            "tx_bp_pos": 100,
        },
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 100,
            # all zeros → INTERGENIC
        },
    ])
    out = _categorize(chunk)
    cnt = out.counts_by_category_strand
    assert cnt.shape == (N_CATEGORIES, N_FRAGMENT_STRANDS)
    assert cnt.sum() == 2
    assert (
        cnt[int(FragmentCategory.INTERGENIC), int(FragmentStrand.SENSE)] == 1
    )
    assert (
        cnt[int(FragmentCategory.EXON_CONTAINED), int(FragmentStrand.SENSE)] == 1
    )


# ---------------------------------------------------------------------------
# Regression: the real-data overhang bug from SRD v1
# (a 192 bp gDNA read where exon_bp = 87, tx_bp = 192 must NOT be exonic)
# ---------------------------------------------------------------------------


def test_overhang_past_transcript_edge_is_exon_incompatible():
    """SRD v1 mis-classified this as exonic (using `min(intron_bp) <= tol`).
    The new geometric rule
    `genomic_footprint - max(exon_bp_pos, exon_bp_neg) <= tol` correctly
    rejects it: 192 - 87 = 105 bp >> tol.
    """
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "genomic_footprint": 192,
            "exon_bp_pos": 87,
            "tx_bp_pos": 87,  # the overhang lies past the transcript edge
        },
    ])
    out = _categorize(chunk, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(FragmentStrand.SENSE)
