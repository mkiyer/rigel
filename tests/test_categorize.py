"""Unit tests for SRD v2 Pass 0 categorization (`rigel.calibration._categorize`)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.buffer import _FinalizedChunk
from rigel.calibration._categorize import (
    FragmentCategory,
    N_CATEGORIES,
    N_STRAND_LABELS,
    StrandLabel,
    categorize_chunk,
)
from rigel.splice import SpliceType


def _make_chunk(fragments: list[dict]) -> _FinalizedChunk:
    """Build a finalized chunk from a list of per-fragment dicts.

    Each dict supplies the SRD v2 strand-bp columns
    (``exon_bp_pos``/``_neg``, ``tx_bp_pos``/``_neg``) plus
    ``splice_type``, ``num_hits``, ``genomic_footprint`` and
    ``read_length``.  The legacy per-candidate arrays
    (``t_indices``/``frag_lengths``/``exon_bp``/``intron_bp``) are
    accepted for fixture convenience but are not consumed by the
    rewritten categorizer.
    """
    n = len(fragments)
    splice_type = np.array([f["splice_type"] for f in fragments], dtype=np.uint8)
    exon_strand = np.array([f.get("exon_strand", 0) for f in fragments], dtype=np.uint8)
    sj_strand = np.array([f.get("sj_strand", 0) for f in fragments], dtype=np.uint8)
    num_hits = np.array([f.get("num_hits", 1) for f in fragments], dtype=np.uint16)
    merge_criteria = np.zeros(n, dtype=np.uint8)
    chimera_type = np.zeros(n, dtype=np.uint8)
    ambig_strand = np.array([f.get("ambig_strand", 0) for f in fragments], dtype=np.uint8)

    t_inds_lists = [f.get("t_indices", []) for f in fragments]
    fl_lists = [f.get("frag_lengths", []) for f in fragments]
    eb_lists = [f.get("exon_bp", []) for f in fragments]
    ib_lists = [f.get("intron_bp", []) for f in fragments]

    t_offsets = np.zeros(n + 1, dtype=np.int32)
    for i, lst in enumerate(t_inds_lists):
        t_offsets[i + 1] = t_offsets[i] + len(lst)

    flat_t = (
        np.concatenate([np.array(x, dtype=np.int32) for x in t_inds_lists])
        if any(t_inds_lists) else np.zeros(0, dtype=np.int32)
    )
    flat_fl = (
        np.concatenate([np.array(x, dtype=np.int32) for x in fl_lists])
        if any(fl_lists) else np.zeros(0, dtype=np.int32)
    )
    flat_eb = (
        np.concatenate([np.array(x, dtype=np.int32) for x in eb_lists])
        if any(eb_lists) else np.zeros(0, dtype=np.int32)
    )
    flat_ib = (
        np.concatenate([np.array(x, dtype=np.int32) for x in ib_lists])
        if any(ib_lists) else np.zeros(0, dtype=np.int32)
    )

    return _FinalizedChunk(
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
        read_length=np.array([f.get("read_length", 100) for f in fragments], dtype=np.uint32),
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
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "tx_bp_pos": 200,
        },
    ])
    out = categorize_chunk(chunk)
    assert not out.keep[0]
    assert out.category[0] == 255  # sentinel
    assert out.counts_by_category_strand.sum() == 0


def test_multimapper_filtered_out():
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "num_hits": 2,
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "tx_bp_pos": 200,
        },
    ])
    out = categorize_chunk(chunk)
    assert not out.keep[0]
    assert out.counts_by_category_strand.sum() == 0


# ---------------------------------------------------------------------------
# Category × strand truth-table cases
# ---------------------------------------------------------------------------


def test_intergenic_zero_strand_bp():
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 175,
            # All four strand-bp columns zero → INTERGENIC.
        },
    ])
    out = categorize_chunk(chunk)
    assert out.keep[0]
    assert out.category[0] == int(FragmentCategory.INTERGENIC)
    assert out.strand[0] == int(StrandLabel.NONE)
    assert out.frag_length[0] == 175


@pytest.mark.parametrize(
    "tx_pos,tx_neg,expected_strand",
    [
        (200, 0, StrandLabel.POS),
        (0, 200, StrandLabel.NEG),
        (200, 200, StrandLabel.AMBIG),
    ],
)
def test_intronic_strand_labels(tx_pos, tx_neg, expected_strand):
    """tx_bp_* > 0 with exon_bp_* == 0 ⇒ INTRONIC, strand from ibp_*."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": 0,
            "exon_bp_neg": 0,
            "tx_bp_pos": tx_pos,
            "tx_bp_neg": tx_neg,
        },
    ])
    out = categorize_chunk(chunk)
    assert out.category[0] == int(FragmentCategory.INTRONIC)
    assert out.strand[0] == int(expected_strand)


@pytest.mark.parametrize(
    "ebp_pos,ebp_neg,expected_strand",
    [
        (200, 0, StrandLabel.POS),
        (0, 200, StrandLabel.NEG),
        (200, 200, StrandLabel.AMBIG),
    ],
)
def test_exon_contained_strand_labels(ebp_pos, ebp_neg, expected_strand):
    """Full exonic coverage (no overhang) ⇒ EXON_CONTAINED, strand from ebp_*."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": ebp_pos,
            "exon_bp_neg": ebp_neg,
            "tx_bp_pos": ebp_pos,
            "tx_bp_neg": ebp_neg,
        },
    ])
    out = categorize_chunk(chunk)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)
    assert out.strand[0] == int(expected_strand)


@pytest.mark.parametrize(
    "ebp_pos,ebp_neg,expected_strand",
    [
        (50, 0, StrandLabel.POS),
        (0, 50, StrandLabel.NEG),
        (50, 50, StrandLabel.AMBIG),
    ],
)
def test_exon_incompatible_strand_labels(ebp_pos, ebp_neg, expected_strand):
    """Some exon overlap but >> tol genomic span outside it ⇒ EXON_INCOMPATIBLE."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": ebp_pos,
            "exon_bp_neg": ebp_neg,
            "tx_bp_pos": max(ebp_pos, 100),
            "tx_bp_neg": max(ebp_neg, 100),
        },
    ])
    out = categorize_chunk(chunk, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(expected_strand)


# ---------------------------------------------------------------------------
# EXON_CONTAINED tolerance edges
# ---------------------------------------------------------------------------


def test_exon_contained_at_tolerance_edge():
    """Exactly tolerance overhang ⇒ EXON_CONTAINED."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": 197,  # 3 bp overhang; tol=3
            "tx_bp_pos": 200,
        },
    ])
    out = categorize_chunk(chunk, exon_fit_tolerance_bp=3)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)


def test_exon_incompatible_one_bp_past_tolerance():
    """One bp past tolerance ⇒ EXON_INCOMPATIBLE."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": 196,  # 4 bp overhang; tol=3
            "tx_bp_pos": 200,
        },
    ])
    out = categorize_chunk(chunk, exon_fit_tolerance_bp=3)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)


# ---------------------------------------------------------------------------
# Strand-overlap zone (a position covered by both + and - transcripts)
# ---------------------------------------------------------------------------


def test_strand_overlap_zone_exon_contained():
    """A 200 bp fragment fully in a region covered by both + and -
    transcripts: ebp_pos == ebp_neg == 200, tbp_* == 200 each.
    Should be EXON_CONTAINED with AMBIG strand label."""
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 200,
            "exon_bp_pos": 200,
            "exon_bp_neg": 200,
            "tx_bp_pos": 200,
            "tx_bp_neg": 200,
        },
    ])
    out = categorize_chunk(chunk)
    assert out.category[0] == int(FragmentCategory.EXON_CONTAINED)
    assert out.strand[0] == int(StrandLabel.AMBIG)


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
            "genomic_footprint": 200,
            "exon_bp_pos": 100,
            "tx_bp_pos": 200,
            "exon_bp_neg": 0,
            "tx_bp_neg": 0,
        },
    ])
    out = categorize_chunk(chunk, exon_fit_tolerance_bp=5)
    # 100 bp exonic + 100 bp intronic on (+) ⇒ EXON_INCOMPATIBLE / POS.
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(StrandLabel.POS)


# ---------------------------------------------------------------------------
# Aggregate output shape & flat tally
# ---------------------------------------------------------------------------


def test_counts_by_category_strand_shape():
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 100,
            "exon_bp_pos": 100,
            "tx_bp_pos": 100,
        },
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "genomic_footprint": 100,
            # all zeros → INTERGENIC
        },
    ])
    out = categorize_chunk(chunk)
    cnt = out.counts_by_category_strand
    assert cnt.shape == (N_CATEGORIES, N_STRAND_LABELS)
    assert cnt.sum() == 2
    assert (
        cnt[int(FragmentCategory.INTERGENIC), int(StrandLabel.NONE)] == 1
    )
    assert (
        cnt[int(FragmentCategory.EXON_CONTAINED), int(StrandLabel.POS)] == 1
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
            "genomic_footprint": 192,
            "exon_bp_pos": 87,
            "tx_bp_pos": 87,  # the overhang lies past the transcript edge
        },
    ])
    out = categorize_chunk(chunk, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)
    assert out.strand[0] == int(StrandLabel.POS)
