"""Unit tests for SRD Pass 0 categorization (`rigel.calibration._categorize`)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.buffer import _FinalizedChunk
from rigel.calibration._categorize import (
    FragmentCategory,
    N_CATEGORIES,
    categorize_chunk,
)
from rigel.splice import SpliceType
from rigel.types import Strand


def _make_chunk(fragments: list[dict]) -> _FinalizedChunk:
    """Build a finalized chunk from a list of per-fragment dicts.

    Each dict supplies: splice_type, exon_strand, sj_strand, ambig_strand,
    num_hits, t_indices (list[int]), frag_lengths (parallel list),
    exon_bp (parallel list), intron_bp (parallel list), genomic_start,
    genomic_footprint, read_length, nm.
    """
    n = len(fragments)
    splice_type = np.array([f["splice_type"] for f in fragments], dtype=np.uint8)
    exon_strand = np.array([f["exon_strand"] for f in fragments], dtype=np.uint8)
    sj_strand = np.array([f.get("sj_strand", 0) for f in fragments], dtype=np.uint8)
    num_hits = np.array([f.get("num_hits", 1) for f in fragments], dtype=np.uint16)
    merge_criteria = np.zeros(n, dtype=np.uint8)
    chimera_type = np.zeros(n, dtype=np.uint8)
    ambig_strand = np.array([f.get("ambig_strand", 0) for f in fragments], dtype=np.uint8)

    t_inds_lists = [f["t_indices"] for f in fragments]
    fl_lists = [f["frag_lengths"] for f in fragments]
    eb_lists = [f["exon_bp"] for f in fragments]
    ib_lists = [f["intron_bp"] for f in fragments]

    t_offsets = np.zeros(n + 1, dtype=np.int32)
    for i, lst in enumerate(t_inds_lists):
        t_offsets[i + 1] = t_offsets[i] + len(lst)

    flat_t = np.concatenate([np.array(x, dtype=np.int32) for x in t_inds_lists]) \
        if any(t_inds_lists) else np.zeros(0, dtype=np.int32)
    flat_fl = np.concatenate([np.array(x, dtype=np.int32) for x in fl_lists]) \
        if any(fl_lists) else np.zeros(0, dtype=np.int32)
    flat_eb = np.concatenate([np.array(x, dtype=np.int32) for x in eb_lists]) \
        if any(eb_lists) else np.zeros(0, dtype=np.int32)
    flat_ib = np.concatenate([np.array(x, dtype=np.int32) for x in ib_lists]) \
        if any(ib_lists) else np.zeros(0, dtype=np.int32)

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
        size=n,
    )


# Fixture: 4 transcripts on 2 references.
#   t0: + on ref 0
#   t1: + on ref 0 (overlaps t0 in places)
#   t2: - on ref 0
#   t3: + on ref 1
@pytest.fixture
def tx_strand_arr() -> np.ndarray:
    return np.array([Strand.POS, Strand.POS, Strand.NEG, Strand.POS], dtype=np.int8)


@pytest.fixture
def tx_ref_arr() -> np.ndarray:
    return np.array([0, 0, 0, 1], dtype=np.int32)


def test_spliced_fragment_classified_first(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.SPLICED_ANNOT),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [200],  # would also "fit" — but SPLICED wins
            "intron_bp": [0],
        },
        {
            "splice_type": int(SpliceType.SPLICED_UNANNOT),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [],
            "frag_lengths": [],
            "exon_bp": [],
            "intron_bp": [],
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.SPLICED)
    assert out.category[1] == int(FragmentCategory.SPLICED)


def test_unspliced_sense_exonic(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),  # matches t0
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [150],
            "exon_bp": [150],
            "intron_bp": [0],  # fits inside an exon
            "genomic_footprint": 150,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_SENSE_EXONIC)
    assert out.frag_length[0] == 150
    assert out.ref_id[0] == 0


def test_unspliced_antisense_exonic(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.NEG),  # opposite of t0 (POS)
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [200],
            "intron_bp": [0],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_ANTISENSE_EXONIC)


def test_unspliced_exonic_ambig_via_mixed_strand(tx_strand_arr, tx_ref_arr):
    # ambig_strand flag set by the BAM scanner (overlapping +/- transcripts).
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 1,  # mixed-strand exon overlap
            "t_indices": [0, 2],  # t0 (+), t2 (-)
            "frag_lengths": [180, 180],
            "exon_bp": [180, 180],
            "intron_bp": [0, 0],
            "read_length": 180,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_EXONIC_AMBIG)


def test_unspliced_exonic_ambig_via_unknown_strand(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.NONE),  # unknown read strand
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [200],
            "intron_bp": [0],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_EXONIC_AMBIG)


def test_exon_incompatible_when_overhang_exceeds_tolerance(tx_strand_arr, tx_ref_arr):
    # Read of 200 bp where the only candidate transcript covers only
    # 150 bp of it as exon (50 bp overhangs into intron/intergenic).
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [150],
            "intron_bp": [50],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)


def test_within_tolerance_overhang_treated_as_exonic(tx_strand_arr, tx_ref_arr):
    # 200 bp read, 197 bp exonic, 3 bp slop — within tolerance.
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [197],
            "intron_bp": [3],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_SENSE_EXONIC)


def test_intronic_when_no_exon_overlap(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [0],
            "intron_bp": [200],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.INTRONIC)


def test_intergenic_when_no_candidates(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.NONE),
            "ambig_strand": 0,
            "t_indices": [],
            "frag_lengths": [],
            "exon_bp": [],
            "intron_bp": [],
            "read_length": 175,
            "genomic_footprint": 175,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.category[0] == int(FragmentCategory.INTERGENIC)
    assert out.frag_length[0] == 175  # frag_length sourced from genomic_footprint
    assert out.ref_id[0] == -1


def test_fits_in_one_of_multiple_candidates(tx_strand_arr, tx_ref_arr):
    # Two candidates: t0 doesn't fit (only 180 bp exonic out of 200),
    # t1 fits perfectly (200 bp exonic). Read should still be exonic.
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0, 1],  # both POS
            "frag_lengths": [200, 200],
            "exon_bp": [180, 200],
            "intron_bp": [20, 0],
            "read_length": 200,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_SENSE_EXONIC)


def test_category_count_vector_len_is_seven(tx_strand_arr, tx_ref_arr):
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [100],
            "exon_bp": [100],
            "intron_bp": [0],
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr)
    assert out.n_per_category.shape == (N_CATEGORIES,)
    assert out.n_per_category.sum() == 1


# ---------------------------------------------------------------------------
# Regression tests for the bugs found on real BAMs (vcap_dna80m, 2026-04-24)
# ---------------------------------------------------------------------------


def test_overhang_past_transcript_edge_is_exon_incompatible(tx_strand_arr, tx_ref_arr):
    """A 192 bp gDNA read with a single transcript candidate that
    explains only 87 bp as exon (the rest overhanging the transcript
    edge into intergenic territory) MUST be classified
    EXON_INCOMPATIBLE.

    The original SRD code computed ``min(intron_bp) <= tol`` and
    accepted this fragment as exonic — because ``intron_bp`` only
    counts within-transcript-span bp and missed the 105 bp overhang.
    Observed in real data at chr11:48244998 (i=23482 in the
    dna80m_subset diagnostic).
    """
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [192],
            "exon_bp": [87],
            "intron_bp": [0],
            "read_length": 192,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)


def test_multiple_candidates_all_intron_heavy_is_exon_incompatible(
    tx_strand_arr, tx_ref_arr,
):
    """A read whose every candidate transcript covers it with mostly
    intronic (or overhanging) bp — none reaching exon-fit — MUST be
    EXON_INCOMPATIBLE. The original SRD rule used ``min(intron_bp)``
    instead of ``max(exon_bp)`` and would mis-class this as exonic
    when one candidate happened to have ``intron_bp = 0`` despite a
    huge transcript-edge overhang.
    """
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0, 1],  # both POS
            "frag_lengths": [300, 300],
            # Cand 0: 150 bp exonic, 100 bp intronic, 50 bp overhang.
            # Cand 1:  50 bp exonic,   0 bp intronic, 250 bp overhang.
            "exon_bp": [150, 50],
            "intron_bp": [100, 0],
            "read_length": 300,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    # max(exon_bp) = 150, read_length - 150 = 150 > 5 → not exon-fit.
    assert out.category[0] == int(FragmentCategory.EXON_INCOMPATIBLE)


def test_one_perfect_candidate_rescues_to_exonic(tx_strand_arr, tx_ref_arr):
    """If at least one candidate covers the entire read in exon
    (within tolerance) the fragment is exonic — even if other
    candidates would call it intron-heavy. This is the same isoform-
    overlap pattern observed at chr1:203798676 (i=19583 in the
    dna80m_subset diagnostic), where 21 candidates split 11/10 between
    intronic and exonic interpretations and the read is genuinely
    indistinguishable from RNA.
    """
    chunk = _make_chunk([
        {
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0, 0, 0, 1],  # 3 cands intronic, 1 cand fully exonic
            "frag_lengths": [250, 250, 250, 250],
            "exon_bp": [0, 0, 0, 250],
            "intron_bp": [250, 250, 250, 0],
            "read_length": 250,
        },
    ])
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    assert out.category[0] == int(FragmentCategory.UNSPLICED_SENSE_EXONIC)


def test_real_data_exon_incompatible_dominates_under_high_gdna(
    tx_strand_arr, tx_ref_arr,
):
    """Bulk regression: a 100-fragment chunk where most reads have a
    50 bp transcript-edge overhang (gDNA splash) MUST surface as
    EXON_INCOMPATIBLE under the corrected rule.

    Pre-fix this entire chunk reported ``EXON_INCOMPATIBLE = 0`` because
    ``intron_bp == 0`` for every fragment.
    """
    fragments = []
    # 80 gDNA-like overhang fragments: 50 bp inside an exon, 150 bp overhang.
    for _ in range(80):
        fragments.append({
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [50],
            "intron_bp": [0],
            "read_length": 200,
        })
    # 20 RNA-like fully-exonic fragments.
    for _ in range(20):
        fragments.append({
            "splice_type": int(SpliceType.UNSPLICED),
            "exon_strand": int(Strand.POS),
            "ambig_strand": 0,
            "t_indices": [0],
            "frag_lengths": [200],
            "exon_bp": [200],
            "intron_bp": [0],
            "read_length": 200,
        })
    chunk = _make_chunk(fragments)
    out = categorize_chunk(chunk, tx_strand_arr, tx_ref_arr, exon_fit_tolerance_bp=5)
    counts = out.n_per_category
    assert counts[int(FragmentCategory.EXON_INCOMPATIBLE)] == 80
    assert counts[int(FragmentCategory.UNSPLICED_SENSE_EXONIC)] == 20
