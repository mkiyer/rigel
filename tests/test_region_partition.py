"""Tests for the v4 fine-grained calibration region builder."""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration import signature as sig
from rigel.calibration.regions import (
    BoundaryKind,
    RegionType,
    build_fine_region_table,
    classify_boundary_kind,
    validate_against_ref_lengths,
)
from rigel.index import build_genomic_intervals, build_index_artifacts
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _mk_tx(
    t_index: int,
    ref: str,
    strand: Strand,
    exons: list[tuple[int, int]],
    is_synthetic: bool = False,
) -> Transcript:
    return Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(start, end) for start, end in exons],
        t_id=f"t{t_index}",
        g_id=f"g{t_index}",
        t_index=t_index,
        g_index=t_index,
        is_synthetic=is_synthetic,
    )


def _regions(
    ref_length: int,
    transcripts: list[Transcript],
    *,
    ref: str = "chr1",
) -> pd.DataFrame:
    return build_fine_region_table(transcripts, {ref: ref_length})


def test_empty_reference_one_intergenic_region():
    region_df = _regions(1000, [])

    assert region_df["signature"].tolist() == [0x0]
    assert region_df.loc[0, "type"] == int(RegionType.INTERGENIC)
    assert (region_df.loc[0, "start"], region_df.loc[0, "end"]) == (0, 1000)
    assert region_df.loc[0, "left_signature"] == 0xFF
    assert region_df.loc[0, "right_signature"] == 0xFF
    assert region_df.loc[0, "boundary_kind_left"] == int(BoundaryKind.NONE)
    assert region_df.loc[0, "boundary_kind_right"] == int(BoundaryKind.NONE)


def test_single_exon_transcript_emits_exon_signature_with_flanks():
    transcript = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    region_df = _regions(1000, [transcript])

    assert region_df["signature"].tolist() == [0x0, 0x2, 0x0]
    assert list(zip(region_df["start"], region_df["end"], strict=True)) == [
        (0, 100),
        (100, 200),
        (200, 1000),
    ]
    assert bool(region_df.loc[1, "exon_pos"])
    assert region_df.loc[1, "strand"] == sig.coarse_strand_from_signature(0x2)
    assert bool(region_df.loc[1, "boundary_flux_left"])
    assert bool(region_df.loc[1, "boundary_flux_right"])


def test_two_exon_transcript_emits_exon_intron_exon():
    transcript = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)])
    region_df = _regions(1000, [transcript])

    assert region_df["signature"].tolist() == [0x0, 0x2, 0x8, 0x2, 0x0]
    assert list(zip(region_df["start"], region_df["end"], strict=True)) == [
        (0, 100),
        (100, 200),
        (200, 300),
        (300, 400),
        (400, 1000),
    ]
    assert region_df.loc[2, "type"] == int(RegionType.INTRON)
    assert region_df.loc[1, "boundary_kind_right"] == int(BoundaryKind.EXON_INTRON)
    assert region_df.loc[2, "boundary_kind_left"] == int(BoundaryKind.EXON_INTRON)


def test_overlapping_same_strand_exons_merge_identical_signature_runs():
    first = _mk_tx(0, "chr1", Strand.POS, [(100, 250)])
    second = _mk_tx(1, "chr1", Strand.POS, [(150, 300)])
    region_df = _regions(1000, [first, second])

    assert region_df["signature"].tolist() == [0x0, 0x2, 0x0]
    assert tuple(region_df.loc[1, ["start", "end"]]) == (100, 300)


def test_overlapping_opposite_strand_exons_emit_ambiguous_exon_signature():
    pos_tx = _mk_tx(0, "chr1", Strand.POS, [(100, 300)])
    neg_tx = _mk_tx(1, "chr1", Strand.NEG, [(100, 300)])
    region_df = _regions(1000, [pos_tx, neg_tx])

    assert region_df["signature"].tolist() == [0x0, 0x3, 0x0]
    assert region_df.loc[1, "strand"] == sig.coarse_strand_from_signature(0x3)


def test_opposite_strand_exon_intron_overlap_sets_both_bits():
    pos_exon = _mk_tx(0, "chr1", Strand.POS, [(250, 350)])
    neg_intron = _mk_tx(1, "chr1", Strand.NEG, [(100, 200), (400, 500)])
    region_df = _regions(1000, [pos_exon, neg_intron])

    assert 0x6 in region_df["signature"].tolist()
    overlap = region_df[region_df["signature"] == 0x6].iloc[0]
    assert (overlap["start"], overlap["end"]) == (250, 350)
    assert overlap["type"] == int(RegionType.EXON)


def test_same_strand_exon_intron_overlap_sets_both_bits():
    intron_tx = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (400, 500)])
    exon_tx = _mk_tx(1, "chr1", Strand.POS, [(250, 350)])
    region_df = _regions(1000, [intron_tx, exon_tx])

    assert 0xA in region_df["signature"].tolist()
    overlap = region_df[region_df["signature"] == 0xA].iloc[0]
    assert (overlap["start"], overlap["end"]) == (250, 350)
    assert overlap["strand"] == sig.coarse_strand_from_signature(0xA)


def test_synthetic_transcripts_excluded_from_regions():
    real = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    synthetic = _mk_tx(1, "chr1", Strand.POS, [(150, 800)], is_synthetic=True)
    region_df = _regions(1000, [real, synthetic])

    assert region_df["signature"].tolist() == [0x0, 0x2, 0x0]
    assert tuple(region_df.loc[1, ["start", "end"]]) == (100, 200)


def test_build_index_artifacts_schema_and_region_id():
    transcripts = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)]),
        _mk_tx(1, "chr2", Strand.NEG, [(50, 150)]),
    ]
    _, region_df = build_index_artifacts(transcripts, {"chr1": 1000, "chr2": 500})

    assert region_df["region_id"].dtype == np.int64
    assert region_df["start"].dtype == np.int64
    assert region_df["end"].dtype == np.int64
    assert region_df["length"].dtype == np.int64
    assert region_df["signature"].dtype == np.uint8
    assert region_df["type"].dtype == np.uint8
    assert region_df["strand"].dtype == np.uint8
    assert region_df["boundary_flux_left"].dtype == np.bool_
    assert region_df["left_signature"].dtype == np.uint8
    assert region_df["boundary_kind_left"].dtype == np.uint8
    assert isinstance(region_df["ref_name"].dtype, pd.StringDtype)
    assert "tx_pos_bp" not in region_df.columns
    assert "exon_pos_bp" not in region_df.columns
    assert region_df["region_id"].tolist() == list(range(len(region_df)))
    validate_against_ref_lengths(region_df, {"chr1": 1000, "chr2": 500})


def test_neighbor_signatures_and_boundary_kinds_match_adjacent_rows():
    transcript = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)])
    region_df = _regions(1000, [transcript])

    signatures = region_df["signature"].tolist()
    for row_index, row in region_df.iterrows():
        expected_left = 0xFF if row_index == 0 else signatures[row_index - 1]
        expected_right = 0xFF if row_index == len(region_df) - 1 else signatures[row_index + 1]
        assert row["left_signature"] == expected_left
        assert row["right_signature"] == expected_right
        assert row["boundary_kind_left"] == classify_boundary_kind(expected_left, row["signature"])
        assert row["boundary_kind_right"] == classify_boundary_kind(
            row["signature"], expected_right
        )


def test_build_index_artifacts_intervals_unchanged_by_refactor():
    transcripts = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)]),
        _mk_tx(1, "chr1", Strand.POS, [(150, 350)]),
        _mk_tx(2, "chr1", Strand.NEG, [(700, 900)]),
        _mk_tx(3, "chr2", Strand.POS, [(50, 150)]),
        _mk_tx(4, "chr1", Strand.POS, [(120, 880)], is_synthetic=True),
    ]
    transcripts.sort(key=lambda transcript: (transcript.ref, transcript.start, transcript.end))
    ref_lengths = {"chr1": 1500, "chr2": 500}

    intervals_new, _ = build_index_artifacts(transcripts, ref_lengths)
    intervals_legacy = build_genomic_intervals(transcripts, ref_lengths)

    pd.testing.assert_frame_equal(intervals_new, intervals_legacy)
