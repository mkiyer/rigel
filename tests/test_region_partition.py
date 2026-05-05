"""Tests for ``rigel.calibration.regions.emit_regions`` and the
``regions.feather`` artifact produced by ``build_index_artifacts``.

Schema and semantics are locked in
``docs/calibration/calibration_v6_plan.md`` §2.2.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.regions import (
    RegionRecord,
    RegionStrand,
    RegionType,
    emit_regions,
)
from rigel.index import (
    _GenicSpan,
    _IntergenicSpan,
    _iter_reference_layout,
    build_index_artifacts,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _mk_tx(t_idx: int, ref: str, strand: Strand, exons: list[tuple[int, int]],
           is_synthetic: bool = False) -> Transcript:
    return Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(s, e) for s, e in exons],
        t_id=f"t{t_idx}",
        g_id=f"g{t_idx}",
        t_index=t_idx,
        g_index=t_idx,
        is_synthetic=is_synthetic,
    )


def _regions(ref_length, transcripts, ref="chr1"):
    layout = _iter_reference_layout(ref_length, transcripts)
    return list(emit_regions(ref, layout))


def test_empty_reference_one_intergenic_region():
    regs = _regions(1000, [])
    assert len(regs) == 1
    r = regs[0]
    assert r.type == RegionType.INTERGENIC
    assert (r.start, r.end) == (0, 1000)
    assert r.strand == RegionStrand.NONE
    assert r.tx_pos_bp == r.tx_neg_bp == 0
    assert r.exon_pos_bp == r.exon_neg_bp == 0
    assert r.boundary_flux_left is False and r.boundary_flux_right is False


def test_single_exon_transcript_no_boundary_flux():
    """Single-exon (+) transcript: one EXON region, both bf flags False."""
    t = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    regs = _regions(1000, [t])
    # Expect: INTERGENIC[0,100), EXON[100,200), INTERGENIC[200,1000)
    assert [r.type for r in regs] == [
        RegionType.INTERGENIC, RegionType.EXON, RegionType.INTERGENIC,
    ]
    exon = regs[1]
    assert (exon.start, exon.end) == (100, 200)
    assert exon.strand == RegionStrand.POS
    assert exon.exon_pos_bp == 100 and exon.exon_neg_bp == 0
    assert exon.tx_pos_bp == 100 and exon.tx_neg_bp == 0
    assert exon.boundary_flux_left is False
    assert exon.boundary_flux_right is False


def test_two_exon_transcript_internal_boundaries_flagged():
    """Two-exon (+) transcript: EXON[0] has bfr=True; EXON[1] has bfl=True."""
    t = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)])
    regs = _regions(1000, [t])
    # Within genic span [100,400): EXON[100,200), INTRON[200,300), EXON[300,400)
    types = [r.type for r in regs]
    assert types == [
        RegionType.INTERGENIC, RegionType.EXON, RegionType.INTRON,
        RegionType.EXON, RegionType.INTERGENIC,
    ]
    e1, intron, e2 = regs[1], regs[2], regs[3]
    assert (e1.boundary_flux_left, e1.boundary_flux_right) == (False, True)
    assert (e2.boundary_flux_left, e2.boundary_flux_right) == (True, False)
    assert (intron.boundary_flux_left, intron.boundary_flux_right) == (False, False)
    assert intron.strand == RegionStrand.POS
    assert intron.tx_pos_bp == 100 and intron.exon_pos_bp == 0


def test_three_exon_transcript_middle_exon_both_flagged():
    """Three-exon (+) transcript: middle EXON has both bf flags True."""
    t = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400), (500, 600)])
    regs = _regions(1000, [t])
    types = [r.type for r in regs]
    assert types == [
        RegionType.INTERGENIC,
        RegionType.EXON, RegionType.INTRON, RegionType.EXON, RegionType.INTRON, RegionType.EXON,
        RegionType.INTERGENIC,
    ]
    middle = regs[3]
    assert (middle.start, middle.end) == (300, 400)
    assert middle.boundary_flux_left is True
    assert middle.boundary_flux_right is True
    # Terminal exons keep one flag False (the TSS / TES side).
    assert regs[1].boundary_flux_left is False  # TSS
    assert regs[5].boundary_flux_right is False  # TES


def test_overlapping_strands_in_exon_marks_ambig():
    """Overlapping (+)-exon and (-)-exon at same coords → strand=AMBIG, both bp counted."""
    tp = _mk_tx(0, "chr1", Strand.POS, [(100, 300)])
    tn = _mk_tx(1, "chr1", Strand.NEG, [(100, 300)])
    regs = _regions(1000, [tp, tn])
    exon = next(r for r in regs if r.type == RegionType.EXON)
    assert exon.strand == RegionStrand.AMBIG
    assert exon.exon_pos_bp == 200
    assert exon.exon_neg_bp == 200
    assert exon.tx_pos_bp == 200
    assert exon.tx_neg_bp == 200


def test_exon_wins_over_intron_at_overlap():
    """Where one transcript's exon overlaps another's intron, EXON wins."""
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (400, 500)])  # intron 200-400
    t2 = _mk_tx(1, "chr1", Strand.POS, [(250, 350)])               # inside the intron
    regs = _regions(1000, sorted([t1, t2], key=lambda x: (x.start, x.end)))
    # Genic span [100,500). Within it:
    #   EXON[100,200) (t1.exon0)
    #   INTRON[200,250) (t1 intron, t2 not yet)
    #   EXON[250,350) (t2.exon)
    #   INTRON[350,400) (t1 intron only)
    #   EXON[400,500) (t1.exon1)
    types = [r.type for r in regs if r.type != RegionType.INTERGENIC]
    assert types == [
        RegionType.EXON, RegionType.INTRON, RegionType.EXON,
        RegionType.INTRON, RegionType.EXON,
    ]


def test_synthetic_transcripts_excluded_from_regions():
    """Synthetics must not contribute to exon or transcript bp counts."""
    real = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    syn = _mk_tx(1, "chr1", Strand.POS, [(150, 800)], is_synthetic=True)
    regs = _regions(1000, sorted([real, syn], key=lambda x: (x.start, x.end)))
    # The synthetic must not extend the genic span beyond [100,200).
    exon_regs = [r for r in regs if r.type == RegionType.EXON]
    assert len(exon_regs) == 1
    assert (exon_regs[0].start, exon_regs[0].end) == (100, 200)
    assert exon_regs[0].exon_pos_bp == 100   # only `real` counted


def test_region_partition_tiles_reference_exactly():
    """Regions must tile [0, ref_length) with no gaps, no overlaps."""
    txs = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (400, 500)]),
        _mk_tx(1, "chr1", Strand.NEG, [(800, 900)]),
    ]
    regs = _regions(2000, sorted(txs, key=lambda x: (x.start, x.end)))
    cursor = 0
    for r in regs:
        assert r.start == cursor
        cursor = r.end
    assert cursor == 2000


def test_build_index_artifacts_schema_and_region_id():
    """End-to-end: build_index_artifacts produces a region_df with locked dtypes
    and a globally monotonic region_id."""
    txs = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)]),
        _mk_tx(1, "chr2", Strand.NEG, [(50, 150)]),
    ]
    iv_df, region_df = build_index_artifacts(txs, {"chr1": 1000, "chr2": 500})

    # dtypes
    assert region_df["region_id"].dtype == np.int64
    assert region_df["start"].dtype == np.int64
    assert region_df["end"].dtype == np.int64
    assert region_df["type"].dtype == np.uint8
    assert region_df["strand"].dtype == np.uint8
    assert region_df["tx_pos_bp"].dtype == np.int64
    assert region_df["exon_pos_bp"].dtype == np.int64
    assert region_df["boundary_flux_left"].dtype == np.bool_
    assert region_df["boundary_flux_right"].dtype == np.bool_
    # ref_name is StringDtype, not object
    assert isinstance(region_df["ref_name"].dtype, pd.StringDtype)

    # region_id is 0..N-1 and matches DataFrame row order
    assert list(region_df["region_id"]) == list(range(len(region_df)))

    # Tiling per ref
    for ref, ref_len in [("chr1", 1000), ("chr2", 500)]:
        sub = region_df[region_df["ref_name"] == ref].sort_values("start")
        cursor = 0
        for _, row in sub.iterrows():
            assert row["start"] == cursor
            cursor = row["end"]
        assert cursor == ref_len


def test_build_index_artifacts_intervals_unchanged_by_refactor(tmp_path_factory):
    """``intervals.feather`` must be byte-identical to the legacy generator output.

    Verified by comparing ``build_index_artifacts(...)[0]`` against
    ``build_genomic_intervals(...)`` on the same inputs.
    """
    from rigel.index import build_genomic_intervals

    txs = [
        _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)]),
        _mk_tx(1, "chr1", Strand.POS, [(150, 350)]),
        _mk_tx(2, "chr1", Strand.NEG, [(700, 900)]),
        _mk_tx(3, "chr2", Strand.POS, [(50, 150)]),
        # Synthetic: must NOT change intervals.feather.
        _mk_tx(4, "chr1", Strand.POS, [(120, 880)], is_synthetic=True),
    ]
    txs.sort(key=lambda t: (t.ref, t.start, t.end))
    ref_lengths = {"chr1": 1500, "chr2": 500}

    iv_new, _ = build_index_artifacts(txs, ref_lengths)
    iv_legacy = build_genomic_intervals(txs, ref_lengths)

    pd.testing.assert_frame_equal(iv_new, iv_legacy)
