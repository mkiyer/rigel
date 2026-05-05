"""Tests for ``rigel.index._iter_reference_layout``.

The layout iterator is the single source of truth for both
``intervals.feather`` and ``regions.feather``. These tests verify it
yields a perfect tiling of every reference: alternating
INTERGENIC / GENIC spans, no overlaps, no gaps, with synthetic
nRNAs excluded.
"""

from __future__ import annotations


from rigel.index import _GenicSpan, _IntergenicSpan, _iter_reference_layout
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _mk_tx(t_idx: int, ref: str, strand: Strand, exons: list[tuple[int, int]],
           is_synthetic: bool = False) -> Transcript:
    """Construct a minimal Transcript for layout tests."""
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


def _assert_tiles(spans, ref_length: int):
    """Assert that ``spans`` exactly tile [0, ref_length) with alternating types."""
    assert spans, "expected at least one span"
    cursor = 0
    last_type = None
    for s in spans:
        assert s.start == cursor, f"gap or overlap before {s} (cursor={cursor})"
        assert s.end > s.start
        cur_type = type(s).__name__
        assert cur_type != last_type, f"two adjacent {cur_type} spans"
        last_type = cur_type
        cursor = s.end
    assert cursor == ref_length, f"layout ends at {cursor}, expected {ref_length}"


def test_empty_reference_yields_single_intergenic():
    spans = list(_iter_reference_layout(1000, []))
    assert spans == [_IntergenicSpan(0, 1000)]


def test_zero_length_reference_yields_nothing():
    assert list(_iter_reference_layout(0, [])) == []


def test_single_transcript_in_middle_three_spans():
    t = _mk_tx(0, "chr1", Strand.POS, [(100, 200), (300, 400)])
    spans = list(_iter_reference_layout(1000, [t]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 3
    assert spans[0] == _IntergenicSpan(0, 100)
    assert isinstance(spans[1], _GenicSpan)
    assert (spans[1].start, spans[1].end) == (100, 400)
    assert spans[1].transcripts == (t,)
    assert spans[2] == _IntergenicSpan(400, 1000)


def test_transcript_at_left_edge_no_left_intergenic():
    t = _mk_tx(0, "chr1", Strand.POS, [(0, 200)])
    spans = list(_iter_reference_layout(500, [t]))
    _assert_tiles(spans, 500)
    assert len(spans) == 2
    assert isinstance(spans[0], _GenicSpan)
    assert spans[1] == _IntergenicSpan(200, 500)


def test_transcript_at_right_edge_no_right_intergenic():
    t = _mk_tx(0, "chr1", Strand.POS, [(800, 1000)])
    spans = list(_iter_reference_layout(1000, [t]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 2
    assert spans[0] == _IntergenicSpan(0, 800)
    assert isinstance(spans[1], _GenicSpan)


def test_overlapping_transcripts_strand_agnostic_coalesce():
    """Two transcripts on opposite strands that overlap → one genic span."""
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 300)])
    t2 = _mk_tx(1, "chr1", Strand.NEG, [(250, 500)])
    spans = list(_iter_reference_layout(1000, sorted([t1, t2], key=lambda x: (x.start, x.end))))
    _assert_tiles(spans, 1000)
    genic = [s for s in spans if isinstance(s, _GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 500)
    assert {t.t_index for t in genic[0].transcripts} == {0, 1}


def test_disjoint_transcripts_yield_separate_genic_spans():
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    t2 = _mk_tx(1, "chr1", Strand.POS, [(500, 700)])
    spans = list(_iter_reference_layout(1000, [t1, t2]))
    _assert_tiles(spans, 1000)
    assert len(spans) == 5  # int, gen, int, gen, int
    assert sum(isinstance(s, _GenicSpan) for s in spans) == 2


def test_synthetic_excluded_from_layout():
    """Synthetic nRNAs must not coalesce or extend genic spans."""
    real = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    syn = _mk_tx(1, "chr1", Strand.POS, [(150, 800)], is_synthetic=True)
    spans = list(_iter_reference_layout(
        1000, sorted([real, syn], key=lambda x: (x.start, x.end))
    ))
    _assert_tiles(spans, 1000)
    # The synthetic [150,800) must NOT extend the real [100,200) genic span.
    genic = [s for s in spans if isinstance(s, _GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 200)
    assert genic[0].transcripts == (real,)


def test_touching_transcripts_coalesce():
    """t1.end == t2.start → coalesced into one genic span (touching boundary)."""
    t1 = _mk_tx(0, "chr1", Strand.POS, [(100, 200)])
    t2 = _mk_tx(1, "chr1", Strand.POS, [(200, 300)])
    spans = list(_iter_reference_layout(1000, [t1, t2]))
    _assert_tiles(spans, 1000)
    genic = [s for s in spans if isinstance(s, _GenicSpan)]
    assert len(genic) == 1
    assert (genic[0].start, genic[0].end) == (100, 300)
