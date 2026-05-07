"""Tests for the M5-B ``Locus`` / ``MultiLocus`` rename.

Verifies that the new ``Locus`` dataclass is the contiguous-interval
type, that ``MultiLocus`` correctly carries one or more ``Locus``
intervals, and that ``gdna_span`` matches the sum of those intervals'
spans.
"""

import dataclasses

import numpy as np
import pytest

from rigel.locus import Locus, MultiLocus, build_multi_loci


# ---------------------------------------------------------------------------
# Locus (contiguous interval) — frozen dataclass invariants
# ---------------------------------------------------------------------------


def test_locus_is_frozen():
    loc = Locus(ref="chr1", ref_id=0, start=10, end=20)
    with pytest.raises(dataclasses.FrozenInstanceError):
        loc.start = 5  # type: ignore[misc]


def test_locus_span_property():
    loc = Locus(ref="chr2", ref_id=1, start=100, end=350)
    assert loc.span == 250


# ---------------------------------------------------------------------------
# MultiLocus invariants
# ---------------------------------------------------------------------------


def _make_multi_locus(loci_seq, *, multi_locus_id=0):
    return MultiLocus(
        multi_locus_id=multi_locus_id,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.array([0], dtype=np.int32),
        gdna_span=sum(loc.span for loc in loci_seq),
        loci=tuple(loci_seq),
    )


def test_multi_locus_single_interval():
    ml = _make_multi_locus([Locus("chr1", 0, 0, 10000)])
    assert len(ml.loci) == 1
    assert ml.gdna_span == ml.loci[0].span == 10000


def test_multi_locus_multi_ref():
    """Paralog cluster spanning two refs produces multiple Locus intervals."""
    loci = (
        Locus("chr1", 0, 100, 600),
        Locus("chr2", 1, 50, 250),
    )
    ml = _make_multi_locus(loci)
    assert len(ml.loci) == 2
    assert ml.gdna_span == sum(loc.span for loc in ml.loci) == 700
    # ref_id ordering assertion: the constructor accepts them sorted by
    # (ref_id, start); confirm we can iterate in declared order.
    assert [loc.ref for loc in ml.loci] == ["chr1", "chr2"]


# ---------------------------------------------------------------------------
# build_multi_loci integration
# ---------------------------------------------------------------------------


def test_build_multi_loci_gdna_span_matches_sum_of_loci(mini_index):
    """``gdna_span`` should equal Σ(end-start) over ``loci`` (within max(_,1))."""
    from rigel.scored_fragments import ScoredFragments

    # Synthetic em_data: one ambiguous unit linking t0 and t1 (both on chr1).
    em_data = ScoredFragments(
        offsets=np.array([0, 2], dtype=np.int64),
        t_indices=np.array([0, 1], dtype=np.int32),
        log_liks=np.zeros(2, dtype=np.float64),
        count_cols=np.zeros(2, dtype=np.uint8),
        coverage_weights=np.ones(2, dtype=np.float64),
        locus_t_indices=np.array([0], dtype=np.int32),
        locus_count_cols=np.array([0], dtype=np.uint8),
        is_spliced=np.zeros(1, dtype=bool),
        gdna_log_liks=np.zeros(1, dtype=np.float64),
        frag_ids=np.array([0], dtype=np.int64),
        frag_class=np.zeros(1, dtype=np.int8),
        splice_type=np.zeros(1, dtype=np.uint8),
        n_units=1,
        n_candidates=2,
    )
    multi_loci = build_multi_loci(em_data, mini_index)
    assert len(multi_loci) > 0
    for ml in multi_loci:
        loci_sum = sum(loc.span for loc in ml.loci)
        # build_multi_loci floors gdna_span at 1 to avoid downstream div-by-zero.
        assert ml.gdna_span == max(loci_sum, 1)
        for loc in ml.loci:
            assert isinstance(loc, Locus)
            assert loc.start <= loc.end
