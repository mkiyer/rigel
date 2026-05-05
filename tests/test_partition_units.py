"""Tests for ``partition_units_to_loci`` (anchor-transcript routing)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._locus_n_obs import (
    build_t_to_local_locus,
    partition_units_to_loci,
)
from rigel.locus import Locus, MultiLocus


class _MiniEM:
    """Minimal stand-in for ScoredFragments."""

    def __init__(self, locus_t_indices: np.ndarray):
        self.locus_t_indices = locus_t_indices


def _make_locus(rid: int, start: int, end: int) -> Locus:
    return Locus(ref=f"chr{rid}", ref_id=rid, start=start, end=end)


def _make_ml(
    multi_locus_id: int,
    transcript_indices: np.ndarray,
    unit_indices: np.ndarray,
    loci: tuple[Locus, ...],
) -> MultiLocus:
    span = sum(loc.span for loc in loci)
    return MultiLocus(
        multi_locus_id=multi_locus_id,
        transcript_indices=transcript_indices.astype(np.int32),
        unit_indices=unit_indices.astype(np.int32),
        gdna_span=int(max(span, 1)),
        loci=loci,
    )


def test_single_locus_fast_path_returns_identity():
    ml = _make_ml(
        0,
        np.array([3, 7]),
        np.array([10, 11, 12]),
        (_make_locus(0, 100, 500),),
    )
    em = _MiniEM(locus_t_indices=np.array([3, 7, 3, 7], dtype=np.int32))
    parts = partition_units_to_loci(ml, em, np.zeros(0, dtype=np.int8))
    assert len(parts) == 1
    # Fast path: returns ml.unit_indices itself (no copy).
    assert parts[0] is ml.unit_indices


def test_two_locus_all_anchored_to_first():
    # Transcripts 3, 7 live in Locus 0; transcript 5 lives in Locus 1.
    ml = _make_ml(
        0,
        np.array([3, 5, 7]),
        np.array([10, 11, 12]),
        (_make_locus(0, 100, 500), _make_locus(0, 1000, 2000)),
    )
    # All 3 units anchored to transcript 3 → Locus 0.
    em = _MiniEM(locus_t_indices=np.array([0, 0, 0, 3, 3, 3, 0, 7, 0, 0, 3, 3, 3], dtype=np.int32))
    t_starts = np.array([0, 0, 0, 200, 0, 1500, 0, 300], dtype=np.int64)
    t_ref = np.array([0] * 8, dtype=np.int32)
    t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
    # transcript 3 (start 200) ⇒ locus 0; transcript 5 (start 1500) ⇒ locus 1; transcript 7 ⇒ locus 0.
    assert t_to_local.tolist() == [0, 1, 0]
    parts = partition_units_to_loci(ml, em, t_to_local)
    assert len(parts) == 2
    assert parts[0].size == 3
    assert parts[1].size == 0


def test_two_locus_mixed_anchors_routed_correctly():
    ml = _make_ml(
        0,
        np.array([3, 5, 7]),
        np.array([10, 11, 12, 13]),
        (_make_locus(0, 100, 500), _make_locus(0, 1000, 2000)),
    )
    # Anchors: unit 0 → t=3 (loc 0); unit 1 → t=5 (loc 1); unit 2 → t=7 (loc 0); unit 3 → t=5 (loc 1).
    em = _MiniEM(
        locus_t_indices=np.array(
            [0, 0, 0, 3, 0, 5, 0, 7, 0, 0, 3, 5, 7, 5], dtype=np.int32
        )
    )
    t_starts = np.array([0, 0, 0, 200, 0, 1500, 0, 300], dtype=np.int64)
    t_ref = np.array([0] * 8, dtype=np.int32)
    t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
    parts = partition_units_to_loci(ml, em, t_to_local)
    assert sorted(parts[0].tolist()) == [10, 12]
    assert sorted(parts[1].tolist()) == [11, 13]


def test_anchor_outside_transcript_indices_raises():
    ml = _make_ml(
        0,
        np.array([3, 5]),
        np.array([10, 11]),
        (_make_locus(0, 100, 500), _make_locus(0, 1000, 2000)),
    )
    # Anchor 99 not in [3, 5].
    em = _MiniEM(locus_t_indices=np.array([0] * 11 + [99], dtype=np.int32))
    em.locus_t_indices[10] = 99
    em.locus_t_indices[11] = 3
    t_starts = np.array([0, 0, 0, 200, 0, 1500], dtype=np.int64)
    t_ref = np.array([0] * 6, dtype=np.int32)
    t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
    with pytest.raises(RuntimeError, match="invariant violation"):
        partition_units_to_loci(ml, em, t_to_local)


def test_slow_path_output_dtype_int32():
    ml = _make_ml(
        0,
        np.array([3, 7]),
        np.array([10, 11], dtype=np.int64),  # input int64
        (_make_locus(0, 100, 500), _make_locus(0, 1000, 2000)),
    )
    em = _MiniEM(locus_t_indices=np.array([0] * 12, dtype=np.int32))
    em.locus_t_indices[10] = 3
    em.locus_t_indices[11] = 7
    t_starts = np.array([0, 0, 0, 200, 0, 0, 0, 1500], dtype=np.int64)
    t_ref = np.array([0] * 8, dtype=np.int32)
    t_to_local = build_t_to_local_locus(ml, t_starts, t_ref)
    parts = partition_units_to_loci(ml, em, t_to_local)
    for arr in parts:
        assert arr.dtype == np.int32
