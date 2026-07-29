"""RegionArrays geometry + boundary↔region index mapping."""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration.region_arrays import (
    RegionArrays,
    boundary_region_indices,
    region_boundary_indices,
)
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    TS_NEG,
    TS_NONE,
    TS_POS,
)


def _scrambled_region_df():
    # Two refs, rows deliberately out of (ref_id, start) order.
    rows = [
        ("chr2", 150, 400, BIT_EXON_NEG),
        ("chr1", 100, 200, BIT_EXON_POS),
        ("chr2", 0, 150, 0),
        ("chr1", 0, 100, 0),
        ("chr1", 200, 300, BIT_INTRON_POS),
    ]
    return pd.DataFrame(
        {
            "region_id": np.arange(len(rows), dtype=np.int64),
            "ref_name": pd.array([r[0] for r in rows], dtype="string"),
            "start": np.array([r[1] for r in rows], dtype=np.int64),
            "end": np.array([r[2] for r in rows], dtype=np.int64),
            "length": np.array([r[2] - r[1] for r in rows], dtype=np.int64),
            "signature": np.array([r[3] for r in rows], dtype=np.uint8),
        }
    )


def test_from_region_df_csr_ordering():
    df = _scrambled_region_df()
    ra = RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})

    # Sorted by (ref_id, start).
    np.testing.assert_array_equal(ra.ref_id, [0, 0, 0, 1, 1])
    np.testing.assert_array_equal(ra.start, [0, 100, 200, 0, 150])
    np.testing.assert_array_equal(ra.end, [100, 200, 300, 150, 400])
    np.testing.assert_array_equal(ra.ref_offsets, [0, 3, 5])
    np.testing.assert_array_equal(ra.region_size_bp, [100, 100, 100, 150, 250])
    np.testing.assert_array_equal(ra.strand_class, [TS_NONE, TS_POS, TS_POS, TS_NONE, TS_NEG])
    assert ra.n_regions == 5
    assert ra.n_refs == 2


def test_from_region_df_requires_signature():
    df = _scrambled_region_df().drop(columns=["signature"])
    try:
        RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})
    except ValueError as exc:
        assert "signature" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for missing signature column")


# Topology: ref0 = 3 regions (4 boundaries), ref1 = 2 regions (3 boundaries),
# ref2 = empty (0 boundaries). Invariant: B == R + nonempty_refs == 5 + 2 == 7.
REF_REGION_OFFSETS = np.array([0, 3, 5, 5], dtype=np.int64)
REF_BOUNDARY_OFFSETS = np.array([0, 4, 7, 7], dtype=np.int64)


def test_region_boundary_indices():
    left, right = region_boundary_indices(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)
    np.testing.assert_array_equal(left, [0, 1, 2, 4, 5])
    np.testing.assert_array_equal(right, [1, 2, 3, 5, 6])


def test_boundary_region_indices_terminals():
    left, right = boundary_region_indices(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)
    np.testing.assert_array_equal(left, [-1, 0, 1, 2, -1, 3, 4])
    np.testing.assert_array_equal(right, [0, 1, 2, -1, 3, 4, -1])


def test_boundary_region_mapping_round_trips():
    lb, rb = region_boundary_indices(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)
    lr, rr = boundary_region_indices(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)

    # Region → boundary → region is the identity (region is right of its left
    # boundary, left of its right boundary).
    r = np.arange(len(lb))
    np.testing.assert_array_equal(rr[lb], r)
    np.testing.assert_array_equal(lr[rb], r)

    # Internal boundary → region → boundary is the identity.
    internal = (lr >= 0) & (rr >= 0)
    np.testing.assert_array_equal(rb[lr[internal]], np.flatnonzero(internal))
    np.testing.assert_array_equal(lb[rr[internal]], np.flatnonzero(internal))
