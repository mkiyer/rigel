"""CalibrationSubstrate channel reductions + D1 per-side boundary projection."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _synthetic import make_synthetic_payload

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.scan_payload import AccumulatorPayload


def test_contained_reductions():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    assert sub.n_regions == 3
    np.testing.assert_array_equal(sub.ts_class, [TS_POS, TS_NEG, TS_NONE])
    np.testing.assert_array_equal(sub.region_len, [100.0, 100.0, 100.0])

    # Unspliced channels are raw genome strand (pos = ch0, neg = ch1).
    np.testing.assert_array_equal(sub.contained.n_unspliced_pos, [10, 1, 7])
    np.testing.assert_array_equal(sub.contained.n_unspliced_neg, [2, 20, 8])
    # Spliced channels are motif-relative sense (ch2) / antisense (ch3).
    np.testing.assert_array_equal(sub.contained.n_spliced_sense, [3, 0, 0])
    np.testing.assert_array_equal(sub.contained.n_spliced_antisense, [0, 5, 0])
    # Totals (properties) + mass (== count for contained).
    np.testing.assert_array_equal(sub.contained.n_unspliced, [12, 21, 15])
    np.testing.assert_array_equal(sub.contained.n_spliced, [3, 5, 0])
    np.testing.assert_array_equal(sub.contained.mass_unspliced, [12.0, 21.0, 15.0])
    np.testing.assert_array_equal(sub.contained.mass_spliced, [3.0, 5.0, 0.0])


def test_boundary_projection_left_view():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # left view of region r ← the RIGHT side of r's LEFT boundary
    # (flux_right / mass_right at left_boundary(r)).
    #   r0: left boundary b0 (terminal) → zero
    #   r1: b1 flux_right [6,2,0,0], mass_right [2.5,0.5,0,0]
    #   r2: b2 flux_right [5,4,0,1], mass_right [0.5,1.0,0,0]
    np.testing.assert_array_equal(sub.left.n_unspliced_pos, [0, 6, 5])
    np.testing.assert_array_equal(sub.left.n_unspliced_neg, [0, 2, 4])
    np.testing.assert_array_equal(sub.left.n_spliced_sense, [0, 0, 0])
    np.testing.assert_array_equal(sub.left.n_spliced_antisense, [0, 0, 1])
    np.testing.assert_allclose(sub.left.mass_unspliced, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(sub.left.mass_spliced, [0.0, 0.0, 0.0])


def test_boundary_projection_right_view():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # right view of region r ← the LEFT side of r's RIGHT boundary
    # (flux_left / mass_left at right_boundary(r)).
    #   r0: b1 flux_left [4,1,0,0], mass_left [1.5,0.5,0,0]
    #   r1: b2 flux_left [2,3,1,0], mass_left [3.0,1.0,0.5,0]
    #   r2: right boundary b3 (terminal) → zero
    np.testing.assert_array_equal(sub.right.n_unspliced_pos, [4, 2, 0])
    np.testing.assert_array_equal(sub.right.n_unspliced_neg, [1, 3, 0])
    np.testing.assert_array_equal(sub.right.n_spliced_sense, [0, 1, 0])
    np.testing.assert_array_equal(sub.right.n_spliced_antisense, [0, 0, 0])
    np.testing.assert_allclose(sub.right.mass_unspliced, [2.0, 4.0, 0.0])
    np.testing.assert_allclose(sub.right.mass_spliced, [0.0, 0.5, 0.0])


def test_flux_is_per_side_not_shared():
    # b1 is shared by r0 (its right boundary) and r1 (its left boundary), but
    # each neighbour sees a DIFFERENT side: r0's right view = flux_left[b1] (4),
    # r1's left view = flux_right[b1] (6). Per-side, not a shared count.
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    assert sub.right.n_unspliced_pos[0] == 4
    assert sub.left.n_unspliced_pos[1] == 6
    assert sub.right.n_unspliced_pos[0] != sub.left.n_unspliced_pos[1]


def test_ambiguous_region_is_flagged_not_silently_oriented():
    # A region with transcripts on BOTH strands (AMBIG) must be flagged via
    # ts_class (distinct from TS_NONE) so downstream excludes its unspliced
    # strand. Strand-agnostic counts stay valid.
    region_contained = np.array([[5, 8, 0, 0]], dtype=np.uint32)
    payload = AccumulatorPayload(
        boundaries=np.array([0, 300], dtype=np.int64),
        ref_pos_offsets=np.array([0, 2], dtype=np.int64),
        ref_region_offsets=np.array([0, 1], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 2], dtype=np.int64),
        region_contained=region_contained,
        boundary_mass_left=np.zeros((2, 4), dtype=np.float32),
        boundary_mass_right=np.zeros((2, 4), dtype=np.float32),
        boundary_flux_left=np.zeros((2, 4), dtype=np.uint32),
        boundary_flux_right=np.zeros((2, 4), dtype=np.uint32),
        n_refs=1,
    )
    region_df = pd.DataFrame(
        {
            "region_id": np.array([0], dtype=np.int64),
            "ref_name": pd.array(["chr1"], dtype="string"),
            "start": np.array([0], dtype=np.int64),
            "end": np.array([300], dtype=np.int64),
            "length": np.array([300], dtype=np.int64),
            "signature": np.array([BIT_EXON_POS | BIT_EXON_NEG], dtype=np.uint8),
        }
    )
    ra = RegionArrays.from_region_df(region_df, {"chr1": 0})
    sub = CalibrationSubstrate.from_payload(payload, ra)

    assert sub.ts_class[0] == TS_AMBIG
    assert sub.ts_class[0] != TS_NONE
    assert sub.contained.n_unspliced[0] == 13  # strand-agnostic density stays valid


def test_misaligned_payload_raises():
    payload, _ = make_synthetic_payload()
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(2, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * 2, dtype="string"),
            "start": np.array([0, 100], dtype=np.int64),
            "end": np.array([100, 200], dtype=np.int64),
            "length": np.array([100, 100], dtype=np.int64),
            "signature": np.array([2, 0], dtype=np.uint8),
        }
    )
    bad = RegionArrays.from_region_df(region_df, {"chr1": 0})
    with pytest.raises(CalibrationSubstrateError):
        CalibrationSubstrate.from_payload(payload, bad)
