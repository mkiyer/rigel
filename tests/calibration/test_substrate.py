"""CalibrationSubstrate channel reductions + D1 boundary→region projection."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _synthetic import make_synthetic_payload

from rigel.calibration.errors import CalibrationSubstrateError
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_NEG, TS_NONE, TS_POS
from rigel.calibration.substrate import CalibrationSubstrate


def test_contained_reductions():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    assert sub.n_regions == 3
    np.testing.assert_array_equal(sub.ts_class, [TS_POS, TS_NEG, TS_NONE])
    np.testing.assert_array_equal(sub.L_eff, [100.0, 100.0, 100.0])

    # n_u = ch0+ch1, n_s = ch2+ch3.
    np.testing.assert_array_equal(sub.contained.n_unspliced, [12, 21, 15])
    np.testing.assert_array_equal(sub.contained.n_spliced, [3, 5, 0])
    # k_plus = sense among unspliced: + → ch0, − → ch1, NONE → ch0 (fixed).
    np.testing.assert_array_equal(sub.contained.k_plus, [10, 20, 7])
    # contained mass == contained count.
    np.testing.assert_array_equal(sub.contained.mass_unspliced, [12.0, 21.0, 15.0])
    np.testing.assert_array_equal(sub.contained.mass_spliced, [3.0, 5.0, 0.0])


def test_boundary_projection_left_view():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # left view of region r ← its LEFT boundary's RIGHT side.
    #   r0: left boundary b0 (terminal) → all zero
    #   r1: b1 flux [4,1,0,0], mass_right [2.5,0.5,0,0]; r1 is − → k_plus = ch1 = 1
    #   r2: b2 flux [2,3,1,0], mass_right [0.5,1.0,0,0]; r2 NONE → k_plus = ch0 = 2
    np.testing.assert_array_equal(sub.left.n_unspliced, [0, 5, 5])
    np.testing.assert_array_equal(sub.left.n_spliced, [0, 0, 1])
    np.testing.assert_array_equal(sub.left.k_plus, [0, 1, 2])
    np.testing.assert_allclose(sub.left.mass_unspliced, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(sub.left.mass_spliced, [0.0, 0.0, 0.0])


def test_boundary_projection_right_view():
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # right view of region r ← its RIGHT boundary's LEFT side.
    #   r0: b1 flux [4,1,0,0], mass_left [1.5,0.5,0,0]; r0 is + → k_plus = ch0 = 4
    #   r1: b2 flux [2,3,1,0], mass_left [3.0,1.0,0.5,0]; r1 is − → k_plus = ch1 = 3
    #   r2: right boundary b3 (terminal) → all zero
    np.testing.assert_array_equal(sub.right.n_unspliced, [5, 5, 0])
    np.testing.assert_array_equal(sub.right.n_spliced, [0, 1, 0])
    np.testing.assert_array_equal(sub.right.k_plus, [4, 3, 0])
    np.testing.assert_allclose(sub.right.mass_unspliced, [2.0, 4.0, 0.0])
    np.testing.assert_allclose(sub.right.mass_spliced, [0.0, 0.5, 0.0])


def test_shared_flux_used_by_both_neighbours():
    # b1's integer flux is consumed by r0's RIGHT view and r1's LEFT view.
    payload, ra = make_synthetic_payload()
    sub = CalibrationSubstrate.from_payload(payload, ra)
    assert sub.right.n_unspliced[0] == sub.left.n_unspliced[1] == 5


def test_misaligned_payload_raises():
    payload, _ = make_synthetic_payload()
    # A 2-region geometry against a 3-region payload must be rejected.
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
