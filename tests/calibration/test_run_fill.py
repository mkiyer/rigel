"""Shared region-geometry helpers (`calibration.run_fill`): same-ref masks + bidirectional run-fill."""

from __future__ import annotations

import numpy as np

from rigel.calibration.run_fill import runfill_bidirectional, same_ref_left_right


def test_same_ref_left_right_edges_and_breaks():
    # Two references: [0,0,0 | 1,1]. Reference edges and cross-ref seams have no same-ref neighbour.
    ref = np.array([0, 0, 0, 1, 1])
    left, right = same_ref_left_right(ref)
    # left_same[i]: region i has a same-ref LEFT neighbour
    np.testing.assert_array_equal(left, [False, True, True, False, True])
    # right_same[i]: region i has a same-ref RIGHT neighbour
    np.testing.assert_array_equal(right, [True, True, False, True, False])


def test_same_ref_left_right_singletons():
    # A length-1 array (both edges) and a per-region singleton run.
    left, right = same_ref_left_right(np.array([7]))
    assert not left[0] and not right[0]
    left, right = same_ref_left_right(np.array([0, 1, 2]))  # every region its own reference
    np.testing.assert_array_equal(left, [False, False, False])
    np.testing.assert_array_equal(right, [False, False, False])


def test_runfill_keeps_set_values():
    ref = np.zeros(4, dtype=int)
    v = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(runfill_bidirectional(v, ref), v)  # nothing to fill


def test_runfill_one_sided_and_mean():
    # Anchors at the ends of one run; interior averages the two carried values; a single trailing
    # anchor carries one-sided.
    ref = np.zeros(5, dtype=int)
    v = np.array([2.0, np.nan, np.nan, np.nan, 6.0])
    out = runfill_bidirectional(v, ref)
    # fwd = [2,2,2,2,6], rev = [2,6,6,6,6] → interior means: (2+6)/2 = 4 each; ends unchanged.
    np.testing.assert_allclose(out, [2.0, 4.0, 4.0, 4.0, 6.0])


def test_runfill_does_not_cross_reference():
    # ref break between index 1 and 2: the right run has no anchor → stays nan (caller's global fallback).
    ref = np.array([0, 0, 1, 1])
    v = np.array([5.0, np.nan, np.nan, np.nan])
    out = runfill_bidirectional(v, ref)
    assert out[1] == 5.0  # filled within ref 0
    assert np.isnan(out[2]) and np.isnan(out[3])  # ref 1 has no anchor — not carried across the seam


def test_runfill_unreachable_run_stays_nan():
    ref = np.zeros(3, dtype=int)
    out = runfill_bidirectional(np.array([np.nan, np.nan, np.nan]), ref)
    assert np.all(np.isnan(out))
