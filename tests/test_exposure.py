"""Tests for :mod:`rigel.calibration._exposure`."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._exposure import (
    boundary_crossing_exposure,
    boundary_side_in_window,
    contained_exposure_clipped,
)
from rigel.frag_length_model import FragmentLengthModel


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _fl_uniform(
    lo: int, hi: int, max_size: int = 1000, weight_per_bin: float = 1e6
) -> FragmentLengthModel:
    """Finalized FL model uniform over fragment lengths ``[lo, hi]`` (inclusive).

    Uses a large ``weight_per_bin`` so Laplace smoothing in
    ``_normalized_probs`` (counts+1)/(total+K) is numerically negligible —
    this keeps tests sensitive to the actual model semantics rather than to
    smoothing artifacts.
    """
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[lo : hi + 1] = weight_per_bin
    fl = FragmentLengthModel.from_counts(counts, max_size=max_size)
    fl.finalize()
    return fl


def _fl_delta(ell: int, max_size: int = 1000, weight: float = 1e9) -> FragmentLengthModel:
    """Finalized FL model with all mass at exactly ``ell``."""
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[ell] = weight
    fl = FragmentLengthModel.from_counts(counts, max_size=max_size)
    fl.finalize()
    return fl


# ---------------------------------------------------------------------------
# contained_exposure_clipped
# ---------------------------------------------------------------------------

class TestContainedExposureClipped:
    def test_clip_equals_full_returns_identical_arrays(self):
        fl = _fl_uniform(50, 200)
        starts = np.array([100, 500, 1000], dtype=np.int64)
        ends = np.array([400, 800, 2000], dtype=np.int64)
        eff_full, eff_clip = contained_exposure_clipped(
            starts, ends, clip_lo=0, clip_hi=10_000, fl=fl
        )
        np.testing.assert_array_equal(eff_full, eff_clip)
        # All spans positive, FL valid → all eff values positive.
        assert np.all(eff_full > 0.0)

    def test_clip_strictly_inside_full_yields_smaller_eff(self):
        fl = _fl_uniform(50, 200)
        starts = np.array([0], dtype=np.int64)
        ends = np.array([1000], dtype=np.int64)
        eff_full, eff_clip = contained_exposure_clipped(
            starts, ends, clip_lo=200, clip_hi=600, fl=fl
        )
        assert eff_clip[0] < eff_full[0]
        assert eff_clip[0] > 0.0

    def test_empty_clip_window_yields_zero_clipped_exposure(self):
        fl = _fl_uniform(50, 200)
        starts = np.array([100, 500], dtype=np.int64)
        ends = np.array([400, 800], dtype=np.int64)
        # clip_lo > clip_hi for one region, disjoint window for both
        eff_full, eff_clip = contained_exposure_clipped(
            starts, ends, clip_lo=900, clip_hi=950, fl=fl
        )
        assert np.all(eff_full > 0.0)
        np.testing.assert_array_equal(eff_clip, np.zeros(2, dtype=np.float64))

    def test_zero_span_region_has_zero_exposure(self):
        fl = _fl_uniform(50, 200)
        starts = np.array([100, 200], dtype=np.int64)
        ends = np.array([100, 250], dtype=np.int64)  # first is empty
        eff_full, eff_clip = contained_exposure_clipped(
            starts, ends, clip_lo=0, clip_hi=10_000, fl=fl
        )
        assert eff_full[0] == 0.0
        assert eff_clip[0] == 0.0
        assert eff_full[1] > 0.0

    def test_ratio_is_monotonic_in_clip_width(self):
        fl = _fl_uniform(50, 200)
        starts = np.array([0], dtype=np.int64)
        ends = np.array([1000], dtype=np.int64)
        widths = [100, 300, 600, 1000]
        ratios = []
        for w in widths:
            _, eff_clip = contained_exposure_clipped(
                starts, ends, clip_lo=0, clip_hi=w, fl=fl
            )
            ratios.append(eff_clip[0])
        # eff_clip monotonically increases as the clip window grows.
        for a, b in zip(ratios, ratios[1:]):
            assert a <= b

    def test_shape_mismatch_raises(self):
        fl = _fl_uniform(50, 200)
        with pytest.raises(ValueError, match="shape"):
            contained_exposure_clipped(
                np.array([0, 100]), np.array([200]), 0, 1000, fl
            )

    def test_empty_input_returns_empty(self):
        fl = _fl_uniform(50, 200)
        eff_full, eff_clip = contained_exposure_clipped(
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            0, 1000, fl,
        )
        assert eff_full.shape == (0,)
        assert eff_clip.shape == (0,)


# ---------------------------------------------------------------------------
# boundary_crossing_exposure
# ---------------------------------------------------------------------------

class TestBoundaryCrossingExposure:
    def test_delta_at_ell_yields_ell_minus_one(self):
        # Pure delta PMF at ℓ=100: B_cross = ℓ - 1 = 99.
        fl = _fl_delta(100, max_size=1000)
        b = boundary_crossing_exposure(fl)
        assert b == pytest.approx(99.0, rel=1e-5)

    def test_uniform_pmf_equals_mean_minus_one_for_min_ell_ge_1(self):
        fl = _fl_uniform(50, 200)
        b = boundary_crossing_exposure(fl)
        # All mass at ℓ ≥ 50, so max(ℓ-1, 0) == ℓ - 1 for the entire support.
        # Hence B_cross == E[ℓ] - 1. Tolerance accounts for the small mass
        # that Laplace smoothing leaks to bins outside [50, 200].
        assert b == pytest.approx(fl.mean - 1.0, rel=1e-3)

    def test_clamps_at_zero_when_pmf_concentrated_at_short_lengths(self):
        # PMF over ℓ ∈ {0, 1}: max(ℓ - 1, 0) == 0 everywhere → B_cross == 0.
        counts = np.zeros(11, dtype=np.float64)
        counts[0] = 5.0
        counts[1] = 5.0
        fl = FragmentLengthModel.from_counts(counts, max_size=10)
        fl.finalize()
        # Laplace smoothing puts a tiny bit of mass at ℓ ≥ 2, so the raw sum
        # is positive but small. The contract is "≥ 0"; here it should be
        # strictly positive due to smoothing tail.
        b = boundary_crossing_exposure(fl)
        assert b >= 0.0

    def test_post_finalize_uses_smoothed_pmf(self):
        # finalize() with a Dirichlet prior should change pmf and thus B_cross.
        counts = np.zeros(101, dtype=np.float64)
        counts[50] = 1.0
        fl_unfin = FragmentLengthModel.from_counts(counts, max_size=100)
        b_unfin = boundary_crossing_exposure(fl_unfin)
        fl_fin = FragmentLengthModel.from_counts(counts, max_size=100)
        prior = np.ones(101, dtype=np.float64)
        fl_fin.finalize(prior_counts=prior, prior_ess=10.0)
        b_fin = boundary_crossing_exposure(fl_fin)
        # Smoothing pulls mass toward the uniform prior → different value.
        assert b_unfin != b_fin


# ---------------------------------------------------------------------------
# boundary_side_in_window
# ---------------------------------------------------------------------------

class TestBoundarySideInWindow:
    def test_basic_inclusion(self):
        starts = np.array([100, 200, 300, 400], dtype=np.int64)
        ends = np.array([150, 250, 350, 450], dtype=np.int64)
        left_in, right_in = boundary_side_in_window(
            starts, ends, clip_lo=200, clip_hi=350
        )
        np.testing.assert_array_equal(left_in, [False, True, True, False])
        np.testing.assert_array_equal(right_in, [False, True, True, False])

    def test_inclusive_endpoints(self):
        starts = np.array([100, 200], dtype=np.int64)
        ends = np.array([150, 350], dtype=np.int64)
        left_in, right_in = boundary_side_in_window(
            starts, ends, clip_lo=100, clip_hi=350
        )
        # 100 sits at clip_lo (inclusive); 350 sits at clip_hi (inclusive).
        np.testing.assert_array_equal(left_in, [True, True])
        np.testing.assert_array_equal(right_in, [True, True])

    def test_left_and_right_can_disagree(self):
        # Region [100, 500): start outside window, end inside.
        starts = np.array([100], dtype=np.int64)
        ends = np.array([500], dtype=np.int64)
        left_in, right_in = boundary_side_in_window(
            starts, ends, clip_lo=400, clip_hi=600
        )
        np.testing.assert_array_equal(left_in, [False])
        np.testing.assert_array_equal(right_in, [True])

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="shape"):
            boundary_side_in_window(
                np.array([0, 100]), np.array([50]), 0, 1000
            )

    def test_empty_input_returns_empty(self):
        left_in, right_in = boundary_side_in_window(
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            0, 1000,
        )
        assert left_in.shape == (0,)
        assert right_in.shape == (0,)
        assert left_in.dtype == np.bool_
        assert right_in.dtype == np.bool_


# ---------------------------------------------------------------------------
# Cross-cutting: pmf accessor on FragmentLengthModel
# ---------------------------------------------------------------------------

class TestFragmentLengthModelPmf:
    def test_pmf_sums_to_one(self):
        fl = _fl_uniform(50, 200)
        pmf = fl.pmf
        assert pmf.shape == (1001,)
        assert pmf.sum() == pytest.approx(1.0, abs=1e-12)

    def test_pmf_matches_internal_normalized_probs(self):
        fl = _fl_uniform(50, 200)
        np.testing.assert_array_equal(fl.pmf, fl._normalized_probs())
