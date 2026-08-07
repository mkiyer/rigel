"""Tests for the fragment length model factory + container serialization.

Covers:
- FragmentLengthModel.from_counts() factory for pre-built distributions

The library-wide gDNA / RNA FL distributions are owned by
``rigel.calibration.fl.FLModels`` (built from these raw counts); their tests
live in ``tests/calibration/test_fl.py``.
"""

import math

import numpy as np
import pytest

from rigel.frag_length_model import (
    FragmentLengthModel,
)


# =====================================================================
# FragmentLengthModel.from_counts() factory
# =====================================================================


class TestFromCounts:
    """Tests for the from_counts() factory method."""

    def test_basic_creation(self):
        """Create a model from a simple histogram."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[250] = 50.0
        model = FragmentLengthModel.from_counts(counts)

        assert model._finalized
        assert model.max_size == 500
        assert model.total_weight == pytest.approx(150.0)
        assert model.n_observations == 150
        assert model._log_prob is not None

    def test_max_size_inferred(self):
        """max_size is inferred as len(counts) - 1."""
        counts = np.ones(101)
        model = FragmentLengthModel.from_counts(counts)
        assert model.max_size == 100

    def test_explicit_max_size(self):
        """Explicit max_size with shorter counts array zero-fills."""
        counts = np.array([0, 0, 10, 20, 10])
        model = FragmentLengthModel.from_counts(counts, max_size=1000)
        assert model.max_size == 1000
        assert model.total_weight == pytest.approx(40.0)
        assert model.counts[2] == pytest.approx(10.0)
        assert model.counts[999] == pytest.approx(0.0)

    def test_log_likelihood_matches_trained(self):
        """from_counts() produces same log_likelihood as observe+finalize."""
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(200)
        for _ in range(50):
            trained.observe(300)
        trained.finalize()

        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 50.0
        factory = FragmentLengthModel.from_counts(counts)

        for length in [0, 100, 200, 250, 300, 500]:
            assert factory.log_likelihood(length) == pytest.approx(trained.log_likelihood(length))

    def test_tail_decay_works(self):
        """Queries beyond max_size use exponential tail decay."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        model = FragmentLengthModel.from_counts(counts)

        ll_500 = model.log_likelihood(500)
        ll_600 = model.log_likelihood(600)
        assert ll_600 < ll_500
        assert ll_600 == pytest.approx(model._tail_base + (600 - 500) * math.log(0.99))

    def test_empty_counts(self):
        """Zero-count histogram yields uniform distribution."""
        counts = np.zeros(101, dtype=np.float64)
        model = FragmentLengthModel.from_counts(counts)
        assert model._finalized
        ll_0 = model.log_likelihood(0)
        ll_50 = model.log_likelihood(50)
        assert ll_0 == pytest.approx(ll_50)
        assert ll_0 == pytest.approx(-np.log(101))

    def test_single_peak(self):
        """Single-peak histogram: mode has highest likelihood."""
        counts = np.zeros(501, dtype=np.float64)
        counts[250] = 1000.0
        model = FragmentLengthModel.from_counts(counts)
        assert model.log_likelihood(250) > model.log_likelihood(100)
        assert model.log_likelihood(250) > model.log_likelihood(400)
        assert model.mode == 250

    def test_statistics(self):
        """Finalized statistics use the same posterior predictive PMF as scoring."""
        counts = np.zeros(501, dtype=np.float64)
        counts[200] = 100.0
        counts[300] = 100.0
        model = FragmentLengthModel.from_counts(counts)
        # In-range mean of the finalized PMF: symmetric spikes at 200/300 → ~250
        # (excluding the empty overflow bin's smoothing mass is a <0.01 shift).
        assert model.mean == pytest.approx(250.0, abs=0.05)
        probs = np.exp(model._log_prob)
        expected_median = float(np.searchsorted(np.cumsum(probs), 0.5))
        assert model.median == pytest.approx(expected_median)
        assert model.mode in (200, 300)

    def test_from_counts_matches_normal_training(self):
        """Manually constructed histogram via from_counts matches
        the same histogram built through observe() calls."""
        trained = FragmentLengthModel(max_size=500)
        for _ in range(100):
            trained.observe(180)
        for _ in range(200):
            trained.observe(220)
        for _ in range(50):
            trained.observe(260)
        trained.finalize()

        counts = np.zeros(501, dtype=np.float64)
        counts[180] = 100
        counts[220] = 200
        counts[260] = 50
        factory = FragmentLengthModel.from_counts(counts)

        for length in range(501):
            assert factory.log_likelihood(length) == pytest.approx(
                trained.log_likelihood(length)
            ), f"Mismatch at length={length}"


# ⛔ `TestToDict` lived here and was DELETED by TRAPS: pure-and-length-censored with the container it tested
# it pinned the SCANNER's raw global + per-splice-category
# histograms, which are gone. The report's fragment-length categories now come from
# `FLModels` alone and are gated in tests/test_summary_report.py.
