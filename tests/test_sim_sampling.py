"""The fragment-length sampler always returns, and draws the window's law: rejection where the window holds enough of
the rounded normal (every panel's law), the window's exact integer masses where it does not, a point mass where the
law has no spread, and a refusal where the window can hold no fragment. A window rejection could never fill looped
forever."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.sim.sampling import truncated_normal_frag_lengths


def _draw(n, mean, std, lo, hi, seed=0):
    return truncated_normal_frag_lengths(np.random.default_rng(seed), n, mean, std, lo, hi)


def test_a_central_window_is_drawn_by_rejection_as_before():
    """A panel's law: the draws are rejection's, the same stream they always were."""
    rng = np.random.default_rng(7)
    raw = np.rint(rng.normal(206.0, 98.0, int(1000 * 1.5) + 10)).astype(int)
    kept = raw[(raw >= 50) & (raw <= 500)][:1000]
    assert np.array_equal(_draw(1000, 206.0, 98.0, 50, 500, seed=7)[: len(kept)], kept)


def test_a_law_with_no_spread_is_its_point():
    assert np.array_equal(_draw(5, 300.0, 0.0, 50, 400), np.full(5, 300))


def test_a_point_outside_the_window_is_refused():
    with pytest.raises(ValueError, match="no spread"):
        _draw(5, 300.0, 0.0, 50, 200)


def test_a_one_length_window_far_in_a_tail_is_that_length():
    """Twenty spreads above the mean: rejection would never keep one; the window's one integer is the law."""
    assert np.array_equal(_draw(4, 0.0, 10.0, 200, 200), np.full(4, 200))


def test_a_tail_window_keeps_its_relative_masses():
    """Ten to twelve spreads up: each length's share is the rounded normal's mass there, renormalised, so it falls
    with length at the normal's rate."""
    draws = _draw(200_000, 0.0, 10.0, 100, 103)
    counts = np.bincount(draws - 100, minlength=4).astype(float)
    k = np.arange(100, 104)
    want = np.exp(-0.5 * (k / 10.0) ** 2)
    np.testing.assert_allclose(counts / counts.sum(), want / want.sum(), atol=5e-3)


def test_no_fragment_requested_is_none_drawn():
    assert _draw(0, 0.0, 10.0, 200, 210).size == 0


@pytest.mark.parametrize("mean, std", [(float("nan"), 10.0), (100.0, float("inf")), (100.0, -1.0)])
def test_a_law_that_is_not_a_law_is_refused(mean, std):
    with pytest.raises(ValueError, match="finite"):
        _draw(5, mean, std, 50, 200)
