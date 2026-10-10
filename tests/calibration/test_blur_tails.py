"""Finite Gaussian convolution retains relative likelihood in deep tails."""

import numpy as np
import pytest
from scipy.special import logsumexp

from rigel.native import transfer_rows as rows


def blur(row, grid, variance):
    row = np.asarray(row, dtype=float)
    if variance <= 0 or len(row) == 1:
        return row - row.max()
    step = grid[1] - grid[0]
    half = max(int(np.ceil(4 * np.sqrt(variance) / step)), 1)
    offsets = np.arange(-half, half + 1) * step
    log_kernel = -0.5 * offsets**2 / variance
    log_kernel -= logsumexp(log_kernel)
    padded = np.pad(row, half, mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, len(offsets))
    result = logsumexp(windows + log_kernel, axis=1)
    return result - result.max()


@pytest.mark.parametrize("slope", [1.0, 100.0, 1000.0])
@pytest.mark.parametrize("variance", [0.01, 0.5, 5.0])
def test_same_convolution_at_every_log_likelihood_scale(slope, variance):
    grid = np.arange(-50, 51) * 0.2
    row = -slope * (grid - 1.3) ** 2
    np.testing.assert_allclose(
        rows.blur_row(row, grid, variance), blur(row, grid, variance), atol=2e-9, rtol=2e-13
    )


@pytest.mark.parametrize("variance", [0.0, 0.1, 10.0])
def test_normalization_constant_and_coordinate_origin_do_not_change_shape(variance):
    grid = np.arange(-40, 41) * 0.2
    row = -150 * (grid - 0.8) ** 2
    expected = blur(row, grid, variance)
    for offset in (-10000.0, 0.0, 10000.0):
        observed = rows.blur_row(row + offset, grid + 32.0, variance)
        np.testing.assert_allclose(observed, expected, atol=2e-9, rtol=2e-13)


@pytest.mark.parametrize("variance", [1e-8, 1e-6])
def test_a_small_kernel_weight_can_be_relevant_to_a_deep_tail(variance):
    grid = np.arange(-3, 4) * 0.2
    row = np.array([-1e8, -1e8, -1e8, 0.0, -1e8, -1e8, -1e8])
    np.testing.assert_allclose(
        rows.blur_row(row, grid, variance), blur(row, grid, variance), atol=2e-8, rtol=2e-13
    )


@pytest.mark.parametrize("count", [100.0, 500.0])
def test_certified_source_keeps_its_low_density_tail(count):
    grid = np.arange(-125, 126) * 0.2
    rate, origin = 1.0, 1.0
    variance = rows.hop_price(count, count / rate, 0.0, 0.0)
    step = grid[1] - grid[0]
    half = max(int(np.ceil(4 * np.sqrt(variance) / step)), 1)
    extended = grid[0] + np.arange(-half, len(grid) + half) * step
    mean = count / rate * origin * np.exp(extended)
    source = count * np.log(mean) - mean
    expected = blur(source, extended, variance)[half : half + len(grid)]
    expected = np.maximum.accumulate(expected)
    expected -= expected.max()
    observed = rows.flux_level(grid, count, rate, origin, variance)
    np.testing.assert_allclose(observed, expected, atol=2e-9, rtol=2e-13)


def test_zero_width_preserves_the_input_exactly():
    grid = np.arange(-50, 51) * 0.2
    row = -1000 * (grid - 1.3) ** 2
    np.testing.assert_array_equal(rows.blur_row(row, grid, 0.0), row - row.max())


def test_zero_likelihood_outside_the_kernel_footprint_stays_zero():
    grid = np.arange(-3, 4) * 0.2
    row = np.full(7, -np.inf)
    row[3] = 0.0
    expected = blur(row, grid, 0.001)
    actual = rows.blur_row(row, grid, 0.001)
    np.testing.assert_array_equal(np.isneginf(actual), np.isneginf(expected))
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)
