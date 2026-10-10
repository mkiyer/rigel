"""The source blur samples its known likelihood beyond the output table."""

import numpy as np
import pytest
from scipy.special import logsumexp

from rigel.native import transfer_rows as rows


def poisson_blur(log_density, count, opportunity, variance, step):
    """Independent discrete convolution, evaluating the source outside any table.

    The four-sd kernel is the existing operator, held fixed in this audit.
    This is not a proposed truncation constant or a different blur law.
    """
    half = max(int(np.ceil(4 * np.sqrt(variance) / step)), 1)
    offsets = np.arange(-half, half + 1) * step
    log_kernel = -0.5 * offsets**2 / variance
    log_kernel -= logsumexp(log_kernel)
    log_mean = np.asarray(log_density)[:, None] + offsets + np.log(opportunity)
    log_likelihood = count * log_mean - np.exp(log_mean)
    return logsumexp(log_likelihood + log_kernel, axis=1)


@pytest.mark.parametrize(
    "count,rate,origin,window",
    [
        (1.0, 0.001, 1.0, 10.0),
        (5.0, 100.0, 0.01, 10.0),
        (1.0, 0.01, 1.0, 5.0),
        (40.0, 0.02, 1.0, 10.0),
    ],
)
def test_source_blur_matches_direct_poisson_evaluation(count, rate, origin, window):
    grid = np.arange(-round(window / 0.2), round(window / 0.2) + 1) * 0.2
    variance = rows.hop_price(count, count / rate, 0.0, 0.0)
    reference = poisson_blur(
        grid + np.log(origin), count, count / rate, variance, float(grid[1] - grid[0])
    )
    reference = np.maximum.accumulate(reference)
    reference -= reference.max()
    observed = rows.flux_level(grid, count, rate, origin, variance)
    np.testing.assert_allclose(observed, reference, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("shift", [-14.0, 14.0])
def test_source_blur_preserves_physical_coordinates(shift):
    grid = np.linspace(-10, 10, 101)
    variance = rows.hop_price(1.0, 1000.0, 0.0, 0.0)
    a = rows.flux_level(grid, 1.0, 0.001, 1.0, variance)
    b = rows.flux_level(grid - shift, 1.0, 0.001, np.exp(shift), variance)
    np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-10)


def test_zero_width_preserves_the_existing_unblurred_path_exactly():
    grid = np.linspace(-10, 10, 101)
    reference = rows.lower_side(rows.poisson_level(grid, 3.0, 300.0, 1.0))
    observed = rows.flux_level(grid, 3.0, 0.01, 1.0, 0.0)
    np.testing.assert_array_equal(observed, reference)


def test_source_extension_preserves_the_original_convolution_spacing():
    grid = np.arange(-50, 51) * 0.03
    step = float(grid[1] - grid[0])
    # Put the kernel radius at an integer: reconstructing spacing from the
    # extended axis must not round it onto a different finite operator.
    variance = (31 * step / 4) ** 2
    reference = poisson_blur(grid, 1.0, 1.0, variance, step)
    reference = np.maximum.accumulate(reference)
    reference -= reference.max()
    observed = rows.flux_level(grid, 1.0, 1.0, 1.0, variance)
    np.testing.assert_allclose(observed, reference, rtol=1e-12, atol=1e-12)
