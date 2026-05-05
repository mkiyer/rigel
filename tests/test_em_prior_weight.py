"""Tests for the M5 ``prior_weight_rna`` ABI extension.

Verifies that:

* ``None`` and an explicit all-ones float32 array produce bit-identical
  EM output (regression guard).
* Per-component weights shift the posterior mass in the expected
  direction.
* The VBEM mode honours the weights identically.
* Degenerate (all-zero) RNA weights fall back to the uniform branch
  without producing NaNs.
"""

import numpy as np

from rigel._em_impl import run_locus_em_native


# ---------------------------------------------------------------------------
# Local helper — almost identical to ``tests/test_em_impl.py::_make_locus``
# but exposes the M5 ``prior_weight_rna`` argument.
# ---------------------------------------------------------------------------


def _run_em(
    *,
    units: list[list[tuple[int, float, float, int, int]]],
    n_components: int,
    alpha_rna: float = 4.0,
    use_vbem: bool = False,
    prior_weight_rna: np.ndarray | None = None,
):
    all_t_indices: list[int] = []
    all_log_liks: list[float] = []
    all_cov_wts: list[float] = []
    all_tx_starts: list[int] = []
    all_tx_ends: list[int] = []
    offsets = [0]
    for unit in units:
        for comp_idx, ll, cw, txs, txe in unit:
            all_t_indices.append(comp_idx)
            all_log_liks.append(ll)
            all_cov_wts.append(cw)
            all_tx_starts.append(txs)
            all_tx_ends.append(txe)
        offsets.append(len(all_t_indices))

    bias_profiles = np.full(n_components, 100_000, dtype=np.int64)
    unambig_totals = np.zeros(n_components, dtype=np.float64)
    prior_eligible = np.ones(n_components, dtype=np.float64)
    eff_lens = np.ones(n_components, dtype=np.float64)

    theta, alpha, em_totals = run_locus_em_native(
        np.asarray(offsets, dtype=np.int64),
        np.asarray(all_t_indices, dtype=np.int32),
        np.asarray(all_log_liks, dtype=np.float64),
        np.asarray(all_cov_wts, dtype=np.float64),
        np.asarray(all_tx_starts, dtype=np.int32),
        np.asarray(all_tx_ends, dtype=np.int32),
        bias_profiles,
        unambig_totals,
        eff_lens,
        prior_eligible,
        n_components,
        alpha_rna,
        1000,           # max_iterations
        1e-9,           # convergence_delta
        use_vbem,
        0,              # n_transcripts (classic/unlinked)
        prior_weight_rna,
    )
    return np.asarray(theta), np.asarray(alpha), np.asarray(em_totals)


def _ambiguous_three_component_units(n_units: int = 200):
    """Build n_units identical 3-way ambiguous units.

    Coverage_wts are equal across components, so the OVR prior alone
    determines the prior asymmetry — making the weight effect easy to
    isolate.
    """
    return [
        [(0, 0.0, 1.0, 0, 1), (1, 0.0, 1.0, 0, 1), (2, 0.0, 1.0, 0, 1)]
        for _ in range(n_units)
    ]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_none_argument_bit_identical_to_baseline():
    """Passing ``prior_weight_rna=None`` matches omission of the kwarg."""
    units = _ambiguous_three_component_units()
    theta_none, alpha_none, em_none = _run_em(
        units=units, n_components=3, prior_weight_rna=None
    )
    theta_kwarg, alpha_kwarg, em_kwarg = _run_em(
        units=units, n_components=3
    )
    np.testing.assert_array_equal(theta_none, theta_kwarg)
    np.testing.assert_array_equal(alpha_none, alpha_kwarg)
    np.testing.assert_array_equal(em_none, em_kwarg)


def test_all_ones_bit_identical_to_none():
    """Explicit all-ones weights produce bit-identical EM output."""
    units = _ambiguous_three_component_units()
    theta_none, alpha_none, em_none = _run_em(
        units=units, n_components=3, prior_weight_rna=None
    )
    ones = np.ones(3, dtype=np.float32)
    theta_ones, alpha_ones, em_ones = _run_em(
        units=units, n_components=3, prior_weight_rna=ones
    )
    np.testing.assert_array_equal(theta_none, theta_ones)
    np.testing.assert_array_equal(alpha_none, alpha_ones)
    np.testing.assert_array_equal(em_none, em_ones)


def test_downweighted_component_gets_smaller_share():
    """Halving one component's weight halves its share of the prior."""
    units = _ambiguous_three_component_units()
    weights = np.array([1.0, 0.5, 1.0], dtype=np.float32)
    theta_w, alpha_w, _ = _run_em(
        units=units, n_components=3, prior_weight_rna=weights
    )
    # Symmetric components 0 and 2 should be equal; component 1 strictly less.
    assert theta_w[0] > theta_w[1]
    assert theta_w[2] > theta_w[1]
    np.testing.assert_allclose(theta_w[0], theta_w[2], rtol=1e-9)
    # The posterior split between {0,2} and {1} should reflect the prior asymmetry.
    # Expected ratio of theta[1] / theta[0] is bounded above by the prior
    # weight ratio (0.5) — equivalence-class evidence is symmetric, so
    # the prior dominates.
    assert theta_w[1] / theta_w[0] < 0.7


def test_vbem_baseline_preserved_with_ones():
    """VBEM mode honours the weights identically (all-ones \u2261 None)."""
    units = _ambiguous_three_component_units()
    theta_none, _, _ = _run_em(
        units=units, n_components=3, use_vbem=True, prior_weight_rna=None
    )
    theta_ones, _, _ = _run_em(
        units=units,
        n_components=3,
        use_vbem=True,
        prior_weight_rna=np.ones(3, dtype=np.float32),
    )
    np.testing.assert_array_equal(theta_none, theta_ones)


def test_degenerate_all_zero_weights_fall_back_to_uniform():
    """All-zero RNA weights trigger the uniform fallback branch (no NaN)."""
    units = _ambiguous_three_component_units(n_units=50)
    zeros = np.zeros(3, dtype=np.float32)
    theta_zero, alpha_zero, em_zero = _run_em(
        units=units, n_components=3, prior_weight_rna=zeros
    )
    assert np.all(np.isfinite(theta_zero))
    assert np.all(np.isfinite(alpha_zero))
    assert np.all(np.isfinite(em_zero))
    # Uniform fallback: all three components receive equal prior \u21d2
    # symmetric posterior on this fully-symmetric input.
    np.testing.assert_allclose(theta_zero[0], theta_zero[1], rtol=1e-6)
    np.testing.assert_allclose(theta_zero[1], theta_zero[2], rtol=1e-6)
