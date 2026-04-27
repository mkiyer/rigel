"""Unit tests for SRD Pass 2 1-D mixture (`rigel.calibration._fl_mixture`)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration._fl_mixture import fit_fl_mixture


def _gaussian_hist(n_bins: int, mean: float, sigma: float, n: int, rng) -> np.ndarray:
    samples = rng.normal(mean, sigma, n).astype(int)
    samples = np.clip(samples, 0, n_bins - 1)
    return np.bincount(samples, minlength=n_bins).astype(np.float64)


def test_recovers_known_pi_and_gdna_distribution():
    rng = np.random.default_rng(0)
    n_bins = 1001
    rna_hist = _gaussian_hist(n_bins, 200, 30, 50_000, rng)
    rna_probs = rna_hist / rna_hist.sum()
    gdna_hist_truth = _gaussian_hist(n_bins, 350, 60, 50_000, rng)
    gdna_truth_probs = gdna_hist_truth / gdna_hist_truth.sum()

    pi_true = 0.3
    n_pool = 20_000
    n_g = int(n_pool * pi_true)
    n_r = n_pool - n_g
    pool = (
        _gaussian_hist(n_bins, 200, 30, n_r, rng)
        + _gaussian_hist(n_bins, 350, 60, n_g, rng)
    )

    out = fit_fl_mixture(pool, rna_probs, max_iter=200, tol=1e-6)
    assert out.converged
    # 1-D mixture EM has some bias when components share support tails.
    assert abs(out.pi - pi_true) < 0.05
    rec_probs = out.gdna_counts / max(out.gdna_counts.sum(), 1e-12)
    assert abs(int(rec_probs.argmax()) - 350) <= 20
    tv = 0.5 * float(np.abs(rec_probs - gdna_truth_probs).sum())
    assert tv < 0.25


def test_pi_zero_when_pool_is_pure_rna():
    rng = np.random.default_rng(1)
    n_bins = 401
    rna_hist = _gaussian_hist(n_bins, 200, 30, 30_000, rng)
    rna_probs = rna_hist / rna_hist.sum()
    pool = _gaussian_hist(n_bins, 200, 30, 5000, rng)

    out = fit_fl_mixture(pool, rna_probs, max_iter=200)
    assert out.pi < 0.05


def test_pi_one_when_pool_disjoint_from_rna():
    rng = np.random.default_rng(2)
    n_bins = 1001
    rna_hist = _gaussian_hist(n_bins, 200, 20, 30_000, rng)
    rna_probs = rna_hist / rna_hist.sum()
    pool = _gaussian_hist(n_bins, 600, 20, 5000, rng)

    out = fit_fl_mixture(pool, rna_probs, max_iter=200)
    assert out.pi > 0.95


def test_empty_pool_returns_zero_pi():
    out = fit_fl_mixture(np.zeros(101), np.full(101, 1.0 / 101))
    assert out.pi == 0.0
    assert out.n_pool == 0
    assert out.gdna_counts.sum() == 0


def test_identical_distributions_unidentifiable():
    # When gDNA_FL ≡ RNA_FL, π is unidentifiable. Fit should still
    # return without crashing; gdna_counts shape is correct.
    rng = np.random.default_rng(3)
    n_bins = 401
    rna_hist = _gaussian_hist(n_bins, 200, 30, 30_000, rng)
    rna_probs = rna_hist / rna_hist.sum()
    pool = _gaussian_hist(n_bins, 200, 30, 5000, rng)
    out = fit_fl_mixture(pool, rna_probs, max_iter=200)
    assert out.gdna_counts.shape == (n_bins,)
    assert 0.0 <= out.pi <= 1.0


def test_mismatched_lengths_raises():
    with pytest.raises(ValueError):
        fit_fl_mixture(np.ones(50), np.ones(40))


def test_gdna_counts_scale_to_pi_n_pool():
    rng = np.random.default_rng(4)
    n_bins = 1001
    rna_hist = _gaussian_hist(n_bins, 200, 30, 50_000, rng)
    rna_probs = rna_hist / rna_hist.sum()
    pi_true = 0.4
    n_pool = 10_000
    pool = (
        _gaussian_hist(n_bins, 200, 30, int(n_pool * (1 - pi_true)), rng)
        + _gaussian_hist(n_bins, 400, 50, int(n_pool * pi_true), rng)
    )
    out = fit_fl_mixture(pool, rna_probs, max_iter=200, tol=1e-6)
    # Total recovered gdna mass must equal π · pool_total.
    assert abs(out.gdna_counts.sum() - out.pi * pool.sum()) < 1.0
