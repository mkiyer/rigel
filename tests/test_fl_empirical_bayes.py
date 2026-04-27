"""Unit tests for SRD Pass 3 EB FL builders (`rigel.calibration._fl_empirical_bayes`)."""

from __future__ import annotations

import numpy as np

from rigel.calibration._fl_empirical_bayes import (
    build_gdna_fl,
    build_global_fl,
    build_rna_fl,
)


def _gaussian_counts(n_bins: int, mean: float, sigma: float, n: int, rng) -> np.ndarray:
    samples = np.clip(rng.normal(mean, sigma, n).astype(int), 0, n_bins - 1)
    return np.bincount(samples, minlength=n_bins).astype(np.float64)


def test_global_fl_finalize_matches_normalized_counts():
    rng = np.random.default_rng(0)
    counts = _gaussian_counts(501, 200, 30, 10_000, rng)
    fl = build_global_fl(counts, max_size=500)
    # finalize sets log_pmf so probs sum to ≈ 1.
    probs = np.exp(fl._log_prob)
    assert abs(probs.sum() - 1.0) < 1e-6


def _shape_matches(a, b) -> None:
    """Assert two log-prob vectors share their mode and are highly correlated
    in probability space.

    The plain-Laplace global FL and a Dirichlet-shrunk fallback FL are both
    shaped by ``global_counts`` but use different smoothing constants, so
    log-probs in empty tail bins differ by a constant — degrading log-domain
    correlation. In probability space (where these models are actually used
    by the scorer), the agreement is essentially perfect.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    assert int(a.argmax()) == int(b.argmax())
    pa = np.exp(a)
    pb = np.exp(b)
    r = float(np.corrcoef(pa, pb)[0, 1])
    assert r > 0.999, f"probability-space correlation {r} below 0.999"


def test_rna_fl_collapses_to_global_shape_when_spliced_empty():
    rng = np.random.default_rng(1)
    global_counts = _gaussian_counts(501, 250, 50, 50_000, rng)
    spliced_counts = np.zeros(501, dtype=np.float64)

    rna_fl = build_rna_fl(spliced_counts, global_counts, max_size=500, prior_ess=500.0)
    global_fl = build_global_fl(global_counts, max_size=500)
    _shape_matches(rna_fl._log_prob, global_fl._log_prob)


def test_gdna_fl_falls_back_to_global_shape_when_counts_none():
    rng = np.random.default_rng(2)
    global_counts = _gaussian_counts(501, 250, 50, 50_000, rng)
    gdna_fl = build_gdna_fl(None, global_counts, max_size=500, prior_ess=500.0)
    global_fl = build_global_fl(global_counts, max_size=500)
    _shape_matches(gdna_fl._log_prob, global_fl._log_prob)


def test_gdna_fl_distinct_from_global_when_signal_strong():
    rng = np.random.default_rng(3)
    n_bins = 501
    # Global FL centred at 200; recovered gDNA counts centred at 400 with
    # large mass — EB shrinkage should not drown it out.
    global_counts = _gaussian_counts(n_bins, 200, 30, 50_000, rng)
    gdna_counts = _gaussian_counts(n_bins, 400, 40, 5_000, rng)

    gdna_fl = build_gdna_fl(gdna_counts, global_counts, max_size=500, prior_ess=500.0)
    global_fl = build_global_fl(global_counts, max_size=500)

    # gDNA_FL mode should sit closer to 400 than to 200.
    gdna_mode = int(np.argmax(gdna_fl._log_prob))
    assert abs(gdna_mode - 400) < abs(gdna_mode - 200)
    # And the two distributions should be measurably different.
    tv = 0.5 * float(np.abs(np.exp(gdna_fl._log_prob) - np.exp(global_fl._log_prob)).sum())
    assert tv > 0.1


def test_higher_ess_pulls_rna_closer_to_global():
    rng = np.random.default_rng(4)
    n_bins = 501
    global_counts = _gaussian_counts(n_bins, 200, 30, 50_000, rng)
    spliced_counts = _gaussian_counts(n_bins, 400, 40, 1_000, rng)

    rna_low = build_rna_fl(spliced_counts, global_counts, max_size=500, prior_ess=10.0)
    rna_high = build_rna_fl(spliced_counts, global_counts, max_size=500, prior_ess=10_000.0)
    global_fl = build_global_fl(global_counts, max_size=500)

    glo = np.exp(global_fl._log_prob)
    tv_low = 0.5 * float(np.abs(np.exp(rna_low._log_prob) - glo).sum())
    tv_high = 0.5 * float(np.abs(np.exp(rna_high._log_prob) - glo).sum())
    assert tv_high < tv_low


def test_three_models_have_correct_max_size():
    rng = np.random.default_rng(5)
    counts = _gaussian_counts(301, 150, 30, 5000, rng)
    glo = build_global_fl(counts, max_size=300)
    rna = build_rna_fl(counts, counts, max_size=300, prior_ess=100.0)
    gdna = build_gdna_fl(counts, counts, max_size=300, prior_ess=100.0)
    assert glo.max_size == rna.max_size == gdna.max_size == 300
    assert glo._log_prob.size == rna._log_prob.size == gdna._log_prob.size == 301
