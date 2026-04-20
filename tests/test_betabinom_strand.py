"""Tests for the Beta-Binomial strand LLR pilot (v3 re-adoption).

See ``docs/calibration/calibration_v3_vs_v4_comparison.md`` and
``docs/calibration/betabinom_strand_pilot_plan.md``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._em import (
    _estimate_kappa_marginal,
    _strand_llr,
    _strand_llr_betabinom,
    run_em,
)


def _make_stats(
    n_pos: np.ndarray,
    n_neg: np.ndarray,
    tx_strand: np.ndarray,
    n_spliced: np.ndarray | None = None,
    mappable_bp: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    n_pos = np.asarray(n_pos, dtype=np.float64)
    n_neg = np.asarray(n_neg, dtype=np.float64)
    tx_strand = np.asarray(tx_strand, dtype=np.int8)
    N = n_pos.shape[0]
    n_unspliced = n_pos + n_neg
    if n_spliced is None:
        n_spliced = np.zeros(N, dtype=np.float64)
    if mappable_bp is None:
        mappable_bp = np.full(N, 1000.0, dtype=np.float64)
    return {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_unspliced": n_unspliced,
        "n_spliced": np.asarray(n_spliced, dtype=np.float64),
        "n_total": n_unspliced + np.asarray(n_spliced, dtype=np.float64),
        "tx_strand": tx_strand,
        "mappable_bp": np.asarray(mappable_bp, dtype=np.float64),
        "region_length": np.asarray(mappable_bp, dtype=np.float64),
        "ref": np.full(N, "chr1", dtype=object),
    }


def test_betabinom_llr_zero_at_ss_one_half():
    """SS = 0.5 ⇒ G and R coincide ⇒ LLR ≡ 0 regardless of κ or counts."""
    rng = np.random.default_rng(0)
    n_pos = rng.integers(0, 50, size=100).astype(np.float64)
    n_neg = rng.integers(0, 50, size=100).astype(np.float64)
    gs = rng.choice([-1, 1], size=100).astype(np.int8)
    stats = _make_stats(n_pos, n_neg, gs)
    for kappa in (0.5, 5.0, 50.0, 500.0):
        llr = _strand_llr_betabinom(stats, 0.5, kappa, enabled=True)
        assert np.allclose(llr, 0.0, atol=1e-10), f"κ={kappa}"


def test_betabinom_llr_disabled_is_zero():
    stats = _make_stats([10, 3], [0, 7], [1, -1])
    assert np.all(_strand_llr_betabinom(stats, 0.99, 10.0, enabled=False) == 0.0)
    # kappa<=0 also zeros out
    assert np.all(_strand_llr_betabinom(stats, 0.99, 0.0, enabled=True) == 0.0)


def test_betabinom_gating_matches_binomial():
    """Same gating as Binomial: n_unspliced < 2 or tx_strand == 0 → LLR = 0."""
    stats = _make_stats(
        n_pos=[0, 1, 5, 10],
        n_neg=[0, 0, 5, 10],
        tx_strand=[1, 1, 0, 1],
    )
    llr_bb = _strand_llr_betabinom(stats, 0.95, 20.0, enabled=True)
    llr_b = _strand_llr(stats, 0.95, enabled=True, noise_floor=1e-3)
    # Regions 0 (n=0), 1 (n=1 < 2), 2 (tx_strand 0) → exactly 0 in both.
    assert llr_bb[0] == 0.0 and llr_bb[1] == 0.0 and llr_bb[2] == 0.0
    assert llr_b[0] == 0.0 and llr_b[1] == 0.0 and llr_b[2] == 0.0
    # Region 3 (n=20, half sense, half anti) carries nonzero LLR in both.
    # 50/50 strongly favours gDNA under SS=0.95 → LLR (G-R) > 0.
    assert llr_bb[3] > 0.0
    assert llr_b[3] > 0.0


def test_betabinom_large_kappa_approaches_binomial():
    """As κ → ∞ the Beta-Binomial LLR approaches the plain Binomial."""
    rng = np.random.default_rng(42)
    # Mid-coverage regions with a mix of sense/antisense counts.
    N = 200
    tot = rng.integers(5, 50, size=N).astype(np.float64)
    ss_true = 0.9
    k_sense = rng.binomial(tot.astype(int), ss_true).astype(np.float64)
    # Map to n_pos, n_neg under tx_strand = +1 convention:
    # k_sense for +strand genes = n_neg (R1-antisense), so n_neg = k_sense.
    gs = np.ones(N, dtype=np.int8)
    stats = _make_stats(n_pos=tot - k_sense, n_neg=k_sense, tx_strand=gs)

    llr_binom = _strand_llr(stats, ss_true, enabled=True, noise_floor=1e-6)
    llr_bb_large = _strand_llr_betabinom(
        stats, ss_true, kappa=1e6, enabled=True, noise_floor=1e-6
    )
    # They should agree to within ~1% on aggregate mass.
    agg_b = float(np.sum(np.abs(llr_binom)))
    agg_bb = float(np.sum(np.abs(llr_bb_large)))
    assert agg_b > 0
    assert abs(agg_bb - agg_b) / agg_b < 0.05, (agg_b, agg_bb)


def test_betabinom_smaller_kappa_tempers_llr():
    """Finite κ < ∞ dampens the per-region LLR vs the Binomial limit."""
    rng = np.random.default_rng(7)
    N = 200
    tot = rng.integers(10, 80, size=N).astype(np.float64)
    ss_true = 0.95
    # Inject overdispersion: per-region p_sense ~ Beta(a, b) with mean 0.95.
    a, b = 9.5, 0.5
    p_i = rng.beta(a, b, size=N)
    k_sense = rng.binomial(tot.astype(int), p_i).astype(np.float64)
    stats = _make_stats(n_pos=tot - k_sense, n_neg=k_sense, tx_strand=np.ones(N))

    llr_binom = _strand_llr(stats, ss_true, enabled=True, noise_floor=1e-6)
    llr_bb = _strand_llr_betabinom(
        stats, ss_true, kappa=5.0, enabled=True, noise_floor=1e-6
    )
    # With κ = 5 the temp ered mixture should carry less total |LLR|.
    assert float(np.sum(np.abs(llr_bb))) < float(np.sum(np.abs(llr_binom)))


def test_kappa_mle_recovers_true_dispersion_roughly():
    """Synthetic Beta-Binomial draw: golden-section recovers the truth order of magnitude."""
    rng = np.random.default_rng(1)
    N = 1000
    tot = rng.integers(20, 200, size=N).astype(np.float64)
    kappa_true = 20.0
    ss = 0.9
    # Pure RNA draws (gamma = 0):
    a_e, b_e = kappa_true * ss, kappa_true * (1.0 - ss)
    p_i = rng.beta(a_e, b_e, size=N)
    k_sense = rng.binomial(tot.astype(int), p_i).astype(np.float64)
    gs = np.ones(N, dtype=np.int8)
    stats = _make_stats(n_pos=tot - k_sense, n_neg=k_sense, tx_strand=gs)

    gamma_all_rna = np.zeros(N, dtype=np.float64)
    kappa_hat = _estimate_kappa_marginal(stats, gamma_all_rna, ss, noise_floor=1e-4)
    assert kappa_hat is not None
    # Allow a wide tolerance: golden-section + single draw noise.
    assert 5.0 < kappa_hat < 80.0, kappa_hat


def test_kappa_mle_returns_none_on_small_sample():
    stats = _make_stats([3, 0], [0, 0], [1, 0])  # only 1 eligible (gs=0 drops the 2nd)
    out = _estimate_kappa_marginal(stats, np.array([0.5, 0.5]), 0.9)
    assert out is None


def test_run_em_betabinom_mode_smoke():
    """End-to-end: mini synthetic stats + run_em in betabinom mode doesn't crash
    and returns kappa > 0 when strand channel is enabled."""
    rng = np.random.default_rng(3)
    N = 500
    mappable = np.full(N, 1000.0, dtype=np.float64)
    # Mix: 80% RNA (sense-biased), 20% gDNA (symmetric).
    is_gdna = rng.random(N) < 0.20
    gs = rng.choice([-1, 1], size=N).astype(np.int8)
    n_total = rng.poisson(5.0, size=N).astype(np.float64)
    n_total[is_gdna] = rng.poisson(3.0, size=int(is_gdna.sum())).astype(np.float64)
    p_sense = np.where(is_gdna, 0.5, rng.beta(90.0, 10.0, size=N))  # κ=100 for RNA
    n_total_int = n_total.astype(int)
    k_sense = np.array([rng.binomial(t, p) for t, p in zip(n_total_int, p_sense)], dtype=np.float64)
    # Convention: for gs=+1, k_sense = n_neg; for gs=-1, k_sense = n_pos.
    n_pos = np.where(gs == 1, n_total - k_sense, k_sense)
    n_neg = n_total - n_pos
    # Add a few hard-spliced anchors so the M-step has something for μ_R/σ_R.
    n_spliced = np.zeros(N, dtype=np.float64)
    rna_ix = np.where(~is_gdna)[0]
    anchor_ix = rng.choice(rna_ix, size=min(50, rna_ix.size), replace=False)
    n_spliced[anchor_ix] = rng.integers(1, 5, size=anchor_ix.size).astype(np.float64)

    stats = _make_stats(n_pos, n_neg, gs, n_spliced=n_spliced, mappable_bp=mappable)
    fit = run_em(
        stats,
        fl_table=None,
        strand_specificity=0.9,
        mean_frag_len=200.0,
        strand_noise_floor=1e-3,
        strand_llr_mode="betabinom",
        max_iterations=20,
        convergence_tol=1e-4,
    )
    assert fit.strand_llr_mode == "betabinom"
    if fit.strand_used:
        assert fit.kappa > 0.0
    assert fit.n_iter > 0
    assert np.all((fit.gamma >= 0.0) & (fit.gamma <= 1.0))


def test_run_em_binomial_is_unchanged_default():
    """Default mode is still 'binomial' and κ stays 0."""
    stats = _make_stats(
        n_pos=np.array([1.0, 10.0, 0.0, 5.0]),
        n_neg=np.array([0.0, 1.0, 5.0, 0.0]),
        tx_strand=np.array([1, 1, -1, 1]),
        n_spliced=np.array([0, 2, 0, 0]),
        mappable_bp=np.array([1000.0, 1000.0, 1000.0, 1000.0]),
    )
    fit = run_em(
        stats,
        fl_table=None,
        strand_specificity=0.9,
        mean_frag_len=200.0,
    )
    assert fit.strand_llr_mode == "binomial"
    assert fit.kappa == 0.0
