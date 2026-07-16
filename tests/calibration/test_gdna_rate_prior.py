"""Unit tests for the pass-0 gDNA-rate prior (`GdnaRatePrior`, the Fixed-Kernel Poisson-lognormal Mixture)."""

from __future__ import annotations

import numpy as np

from rigel.calibration.gdna_rate_prior import GdnaRatePrior


def _bimodal_counts(seed=0):
    """8000 depleted (rate 1e-5, long E) + 2000 enriched (rate 1e-2, short E) Poisson counts — zero-inflated."""
    rng = np.random.default_rng(seed)
    eff = np.concatenate([rng.uniform(1e3, 3e3, 8000), rng.uniform(80, 200, 2000)])
    rate = np.concatenate([np.full(8000, 1e-5), np.full(2000, 1e-2)])
    count = rng.poisson(rate * eff).astype(float)
    return count, eff


def test_fit_is_deterministic():
    count, eff = _bimodal_counts()
    a = GdnaRatePrior.fit(count, eff)
    b = GdnaRatePrior.fit(count, eff)
    assert np.array_equal(a.log_rho, b.log_rho)
    assert np.array_equal(a.logP, b.logP)  # pure EM + arithmetic ⇒ bit-identical, no RNG


def test_recovers_two_modes():
    count, eff = _bimodal_counts()
    pr = GdnaRatePrior.fit(count, eff, bandwidth=0.15)
    p = np.exp(pr.logP - pr.logP.max())
    modes = pr.log_rho[1:-1][(p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:])] / np.log(10.0)
    assert modes.size >= 2  # a depleted + an enriched mode
    # the enriched mode sits near the injected 1e-2 (log10 ≈ −2); the depleted well below it.
    assert modes.max() > -2.6
    assert modes.min() < modes.max() - 1.0


def test_zero_inflated_native():
    count, eff = _bimodal_counts()
    assert (count == 0).mean() > 0.5  # the substrate really is zero-dominated
    pr = GdnaRatePrior.fit(count, eff)  # count-0 handled natively (Poisson e^{−ρE}); must not error
    assert np.isfinite(pr.logP).all()


def test_projection_is_bare_no_reference_prior():
    """``logprior`` is ``log P(log ρ_g)`` and NOTHING else — no reference prior, no measure, no Jacobian
    (`prior_ramp_and_bp_roadmap.md` §2). Regression guard for the ``+0.5·λ`` ramp: any residual linear term
    in λ is an improper, curvature-free pull whose strength is set only by the grid width."""
    count, eff = _bimodal_counts()
    pr = GdnaRatePrior.fit(count, eff)
    lam = np.linspace(-10, 10, 60)
    fg = 1.0 / (1.0 + np.exp(-lam))
    mass, e = count[:50] + 1.0, eff[:50]
    term = pr.logprior(fg, mass, e)
    assert term.shape == (50, 60)
    assert np.isfinite(term).all()
    # It must equal the interpolated log-density at ρ_g = f_g·M/E exactly — no additive f_g-only term.
    log_rho = np.log(fg)[None, :] + (np.log(mass) - np.log(e))[:, None]
    expect = np.interp(log_rho.ravel(), pr.log_rho, pr.logP, left=pr.logP[0], right=pr.logP[-1])
    assert np.allclose(term, expect.reshape(term.shape), atol=1e-12)
    # A node whose density is far ABOVE the fitted support must be penalised at f_g→1, not rewarded.
    # (Under the ramp this was inverted: the +0.5·λ pull beat the kde's own tail.)
    dense = pr.logprior(fg, np.array([mass.max() * 50.0]), np.array([e[0]]))[0]
    assert dense[-1] <= dense[len(fg) // 2] + 1e-9


def test_prior_is_weak():
    """The projected prior must stay weak — worth ~a read, so the strand can dominate wherever it has any
    information at all (`CALIBRATION_ARCHITECTURE.md`, count-zero-information).

    The bound is **2** pseudo-observations, not 1. The old <1 threshold was measured against a prior that
    silently carried the ``+0.5·λ`` ramp, and **a ramp has zero curvature by construction** — it widened
    Var(λ) and dragged the mass to the vertex, so the prior *scored* weak while in fact overwhelming every
    unstranded node (it alone returned f_g = 0.9994). Stripping it reveals the kde's true curvature:
    n_eff ≈ 1.2. That is still weak — the strand at κ=0.99, n=100 is worth ~5.3 — but it is honest.
    See prior_ramp_and_bp_roadmap.md §2.
    """
    count, eff = _bimodal_counts()
    pr = GdnaRatePrior.fit(count, eff)
    lam = np.linspace(-10, 10, 400)
    fg = 1.0 / (1.0 + np.exp(-lam))
    nz = count > 0
    t = pr.logprior(fg, count[nz], eff[nz])  # bare — nothing to strip
    w = np.exp(t - t.max(axis=1, keepdims=True))
    w /= w.sum(axis=1, keepdims=True)
    m = (w * lam[None, :]).sum(1)
    var = (w * (lam[None, :] - m[:, None]) ** 2).sum(1)
    n_eff = (1.0 / np.maximum(var, 1e-9)) / 0.25  # pseudo-observations
    assert np.median(n_eff) < 2.0
    assert (n_eff > 4.0).mean() < 0.05  # no node may rival a well-stranded count


def test_lognormal_belief_width_path():
    """The refit path (nonzero belief width τ²=var_g, via Gauss-Hermite) fits a valid, deterministic P(ρ).
    Note: unlike a non-deconvolving average, the fixed-kernel EM DECONVOLVES the belief width out — a wide τ
    is attributed to observation noise, not true rate spread — so uncertain nodes are down-weighted rather
    than broadening P(ρ). Here we only assert the path is valid + reproducible; the deconvolution behaviour is
    exercised end-to-end by the calibration/oracle benchmarks."""
    count, eff = _bimodal_counts()
    var_g = np.full_like(count, 1.5)
    a = GdnaRatePrior.fit(count, eff, var_g=var_g)
    b = GdnaRatePrior.fit(count, eff, var_g=var_g)
    assert np.isfinite(a.logP).all()
    assert np.array_equal(a.logP, b.logP)  # deterministic
    # differs from the pure-Poisson (τ=0) fit — the belief width genuinely enters the likelihood.
    assert not np.allclose(a.logP, GdnaRatePrior.fit(count, eff).logP)


def test_degenerate_all_zero():
    """All count-0 (a truly gDNA-free library) must fit without error and place mass at the low end."""
    eff = np.full(1000, 2000.0)
    pr = GdnaRatePrior.fit(np.zeros(1000), eff)
    assert np.isfinite(pr.logP).all()
    # mass concentrates in the lower half of the grid (ρ→0), not the top.
    p = np.exp(pr.logP - pr.logP.max())
    assert (p[: len(p) // 2].sum()) > (p[len(p) // 2 :].sum())
