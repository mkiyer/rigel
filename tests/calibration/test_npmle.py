"""Unit tests for the enrichment NPMLE (`DensityNPMLE`, the Fixed-Kernel Poisson-lognormal Mixture)."""

from __future__ import annotations


import numpy as np

from rigel.calibration.background_reference import BackgroundReference
from rigel.calibration.npmle import DensityNPMLE


def _bimodal_counts(seed=0):
    """8000 depleted (rate 1e-5, long E) + 2000 enriched (rate 1e-2, short E) Poisson counts — zero-inflated."""
    rng = np.random.default_rng(seed)
    eff = np.concatenate([rng.uniform(1e3, 3e3, 8000), rng.uniform(80, 200, 2000)])
    rate = np.concatenate([np.full(8000, 1e-5), np.full(2000, 1e-2)])
    count = rng.poisson(rate * eff).astype(float)
    return count, eff


def test_fit_is_deterministic():
    count, eff = _bimodal_counts()
    a = DensityNPMLE.fit(count, eff)
    b = DensityNPMLE.fit(count, eff)
    assert np.array_equal(a.log_rho, b.log_rho)
    assert np.array_equal(a.logP, b.logP)  # pure EM + arithmetic ⇒ bit-identical, no RNG


def test_recovers_two_modes():
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff, bandwidth=0.15)
    p = np.exp(pr.logP - pr.logP.max())
    modes = pr.log_rho[1:-1][(p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:])] / np.log(10.0)
    assert modes.size >= 2  # a depleted + an enriched mode
    # the enriched mode sits near the injected 1e-2 (log10 ≈ −2); the depleted well below it.
    assert modes.max() > -2.6
    assert modes.min() < modes.max() - 1.0


def test_zero_inflated_native():
    count, eff = _bimodal_counts()
    assert (count == 0).mean() > 0.5  # the substrate really is zero-dominated
    pr = DensityNPMLE.fit(count, eff)  # count-0 handled natively (Poisson e^{−ρE}); must not error
    assert np.isfinite(pr.logP).all()


def test_projection_is_bare_no_reference_prior():
    """``logprior`` is ``log P(log ρ_g)`` and NOTHING else — no reference prior, no measure, no Jacobian
    Regression guard for the ``+0.5·λ`` ramp: any residual linear term
     in λ is an improper, curvature-free pull whose strength is set only by the grid width."""
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff)
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
    # A region whose density is far ABOVE the fitted support must be penalised at f_g→1, not rewarded.
    # (Under the ramp this was inverted: the +0.5·λ pull beat the kde's own tail.)
    dense = pr.logprior(fg, np.array([mass.max() * 50.0]), np.array([e[0]]))[0]
    assert dense[-1] <= dense[len(fg) // 2] + 1e-9


def test_prior_is_weak():
    """The projected prior must stay weak — worth ~a read, so the strand can dominate wherever it has any
    information at all (count-zero-information).

    The bound is **2** pseudo-observations, not 1. The old <1 threshold was measured against a prior that
    silently carried the ``+0.5·λ`` ramp, and **a ramp has zero curvature by construction** — it widened
    Var(λ) and dragged the mass to the vertex, so the prior *scored* weak while in fact overwhelming every
    unstranded region (it alone returned f_g = 0.9994). Stripping it reveals the kde's true curvature:
    n_eff ≈ 1.2. That is still weak — the strand at κ=0.99, n=100 is worth ~5.3 — but it is honest.
    See
    """
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff)
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
    assert (n_eff > 4.0).mean() < 0.05  # no region may rival a well-stranded count


def test_lognormal_belief_width_path():
    """The refit path (nonzero belief width τ²=var_g, via Gauss-Hermite) fits a valid, deterministic P(ρ).
    Note: unlike a non-deconvolving average, the fixed-kernel EM DECONVOLVES the belief width out — a wide τ
    is attributed to observation noise, not true rate spread — so uncertain regions are down-weighted rather
    than broadening P(ρ). Here we only assert the path is valid + reproducible; the deconvolution behaviour is
    exercised end-to-end by the calibration/oracle benchmarks."""
    count, eff = _bimodal_counts()
    var_g = np.full_like(count, 1.5)
    a = DensityNPMLE.fit(count, eff, var_g=var_g)
    b = DensityNPMLE.fit(count, eff, var_g=var_g)
    assert np.isfinite(a.logP).all()
    assert np.array_equal(a.logP, b.logP)  # deterministic
    # differs from the pure-Poisson (τ=0) fit — the belief width genuinely enters the likelihood.
    assert not np.allclose(a.logP, DensityNPMLE.fit(count, eff).logP)


def test_degenerate_all_zero():
    """All count-0 (a truly gDNA-free library) must fit without error and place mass at the low end."""
    eff = np.full(1000, 2000.0)
    pr = DensityNPMLE.fit(np.zeros(1000), eff)
    assert np.isfinite(pr.logP).all()
    # mass concentrates in the lower half of the grid (ρ→0), not the top.
    p = np.exp(pr.logP - pr.logP.max())
    assert (p[: len(p) // 2].sum()) > (p[len(p) // 2 :].sum())


# the projection message-transfer-variance (design) ----


def test_projection_floor_is_h_squared():
    """Deep inside a single mode the responsibilities concentrate ⇒ ``var_proj`` sits at the within-mode floor
    ``h²`` (the count-zero-information max-precision cap). ``var_proj ≥ h²`` everywhere, and it scales with h."""
    count, eff = _bimodal_counts()
    for h in (0.15, 0.25):
        pr = DensityNPMLE.fit(count, eff, bandwidth=h)
        h2 = (h * np.log(10.0)) ** 2
        _mu, var = pr.project(count + 1.0, eff)
        assert np.isfinite(var).all()
        assert (var >= h2 - 1e-9).all()  # the within-mode floor is the hard lower bound
        # a region deep in the enriched mode (rate 1e-2) sits at ~the floor (unambiguous membership)
        assert pr.project(np.array([1e-2 * 150.0]), np.array([150.0]))[1][0] < 3.0 * h2


def test_projection_crossing_dominates_within_mode():
    """The per-boundary transfer variance TRAPS: specificity-and-sense-are-complements = ``var_dst + (mu_dst − mu_src)²``: a depleted↔enriched crossing is
    orders of magnitude larger than a within-mode boundary — this is the enrichment-crossing damping that gags a
    message across a capture boundary while letting same-mode messages flow."""
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff, bandwidth=0.15)
    E = np.array([1500.0, 150.0])
    dep = np.array([1e-5 * 1500.0, 1e-5 * 150.0])  # two depleted (different E)
    mix = np.array([1e-5 * 1500.0, 1e-2 * 150.0])  # depleted ↔ enriched
    mu_d, var_d = pr.project(dep, E)
    mu_m, var_m = pr.project(mix, E)
    within = var_d[1] + (mu_d[1] - mu_d[0]) ** 2
    crossing = var_m[1] + (mu_m[1] - mu_m[0]) ** 2
    assert crossing > 20.0 * within


def test_projection_deterministic():
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff)
    a0, a1 = pr.project(count[:200] + 1.0, eff[:200])
    b0, b1 = pr.project(count[:200] + 1.0, eff[:200])
    assert np.array_equal(a0, b0) and np.array_equal(a1, b1)


# the one-sided background floor (Phase 2) ----


# ---- the aggregate background CELL in the fit (the pooled scalar as one genome-length Poisson observation) ----


def test_aggregate_cell_concentrates_mass_at_background():
    """A present aggregate (count Σg over a genome-scale ΣE) injected as a Poisson cell concentrates prior mass
    at ``log ρ_bg = log(Σg/ΣE)`` — the smooth low-density anchor, sharp by construction (the huge ΣE)."""
    rng = np.random.default_rng(1)
    eff = rng.uniform(1000.0, 3000.0, 5000)
    count = rng.poisson(1e-2 * eff).astype(float)  # enriched individual regions at ρ ≈ 1e-2
    target = float(np.log(1e-3))  # the background sits a decade BELOW the individual regions
    bg = BackgroundReference(
        log_rho_bg=target,
        sigma_bg=0.001,
        n_counts=float(1e-3 * 6.0e6),
        eff_total=6.0e6,
        n_regions=2000,
    )
    pr_bg = DensityNPMLE.fit(count, eff, background=bg)
    pr_0 = DensityNPMLE.fit(count, eff)
    assert np.isfinite(pr_bg.logP).all()
    assert np.array_equal(
        pr_bg.logP, DensityNPMLE.fit(count, eff, background=bg).logP
    )  # deterministic
    near = np.abs(pr_bg.log_rho - target) < (
        1.5 * 0.15 * np.log(10.0)
    )  # within ~1 bandwidth of ρ_bg
    m_bg = float(np.exp(pr_bg.logP)[near].sum() / np.exp(pr_bg.logP).sum())
    # the plain fit's grid may not even reach ρ_bg — the aggregate both EXTENDS the grid there and puts mass on it
    m_0 = float(
        np.exp(pr_0.logP)[np.abs(pr_0.log_rho - target) < (1.5 * 0.15 * np.log(10.0))].sum()
        / np.exp(pr_0.logP).sum()
    )
    assert m_bg > m_0


def test_aggregate_cell_zero_counts_anchors_the_derived_floor():
    """The critical zero-DNA case: an aggregate with Σg=0 must anchor prior mass at the DERIVED resolution floor
    ``log_rho_floor`` — the per-region resolution wall ``ρ_res = 1/harmmean(E of zero-count regions)`` — NOT the
    old ``1/ΣE`` (which pools the genome as one region ⇒ ~3 logs too low = the confident-FP seed). No single
    region can resolve below its own ``1/E``; the honest floor is where a TYPICAL region still reads ~zero.
    (.)"""
    rng = np.random.default_rng(2)
    eff = rng.uniform(1000.0, 3000.0, 3000)
    count = rng.poisson(1e-2 * eff).astype(float)  # enriched individual regions
    rho_res = float(np.mean(1.0 / eff))  # mean(1/E_i) = 1/harmmean(E) — the resolution wall
    lrf = float(np.log(rho_res))
    bg0 = BackgroundReference(
        log_rho_bg=-np.inf,
        sigma_bg=np.inf,
        log_rho_floor=lrf,
        n_counts=0.0,
        eff_total=6.0e6,
        n_regions=3000,
    )
    pr = DensityNPMLE.fit(count, eff, background=bg0)
    assert np.isfinite(pr.logP).all()
    assert pr.log_rho[0] <= lrf + 3.0 * 0.15 * np.log(
        10.0
    )  # the grid reached down to the DERIVED floor
    p = np.exp(pr.logP - pr.logP.max())
    lowmass = float(p[pr.log_rho < lrf + 0.15 * np.log(10.0)].sum() / p.sum())
    assert (
        lowmass > 0.10
    )  # real prior mass sits at the derived floor (not all at the enriched individuals)
    # …and the derived floor is well ABOVE the old 1/ΣE collapse (the whole point of the fix):
    assert lrf > np.log(1.0 / 6.0e6) + 3.0  # >3 nats above 1/ΣE


def test_background_none_is_byte_identical():
    """``background=None`` must leave the fit EXACTLY as before the aggregate cell existed (safe default)."""
    count, eff = _bimodal_counts()
    a = DensityNPMLE.fit(count, eff)
    b = DensityNPMLE.fit(count, eff, background=None)
    assert np.array_equal(a.logP, b.logP) and np.array_equal(a.weights, b.weights)


def test_projection_is_continuous_across_the_valley():
    """Regression guard for the external reviewer's Risk 1 (step-function discontinuities → refit oscillation).

    The projection is a softmax over the fixed mixture, so ``mu_proj``/``var_proj`` are SMOOTH (Lipschitz)
    functions of the observed density — there are no regime strata, no thresholds. The discriminating property,
    with no magic threshold: for a Lipschitz function the max adjacent jump ``∝ grid step``, so HALVING the step
    roughly halves the max jump; a step-function's jump (≈ the mode gap, several nats) is resolution-INVARIANT.
    We sweep a region's density straight through the depleted↔enriched valley — the one place a stratified design
    would step — at two resolutions and assert the max jump shrinks like the step. Reintroducing discrete regime
    bins would make the ratio ≈ 1, failing this."""
    count, eff = _bimodal_counts()
    pr = DensityNPMLE.fit(count, eff, bandwidth=0.15)
    lo, hi = pr.log_rho[0] + 0.5, pr.log_rho[-1] - 0.5  # inside the fitted support

    def max_jump(n):
        ld = np.linspace(lo, hi, n)  # sweep across BOTH modes and the valley between
        mu, var = pr.project(np.exp(ld) * 1000.0, np.full(n, 1000.0))  # density = exp(ld)
        assert np.isfinite(mu).all() and np.isfinite(var).all()
        return np.abs(np.diff(mu)).max(), np.abs(np.diff(var)).max()

    mu_coarse, var_coarse = max_jump(2000)
    mu_fine, var_fine = max_jump(4000)  # half the step
    # Lipschitz ⇒ refining halves the jump (ratio → 0.5); a step-function ⇒ ratio → 1. Allow slack for the
    # discrete grid, but demand the jump clearly shrinks — anything ≤ 0.7 excludes a resolution-invariant step.
    assert mu_fine / max(mu_coarse, 1e-30) < 0.7
    assert var_fine / max(var_coarse, 1e-30) < 0.7


# ---------------------------------------------------------------------------
# The ADDITIVE Role-B representation: occupancy-weighted, fixed-bandwidth KDE +
# a weak 1-pseudo-observation floor. These pin the three properties the design exists to guarantee.
# ---------------------------------------------------------------------------


def _two_pop(n_dep=400, n_enr=100, tau_dep=0.40, tau_enr=0.80):
    """A depleted majority (ρ=1e-3) + an ENRICHED minority (ρ=1e-2) that is MORE imprecise (larger τ) — the
    regime where EM starves the minority. Returns ``(g_hat, eff, var_g)`` with ρ = g_hat/eff."""
    eff = np.full(n_dep + n_enr, 1000.0)
    g_hat = np.concatenate([np.full(n_dep, 1.0), np.full(n_enr, 10.0)])  # 1e-3 vs 1e-2
    var_g = np.concatenate([np.full(n_dep, tau_dep**2), np.full(n_enr, tau_enr**2)])
    return g_hat, eff, var_g


def _height(fit, log10_rho):
    x = np.asarray(fit.log_rho) / np.log(10.0)
    p = np.exp(fit.logP - fit.logP.max())
    return float(p[np.argmin(np.abs(x - log10_rho))])


def test_additive_default_off_is_the_em_path():
    """``additive=False`` must be the unchanged EM NPMLE (Role A / σ²_transfer relies on this)."""
    count, eff = _bimodal_counts()
    a = DensityNPMLE.fit(count, eff, bandwidth=0.15)
    b = DensityNPMLE.fit(count, eff, bandwidth=0.15, additive=False)
    assert np.array_equal(a.logP, b.logP) and np.array_equal(a.weights, b.weights)


def test_additive_preserves_a_starved_minority_mode():
    """The core guarantee: EM competes the enriched minority to ~zero height; the additive KDE keeps it."""
    g_hat, eff, var_g = _two_pop()
    em = DensityNPMLE.fit(g_hat, eff, var_g=var_g, bandwidth=0.15, additive=False)
    kde = DensityNPMLE.fit(g_hat, eff, var_g=var_g, bandwidth=0.15, additive=True)
    assert _height(em, -2.0) < 0.05  # EM starves the minority
    assert _height(kde, -2.0) > 0.15  # additive keeps it — the enriched mode survives


def test_additive_occupancy_equals_height_despite_imprecision():
    """Occupancy ratio ≈ rendered-height ratio (fixed bandwidth ⇒ no τ discounting), even though the minority
    is the more imprecise population. This is why TRAPS: perturb-every-gate (additive) needs TRAPS: a-gate-that-already-passed (fixed h)."""
    g_hat, eff, var_g = _two_pop(n_dep=400, n_enr=100, tau_enr=0.80)
    kde = DensityNPMLE.fit(g_hat, eff, var_g=var_g, bandwidth=0.15, additive=True)
    ratio = _height(kde, -2.0) / _height(kde, -3.0)
    assert 0.18 < ratio < 0.33  # ≈ the 100/400 = 0.25 occupancy ceiling (not τ-discounted toward 0)


def test_additive_weak_floor_cannot_dominate_the_flood():
    """TRAPS: self-checking-validator: the background floor is ONE pseudo-observation regardless of ``n_regions`` — a 100k-region intergenic
    flood must NOT crush the enriched mode (the exact real-data failure the EM aggregate cell would cause)."""
    g_hat, eff, var_g = _two_pop()
    bg = BackgroundReference(
        log_rho_bg=float(np.log(1e-3)),
        sigma_bg=float(1 / np.sqrt(50.0)),
        n_counts=50.0,
        eff_total=5e4,
        n_regions=100_000,
        log_rho_floor=float(np.log(1e-3)),
    )
    kde = DensityNPMLE.fit(g_hat, eff, var_g=var_g, bandwidth=0.15, additive=True, background=bg)
    assert _height(kde, -2.0) > 0.15  # enriched mode still present under a 100k-region flood


def test_additive_pure_rna_concentrates_low_no_phantom_mode():
    """FP safety: pure-RNA regions (ĝ≈0) read as gDNA ≤ 1/E ⇒ the density concentrates at the low resolution
    wall, with no manufactured high-gDNA mode."""
    eff = np.full(300, 1000.0)
    g_hat = np.zeros(300)  # pure RNA: no deconvolved gDNA
    kde = DensityNPMLE.fit(g_hat, eff, bandwidth=0.15, additive=True)
    x = np.asarray(kde.log_rho) / np.log(10.0)
    p = np.exp(kde.logP - kde.logP.max())
    assert (
        x[int(np.argmax(p))] < -2.5
    )  # the mode sits low (≈ log10(1/E) = −3), not at a phantom high rate


def test_additive_sub_one_counts_normalize_and_stay_low():
    """REGRESSION (adversarial review, 2026-07-19): the common mostly-RNA case has deconvolved gDNA counts
    ĝ = f_g·M that are positive-but-sub-1, so the resolution-floored kernel centres (log(1/E)) sit ABOVE a
    raw-density grid — leaving every kernel off-grid, underflowing the density (phantom high mode + an
    unclamped-blo IndexError). The grid must span the floored centres. Assert: no crash, normalized, low mode."""
    for g in (
        np.full(300, 0.02),
        np.concatenate([np.full(120, 0.02), np.zeros(180)]),
        np.full(300, 0.001),
    ):
        kde = DensityNPMLE.fit(g, np.full(300, 1000.0), bandwidth=0.15, additive=True)
        p = np.exp(kde.logP - kde.logP.max())
        assert (
            abs(float(kde.weights.sum()) - 1.0) < 1e-9
        )  # normalized density (the bug collapsed it to ~2e-4)
        assert np.isfinite(kde.logP).all()
        x = np.asarray(kde.log_rho) / np.log(10.0)
        assert (
            x[int(np.argmax(p))] < -2.0
        )  # mode LOW (≈ 1/E resolution wall), not a phantom high-gDNA mode


def test_additive_sub_one_counts_with_resolution_wall_floor():
    """The Σg=0 / gdna_none background (floor at the resolution wall, sigma_bg=∞) must also handle sub-1 ĝ
    without crashing — the floor location must land on the grid and h_floor must be finite."""
    bg = BackgroundReference(
        log_rho_bg=-np.inf,
        sigma_bg=np.inf,
        n_counts=0.0,
        eff_total=1e6,
        n_regions=5000,
        log_rho_floor=float(np.log(1.0 / 3000.0)),  # the resolution wall
    )
    kde = DensityNPMLE.fit(
        np.full(300, 0.001), np.full(300, 1000.0), bandwidth=0.15, additive=True, background=bg
    )
    assert abs(float(kde.weights.sum()) - 1.0) < 1e-9
    assert np.isfinite(kde.logP).all()
