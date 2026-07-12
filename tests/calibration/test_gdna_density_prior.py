"""Phase-2 gDNA-density mixture prior — the substrate extractor + the KDE fit
(`calibration.gdna_density_prior`; design `docs/calibration/PHASE2_gdna_mixture_fit_design.md`)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.bp_solver import BOUNDARY, REGION, NodeBelief, NodeGeometry
from rigel.calibration.gdna_density_prior import (
    KIND_EXON,
    GdnaDensityPrior,
    TrainingSubstrate,
    _weighted_kde_logpdf,
    _weighted_median,
    build_training_substrate,
)
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, TS_AMBIG, TS_POS


# --------------------------------------------------------------------------- KDE (P2.1)


def _synthetic(x, *, w=None, std=None, kind=None):
    x = np.asarray(x, dtype=np.float64)
    w = np.ones_like(x) if w is None else np.asarray(w, dtype=np.float64)
    std = np.full_like(x, 0.05) if std is None else np.asarray(std, dtype=np.float64)
    kind = np.zeros(x.shape[0], dtype=np.int64) if kind is None else np.asarray(kind)
    return TrainingSubstrate(
        log_rho=x, weight=w, node_kind=kind, node_index=np.arange(x.shape[0]), log_rho_std=std
    )


def test_kde_recovers_bimodal_modes():
    """Two well-separated clusters (depleted ⊕ enriched) → the fit is bimodal with modes at both levels."""
    x = np.concatenate([np.linspace(-8.2, -7.8, 20), np.linspace(-1.2, -0.8, 20)])
    pr = GdnaDensityPrior.fit(_synthetic(x, std=np.full(x.size, 0.05)), bandwidth="silverman")
    top2 = sorted(m[0] for m in pr.modes[:2])
    assert len(pr.modes) >= 2, pr.modes
    assert -8.6 < top2[0] < -7.4, top2
    assert -1.6 < top2[1] < -0.4, top2


def test_kde_unimodal_when_single_cluster():
    """One tight cluster (uniform gDNA) → a single dominant mode near the cluster centre."""
    x = np.linspace(-2.3, -1.7, 40)
    pr = GdnaDensityPrior.fit(_synthetic(x, std=np.full(x.size, 0.2)), bandwidth="silverman")
    assert -2.2 < pr.modes[0][0] < -1.8, pr.modes


def test_logpdf_matches_direct_kde():
    """`logpdf` (grid interpolation) agrees with a direct weighted-KDE evaluation."""
    x = np.array([-3.0, -2.0, -1.0, 0.0, 1.0])
    w = np.array([1.0, 2.0, 3.0, 2.0, 1.0])
    pr = GdnaDensityPrior.fit(_synthetic(x, w=w, std=np.full(5, 0.1)), bandwidth=0.5, n_grid=4096)
    q = np.array([-2.5, -1.5, -0.5, 0.5])
    direct = _weighted_kde_logpdf(q, pr.train_x, pr.train_w, pr.bandwidth)
    assert np.allclose(pr.logpdf(q), direct, atol=0.02), (pr.logpdf(q), direct)


def test_bandwidth_noise_floor():
    """An ultra-tight cluster with LARGE per-node sampling std must NOT be resolved into spurious spikes —
    the bandwidth is floored at the density-noise scale (the fix for the uniform-gDNA fracturing)."""
    x = np.linspace(-2.01, -1.99, 60)  # spread 0.02 — far below the noise
    pr = GdnaDensityPrior.fit(_synthetic(x, std=np.full(x.size, 0.3)), bandwidth="silverman")
    assert pr.bandwidth >= 0.25, pr.bandwidth  # floored at ~median std 0.3
    assert len(pr.modes) == 1, pr.modes  # one smooth mode, no noise spikes


def test_fixed_bandwidth_and_lscv_run():
    x = np.concatenate([np.linspace(-6, -5.5, 15), np.linspace(-1.5, -1.0, 15)])
    sub = _synthetic(x, std=np.full(x.size, 0.05))
    assert abs(GdnaDensityPrior.fit(sub, bandwidth=0.4).bandwidth - 0.4) < 1e-6
    pr = GdnaDensityPrior.fit(
        sub, bandwidth="lscv"
    )  # LSCV must at least run + separate the two clusters
    assert len(pr.modes) >= 2, pr.modes


def test_floor_anchor_adds_depleted_mass():
    """The optional ρ_floor virtual sample seeds the depleted mode when depleted training nodes are absent."""
    x = np.linspace(-1.2, -0.8, 30)  # only enriched training nodes
    sub = _synthetic(x, std=np.full(x.size, 0.05))
    pr = GdnaDensityPrior.fit(
        sub, bandwidth=0.2, floor_log_rho=-8.0, floor_weight=float(np.sum(sub.weight))
    )
    assert (
        pr.logpdf(np.array([-8.0]))[0] > pr.logpdf(np.array([-5.0]))[0]
    )  # mass appears at the floor


def test_weighted_median_is_continuous():
    """The bandwidth-floor weighted median must be CONTINUOUS in its inputs (interpolated, not the old
    ``searchsorted`` step) — an ε change in a weight/value moves the result by ε, not discretely. This is the
    calibration cross-process determinism fix (a 1e-15 nudge must not flip the KDE bandwidth)."""
    rng = np.random.default_rng(0)
    v = np.sort(rng.normal(size=201))
    w = np.ones_like(v)
    m0 = _weighted_median(v, w)
    w2 = w.copy()
    w2[100] += 1e-12
    assert abs(_weighted_median(v, w2) - m0) < 1e-9  # continuous in weight
    v2 = v.copy()
    v2[100] += 1e-12
    assert abs(_weighted_median(v2, w) - m0) < 1e-9  # continuous in value
    assert (
        abs(_weighted_median(np.array([0.0, 1.0]), np.array([1.0, 1.0])) - 0.5) < 1e-9
    )  # even split → midpoint


def test_empty_substrate_raises():
    with pytest.raises(ValueError):
        GdnaDensityPrior.fit(_synthetic(np.zeros(0)))


# --------------------------------------------------------------------------- genome-scale memory fix


def test_weighted_kde_logpdf_query_tiling_is_exact():
    """Query-axis tiling (``eval_chunk``) is a pure memory knob — it must be bit-identical to a single
    untiled block. Regression for the genome-scale OOM: a ~90M-point query (nodes × grid) previously built
    a single (n_eval, n_samp) matrix = 2.69 TiB and crashed; tiling bounds it without changing the result."""
    rng = np.random.default_rng(3)
    xe = rng.normal(size=40_000)
    x = rng.normal(size=5_000)
    w = rng.random(5_000) + 0.1
    tiled = _weighted_kde_logpdf(xe, x, w, 0.4, eval_chunk=4096)
    untiled = _weighted_kde_logpdf(xe, x, w, 0.4, eval_chunk=10**9)
    assert np.array_equal(tiled, untiled)


def test_kde_logprior_tabulation_matches_direct_at_scale():
    """At genome scale ``_kde_logprior`` tabulates the exact kernel on a lattice spanning the query range
    and interpolates (instead of evaluating the KDE at every one of m·K points). For a large query it must
    agree with a direct per-point kernel evaluation to well within the kernel bandwidth — the real quadratic
    tails are preserved because the lattice covers the full range (no clamping)."""
    from rigel.calibration.bp_solver import _kde_logprior

    xs = np.concatenate([np.linspace(-9.0, -8.5, 40), np.linspace(-1.5, -1.0, 40)])
    pr = GdnaDensityPrior.fit(_synthetic(xs, std=np.full(xs.size, 0.1)), bandwidth="silverman")
    rng = np.random.default_rng(4)
    m = 50_000
    mass = rng.uniform(1.0, 500.0, size=m)
    eff = rng.uniform(50.0, 5_000.0, size=m)
    fgg = np.linspace(1e-3, 1.0 - 1e-3, 60)  # m·K = 3M ≫ lattice ⇒ tabulation path
    out = _kde_logprior(fgg, mass, eff, pr)
    assert out.shape == (m, 60)
    assert np.all(np.isfinite(out))
    # direct reference on a subset (exact per-point kernel) — must match the interpolated result
    sub = slice(0, 400)
    log_me = np.log(mass[sub]) - np.log(eff[sub])
    fg = np.clip(fgg, 1e-9, 1.0 - 1e-9)
    log_rho = np.log(fg)[None, :] + log_me[:, None]
    ref = pr.logpdf_kernel(log_rho.ravel()).reshape(log_rho.shape) - np.log1p(-fg)[None, :]
    assert (
        np.max(np.abs(out[sub] - ref)) < 5e-3
    )  # interp error ≪ bandwidth; a broken lattice would be O(1)


# --------------------------------------------------------------------------- substrate (P2.0)


def _mock_chain_beliefs(kappa_marks_ambig: bool):
    """A 5-node chain B0 R0(SS exon) B1 R1(AMBIG exon) B2 with known densities. include_boundaries=False in
    the tests so only the two region nodes can teach — R0 (single-strand) does, R1 (AMBIG) must not."""
    n = 5
    kind = np.array([BOUNDARY, REGION, BOUNDARY, REGION, BOUNDARY])
    ref_idx = np.array([0, 0, 1, 1, 2])
    chain = SimpleNamespace(
        kind=kind,
        ref_idx=ref_idx,
        n_nodes=n,
        left=np.array([-1, 0, 1, 2, 3]),
        right=np.array([1, 2, 3, 4, -1]),
    )
    # region signatures: R0 exon(+), R1 exon(+/−) = AMBIG
    region_arrays = SimpleNamespace(
        signature=np.array([BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG], dtype=np.int64),
        strand_class=np.array([TS_POS, TS_AMBIG], dtype=np.int8),
    )
    # per NODE free masks: R0 single-strand (pos only), R1 AMBIG (both)
    statics = SimpleNamespace(
        free_pos=np.array([False, True, False, True, False]),
        free_neg=np.array([False, False, False, True, False]),
    )
    z = np.zeros(n)
    M = np.array([0.0, 500.0, 0.0, 800.0, 0.0])  # mass on R0, R1
    Eg = np.array([1.0, 100.0, 1.0, 100.0, 1.0])  # gDNA eff-len
    geometry = NodeGeometry(
        n_nodes=n,
        mass_left=M.copy(),
        mass_right=M.copy(),
        eff_gdna_left=Eg.copy(),
        eff_gdna_right=Eg.copy(),
        eff_rna_left=Eg.copy(),
        eff_rna_right=Eg.copy(),
        eff_spl_left=np.ones(n),
        eff_spl_right=np.ones(n),
        spliced_pos_left=z.copy(),
        spliced_pos_right=z.copy(),
        spliced_neg_left=z.copy(),
        spliced_neg_right=z.copy(),
    )
    belief = NodeBelief(
        f_pos=np.array([0.0, 0.8, 0.0, 0.5, 0.0]),
        f_neg=z.copy(),
        f_g=np.array([1.0, 0.2, 1.0, 0.5, 1.0]),  # R0 f_g=0.2, R1 f_g=0.5
        var_pos=z.copy(),
        var_neg=z.copy(),
        var_gdna=np.array([0.0, 0.01, 0.0, 0.01, 0.0]),
    )
    bsub = SimpleNamespace(
        left_region=np.array([-1, 0, 1, 2, 3]), right_region=np.array([0, 1, 2, 3, -1])
    )
    return chain, belief, geometry, statics, region_arrays, bsub


def test_substrate_excludes_ambig_and_matches_density():
    chain, belief, geom, st, ra, bsub = _mock_chain_beliefs(True)
    sub = build_training_substrate(chain, belief, geom, st, ra, bsub, include_boundaries=False)
    assert sub.n == 1  # only R0 (SS exon); R1 (AMBIG) excluded
    assert sub.node_index[0] == 1  # chain id of R0
    assert sub.node_kind[0] == KIND_EXON
    # log ρ_g = log(f_g·M/E) = log(0.2·500/100) = log(1.0) = 0
    assert abs(sub.log_rho[0] - np.log(0.2 * 500.0 / 100.0)) < 1e-9


def test_substrate_weight_is_unit():
    """Every training node carries UNIT weight — precision is NOT used as a weight (it correlates with the
    density and would bias the KDE shape; design §8e)."""
    chain, belief, geom, st, ra, bsub = _mock_chain_beliefs(True)
    sub = build_training_substrate(chain, belief, geom, st, ra, bsub, include_boundaries=False)
    assert np.all(sub.weight == 1.0), sub.weight


def test_substrate_excludes_tiny_regions():
    """A region shorter than the fragment (E_gdna < min_eff_length) can't contain one → excluded (the
    1/E-blowup fix). R0 has E_gdna=100; excluding below 200 drops it → empty substrate."""
    chain, belief, geom, st, ra, bsub = _mock_chain_beliefs(True)
    sub = build_training_substrate(
        chain, belief, geom, st, ra, bsub, min_eff_length=200.0, include_boundaries=False
    )
    assert sub.n == 0, sub.n
