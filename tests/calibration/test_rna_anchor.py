"""Gates for the RNA-anchored evidence factor (`calibration.rna_anchor`), round-2 form.

The factor anchors the RNA side of the unspliced count on quantities hybrid capture cannot
mis-scale. Round 2 (owner-approved, 2026-08-24, after the adversarial design review): the exon flank rate is the SUM of per-route junction rates
(the review-confirmed pooling bug fix), scoring is the NegBinomial marginal quadrature over the
transport scatter (median-preserving), and the nascent term is MARGINALIZED into the quadrature
via truncated-posterior quantile nodes (the plug-in and the truncated posterior MEAN are both
positively biased at nascent-free truth — Jensen; measured as the round-1 arm-B g98 regression).

Written against the new API BEFORE the module was rewritten; every gate was watched firing
against a deliberately broken build.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import rna_anchor as RA
from rigel.calibration.simplex_logodds import _logodds_grid

K = 61
WINDOW = 10.0
_, FG = _logodds_grid(K, WINDOW)


def _argbest_fg(row):
    return float(FG[int(np.argmax(row))])


# ── the vectorized NegBinomial and the s = 0 identity ────────────────────────────────────────────


def test_the_vectorized_nb_matches_the_shipped_one_exactly():
    """`_log_negbinom_rows` must be the shipped `density_deconv._log_negbinom` elementwise — the
    quadrature is that same marginal averaged over scatter, and any drift here is silent bias."""
    from rigel.calibration.density_deconv import _log_negbinom

    g = np.linspace(0.0, 200.0, 23)[None, :]
    for mu, size in ((20.0, 5.0), (0.5, 0.5), (150.0, 300.5)):
        a = _log_negbinom(g, mu, size)
        b = RA._log_negbinom_rows(g, np.array([[mu]]), np.array([[size]]))
        # scalar-vs-array numpy paths differ at ULP association order (~4e-14); the gate demands
        # identity to within that, far below any consumer's sensitivity
        np.testing.assert_allclose(a, b, rtol=0.0, atol=1e-10)


def test_zero_scatter_quadrature_reduces_to_the_pure_nb():
    """s = 0 ⇒ every node multiplier is 1 ⇒ the marginal IS the NegBinomial (the intron factory's
    own Gamma⊗Poisson pattern). Row-offset invariance is allowed; shape must be exact."""
    unspl = np.array([300.0, 40.0])
    mu = np.array([120.0, 35.0])
    size = np.array([80.5, 12.5])
    rows = RA._quadrature_rows(
        unspl, mu, size, scatter_log_variance=0.0, nascent_count_nodes=None, fg_grid=FG
    )
    ref = RA._log_negbinom_rows(
        np.clip(1.0 - FG, 1e-12, 1.0)[None, :] * unspl[:, None], mu[:, None], size[:, None]
    )
    ref = ref - ref.max(axis=1, keepdims=True)
    np.testing.assert_allclose(rows, ref, atol=1e-10)


def test_the_marginal_is_bulk_soft_and_tail_bounded():
    """The SHAPE property the panel actually rewarded (measured, not the design prose): the
    stratified equal-mass marginal is SOFTER than the count-scale Gaussian in the BULK — 1–2×
    misses, where the clean-library damage lives (claimed exons 5.5k vs 45k under a bulk-tight
    Gauss–Hermite weighting, measured) — while the tail stays within a bounded multiple of the
    Gaussian and never goes over-generous. A heavier-tailed scatter law is the recorded
    refinement for the deep tail."""
    size = np.array([200.5])
    V = 0.16  # the measured post-route-fix scatter scale
    mu = np.array([300.0])

    def costs(miss):
        C = np.array([mu[0] * miss])
        rows = RA._quadrature_rows(
            C, mu, size, scatter_log_variance=V, nascent_count_nodes=None, fg_grid=FG
        )
        quad = float(rows[0].max() - rows[0][0])
        var = float(mu[0] + np.expm1(V) * mu[0] ** 2 + mu[0] ** 2 / size[0])
        gauss = 0.5 * float(C[0] - mu[0]) ** 2 / var
        return quad, gauss

    q2, g2 = costs(2.0)
    assert q2 < 0.75 * g2, (
        f"2× miss: the marginal must be measurably SOFTER than the Gaussian in the bulk "
        f"({q2:.2f} vs {g2:.2f} nats) — bulk-tightness is the measured clean-library killer"
    )
    q3, g3 = costs(3.0)
    assert 0.25 * g3 < q3 < 2.5 * g3, (
        f"3× miss: the tail must stay within a bounded band of the Gaussian "
        f"({q3:.1f} vs {g3:.1f} nats) — neither fanatical nor free"
    )


# ── the route-sum pooling (the confirmed bug fix) ────────────────────────────────────────────────


def test_the_flank_rate_is_the_sum_of_route_rates_not_their_mean():
    """Two disjoint routes at one flank each carry molecules the exon contains, so the flank's
    molecule rate is the SUM of per-route rates. The round-1 ratio-of-sums (the opportunity-
    weighted MEAN) under-predicted k-route exons ~k× — measured median 2.0–2.3 at multi-route
    exons, 0.99–1.15 after this fix."""
    flux = np.array([100.0, 100.0])
    eff = np.array([200.0, 200.0])
    rate, flux_total = RA._flank_route_rate(flux, eff)
    assert np.isclose(rate, (100.0 + 0.25) / 200.0 * 2.0), (
        "two equal routes ⇒ twice one route's rate"
    )
    assert flux_total == 200.0
    # ratio-of-sums would give (200+0.5)/400 ≈ 0.501 — half the correct 1.0025
    assert rate > 1.9 * (200.5 / 400.0)


def test_single_route_flanks_are_unchanged_by_the_fix():
    flux = np.array([100.0])
    eff = np.array([200.0])
    rate, _ = RA._flank_route_rate(flux, eff)
    assert np.isclose(rate, 100.5 / 200.0)


# ── the marginalized nascent term ────────────────────────────────────────────────────────────────


def test_nascent_nodes_carry_an_atom_at_zero_for_a_clean_deep_intron():
    """A deep intron AT the background rate must contribute ~zero nascent at its lower nodes — the
    marginalization's whole point: the null keeps mass at exactly zero, where both the plug-in and
    the truncated posterior MEAN are strictly positive (the Jensen refutation, measured as the
    round-1 g98 regression)."""
    nodes = RA.nascent_rate_nodes(
        counts=np.array([10000.0]),
        gdna_opportunity=np.array([10000.0]),
        rna_opportunity=np.array([10000.0]),
        background_rate=1.0,
    )
    assert nodes.shape[1] >= 3
    assert nodes[0, 0] == 0.0, "the low node of an at-background intron must be exactly zero"
    assert nodes[0, -1] < 0.05, "even the high node must be within counting noise of zero"


def test_nascent_nodes_detect_a_real_excess_and_spread_with_noise():
    nodes_rich = RA.nascent_rate_nodes(
        counts=np.array([9000.0]),
        gdna_opportunity=np.array([10000.0]),
        rna_opportunity=np.array([10000.0]),
        background_rate=0.1,
    )
    assert nodes_rich[0, 0] > 0.5, "a 9x-over-background intron must read decisive nascent"
    shallow = RA.nascent_rate_nodes(
        counts=np.array([18.0]),
        gdna_opportunity=np.array([100.0]),
        rna_opportunity=np.array([100.0]),
        background_rate=0.1,
    )
    spread = shallow[0, -1] - shallow[0, 0]
    assert spread > 0.05, "a shallow intron's nascent must be a SPREAD, not a point"


def test_the_nascent_frame_uses_the_rna_opportunity_divisor():
    """The excess COUNTS convert to a rate with the intron's RNA opportunity (nascent fragments
    carry the RNA length distribution) — identical when the two opportunities coincide (this
    panel, BY DESIGN) and scaled by their ratio when they do not (real data; the fl-gap side
    panels are the falsification substrate). Review-refuted round 1 used the gDNA divisor."""
    a = RA.nascent_rate_nodes(
        counts=np.array([900.0]),
        gdna_opportunity=np.array([1000.0]),
        rna_opportunity=np.array([1000.0]),
        background_rate=0.1,
    )
    b = RA.nascent_rate_nodes(
        counts=np.array([900.0]),
        gdna_opportunity=np.array([1000.0]),
        rna_opportunity=np.array([500.0]),
        background_rate=0.1,
    )
    np.testing.assert_allclose(b, 2.0 * a, rtol=1e-9)


def test_marginalized_nascent_reaches_the_exon_factor():
    """With a nascent-rich adjacent intron the exon factor must tolerate un-fluxed RNA — nascent
    nodes must move the factor's peak toward RNA (TRAPS: could-the-arm-have-fired for the whole
    nascent path)."""
    C = np.array([1000.0])
    mu = np.array([300.0])
    size = np.array([200.5])
    no_nas = RA._quadrature_rows(
        C, mu, size, scatter_log_variance=0.04, nascent_count_nodes=None, fg_grid=FG
    )
    with_nas = RA._quadrature_rows(
        C,
        mu,
        size,
        scatter_log_variance=0.04,
        nascent_count_nodes=np.array([[600.0, 700.0, 800.0]]),
        fg_grid=FG,
    )
    assert _argbest_fg(with_nas[0]) < _argbest_fg(no_nas[0]) - 0.2, (
        "nascent nodes must move the factor's peak toward RNA"
    )


def test_the_boundary_nascent_frame_uses_the_rna_opportunity_divisor():
    """The boundary factor's divisor swap: halving the intron's RNA opportunity doubles the
    implied nascent rate, so the factor must tolerate more RNA at the boundary (peak at lower
    f_g). Identical when the opportunities coincide (this panel BY DESIGN); the fl-gap side
    panels are the real-data falsification substrate."""
    # sized so the base prediction sits at ~half the crossing count: the halved-opportunity arm
    # then doubles it to ~the full count, moving the peak measurably (an over-sized intron rate
    # saturates both arms at f≈0 and the gate reads nothing)
    kw = dict(
        unspl_boundary=np.array([500.0]),
        unspl_intron=np.array([520.0]),
        gdna_opportunity_intron=np.array([1000.0]),
        rna_opportunity_boundary=np.array([600.0]),
        background_rate=0.1,
        fg_grid=FG,
        pair_log_variance=0.02,
    )
    base = RA.boundary_anchor_rows(rna_opportunity_intron=np.array([1000.0]), **kw)
    halved = RA.boundary_anchor_rows(rna_opportunity_intron=np.array([500.0]), **kw)
    assert _argbest_fg(halved[0]) < _argbest_fg(base[0]) - 0.1, (
        "halving the intron RNA opportunity must move the boundary factor toward RNA"
    )


def test_the_boundary_goes_flat_on_a_capture_cliff_where_the_quadrature_asserts():
    """The capture cliff: boundary crossings far exceed the intron-derived prediction with
    decade-scale heterogeneity (capture enriches the crossing, depletes the intron the anchor
    reads, differently per exon). The Gaussian family's guarded fit reads that as an enormous
    spread and must go NEAR-FLAT; the quadrature on the same inputs confidently delivers gDNA —
    the priced 734 → 27.5k zero-control leak this gate exists to keep dead."""
    C_b = np.array([1.0, 3.0, 8.0, 40.0, 200.0, 15.0, 60.0, 2.0, 110.0, 25.0])
    C_i = np.ones(10)
    E_gi = np.full(10, 1000.0)
    rows = RA.boundary_anchor_rows(
        unspl_boundary=C_b,
        unspl_intron=C_i,
        gdna_opportunity_intron=E_gi,
        rna_opportunity_intron=np.full(10, 1000.0),
        rna_opportunity_boundary=np.full(10, 300.0),
        background_rate=0.0005,
        fg_grid=FG,
        pair_log_variance=None,
    )
    assert float(np.ptp(rows, axis=1).max()) < 2.0, (
        "a heterogeneous cliff must send the boundary factor near-flat, never a confident call"
    )
    # the contrast that makes this gate self-falsifying: the quadrature form on the SAME slots
    nodes = RA.nascent_rate_nodes(
        counts=C_i,
        gdna_opportunity=E_gi,
        rna_opportunity=np.full(10, 1000.0),
        background_rate=0.0005,
    )
    quad = RA._quadrature_rows(
        C_b,
        np.zeros(10),
        C_i + 0.5,
        scatter_log_variance=0.02,
        nascent_count_nodes=nodes * 300.0,
        fg_grid=FG,
    )
    assert float((quad[:, -1] - quad[:, 0]).max()) > 10.0, (
        "the quadrature must be visibly asserting here, or the flat assertion above is vacuous"
    )


def test_the_builder_routes_boundaries_through_the_gaussian_family(anchor_toy, monkeypatch):
    """The build path itself must call `boundary_anchor_rows` for boundary slots — a correct
    standalone function is DEAD CODE unless the builder dispatches to it (this exact hole shipped
    once: the reverted family existed while the builder still scored boundaries with the
    quadrature, and 17 green gates could not see it)."""
    t = anchor_toy
    sentinel = 7.25

    def stamped(**kw):
        rows = RA._gaussian_rows(
            np.asarray(kw["unspl_boundary"], np.float64),
            np.zeros(np.asarray(kw["unspl_boundary"]).shape[0]),
            np.ones(np.asarray(kw["unspl_boundary"]).shape[0]),
            kw["fg_grid"],
        )
        return np.full_like(rows, sentinel)

    monkeypatch.setattr(RA, "boundary_anchor_rows", stamped)
    out = RA.build_rna_anchor_factor(
        t["chain"],
        t["statics"],
        t["geometry"],
        t["region_arrays"],
        t["routes"],
        n_grid=K,
        logodds_window=WINDOW,
    )
    eligible = RA.eligible_slots(t["chain"], t["statics"], t["geometry"], t["region_arrays"])
    is_boundary = np.asarray(t["chain"].kind) == RA.BOUNDARY
    b_live = eligible & is_boundary
    assert np.any(b_live), "the toy must expose at least one eligible boundary"
    assert np.allclose(out[b_live], sentinel), (
        "boundary slots must be scored by boundary_anchor_rows, not an inline family"
    )


# ── the dispersion estimators ────────────────────────────────────────────────────────────────────


def test_route_pair_variance_is_conservative_and_refuses_small_populations():
    """The per-route two-flank estimator: MAD of log(rate_l/rate_r) halved, NO counting-noise
    subtraction — deliberately an over-estimate (the campaign's measured lesson: over-confidence
    is the killer, over-width is cheap). It must bracket an injected dispersion from above and
    refuse below the population minimum."""
    rng = np.random.default_rng(5)
    n = 800
    s_true = 0.3
    base = rng.uniform(0.5, 3.0, n)
    rl = base * np.exp(rng.normal(0, s_true, n))
    rr = base * np.exp(rng.normal(0, s_true, n))
    V = RA.route_pair_log_variance(rl, rr)
    assert V is not None
    assert s_true**2 * 0.7 < V < s_true**2 * 2.0, f"must bracket the injected s² ({V:.3f} vs 0.09)"
    assert RA.route_pair_log_variance(rl[:5], rr[:5]) is None


def test_left_fit_center_and_guard_survive_round_two():
    """The guarded lower-quantile center fit is unchanged machinery — re-pinned here: accepts a
    clean offset population, refuses a gDNA-saturated tail."""
    rng = np.random.default_rng(29)
    n = 2000
    mu = np.full(n, 600.0)
    clean = mu * 1.2 * np.exp(rng.normal(0.0, 0.3, n))
    fit = RA.left_fit_center_spread(clean, mu)
    assert fit is not None and 0.05 < fit[0] < 0.35
    mask = rng.random(n) < 0.97
    contaminated = mu * np.exp(rng.normal(0.0, 0.2, n))
    contaminated[mask] *= rng.uniform(2.0, 50.0, int(mask.sum()))
    assert RA.left_fit_center_spread(contaminated, mu) is None


# ── deliver / refute through the quadrature ──────────────────────────────────────────────────────


def test_zero_flux_delivers_and_full_flux_refutes():
    """The two ends of the mechanism: MU ≈ 0 ⇒ the count is gDNA (robust — no multiplicative
    scatter un-zeroes a zero prediction); MU ≈ C ⇒ the count is RNA."""
    C = np.array([500.0, 500.0])
    mu = np.array([0.5, 495.0])
    size = np.array([0.5, 400.5])
    rows = RA._quadrature_rows(
        C, mu, size, scatter_log_variance=0.36, nascent_count_nodes=None, fg_grid=FG
    )
    assert _argbest_fg(rows[0]) > 0.95, "a near-zero prediction must deliver gDNA"
    assert float(rows[0][-1] - rows[0][0]) > 10.0, "the deliver margin must be decisive"
    # refute: the f≈0 end must be decisively preferred over f≈1 AND near-optimal — the exact peak
    # may sit mid-grid because sharper (smaller-mu) mixture components have taller densities, and
    # a position-only assertion would read that artifact as a failure
    assert float(rows[1][0] - rows[1][-1]) > 10.0, "a full-count prediction must refute gDNA"
    assert float(rows[1].max() - rows[1][0]) < 2.0, "the RNA end must be near-optimal"


# ── assembly on a real chain ─────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def anchor_toy(tmp_path_factory):
    """The prior_vs_oracle toy scanned once, with chain pieces, the ROUTE TABLE, and full
    `calibrate` kwargs built the production way."""
    import dataclasses
    import importlib.util
    import sys
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "tpo_for_anchor", Path(__file__).parent / "test_prior_vs_oracle.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tpo_for_anchor"] = m
    spec.loader.exec_module(m)
    toy = m.toy.__wrapped__(tmp_path_factory)

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.region_chain import build_region_chain
    from rigel.calibration.region_geometry import build_region_geometry, build_region_statics
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )
    from rigel.calibration.structural_claims import build_structural_claims
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

    index = toy.index
    scan = dataclasses.replace(
        PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(str(toy.bam_path))
    )
    _stats, strand_model, _buf, payload = scan_and_buffer(str(toy.bam_path), index, scan)
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    geometry = build_region_geometry(chain, substrate, ra, sj, fl.gdna_pmf, fl.rna_pmf, None)
    statics = build_region_statics(chain, ra, bflags)
    claims = build_structural_claims(chain, statics)
    routes = RA.build_route_table(sj, substrate, fl.rna_pmf)
    kwargs = dict(
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        sj=sj,
        boundary_flags=bflags,
    )
    return dict(
        payload=payload,
        chain=chain,
        statics=statics,
        geometry=geometry,
        region_arrays=ra,
        claims=claims,
        routes=routes,
        calibrate_kwargs=kwargs,
    )


def test_the_factor_touches_only_the_claimed_populations(anchor_toy, monkeypatch):
    """No leakage, with the factor forced live (the toy is below the estimator minimums and
    gDNA-heavy, so the honest factor is near-flat there — a flat factor cannot leak, so the
    estimators are pinned to known values; TRAPS: could-the-arm-have-fired)."""
    t = anchor_toy
    monkeypatch.setattr(RA, "left_fit_center_spread", lambda o, p: (0.0, 0.01))
    monkeypatch.setattr(RA, "route_pair_log_variance", lambda a, b: 0.01)
    out = RA.build_rna_anchor_factor(
        t["chain"],
        t["statics"],
        t["geometry"],
        t["region_arrays"],
        t["routes"],
        n_grid=K,
        logodds_window=WINDOW,
    )
    touched = np.ptp(out, axis=1) > 0.0
    assert np.any(touched), "eligible slots must go live under pinned estimators"
    eligible = RA.eligible_slots(t["chain"], t["statics"], t["geometry"], t["region_arrays"])
    assert not np.any(touched & ~eligible), "the factor moved a slot outside its populations"


def test_the_anchor_is_a_message_silent_stays_silent(anchor_toy, monkeypatch):
    """THE CITIZENSHIP CONTRACT (owner ruling, 2026-08-25): the anchor is a MESSAGE. With message
    propagation OFF (the silent control) the flag must be inert — the row arithmetic never runs
    and the result is byte-identical to `rna_anchor=False`. The control is a control again."""
    import sys

    import rigel.calibration.calibrate  # noqa: F401

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    from rigel.config import CalibrationConfig

    class Boom(Exception):
        pass

    def boom(*a, **kw):
        raise Boom

    t = anchor_toy
    res_off = calibrate_mod.calibrate(
        payload=t["payload"],
        config=CalibrationConfig(rna_anchor=False, message_propagation=False),
        **t["calibrate_kwargs"],
    )
    monkeypatch.setattr(calibrate_mod, "prepare_flux_evidence", boom)
    res_on = calibrate_mod.calibrate(
        payload=t["payload"],
        config=CalibrationConfig(rna_anchor=True, message_propagation=False),
        **t["calibrate_kwargs"],
    )
    np.testing.assert_array_equal(
        np.asarray(res_on.mass_gdna_region), np.asarray(res_off.mass_gdna_region)
    )
    np.testing.assert_array_equal(
        np.asarray(res_on.mass_gdna_boundary), np.asarray(res_off.mass_gdna_boundary)
    )


def test_the_relay_carries_the_flux_stream_both_ways(anchor_toy, monkeypatch):
    """With the relay, `rna_anchor=True` must reach the evidence builder (explosive stub), and
    `rna_anchor=False` must not — the flag gates the STREAM, inside the policy path."""
    import sys

    import rigel.calibration.calibrate  # noqa: F401

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    from rigel.config import CalibrationConfig

    class Boom(Exception):
        pass

    def boom(*a, **kw):
        raise Boom

    monkeypatch.setattr(calibrate_mod, "prepare_flux_evidence", boom)
    t = anchor_toy
    calibrate_mod.calibrate(
        payload=t["payload"],
        config=CalibrationConfig(rna_anchor=False, message_propagation=True),
        **t["calibrate_kwargs"],
    )
    with pytest.raises(Boom):
        calibrate_mod.calibrate(
            payload=t["payload"],
            config=CalibrationConfig(rna_anchor=True, message_propagation=True),
            **t["calibrate_kwargs"],
        )


def test_the_stream_rows_equal_the_reference_arithmetic(anchor_toy, monkeypatch):
    """ARITHMETIC PARITY: the rows the relay's certified-flux stream delivers must be byte-equal
    to `build_rna_anchor_factor` — the citizenship moved, the arithmetic did not. A drift here is
    a bug in one of them, full stop. The estimators are PINNED so the rows are LIVE on this toy
    (below the minimums the honest rows are all-zero and any corruption of the policy path would
    multiply zeros — TRAPS: could-the-arm-have-fired)."""
    from rigel.calibration.messages.relay import RelayPolicy

    monkeypatch.setattr(RA, "left_fit_center_spread", lambda o, p: (0.0, 0.01))
    monkeypatch.setattr(RA, "route_pair_log_variance", lambda a, b: 0.01)
    t = anchor_toy
    ev = RA.prepare_flux_evidence(
        t["chain"], t["statics"], t["geometry"], t["region_arrays"], t["routes"]
    )
    got = RelayPolicy(flux=ev)._rows_at(K, WINDOW)
    assert got is not None and float(np.ptp(got)) > 0.0, "the pinned rows must be live"
    want = RA.build_rna_anchor_factor(
        t["chain"],
        t["statics"],
        t["geometry"],
        t["region_arrays"],
        t["routes"],
        n_grid=K,
        logodds_window=WINDOW,
    )
    np.testing.assert_array_equal(got, want)


def test_a_flux_message_is_not_silent_and_silent_carries_none():
    """`PsiMessage.lam_rows` is a channel: present ⇒ not silent; `silent()` ⇒ None."""
    from rigel.calibration.messages import PsiMessage

    assert PsiMessage.silent().lam_rows is None
    assert PsiMessage.silent().is_silent
    assert not PsiMessage(lam_rows=np.zeros((3, 5))).is_silent


def test_the_real_relay_delivers_the_stream_through_calibrate(anchor_toy, monkeypatch):
    """The REAL integration: `calibrate` with the relay and the flag ON must complete through the
    genuine path (no stubs), and the stream must REACH ψ — proven by pinning the recipient
    arithmetic to a biased claim and watching the result move (on this toy the honest rows are
    exactly zero — below the estimator minimums — so completion alone cannot distinguish a wired
    stream from a dead one; TRAPS: could-the-arm-have-fired)."""
    import sys

    import rigel.calibration.calibrate  # noqa: F401
    import rigel.calibration.rna_anchor  # noqa: F401

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    anchor_mod = sys.modules["rigel.calibration.rna_anchor"]
    from rigel.config import CalibrationConfig

    t = anchor_toy
    cfg = CalibrationConfig(rna_anchor=True, message_propagation=True)
    res_real = calibrate_mod.calibrate(payload=t["payload"], config=cfg, **t["calibrate_kwargs"])
    assert res_real is not None

    def biased_rows(evidence, *, n_grid, logodds_window):
        return np.tile(np.linspace(0.0, 60.0, int(n_grid)), (int(evidence["n_slots"]), 1))

    monkeypatch.setattr(anchor_mod, "flux_rows", biased_rows)
    res_biased = calibrate_mod.calibrate(payload=t["payload"], config=cfg, **t["calibrate_kwargs"])
    real = np.asarray(res_real.mass_gdna_region, np.float64)
    biased = np.asarray(res_biased.mass_gdna_region, np.float64)
    assert not np.array_equal(real, biased), (
        "a pinned biased claim must move the result — identical outputs mean the stream never "
        "reached psi"
    )


def test_the_flag_exists_and_defaults_on():
    from rigel.config import CalibrationConfig

    assert CalibrationConfig().rna_anchor is True
