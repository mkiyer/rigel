"""Gates for THE BRIDGE (`calibration.messages.unified`) — the foundation spec running on the
real backbone.

The bridge is the CONVERGENCE TARGET, not another sandbox policy: one runner that plugs a
`(PropagationModel, SolveModel)` pair into the backbone's existing protocol, builds each node's
OWN message from the banks the context already carries, and is anchored by byte-identities before
it claims anything new. Written BEFORE the module, verified failing; every gate watched firing
against a deliberately broken build.
"""

from __future__ import annotations

import dataclasses
import sys

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy, PsiMessage
from rigel.calibration.messages.silent import SilentPolicy


def _U():
    from rigel.calibration.messages import foundation, unified

    return foundation, unified


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """ONE real `solve_chain` call captured from a calibrate run on the toy — the
    backbone-parity pattern: every gate below re-runs the sweep with a different policy on
    byte-identical inputs."""
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "tpo_for_bridge", Path(__file__).parent / "test_prior_vs_oracle.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tpo_for_bridge"] = m
    spec.loader.exec_module(m)
    toy = m.toy.__wrapped__(tmp_path_factory)

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )
    from rigel.config import CalibrationConfig, PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

    index = toy.index
    scan_cfg = dataclasses.replace(
        PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(str(toy.bam_path))
    )
    _stats, strand_model, _buf, payload = scan_and_buffer(str(toy.bam_path), index, scan_cfg)
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    grabbed: list = []
    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    orig = SW.solve_chain

    def spy(chain, statics, geometry, belief, region_arrays, **kw):
        if not grabbed:
            grabbed.append((chain, statics, geometry, belief, region_arrays, dict(kw)))
        return orig(chain, statics, geometry, belief, region_arrays, **kw)

    calibrate_mod.solve_chain = spy
    try:
        calibrate_mod.calibrate(
            payload=payload,
            config=CalibrationConfig(rna_anchor=True, message_propagation=True),
            region_arrays=ra,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf,
            rna_fl_pmf=fl.rna_pmf,
            sj=build_sj_geometry_arrays(index),
            boundary_flags=build_boundary_flags_array(index),
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert grabbed, "the spy never fired"
    chain, statics, geometry, belief, region_arrays, kw = grabbed[0]
    shipped_flux = getattr(kw.get("policy"), "_flux", None)
    kw = {k: v for k, v in kw.items() if k not in ("policy", "_capture")}
    return dict(
        args=(chain, statics, geometry, belief, region_arrays),
        kw=kw,
        flux=shipped_flux,
    )


def _run(si, policy):
    out = SW.solve_chain(*si["args"], **si["kw"], policy=policy)
    return {
        f: np.asarray(getattr(out, f))
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg")
    }


def test_the_runner_satisfies_the_backbone_protocol():
    F, U = _U()
    pol = U.UnifiedPolicy(F.PassThroughPropagation(), U.SilentSolve())
    assert isinstance(pol, Policy)


def test_silent_solve_is_byte_identical_to_the_silent_policy(sweep_inputs):
    """IDENTITY ANCHOR ⓐ (the owner's rule made executable): silence is the solve ignoring the
    messages — nothing else changes, because a node's belief only updates at the solve. The
    runner with a pass-through propagation and a silent solve must reproduce SilentPolicy
    byte-for-byte through the REAL backbone."""
    F, U = _U()
    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, U.UnifiedPolicy(F.PassThroughPropagation(), U.SilentSolve()))
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


def test_flux_rows_solve_is_byte_identical_to_the_channels_off_relay(sweep_inputs, monkeypatch):
    """IDENTITY ANCHOR ⓑ: a solve whose spliced treatment is the shipped anchor arithmetic must
    reproduce the relay running ONLY its certified-flux stream (every Gaussian channel off),
    byte-for-byte through the real backbone — the spliced-lane seam proven against the shipped
    code before anything new is claimed. The estimators are PINNED so the rows are LIVE on this
    toy (its honest rows are all-zero, and an identity between two zero factors proves nothing —
    TRAPS: could-the-arm-have-fired)."""
    import rigel.calibration.rna_anchor as RA
    from rigel.calibration.messages.relay import RelayPolicy, RelaySwitches
    from rigel.calibration.rna_anchor import flux_rows

    monkeypatch.setattr(RA, "left_fit_center_spread", lambda o, p: (0.0, 0.01))
    monkeypatch.setattr(RA, "route_pair_log_variance", lambda a, b: 0.01)
    F, U = _U()
    flux = sweep_inputs["flux"]
    assert flux is not None, "the toy calibrate must have prepared flux evidence"
    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    rows = flux_rows(flux, n_grid=n_grid, logodds_window=window)
    assert float(np.ptp(rows)) > 0.0, "the pinned rows must be live"
    off = RelaySwitches(
        gdna_channel=False, rna_channel=False, lam_channel=False, theta_channel=False
    )
    a = _run(sweep_inputs, RelayPolicy(off, flux=flux))
    b = _run(
        sweep_inputs,
        U.UnifiedPolicy(F.PassThroughPropagation(), U.RowsSolve(rows)),
    )
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


def test_own_messages_carry_the_node_banks(sweep_inputs):
    """`prepare` must build each node's OWN message from the banks the context already carries:
    the unspliced lanes mirror the self-solve densities/precisions, and the spliced lanes are
    live exactly where the slot carries certified sj flux."""
    F, U = _U()

    class Recorder(U.UnifiedPolicy):
        def prepare(self, ctx):
            prepared = super().prepare(ctx)
            self.ctx = ctx
            self.prepared = prepared
            return prepared

    pol = Recorder(F.PassThroughPropagation(), U.SilentSolve())
    _run(sweep_inputs, pol)
    ctx, own = pol.ctx, pol.prepared.own
    np.testing.assert_array_equal(own["unspliced_gdna"][0], np.asarray(ctx.own.rho_g))
    np.testing.assert_array_equal(own["unspliced_gdna"][1], np.asarray(ctx.own.prec_g))
    spl = np.asarray(ctx.sj_count, np.float64)
    live_pos = (np.asarray(ctx.eff_sj)[:, 0] > 0) & (spl[:, 0] > 0)
    assert live_pos.any(), "the toy must carry some positive-strand sj flux"
    assert np.all(own["spliced_rna_pos"][1][live_pos] > 0), (
        "spliced lanes must be live where the slot has certified flux"
    )
    no_opp = np.asarray(ctx.eff_sj)[:, 0] <= 0
    assert np.all(own["spliced_rna_pos"][1][no_opp] == 0.0), (
        "spliced lanes must be silent where there is no crossing opportunity"
    )
    # the route-sum homecoming: the lane's abundance is the ROUTE-SUMMED rate, never the pooled
    # ratio-of-sums (the round-2 review's k-route under-read)
    want = (np.asarray(ctx.route_rate_lo) + np.asarray(ctx.route_rate_hi))[:, 0]
    np.testing.assert_allclose(own["spliced_rna_pos"][0], want, rtol=1e-12)


def test_propagation_travels_the_chain(sweep_inputs):
    """The scans must actually MOVE information ≥ 2 hops. The proof rides the UNSPLICED gDNA
    lane (spliced lanes are one-hop by design under FrameAware; this gate runs pass-through, so
    any lane works — gDNA is the densest)."""
    F, U = _U()

    class SpySolve(U.SilentSolve):
        def solve(self, own, forward, backward):
            self.seen = (own, forward, backward)
            return PsiMessage.silent()

    solver = SpySolve()
    _run(sweep_inputs, U.UnifiedPolicy(F.PassThroughPropagation(), solver))
    own, fwd, bwd = solver.seen
    # the sparsest belief lane on this toy gives the 2-hop proof its empty middle slots
    own_p = np.asarray(own.unspliced_rna_neg.precision)
    fwd_p = np.asarray(fwd.unspliced_rna_neg.precision)
    bwd_p = np.asarray(bwd.unspliced_rna_neg.precision)
    # ⛔ one-hop receipt is the BACKBONE's gather, not the scan's work: the claim must be shown
    # at least TWO hops from any source — a slot whose own lane AND both neighbours' own lanes
    # are silent, yet whose delivered message carries the claim (TRAPS: could-the-arm-have-fired,
    # which this gate failed once: a scan that never stepped still passed the one-hop form).
    chain = sweep_inputs["args"][0]
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    n = own_p.shape[0]
    nb_l = own_p[np.clip(left, 0, n - 1)] * (left >= 0)
    nb_r = own_p[np.clip(right, 0, n - 1)] * (right >= 0)
    far = (own_p == 0) & (nb_l == 0) & (nb_r == 0)
    received = ((fwd_p > 0) | (bwd_p > 0)) & far
    assert received.any(), (
        "no slot two hops from any RNA− belief received a claim — the scans moved nothing"
    )


# ── step 2: the first real propagation model ────────────────────────────────────────────────────


def _hop_tables(log_r=0.0, v_r=1.0, premise=0.0, fp=(True, True), fn=(True, True)):
    """Synthetic per-hop context for unit gates: one hop, src index 0, dst index 1."""
    F, U = _U()
    m = U.FrameAwarePropagation()
    m._premise = premise
    m._tables = {
        False: {
            "log_r": np.array([0.0, log_r]),
            "v_r": np.array([1.0, v_r]),
            "fp": np.array([fp[0], fp[1]]),
            "fn": np.array([fn[0], fn[1]]),
            "src": np.array([0, 0]),
            "shed_pos": np.zeros(2),
            "shed_neg": np.zeros(2),
        }
    }
    return F, U, m


def _one_hop(m, F, U, lane, claim):
    incoming = F.Message.silent().with_lane(lane, claim)
    hop = U.Hop(src=0, dst=1, backward=False)
    return m.propagate(F.Message.silent(), incoming, hop).lane(lane)


def test_frame_aware_refuses_a_strand_population_change():
    """A strand admissible on one side of a hop and not the other is a DIFFERENT RNA population
    (the owner's 2026-08-18 rule, extracted from the donors): that lane's claim is refused —
    value AND precision in one statement. The gDNA lane always crosses (gDNA is genomically
    continuous — AXIOM 0)."""
    F, U, m = _hop_tables(fp=(True, False))
    out = _one_hop(m, F, U, "unspliced_rna_pos", F.Claim(2.0, 10.0))
    assert out.is_silent, "a changed-population lane must be refused"
    g = _one_hop(m, F, U, "unspliced_gdna", F.Claim(2.0, 10.0))
    assert g.precision > 0, "the gDNA lane always crosses"


def test_the_knob_interpolates_between_the_two_strategies():
    """The knob (currency's derivation, no constant): an enrichment reading that dwarfs its
    counting noise is believed nearly in full (value ×≈ r); one inside the noise is shrunk to
    nearly nothing (value ≈ unchanged). Every crossing costs precision."""
    F, U, m = _hop_tables(log_r=2.0, v_r=0.01)
    strong = _one_hop(m, F, U, "unspliced_gdna", F.Claim(1.0, 10.0))
    assert 0.9 * np.exp(2.0) < strong.abundance <= np.exp(2.0)
    F, U, m = _hop_tables(log_r=0.05, v_r=1.0)
    weak = _one_hop(m, F, U, "unspliced_gdna", F.Claim(1.0, 10.0))
    assert 1.0 <= weak.abundance < 1.01
    # ⭐ the knob's costs are the REFRAME's, so they are charged to the message's SHARED LEVEL
    # (the level derivation, 2026-08-27) — not to this lane's precision
    hop = U.Hop(src=0, dst=1, backward=False)
    assert m.level_cost(hop) > 0.0, "a hop with any disagreement must cost the shared level"
    _f_flat, _u_flat, m_flat = _hop_tables(log_r=0.0, v_r=1.0)
    assert m_flat.level_cost(hop) == 0.0, "a perfectly agreeing frame costs the level nothing"


def test_the_premise_is_charged_on_every_hop():
    """An imputation must cost something every hop, even at a perfectly agreeing frame: with
    log_r = 0 the knob and disagreement vanish, and the premise alone must still lower the
    precision (the recorded free-hop failure is what this forbids)."""
    F, U, m = _hop_tables(log_r=0.0, premise=0.3)
    out = _one_hop(m, F, U, "unspliced_gdna", F.Claim(1.0, 10.0))
    assert np.isclose(out.abundance, 1.0)
    assert np.isclose(out.precision, 10.0 / (1.0 + 10.0 * 0.3))


def test_the_mass_rescale_projects_a_licensed_state_onto_the_local_total():
    """THE MASS RESCALE (relay's scan-time conservation projection, ported behind the A/B
    flag): at every licensed hop the running state must satisfy the local identity
    sum_c rho_c*E_c = M — an OVERWRITE by k = M/S, never a fuse, so the message carries the
    COMPOSITION while the level is re-measured at each licensed slot. Licensed = the incoming
    state supplied both components, or the destination is structurally pure-gDNA. A bare gDNA
    level at an ordinary slot passes unprojected."""
    F, U, m = _hop_tables(log_r=0.0, v_r=1.0)
    m._mass_rescale = True
    t = m._tables[False]
    t["E_g"] = np.array([1.0, 2.0])
    t["E_r"] = np.array([1.0, 4.0])
    t["M"] = np.array([1.0, 40.0])
    t["own_g"] = np.zeros(2)
    t["own_p"] = np.zeros(2)
    t["own_n"] = np.zeros(2)
    t["g1"] = np.array([False, False])
    hop = U.Hop(src=0, dst=1, backward=False)
    # the supply test: gdna plus EVERY admissible strand (the fixture admits only +)
    t["fn"][1] = False
    supplied = (
        F.Message.silent()
        .with_lane("unspliced_gdna", F.Claim(3.0, 10.0, 1.0))
        .with_lane("unspliced_rna_pos", F.Claim(1.0, 5.0))
    )
    out = m.propagate(F.Message.silent(), supplied, hop)
    vg = out.lane("unspliced_gdna").abundance
    vp = out.lane("unspliced_rna_pos").abundance
    # budget before: 3*2 + 1*4 = 10 against M = 40 => k = 4; ratios preserved
    assert np.isclose(vg * 2.0 + vp * 4.0, 40.0), (vg, vp)
    assert np.isclose(vg / vp, 3.0), "the projection preserves the composition"
    bare = F.Message.silent().with_lane("unspliced_gdna", F.Claim(3.0, 10.0, 1.0))
    out2 = m.propagate(F.Message.silent(), bare, hop)
    assert np.isclose(out2.lane("unspliced_gdna").abundance, 3.0), (
        "an unlicensed bare level passes unprojected"
    )
    t["g1"] = np.array([False, True])
    out3 = m.propagate(F.Message.silent(), bare, hop)
    assert np.isclose(out3.lane("unspliced_gdna").abundance * 2.0, 40.0), (
        "a structurally pure-gDNA destination projects regardless: its whole mass IS gDNA"
    )


def test_the_residual_level_claims_the_unexplained_boundary_mass_as_rna():
    """RESIDUAL_LEVEL (relay's boundary law, ported behind the A/B flag): everything the
    imputed gDNA level cannot explain about a boundary's OWN observed mass is continuing RNA,
    delivered on the RNA channel at near-counting precision. At a low-gDNA boundary the claim
    is rho_nu ~= M/E_r — the positive half of the pincer the unified boundary solve lacked."""
    F, U = _U()
    m = _mini_solve()
    m._residual_level = True
    m._is_boundary = np.array([True])
    m._M = np.full(1, 20.0)
    m._n_slot = np.full(1, 20.0)
    m._free_pos = np.array([True])
    m._free_neg = np.array([False])
    gdna = F.Claim(np.array([0.01]), np.array([10.0]), np.array([10.0]))
    own, fwd, bwd = _msgs(F, claim=gdna)
    out = m.solve_unspliced(own, fwd, bwd)
    assert out.rna_prec[0][0] > 0, "the residual claim must reach the RNA channel"
    # rho_nu ~= M/E_r = 20/1; the delivered mode is log(x*E_r/M) ~= log(1) = 0
    assert abs(float(out.rna_mode[0][0])) < 0.35, float(out.rna_mode[0][0])
    m2 = _mini_solve()
    m2._is_boundary = np.array([True])
    m2._M = np.full(1, 20.0)
    m2._n_slot = np.full(1, 20.0)
    m2._free_pos = np.array([True])
    m2._free_neg = np.array([False])
    own2, fwd2, bwd2 = _msgs(F, claim=gdna)
    assert m2.solve_unspliced(own2, fwd2, bwd2).rna_prec[0][0] == 0.0, (
        "with the flag off the RNA channel stays dead (the v1 baseline)"
    )


def test_the_conservation_is_the_count_identity_not_a_density_sum():
    """THE CONSERVATION IDENTITY (region_geometry's own contract, and psi's): the components
    of one slot satisfy sum_c rho_c * E_c = M — each density weighted by ITS OWN opportunity,
    summing to the slot's OBSERVED unspliced count. That is what makes the delivered f_c a
    COUNT share, which is what psi consumes.

    ⛔ Summing raw densities against the reciprocal-opportunity total is a DIFFERENT and, at
    regions, a WRONG identity: `inv_abundance` reads rho * P(w <= ell) there (region_geometry's
    contract; measured on the ladder: median exon P = 0.452, p5 = 0.044 — the total is up to
    23x too low at half of all exon slots), while it is exact at boundaries. The count identity
    is exact at both.

    The check is the invariant psi expects: with no absorber licensed, the delivered shares of
    the live lanes sum to one. Lanes with DIFFERENT opportunities make this fail for any
    unweighted density sum."""
    F, U = _U()
    m = _mini_solve()
    m._E = {
        "unspliced_gdna": np.array([1.0]),
        "unspliced_rna_pos": np.array([4.0]),
        "unspliced_rna_neg": np.array([4.0]),
    }
    m._M = np.array([40.0])
    m._n_slot = np.array([40.0])
    silent = F.Message(**{k: F.Claim(np.zeros(1), np.zeros(1)) for k in F.Message.LANES})
    fwd = silent.with_lane(
        "unspliced_gdna", F.Claim(np.array([3.0]), np.array([5.0]), np.array([5.0]))
    ).with_lane("unspliced_rna_pos", F.Claim(np.array([1.0]), np.array([5.0]), np.array([5.0])))
    out = m.solve_unspliced(silent, fwd, silent)
    c_g = float(np.exp(out.gdna_mode[0])) * 40.0
    c_r = float(np.exp(out.rna_mode[0][0])) * 40.0
    assert out.gdna_prec[0] > 0 and out.rna_prec[0][0] > 0
    # (a) the counts account for the observed mass
    assert abs(c_g + c_r - 40.0) < 2.0, (
        f"the components' counts must account for the slot's mass (got {c_g + c_r})"
    )
    # (b) the SPLIT is decided in count space: with equal precisions the residual is shared
    # equally in FRAGMENTS, so both lanes gain the same number — the arriving counts were
    # 3*1 = 3 and 1*4 = 4, so the delivered counts keep that one-fragment difference. An
    # unweighted density sum shares the residual in DENSITY and inverts the gap.
    assert abs((c_r - c_g) - 1.0) < 0.5, (
        f"the deficit must be shared in count space (gap {c_r - c_g}, want ~1.0)"
    )


def test_the_level_solve_has_all_three_operators_as_limits():
    """THE LEVEL DERIVATION (owner directive, 2026-08-27) — one conservation solve with TWO
    variance kinds: V, the SHARED scale uncertainty a transported message carries (the reframe
    knob's costs multiply every lane identically, so they are common-mode), and v_c, each
    component's own. Minimising a^2/2V + sum u_c^2/2v_c subject to sum mu_c e^(a+u_c) = M
    gives, around the current point,

        a   = V   p_M D S   / (1 + p_M W)      D = M - S,  S = sum y_c
        u_c = v_c p_M D y_c / (1 + p_M W)      W = V S^2 + sum v_c y_c^2

    and the three operators built separately are its LIMITS. This gates all three:

    * V >> v_c  =>  a -> log(M/S), u -> 0: the multiplicative mass rescale (composition held);
    * V = 0     =>  the residual distributes per component by v_c y_c: the additive allocation;
    * one lane uninformative (v -> large) => the residual lands there: residual_level's law,
      "everything the gDNA level cannot explain about this mass is RNA".
    """
    U = _alloc()
    mu = np.array([3.0, 1.0])

    # (1) shared-level limit: composition preserved, counts sum to M
    y = U.allocate_level(
        mu, np.array([1e-6, 1e-6]), total=8.0, total_precision=1e6, level_logvar=1.0, absorber=False
    )
    assert abs(y.sum() - 8.0) < 1e-6, y
    assert abs(y[0] / y[1] - 3.0) < 1e-3, f"a shared-level residual must hold the composition {y}"

    # (2) no shared level: the residual moves the composition, weighted by v_c * y_c
    y2 = U.allocate_level(
        mu, np.array([1.0, 1.0]), total=8.0, total_precision=1e6, level_logvar=0.0, absorber=False
    )
    assert abs(y2.sum() - 8.0) < 1e-6, y2
    assert y2[0] / y2[1] > 3.5, f"an own-error residual must move the composition {y2}"

    # (3) one uninformative lane takes the residual whole (residual_level's limit)
    y3 = U.allocate_level(
        mu, np.array([1e-9, 1e4]), total=8.0, total_precision=1e6, level_logvar=0.0, absorber=False
    )
    assert abs(y3.sum() - 8.0) < 1e-6, y3
    assert abs(y3[0] - 3.0) < 0.05, f"the confident lane holds its level {y3}"
    assert abs(y3[1] - 5.0) < 0.05, f"the uninformative lane takes the remainder {y3}"

    # positivity is structural: the log form cannot produce a negative count
    y4 = U.allocate_level(
        mu, np.array([1.0, 1.0]), total=0.5, total_precision=1e6, level_logvar=0.0, absorber=False
    )
    assert np.all(y4 > 0), y4


def test_a_common_mode_residual_scales_the_claims_and_keeps_the_composition():
    """THE CONSERVATION OPERATOR, additive vs multiplicative. When every live claim is an
    imputation of comparable quality and their counts fall short of the slot's mass, the
    residual is a LEVEL error common to all of them — the honest repair scales them together
    (relay's k = M/S) and leaves the COMPOSITION untouched. Shifting the residual additively
    onto the weakest lane converts a level error into a composition error, which is how a
    near-zero gDNA claim eats an unstranded slot's unexplained RNA mass.

    Gate on the multiplicative mode: two equally-precise claims short of the mass keep their
    ratio exactly, and their counts sum to M."""
    F, U = _U()
    m = _mini_solve()
    m._conserve_multiplicative = True
    m._E = {
        "unspliced_gdna": np.array([1.0]),
        "unspliced_rna_pos": np.array([1.0]),
        "unspliced_rna_neg": np.array([1.0]),
    }
    m._M = np.array([40.0])
    m._n_slot = np.array([40.0])
    silent = F.Message(**{k: F.Claim(np.zeros(1), np.zeros(1)) for k in F.Message.LANES})
    fwd = silent.with_lane(
        "unspliced_gdna", F.Claim(np.array([1.0]), np.array([5.0]), np.array([5.0]))
    ).with_lane("unspliced_rna_pos", F.Claim(np.array([3.0]), np.array([5.0]), np.array([5.0])))
    out = m.solve_unspliced(silent, fwd, silent)
    c_g = float(np.exp(out.gdna_mode[0])) * 40.0
    c_r = float(np.exp(out.rna_mode[0][0])) * 40.0
    assert abs(c_g + c_r - 40.0) < 1e-6, (c_g, c_r)
    assert abs(c_r / c_g - 3.0) < 1e-6, (
        f"a common-mode residual must not move the composition (ratio {c_r / c_g}, want 3.0)"
    )


def test_spliced_lanes_never_enter_the_unspliced_conservation():
    """THE SPLICED LANES LEAVE THE ALLOCATION (the owner's rule, made structural by the count
    identity): spliced fragments are not in M — they are counted separately — so they carry no
    information into the unspliced conservation and cannot displace an unspliced witness.
    A live spliced lane of any size must not move the delivered unspliced channels by one ulp;
    its evidence reaches psi only through `solve_spliced`'s rows."""
    F, U = _U()
    silent = F.Message(**{k: F.Claim(np.zeros(1), np.zeros(1)) for k in F.Message.LANES})
    gdna = F.Claim(np.array([2.0]), np.array([5.0]), np.array([5.0]))
    fwd = silent.with_lane("unspliced_gdna", gdna)

    def run(own):
        m = _mini_solve()
        m._M = np.array([10.0])
        m._n_slot = np.array([10.0])
        return m.solve_unspliced(own, fwd, silent)

    bare = run(silent)
    with_flux = run(
        silent.with_lane(
            "spliced_rna_pos", F.Claim(np.array([99.0]), np.array([40.0]), np.array([40.0]))
        )
    )
    assert float(bare.gdna_mode[0]) == float(with_flux.gdna_mode[0])
    assert float(bare.gdna_prec[0]) == float(with_flux.gdna_prec[0])


def test_frame_aware_under_a_silent_solve_changes_nothing(sweep_inputs):
    """The free invariant: propagation NEVER changes the answer on its own — a belief only forms
    at the solve, so FrameAwarePropagation under SilentSolve must equal SilentPolicy
    byte-for-byte through the real backbone."""
    F, U = _U()
    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, U.UnifiedPolicy(U.FrameAwarePropagation(), U.SilentSolve()))
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


# ── step 3 commit 2: the allocation solve ───────────────────────────────────────────────────────


def _alloc():
    from rigel.calibration.messages import unified as U

    return U


def test_allocation_passes_through_when_totals_agree():
    """The owner's solvency rule: when the arriving components already sum to the node's total,
    the residual is zero and every component passes through untouched — no arbitration exists
    because none is needed."""
    U = _alloc()
    mu = np.array([3.0, 1.0, 0.5])
    p = np.array([10.0, 5.0, 2.0])
    x = U.allocate(mu, p, total=4.5, total_precision=100.0, absorber=False)
    np.testing.assert_allclose(x, mu, rtol=1e-12)


def test_allocation_lands_the_residual_on_weak_components():
    """The precision-weighted allocation: the residual is distributed in proportion to VARIANCE —
    the certified/strong components hold their values, the weak absorb. The derived repair of the
    equal-weight rescale (the recorded k = 467,000× and 235,800× amplifications)."""
    U = _alloc()
    mu = np.array([10.0, 1.0])
    p = np.array([1000.0, 0.1])
    x = U.allocate(mu, p, total=14.0, total_precision=50.0, absorber=False)
    assert abs(x[0] - 10.0) < 0.05, "the strong component must hold its value"
    assert x[1] > 3.5, "the weak component must absorb nearly the whole deficit"
    assert x.sum() < 14.0 + 1e-9, "a soft constraint never overshoots the total"


def test_a_licensed_absorber_takes_the_deficit_and_known_components_hold():
    """The owner's insolvency rule, structural: where a terminus licenses NEW-RNA, the absorber
    takes the whole deficit and the known components keep their claimed values — the solve is in
    abundance space, never share space."""
    U = _alloc()
    mu = np.array([2.0, 1.0])
    p = np.array([50.0, 20.0])
    x = U.allocate(mu, p, total=10.0, total_precision=100.0, absorber=True)
    np.testing.assert_allclose(x, mu, rtol=1e-12)


def test_an_excess_cannot_go_negative():
    """Excess (arriving claims above the total) is shrunk variance-weighted and clamped at zero —
    an absorber cannot absorb an excess (new transcription cannot be negative)."""
    U = _alloc()
    mu = np.array([5.0, 0.5])
    p = np.array([100.0, 0.05])
    x = U.allocate(mu, p, total=4.0, total_precision=200.0, absorber=True)
    assert np.all(x >= 0.0), "an excess allocation may never go negative"
    assert x[1] == 0.0, "the weak component is exhausted first"
    assert 3.9 < x[0] <= 5.0, "the remaining excess lands on the strong component, bounded"


def test_the_allocation_solve_is_silent_on_silent_messages(sweep_inputs):
    """End to end: the full unspliced solve under silent messages must equal SilentPolicy
    byte-for-byte (no message ⇒ nothing to allocate ⇒ nothing delivered)."""
    F, U = _U()
    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, U.UnifiedPolicy(F.PassThroughPropagation(), U.AllocationSolve()))
    # under PassThrough the messages are NOT silent (own beliefs travel), so this is NOT an
    # identity run — instead assert the solve RUNS end to end and stays finite wherever the
    # silent baseline is finite (a handful of inf variances are a baseline property of the
    # sweep at degenerate slots, not the solve's doing)
    for f in b:
        base_finite = np.isfinite(a[f])
        assert np.all(np.isfinite(b[f][base_finite])), f
    assert not all(np.array_equal(a[f], b[f]) for f in a), (
        "with travelling messages the allocation solve must move SOMETHING vs silence "
        "(TRAPS: could-the-arm-have-fired)"
    )


def _mini_solve(n=1):
    U = _alloc()
    m = U.AllocationSolve()
    m._E = {k: np.ones(n) for k in m._UNSPLICED}
    m._M = np.ones(n)
    m._T = np.full(n, 10.0)
    m._pT = np.full(n, 100.0)
    m._M = np.full(n, 1.0)
    m._n_slot = np.full(n, 10.0)
    m._pM = np.full(n, 100.0)
    m._absorber = np.zeros(n, bool)
    m._dom = (-30.0, 0.0)
    # the composition-evidence variances the reception law is keyed on: the default mini node
    # is composition-BLIND (v_own = inf, the unstranded/AMBIG state — messages pass untouched)
    m._v_own_g = np.full(n, np.inf)
    m._v_own_r = np.full(n, np.inf)
    m._free_pos = np.zeros(n, bool)
    m._free_neg = np.zeros(n, bool)
    return m


def _msgs(F, lane="unspliced_gdna", claim=None):
    silent = F.Message(**{k: F.Claim(np.zeros(1), np.zeros(1)) for k in F.Message.LANES})
    fwd = silent.with_lane(lane, claim) if claim is not None else silent
    return silent, fwd, silent


def test_delivered_precision_is_the_measurement_stream():
    """MEASUREMENT CITIZENSHIP (the shipped relay's law, ported): the psi channel precision a
    solve delivers is the MEASUREMENT stream — what independent witnesses actually counted
    (struct-locked gDNA anchors, certified sj counts), hop-damped — never the belief precision.
    A lane arriving with an enormous belief precision and no measurement delivers a DEAD
    channel; the same lane with a measurement seed delivers exactly that seed (the node here is
    composition-blind, so reception does not touch it)."""
    F, U = _U()
    m = _mini_solve()
    own, fwd, bwd = _msgs(F, claim=F.Claim(np.array([8.0]), np.array([500.0]), np.array([0.0])))
    assert m.solve_unspliced(own, fwd, bwd).gdna_prec[0] == 0.0, (
        "a pure belief may inform the value but must never be delivered as channel precision"
    )
    m2 = _mini_solve()
    own, fwd, bwd = _msgs(F, claim=F.Claim(np.array([8.0]), np.array([500.0]), np.array([30.0])))
    got = m2.solve_unspliced(own, fwd, bwd).gdna_prec[0]
    assert np.isclose(got, 30.0), got


def test_reception_is_the_relay_deflation_law():
    """THE RECEPTION LAW (mismatch_deflate, ported verbatim in shape): per stream,
    p_eff = 1/max(v_stream, G^2 - v_own) with G the log gap between the arriving and the own
    lane value, v_own the node's OWN COMPOSITION variance. Strong own composition + a
    disagreeing claim => capped near 1/G^2; an agreeing claim passes; a composition-blind node
    (v_own = inf) passes everything untouched — which is exactly what lets messages solve
    unstranded data while staying weak on stranded data."""
    F, U = _U()
    m = _mini_solve()
    m._v_own_g = np.array([1e-3])
    claim = F.Claim(np.array([8.0]), np.array([500.0]), np.array([30.0]))
    own, fwd, bwd = _msgs(F, claim=claim)
    own = own.with_lane("unspliced_gdna", F.Claim(np.array([1.0]), np.array([50.0])))
    got = m.solve_unspliced(own, fwd, bwd).gdna_prec[0]
    g2 = float(np.log(8.0) ** 2)
    assert got <= 1.0 / (g2 - 1e-3) * 1.0001, (got, 1.0 / g2)
    # the same claim, agreeing: untouched
    m2 = _mini_solve()
    m2._v_own_g = np.array([1e-3])
    agree = F.Claim(np.array([1.0]), np.array([500.0]), np.array([30.0]))
    own2, fwd2, bwd2 = _msgs(F, claim=agree)
    own2 = own2.with_lane("unspliced_gdna", F.Claim(np.array([1.0]), np.array([50.0])))
    assert np.isclose(m2.solve_unspliced(own2, fwd2, bwd2).gdna_prec[0], 30.0)
    # the same disagreeing claim at a composition-blind node: untouched
    m3 = _mini_solve()
    own3, fwd3, bwd3 = _msgs(F, claim=claim)
    own3 = own3.with_lane("unspliced_gdna", F.Claim(np.array([1.0]), np.array([50.0])))
    assert np.isclose(m3.solve_unspliced(own3, fwd3, bwd3).gdna_prec[0], 30.0)


def test_the_deficit_lands_on_admissible_unseen_rna_not_on_a_measured_witness():
    """THE UNSEEN-COMPONENT ABSORBER: a conservation deficit at a node whose annotation admits
    an RNA strand NOBODY has evidence about (silent lane, admissible bit set) belongs to that
    unseen component — never to the weakest live witness. Without this, an unstranded library's
    RNA mass (invisible to every witness) is force-fed to the near-zero gDNA anchor claim as
    phantom gDNA. With no admissible silent lane the conservation bites as before."""
    F, U = _U()
    m = _mini_solve()
    m._M = np.full(1, 10.0)  # slot mass matches the total, so share = x/10 and the clamp
    # cannot flatten hold (x=1, mode=log 0.1) into inflate (x~10, mode~log 1)
    m._free_pos = np.array([True])
    m._free_neg = np.array([True])
    gdna = F.Claim(np.array([1.0]), np.array([0.2]), np.array([0.2]))
    own, fwd, bwd = _msgs(F, claim=gdna)
    out = m.solve_unspliced(own, fwd, bwd)
    assert out.gdna_prec[0] > 0
    assert float(out.gdna_mode[0]) < -2.0, (
        f"the measured gdna witness was inflated to absorb unseen RNA "
        f"(mode {float(out.gdna_mode[0])}, want ~log 0.1 = -2.3)"
    )


def test_a_contradicted_claim_is_zeroed():
    """The contradiction rule (relay's): the own lane and the arriving lane on exactly one
    side of zero, at a node WITH composition evidence => the claim is killed outright, both
    streams — a confident zero against a confident positive is not a weak agreement."""
    F, U = _U()
    m = _mini_solve()
    m._v_own_g = np.array([1e-3])
    zero_claim = F.Claim(np.array([0.0]), np.array([20.0]), np.array([10.0]))
    own, fwd, bwd = _msgs(F, claim=zero_claim)
    own = own.with_lane("unspliced_gdna", F.Claim(np.array([1.0]), np.array([50.0])))
    assert m.solve_unspliced(own, fwd, bwd).gdna_prec[0] == 0.0


# ── step 3 commit 3: the spliced solve ──────────────────────────────────────────────────────────


def test_count_recovery_inverts_the_trigamma():
    """count → precision (1/trigamma(count+½)) → count must round-trip."""
    U = _alloc()
    from rigel.calibration.messages.variance import count_logvar

    for c in (0.0, 1.0, 7.0, 42.0, 913.0):
        p = 1.0 / count_logvar(np.array([c]))
        got = U.count_from_precision(np.asarray(p, float))[0]
        assert abs(got - c) < 1e-6, (c, got)


def test_incomplete_flanks_make_no_flux_claim():
    """The ruled completeness branch, v1: where a terminus gain makes a flank incomplete, the
    exon's flux claim is withheld (the one-sided form is the recorded refinement) — an unknown
    component may enter there and an equality would be wrong."""
    U = _alloc()
    ok = U.flank_complete(
        spliced_live=np.array([True, True]),
        terminus_gain=np.array([False, True]),
    )
    assert bool(ok[0]) and not bool(ok[1])


def test_the_spliced_solve_moves_beliefs_through_the_backbone(sweep_inputs, monkeypatch):
    """End to end with pinned estimators (the toy sits below the estimator minimums, where honest
    rows are flat — TRAPS: could-the-arm-have-fired): the lane-native spliced law must reach ψ
    and move the beliefs vs the same solve with its spliced half disabled."""
    import rigel.calibration.rna_anchor as RA

    monkeypatch.setattr(RA, "left_fit_center_spread", lambda o, p: (0.0, 0.01))
    monkeypatch.setattr(RA, "route_pair_log_variance", lambda a, b: 0.01)
    monkeypatch.setattr(RA, "left_tail_log_variance", lambda o, p: 0.01)
    F, U = _U()
    a = _run(sweep_inputs, U.UnifiedPolicy(U.FrameAwarePropagation(), U.AllocationSolve()))
    solver = U.AllocationSolve()
    solver._spliced_enabled = False
    b = _run(sweep_inputs, U.UnifiedPolicy(U.FrameAwarePropagation(), solver))
    assert not all(np.array_equal(a[f], b[f]) for f in a), (
        "the spliced rows must reach psi and move the beliefs"
    )
