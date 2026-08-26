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


def test_propagation_travels_the_chain(sweep_inputs):
    """The scans must actually MOVE information: a slot whose own spliced lane is silent must
    still RECEIVE a spliced claim from a fluxed neighbour through the delivered directional
    message (pass-through propagation; the backbone gathers at the source)."""
    F, U = _U()

    class SpySolve(U.SilentSolve):
        def solve(self, own, forward, backward):
            self.seen = (own, forward, backward)
            return PsiMessage.silent()

    solver = SpySolve()
    _run(sweep_inputs, U.UnifiedPolicy(F.PassThroughPropagation(), solver))
    own, fwd, bwd = solver.seen
    own_p = np.asarray(own.spliced_rna_pos.precision)
    fwd_p = np.asarray(fwd.spliced_rna_pos.precision)
    bwd_p = np.asarray(bwd.spliced_rna_pos.precision)
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
        "no slot two hops from any flux received a spliced claim — the scans moved nothing"
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
    assert weak.precision < 10.0, "a hop with any disagreement must cost precision"


def test_the_premise_is_charged_on_every_hop():
    """An imputation must cost something every hop, even at a perfectly agreeing frame: with
    log_r = 0 the knob and disagreement vanish, and the premise alone must still lower the
    precision (the recorded free-hop failure is what this forbids)."""
    F, U, m = _hop_tables(log_r=0.0, premise=0.3)
    out = _one_hop(m, F, U, "unspliced_gdna", F.Claim(1.0, 10.0))
    assert np.isclose(out.abundance, 1.0)
    assert np.isclose(out.precision, 10.0 / (1.0 + 10.0 * 0.3))


def test_spliced_lanes_cross_unreframed():
    """Certified counts are capture-invariant (the anchor's measured property): a spliced claim
    crosses an enrichment frame with its VALUE untouched — no knob, no reframe."""
    F, U, m = _hop_tables(log_r=2.0, v_r=0.01)
    out = _one_hop(m, F, U, "spliced_rna_pos", F.Claim(3.0, 25.0))
    assert np.isclose(out.abundance, 3.0), "a measurement's value never reframes"


def test_frame_aware_under_a_silent_solve_changes_nothing(sweep_inputs):
    """The free invariant: propagation NEVER changes the answer on its own — a belief only forms
    at the solve, so FrameAwarePropagation under SilentSolve must equal SilentPolicy
    byte-for-byte through the real backbone."""
    F, U = _U()
    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, U.UnifiedPolicy(U.FrameAwarePropagation(), U.SilentSolve()))
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)
