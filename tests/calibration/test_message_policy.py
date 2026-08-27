"""Gates for THE MESSAGE POLICY (`calibration.messages.policy`) — the skeleton of the
ground-up rebuild (owner directive, 2026-08-27).

The skeleton is anchored before it claims anything: MessagePolicy with a pass-through
propagation and a silent solve reproduces SilentPolicy byte-for-byte through the REAL backbone
(rung 0 IS the silent floor), and the spliced-lane seam reproduces the relay running only its
certified-flux stream. Every rung built on top must arrive with its own fail-first gate here.
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
    from rigel.calibration.messages import foundation, policy as unified

    return foundation, unified


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """ONE real `solve_chain` call captured from a calibrate run on the toy — the
    backbone-parity pattern: every gate below re-runs the sweep with a different policy on
    byte-identical inputs."""
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "tpo_for_message_policy", Path(__file__).parent / "test_prior_vs_oracle.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tpo_for_message_policy"] = m
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
    pol = U.MessagePolicy(F.PassThroughPropagation(), U.SilentSolve())
    assert isinstance(pol, Policy)


def test_silent_solve_is_byte_identical_to_the_silent_policy(sweep_inputs):
    """IDENTITY ANCHOR ⓐ (the owner's rule made executable): silence is the solve ignoring the
    messages — nothing else changes, because a node's belief only updates at the solve. The
    runner with a pass-through propagation and a silent solve must reproduce SilentPolicy
    byte-for-byte through the REAL backbone."""
    F, U = _U()
    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, U.MessagePolicy(F.PassThroughPropagation(), U.SilentSolve()))
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
        U.MessagePolicy(F.PassThroughPropagation(), U.RowsSolve(rows)),
    )
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


def test_own_messages_carry_the_node_banks(sweep_inputs):
    """`prepare` must build each node's OWN message from the banks the context already carries:
    the unspliced lanes mirror the self-solve densities/precisions, and the spliced lanes are
    live exactly where the slot carries certified sj flux."""
    F, U = _U()

    class Recorder(U.MessagePolicy):
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
    lane (spliced lanes are one-hop by the SKELETON's law; this gate runs pass-through, so any
    unspliced lane works)."""
    F, U = _U()

    class SpySolve(U.SilentSolve):
        def solve(self, own, forward, backward):
            self.seen = (own, forward, backward)
            return PsiMessage.silent()

    solver = SpySolve()
    _run(sweep_inputs, U.MessagePolicy(F.PassThroughPropagation(), solver))
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
