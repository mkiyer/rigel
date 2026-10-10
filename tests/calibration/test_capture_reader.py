"""The capture reader's publication: the block solve reads every slot's posterior mode of gDNA density on
the LAST refit's sweep, and ``calibrate`` publishes the weights relative to the typical read object, 1 where nothing was
read. The reader itself is held to its specification in ``tests/native/test_honest_reader.py``."""

import importlib

import numpy as np
import pytest

import rigel.calibration.sweep as SWEEP
from rigel.calibration.region_chain import REGION
from rigel.config import CalibrationConfig

# the package exports the function `calibrate`, which shadows the module of the same name
CAL = importlib.import_module("rigel.calibration.calibrate")


def _run(sweep_inputs, iters):
    debug = {}
    res = CAL.calibrate(
        payload=sweep_inputs["payload"],
        config=CalibrationConfig(calib_refit_iters=iters),
        _debug=debug,
        **sweep_inputs["calibrate_kw"],
    )
    return res, debug


def test_no_refit_means_no_reader_and_every_weight_is_exactly_one(sweep_inputs):
    """With no refit there is no landscape to read under, so nothing is read and nothing contracts."""
    res, debug = _run(sweep_inputs, 0)
    assert debug["reader_mode"] is None
    np.testing.assert_array_equal(res.gdna_capture_efficiency_region, 1.0)
    np.testing.assert_array_equal(res.gdna_capture_efficiency_boundary, 1.0)


def test_the_reader_runs_once_on_the_last_refit_and_the_weights_are_its_modes_relative_to_the_median(
    sweep_inputs, monkeypatch
):
    """Three sweeps run for two refits; exactly the last one carries the reader. Its modes, published:
    ``exp(mode − median)`` on every read slot — the typical read object exactly 1 — and 1 on every unread
    slot, split onto the two object axes by the chain's own kind and index."""
    seen = []
    orig = SWEEP.solve_blocks

    def spy(*args, **kw):
        seen.append(kw.get("reader") is not None)
        return orig(*args, **kw)

    monkeypatch.setattr(SWEEP, "solve_blocks", spy)
    res, debug = _run(sweep_inputs, 2)
    assert seen == [False, False, True], seen
    mode = debug["reader_mode"]
    finite = np.isfinite(mode)
    assert finite.sum() > 0, "the reader read no slot"
    w = np.where(finite, np.exp(mode - np.median(mode[finite])), 1.0)
    assert np.median(w[finite]) == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.isfinite(w)) and np.all(w > 0.0)
    chain = debug["chain"]
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)
    np.testing.assert_array_equal(
        res.gdna_capture_efficiency_region[obj[kind == REGION]], w[kind == REGION]
    )
    np.testing.assert_array_equal(
        res.gdna_capture_efficiency_boundary[obj[kind != REGION]], w[kind != REGION]
    )


def test_a_forced_set_of_modes_is_published_exactly(sweep_inputs, monkeypatch):
    """PERTURBATION: with the solve's modes replaced by stated ones, the result carries exactly their
    relative weights — a reader whose answer did not reach the result would publish 1 everywhere."""
    orig = CAL._solve

    def forced(s, _debug):
        belief, hyperprior, mode = orig(s, _debug)
        n = int(s.chain.n_slots)
        mode = np.full(n, np.nan)
        mode[0], mode[1], mode[n - 1] = -3.0, -1.0, -2.0
        return belief, hyperprior, mode

    monkeypatch.setattr(CAL, "_solve", forced)
    res, debug = _run(sweep_inputs, 1)
    chain = debug["chain"]
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)
    both = np.concatenate(
        [res.gdna_capture_efficiency_region, res.gdna_capture_efficiency_boundary]
    )
    n_read = 3
    assert int((both != 1.0).sum()) == n_read - 1  # the median read slot is exactly 1
    for slot, want in ((0, np.exp(-1.0)), (1, np.exp(1.0)), (int(chain.n_slots) - 1, 1.0)):
        arr = (
            res.gdna_capture_efficiency_region
            if kind[slot] == REGION
            else res.gdna_capture_efficiency_boundary
        )
        assert arr[obj[slot]] == pytest.approx(want, rel=1e-15, abs=0.0)


def test_capture_weights_split_the_modes_onto_the_two_axes():
    """The publication function alone, on a hand-made chain: NaN is unread (1), the median finite mode is 1,
    regions and boundaries are addressed by the chain's kind and index."""
    from types import SimpleNamespace

    chain = SimpleNamespace(
        kind=np.array([REGION, REGION + 1, REGION, REGION + 1]), obj_idx=np.array([0, 0, 1, 1])
    )
    reg, bnd = CAL.capture_weights(
        np.array([np.log(2.0), np.nan, np.log(8.0), np.log(4.0)]), chain, 2, 2
    )
    np.testing.assert_allclose(reg, [0.5, 2.0], rtol=1e-15, atol=0.0)
    np.testing.assert_allclose(bnd, [1.0, 1.0], rtol=1e-15, atol=0.0)
    reg0, bnd0 = CAL.capture_weights(None, chain, 2, 2)
    np.testing.assert_array_equal(reg0, 1.0)
    np.testing.assert_array_equal(bnd0, 1.0)
