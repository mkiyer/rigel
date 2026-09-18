"""``native.transfer_pass`` — the composition transfer's directional pass as ONE native call — against
the Python kernel it replaces (`sweep._pass` running ``propagate``'s ``receive(source, destination)``
hop by hop), on the toy's captured sweep: both passes of the shipped policy, every table compared field
by field — a row matrix under its presence bits, since the matrices are allocated unfilled and an absent
row is no row. The port is the same arithmetic in another summation order, so profiles are held to 1e-12 and
every boolean, count and witness to equality; ``trigamma`` is held to SciPy's Hurwitz zeta. The wiring
gate proves the backbone reaches the native kernel for the shipped policy and still runs the per-hop
path for a policy that offers no whole-pass kernel.
"""

from __future__ import annotations

import numpy as np
import pytest
from _transfer_harness import (
    _ctx_of,
    _full_policy,
    _leaves,
    _native_passes,
    _passes,
    _prepared,
    capture_sweep_inputs,
)

import rigel.calibration.sweep as SW
from rigel.native import trigamma


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


def test_trigamma_is_scipys_hurwitz_zeta():
    """The counting variance's one home in C++ agrees with ``zeta(2, x)`` from a half to ten million."""
    from scipy.special import zeta

    x = np.concatenate([np.linspace(0.5, 12.0, 400), np.logspace(1.0, 7.0, 200)])
    np.testing.assert_allclose(np.array([trigamma(float(v)) for v in x]), zeta(2, x), rtol=1e-13)


def test_the_native_pass_is_the_python_kernels_pass_on_every_table(sweep_inputs):
    """Both passes of the shipped policy on the toy's captured sweep, the per-hop kernel against the
    native call: every presence bit, count, opportunity and witness equal, every present row within
    1e-12, and the tables not trivially silent."""
    pol, _rows, _n_grid, _window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(pol, ctx)
    python = _passes(prepared, ctx)
    native = _native_passes(prepared, ctx)
    assert python[0].heard.any() and python[1].heard.any(), "the toy's sweep must carry messages"
    for py, nat in zip(python, native):
        for name, x, y in _leaves(py, nat):
            if x.dtype == bool:
                np.testing.assert_array_equal(x, y, err_msg=name)
            else:
                both_nan = np.isnan(x) & np.isnan(y)
                np.testing.assert_array_equal(np.isnan(x), np.isnan(y), err_msg=name)
                np.testing.assert_allclose(
                    np.where(both_nan, 0.0, x),
                    np.where(both_nan, 0.0, y),
                    rtol=0.0,
                    atol=1e-12,
                    err_msg=name,
                )


def test_the_backbone_runs_the_shipped_policys_pass_natively(sweep_inputs, monkeypatch):
    """The wiring: `sweep._pass` reaches `transfer_pass` for the transfer policy — a spy counts the
    calls and the table it fills carries messages — while the silent policy, which offers no whole-pass
    kernel, runs the per-hop path and yields an all-silence table with its neighbours marked."""
    import rigel.calibration.messages.transfer as TP
    from rigel.calibration.messages.silent import SilentPolicy

    calls = []
    real = TP.transfer_pass

    def spy(**kw):
        calls.append(1)
        return real(**kw)

    monkeypatch.setattr(TP, "transfer_pass", spy)
    pol, _rows, n_grid, _window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    order = list(range(int(ctx.n_slots)))
    from_left = SW._pass(order, list(ctx.left), _prepared(pol, ctx), n_grid, backward=False)
    assert calls == [1] and from_left.heard.any()
    silent = SilentPolicy().prepare(ctx, SilentPolicy().library(ctx))
    table = SW._pass(order, list(ctx.left), silent, n_grid, backward=False)
    assert calls == [1] and table.has_neighbour.any() and not table.heard.any()


def test_the_pass_reads_the_builders_own_tables_without_a_copy(sweep_inputs):
    """The tables `run_pass` hands the native kernel are the builders' own arrays — the claims' matrix
    and mask, the faces' arrays and the written prefix of their row store, each lane's arrays — every
    one C-contiguous (the binding would copy a strided one silently), none built at the call.
    PERTURBATION: a builder that keeps a strided column, or a table packed at the call, fails here."""
    pol, _rows, _n_grid, _window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = _prepared(pol, ctx)
    t = prepared.tables()
    assert t["own"] is prepared.own.rows and t["own_mask"] is prepared.own.mask
    assert t["f_rows"].shape[0] == prepared.faces.n_rows > 0
    assert np.shares_memory(t["f_rows"], prepared.faces.rows)
    assert {lane[0] for lane in t["lanes"]} == {0, 1, 2}
    for lane, ln in zip(t["lanes"], prepared.lanes.values()):
        assert lane[4] is ln.own_level.rows and lane[9] is ln.flux_witness.rows
        assert lane[6] is ln.count and lane[7] is ln.a
    arrays = [(k, v) for k, v in t.items() if k != "lanes"] + [
        (f"lane[{j}][{e}]", x)
        for j, lane in enumerate(t["lanes"])
        for e, x in enumerate(lane)
        if isinstance(x, np.ndarray)
    ]
    for name, a in arrays:
        assert a.flags.c_contiguous, f"{name} is not C-contiguous"
