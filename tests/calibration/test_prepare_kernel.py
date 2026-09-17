"""``native.transfer_prepare`` — the composition transfer's builders for one block as ONE native call —
against the Python builders they replace (`transfer._prepare_reference`, the executable specification), on
the toy's captured sweep: every table of the two prepared objects compared — the claims, the face tables and
their row store, each lane's faces, own levels, junction flux levels and flux witnesses. The port is the same
arithmetic under the compiler's fused multiply-add, so every profile is held to 1e-9 (a max-normalised
log-row whose raw terms reach 1e4 on a deep node moves by a few ulps of those terms) and every boolean, index,
count, opportunity and witness to equality. Each case first proves its mechanism is live in the reference, so
a comparison of two silent tables cannot pass it. The wiring gate proves `TransferPolicy.prepare` reaches the
native builders.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from _transfer_harness import (
    _ctx_of,
    _full_policy,
    _with_alt_splice_sites,
    _with_populated_inside,
    capture_sweep_inputs,
)

import rigel.calibration.messages.transfer as TP
from rigel.calibration.messages.faces import LEVEL, SPLICE_OUT


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


#: the tables that are copies or exact sums of the context's numbers — held to equality
_EXACT = ("count", "a", "total", "other", "flux_witness.rows", "n_u", "n_s", "a_b", "a_x")


def _leaves(p):
    """Every table of a prepared object, by name."""
    yield "own.mask", p.own.mask
    yield "own.rows", p.own.rows
    f = p.faces
    for name in ("kind", "row", "row2", "n_u", "n_s", "a_b", "a_x", "width", "var"):
        yield f"faces.{name}", getattr(f, name)
    yield "faces.rows", f.rows[: f.n_rows]
    for k in sorted(p.lanes):
        ln = p.lanes[k]
        for name in ("face", "two_sided", "empty", "count", "a", "total"):
            yield f"{k}.{name}", getattr(ln, name)
        yield f"{k}.other", np.zeros(0) if ln.other is None else ln.other
        yield f"{k}.rho_ref", np.array(ln.rho_ref)
        for t in ("own_level", "flux", "flux_witness"):
            yield f"{k}.{t}.mask", getattr(ln, t).mask
            yield f"{k}.{t}.rows", getattr(ln, t).rows
    yield "site.ambig", p.site.ambig


def _both(pol, ctx, library=None):
    src = ctx.factory_rows
    if src is None:
        src = np.zeros((int(ctx.n_slots), int(ctx.n_grid)))
    chain = TP._Chain(ctx, np.asarray(src, np.float64), pol._strand)
    library = pol.library(ctx) if library is None else library
    return TP._prepare_reference(chain, library), TP._prepare_native(chain, library)


def _assert_same_tables(ref, nat):
    assert set(ref.lanes) == set(nat.lanes)
    assert ref.faces.n_rows == nat.faces.n_rows > 0
    for (name, x), (_n, y) in zip(_leaves(ref), _leaves(nat)):
        x, y = np.asarray(x), np.asarray(y)
        assert x.shape == y.shape and x.dtype == y.dtype, name
        if x.dtype == bool or x.dtype.kind in "iu" or name.split(".", 1)[-1] in _EXACT:
            np.testing.assert_array_equal(x, y, err_msg=name)
        else:
            np.testing.assert_allclose(x, y, rtol=0.0, atol=1e-9, err_msg=name)


def test_the_native_builders_write_the_python_builders_tables(sweep_inputs):
    """The toy's captured sweep, the shipped policy: claims at introns, exons and boundaries; forward,
    transport and splice-out rules; all three lanes with own levels and flux levels — every table the same."""
    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ref, nat = _both(pol, _ctx_of(sweep_inputs))
    assert ref.own.mask.any() and ref.faces.any() and set(ref.lanes) == {"gdna", "pos", "neg"}
    assert all(ref.lanes[k].own_level.mask.any() for k in ref.lanes)
    assert any(ref.lanes[k].flux.mask.any() for k in ("pos", "neg"))
    _assert_same_tables(ref, nat)


def test_the_native_terminus_level_rule_is_the_pythons(sweep_inputs):
    """The inside pieces populated (`_with_populated_inside`): the level rule fires, with the pair's own
    discrepancies — its map, its bound and its width the same."""
    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ref, nat = _both(pol, _with_populated_inside(_ctx_of(sweep_inputs)))
    assert (ref.faces.kind == LEVEL).any(), "the level rule must fire or this gate proves nothing"
    assert (ref.faces.var[ref.faces.kind == LEVEL] > 0.0).all()
    _assert_same_tables(ref, nat)


def test_the_native_alternative_splice_site_is_the_pythons(sweep_inputs):
    """Two alternative splice sites with route flux (`_with_alt_splice_sites`): both flanks served, the
    pair's width beyond counting the same."""
    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ref, nat = _both(pol, _with_alt_splice_sites(_ctx_of(sweep_inputs)))
    served = (ref.faces.kind == SPLICE_OUT) & (ref.faces.width > 0.0)
    assert served.any(), (
        "a splice-out rule with a pair width must exist or this gate proves nothing"
    )
    _assert_same_tables(ref, nat)


def test_the_native_builders_without_a_strand_model_and_without_a_gdna_coordinate(sweep_inputs):
    """No strand model: no strand claim, no alternative splice site, the lanes still built. No gDNA
    coordinate (a gDNA-free library): no gDNA lane, the RNA lanes still built — the same in both."""
    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    ref, nat = _both(TP.TransferPolicy(), ctx)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    assert not ref.own.mask[is_exon].any() and ref.faces.any()
    _assert_same_tables(ref, nat)
    library = dataclasses.replace(pol.library(ctx), rho_gdna=0.0)
    ref, nat = _both(pol, ctx, library)
    assert set(ref.lanes) == {"pos", "neg"} and ref.lanes["pos"].own_level.mask.any()
    _assert_same_tables(ref, nat)


def test_the_policy_prepares_natively(sweep_inputs, monkeypatch):
    """The wiring: `TransferPolicy.prepare` reaches `transfer_prepare` once per block — a spy counts the
    calls — and the object it returns carries claims, rules and lanes."""
    calls = []
    real = TP.transfer_prepare

    def spy(**kw):
        calls.append(1)
        return real(**kw)

    monkeypatch.setattr(TP, "transfer_prepare", spy)
    pol, _rows, _g, _w = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = pol.prepare(ctx, pol.library(ctx))
    assert calls == [1]
    assert prepared.own.mask.any() and prepared.faces.any() and prepared.faces.n_rows > 0
    assert set(prepared.lanes) == {"gdna", "pos", "neg"}
