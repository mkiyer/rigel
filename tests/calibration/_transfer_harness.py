"""The shared harness of the transfer policy's gate files (`test_transfer_*.py`): ONE real
`solve_chain` call captured from a calibrate run on the toy — the backbone-parity pattern, so every
gate re-runs the sweep with a different policy on byte-identical inputs — plus the policy, context
and pass builders they share. ⚠ Not a test module: `capture_sweep_inputs` is what each gate file's
``sweep_inputs`` fixture returns, and nothing here asserts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

import rigel.calibration.sweep as SW
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION


def capture_sweep_inputs(tmp_path_factory):
    """ONE real `solve_chain` call captured from a calibrate run on the toy — the backbone-parity
    pattern: every gate re-runs the sweep with a different policy on byte-identical inputs. Each gate
    file wraps this in its own module-scoped ``sweep_inputs`` fixture."""
    import dataclasses

    spec = importlib.util.spec_from_file_location(
        "tpo_for_transfer_policy", Path(__file__).parent / "test_prior_vs_oracle.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tpo_for_transfer_policy"] = m
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
            config=CalibrationConfig(message_propagation=True),
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
        payload=payload,
        calibrate_kw=dict(
            region_arrays=ra,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf,
            rna_fl_pmf=fl.rna_pmf,
            sj=build_sj_geometry_arrays(index),
            boundary_flags=build_boundary_flags_array(index),
        ),
    )


def _run(si, policy, capture=None):
    kw = dict(si["kw"])
    if capture is not None:
        kw["_capture"] = capture
    out = SW.solve_chain(*si["args"], **kw, policy=policy)
    return {
        f: np.asarray(getattr(out, f))
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg")
    }


def _expected_pairs(si):
    """The intron|exon pairs derived INDEPENDENTLY of the policy (its falsification power):
    a BOUNDARY whose one flank is an exon REGION and whose other flank is an intron REGION
    that admits RNA."""
    from rigel.calibration.signature import coarse_type_array

    chain, statics, _geometry, _belief, region_arrays = si["args"]
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_exon = is_reg & (rtype[np.clip(obj, 0, rtype.shape[0] - 1)] == 2)
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    is_intron = is_reg & ~is_exon & (fp | fn)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    pairs = []
    for i in np.flatnonzero(~is_reg):
        lo, hi = left[i], right[i]
        if lo < 0 or hi < 0:
            continue
        if is_exon[lo] and is_intron[hi]:
            pairs.append((int(i), int(hi)))
        elif is_exon[hi] and is_intron[lo]:
            pairs.append((int(i), int(lo)))
    return pairs, int(chain.n_slots)


def _live_provider(si, n_grid, window):
    """A synthetic factory-row provider: a distinct non-flat row at every intron REGION slot."""
    from rigel.calibration.simplex_logodds import _logodds_grid

    pairs, n_slots = _expected_pairs(si)
    lam, _ = _logodds_grid(n_grid, window)
    rows = np.zeros((n_slots, lam.shape[0]))
    for _b, j in pairs:
        rows[j] = -0.05 * (lam - (0.1 * (j % 7) - 0.3)) ** 2  # non-flat, slot-distinct
    return lambda g, w: rows if (int(g), float(w)) == (int(n_grid), float(window)) else None


def _ctx_of(si):
    """Rebuild the StepContext exactly as the backbone would, by running a silent sweep with a
    capture and reading nothing — instead we call solve_chain's own construction path via a spy
    on the policy prepare."""
    grabbed = []

    class _Spy:
        name = "ctx-spy"

        def prepare(self, ctx):
            grabbed.append(ctx)
            return SilentPolicy().prepare(ctx)

    _run(si, _Spy())
    assert grabbed, "the spy never fired"
    return grabbed[0]


def _nothing_held(n_slots):
    """The two held lists of a policy that sent nothing: SILENCE at every interior node."""
    from rigel.calibration.messages import SILENCE

    return [SILENCE] * int(n_slots)


def _strand_of(si):
    kw = si["kw"]
    return (
        float(kw["rna_sense_frac"]),
        float(kw.get("gdna_strand_overdispersion", 0.0)),
        float(kw.get("rna_strand_overdispersion", 0.0)),
    )


def _intron_mask(ctx):
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    return ~is_bnd & ~is_exon & (fp | fn)


def _dead_boundaries(ctx):
    """The context with every BOUNDARY's strand channel declared dead and every region's intact."""
    import dataclasses as _dc

    tau = np.asarray(ctx.own.tau_lam, np.float64).copy()
    tau[np.asarray(ctx.is_boundary, bool)] = 0.0
    return _dc.replace(ctx, own=_dc.replace(ctx.own, tau_lam=tau))


def _drive_the_backbone(prepared, ctx):
    """The backbone's own contract, reproduced (`sweep.solve_chain`'s two directional passes and its
    solve): each pass calls ``receive(source, destination)`` over the chain order — the forward pass
    reading each slot's LOW neighbour, the backward pass its HIGH one — and holds the result at the
    destination; ``solve`` receives the two held lists. Returns the delivered rows (zeros when the
    policy is silent)."""
    from rigel.calibration.messages import SILENCE

    order = list(range(ctx.n_slots))
    held = []
    for nbr, seq, backward in (
        (np.asarray(ctx.left, np.int64), order, False),
        (np.asarray(ctx.right, np.int64), order[::-1], True),
    ):
        receive = prepared.propagate(backward=backward)
        got = [None] * len(order)
        for i in seq:
            s = int(nbr[i])
            if s >= 0:
                got[i] = SILENCE if receive is None else receive(s, i)
        held.append(got)
    msg = prepared.solve(*held)
    if msg.lam_rows is None:
        return np.zeros((len(order), int(ctx.n_grid)))
    return np.asarray(msg.lam_rows)


def _rows_of(pol, ctx):
    return _drive_the_backbone(pol.prepare(ctx), ctx)


def _strand_row_of(ctx, strand, lam, x):
    """A slot's own strand profile, recomputed independently (the frozen-variance count form)."""
    kappa, od_g, od_r = strand
    fp = np.asarray(ctx.free_pos, bool)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    n = cnt[x].sum()
    ks = kappa if fp[x] else 1.0 - kappa
    f_ref = float(np.clip(belief[x], 1e-9, 1 - 1e-9))
    p = 0.5 * fg + ks * (1 - fg)
    p_ref = 0.5 * f_ref + ks * (1 - f_ref)
    var = max(
        n * p_ref * (1 - p_ref)
        + (n * f_ref) ** 2 * 0.25 * od_g
        + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
        1e-9,
    )
    row = -0.5 * (cnt[x, 0] - n * p) ** 2 / var
    return row - row.max()


def _full_policy(sweep_inputs):
    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    return TransferPolicy(provider, strand=_strand_of(sweep_inputs)), provider, n_grid, window


def _rna_lanes_of(sweep_inputs, ctx=None):
    pol = _full_policy(sweep_inputs)[0]
    ctx = _ctx_of(sweep_inputs) if ctx is None else ctx
    prepared = pol.prepare(ctx)
    assert set(prepared.lanes) == {"gdna", "pos", "neg"}
    return ctx, prepared


def _rna(prepared) -> dict:
    """The two RNA lanes by strand."""
    return {name: prepared.lanes[name] for name in ("pos", "neg")}


def _strand_intron(ctx, name):
    """The PER-STRAND intron test the lane uses: a region that admits ``s`` and carries no exon of ``s``."""
    is_bnd = np.asarray(ctx.is_boundary, bool)
    free = np.asarray(ctx.free_pos if name == "pos" else ctx.free_neg, bool)
    exon_s = np.asarray(ctx.exon_pos if name == "pos" else ctx.exon_neg, bool)
    return ~is_bnd & free & ~exon_s


def _held_rna(held, x, field):
    m = held[x]
    return None if (m is None or getattr(m, field) is None) else getattr(m, field).profile.copy()


def _with_populated_inside(ctx):
    """The toy's inside pieces are EMPTY (shorter than a fragment: no contained fragment, a contained
    opportunity below one base), so the level rule serves nothing there and a gate on the bare toy proves
    nothing. A context is data: populate the two inside slots with consistent counts (both strand
    channels live, so the step's spread is fitted from two pairs) and a plausible contained
    opportunity, for the policy and the independent recompute alike."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer_rows import outside_flank

    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    cnt = np.asarray(ctx.unspliced_count, np.float64).copy()
    a_g = np.asarray(ctx.eff_gdna, np.float64).copy()
    a_r = np.asarray(ctx.eff_rna, np.float64).copy()
    tau = np.asarray(ctx.own.tau_lam, np.float64).copy()
    fills = iter(
        [(20.0, 180.0), (80.0, 120.0), (25.0, 175.0), (70.0, 130.0)]
    )  # two pairs that DISAGREE
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        _o, i = outside_flank(flags[b], lo, hi)
        if i is None:
            continue
        cnt[i] = next(fills)
        a_g[i] = a_r[i] = 150.0
        tau[i] = 1.0
    return _dc.replace(
        ctx,
        unspliced_count=cnt,
        n_slot=cnt.sum(axis=1),
        eff_gdna=a_g,
        eff_rna=a_r,
        own=_dc.replace(ctx.own, tau_lam=tau),
    )


def _with_alt_splice_sites(ctx):
    """The toy carries no alternative splice site: turn its two exon|exon terminus boundaries into a
    DONOR (intron to the right) and an ACCEPTOR (intron to the left) with a route flux each, and
    populate the pieces beyond them — a context is data; the policy and the recompute read the same."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer_rows import TERMINUS
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_NEG, FLAG_DONOR_POS

    base = _with_populated_inside(ctx)
    is_bnd = np.asarray(base.is_boundary, bool)
    is_exon = np.asarray(base.is_exon_region, bool)
    left, right = np.asarray(base.left, np.int64), np.asarray(base.right, np.int64)
    flags = np.asarray(base.boundary_flags, np.uint16).copy()
    sjc = np.asarray(base.sj_count, np.float64).copy()
    kinds = iter([(FLAG_DONOR_POS, (40.0, 0.0)), (FLAG_ACCEPTOR_NEG, (0.0, 50.0))])
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]) or not (flags[b] & TERMINUS):
            continue
        kind, fl = next(kinds)
        flags[b] = kind
        sjc[b] = fl
    return _dc.replace(base, boundary_flags=flags, sj_count=sjc)
