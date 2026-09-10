"""THE ALTERNATIVE SPLICE SITE (`transfer._alternative_splice_site`): the junction flanks read
from the flag kind alone, and on a toy given alternative sites both flanks' rules — the intron-side
flank on the full crossing, the exon-of-both flank on the crossing plus the leaving isoform's flux —
each pair at its own width beyond counting, recomputed independently."""

from __future__ import annotations

import numpy as np
import pytest

from _transfer_harness import (
    _ctx_of,
    _full_policy,
    _strand_of,
    _with_alt_splice_sites,
    capture_sweep_inputs,
)


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


def test_the_junction_flanks_read_the_flag_kind_alone():
    """The flank predicate, gated DIRECTLY: a DONOR bit puts the intron to the right on either strand
    (the flags are genomic-order), an ACCEPTOR bit to the left; a terminus on the boundary, no junction,
    or junctions both ways give no flanks."""
    from rigel.calibration.messages.transfer_rows import junction_flanks
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
        FLAG_TSS_POS,
    )

    assert junction_flanks(FLAG_DONOR_POS, 7, 9) == (9, 7)
    assert junction_flanks(FLAG_DONOR_NEG, 7, 9) == (9, 7)
    assert junction_flanks(FLAG_ACCEPTOR_POS, 7, 9) == (7, 9)
    assert junction_flanks(FLAG_ACCEPTOR_NEG, 7, 9) == (7, 9)
    assert junction_flanks(FLAG_DONOR_POS | FLAG_TSS_POS, 7, 9) == (None, None), (
        "sj+terminus: its own item"
    )
    assert junction_flanks(0, 7, 9) == (None, None)
    assert junction_flanks(FLAG_DONOR_POS | FLAG_ACCEPTOR_NEG, 7, 9) == (None, None), (
        "junctions both ways"
    )


def _expected_alt_splice_rows(si, ctx, strand, lam):
    """Item 7 recomputed INDEPENDENTLY of the policy on the patched context: at each junction boundary
    the two flanks' own strand rows through the splice-out map (E with S_b + F, C with S_b) into the
    boundary and the boundary's own row through the face map into each flank, each blurred by the
    pair's OWN disagreement beyond counting (the owner's discrepancy rule, nothing pooled). Keyed by
    slot; the pair widths returned beside."""
    from rigel.calibration.messages.transfer_rows import (
        blur_row,
        boundary_shares_strand,
        face_map_lambda,
        junction_flanks,
        splice_out_row,
        transport_row,
    )

    kappa, od_g, od_r = strand
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    flux = np.asarray(ctx.sj_count, np.float64).sum(axis=1)
    A_g = np.asarray(ctx.eff_gdna, np.float64)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))

    def strand_row(x):
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

    served = []
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        c_side, e_side = junction_flanks(flags[b], lo, hi)
        if c_side is None or not (n_u[b] > 0 and A_g[b] > 0):
            continue
        for x, s_out, kind in ((e_side, n_s[b] + flux[b], "E"), (c_side, n_s[b], "C")):
            if boundary_shares_strand(fp[b], fn[b], fp[x], fn[x]) and A_g[x] > 0:
                served.append((int(b), int(x), float(s_out), kind))
    width = {}
    for b, x, s_out, _kind in served:
        if not (tau[b] > 0.0 and tau[x] > 0.0):
            continue
        lo_v = []
        for y in (b, x):
            n = cnt[y].sum()
            ks = kappa if fp[y] else 1.0 - kappa
            p = cnt[y, 0] / n
            f = (p - ks) / (0.5 - ks)
            if not 0.0 < f < 1.0:
                break
            lo_v.append((np.log(f / (1 - f)), p * (1 - p) / n / (p - ks) ** 2 / (1 - f) ** 2))
        if len(lo_v) < 2:
            continue
        v_ratio = s_out / (n_u[b] * (n_u[b] + s_out)) if s_out > 0 else 0.0
        d = lo_v[0][0] - lo_v[1][0] - np.log((n_u[b] + s_out) / n_u[b])
        width[(b, x)] = max(0.0, d * d - (lo_v[0][1] + lo_v[1][1] + v_ratio))
    out = {}
    for b, x, s_out, _kind in served:
        w = width.get((b, x), 0.0)
        if tau[x] > 0.0 and fp[x] != fn[x]:
            row = splice_out_row(strand_row(x), lam, n_u[b], s_out, A_g[b], A_g[x])
            if np.ptp(row) > 1e-9:
                out.setdefault(b, np.zeros(lam.shape[0]))
                out[b] += blur_row(row, lam, w)
        if tau[b] > 0.0:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[x], A_g[x], s_out / A_g[b])
            row = transport_row(strand_row(b), lam, le, n_u[b], s_out)
            if np.ptp(row) > 1e-9:
                out.setdefault(x, np.zeros(lam.shape[0]))
                out[x] += blur_row(row, lam, w)
    return out, width


def test_the_alt_splice_rules_carry_both_flanks_with_the_pair_width(sweep_inputs):
    """ITEM 7 as claims and rules on the patched toy: the rule flank → junction boundary applied to the
    flank's claim, summed over the two flanks, equals the independently recomputed rows at the
    boundary; the rule boundary → flank applied to the boundary's claim equals the recomputed row at
    each flank — every message blurred by its own pair's disagreement beyond counting."""
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _with_alt_splice_sites(_ctx_of(sweep_inputs))
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    expected, widths = _expected_alt_splice_rows(sweep_inputs, ctx, strand, lam)
    assert expected and any(w > 0.0 for w in widths.values()), (
        "the patched toy must carry served junctions with a live pair width or this gate proves nothing"
    )
    is_bnd = np.asarray(ctx.is_boundary, bool)
    got = {}
    for (s, d), fn in prepared.rule.items():
        if not (
            (is_bnd[s] and s in widths_keys(widths)) or (is_bnd[d] and d in widths_keys(widths))
        ):
            continue
        if prepared.own[s] is None:
            continue
        r = fn(prepared.own[s], None)
        if r is not None and np.ptp(r) > 1e-9:
            got[int(d)] = got.get(int(d), 0.0) + (r - r.max())
    for slot, want in expected.items():
        assert slot in got, f"slot {slot} has no live item-7 rule"
        np.testing.assert_allclose(
            got[slot] - got[slot].max(), want - want.max(), rtol=0, atol=1e-10
        )


def widths_keys(widths):
    """The junction boundaries the recompute served (its width keys are ``(boundary, flank)`` pairs)."""
    return {int(k[0]) if isinstance(k, tuple) else int(k) for k in widths}
