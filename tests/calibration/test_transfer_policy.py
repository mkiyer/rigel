"""Gates for THE TRANSFER POLICY's shape (`calibration.messages.transfer`): the backbone protocol,
the silence identity (an evidence-free transfer is byte-identical to `SilentPolicy` through the real
backbone), the installer, the two passes against an independent recursive reference and the no-echo
law, the own claims' liveness under the strand deadband, and the completion contract — every
directed face carries a rule or a lane face, or leads into structural pure gDNA.

The messages themselves are gated by subject in the sibling files, one per builder:
`test_transfer_splice_faces.py` (the intron|exon face), `test_transfer_edge_and_terminus.py` (the
edge's level and the terminus), `test_transfer_alt_splice.py`, `test_transfer_gdna_lane.py`,
`test_transfer_rna_lanes.py` (the lanes and the cube) and `test_transfer_ceilings.py`; the shared
fixture and builders are `_transfer_harness.py`.

⚠ Falsification note (2026-09-01, watched): dropping the exon-flank requirement from the pair
predicate is NOT catchable on this toy — every intron-flanking boundary here also has an exon
flank — so the perturbation that must fire instead is the SOURCE-SIDE flip (exon as source),
which the pair gate catches.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy
from rigel.calibration.messages.silent import SilentPolicy
from _transfer_harness import (
    _ctx_of,
    _dead_boundaries,
    _drive_the_backbone,
    _full_policy,
    _run,
    _with_alt_splice_sites,
    capture_sweep_inputs,
)


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    return capture_sweep_inputs(tmp_path_factory)


def test_the_transfer_policy_satisfies_the_backbone_protocol():
    from rigel.calibration.messages.transfer import TransferPolicy

    assert isinstance(TransferPolicy(lambda g, w: None), Policy)


def test_an_evidence_free_transfer_is_byte_identical_to_silence(sweep_inputs):
    """THE RUNG-0 IDENTITY: with no factory rows to transfer, the policy must reproduce
    `SilentPolicy` byte-for-byte through the real backbone — and it must do so by delivering a
    SILENT message, never zero-filled channel arrays (the measured 1-ULP non-identity of a
    present-but-zero channel)."""
    from rigel.calibration.messages.transfer import TransferPolicy

    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, TransferPolicy(lambda g, w: None))
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


def test_the_policy_name_installs_the_transfer_policy(sweep_inputs):
    """`message_policy = "transfer"` must install `TransferPolicy` (an unreadable knob is worse
    than no knob), and the unknown-name refusal must survive the new branch."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.config import CalibrationConfig

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    seen: list = []
    orig = SW.solve_chain

    def spy(*a, **kw):
        seen.append(kw.get("policy"))
        return orig(*a, **kw)

    calibrate_mod.solve_chain = spy
    try:
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(
                CalibrationConfig(), message_propagation=True, message_policy="transfer"
            ),
            **sweep_inputs["calibrate_kw"],
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert seen and all(isinstance(p, TransferPolicy) for p in seen)
    assert all(p._strand is not None for p in seen), (
        "the exon -> boundary message must be ON in production"
    )
    with pytest.raises(ValueError, match="unknown message_policy"):
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(
                CalibrationConfig(), message_propagation=True, message_policy="no-such-policy"
            ),
            **sweep_inputs["calibrate_kw"],
        )


def _reference_rows(prepared, ctx):
    """AN INDEPENDENT IMPLEMENTATION OF THE TWO PASSES — recursive rather than sequential: the message
    into ``i`` from its neighbour ``s`` is the face's rule applied to ``s``'s own claim composed with
    the message into ``s`` from ``s``'s OTHER neighbour, and the rows are the two messages composed.
    Nothing here reads the policy's pass state; only its claims and rules."""
    from rigel.calibration.messages.transfer_rows import EPS

    left, right = list(ctx.left), list(ctx.right)
    own, rule = prepared.own, prepared.rule
    memo = {}

    def norm(r):
        return r - r.max()

    def into(s, i):
        key = (s, i)
        if key in memo:
            return memo[key]
        fn = rule.get((s, i))
        out = None
        if fn is not None:
            far = left[s] if right[s] == i else right[s]
            m = into(far, s) if far >= 0 else None
            r = fn(own[s], m)
            if r is not None and np.ptp(r) > EPS:
                out = norm(r)
        memo[key] = out
        return out

    lane = prepared.lanes.get("gdna")
    lmemo = {}

    def level_into(s, i):
        """THE LANE, recursively: nothing unless the face is a lane face; an EMPTY sender forwards the
        level that reaches it from its far side unchanged; a full sender sends the product of its own
        level and the level it holds; a full recipient prices the hop (both totals' counting plus the
        discrepancy beyond it) and takes the level as a lower bound."""
        from rigel.calibration.messages.transfer_rows import blur_row, hop_price, lower_side

        key = (s, i)
        if key in lmemo:
            return lmemo[key]
        out = None
        if lane is not None and (s, i) in lane.faces:
            far = left[s] if right[s] == i else right[s]
            held = level_into(far, s) if far >= 0 else None
            if lane.empty[s]:
                out = held
            else:
                own_s = lane.own_level[s]
                if own_s is not None:
                    own_s = np.maximum.accumulate(own_s)  # its lower side
                parts = [q for q in (own_s, None if held is None else held[0]) if q is not None]
                if parts:
                    tight = parts[0] if len(parts) == 1 else np.minimum(parts[0], parts[1])
                    out = (norm(tight), float(lane.count[s]), float(lane.a[s]))
            if out is not None and not lane.empty[i]:
                v = hop_price(out[1], out[2], lane.count[i], lane.a[i])
                priced = lower_side(out[0])
                priced = blur_row(priced, lane.u, v) if v > 0.0 else priced
                out = (priced, float(lane.count[i]), float(lane.a[i]))
        lmemo[key] = out
        return out

    n = len(own)
    rows = np.zeros((n, int(ctx.n_grid)))
    for i in range(n):
        parts = [into(s, i) for s in (left[i], right[i]) if s >= 0]
        if lane is not None and not lane.empty[i]:
            from rigel.calibration.messages.transfer_rows import profile_of_level

            bounds = []
            for s in (left[i], right[i]):
                lv = level_into(s, i) if s >= 0 else None
                if lv is not None:
                    bounds.append(
                        profile_of_level(
                            lv[0], lane.u, lane.lam, lane.total[i], lane.a[i], lane.rho_ref
                        )
                    )
            if bounds:  # two bounds on one density intersect: the pointwise tighter
                parts.append(
                    norm(bounds[0] if len(bounds) == 1 else np.minimum(bounds[0], bounds[1]))
                )
        parts = [p for p in parts if p is not None]
        if parts:
            rows[i] = norm(sum(parts))
    return rows


def test_the_backbone_rows_equal_an_independent_recursive_reference_of_the_passes(sweep_inputs):
    """THE PASSES, gated against a second implementation: what the backbone's two sequential passes
    deliver equals the recursive definition — the message into a node is the face's rule applied to
    the sender's own claim composed with what reached the sender from ITS far side — on the live toy
    and on the toy with populated inside pieces and alternative splice sites (every rule family live)."""
    pol, _p, _g, _w = _full_policy(sweep_inputs)
    for ctx in (_ctx_of(sweep_inputs), _with_alt_splice_sites(_ctx_of(sweep_inputs))):
        prepared = pol.prepare(ctx)
        rows = _drive_the_backbone(prepared, ctx)
        assert rows.any(), "the toy delivered nothing — this gate would prove nothing"
        np.testing.assert_allclose(rows, _reference_rows(prepared, ctx), rtol=0, atol=1e-10)


def test_PERTURBATION_no_node_ever_hears_its_own_claim_back(sweep_inputs):
    """THE NO-ECHO LAW, watched: replace ONE node's own claim by a distinctive profile and re-run the
    passes — what that node HOLDS from either side must not move (its claim never returns to it),
    while some other node's rows must (the claim did travel). Checked at an intron with two live faces
    and at a strand-live exon."""
    from rigel.calibration.messages.transfer_rows import EPS

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    prepared = pol.prepare(ctx)
    lam = np.linspace(-window, window, n_grid)
    probes = [
        i
        for i in range(ctx.n_slots)
        if prepared.own[i] is not None and np.ptp(prepared.own[i]) > EPS
    ]
    assert probes, "no node carries a claim — this gate would prove nothing"
    checked = 0
    for i in probes[:12]:
        base_rows = _drive_the_backbone(prepared, ctx)
        held_before = [
            None if h is None else h.composition
            for h in (prepared.held[False][i], prepared.held[True][i])
        ]
        saved = prepared.own[i]
        prepared.own[i] = -0.5 * ((lam - 3.3) / 0.2) ** 2  # a spike nowhere near any real claim
        rows = _drive_the_backbone(prepared, ctx)
        held_after = [
            None if h is None else h.composition
            for h in (prepared.held[False][i], prepared.held[True][i])
        ]
        prepared.own[i] = saved
        for a, b in zip(held_before, held_after):
            if a is None or b is None:
                assert a is None and b is None
            else:
                np.testing.assert_array_equal(a, b, err_msg=f"slot {i} heard its own claim back")
        moved = [
            j for j in range(ctx.n_slots) if j != i and not np.array_equal(rows[j], base_rows[j])
        ]
        if moved:
            checked += 1
    assert checked >= 1, "no perturbed claim travelled anywhere — the gate could not have fired"


def test_a_dead_strand_channel_carries_no_own_claim(sweep_inputs):
    """THE VACUITY LAW for every strand-borne claim (items 1, 2, 5, 6, 7): where the solver's derived
    deadband declares a node's strand channel dead (``own.tau_lam == 0``), that node's OWN CLAIM is
    absent — an exon's and a boundary's alike — so nothing of its own can travel; and a policy built
    without strand parameters carries no strand claim anywhere."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy

    pol, provider, _g, _w = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    live = pol.prepare(ctx)
    assert any(live.own[e] is not None for e in np.flatnonzero(is_exon)), "no live exon claim"
    dead = pol.prepare(_dc.replace(ctx, own=_dc.replace(ctx.own, tau_lam=np.zeros(n))))
    for x in range(n):
        if is_exon[x]:
            assert dead.own[x] is None, f"a dead exon {x} carries a claim"
        elif is_bnd[x] and dead.own[x] is not None:
            assert not dead.own[x].any(), f"a dead boundary {x} carries a strand claim"
    no_strand = TransferPolicy(provider).prepare(ctx)
    for x in range(n):
        if is_exon[x]:
            assert no_strand.own[x] is None
    # boundaries dead, exons live: the boundary's own claim is absent at every boundary
    dead_b = pol.prepare(_dead_boundaries(ctx))
    for b in np.flatnonzero(is_bnd):
        assert dead_b.own[b] is None or not dead_b.own[b].any()


def test_every_directed_face_is_served_or_faces_structural_pure_gdna(sweep_inputs):
    """THE COMPLETION CONTRACT, structurally: after `prepare`, every directed face (x → y) of the chain
    carries a composition rule or is a lane face, unless y is an intergenic region (structurally pure
    gDNA: nothing to impute there). STOP by omission is impossible by construction."""
    ctx = _ctx_of(sweep_inputs)
    prepared = _full_policy(sweep_inputs)[0].prepare(ctx)
    lane = prepared.lanes.get("gdna")
    assert lane is not None and lane.faces
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    unserved = [
        (int(x), int(y))
        for x in range(int(ctx.n_slots))
        for y in (left[x], right[x])
        if y >= 0
        and not intergenic[y]
        and (int(x), int(y)) not in prepared.rule
        and (int(x), int(y)) not in lane.faces
    ]
    assert not unserved, f"faces with no rule and no lane: {unserved[:8]}"
    # and a lane face never doubles a composition rule (no witness counted twice)
    assert not (set(prepared.rule) & lane.faces)
