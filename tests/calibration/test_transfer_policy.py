"""Gates for THE TRANSFER POLICY (`calibration.messages.transfer`) — rungs 1 and 2 of the
ground-up message rebuild (owner rulings, 2026-09-01): COMPOSITION TRANSFER, strictly one hop —
the intron's factory row verbatim at intron|exon boundaries (rung 1) and transported through
the face map into exons (rung 2).

The rung's anchors, in the fail-first order they were built:

* an evidence-free transfer is byte-identical to `SilentPolicy` through the real backbone —
  the rung-0 identity every new rung must keep;
* a live transfer delivers `lam_rows` at exactly the intron|exon PAIR boundaries (rung 1:
  each row its intron source's row VERBATIM — the hop cost was measured at zero beyond the
  row's own ``alpha_eff`` width) and at the LICENSED exon faces (rung 2: the intron row
  transported through the face map and widened by the face's own counting variance), zero
  everywhere else, and the policy relays NOTHING (one hop is structural);
* the face map is MONOTONE with the certified-flux CEILING (a measured spliced density caps the
  exon's claimable f_g; s = 0 degenerates to a pure opportunity shift);
* the ingredient width adds EXACTLY the face's counting variance (trigamma(n_u+1/2) +
  trigamma(n_s+1/2)) to a transported Gaussian factor — an inverted or mis-scaled kernel fails;
* `message_policy = "transfer"` actually installs the policy, and the unknown-name refusal
  stays intact.

⚠ Falsification note (2026-09-01, watched): dropping the exon-flank requirement from the pair
predicate is NOT catchable on this toy — every intron-flanking boundary here also has an exon
flank — so the perturbation that must fire instead is the SOURCE-SIDE flip (exon as source),
which the pair gate catches.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION


def _mp():
    spec = importlib.util.spec_from_file_location(
        "tmp_for_transfer_policy", Path(__file__).parent / "test_message_policy.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tmp_for_transfer_policy"] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """The message-policy gate file's captured `solve_chain` inputs, reused verbatim so every
    policy gate in this package runs on byte-identical inputs."""
    return _mp().sweep_inputs.__wrapped__(tmp_path_factory)


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


def _expected_exon_rows(si, ctx, src, lam):
    """The rung-2 exon rows recomputed INDEPENDENTLY of the policy: for every licensed
    intron|exon face of every exon, the intron row pushed through the face map and widened by
    the face's counting variance — a second implementation, so a policy bug cannot hide."""
    from scipy.special import polygamma

    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    term = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left = np.asarray(ctx.left, np.int64)
    right = np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
    A_r = np.asarray(ctx.eff_rna, np.float64)
    rr_lo = np.asarray(ctx.route_rate_lo, np.float64).sum(axis=1)
    rr_hi = np.asarray(ctx.route_rate_hi, np.float64).sum(axis=1)
    sc_lo = np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1)
    sc_hi = np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1)
    sig = 1.0 / (1.0 + np.exp(-lam))
    dlam = float(lam[1] - lam[0])
    out = {}
    for e in np.flatnonzero(is_exon):
        add = np.zeros(lam.shape[0])
        live = False
        for b, hi in ((left[e], True), (right[e], False)):
            if b < 0 or not is_bnd[b]:
                continue
            i = left[b] if right[b] == e else right[b]
            if i < 0 or not is_intron[i] or (flags[b] & term):
                continue
            if fp[e] != fp[i] or fn[e] != fn[i]:
                continue
            row_i = src[i]
            if np.ptp(row_i) <= 1e-9 or not (
                n_u[b] > 0 and A_g[b] > 0 and A_r[b] > 0 and A_g[e] > 0 and A_r[e] > 0
            ):
                continue
            s = float((rr_hi if hi else rr_lo)[b])
            n_s = float((sc_hi if hi else sc_lo)[b])
            le = np.log(n_u[b] * sig / A_g[b] * A_g[e]) - np.log(
                (n_u[b] * (1.0 - sig) / A_r[b] + s) * A_r[e]
            )
            r = np.interp(np.interp(lam, le, lam, left=lam[0], right=lam[-1]), lam, row_i - row_i.max())
            v = float(polygamma(1, n_u[b] + 0.5) + polygamma(1, n_s + 0.5))
            half = max(int(np.ceil(4.0 * np.sqrt(v) / dlam)), 1)
            x = np.arange(-half, half + 1) * dlam
            kern = np.exp(-0.5 * x * x / v)
            kern /= kern.sum()
            pr = np.convolve(np.pad(np.exp(r - r.max()), half, mode="edge"), kern, mode="valid")
            r = np.log(np.maximum(pr, 1e-300))
            add += r - r.max()
            live = True
        if live:
            out[int(e)] = add
    return out


def test_the_transfer_delivers_at_pairs_and_licensed_faces_and_relays_nothing(sweep_inputs):
    """The rungs-1+2 contract: rung-1 rows VERBATIM at exactly the independently-derived
    intron|exon pair boundaries; rung-2 rows at exactly the licensed exon faces, equal to the
    independently recomputed map+width transport; zero everywhere else — and the policy's scan
    is None, so nothing can travel a second hop by construction."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    pairs, n_slots = _expected_pairs(sweep_inputs)
    assert pairs, "the toy must carry at least one intron|exon pair or this gate proves nothing"

    ctx = _ctx_of(sweep_inputs)
    pol = TransferPolicy(provider)
    prepared = pol.prepare(ctx)
    assert prepared.scan(backward=False) is None and prepared.scan(backward=True) is None
    import rigel.calibration.messages as M

    msg = prepared.deliver(_dummy_nb(n_slots), _dummy_nb(n_slots))
    assert isinstance(msg, M.PsiMessage) and msg.lam_rows is not None
    rows = np.asarray(msg.lam_rows)
    src = provider(n_grid, window)
    lam, _ = _logodds_grid(n_grid, window)
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, src, lam)
    assert exon_rows, "the toy must license at least one exon face or this gate proves nothing"
    expected_b = {b for b, _j in pairs}
    for i in range(n_slots):
        if i in expected_b:
            j = dict(pairs)[i]
            np.testing.assert_array_equal(rows[i], src[j] - src[j].max(), err_msg=f"slot {i}")
        elif int(i) in exon_rows:
            np.testing.assert_allclose(
                rows[i], exon_rows[int(i)], rtol=0, atol=1e-12, err_msg=f"exon slot {i}"
            )
        else:
            assert not rows[i].any(), f"slot {i} received a transfer it is not licensed for"


def test_the_face_licence_refuses_unmeasured_population_changes():
    """The licence predicate, gated DIRECTLY — the integration gate cannot falsify the terminus
    branch on this toy (no terminus-flagged intron|exon face exists there; watched 2026-09-01),
    so the pure function carries the falsification instead."""
    from rigel.calibration.messages.transfer import face_is_licensed
    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    assert face_is_licensed(0, True, False, True, False)
    for flag in (FLAG_TSS_POS, FLAG_TSS_NEG, FLAG_TES_POS, FLAG_TES_NEG):
        assert not face_is_licensed(flag, True, False, True, False), f"terminus {flag} must refuse"
    assert not face_is_licensed(0, True, False, True, True), "a strand-set change must refuse"
    assert not face_is_licensed(0, True, True, False, True), "a strand-set change must refuse"
    assert face_is_licensed(0, True, True, True, True), "matching AMBIG sides are licensed"


def test_the_face_map_is_monotone_and_flux_capped():
    """The map lam_e(lam_u) is monotone nondecreasing; a measured spliced density CAPS the
    exon's claimable f_g at the closed-form ceiling; s = 0 degenerates to a pure shift by the
    opportunity ratio."""
    from rigel.calibration.messages.transfer import face_map_lambda

    lam = np.linspace(-10, 10, 401)
    n_u, A_g_b, A_r_b, A_g_e, A_r_e, s = 40.0, 200.0, 210.0, 800.0, 790.0, 3.0
    le = face_map_lambda(lam, n_u, A_g_b, A_r_b, A_g_e, A_r_e, s)
    assert np.all(np.diff(le) >= -1e-12), "the map must be monotone nondecreasing"
    ceiling = np.log(n_u / A_g_b * A_g_e) - np.log(s * A_r_e)
    # the ceiling is a SUPREMUM: the grid's last point sits sigma(+L)-short of it (~5e-5 nats
    # at L = 10), so approach is asserted at that scale and the CAP is asserted hard.
    assert le[-1] == pytest.approx(ceiling, abs=1e-3)
    assert le.max() <= ceiling + 1e-9, "certified flux must CAP the claimable f_g"
    le0 = face_map_lambda(lam, n_u, A_g_b, A_r_b, A_g_e, A_r_e, 0.0)
    shift = np.log((A_g_e / A_g_b) / (A_r_e / A_r_b))
    np.testing.assert_allclose(le0, lam + shift, rtol=0, atol=1e-9)


def test_the_ingredient_width_adds_the_counting_variance():
    """A Gaussian factor transported through an identity-shaped map must widen by EXACTLY the
    face's counting variance trigamma(n_u+1/2) + trigamma(n_s+1/2) (the delta-method width) —
    the analytic property an inverted or mis-scaled kernel cannot fake (rung 1's lesson,
    watched 2026-09-01)."""
    from scipy.special import polygamma

    from rigel.calibration.messages.transfer import transport_row

    lam = np.linspace(-10, 10, 401)
    row = -0.5 * (lam - 0.8) ** 2  # unit-variance Gaussian factor
    n_u, n_s = 9.0, 3.0
    v = float(polygamma(1, n_u + 0.5) + polygamma(1, n_s + 0.5))
    out = transport_row(row, lam, lam.copy(), n_u, n_s)  # identity map isolates the width

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        m = w @ lam
        return float(w @ (lam * lam) - m * m)

    assert _var(row) == pytest.approx(1.0, rel=0.02)
    assert _var(out) == pytest.approx(1.0 + v, rel=0.02), "the width must add exactly the variance"
    assert abs(float(out.max())) < 1e-12


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


def _dummy_nb(n_slots):
    from rigel.calibration.messages import NeighbourState

    idx = np.zeros(n_slots, np.int64)
    return NeighbourState(state=(), valid=np.zeros(n_slots, bool), src=idx)


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
                CalibrationConfig(),
                message_propagation=True,
                message_policy="transfer",
                rna_anchor=False,
            ),
            **sweep_inputs["calibrate_kw"],
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert seen and all(isinstance(p, TransferPolicy) for p in seen)
    with pytest.raises(ValueError, match="unknown message_policy"):
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(
                CalibrationConfig(), message_propagation=True, message_policy="no-such-policy"
            ),
            **sweep_inputs["calibrate_kw"],
        )
