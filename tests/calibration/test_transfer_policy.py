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
  stays intact;
* ITEM 1 (the exon → boundary message) and ITEM 2 (the boundary → intron message) each carry a
  deadband-silence gate, an independent recompute of their rows at exactly the licensed slots
  with every other slot untouched, and a directly-gated pure predicate or row law.

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
    is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
    sig = 1.0 / (1.0 + np.exp(-lam))
    dlam = float(lam[1] - lam[0])
    out = {}
    for e in np.flatnonzero(is_exon):
        add = np.zeros(lam.shape[0])
        live = False
        for b in (left[e], right[e]):
            # rung 3: the intergenic|exon EDGE, lower bound only (owner ruling 2026-09-02)
            if b < 0 or not is_bnd[b]:
                continue
            o = left[b] if right[b] == e else right[b]
            if o < 0 or not is_intergenic[o]:
                continue
            if not (n_u[b] > 0 and A_g[b] > 0 and A_g[e] > 0):
                continue  # a zero edge is VACUOUS by the profile (no false claim possible)
            c = sig * n_u[e] * A_g[b] / A_g[e]
            n_b = float(n_u[b])
            r = np.where(c >= n_b, 0.0, n_b * np.log(np.maximum(c, 1e-300) / n_b) - (c - n_b))
            if np.ptp(r) <= 1e-9:
                continue
            add += r - r.max()
            live = True
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
            r = np.interp(
                np.interp(lam, le, lam, left=lam[0], right=lam[-1]), lam, row_i - row_i.max()
            )
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
            out[int(e)] = add - add.max()
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
    item5 = _item5_slots(ctx)  # the terminus boundaries and outside exons: item 5's own gate
    for i in range(n_slots):
        if i in expected_b:
            j = dict(pairs)[i]
            np.testing.assert_array_equal(rows[i], src[j] - src[j].max(), err_msg=f"slot {i}")
        elif int(i) in exon_rows and i not in item5:
            np.testing.assert_allclose(
                rows[i], exon_rows[int(i)], rtol=0, atol=1e-12, err_msg=f"exon slot {i}"
            )
        elif i not in item5:
            assert not rows[i].any(), f"slot {i} received a transfer it is not licensed for"


def test_the_edge_bound_row_is_one_sided_and_vacuous_at_zero():
    """RUNG 3 (owner ruling 2026-09-02: LOWER BOUND ONLY — the upper side is refused as
    over-engineering and the g00-edge residual is an ACCEPTED error): the profile-likelihood
    row sup_{s>=1} Pois(n_b; c/s) is 0 wherever the exon's implied gDNA count c >= n_b (any
    enrichment explains an excess), the edge count's own Poisson tail below it (the sign
    certificate makes the bias direction structural), and IDENTICALLY vacuous at n_b = 0 —
    a zero edge can never manufacture a claim."""
    from rigel.calibration.messages.transfer_rows import edge_bound_row

    lam = np.linspace(-10, 10, 401)
    n_b, n_e, a_g_b, a_g_e = 30.0, 800.0, 200.0, 780.0
    row = edge_bound_row(lam, n_b, n_e, a_g_b, a_g_e)
    sig = 1.0 / (1.0 + np.exp(-lam))
    c = sig * n_e * a_g_b / a_g_e
    assert np.all(row[c >= n_b] == 0.0), "at or above the bound the row must be exactly flat"
    below = c < n_b
    assert np.all(row[below] < 0.0) and np.all(np.diff(row[below]) > 0), (
        "below the bound the penalty must be the one-sided increasing Poisson tail"
    )
    np.testing.assert_allclose(
        row[below], n_b * np.log(c[below] / n_b) - (c[below] - n_b), rtol=0, atol=1e-9
    )
    assert not edge_bound_row(lam, 0.0, n_e, a_g_b, a_g_e).any(), "n_b = 0 must be vacuous"


def test_the_face_licence_refuses_unmeasured_population_changes():
    """The licence predicate, gated DIRECTLY — the integration gate cannot falsify the terminus
    branch on this toy (no terminus-flagged intron|exon face exists there; watched 2026-09-01),
    so the pure function carries the falsification instead."""
    from rigel.calibration.messages.transfer_rows import face_is_licensed
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
    from rigel.calibration.messages.transfer_rows import face_map_lambda

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

    from rigel.calibration.messages.transfer_rows import transport_row

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


# ── ITEM 1 (owner design 2026-09-02): the exon -> intron|exon boundary message ────────────────────


def _expected_splice_out_rows(si, ctx, strand, lam):
    """The item-1 boundary rows recomputed INDEPENDENTLY of the policy: for every licensed
    intron|exon face of every exon whose strand channel is live, the exon's own frozen-variance
    strand row read at the GEOMETRIC splice-in map (both components on ``eff_gdna_global``, the
    spliced density ``S / A_g^b``), marginalised over ``log rho ~ N(log S/U, trigamma(S+1/2) + trigamma(U+1/2))`` on nine
    equal-probability nodes — a second implementation, so a policy bug cannot hide."""
    from scipy.special import polygamma
    from scipy.stats import norm

    from rigel.calibration.splice_graph import (
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    kappa, od_g, od_r = strand
    term = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
    sc_lo = np.asarray(ctx.sj_count_lo, np.float64).sum(axis=1)
    sc_hi = np.asarray(ctx.sj_count_hi, np.float64).sum(axis=1)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    nodes = norm.ppf((np.arange(9) + 0.5) / 9.0)
    out = {}
    for e in np.flatnonzero(is_exon):
        if fp[e] == fn[e] or not tau[e] > 0.0:
            continue
        n = cnt[e].sum()
        ks = kappa if fp[e] else 1.0 - kappa
        f_ref = float(np.clip(belief[e], 1e-3, 1 - 1e-3))
        p = 0.5 * fg + ks * (1 - fg)
        p_ref = 0.5 * f_ref + ks * (1 - f_ref)
        var = max(
            n * p_ref * (1 - p_ref)
            + (n * f_ref) ** 2 * 0.25 * od_g
            + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
            1e-9,
        )
        row_e = -0.5 * (cnt[e, 0] - n * p) ** 2 / var
        row_e -= row_e.max()
        if np.ptp(row_e) <= 1e-9:
            continue
        for b, hi in ((left[e], True), (right[e], False)):
            if b < 0 or not is_bnd[b]:
                continue
            i = left[b] if right[b] == e else right[b]
            if i < 0 or not is_intron[i] or (flags[b] & term) or fp[e] != fp[i] or fn[e] != fn[i]:
                continue
            n_s = float((sc_hi if hi else sc_lo)[b])
            if not (n_u[b] > 0 and A_g[b] > 0 and A_g[e] > 0):
                continue
            v = float(polygamma(1, n_s + 0.5) + polygamma(1, n_u[b] + 0.5))
            acc = np.zeros(lam.shape[0])
            for z in nodes:
                s = n_s / A_g[b] * np.exp(z * np.sqrt(v))
                m = np.log(n_u[b] * fg / A_g[b] * A_g[e]) - np.log(
                    (n_u[b] * (1 - fg) / A_g[b] + s) * A_g[e]
                )
                acc += np.exp(np.interp(m, lam, row_e, left=row_e[0], right=row_e[-1]))
            r = np.log(np.maximum(acc / nodes.size, 1e-300))
            r -= r.max()
            if np.ptp(r) <= 1e-9:
                continue
            out.setdefault(int(b), np.zeros(lam.shape[0]))
            out[int(b)] += r
    return out


def _strand_of(si):
    kw = si["kw"]
    return (
        float(kw["rna_sense_frac"]),
        float(kw.get("gdna_strand_overdispersion", 0.0)),
        float(kw.get("rna_strand_overdispersion", 0.0)),
    )


def test_the_exon_message_is_silence_inside_the_strand_deadband(sweep_inputs):
    """ITEM 1's vacuity law: an exon whose strand channel the solver declares dead (its own
    ``tau_lam`` is 0 — the derived noise-floor deadband, no constant) sends NOTHING, so the
    policy's rows equal the rungs-1–3 rows byte for byte; and a policy built without strand
    parameters is the rungs-1–3 policy exactly."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    base = np.asarray(
        TransferPolicy(provider).prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    dead = _dc.replace(ctx, own=_dc.replace(ctx.own, tau_lam=np.zeros(n)))
    gated = TransferPolicy(provider, strand=_strand_of(sweep_inputs)).prepare(dead)
    np.testing.assert_array_equal(
        np.asarray(gated.deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows), base
    )


def test_the_exon_message_lands_at_licensed_faces_beside_the_intron_row(sweep_inputs):
    """ITEM 1's contract: at every licensed intron|exon face of a strand-live exon the boundary
    receives the exon's splice-out row (independently recomputed) SUMMED with rung 1's intron row —
    two witnesses; every other slot is exactly as rungs 1–3 leave it, except the INTRON slots, which
    item 2 serves off the same ``strand`` parameter and whose contract is item 2's own gate."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    base = np.asarray(
        TransferPolicy(provider).prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    rows = np.asarray(
        TransferPolicy(provider, strand=strand)
        .prepare(ctx)
        .deliver(_dummy_nb(n), _dummy_nb(n))
        .lam_rows
    )
    expected = _expected_splice_out_rows(sweep_inputs, ctx, strand, lam)
    assert expected, (
        "the toy must carry a strand-live exon with a licensed face or this gate proves nothing"
    )
    intron = _intron_mask(ctx)
    item5 = _item5_slots(ctx)
    for i in range(n):
        if i in expected:
            want = base[i] + expected[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"slot {i}"
            )
        elif not intron[i] and i not in item5:
            np.testing.assert_array_equal(
                rows[i], base[i], err_msg=f"slot {i} moved without a licence"
            )


def test_the_splice_out_row_is_the_count_form_widened_by_the_marginal():
    """The splice-out row's three analytic properties: (i) with equal face and region geometry its
    mode sits at f_b = f_E (U+S)/U — the owner's arithmetic with the enrichment ratio cancelled;
    (ii) more gDNA at the boundary than in the exon whenever S > 0 (the reversed subtraction cannot
    pass); (iii) the marginal over log rho is WIDER at a thin face than at a deep one."""
    from rigel.calibration.messages.transfer_rows import splice_out_row

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    f_e = 0.20
    row_e = -0.5 * ((lam - np.log(f_e / (1 - f_e))) / 0.05) ** 2  # a sharp exon belief at f_E

    def _mode_f(r):
        return float(sig[np.argmax(r)])

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        m = w @ lam
        return float(w @ (lam * lam) - m * m)

    deep = splice_out_row(row_e, lam, n_u=2000.0, n_s=3000.0, a_g_b=210.0, a_g_e=210.0)
    assert _mode_f(deep) == pytest.approx(f_e * (2000 + 3000) / 2000, abs=0.01)
    assert _mode_f(deep) > f_e
    thin = splice_out_row(row_e, lam, n_u=6.0, n_s=9.0, a_g_b=210.0, a_g_e=210.0)
    assert _var(thin) > 3 * _var(deep), "a thin face must deliver a much wider claim"
    assert abs(float(deep.max())) < 1e-12 and abs(float(thin.max())) < 1e-12
    assert not splice_out_row(row_e, lam, n_u=0.0, n_s=9.0, a_g_b=210.0, a_g_e=210.0).any()


# ── ITEM 2 (owner design 2026-09-02): the intron|exon boundary -> intron message ──────────────────


def _expected_boundary_rows(si, ctx, strand, lam):
    """The item-2 intron rows recomputed INDEPENDENTLY of the policy: for every intron|exon pair whose
    boundary admits the intron's SINGLE strand set and whose strand channel the solver declares live
    (``own.tau_lam > 0`` at the boundary — the strand Fisher information itself, there being no
    factory at a boundary), the boundary's own frozen-variance strand row VERBATIM, summed at the
    intron — a second implementation, so a policy bug cannot hide."""
    kappa, od_g, od_r = strand
    pairs, _n = _expected_pairs(si)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    out = {}
    for b, i in pairs:
        if fp[b] != fp[i] or fn[b] != fn[i] or fp[b] == fn[b] or not tau[b] > 0.0:
            continue
        n = cnt[b].sum()
        ks = kappa if fp[b] else 1.0 - kappa
        f_ref = float(np.clip(belief[b], 1e-9, 1 - 1e-9))
        p = 0.5 * fg + ks * (1 - fg)
        p_ref = 0.5 * f_ref + ks * (1 - f_ref)
        var = max(
            n * p_ref * (1 - p_ref)
            + (n * f_ref) ** 2 * 0.25 * od_g
            + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
            1e-9,
        )
        row = -0.5 * (cnt[b, 0] - n * p) ** 2 / var
        row -= row.max()
        if np.ptp(row) <= 1e-9:
            continue
        out.setdefault(int(i), np.zeros(lam.shape[0]))
        out[int(i)] += row
    return out


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


def test_the_boundary_message_is_silence_inside_the_strand_deadband(sweep_inputs):
    """ITEM 2's vacuity law: a boundary whose strand channel the solver declares dead (its own
    ``tau_lam`` is 0 — at a boundary that IS the strand Fisher information, no constant) sends
    NOTHING to its intron. With every boundary dead and every exon live, no intron receives a row and
    every other slot is exactly as the fully-live policy leaves it (item 1 reads the EXON's deadband,
    which is intact). ⚠ This gate cannot fail before the mechanism exists — nothing delivered to an
    intron before item 2 — so its falsification is the watched perturbation (the deadband dropped)."""
    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    pol = TransferPolicy(provider, strand=_strand_of(sweep_inputs))
    live = np.asarray(pol.prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    dead = np.asarray(
        pol.prepare(_dead_boundaries(ctx)).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    intron = _intron_mask(ctx)
    keep = ~intron
    keep[sorted(_item5_slots(ctx))] = False  # a dead boundary also stops item 5's row to its exon
    assert not dead[intron].any(), "a dead boundary delivered to its intron"
    np.testing.assert_array_equal(dead[keep], live[keep])


def test_the_boundary_message_lands_at_introns_beside_nothing_else(sweep_inputs):
    """ITEM 2's contract: every intron whose licensed, strand-live boundaries send receives the SUM of
    their own strand rows (independently recomputed — two witnesses at a two-faced intron); an intron
    with no such boundary receives nothing; every other slot is exactly as rungs 1–3 + item 1 leave
    it (the boundaries-dead rows, which the vacuity gate proves item-2-free)."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    pol = TransferPolicy(provider, strand=strand)
    rows = np.asarray(pol.prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    base = np.asarray(
        pol.prepare(_dead_boundaries(ctx)).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    expected = _expected_boundary_rows(sweep_inputs, ctx, strand, lam)
    assert expected, "the toy must carry a strand-live intron|exon pair or this gate proves nothing"
    intron = _intron_mask(ctx)
    item5 = _item5_slots(ctx)
    for i in range(n):
        if i in expected:
            want = expected[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"intron slot {i}"
            )
        elif intron[i]:
            assert not rows[i].any(), f"intron slot {i} received a row without a licensed live face"
        elif i not in item5:
            np.testing.assert_array_equal(rows[i], base[i], err_msg=f"slot {i} moved with item 2")


def test_the_boundary_strand_licence_requires_one_shared_strand():
    """The pair predicate, gated DIRECTLY — the integration toy cannot falsify it (every pair there
    shares one strand), so the pure function carries the falsification: one strand, shared by both
    flanks, is licensed; an AMBIG pair (whose strand split constrains the tilt, never f_g), a
    strand-set change and an empty boundary refuse."""
    from rigel.calibration.messages.transfer_rows import boundary_shares_strand

    assert boundary_shares_strand(True, False, True, False)
    assert boundary_shares_strand(False, True, False, True)
    assert not boundary_shares_strand(True, True, True, True), "an AMBIG pair must refuse"
    assert not boundary_shares_strand(True, False, True, True), "a strand-set change must refuse"
    assert not boundary_shares_strand(False, False, True, False), "an empty boundary must refuse"


# ── ITEM 5 (rung 4, 2026-09-02): the exon|exon TERMINUS boundary and its OUTSIDE exon ────────────


def _expected_terminus_rows(si, ctx, strand, lam, exon_rows):
    """Item 5 recomputed INDEPENDENTLY of the policy. At every exon|exon boundary carrying exactly one
    terminus direction and no splice junction, the OUTSIDE flank is read off the flag alone (TSS+ and
    TES− bodies extend right, outside = LEFT; TES+ and TSS− extend left, outside = RIGHT). The boundary
    receives (i) the rungs-2/3 rows of its outside exon (``exon_rows``, the independent recompute) and
    (ii) the exon's own frozen-variance strand row when its channel is live, each through the
    splice-out map with S = the boundary's SPLICED crossing; the outside exon receives the boundary's
    own strand row through the face map with that spliced density, widened by the counting variance.
    Returns ``(rows_at_boundaries, rows_at_outside_exons)`` keyed by slot."""
    from rigel.calibration.messages.transfer_rows import (
        face_map_lambda,
        splice_out_row,
        transport_row,
    )
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    kappa, od_g, od_r = strand
    sj = FLAG_DONOR_POS | FLAG_DONOR_NEG | FLAG_ACCEPTOR_POS | FLAG_ACCEPTOR_NEG
    body_right, body_left = FLAG_TSS_POS | FLAG_TES_NEG, FLAG_TES_POS | FLAG_TSS_NEG
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
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

    at_b, at_o = {}, {}
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        f = int(flags[b])
        to_right, to_left = bool(f & body_right), bool(f & body_left)
        if (f & sj) or to_right == to_left:
            continue
        o = lo if to_right else hi
        if fp[b] != fp[o] or fn[b] != fn[o] or fp[b] == fn[b]:
            continue
        if not (n_u[b] > 0 and A_g[b] > 0 and A_g[o] > 0):
            continue
        acc = np.zeros(lam.shape[0])
        if int(o) in exon_rows:
            acc += splice_out_row(exon_rows[int(o)], lam, n_u[b], n_s[b], A_g[b], A_g[o])
        if tau[o] > 0.0:
            acc += splice_out_row(strand_row(o), lam, n_u[b], n_s[b], A_g[b], A_g[o])
        if np.ptp(acc) > 1e-9:
            at_b[int(b)] = acc
        if tau[b] > 0.0:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            at_o.setdefault(int(o), np.zeros(lam.shape[0]))
            at_o[int(o)] += transport_row(strand_row(b), lam, le, n_u[b], n_s[b])
    return at_b, at_o


def _item5_slots(ctx):
    """Item 5's delivery sites — every exon|exon boundary carrying one terminus direction and no sj
    whose outside flank shares its single strand, and that outside exon — derived from the flags
    alone, so the earlier gates can hand these slots to item 5's own gate."""
    from rigel.calibration.messages.transfer_rows import boundary_shares_strand, outside_flank

    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    out = set()
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        o, _i = outside_flank(flags[b], lo, hi)
        if o is not None and boundary_shares_strand(fp[b], fn[b], fp[o], fn[o]):
            out |= {int(b), int(o)}
    return out


def _terminus_flags_cleared(ctx):
    """The context with every exon|exon boundary's terminus bits cleared — item 5 sees no terminus,
    every other message is unchanged (they read flags at intron|exon faces only)."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer_rows import TERMINUS

    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16).copy()
    for b in np.flatnonzero(is_bnd):
        if left[b] >= 0 and right[b] >= 0 and is_exon[left[b]] and is_exon[right[b]]:
            flags[b] &= ~TERMINUS
    return _dc.replace(ctx, boundary_flags=flags)


def test_the_terminus_orientation_reads_the_flag_alone():
    """The orientation table, gated DIRECTLY: a + start or a − end extends genomic-right, so the
    OUTSIDE flank is the left one; a + end or a − start extends left, so it is the right one; termini
    pointing both ways, no terminus at all, or a splice junction sharing the boundary give no side."""
    from rigel.calibration.messages.transfer_rows import outside_flank
    from rigel.calibration.splice_graph import (
        FLAG_DONOR_POS,
        FLAG_TES_NEG,
        FLAG_TES_POS,
        FLAG_TSS_NEG,
        FLAG_TSS_POS,
    )

    assert outside_flank(FLAG_TSS_POS, 7, 9) == (7, 9)
    assert outside_flank(FLAG_TES_NEG, 7, 9) == (7, 9)
    assert outside_flank(FLAG_TES_POS, 7, 9) == (9, 7)
    assert outside_flank(FLAG_TSS_NEG, 7, 9) == (9, 7)
    assert outside_flank(FLAG_TSS_POS | FLAG_TES_POS, 7, 9) == (None, None), "both ways: no side"
    assert outside_flank(0, 7, 9) == (None, None), "no terminus: no side"
    assert outside_flank(FLAG_TSS_POS | FLAG_DONOR_POS, 7, 9) == (None, None), (
        "sj+terminus: its own item"
    )


def test_the_terminus_messages_land_at_the_outside_pair_beside_nothing_else(sweep_inputs):
    """ITEM 5's contract on the toy's two exon|exon terminus boundaries (a + start and a + end): the
    boundary receives exactly the independently recomputed composed transport plus the outside
    exon's mapped strand row; the outside exon receives exactly the boundary's mapped strand row on
    top of its rungs-2/3 rows; the INSIDE exon and every other slot are exactly as the policy leaves
    them when the terminus bits are cleared — so a reversed orientation, a dropped spliced crossing,
    or a delivery to the inside flank fails (each watched firing, 2026-09-02). ⚠ An ECHO — the
    composed transport reading a row item 5 itself delivered to the exon — is NOT catchable on this
    toy (each outside exon serves one boundary and is read before anything reaches it; watched
    silent), so the policy keeps that law STRUCTURALLY: item 5 accumulates apart and reads only the
    rows the earlier messages left."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _ctx_of(sweep_inputs)
    n = int(ctx.n_slots)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    pol = TransferPolicy(provider, strand=strand)
    rows = np.asarray(pol.prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    base = np.asarray(
        pol.prepare(_terminus_flags_cleared(ctx)).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, provider(n_grid, window), lam)
    at_b, at_o = _expected_terminus_rows(sweep_inputs, ctx, strand, lam, exon_rows)
    assert len(at_b) >= 2 and at_o, (
        "the toy must carry two served terminus boundaries or this gate proves nothing"
    )
    for i in range(n):
        if i in at_b:
            assert not base[i].any(), f"boundary {i} receives something other than item 5"
            want = at_b[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"terminus boundary {i}"
            )
        elif i in at_o:
            want = base[i] + at_o[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"outside exon {i}"
            )
        else:
            np.testing.assert_array_equal(rows[i], base[i], err_msg=f"slot {i} moved with item 5")


# ── ITEM 6 (owner design 2026-09-02): the abundance-discrepancy message into the inside flank ────


def _expected_abundance_rows(si, ctx, strand, lam):
    """Item 6 recomputed INDEPENDENTLY of the policy: at every served terminus boundary (one direction,
    no sj, exon flanks, a live strand channel), the boundary's own strand row carried into the INSIDE
    flank through the abundance map, the step's spread refit here from the served pairs' two modes,
    the hard cap and the delta-method width — keyed by the inside slot."""
    from rigel.calibration.messages.transfer_rows import (
        abundance_row,
        boundary_shares_strand,
        outside_flank,
    )

    kappa, od_g, od_r = strand
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
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
        _o, i = outside_flank(flags[b], lo, hi)
        if i is None or not boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
            continue
        if not (tau[b] > 0.0 and n_u[b] > 0 and n_u[i] > 0 and A_g[b] > 0 and A_g[i] > 0):
            continue
        served.append((int(b), int(i), (n_u[i] / A_g[i]) / ((n_u[b] + n_s[b]) / A_g[b])))
    ls, vs = [], []
    for b, i, r in served:
        if not tau[i] > 0.0:
            continue
        ks = kappa if fp[b] else 1.0 - kappa
        fs, vv = [], 0.0
        for x in (b, i):
            n = cnt[x].sum()
            p = cnt[x, 0] / n
            f = (p - ks) / (0.5 - ks)
            fs.append(f)
            vv += p * (1 - p) / n / (p - ks) ** 2
        if not all(0.0 < f < 1.0 for f in fs):
            continue
        ls.append(np.log(fs[1] / (fs[0] * n_u[b] / (n_u[b] + n_s[b]) / r)))
        vs.append(vv + 1.0 / n_u[i] + 1.0 / (n_u[b] + n_s[b]))
    v_step = 0.0
    if len(ls) >= 2:
        ls, vs = np.asarray(ls), np.asarray(vs)
        w = 1.0 / vs
        w /= w.sum()
        mu = float(w @ ls)
        v_step = max(0.0, float(w @ (ls - mu) ** 2) - float(w @ vs))
    out = {}
    for b, i, r in served:
        row = abundance_row(strand_row(b), lam, n_u[b], n_s[b], r, n_u[i], v_step)
        if np.ptp(row) > 1e-9:
            out.setdefault(i, np.zeros(lam.shape[0]))
            out[i] += row
    return out, v_step


def _with_populated_inside(ctx):
    """The toy's inside pieces are EMPTY (shorter than a fragment: no contained fragment, a contained
    opportunity below one base), so item 6 serves nothing there and a gate on the bare toy proves
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
    a_g = np.asarray(ctx.eff_gdna_global, np.float64).copy()
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
        eff_gdna_global=a_g,
        eff_rna=a_r,
        own=_dc.replace(ctx.own, tau_lam=tau),
    )


def test_the_abundance_map_holds_the_two_hypotheses_and_the_cap():
    """The map's three analytic properties: at s = r the flank's share is the composition-transfer
    value f_c = f_b U/(U+S) (enrichment); at s = 1 it is f_c / r (new RNA); it is monotone in the
    boundary's log-odds; and the row never claims a share above f_c — a source peaked at f_b arrives
    at or below the composition-transfer value whatever the step prior, the hard cap."""
    from rigel.calibration.messages.transfer_rows import abundance_map, abundance_row

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    U, S, r = 300.0, 60.0, 2.5
    fc = U * sig / (U + S)
    m_enrich = abundance_map(lam, U, S, r, r)
    m_new = abundance_map(lam, U, S, r, 1.0)
    np.testing.assert_allclose(1 / (1 + np.exp(-m_enrich)), fc, rtol=0, atol=1e-9)
    np.testing.assert_allclose(1 / (1 + np.exp(-m_new)), fc / r, rtol=0, atol=1e-9)
    assert np.all(np.diff(m_enrich) >= -1e-12) and np.all(np.diff(m_new) >= -1e-12)

    def _var(rw):
        p = np.exp(rw - rw.max())
        p /= p.sum()
        m = p @ lam
        return float(p @ (lam * lam) - m * m)

    f_b = 0.6
    row_b = -0.5 * ((lam - np.log(f_b / (1 - f_b))) / 0.05) ** 2
    # the cap is exact in the map and softened only by the ingredients' counting width, so it is
    # asserted on DEEP counts (the blur's sd ~0.03 nats there); the shallow case is asserted by mode
    Ud, Sd, nd = 3000.0, 600.0, 4000.0
    fc_b = f_b * Ud / (Ud + Sd)
    for v in (0.0, 0.3):
        row = abundance_row(row_b, lam, Ud, Sd, r, nd, v)
        w = np.exp(row - row.max())
        assert float(sig[np.argmax(row)]) <= fc_b + 0.02, (
            "the cap: never above the composition-transfer value"
        )
        assert w[sig > fc_b + 0.05].max() < 1e-3, "no mass above the cap"
    sharp = abundance_row(row_b, lam, Ud, Sd, r, nd, 0.0)
    wide = abundance_row(row_b, lam, Ud, Sd, r, nd, 0.3)
    assert abs(float(sig[np.argmax(sharp)]) - fc_b / r) < 0.02, (
        "v_step = 0 is the new-RNA point when the totals rise"
    )
    down = abundance_row(row_b, lam, Ud, Sd, 0.4, nd, 0.0)
    assert abs(float(sig[np.argmax(down)]) - fc_b) < 0.02, (
        "v_step = 0 is the composition point when the totals fall"
    )
    down_wide = abundance_row(row_b, lam, Ud, Sd, 0.4, nd, 0.3)
    assert float(sig[np.argmax(down_wide)]) <= fc_b + 0.02 and _var(down_wide) > _var(down), (
        "below the cap, and wider with a fitted spread"
    )
    assert _var(wide) > _var(sharp), "a fitted spread widens the claim"
    assert not abundance_row(row_b, lam, 0.0, S, r, 400.0, 0.0).any(), (
        "a depleted boundary is vacuous"
    )


def test_the_abundance_message_lands_at_inside_exons_beside_their_face_rows(sweep_inputs):
    """ITEM 6's contract on the toy's two terminus boundaries: each INSIDE exon receives exactly the
    independently recomputed abundance row (the step's spread refit here) on top of its rungs-2/3
    face rows; every slot that is neither an item-5 site nor an inside exon is exactly as the policy
    leaves it with the terminus bits cleared — so a dropped cap, a dropped fit, an inverted ratio or a
    delivery to the outside flank fails."""
    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _with_populated_inside(_ctx_of(sweep_inputs))
    n = int(ctx.n_slots)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    pol = TransferPolicy(provider, strand=strand)
    rows = np.asarray(pol.prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    base = np.asarray(
        pol.prepare(_terminus_flags_cleared(ctx)).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, provider(n_grid, window), lam)
    expected, v_step = _expected_abundance_rows(sweep_inputs, ctx, strand, lam)
    assert len(expected) >= 2 and v_step > 0.0, (
        "the toy must carry two served inside exons or this gate proves nothing"
    )
    item5 = _item5_slots(ctx)
    for i in range(n):
        if i in expected:
            want = exon_rows.get(i, np.zeros(lam.shape[0])) + expected[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"inside exon {i}"
            )
        elif i not in item5:
            np.testing.assert_array_equal(rows[i], base[i], err_msg=f"slot {i} moved with item 6")


def test_the_abundance_message_is_silence_inside_the_strand_deadband(sweep_inputs):
    """ITEM 6's vacuity law: a terminus boundary whose strand channel the solver declares dead sends
    nothing into its inside exon — with every boundary dead, each inside exon's row equals its row
    under cleared terminus bits (rungs 2/3 alone). ⚠ Cannot fail before the mechanism exists; its
    falsification is the watched perturbation (the deadband dropped)."""
    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _with_populated_inside(_ctx_of(sweep_inputs))
    n = int(ctx.n_slots)
    pol = TransferPolicy(provider, strand=_strand_of(sweep_inputs))
    dead = np.asarray(
        pol.prepare(_dead_boundaries(ctx)).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows
    )
    cleared = np.asarray(
        pol.prepare(_terminus_flags_cleared(_dead_boundaries(ctx)))
        .deliver(_dummy_nb(n), _dummy_nb(n))
        .lam_rows
    )
    from rigel.calibration.simplex_logodds import _logodds_grid

    lam, _ = _logodds_grid(n_grid, window)
    expected, _v = _expected_abundance_rows(sweep_inputs, ctx, _strand_of(sweep_inputs), lam)
    for i in expected:
        np.testing.assert_array_equal(
            dead[i], cleared[i], err_msg=f"inside exon {i} received a row from a dead boundary"
        )


# ── ITEM 7 (2026-09-02): the alternative splice site ──────────────────────────────────────────────


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
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
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


def test_the_alt_splice_messages_land_at_the_junction_and_both_flanks(sweep_inputs):
    """ITEM 7's contract on the patched toy (one donor, one acceptor, both flanks populated): the
    junction boundary receives exactly both flanks' mapped rows and each flank exactly the boundary's
    mapped row — each blurred by the pair's own disagreement beyond counting, the fills chosen to
    DISAGREE so at least one width is non-trivial — on top of its rungs-2/3 rows; every other slot is
    exactly as the policy leaves it with the junction bits cleared. So swapped flanks, a dropped flux,
    a dropped width, or a delivery to the wrong flank fails."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.messages.transfer_rows import SJ_FLAGS
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _with_alt_splice_sites(_ctx_of(sweep_inputs))
    n = int(ctx.n_slots)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    pol = TransferPolicy(provider, strand=strand)
    rows = np.asarray(pol.prepare(ctx).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    flags = np.asarray(ctx.boundary_flags, np.uint16).copy()
    flags[np.asarray(ctx.is_boundary, bool)] &= ~SJ_FLAGS
    base = np.asarray(
        pol.prepare(_dc.replace(ctx, boundary_flags=flags))
        .deliver(_dummy_nb(n), _dummy_nb(n))
        .lam_rows
    )
    expected, widths = _expected_alt_splice_rows(sweep_inputs, ctx, strand, lam)
    assert len(expected) >= 4, (
        "the patched toy must serve two junctions and their flanks or this gate proves nothing"
    )
    assert any(w > 0.0 for w in widths.values()), (
        "the fills must disagree beyond counting somewhere or the width is untested"
    )
    for i in range(n):
        if i in expected:
            want = base[i] + expected[i]
            np.testing.assert_allclose(
                rows[i], want - want.max(), rtol=0, atol=1e-10, err_msg=f"slot {i}"
            )
        else:
            np.testing.assert_array_equal(rows[i], base[i], err_msg=f"slot {i} moved with item 7")


def test_the_alt_splice_messages_are_silence_inside_the_strand_deadband(sweep_inputs):
    """ITEM 7's vacuity law: with every boundary's strand channel dead and the flanks' too, the junction
    bits change nothing — each slot's row equals its row with the junction bits cleared. ⚠ Cannot fail
    before the mechanism exists; its falsification is the watched perturbation (a deadband dropped)."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.calibration.messages.transfer_rows import SJ_FLAGS

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    ctx = _with_alt_splice_sites(_ctx_of(sweep_inputs))
    n = int(ctx.n_slots)
    tau = np.zeros(n)
    dead = _dc.replace(ctx, own=_dc.replace(ctx.own, tau_lam=tau))
    pol = TransferPolicy(provider, strand=_strand_of(sweep_inputs))
    rows = np.asarray(pol.prepare(dead).deliver(_dummy_nb(n), _dummy_nb(n)).lam_rows)
    flags = np.asarray(dead.boundary_flags, np.uint16).copy()
    flags[np.asarray(dead.is_boundary, bool)] &= ~SJ_FLAGS
    cleared = np.asarray(
        pol.prepare(_dc.replace(dead, boundary_flags=flags))
        .deliver(_dummy_nb(n), _dummy_nb(n))
        .lam_rows
    )
    np.testing.assert_array_equal(rows, cleared)
