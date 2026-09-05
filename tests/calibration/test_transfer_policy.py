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
from scipy.special import polygamma

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """ONE real `solve_chain` call captured from a calibrate run on the toy — the backbone-parity
    pattern: every gate below re-runs the sweep with a different policy on byte-identical inputs.
    (Moved here from the retired foundation-spec gate file, 2026-09-04.)"""
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


def _nothing_held(n_slots):
    """The two held lists of a policy that sent nothing: SILENCE at every interior node."""
    from rigel.calibration.messages import SILENCE

    return [SILENCE] * int(n_slots)


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


# ── THE LEVEL RULE (MESSAGE_PLAN.md step A, 2026-09-04): the terminus boundary into the region inside it ──


def _expected_level_rows(si, ctx, strand, lam):
    """THE LEVEL RULE recomputed INDEPENDENTLY of the policy: at every terminus boundary with an inside
    EXON (the outside an exon or an intron), the boundary's own strand row through the level-kept map
    (the boundary's share times its crossing density, times the inside's opportunity, over the
    inside's total), blurred by both totals' counting plus the pair's discrepancies — the totals'
    disagreement beyond counting and, where both strand channels are live, the two strand modes'
    disagreement beyond counting. With no own row the crossing total's one-sided upper bound. Keyed by
    the inside slot; ``(rows, served pairs)``."""
    from rigel.calibration.messages.transfer_rows import (
        boundary_shares_strand,
        level_bound_row,
        level_map_lambda,
        level_row,
        outside_flank,
    )

    kappa, od_g, od_r = strand
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    is_intron = ~is_bnd & ~is_exon & (fp | fn)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    out, served = {}, []
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0:
            continue
        o, i = outside_flank(flags[b], lo, hi)
        if i is None or not is_exon[i] or not (is_exon[o] or is_intron[o]):
            continue
        if not boundary_shares_strand(fp[b], fn[b], fp[i], fn[i]):
            continue
        if not (n_u[b] > 0 and n_u[i] > 0 and A_g[b] > 0 and A_g[i] > 0):
            continue
        served.append((int(b), int(i), "exon|exon" if is_exon[o] else "exon|intron"))
        d_b, T_b = n_u[b] / A_g[b], n_u[b] + n_s[b]
        r = (n_u[i] / A_g[i]) / (T_b / A_g[b])
        v = max(0.0, np.log(r) ** 2 - (1.0 / n_u[i] + 1.0 / T_b))
        if tau[b] > 0.0 and tau[i] > 0.0:
            ok, modes = True, []
            for y in (b, i):
                n = cnt[y].sum()
                p = cnt[y, 0] / n
                ks = kappa if fp[y] else 1.0 - kappa
                f = (p - ks) / (0.5 - ks)
                if not 0.0 < f < 1.0:
                    ok = False
                    break
                modes.append((f, p * (1 - p) / n / (p - ks) ** 2 / (1 - f) ** 2))
            if ok:
                (f_b, v_b), (f_i, v_i) = modes
                f_pred = min(f_b / r, 1 - 1e-9)
                dd = np.log(f_i / (1 - f_i)) - np.log(f_pred / (1 - f_pred))
                v += max(0.0, dd * dd - (v_b + v_i + 1.0 / n_u[i] + 1.0 / T_b))
        v += float(polygamma(1, n_u[b] + 0.5) + polygamma(1, n_u[i] + 0.5))
        m = level_map_lambda(lam, d_b, A_g[i], n_u[i])
        if tau[b] > 0.0 and fp[b] != fn[b]:
            row = level_row(_strand_row_of(ctx, strand, lam, b), lam, m, v)
        else:
            row = level_bound_row(lam, d_b, A_g[i], n_u[i], v)
        if np.ptp(row) > 1e-9:
            out.setdefault(int(i), np.zeros(lam.shape[0]))
            out[int(i)] += row
    return out, served


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


def test_the_level_map_keeps_the_level_and_the_bound_is_vacuous_below_the_total():
    """THE LEVEL RULE's arithmetic, gated directly: through the level-kept map a boundary whose share is
    f_b and whose crossing density is d lands the inside at f_b * d * E_i / T_i (the gDNA level kept,
    the inside's own total supplying the rest); the map is monotone; a profile read through it keeps
    its peak there and its width grows with the dampening; the total's upper bound charges nothing
    below the crossing's density and charges above it."""
    from rigel.calibration.messages.transfer_rows import (
        level_bound_row,
        level_map_lambda,
        level_row,
    )

    lam = np.linspace(-8, 8, 801)
    sig = 1 / (1 + np.exp(-lam))
    d_b, E_i, T_i = (
        0.2,
        500.0,
        200.0,
    )  # the boundary crosses 0.2 fragments/base; the inside holds 200
    m = level_map_lambda(lam, d_b, E_i, T_i)
    assert np.all(np.diff(m) >= 0.0)
    f_b = 0.4
    row_b = -0.5 * ((lam - np.log(f_b / (1 - f_b))) / 0.05) ** 2
    want = f_b * d_b * E_i / T_i  # 0.2: the level kept
    tight = level_row(row_b, lam, m, 0.001)
    assert float(sig[np.argmax(tight)]) == pytest.approx(want, abs=0.01)

    def _var(r):
        w = np.exp(r - r.max())
        w /= w.sum()
        mu = w @ lam
        return float(w @ (lam * lam) - mu * mu)

    wide = level_row(row_b, lam, m, 0.5)
    assert _var(wide) > 3 * _var(tight), "the dampening must widen the delivered profile"
    assert not level_row(np.zeros_like(lam), lam, m, 0.1).any(), "a flat profile is vacuous"
    ub = level_bound_row(lam, d_b, E_i, T_i, 0.01)
    below = sig * T_i / E_i < d_b  # inside gDNA density below the crossing's total density
    assert np.all(ub[below] == 0.0) and np.all(ub[~below] <= 0.0) and np.any(ub[~below] < 0.0)


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


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE PASS FORM (owner ruling 2026-09-04, `DESIGN.md` §6b.12): every node's OWN CLAIM, a RULE per
# directed face, and the two passes. Each family below is gated on its claim and its rule directly
# (the independent recomputations above supply the expected rows); the passes are gated once, against
# an independent recursive reference, and once for the no-echo law.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _drive_the_backbone(prepared, ctx):
    """The backbone's own contract, reproduced (`sweep.solve_chain`'s two directional passes and its
    solve): each pass calls ``receive(source, destination)`` over the chain order — the forward pass
    reading each slot's LOW neighbour, the backward pass its HIGH one — and holds the result at the
    destination; ``solve`` receives the two held lists. Returns the delivered rows (zeros when the
    policy is silent)."""
    from rigel.calibration.messages import SILENCE

    order = list(ctx.order)
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

    n = len(own)
    rows = np.zeros((n, int(ctx.n_grid)))
    for i in range(n):
        parts = [into(s, i) for s in (left[i], right[i]) if s >= 0]
        parts = [p for p in parts if p is not None]
        if parts:
            rows[i] = norm(sum(parts))
    return rows


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


def test_the_intron_face_carries_the_pair_identity_and_the_face_map(sweep_inputs):
    """RUNGS 1–3 as claims and rules: an intron's own claim is its factory profile and its rule into
    each boundary is the identity (FORWARD); the rule from a licensed face into the exon, applied to
    the intron's claim and summed with the edge's level rule, equals the independently recomputed
    rung-2 + rung-3 rows of that exon; an unlicensed face has no rule into the exon."""
    from rigel.calibration.messages.transfer import TransferPolicy, _forward
    from rigel.calibration.simplex_logodds import _logodds_grid

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    pairs, _n = _expected_pairs(sweep_inputs)
    assert pairs, "the toy must carry at least one intron|exon pair or this gate proves nothing"
    ctx = _ctx_of(sweep_inputs)
    prepared = TransferPolicy(provider).prepare(ctx)
    src = provider(n_grid, window)
    lam, _ = _logodds_grid(n_grid, window)
    for b, j in pairs:
        assert prepared.rule.get((int(j), int(b))) is _forward, f"pair ({j}, {b}) is not FORWARD"
        np.testing.assert_array_equal(prepared.own[j], src[j] - src[j].max())
    exon_rows = _expected_exon_rows(sweep_inputs, ctx, src, lam)
    assert exon_rows, "the toy must license at least one exon face or this gate proves nothing"
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    intron = _intron_mask(ctx)
    for e in np.flatnonzero(is_exon):
        acc = None
        for b in (left[e], right[e]):
            fn = prepared.rule.get((int(b), int(e)))
            if fn is None:
                continue
            i = left[b] if right[b] == e else right[b]
            claim = prepared.own[i] if (i >= 0 and intron[i]) else prepared.own[b]
            if claim is None:
                continue
            r = fn(claim, None)
            if r is not None and np.ptp(r) > 1e-9:
                acc = (r - r.max()) if acc is None else acc + (r - r.max())
        if int(e) in exon_rows:
            assert acc is not None, f"exon {e} has no rule where the recompute expects rows"
            np.testing.assert_allclose(acc - acc.max(), exon_rows[int(e)], rtol=0, atol=1e-12)
        else:
            assert acc is None, f"exon {e} received a rule it is not licensed for"


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
    # boundaries dead, exons live: item 2's claim is absent at every boundary
    dead_b = pol.prepare(_dead_boundaries(ctx))
    for b in np.flatnonzero(is_bnd):
        assert dead_b.own[b] is None or not dead_b.own[b].any()


def test_the_exon_and_boundary_own_claims_are_the_strand_rows_and_their_rules_the_maps(
    sweep_inputs,
):
    """ITEMS 1 and 2 as claims and rules: at every licensed face of a strand-live exon the rule
    exon → boundary applied to the exon's claim, summed per boundary, equals the independently
    recomputed splice-out rows; at every intron|exon pair sharing one strand with a live boundary,
    the rule boundary → intron is the identity and the boundary's claim, summed per intron, equals
    the independently recomputed strand rows."""
    from rigel.calibration.messages.transfer import _forward
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    expected = _expected_splice_out_rows(sweep_inputs, ctx, strand, lam)
    assert expected, (
        "the toy must carry a strand-live exon with a licensed face or this gate proves nothing"
    )
    intron = _intron_mask(ctx)
    got = {}
    for e in np.flatnonzero(is_exon):
        if prepared.own[e] is None:
            continue
        for b in (left[e], right[e]):
            fn = prepared.rule.get((int(e), int(b)))
            if fn is None or b < 0:
                continue
            other = left[b] if right[b] == e else right[b]
            if other < 0 or not intron[other]:
                continue  # item 5 / item 7 rules leave through exon|exon faces: their own gates
            r = fn(prepared.own[e], None)
            if r is not None and np.ptp(r) > 1e-9:
                got[int(b)] = got.get(int(b), 0.0) + (r - r.max())
    assert set(got) == set(expected), set(got) ^ set(expected)
    for b in expected:
        np.testing.assert_allclose(
            got[b] - got[b].max(), expected[b] - expected[b].max(), rtol=0, atol=1e-10
        )
    expected_i = _expected_boundary_rows(sweep_inputs, ctx, strand, lam)
    assert expected_i, (
        "the toy must carry a strand-live intron|exon pair or this gate proves nothing"
    )
    pairs, _n = _expected_pairs(sweep_inputs)
    got_i = {}
    for b, i in pairs:
        fn = prepared.rule.get((int(b), int(i)))
        if fn is None or prepared.own[b] is None:
            continue
        assert fn is _forward, f"the boundary → intron rule at ({b}, {i}) is not the identity"
        got_i[int(i)] = got_i.get(int(i), 0.0) + prepared.own[b]
    assert set(got_i) == set(expected_i), set(got_i) ^ set(expected_i)
    for i in expected_i:
        np.testing.assert_allclose(
            got_i[i] - got_i[i].max(), expected_i[i] - expected_i[i].max(), rtol=0, atol=1e-10
        )


def test_the_terminus_rules_land_at_the_outside_pair_and_nowhere_when_the_flags_clear(sweep_inputs):
    """ITEM 5 as claims and rules: at every exon|exon terminus boundary the rule outside exon → boundary
    applied to the exon's claim is the independently recomputed splice-out row with the SPLICED
    crossing, the rule boundary → outside exon applied to the boundary's claim the recomputed face-map
    row; the rules exist exactly at the outside pairs and vanish when the terminus bits are cleared."""
    from rigel.calibration.messages.transfer_rows import (
        face_map_lambda,
        splice_out_row,
        transport_row,
    )
    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _ctx_of(sweep_inputs)
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    n_u = np.asarray(ctx.n_slot, np.float64)
    n_s = np.asarray(ctx.spliced_slot, np.float64)
    A_g = np.asarray(ctx.eff_gdna_global, np.float64)
    tau = np.asarray(ctx.own.tau_lam, np.float64)
    sites = _item5_slots(ctx)
    served = 0
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if not (is_exon[left[b]] and is_exon[right[b]]):
            continue
        outs = [
            (s, i) for (s, i) in prepared.rule if (s == b and is_exon[i]) or (i == b and is_exon[s])
        ]
        if int(b) not in sites:
            assert not outs, f"boundary {b} carries exon|exon rules without a served terminus"
            continue
        served += 1
        (o,) = (
            {s if s != b else i for (s, i) in outs} & {int(left[b]), int(right[b])}
            if outs
            else (None,)
        )
        assert o is not None
        want = splice_out_row(
            _strand_row_of(ctx, strand, lam, o), lam, n_u[b], n_s[b], A_g[b], A_g[o]
        )
        if tau[o] > 0.0 and prepared.own[o] is not None:
            got = prepared.rule[(int(o), int(b))](prepared.own[o], None)
            np.testing.assert_allclose(got - got.max(), want - want.max(), rtol=0, atol=1e-10)
        if tau[b] > 0.0:
            le = face_map_lambda(lam, n_u[b], A_g[b], A_g[b], A_g[o], A_g[o], n_s[b] / A_g[b])
            want_o = transport_row(_strand_row_of(ctx, strand, lam, b), lam, le, n_u[b], n_s[b])
            got_o = prepared.rule[(int(b), int(o))](prepared.own[b], None)
            np.testing.assert_allclose(
                got_o - got_o.max(), want_o - want_o.max(), rtol=0, atol=1e-10
            )
    assert served >= 2, (
        "the toy must carry two served terminus boundaries or this gate proves nothing"
    )
    cleared = pol.prepare(_terminus_flags_cleared(ctx))
    for b in np.flatnonzero(is_bnd & (left >= 0) & (right >= 0)):
        if is_exon[left[b]] and is_exon[right[b]]:
            assert not any(s == b or i == b for (s, i) in cleared.rule), (
                f"rules survive at {b} with no terminus"
            )


def test_the_level_rule_serves_every_terminus_inside_from_the_measurement_alone(sweep_inputs):
    """THE LEVEL RULE (MESSAGE_PLAN.md step A) as a rule: at every terminus boundary with an inside exon —
    exon|exon and exon|intron alike — the rule applied to the boundary's OWN claim equals the
    independent recompute; what the boundary HOLDS changes nothing (an imputation never becomes a
    level); with no own claim the crossing total's upper bound is what crosses; the rule vanishes when
    the terminus bits are cleared; and no library-wide parameter exists — changing another pair's
    counts leaves this pair's message unchanged."""
    import dataclasses as _dc

    from rigel.calibration.simplex_logodds import _logodds_grid

    pol, _p, n_grid, window = _full_policy(sweep_inputs)
    ctx = _with_populated_inside(_ctx_of(sweep_inputs))
    strand = _strand_of(sweep_inputs)
    lam, _ = _logodds_grid(n_grid, window)
    prepared = pol.prepare(ctx)
    expected, served = _expected_level_rows(sweep_inputs, ctx, strand, lam)
    assert len(served) >= 2, (
        "the toy must carry two served terminus pairs or this gate proves nothing"
    )
    spike = -0.5 * ((lam - 2.0) / 0.1) ** 2  # a held imputation that must not cross a level face
    for b, i, _kind in served:
        fn = prepared.rule.get((b, i))
        assert fn is not None, f"no level rule at terminus pair ({b}, {i})"
        got = fn(prepared.own[b], None)
        np.testing.assert_allclose(
            got - got.max(), expected[i] - expected[i].max(), rtol=0, atol=1e-10
        )
        np.testing.assert_array_equal(
            fn(prepared.own[b], spike), got, err_msg="what is held crossed a level face"
        )
        bound = fn(None, spike)
        assert bound is not None and np.ptp(bound) > 0.0 and np.all(bound <= 0.0), (
            "no upper bound without a claim"
        )
    cleared = pol.prepare(_terminus_flags_cleared(ctx))
    for b, i, _kind in served:
        assert (b, i) not in cleared.rule, f"a level rule survives at ({b}, {i}) with no terminus"
    b0, i0, _k = served[0]
    b1, i1, _k = served[-1]
    cnt = np.asarray(ctx.unspliced_count, np.float64).copy()
    cnt[i1] = cnt[i1] * 3.0 + 7.0
    other = pol.prepare(_dc.replace(ctx, unspliced_count=cnt, n_slot=cnt.sum(axis=1)))
    np.testing.assert_array_equal(
        other.rule[(b0, i0)](other.own[b0], None), prepared.rule[(b0, i0)](prepared.own[b0], None)
    )


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
