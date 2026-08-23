"""STAGES 3 AND 4 OF THE FIRST-PASS REDESIGN — the pass-0 FAN-OUT policy, falsified clause by clause.

**Stage 3**: solved single-stranded introns speak to their flanking ``intron|exon`` boundaries on the
λ channel — the intron's own local solve, clipped into the grid domain, at ``tau_lam``, with NO drift
allowance (measured: `hop_currency.py` excess 0.2–0.4 % on both panels, capture ON and OFF). The
licence lives on the SOURCE alone; the stage-3 destination is implied by adjacency.

**Stage 4**: each pass-through boundary COMPOSES the intron claim with its own strand-λ evidence, and
the licensed exons receive the jointly-derived transfer (:func:`flank_to_exon_lambda`) from the sj
FACE on their side of the flank, routed by the completeness bit — two-sided λ where the flank's
account is complete (two complete flanks fuse, precisions adding), the one-sided at-least-this-much-
RNA bound on the exon's live strand where it is not (two incomplete flanks keep the BINDING bound).
Nothing echoes back into the introns — the no-echo gate is what makes the returned destination mask
real.

The unit tests drive the policy through its own interface on `make_chain_parts` fixtures, emulating
the backbone's source-indexed gather verbatim, with every stage-4 expectation composed from the
SEPARATELY-GATED primitives so these tests gate the WIRING (inputs, faces, masks, routing) alone; the
integration tests run the REAL `solve_chain` and assert exactly the claimed populations move.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.messages import NeighbourState
from rigel.calibration.messages.fanout import FanOutPolicy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_geometry import init_beliefs
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.structural_claims import build_structural_claims
from rigel.calibration.sweep import solve_chain

from _synthetic import make_chain_parts

_L = 10.0
_INTERGENIC = 0


def _logit(f: float) -> float:
    return float(np.log(f) - np.log1p(-f))


def _fixture(signatures, f_g, tau_lam, **counts):
    """parts + claims + a duck ctx carrying what the policy reads — the stage-3 view: no flags, no
    counts at the flanks, so no exon is licensed and the stage-4 machinery stays provably inert."""
    return _rich_fixture(signatures, f_g, tau_lam, **counts)


def _deliver(parts, relay):
    """The backbone's gather, verbatim: publish each direction, index the state at the SOURCE."""
    n = int(parts.chain.n_slots)
    left = np.asarray(parts.chain.left, np.int64)
    right = np.asarray(parts.chain.right, np.int64)
    out = {}
    for backward, nbr in ((False, left), (True, right)):
        step, publish = relay.scan(backward=backward)
        for i in range(n):
            s = int(nbr[i])
            if s >= 0:
                step(s, i)
        src = np.clip(nbr, 0, n - 1)
        out[backward] = NeighbourState(
            state=tuple(np.asarray(a)[src] for a in publish()), valid=nbr >= 0, src=src
        )
    return relay.deliver(out[False], out[True])


#: R0 intergenic · R1 exon+ · R2 intron+ · R3 exon+ · R4 intergenic — slots R B R B R B R B R,
#: so slot 4 is the intron and slots 3/5 are its two ``intron|exon`` boundaries.
SIGS = [_INTERGENIC, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, _INTERGENIC]
#: a local solve everywhere NON-ZERO on purpose: only the LICENCE may silence a slot, so a licence
#: clause that leaks shows up as a delivered claim rather than hiding behind an empty input.
#: Slot order R B R B R B R B R — the intron REGION is slot 4, its boundaries slots 3 and 5.
F_G = [1.0, 0.9, 0.3, 0.45, 0.8, 0.55, 0.4, 0.9, 1.0]
TAU = [9.0, 7.0, 25.0, 30.0, 50.0, 30.0, 25.0, 7.0, 9.0]


def test_the_boundaries_receive_the_introns_claim():
    """Each flanking boundary gets the intron's λ mode at the intron's own tau_lam — from the intron
    side only (the other flank is an exon, which is not a licensed source)."""
    parts, claims, ctx = _fixture(SIGS, F_G, TAU)
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert msg.lam_mode is not None and msg.lam_prec is not None
    for b in (3, 5):
        assert msg.lam_prec[b] == pytest.approx(50.0)
        assert msg.lam_mode[b] == pytest.approx(_logit(0.8))


def test_the_licence_silences_everything_else():
    """Zero delivered precision at every slot the stage does not claim: regions of every kind, the
    gene-edge boundaries — even though every slot's own solve is non-zero."""
    parts, claims, ctx = _fixture(SIGS, F_G, TAU)
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    quiet = np.ones(9, bool)
    quiet[[3, 5]] = False
    assert not np.asarray(msg.lam_prec)[quiet].any()
    assert msg.gdna_mode is None and msg.rna_mode is None and msg.theta_mode is None


def test_two_intron_flanks_fuse_by_precision():
    """A boundary between TWO single-stranded introns receives both claims, fused as a precision-
    weighted mean with the precisions ADDED — two independent measurements of one shared composition."""
    sigs = [_INTERGENIC, BIT_INTRON_POS, BIT_INTRON_POS, _INTERGENIC]  # slots R B R B R B R
    f_g = [1.0, 0.5, 0.9, 0.5, 0.6, 0.5, 1.0]
    tau = [9.0, 9.0, 30.0, 9.0, 10.0, 9.0, 9.0]
    parts, claims, ctx = _fixture(sigs, f_g, tau)
    assert claims.ss_intron_boundary[3], "fixture: slot 3 must be an intron|intron boundary"
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    want = (30.0 * _logit(0.9) + 10.0 * _logit(0.6)) / 40.0
    assert msg.lam_prec[3] == pytest.approx(40.0)
    assert msg.lam_mode[3] == pytest.approx(want)


def test_an_empty_intron_transports_nothing():
    """``tau_lam = 0`` at the source ⇒ no claim exists — a licence on the CLAIM, not a floored knob."""
    tau = list(TAU)
    tau[4] = 0.0
    parts, claims, ctx = _fixture(SIGS, F_G, tau)
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert not np.asarray(msg.lam_prec).any()


def test_the_mode_is_clipped_into_the_grid_domain():
    """A vertex local solve (f_g exactly 1) must arrive at λ = +L, never off-grid
    (`TRAPS: off-grid-message-mode`) — and the clip must follow the window the solve runs on."""
    f_g = list(F_G)
    f_g[4] = 1.0
    parts, claims, ctx = _fixture(SIGS, f_g, TAU)
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert msg.lam_mode[3] == pytest.approx(_L)
    ctx.logodds_window = 5.0
    msg5 = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert msg5.lam_mode[3] == pytest.approx(5.0)


def test_an_ambig_intron_is_not_a_source():
    """Both strands admissible ⇒ outside the stage-0 substrate ⇒ silent, however solved it looks."""
    sigs = [
        _INTERGENIC,
        BIT_EXON_POS | BIT_EXON_NEG,
        BIT_INTRON_POS | 0x4,
        BIT_EXON_POS,
        _INTERGENIC,
    ]
    parts, claims, ctx = _fixture(sigs, F_G, TAU)
    assert not claims.ss_intron_region.any(), "fixture: the AMBIG intron must not be claimed"
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert not np.asarray(msg.lam_prec).any()


def test_through_the_real_backbone_only_the_claimed_boundaries_move():
    """Integration: `solve_chain` under FanOutPolicy vs SilentPolicy on the same fixture — the two
    solves differ AT the intron|exon boundaries and are BIT-IDENTICAL everywhere else."""
    parts = make_chain_parts(
        SIGS,
        region_pos=np.array([0.0, 50.0, 40.0, 50.0, 0.0]),
        region_neg=np.array([0.0, 2.0, 2.0, 2.0, 0.0]),
        boundary_pos=np.array([1.0, 30.0, 30.0, 1.0]),
        boundary_neg=np.array([1.0, 20.0, 20.0, 1.0]),
    )
    claims = build_structural_claims(parts.chain, parts.statics)

    def run(policy):
        belief = init_beliefs(
            parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60
        )
        return solve_chain(
            parts.chain,
            parts.statics,
            parts.geometry,
            belief,
            parts.region_arrays,
            rna_sense_frac=0.95,
            n_rna_obs=10000.0,
            n_gdna_obs=10000.0,
            n_grid=60,
            logodds_window=_L,
            policy=policy,
        )

    fan, silent = run(FanOutPolicy()), run(SilentPolicy())
    moved = np.asarray(fan.f_g) != np.asarray(silent.f_g)
    assert moved[claims.ss_intron_boundary].all(), "the claimed boundaries must hear the message"
    assert not moved[~claims.ss_intron_boundary].any(), "nothing else may move"


# ── stage 4's transfer function: route-merge ∘ reframe, jointly derived ────────────────────────────


def _f2e(lam_b, tau_b, U, S, E_sj, gx, rx, gc, rc, L=_L):
    from rigel.calibration.messages.fanout import flank_to_exon_lambda

    return flank_to_exon_lambda(
        np.atleast_1d(np.asarray(lam_b, np.float64)),
        np.atleast_1d(np.asarray(tau_b, np.float64)),
        np.atleast_1d(np.asarray(U, np.float64)),
        np.atleast_1d(np.asarray(S, np.float64)),
        np.atleast_1d(np.asarray(E_sj, np.float64)),
        np.atleast_1d(np.asarray(gx, np.float64)),
        np.atleast_1d(np.asarray(rx, np.float64)),
        np.atleast_1d(np.asarray(gc, np.float64)),
        np.atleast_1d(np.asarray(rc, np.float64)),
        L,
    )


def test_flank_to_exon_no_flux_is_a_pure_reparameterization():
    """With no spliced flux the transfer is λ_b plus an opportunity constant, at EXACTLY the source's
    precision: the count noise cancels in the ratio and the composition terms collapse to 1/τ_b.
    ⛔ This is the covariance gate — a split (density, variance) interface between route-merge and
    reframe double-counts the shared count and λ (`TRAPS: two-gaussians-one-latent`) and fails here."""
    lam, tau = _f2e(0.7, 25.0, U=40.0, S=0.0, E_sj=100.0, gx=180.0, rx=220.0, gc=150.0, rc=90.0)
    assert lam[0] == pytest.approx(0.7 + np.log(220.0 * 150.0 / (180.0 * 90.0)), abs=1e-12)
    assert tau[0] == pytest.approx(25.0, abs=1e-9)


def test_flank_to_exon_capture_invariance():
    """A uniform capture pull multiplies every COUNT at the locus by one factor; the density RATIO —
    and therefore the delivered λ — must not move at all (the reframe's whole reason to exist)."""
    base_lam, _ = _f2e(
        0.3, 12.0, U=50.0, S=30.0, E_sj=100.0, gx=200.0, rx=200.0, gc=120.0, rc=100.0
    )
    pulled_lam, _ = _f2e(
        0.3, 12.0, U=50.0 * 37.0, S=30.0 * 37.0, E_sj=100.0, gx=200.0, rx=200.0, gc=120.0, rc=100.0
    )
    assert pulled_lam[0] == base_lam[0]


def test_flank_to_exon_variance_matches_monte_carlo():
    """The delta-method precision against brute force: sample the three raw statistics from their own
    sampling models, push each through the same arithmetic, and the empirical variance of λ must match
    ``1/τ`` — the working rule's enumeration check on the covariance handling."""
    rng = np.random.default_rng(7)
    lam0, tau0, U0, S0, E_sj = 0.4, 30.0, 400.0, 250.0, 100.0
    gx, rx, gc, rc = 200.0, 200.0, 150.0, 120.0
    lam, tau = _f2e(lam0, tau0, U0, S0, E_sj, gx, rx, gc, rc)

    n = 40_000
    lam_s = rng.normal(lam0, 1.0 / np.sqrt(tau0), n)
    U_s = rng.poisson(U0, n).astype(float)
    S_s = rng.poisson(S0, n).astype(float)
    f = 1.0 / (1.0 + np.exp(-lam_s))
    rho_g = f * U_s / gx
    rho_r = (1.0 - f) * U_s / rx + S_s / E_sj
    ok = (rho_g > 0) & (rho_r > 0)
    emp = np.log(rho_g[ok] * gc) - np.log(rho_r[ok] * rc)
    assert emp.mean() == pytest.approx(float(lam[0]), abs=0.02)
    assert emp.var() == pytest.approx(1.0 / float(tau[0]), rel=0.10)


def test_flank_to_exon_guards_and_direction():
    """No unspliced measurement, or no composed claim, or a dead opportunity ⇒ NO claim (τ = 0, one
    statement); more spliced flux moves the claim toward RNA; extreme inputs clip into the window."""
    _, tau = _f2e(0.5, 20.0, U=0.0, S=10.0, E_sj=100.0, gx=200.0, rx=200.0, gc=100.0, rc=100.0)
    assert tau[0] == 0.0
    _, tau = _f2e(0.5, 0.0, U=10.0, S=0.0, E_sj=100.0, gx=200.0, rx=200.0, gc=100.0, rc=100.0)
    assert tau[0] == 0.0
    _, tau = _f2e(0.5, 20.0, U=10.0, S=0.0, E_sj=100.0, gx=0.0, rx=200.0, gc=100.0, rc=100.0)
    assert tau[0] == 0.0

    low, _ = _f2e(0.2, 20.0, U=50.0, S=0.0, E_sj=100.0, gx=200.0, rx=200.0, gc=100.0, rc=100.0)
    high, _ = _f2e(0.2, 20.0, U=50.0, S=80.0, E_sj=100.0, gx=200.0, rx=200.0, gc=100.0, rc=100.0)
    assert high[0] < low[0], "spliced flux is RNA — it must push the claim toward RNA"

    lam, _ = _f2e(9.9, 20.0, U=1e12, S=0.0, E_sj=100.0, gx=1e9, rx=1e-6, gc=1e9, rc=1e-9)
    assert -_L <= lam[0] <= _L


# ── stage 4e: the boundary→exon wiring — compose, face selection, sidedness routing, no echo ───────


def _rich_fixture(
    signatures,
    f_g,
    tau_lam,
    *,
    flags=None,
    n_slot=None,
    sj_lo=None,
    sj_hi=None,
    eff_g=None,
    eff_r=None,
    e_sj=100.0,
    **counts,
):
    """The stage-4 duck ctx: everything the depth-2 policy reads, every array hand-settable."""
    parts = make_chain_parts(signatures, boundary_flags=flags, **counts)
    claims = build_structural_claims(parts.chain, parts.statics)
    n = int(parts.chain.n_slots)

    def arr(x, default):
        return np.full(n, default, np.float64) if x is None else np.asarray(x, np.float64)

    def strands(x):
        a = np.zeros((n, 2), np.float64)
        if x is not None:
            a[:, 0] = np.asarray(x, np.float64)
        return a

    geometry = SimpleNamespace(
        sj_count_lo=strands(sj_lo),
        sj_count_hi=strands(sj_hi),
        eff_sj_lo=strands([e_sj] * n),
        eff_sj_hi=strands([e_sj] * n),
    )
    ctx = SimpleNamespace(
        claims=claims,
        own=SimpleNamespace(
            f_g=np.asarray(f_g, np.float64), tau_lam=np.asarray(tau_lam, np.float64)
        ),
        logodds_window=_L,
        n_slot=arr(n_slot, 0.0),
        eff_gdna=arr(eff_g, 100.0),
        eff_rna=arr(eff_r, 100.0),
        free_pos=np.asarray(parts.statics.free_pos, bool),
        free_neg=np.asarray(parts.statics.free_neg, bool),
        geometry=geometry,
    )
    return parts, claims, ctx


def _stage4_setup(extra_flag_at_b1=0):
    """SIGS with a donor bit at slot 3 and an acceptor-role donor bit at slot 5 — both exons licensed; asymmetric
    sj faces so a wrong-face implementation cannot match the expectations."""
    from rigel.calibration.splice_graph import FLAG_DONOR_POS

    flags = np.array(
        [0, FLAG_DONOR_POS | extra_flag_at_b1, FLAG_DONOR_POS, 0], np.uint16
    )  # boundary axis: slots 1, 3, 5, 7
    n_slot = [0.0, 0.0, 0.0, 100.0, 60.0, 90.0, 0.0, 0.0, 0.0]
    sj_lo = [0.0, 0.0, 0.0, 12.0, 0.0, 55.0, 0.0, 0.0, 0.0]
    sj_hi = [0.0, 0.0, 0.0, 99.0, 0.0, 7.0, 0.0, 0.0, 0.0]
    return _rich_fixture(SIGS, F_G, TAU, flags=flags, n_slot=n_slot, sj_lo=sj_lo, sj_hi=sj_hi)


def _expected_exon_claim(ctx, flank, exon, face_count, tau0_intron, lam0_intron):
    """The expectation composed from the SEPARATELY-GATED primitives: fuse at the flank, then the
    joint transfer — so this test gates the WIRING (inputs, faces, masks), nothing else."""
    from rigel.calibration.messages.fanout import flank_to_exon_lambda

    own_l = _logit(float(ctx.own.f_g[flank]))
    own_t = float(ctx.own.tau_lam[flank])
    p = tau0_intron + own_t
    lam_c = (tau0_intron * lam0_intron + own_t * own_l) / p
    return flank_to_exon_lambda(
        np.array([lam_c]),
        np.array([p]),
        np.array([float(ctx.n_slot[flank])]),
        np.array([face_count]),
        np.array([100.0]),
        np.array([float(ctx.eff_gdna[flank])]),
        np.array([float(ctx.eff_rna[flank])]),
        np.array([float(ctx.eff_gdna[exon])]),
        np.array([float(ctx.eff_rna[exon])]),
        _L,
    )


def test_stage4_the_exon_receives_the_composed_transfer_from_the_correct_face():
    """Both exons (complete flanks) receive the two-sided λ claim: the intron's claim fused with the
    flank's OWN strand evidence, then transferred through the joint function — and each exon takes the
    sj face on ITS side of the flank (R1 is LEFT of its flank ⇒ the _lo face, 12; R3 is RIGHT of its
    flank ⇒ the _hi face, 7). A wrong face cannot reproduce these numbers."""
    parts, claims, ctx = _stage4_setup()
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    want_r1 = _expected_exon_claim(
        ctx, flank=3, exon=2, face_count=12.0, tau0_intron=50.0, lam0_intron=_logit(0.8)
    )
    want_r3 = _expected_exon_claim(
        ctx, flank=5, exon=6, face_count=7.0, tau0_intron=50.0, lam0_intron=_logit(0.8)
    )
    assert msg.lam_prec[2] == pytest.approx(float(want_r1[1][0]))
    assert msg.lam_mode[2] == pytest.approx(float(want_r1[0][0]))
    assert msg.lam_prec[6] == pytest.approx(float(want_r3[1][0]))
    assert msg.lam_mode[6] == pytest.approx(float(want_r3[0][0]))


def test_stage4_the_intron_receives_no_echo():
    """The composed boundary states must NOT flow back into the intron: stage 2 solved it, and an
    echo would hand its own strand evidence back at message precision. Zero on every channel."""
    parts, claims, ctx = _stage4_setup()
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert msg.lam_prec[4] == 0.0
    assert msg.rna_prec is None or (msg.rna_prec[0][4] == 0.0 and msg.rna_prec[1][4] == 0.0)


def test_stage4_an_incomplete_flank_rides_the_one_sided_rna_channel():
    """With a transcript-end bit facing R1 (TES at its flank), R1's claim becomes the one-sided
    at-least-this-much-RNA bound: the RNA-share channel on the exon's live strand, flagged one-sided,
    at the Jacobian precision — and the λ channel stays silent there. R3 is untouched."""
    from rigel.calibration.splice_graph import FLAG_TES_POS

    parts, claims, ctx = _stage4_setup(extra_flag_at_b1=FLAG_TES_POS)
    assert not claims.exon_flank_right_complete[2] and claims.exon_flank_right[2]
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    want = _expected_exon_claim(
        ctx, flank=3, exon=2, face_count=12.0, tau0_intron=50.0, lam0_intron=_logit(0.8)
    )
    lam_e, tau_e = float(want[0][0]), float(want[1][0])
    f_e = 1.0 / (1.0 + np.exp(-lam_e))
    assert msg.lam_prec[2] == 0.0
    assert msg.rna_one_sided is not None and bool(msg.rna_one_sided[2])
    assert msg.rna_mode[0][2] == pytest.approx(float(np.log1p(-f_e)))
    assert msg.rna_prec[0][2] == pytest.approx(tau_e / f_e**2)
    assert msg.rna_prec[1][2] == 0.0, "the dead strand says nothing"
    assert msg.lam_prec[6] > 0.0 and not msg.rna_one_sided[6]


def test_stage4_two_complete_flanks_fuse_on_lambda():
    """An exon between two solved introns, both flanks licensed and complete: the two independent
    transfers fuse precision-weighted, precisions adding."""
    from rigel.calibration.splice_graph import FLAG_DONOR_POS

    sigs = [_INTERGENIC, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, _INTERGENIC]
    f_g = [1.0, 0.9, 0.85, 0.5, 0.4, 0.6, 0.75, 0.9, 1.0]
    tau = [9.0, 7.0, 40.0, 20.0, 0.0, 25.0, 35.0, 7.0, 9.0]
    flags = np.array([0, FLAG_DONOR_POS, FLAG_DONOR_POS, 0], np.uint16)
    n_slot = [0.0, 0.0, 0.0, 80.0, 0.0, 70.0, 0.0, 0.0, 0.0]
    parts, claims, ctx = _rich_fixture(sigs, f_g, tau, flags=flags, n_slot=n_slot)
    assert claims.exon_flank_left_complete[4] and claims.exon_flank_right_complete[4]
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    left = _expected_exon_claim(
        ctx, flank=3, exon=4, face_count=0.0, tau0_intron=40.0, lam0_intron=_logit(0.85)
    )
    right = _expected_exon_claim(
        ctx, flank=5, exon=4, face_count=0.0, tau0_intron=35.0, lam0_intron=_logit(0.75)
    )
    lt, rt = float(left[1][0]), float(right[1][0])
    want_prec = lt + rt
    want_mode = (lt * float(left[0][0]) + rt * float(right[0][0])) / want_prec
    assert msg.lam_prec[4] == pytest.approx(want_prec)
    assert msg.lam_mode[4] == pytest.approx(want_mode)


def test_stage4_two_incomplete_flanks_keep_the_binding_bound():
    """Two one-sided bounds do not fuse (a precision mix of bounds is not a bound): the BINDING one —
    the smaller λ, the larger RNA floor — is kept whole, on the RNA channel."""
    from rigel.calibration.splice_graph import FLAG_DONOR_POS, FLAG_TES_POS, FLAG_TSS_POS

    sigs = [_INTERGENIC, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, _INTERGENIC]
    f_g = [1.0, 0.9, 0.85, 0.5, 0.4, 0.6, 0.75, 0.9, 1.0]
    tau = [9.0, 7.0, 40.0, 20.0, 0.0, 25.0, 35.0, 7.0, 9.0]
    flags = np.array(
        [0, FLAG_DONOR_POS | FLAG_TSS_POS, FLAG_DONOR_POS | FLAG_TES_POS, 0], np.uint16
    )
    n_slot = [0.0, 0.0, 0.0, 80.0, 0.0, 70.0, 0.0, 0.0, 0.0]
    parts, claims, ctx = _rich_fixture(sigs, f_g, tau, flags=flags, n_slot=n_slot)
    assert claims.exon_flank_left[4] and not claims.exon_flank_left_complete[4]
    assert claims.exon_flank_right[4] and not claims.exon_flank_right_complete[4]
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    left = _expected_exon_claim(
        ctx, flank=3, exon=4, face_count=0.0, tau0_intron=40.0, lam0_intron=_logit(0.85)
    )
    right = _expected_exon_claim(
        ctx, flank=5, exon=4, face_count=0.0, tau0_intron=35.0, lam0_intron=_logit(0.75)
    )
    (lam_w, tau_w) = min(
        (float(left[0][0]), float(left[1][0])), (float(right[0][0]), float(right[1][0]))
    )
    f_w = 1.0 / (1.0 + np.exp(-lam_w))
    assert msg.lam_prec[4] == 0.0
    assert bool(msg.rna_one_sided[4])
    assert msg.rna_mode[0][4] == pytest.approx(float(np.log1p(-f_w)))
    assert msg.rna_prec[0][4] == pytest.approx(tau_w / f_w**2)


def test_stage4_through_the_real_backbone():
    """Integration with flags: `solve_chain` under FanOutPolicy vs SilentPolicy — the claimed
    boundaries AND the licensed exons move, and nothing else does."""
    from rigel.calibration.splice_graph import FLAG_DONOR_POS

    parts = make_chain_parts(
        SIGS,
        region_pos=np.array([0.0, 50.0, 40.0, 50.0, 0.0]),
        region_neg=np.array([0.0, 2.0, 2.0, 2.0, 0.0]),
        boundary_pos=np.array([1.0, 30.0, 30.0, 1.0]),
        boundary_neg=np.array([1.0, 20.0, 20.0, 1.0]),
        boundary_flags=np.array([0, FLAG_DONOR_POS, FLAG_DONOR_POS, 0], np.uint16),
    )
    claims = build_structural_claims(parts.chain, parts.statics)
    assert claims.solvable_exon.any(), "fixture: the flags must license the exons"

    def run(policy):
        belief = init_beliefs(
            parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60
        )
        return solve_chain(
            parts.chain,
            parts.statics,
            parts.geometry,
            belief,
            parts.region_arrays,
            rna_sense_frac=0.95,
            n_rna_obs=10000.0,
            n_gdna_obs=10000.0,
            n_grid=60,
            logodds_window=_L,
            policy=policy,
        )

    fan, silent = run(FanOutPolicy()), run(SilentPolicy())
    moved = np.asarray(fan.f_g) != np.asarray(silent.f_g)
    allowed = claims.ss_intron_boundary | claims.solvable_exon
    assert moved[claims.ss_intron_boundary].all(), "the boundaries must hear stage 3"
    assert moved[claims.solvable_exon].all(), "the licensed exons must hear stage 4"
    assert not moved[~allowed].any(), "nothing else may move — the intron must hear NO echo"


def test_stage4_a_minus_strand_exon_rides_the_neg_channel():
    """The locus-F lesson replayed at stage 4: on a − exon the one-sided bound must land on the
    NEGATIVE strand's RNA channel, with the positive channel silent — a fixture family that is all
    ``+`` cannot see a hard-coded column."""
    from rigel.calibration.signature import BIT_EXON_NEG, BIT_INTRON_NEG
    from rigel.calibration.splice_graph import FLAG_DONOR_NEG, FLAG_TSS_NEG

    sigs = [_INTERGENIC, BIT_EXON_NEG, BIT_INTRON_NEG, BIT_EXON_NEG, _INTERGENIC]
    flags = np.array([0, FLAG_DONOR_NEG | FLAG_TSS_NEG, FLAG_DONOR_NEG, 0], np.uint16)
    n_slot = [0.0, 0.0, 0.0, 100.0, 60.0, 90.0, 0.0, 0.0, 0.0]
    parts, claims, ctx = _rich_fixture(sigs, F_G, TAU, flags=flags, n_slot=n_slot)
    assert claims.exon_flank_right[2] and not claims.exon_flank_right_complete[2], (
        "fixture: R1's right flank must be licensed and incomplete (TSS_NEG is its HIGH-end bit)"
    )
    # the sj face arrays: the − strand column must carry the flux, so a hard-coded column 0 reads 0
    ctx.geometry.sj_count_lo[3, 1] = 12.0
    ctx.geometry.eff_sj_lo[3, 1] = 100.0
    msg = _deliver(parts, FanOutPolicy().prepare(ctx))
    assert msg.rna_one_sided is not None and bool(msg.rna_one_sided[2])
    assert msg.rna_prec[1][2] > 0.0, "the − strand carries the claim"
    assert msg.rna_prec[0][2] == 0.0, "the + strand says nothing"
    want = _expected_exon_claim(
        ctx, flank=3, exon=2, face_count=12.0, tau0_intron=50.0, lam0_intron=_logit(0.8)
    )
    f_e = 1.0 / (1.0 + np.exp(-float(want[0][0])))
    assert msg.rna_mode[1][2] == pytest.approx(float(np.log1p(-f_e)))
