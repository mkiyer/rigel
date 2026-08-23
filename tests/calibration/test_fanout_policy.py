"""STAGE 3 OF THE FIRST-PASS REDESIGN — the pass-0 FAN-OUT policy, falsified clause by clause.

`rigel.calibration.messages.fanout.FanOutPolicy`: solved single-stranded introns speak to their
flanking ``intron|exon`` boundaries on the λ channel, and nothing else says anything. The message is
the intron's own local solve — mode ``logit(f_g)`` clipped into the grid domain, precision
``tau_lam`` (the strand + intron-factory Fisher; the reference's location term does not enter, the
§9c.2 ruling) — with NO drift allowance, which is measured rather than assumed: `hop_currency.py`'s
excess-over-floor for this hop is 0.2–0.4 % on both panels, capture ON and OFF (2026-08-23).

The licence is STRUCTURAL and lives on the SOURCE alone (`structural_claims.ss_intron_region`); the
destination end is IMPLIED by adjacency and an empty intron transports nothing by arithmetic — the
policy docstring carries both derivations, and the perturbation sweep is what proved the redundant
masks did no work.

The unit tests drive the policy through its own interface on a `make_chain_parts` fixture, emulating
the backbone's source-indexed gather verbatim; the integration test runs the REAL `solve_chain` and
asserts the message moves ONLY the claimed boundaries.
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
    """parts + claims + a duck ctx carrying exactly what the policy reads. ``f_g``/``tau_lam`` are the
    per-SLOT local solve the policy transports (chain order: R B R B … R)."""
    parts = make_chain_parts(signatures, **counts)
    claims = build_structural_claims(parts.chain, parts.statics)
    ctx = SimpleNamespace(
        claims=claims,
        own=SimpleNamespace(
            f_g=np.asarray(f_g, np.float64), tau_lam=np.asarray(tau_lam, np.float64)
        ),
        logodds_window=_L,
    )
    return parts, claims, ctx


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
