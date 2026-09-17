"""Two properties of ``simplex_logodds``: how ψ prices the simplex vertex, and that the composition
it publishes is a point on the simplex by construction.

ψ's reference is one constant used twice with opposite signs, which on the ``λ = logit f_g`` axis is
Beta(½,½) with a tail of ``−C·|λ|``. A Gaussian claim of precision ``p`` written on ``log f_c`` has
an exponentially decaying gradient there, so it balances the reference at ``λ* = ½·log(p/C)``: the
vertex is priced rather than forbidden, and no achievable depth buys it. G1 and G1b measure that
law at the two vertices, G2 the tail slope, G3 the halves' orthogonality, G4 how weak a lever
softening the reference is, G5 the control that the same claim on ``λ`` is honoured at par, and G6
ψ's blindness to the certified-splice channel. The second half gates the read-out, where
``_compose`` maps ψ's level and tilt onto the simplex so that closure is an identity.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from _psi_reference import jeffreys_arms, lse, row_at, strand_loglik_mixture
from scipy.special import expit, log_expit

from rigel.calibration import simplex_logodds as SL
from rigel.calibration.simplex_logodds import (
    _logodds_grid,
    _solve_regions_logodds_all,
    compose,
    posterior_median_fg,
    psi_cube,
)

#: κ = ½ EXACTLY, no overdispersion, no fitted priors. On that substrate the strand term is bit-flat, so
#: the only things speaking about λ are the message under test and ψ's reference — which isolates the
#: reference without ablating anything.
#: the lattice is deliberately FINE: the price law is read in λ, and a coarse lattice quantises the
#: very quantity being measured. The window ``L`` is the one real hard limit and G1 stays clear of it.
_BASE = dict(kappa=0.5, od_g=0.0, od_r=0.0, n_grid=4096, L=10.0)

#: one decade of precision on a ``log f_c`` message buys this many nats of log-odds, by derivation
#: (``λ* = ½·log(p/C)`` ⇒ ``dλ/dlog₁₀p = ½·ln 10``). Not a tuned tolerance — the prediction itself.
_NATS_PER_DECADE = 0.5 * float(np.log(10.0))


def _regions(n: int = 2):
    """``n`` single-strand (+) regions carrying data. Single-strand keeps this on the exact 1-D λ solve."""
    u_pos = np.full(n, 200.0)
    u_neg = np.zeros(n)
    return u_pos, u_neg, np.ones(n, bool), np.zeros(n, bool), u_pos + u_neg, np.zeros(n)


def _solve(**over):
    n = int(over.pop("n", 2))
    u_pos, u_neg, ap, an, mu, ms = _regions(n)
    ms = np.asarray(over.pop("mass_spliced", ms), np.float64)
    kw = dict(_BASE)
    kw.update(over)
    return np.asarray(
        _solve_regions_logodds_all(u_pos, u_neg, ap, an, mu, ms, **kw).gdna_frac, np.float64
    )


def _lam(f):
    f = np.clip(np.asarray(f, np.float64), 1e-300, 1.0 - 1e-16)
    return np.log(f / (1.0 - f))


def _rows(fn, n: int = 2):
    """``n`` copies of one claim as a λ-row on the coarse solve grid — the socket every message uses."""
    lam, _ = _logodds_grid(int(_BASE["n_grid"]), float(_BASE["L"]))
    row = np.asarray(fn(lam), np.float64)
    return dict(lam_logprior=np.tile(row[None, :], (n, 1)))


def _msg_up(p):
    """A claim written on ``log f_g`` that ``f_g = 1`` — mode ``log f_g = 0`` — at precision ``p``."""
    return _rows(lambda lam: -0.5 * float(p) * log_expit(lam) ** 2)


def _msg_dn(p):
    """A claim written on ``log(1 − f_g)`` that ``f_g = 0`` — mode ``0`` — at precision ``p``."""
    return _rows(lambda lam: -0.5 * float(p) * log_expit(-lam) ** 2)


def _msg_lam(lam_t, p):
    """The same kind of claim written on ``λ`` itself: a Gaussian at ``lam_t`` with precision ``p``."""
    return _rows(lambda lam: -0.5 * float(p) * (lam - float(lam_t)) ** 2)


# ── G1 — the price law ──────────────────────────────────────────────────────────────────────────────


def test_G1_a_log_f_g_message_buys_log_odds_only_LOGARITHMICALLY_in_its_precision():
    """The price law. Sweep the message precision over four decades and regress ``λ`` on
    ``log₁₀ p``: the slope must be ``½·ln 10`` — the vertex costs ``e²`` in precision per nat.

    The consequence is the finding: an object with the precision the toy actually measures
    (``p ≈ 10``) sits at ``λ ≈ 1.5``, i.e. ``f_g ≈ 0.82``, and no achievable depth closes that, because
    ``p`` is at best linear in the count while ``λ`` is logarithmic in ``p``.

    Kept clear of ``p ≥ 1e8``, where ``λ`` hits the ``L = 10`` window and the law stops being the
    binding constraint — that ceiling is a separate fact and G5 is what shows it is not the cause."""
    lam_of = {e: float(_lam(_solve(**_msg_up(10.0**e))[0])) for e in (2, 3, 4, 5, 6)}
    es = sorted(lam_of)
    slope = float(np.polyfit(es, [lam_of[e] for e in es], 1)[0])
    assert slope == pytest.approx(_NATS_PER_DECADE, rel=0.15), (slope, _NATS_PER_DECADE, lam_of)
    # and the practical statement: at the precision real objects carry, the answer is far from 1.
    assert float(_solve(**_msg_up(10.0))[0]) < 0.95, _solve(**_msg_up(10.0))[0]


def test_G1b_the_price_law_is_the_SAME_at_the_other_vertex():
    """The mirror. ``f_g → 0`` is bounded by the gDNA half of the same constant, so the slope in
    ``|λ|`` must match G1's — one constant, two signs (which G3 then pins as orthogonal)."""
    lam_of = {e: float(_lam(_solve(**_msg_dn(10.0**e))[0])) for e in (2, 3, 4, 5, 6)}
    es = sorted(lam_of)
    slope = float(np.polyfit(es, [-lam_of[e] for e in es], 1)[0])
    assert slope == pytest.approx(_NATS_PER_DECADE, rel=0.15), (slope, _NATS_PER_DECADE, lam_of)


# ── G2 — the reference is the only slope in the tail ─────────────────────────────────────────────────


def test_G2_psi_slope_in_the_vertex_tail_is_exactly_minus_the_reference_exponent():
    """ψ built directly on the λ grid with a claim written on ``log f_g``: every non-reference term is flat in
    the far tail (κ=½ ⇒ bit-flat strand; the message's gradient decays as ``e^{−2λ}``), so ``dψ/dλ``
    must equal ``−_JEFFREYS_REF``. The whole mechanism in one number.

    ψ is CALLED (`psi_cube`, the solver's own cube), not reimplemented, so this cannot drift from
    what the solver computes (TRAPS: self-checking-validator)."""
    lam, _fg = _logodds_grid(1024, 10.0)
    u_pos, u_neg, ap, an, _mu, _ms = _regions(1)
    psi, _fp, _fn, _tau = psi_cube(
        u_pos,
        u_neg,
        ap,
        an,
        np.full(1, 0.5),
        np.full(1, 0.5),
        np.zeros(1),
        kappa=0.5,
        od_g=0.0,
        od_r=0.0,
        lam=lam,
        ambig=False,
        lam_logprior=(-0.5 * 1e3 * log_expit(lam) ** 2)[None, :],
    )
    p = psi[0, :, 0]
    lo, hi = int(0.90 * lam.size), lam.size - 1
    slope = (p[hi] - p[lo]) / (lam[hi] - lam[lo])
    assert abs(slope + SL._JEFFREYS_REF) < 5e-3, (slope, -SL._JEFFREYS_REF)


# ── G3 — the two halves are orthogonal ──────────────────────────────────────────────────────────────


def _with_exponents(msg: dict, c_g: float, c_r: float) -> dict:
    """The message with ψ's two arms re-written to independent exponents. The arms are ADDITIVE, so an
    exponent change is a λ-row — ``(c − ½)·log f`` per arm — delivered on the same channel as the message
    and reaching the solver's cube exactly as the arm itself does; no patching of the solver, which is
    native. The assertions below demand a MOVE rather than merely a difference, so an ablation that never
    reached ψ would fail them."""
    lam, _ = _logodds_grid(int(_BASE["n_grid"]), float(_BASE["L"]))
    shift = jeffreys_arms(lam, c_g - 0.5, c_r - 0.5)
    return dict(lam_logprior=np.asarray(msg["lam_logprior"], np.float64) + shift[None, :])


def test_G3_each_half_of_the_constant_holds_ONE_vertex_and_is_NEGLIGIBLE_at_the_other():
    """The orthogonality — what makes one fix cover both vertices.

    ``C·log(1−f_g)`` is the only term bounding ``f_g → 1``. Deleting it must move that answer by a full
    nat of log-odds while barely touching ``f_g → 0``, and symmetrically for the other half. Two claims
    per half, so neither can pass by accident.

    Scored in λ, not in ``f_g``: near a vertex ``f_g`` compresses a nat into ~1e-3, and a threshold on
    ``f_g`` would read a large move as a small one.

    Stated as a RATIO, not as bit-identity, and that distinction is a measurement rather than
    caution. On a 256-point lattice (``Δλ = 0.078``) the off-vertex half IS
    bit-identical — but at the 16× finer lattice used here it moves 0.015 nats, so the identity was
    lattice quantisation and not orthogonality (TRAPS: byte-identity-gate: a bit-identity gate has lied in both
    directions). The real claim is that each half is worth two orders of magnitude more at its own
    vertex than at the other, and that survives any resolution."""
    up, dn = _msg_up(1e3), _msg_dn(1e3)
    base_up = float(_lam(_solve(**up)[0]))
    base_dn = float(_lam(_solve(**dn)[0]))

    # delete the f_g→1 bound
    own = float(_lam(_solve(**_with_exponents(up, 0.5, 0.0))[0])) - base_up
    other = abs(float(_lam(_solve(**_with_exponents(dn, 0.5, 0.0))[0])) - base_dn)
    assert own > 1.0, own  # it MOVED, toward its own vertex
    assert other < 0.1, other  # and barely registered at the other
    assert own > 20.0 * max(other, 1e-12), (own, other)

    # delete the f_g→0 bound instead
    own = base_dn - float(_lam(_solve(**_with_exponents(dn, 0.0, 0.5))[0]))
    other = abs(float(_lam(_solve(**_with_exponents(up, 0.0, 0.5))[0])) - base_up)
    assert own > 1.0, own
    assert other < 0.1, other
    assert own > 20.0 * max(other, 1e-12), (own, other)


# ── G4 — the perturbation that makes G2 non-vacuous ─────────────────────────────────────────────────


def test_G4_softening_the_reference_moves_lambda_MONOTONICALLY_and_by_a_BOUNDED_step():
    """The perturbation. ``λ* = ½·log(p/C)`` predicts that halving ``C`` buys only ``½·ln 2 = 0.35``
    nats — so softening the reference is a *weak* lever, and that is why "just lower the exponent" is
    not the fix (and why ``C → 0``, which would be, makes ψ improper — TRAPS: no-prior-means-haldane).

    Gated as monotone-and-bounded rather than as an exact step: the message's residual curvature and the
    other half's own term both contribute an offset that does not cancel across ``C``."""
    msg = _msg_up(1e4)
    seen = []
    for c in (0.5, 0.25, 0.125, 0.0625):
        seen.append(float(_lam(_solve(**_with_exponents(msg, 0.5, c))[0])))
    steps = np.diff(seen)
    assert np.all(steps > 0.0), seen  # monotone toward the vertex
    assert np.all(steps < 1.0), seen  # and each halving buys well under a nat
    # four halvings — a 16x softening — buy less than the ONE decade of precision G1 measured.
    assert (seen[-1] - seen[0]) < 4.0 * _NATS_PER_DECADE, seen


# ── G5 — the CONTROL: on λ, the same claim is honoured at par ────────────────────────────────────────


def test_G5_the_SAME_claim_delivered_on_lambda_is_honoured_AT_PAR():
    """The control that makes G1 a falsification (TRAPS: a-gate-that-already-passed).

    ``f_g = 1 − 1e-4`` is ``λ = +9.21``, which the ``L = 10`` grid represents exactly. Delivered on the λ
    axis the claim is honoured; delivered on ``log f_g`` — the SAME claim, the same precision — it is
    not. So the vertex is reachable and the lattice is not the limit: the coordinate is the price.

    This is what forbids "the grid cannot represent a vertex" as a diagnosis."""
    target = 1.0 - 1e-4
    lam_t = float(_lam(np.array([target]))[0])
    on_lambda = float(_solve(**_msg_lam(lam_t, 1e3))[0])
    on_log_fg = float(_solve(**_msg_up(1e3))[0])
    assert on_lambda > 0.999, on_lambda
    assert on_log_fg < 0.995, on_log_fg
    assert float(_lam(np.array([on_lambda]))[0] - _lam(np.array([on_log_fg]))[0]) > 2.0


# ── G6 — the root cause: the certified-RNA channel is not read ───────────────────────────────────────


def test_G6_psi_is_BLIND_to_the_certified_RNA_channel():
    """The root cause, as an observable — and the one a fix must move
    (TRAPS: name-the-observable-per-site).

    ``f_g = 1`` is a claim about ``ρ_r = 0``, and no gDNA-side message can establish it at any
    precision: G1 prices exactly that attempt. The channel that speaks about RNA directly is the
    certified splice — ``boundary_spliced``, the fragments that crossed this boundary contiguously
    having spliced elsewhere — which has no gDNA term at all, because a spliced fragment cannot be
    gDNA. (Spliced-vs-unspliced is a property of the fragment, never a second species of RNA.)

    ψ does not read it: sweeping ``mass_spliced`` from 0 to 1e5 at fixed strand counts leaves
    ``f_g`` bit-identical. The exclusion is deliberate, and its reason is sound as far as it goes —
    the spliced mass must not be double-counted into the unspliced likelihood — but "do not add it
    to that likelihood" is not "discard the information in it", and this gate is what fires if the
    second is ever addressed. ``boundary_spliced`` is boundary-only and structurally zero on every
    toy geometry, so the sweep uses ``mass_spliced`` purely as a handle; it is not a claim about
    the channel's size.

    Read ``test_certified_rna_licence.py`` before acting on this: the obvious fix is refuted. A
    zero spliced count is not evidence for the vertex, because the chance a crossing RNA fragment
    shows a visible splice is a free nuisance and ``S = 0`` is explained at ``f_g = 1`` too; the
    information lives entirely in ``S > 0`` and points away from this vertex. What survives is the
    algebra — reference plus term is exactly Beta(½, ½+S) — and not the coefficient, which scores
    worse than the uninformative reference on a third of the ladder."""
    quiet = float(_solve(mass_spliced=np.zeros(2), **_msg_up(10.0))[0])
    loud = float(_solve(mass_spliced=np.full(2, 1e5), **_msg_up(10.0))[0])
    assert loud == pytest.approx(quiet, abs=0.0), (quiet, loud)
    # and the object it would matter for is nowhere near the vertex without it.
    assert quiet < 0.95, quiet


# ── the read-out: ψ's composition is a point on the simplex by construction ────────────────────


_L = 10.0


def _solve_composition(u_pos, u_neg, *, kappa, allow_pos, allow_neg, n_grid=60):
    d = _solve_regions_logodds_all(
        np.asarray(u_pos, np.float64),
        np.asarray(u_neg, np.float64),
        np.asarray(allow_pos, bool),
        np.asarray(allow_neg, bool),
        np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64),
        np.zeros_like(np.asarray(u_pos, np.float64)),
        kappa=kappa,
        od_g=0.0,
        od_r=0.0,
        n_grid=n_grid,
        L=_L,
    )
    return (
        np.asarray(d.gdna_frac, np.float64),
        np.asarray(d.rna_pos_frac, np.float64),
        np.asarray(d.rna_neg_frac, np.float64),
    )


_DEPTHS = (1.0, 3.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
_KAPPAS = (0.5, 0.75, 0.9, 0.99)
_SPLITS = (0.0, 0.01, 0.2, 0.5)


@pytest.mark.parametrize("kappa", _KAPPAS)
@pytest.mark.parametrize(
    "allow_pos,allow_neg", ((True, False), (False, True), (True, True)), ids=("ss+", "ss-", "ambig")
)
def test_the_composition_closes_exactly_on_every_solved_slot(kappa, allow_pos, allow_neg):
    """The property the parametrisation buys. Not "closes to a tolerance" — closes to float64
    round-off, on both ψ paths, at every depth and strand split, because the read-out is a map onto the
    simplex rather than three independent summaries.

    The perturbation is below: a composition assembled from three independent read-outs on the same
    posterior does not close, so this gate is not passing because the quantity is trivial.
    """
    u_pos = np.array([n * (1.0 - f) for n, f in itertools.product(_DEPTHS, _SPLITS)])
    u_neg = np.array([n * f for n, f in itertools.product(_DEPTHS, _SPLITS)])
    ap = np.full(u_pos.shape, allow_pos)
    an = np.full(u_pos.shape, allow_neg)
    f_g, f_pos, f_neg = _solve_composition(u_pos, u_neg, kappa=kappa, allow_pos=ap, allow_neg=an)
    total = f_g + f_pos + f_neg
    assert np.max(np.abs(total - 1.0)) < 1e-12, (kappa, np.max(np.abs(total - 1.0)))
    # and the dead strand really is dead, so "closes" is not being bought by leaking RNA onto it
    if not allow_pos:
        assert np.all(f_pos == 0.0)
    if not allow_neg:
        assert np.all(f_neg == 0.0)
    # non-vacuous: BOTH components are actually being placed, so closure is not holding because one
    #   side is identically zero. Deliberately not a span check on ``f_g`` — at κ = ½ the strand
    #   channel is dead by derivation and every slot sits near the reference, which is correct.
    assert f_g.max() > 0.0 and (f_pos + f_neg).max() > 0.0, (f_g.max(), (f_pos + f_neg).max())


def test_the_old_read_out_did_not_close_and_the_gap_was_the_skew():
    """The perturbation for the gate above. Reconstruct a composition assembled from three
    independent read-outs — ``f_g`` from the median, the RNA total from the posterior mean — and
    show it misses the simplex by exactly ``median − mean``.

    The posterior mean is recoverable from the shipped output without re-solving: the RNA total is
    ``1 − f_g``, so the three-read-out RNA total is ``1 − E[f_g]`` in its place."""
    lam, fg = _logodds_grid(60, _L)
    u_pos, u_neg = np.array([27.0, 9.0, 3.0]), np.array([3.0, 1.0, 0.0])
    psi = psi_cube(
        u_pos,
        u_neg,
        np.array([True, True, True]),
        np.array([False, False, False]),
        np.full(3, 0.5),
        np.full(3, 0.5),
        np.zeros(3),
        kappa=0.99,
        od_g=0.0,
        od_r=0.0,
        lam=lam,
        ambig=False,
    )[0][:, :, 0]
    post = np.exp(psi - lse(psi, axis=1, keepdims=True))
    median = posterior_median_fg(post, lam)
    mean = np.sum(post * fg[None, :], axis=1)
    old_sum = median + (1.0 - mean)  # a three-read-out composition
    assert np.max(np.abs(old_sum - (1.0 + median - mean))) < 1e-12
    # it genuinely misses the simplex on this fixture — otherwise the gate above is vacuous
    assert np.max(np.abs(old_sum - 1.0)) > 0.02, old_sum
    # and the shipped read-out on the same counts closes
    f_g, f_pos, f_neg = _solve_composition(
        u_pos, u_neg, kappa=0.99, allow_pos=np.full(3, True), allow_neg=np.full(3, False)
    )
    assert np.max(np.abs(f_g + f_pos + f_neg - 1.0)) < 1e-12


def test_compose_enforces_admissibility_so_rna_cannot_leak_onto_a_forbidden_strand():
    """The one way a map onto the simplex can still be wrong. If the
    tilt share were applied blind, a slot with only the + strand admissible and ``w_pos = ½`` would place
    half its RNA on the forbidden − strand, where it is zeroed — and the RNA would simply VANISH, giving
    ``SUM = f_g + (1−f_g)/2``. `compose` therefore restricts the share to the admissible strands itself
    rather than trusting the caller, so the single-strand path has no tilt to supply at all.

    Closure must hold for EVERY ``w_pos``, including ones that are wrong for the slot."""
    f_g = np.array([0.0, 0.05, 0.5, 0.95, 1.0])
    T, F = np.full(5, True), np.full(5, False)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        for ap, an in ((T, T), (T, F), (F, T)):
            p, n = compose(f_g, np.full(5, w), ap, an)
            assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, (w, ap[0], an[0])
            assert np.all(p >= 0.0) and np.all(n >= 0.0)
            assert np.all(p[~ap] == 0.0) and np.all(n[~an] == 0.0)
    # a single-strand slot takes the WHOLE RNA total on its admissible strand, whatever w_pos says
    for w in (0.0, 0.5, 1.0):
        p, n = compose(f_g, np.full(5, w), T, F)
        assert np.allclose(p, 1.0 - f_g) and np.all(n == 0.0), w
        p, n = compose(f_g, np.full(5, w), F, T)
        assert np.allclose(n, 1.0 - f_g) and np.all(p == 0.0), w
    # neither strand admissible ⇒ no RNA to place; the composition is f_g alone, which is a
    #   statement about opportunity and not a closure failure
    p, n = compose(f_g, np.full(5, 0.5), F, F)
    assert np.all(p == 0.0) and np.all(n == 0.0)


def test_a_share_outside_the_unit_interval_cannot_produce_a_negative_fraction():
    """The one way a closing composition can still be nonsense. An unclamped share yields a negative
    fraction that still sums to 1 and would pass every closure gate in this file: at ``w_pos = −1``
    it gives ``f_pos = −0.8, f_neg = 1.6``, summing to 1.

    Not reachable from the shipped AMBIG caller — ``w_pos = m_pos/(m_pos+m_neg)`` is a ratio of two
    non-negative expectations — which is precisely why the constraint is asserted here rather than assumed
    at the call site."""
    f_g = np.array([0.0, 0.3, 0.9])
    T = np.full(3, True)
    for w in (-1.0, -1e-9, 1.0 + 1e-9, 2.0, 1e9):
        p, n = compose(f_g, np.full(3, w), T, T)
        assert np.all(p >= 0.0) and np.all(n >= 0.0), (w, p, n)
        assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, w
    # the perturbation: an IN-range share is passed through untouched, so the clamp is not flattening
    #   the tilt into a constant
    p_lo, _ = compose(f_g, np.full(3, 0.25), T, T)
    p_hi, _ = compose(f_g, np.full(3, 0.75), T, T)
    assert np.all(p_hi > p_lo)


def test_the_ambig_tilt_share_stays_a_share():
    """The AMBIG caller derives ``w_pos`` from two cube expectations of non-negative grid quantities, so
    it is in [0,1] by construction. This asserts that on the SHIPPED grid rather than trusting the
    algebra — an f32 cube reduction is where such a guarantee would quietly fail."""
    for kappa in (0.5, 0.9, 0.99):
        for n in (1.0, 10.0, 1e3, 1e6):
            for frac in (0.0, 0.1, 0.5, 0.9, 1.0):
                f_g, f_pos, f_neg = _solve_composition(
                    np.array([n * (1 - frac)]),
                    np.array([n * frac]),
                    kappa=kappa,
                    allow_pos=np.array([True]),
                    allow_neg=np.array([True]),
                )
                rna = f_pos + f_neg
                assert np.all(rna >= -1e-15) and np.all(f_pos >= -1e-15) and np.all(f_neg >= -1e-15)
                assert np.allclose(rna, 1.0 - f_g, atol=1e-15), (kappa, n, frac)


def test_a_zero_count_slot_publishes_no_data_rather_than_a_composition():
    """The one population that does not close, and it is deliberate. A slot with no fragments, or
    with neither strand admissible, is not dispatched to either ψ solve; it publishes ``(0, 0, 0)``,
    which is "no data" and not a composition claim. `region_init` replaces it with the signature-binary
    init, which IS a simplex point.

    Recorded so the closure gate above is read with its true scope: closure is asserted over SOLVED
    slots, and this is what the complement looks like."""
    f_g, f_pos, f_neg = _solve_composition(
        np.array([0.0, 30.0]),
        np.array([0.0, 1.0]),
        kappa=0.99,
        allow_pos=np.array([True, False]),
        allow_neg=np.array([False, False]),
    )
    assert np.all(f_g == 0.0) and np.all(f_pos == 0.0) and np.all(f_neg == 0.0)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE θ QUADRATURE — the nodes follow the strand term's peak, so the marginal is exact at every depth
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
#
# At fixed λ the strand term is an exact Gaussian in τ (`EQUATIONS.md`), of width ∝ n^{−½}: at 50k
# fragments its θ peak is 0.005 rad against a 60-node lattice's 0.053 step, so a fixed lattice's sum is a
# comb across λ (the ladder's recorded K_t 30 failure was ONE 25k-fragment slot whose peak fell between
# nodes). The rule places K_t nodes across the window where the term is within T nats of its maximum on
# the domain and weights them as the trapezoid rule; these gates hold it to an adaptive reference.

_THETA_KAPPA = 0.99


def _theta_case(n: float, fg_true: float, tau_true: float, K: int = 41):
    """One AMBIG slot with the EXPECTED counts of composition ``(f_g, τ)`` and the variance frozen there."""
    fpos = (1 - fg_true) * (1 + tau_true) / 2
    fneg = (1 - fg_true) * (1 - tau_true) / 2
    p = 0.5 * fg_true + _THETA_KAPPA * fpos + (1 - _THETA_KAPPA) * fneg
    lam, _fg = _logodds_grid(K, 10.0)
    args = (
        np.array([n * p]),
        np.array([n * (1 - p)]),
        np.array([True]),
        np.array([True]),
        np.array([fg_true]),
        np.array([fpos]),
        np.array([fneg]),
    )
    kw = dict(kappa=_THETA_KAPPA, od_g=0.0, od_r=0.0, lam=lam)
    return args, kw


def _theta_reference(args, kw):
    """log M(λ) by adaptive quadrature in θ per λ, the peak located analytically — a second
    implementation of the integral, not of the integrand (the readable strand term of `_psi_reference`)."""
    from scipy.integrate import quad

    u_pos, u_neg, _ap, _an, fgr, fpr, fnr = args
    lam = kw["lam"]
    fg = expit(lam)
    n = float(u_pos[0] + u_neg[0])
    d = float(u_pos[0]) / n - 0.5
    out = np.empty(lam.shape[0])
    for j in range(lam.shape[0]):

        def g(theta, j=j):
            tau = np.sin(theta)
            f_act = 1.0 - fg[j]
            return strand_loglik_mixture(
                float(u_pos[0]),
                n,
                fg[j],
                f_act * (1 + tau) / 2,
                f_act * (1 - tau) / 2,
                _THETA_KAPPA,
                0.0,
                0.0,
                float(fgr[0]),
                float(fpr[0]),
                float(fnr[0]),
            )

        a = (1 - fg[j]) * (_THETA_KAPPA - 0.5)
        th = float(np.arcsin(np.clip(d / a, -1, 1)))
        peak = g(th)
        val, _err = quad(
            lambda t: np.exp(g(t) - peak),
            -np.pi / 2,
            np.pi / 2,
            points=[th],
            limit=400,
            epsabs=0,
            epsrel=1e-11,
        )
        # the tilt atom's measure (EQUATIONS §9f): the continuum normalised to the domain plus the two
        # pure hypotheses at τ = ±1, each at the same reference weight
        atoms = np.logaddexp(g(0.5 * np.pi), g(-0.5 * np.pi))
        out[j] = np.logaddexp(peak + np.log(val / np.pi), atoms)
    return out + jeffreys_arms(lam)


@pytest.mark.parametrize(
    "n, fg_true, tau_true",
    [(500.0, 0.0, 0.0), (50_000.0, 0.0, 0.5), (500_000.0, 0.3, 0.9), (500_000.0, 0.0, 1.0)],
)
def test_the_theta_marginal_matches_adaptive_quadrature_at_every_depth(n, fg_true, tau_true):
    """ψ's θ-marginal, log Σ_k exp ψ(λ, θ_k) over the continuum's nodes and the two atoms, against the
    adaptive reference under the same measure: the λ-SHAPE of the error (a constant offset is harmless;
    a λ-dependent one is a gDNA bias) is below 1e−5 nats at 500 fragments and at 500k, interior tilt,
    near-pure tilt and the strand-pure boundary alike. The fixed 60-node lattice read 0.01–0.5 nats at
    500 fragments and 90–130 at 500k."""
    args, kw = _theta_case(n, fg_true, tau_true)
    psi, _fp, _fn, _tau = psi_cube(*args, ambig=True, **kw)
    delta = lse(psi, axis=2)[0] - _theta_reference(args, kw)
    shape = float(np.max(np.abs(delta - delta.mean())))
    assert shape < 1e-5, (n, fg_true, tau_true, shape)


def _solve_case(args, cube_rows=None, n_grid: int = 101, n_tilt=None):
    """One `_theta_case` slot through the dispatcher, the reference composition its own."""
    u_pos, u_neg, ap, an, fgr, fpr, fnr = args
    return _solve_regions_logodds_all(
        u_pos,
        u_neg,
        ap,
        an,
        u_pos + u_neg,
        np.zeros(1),
        kappa=_THETA_KAPPA,
        od_g=0.0,
        od_r=0.0,
        n_grid=n_grid,
        L=10.0,
        fg_ref=fgr,
        fpos_ref=fpr,
        fneg_ref=fnr,
        cube_rows=cube_rows,
        n_tilt=n_tilt,
    )


def test_the_derived_node_count_is_converged():
    """K_t = 2T/π + 1 = 24 is the node count at which the trapezoid error on the window falls below
    e^{−T}: the read-out at 24 nodes equals the read-out at 60 to 1e−6 on a 500k-fragment slot at an
    interior tilt, where a fixed lattice at 24 and 60 disagreed by tenths. There is no other θ count
    anywhere: the config has no tilt knob and ψ's module has no tilt lattice; the dispatcher's ``n_tilt``
    exists for this gate."""
    from rigel.config import CalibrationConfig

    assert SL._TILT_NODES == 24
    assert not hasattr(CalibrationConfig(), "sweep_n_tilt")
    assert not hasattr(SL, "_tilt_grid")
    args, _kw = _theta_case(500_000.0, 0.1, 0.6)
    a = _solve_case(args)
    b = _solve_case(args, n_tilt=60)
    assert abs(a.gdna_frac[0] - b.gdna_frac[0]) < 1e-6
    assert abs(a.rna_pos_frac[0] - b.rna_pos_frac[0]) < 1e-6


def test_a_delivered_row_is_evaluated_at_the_nodes_exactly():
    """The RNA level lanes deliver a row's INGREDIENTS (`CubeRow`: the held profiles, the slot's total
    and RNA opportunity, the lanes' coordinates) and ψ evaluates them at its own nodes — no lattice, no
    interpolation in θ: ψ with the row minus ψ without it IS the row's map (`_psi_reference.row_at`) on
    ψ's own tilt. A row with no profile changes nothing."""
    args, kw = _theta_case(50_000.0, 0.0, 0.5)
    u = kw["lam"]
    floor = -0.5 * np.maximum(0.0, (0.0 - u) / 0.3) ** 2
    row = SL.CubeRow(
        profile_pos=floor,
        profile_neg=None,
        u=u,
        total=400.0,
        opportunity=100.0,
        rho_ref=0.5,
    )
    bare, _fp, _fn, tau = psi_cube(*args, ambig=True, **kw)
    with_row, _fp, _fn, tau2 = psi_cube(*args, ambig=True, cube_rows={0: row}, **kw)
    assert np.array_equal(tau, tau2)
    # to rounding: the row is added before the quadrature's log-weights, so the difference of two
    # sums is not the row to the bit; an interpolated row would miss by 1e-2. The row's + profile is
    # also the witness that rules the pure − atom out (the last column), so that column is compared
    # apart: −∞ with the row, finite without
    assert np.allclose(
        (with_row - bare)[..., :-1], row_at(row, expit(u), tau[0])[:, :-1], atol=1e-9, rtol=0.0
    )
    assert np.all(np.isneginf(with_row[0, :, -1])) and np.all(np.isfinite(bare[0, :, -1]))
    assert not np.array_equal(with_row, bare), "the row must do something"
    empty = SL.CubeRow(None, None, u, 400.0, 100.0, 0.5)
    nothing, *_ = psi_cube(*args, ambig=True, cube_rows={0: empty}, **kw)
    assert np.array_equal(nothing, bare)


def test_without_strand_information_the_nodes_are_the_whole_domain():
    """κ = ½: the strand term is flat in the tilt, the window is the whole domain, and the nodes are
    uniform in θ over ``[−π/2, π/2]`` — the rule degrades to a lattice exactly where a lattice was right."""
    args, kw = _theta_case(1000.0, 0.2, 0.3)
    kw["kappa"] = 0.5
    _psi_, _fp, _fn, tau = psi_cube(*args, ambig=True, **kw)
    lattice = np.sin(np.linspace(-0.5 * np.pi, 0.5 * np.pi, SL._TILT_NODES))
    mixed = tau[..., : SL._TILT_NODES]  # the continuum's columns; the two atoms follow
    assert np.allclose(np.broadcast_to(lattice, mixed.shape), mixed, atol=1e-14)


# ── the tilt atom: the AMBIG tilt's hypothesis space is {pure +, pure −, mixed} ──────────────────────


def _pure_plus_case(n: float, fg_true: float, K: int = 101):
    """One strand-pure AMBIG slot — all of its RNA on + — with the expected counts and the variance
    frozen at the truth; prior-free."""
    args, _kw = _theta_case(n, fg_true, 1.0, K=K)
    return args


@pytest.mark.parametrize("fg_true", [0.5, 0.85, 0.97])
@pytest.mark.parametrize("n", [30.0, 300.0, 3000.0, 30000.0])
def test_a_strand_pure_slot_reads_its_gdna_at_the_strand_cap(n, fg_true):
    """At a slot whose RNA is all on one strand the truth sits AT the strand cap, and the pure
    hypothesis explains the split with no tilt parameter: the prior-free read-out lands within 0.07 of
    the truth at every depth from 30 to 30k fragments (the tilt continuum alone read 0.31–0.37 for a
    truth of 0.50, 0.69–0.80 for 0.85 and 0.85–0.95 for 0.97: every f_g below the cap fits the split
    with a slightly impure tilt, and the marginal's median sat below the cap)."""
    out = _solve_case(_pure_plus_case(n, fg_true))
    assert abs(float(out.gdna_frac[0]) - fg_true) < 0.07, (n, fg_true, float(out.gdna_frac[0]))


def test_the_three_tilt_hypotheses_carry_equal_reference_weight():
    """With no strand information (κ = ½) the strand term is flat in the tilt, so the three hypotheses'
    marginal masses must be equal at every λ: the mixed continuum's — its nodes across the whole domain
    with the trapezoid weights, normalised to the domain — equals each atom's. A continuum weighted by
    the window's length rather than its share of the domain would carry π times an atom's mass."""
    args, kw = _theta_case(1000.0, 0.2, 0.3)
    kw["kappa"] = 0.5
    psi, _fp, _fn, tau = psi_cube(*args, ambig=True, **kw)
    assert psi.shape[2] == SL._TILT_NODES + 2 and tau.shape[2] == SL._TILT_NODES + 2
    assert np.all(tau[..., -2] == 1.0) and np.all(tau[..., -1] == -1.0)
    mixed = lse(psi[..., :-2], axis=2)
    assert np.allclose(mixed, psi[..., -2], atol=1e-9, rtol=0.0)
    assert np.allclose(mixed, psi[..., -1], atol=1e-9, rtol=0.0)


def test_a_delivered_level_on_a_strand_rules_the_other_strands_pure_hypothesis_out():
    """A held RNA level on strand s is a certified witness that s carries RNA, so the hypothesis that ALL
    the slot's RNA is on the other strand is out (−∞ in ψ); a level on s says nothing against "pure s";
    with nothing delivered both atoms stand. And the witness reaches the read-out: at a strand-pure +
    slot with a truth of 0.50 a − level pulls f_g back below the cap the atom had recovered."""
    args, kw = _theta_case(3000.0, 0.5, 1.0)
    u = kw["lam"]
    floor = -0.5 * np.maximum(0.0, (0.0 - u) / 0.3) ** 2
    neg_level = SL.CubeRow(None, floor, u, 3000.0, 1000.0, 0.5)
    pos_level = SL.CubeRow(floor, None, u, 3000.0, 1000.0, 0.5)
    bare, *_ = psi_cube(*args, ambig=True, **kw)
    with_neg, *_ = psi_cube(*args, ambig=True, cube_rows={0: neg_level}, **kw)
    with_pos, *_ = psi_cube(*args, ambig=True, cube_rows={0: pos_level}, **kw)
    assert np.all(np.isfinite(bare[0, :, -2:]))
    assert np.all(np.isneginf(with_neg[0, :, -2])) and np.all(np.isfinite(with_neg[0, :, -1]))
    assert np.all(np.isneginf(with_pos[0, :, -1])) and np.all(np.isfinite(with_pos[0, :, -2]))
    unwitnessed = _solve_case(args, n_grid=41)
    witnessed = _solve_case(args, cube_rows={0: neg_level}, n_grid=41)
    assert float(witnessed.gdna_frac[0]) < float(unwitnessed.gdna_frac[0]) - 0.05
