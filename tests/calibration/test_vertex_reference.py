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

from rigel.calibration import simplex_logodds as SL
from rigel.calibration.simplex_logodds import (
    _compose,
    _logodds_grid,
    _solve_regions_logodds_all,
)

#: κ = ½ EXACTLY, no overdispersion, no fitted priors. On that substrate the strand term is bit-flat, so
#: the only things speaking about λ are the message under test and ψ's reference — which isolates the
#: reference without ablating anything.
#: ``n_grid_ss`` is deliberately FINE: the price law is read in λ, and a coarse lattice quantises the
#: very quantity being measured. The window ``L`` is the one real hard limit and G1 stays clear of it.
_BASE = dict(kappa=0.5, od_g=0.0, od_r=0.0, n_grid=128, L=10.0, n_tilt=64, n_grid_ss=4096)

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
    return _rows(lambda lam: -0.5 * float(p) * SL._log_fg(lam) ** 2)


def _msg_dn(p):
    """A claim written on ``log(1 − f_g)`` that ``f_g = 0`` — mode ``0`` — at precision ``p``."""
    return _rows(lambda lam: -0.5 * float(p) * SL._log1m_fg(lam) ** 2)


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

    ``_psi`` is CALLED, not reimplemented, so this cannot drift from what the solver
    computes (TRAPS: self-checking-validator)."""
    lam, fg = _logodds_grid(1024, 10.0)
    u_pos, u_neg, ap, an, _mu, _ms = _regions(1)
    psi, _fp, _fn, _tau = SL._psi(
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
        fg=fg,
        lam_logprior=(-0.5 * 1e3 * SL._log_fg(lam) ** 2)[None, :],
    )
    p = psi[0, :, 0]
    lo, hi = int(0.90 * lam.size), lam.size - 1
    slope = (p[hi] - p[lo]) / (lam[hi] - lam[lo])
    assert abs(slope + SL._JEFFREYS_REF) < 5e-3, (slope, -SL._JEFFREYS_REF)


# ── G3 — the two halves are orthogonal ──────────────────────────────────────────────────────────────


def _set_exponents(monkeypatch, c_g: float, c_r: float):
    """Re-write ψ's two arms with independent exponents.

    ``_gdna_arm`` / ``_rna_arm`` are looked up on the module at CALL time and nothing else imports the
    names, so patching the functions reaches BOTH the 1-D and the 2-D solver — TRAPS: an-ablation-that-never-ran's
    import-binding trap does not apply, and the assertions below would catch it if it did, because each
    one demands a MOVE rather than merely a difference."""
    monkeypatch.setattr(
        SL,
        "_gdna_arm",
        lambda lam, glp: (
            c_g * SL._log_fg(lam)[None, :]
            if glp is None
            else c_g * SL._log_fg(lam)[None, :] + np.asarray(glp, np.float64)
        ),
    )
    # The RNA replacement mirrors the gDNA one, fitted-prior socket included. Patching it with the
    # wrong arity raises rather than silently ignoring the arm, which is the outcome to want.
    monkeypatch.setattr(
        SL,
        "_rna_arm",
        lambda lam, rlp=None: (
            c_r * SL._log1m_fg(lam)[None, :]
            if rlp is None
            else c_r * SL._log1m_fg(lam)[None, :] + np.asarray(rlp, np.float64)
        ),
    )


def test_G3_each_half_of_the_constant_holds_ONE_vertex_and_is_NEGLIGIBLE_at_the_other(monkeypatch):
    """The orthogonality — what makes one fix cover both vertices.

    ``C·log(1−f_g)`` is the only term bounding ``f_g → 1``. Deleting it must move that answer by a full
    nat of log-odds while barely touching ``f_g → 0``, and symmetrically for the other half. Two claims
    per half, so neither can pass by accident.

    Scored in λ, not in ``f_g``: near a vertex ``f_g`` compresses a nat into ~1e-3, and a threshold on
    ``f_g`` would read a large move as a small one.

    Stated as a RATIO, not as bit-identity, and that distinction is a measurement rather than
    caution. On the toy's shipped lattice (``n_grid_ss = 256``, ``Δλ = 0.078``) the off-vertex half IS
    bit-identical — but at the 16× finer lattice used here it moves 0.015 nats, so the identity was
    lattice quantisation and not orthogonality (TRAPS: byte-identity-gate: a bit-identity gate has lied in both
    directions). The real claim is that each half is worth two orders of magnitude more at its own
    vertex than at the other, and that survives any resolution."""
    up, dn = _msg_up(1e3), _msg_dn(1e3)
    base_up = float(_lam(_solve(**up)[0]))
    base_dn = float(_lam(_solve(**dn)[0]))

    _set_exponents(monkeypatch, 0.5, 0.0)  # delete the f_g→1 bound
    own = float(_lam(_solve(**up)[0])) - base_up
    other = abs(float(_lam(_solve(**dn)[0])) - base_dn)
    assert own > 1.0, own  # it MOVED, toward its own vertex
    assert other < 0.1, other  # and barely registered at the other
    assert own > 20.0 * max(other, 1e-12), (own, other)

    _set_exponents(monkeypatch, 0.0, 0.5)  # delete the f_g→0 bound instead
    own = base_dn - float(_lam(_solve(**dn)[0]))
    other = abs(float(_lam(_solve(**up)[0])) - base_up)
    assert own > 1.0, own
    assert other < 0.1, other
    assert own > 20.0 * max(other, 1e-12), (own, other)


# ── G4 — the perturbation that makes G2 non-vacuous ─────────────────────────────────────────────────


def test_G4_softening_the_reference_moves_lambda_MONOTONICALLY_and_by_a_BOUNDED_step(monkeypatch):
    """The perturbation. ``λ* = ½·log(p/C)`` predicts that halving ``C`` buys only ``½·ln 2 = 0.35``
    nats — so softening the reference is a *weak* lever, and that is why "just lower the exponent" is
    not the fix (and why ``C → 0``, which would be, makes ψ improper — TRAPS: no-prior-means-haldane).

    Gated as monotone-and-bounded rather than as an exact step: the message's residual curvature and the
    other half's own term both contribute an offset that does not cancel across ``C``."""
    msg = _msg_up(1e4)
    seen = []
    for c in (0.5, 0.25, 0.125, 0.0625):
        _set_exponents(monkeypatch, 0.5, c)
        seen.append(float(_lam(_solve(**msg)[0])))
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


def _solve_composition(
    u_pos, u_neg, *, kappa, allow_pos, allow_neg, n_grid=60, n_grid_ss=256, n_tilt=None
):
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
        n_grid_ss=n_grid_ss,
        n_tilt=n_tilt,
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
    from rigel.calibration.simplex_logodds import (
        _psi,
        _logodds_grid,
        _lse,
        _posterior_median_fg,
    )

    lam, fg = _logodds_grid(60, _L)
    u_pos, u_neg = np.array([27.0, 9.0, 3.0]), np.array([3.0, 1.0, 0.0])
    psi = _psi(
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
        fg=fg,
    )[0][:, :, 0]
    post = np.exp(psi - _lse(psi, axis=1, keepdims=True))
    median = _posterior_median_fg(post, lam, fg)
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
    ``SUM = f_g + (1−f_g)/2``. `_compose` therefore restricts the share to the admissible strands itself
    rather than trusting the caller, so the single-strand path has no tilt to supply at all.

    Closure must hold for EVERY ``w_pos``, including ones that are wrong for the slot."""
    f_g = np.array([0.0, 0.05, 0.5, 0.95, 1.0])
    T, F = np.full(5, True), np.full(5, False)
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        for ap, an in ((T, T), (T, F), (F, T)):
            p, n = _compose(f_g, np.full(5, w), ap, an)
            assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, (w, ap[0], an[0])
            assert np.all(p >= 0.0) and np.all(n >= 0.0)
            assert np.all(p[~ap] == 0.0) and np.all(n[~an] == 0.0)
    # a single-strand slot takes the WHOLE RNA total on its admissible strand, whatever w_pos says
    for w in (0.0, 0.5, 1.0):
        p, n = _compose(f_g, np.full(5, w), T, F)
        assert np.allclose(p, 1.0 - f_g) and np.all(n == 0.0), w
        p, n = _compose(f_g, np.full(5, w), F, T)
        assert np.allclose(n, 1.0 - f_g) and np.all(p == 0.0), w
    # neither strand admissible ⇒ no RNA to place; the composition is f_g alone, which is a
    #   statement about opportunity and not a closure failure
    p, n = _compose(f_g, np.full(5, 0.5), F, F)
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
        p, n = _compose(f_g, np.full(3, w), T, T)
        assert np.all(p >= 0.0) and np.all(n >= 0.0), (w, p, n)
        assert np.max(np.abs(f_g + p + n - 1.0)) < 1e-15, w
    # the perturbation: an IN-range share is passed through untouched, so the clamp is not flattening
    #   the tilt into a constant
    p_lo, _ = _compose(f_g, np.full(3, 0.25), T, T)
    p_hi, _ = _compose(f_g, np.full(3, 0.75), T, T)
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
