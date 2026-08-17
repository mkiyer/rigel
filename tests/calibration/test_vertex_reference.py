"""ψ PRICES THE SIMPLEX VERTEX AT HALF A NAT PER NAT OF PRECISION — the mechanism, pinned as numbers.

The composition posterior lives on ``λ = logit f_g``. ψ's reference is ONE constant used twice with
opposite signs (``simplex_logodds._JEFFREYS_REF``): ``+C·log f_g`` bounds the ``f_g → 0`` vertex and
``+C·log(1−f_g)`` bounds ``f_g → 1``. On the λ axis that pair is exactly Beta(½,½), whose tail is
``−C·|λ|``.

⭐⭐ Every composition MESSAGE is written on ``log f_c``, and near the vertex ``log f_c ≈ −e^{−|λ|}``. So a
Gaussian message of precision ``p`` contributes ``−½p·e^{−2|λ|}``, whose gradient decays exponentially
while the reference's stays flat at ``C``. Setting them equal:

    λ*  =  ½ · log(p / C)

⭐⭐⭐ **The vertex is not forbidden — it is PRICED, and the price is exponential.** One extra nat of
log-odds costs ``e²`` ≈ 7.4× the precision; one decade of precision buys ``½·ln 10 ≈ 1.15`` nats. Real
objects on the toy carry ``p ≈ 6–14``, i.e. ``λ* ≈ 1.7``, so ``f_g ≈ 0.85`` before the λ-channel tops it
up — and reaching ``f_g = 0.999`` would need ``p ≈ 5e5``, some 40,000× more precision than exists.
⭐ That is why DEPTH cannot fix it: ``p`` grows at best linearly in ``n``, so ``λ`` grows as ``½ log n``.

⛔ **AND THE ROOT CAUSE IS NOT THE REFERENCE, IT IS A MISSING CHANNEL** — G6. ``f_g = 1`` is a statement
about ``ρ_r = 0``, and no gDNA-side channel can make it however precise it is. The one channel that CAN
is the CERTIFIED-SPLICE mass (``boundary_spliced``) — a spliced fragment cannot be gDNA — whose value at a
silent gene is exactly zero, the strongest possible evidence for the vertex. ψ never reads it.
⛔ There are TWO components, gDNA and RNA; spliced-vs-unspliced is a property of the FRAGMENT, never a
second RNA species.

===  =========================================================================================
TRAPS: no-magic-numbers   the price law — ``λ`` buys ~½·ln10 nats per DECADE of precision on ``log f_c``
G1b  the mirror at the other vertex
TRAPS: one-thing-varied   ψ's slope in the vertex tail is exactly ``−C``
TRAPS: converge-and-delete   the two halves are ORTHOGONAL — each is worth >20x more at its own vertex than at the other
TRAPS: the-source-does-not-cite-docs   halving ``C`` moves ``λ`` by a bounded, sub-nat step  ⭐ the perturbation for TRAPS: one-thing-varied
TRAPS: no-enumeration-without-a-census   the SAME claim delivered on ``λ`` is honoured AT PAR  ⭐ the control for TRAPS: no-magic-numbers
G6   ⛔ ψ is BLIND to the certified-RNA channel — the observable a fix must move
===  =========================================================================================

⭐ TRAPS: no-enumeration-without-a-census is why G1 is a falsification and not a truism: without it, G1 reads as "a grid cannot represent a
vertex", which is false — the grid reaches ``σ(±10) = 4.5e-5`` and a λ-message gets there at par. The
COORDINATE is the price, not the lattice.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration import simplex_logodds as SL
from rigel.calibration.simplex_logodds import _logodds_grid, _solve_regions_logodds_all

#: κ = ½ EXACTLY, no overdispersion, no fitted priors. On that substrate the strand term is bit-flat, so
#: the only things speaking about λ are the message under test and ψ's reference — which isolates the
#: reference without ablating anything.
#: ⚠ ``n_grid_ss`` is deliberately FINE: the price law is read in λ, and a coarse lattice quantises the
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


def _msg_up(p):
    """A gDNA message claiming ``f_g = 1`` — mode ``log f_g = 0`` — at precision ``p``."""
    return dict(gdna_imp_mode=np.zeros(2), gdna_imp_prec=np.full(2, float(p)))


def _msg_dn(p):
    """An RNA message claiming ``f_g = 0`` — mode ``log(1−f_g) = 0`` — at precision ``p``."""
    return dict(
        rna_imp_mode=(np.zeros(2), np.zeros(2)),
        rna_imp_prec=(np.full(2, float(p)), np.zeros(2)),
    )


# ── TRAPS: no-magic-numbers — the price law ──────────────────────────────────────────────────────────────────────────────


def test_G1_a_log_f_g_message_buys_log_odds_only_LOGARITHMICALLY_in_its_precision():
    """⭐⭐⭐ THE PRICE LAW. Sweep the message precision over four decades and regress ``λ`` on
    ``log₁₀ p``: the slope must be ``½·ln 10`` — the vertex costs ``e²`` in precision per nat.

    ⛔ The consequence is the finding: an object with the precision the toy actually measures
    (``p ≈ 10``) sits at ``λ ≈ 1.5``, i.e. ``f_g ≈ 0.82``, and no achievable depth closes that, because
    ``p`` is at best linear in the count while ``λ`` is logarithmic in ``p``.

    ⚠ Kept clear of ``p ≥ 1e8``, where ``λ`` hits the ``L = 10`` window and the law stops being the
    binding constraint — that ceiling is a separate fact and TRAPS: no-enumeration-without-a-census is what shows it is not the cause."""
    lam_of = {e: float(_lam(_solve(**_msg_up(10.0**e))[0])) for e in (2, 3, 4, 5, 6)}
    es = sorted(lam_of)
    slope = float(np.polyfit(es, [lam_of[e] for e in es], 1)[0])
    assert slope == pytest.approx(_NATS_PER_DECADE, rel=0.15), (slope, _NATS_PER_DECADE, lam_of)
    # ⭐ and the practical statement: at the precision real objects carry, the answer is far from 1.
    assert float(_solve(**_msg_up(10.0))[0]) < 0.95, _solve(**_msg_up(10.0))[0]


def test_G1b_the_price_law_is_the_SAME_at_the_other_vertex():
    """The mirror. ``f_g → 0`` is bounded by the gDNA half of the same constant, so the slope in
    ``|λ|`` must match TRAPS: no-magic-numbers' — one constant, two signs (which TRAPS: converge-and-delete then pins as orthogonal)."""
    lam_of = {e: float(_lam(_solve(**_msg_dn(10.0**e))[0])) for e in (2, 3, 4, 5, 6)}
    es = sorted(lam_of)
    slope = float(np.polyfit(es, [-lam_of[e] for e in es], 1)[0])
    assert slope == pytest.approx(_NATS_PER_DECADE, rel=0.15), (slope, _NATS_PER_DECADE, lam_of)


# ── TRAPS: one-thing-varied — the reference is the only slope in the tail ─────────────────────────────────────────────────


def test_G2_psi_slope_in_the_vertex_tail_is_exactly_minus_the_reference_exponent():
    """⭐⭐ ψ built directly on the λ grid with a ``log f_g`` message: every non-reference term is flat in
    the far tail (κ=½ ⇒ bit-flat strand; the message's gradient decays as ``e^{−2λ}``), so ``dψ/dλ``
    must equal ``−_JEFFREYS_REF``. The whole mechanism in one number.

    ⛔ ``_local_loglik_logodds`` is CALLED, not reimplemented, so this cannot drift from what the solver
    computes (TRAPS: self-checking-validator)."""
    lam, fg = _logodds_grid(1024, 10.0)
    u_pos, u_neg, ap, an, _mu, _ms = _regions(1)
    psi, _fp, _fn = SL._local_loglik_logodds(
        u_pos,
        u_neg,
        ap,
        an,
        0.5,
        0.0,
        0.0,
        lam,
        fg,
        np.full(1, 0.5),
        np.full(1, 0.5),
        np.zeros(1),
        gdna_imp_mode=np.zeros(1),
        gdna_imp_prec=np.full(1, 1e3),
    )
    p = psi[0]
    lo, hi = int(0.90 * lam.size), lam.size - 1
    slope = (p[hi] - p[lo]) / (lam[hi] - lam[lo])
    assert abs(slope + SL._JEFFREYS_REF) < 5e-3, (slope, -SL._JEFFREYS_REF)


# ── TRAPS: converge-and-delete — the two halves are orthogonal ──────────────────────────────────────────────────────────────


def _set_exponents(monkeypatch, c_g: float, c_r: float):
    """Re-write ψ's two arms with independent exponents.

    ⚠ ``_gdna_arm`` / ``_rna_arm`` are looked up on the module at CALL time and nothing else imports the
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
    # ⭐ The RNA replacement mirrors the gDNA one, fitted-prior socket included. It took no prior until
    # 2026-08-15, when `_rna_arm` gained the parameter its twin always had; patching it with the old
    # one-argument shape raised rather than silently ignoring the arm, which is the outcome to want.
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
    """⭐⭐⭐ THE ORTHOGONALITY — what makes ONE fix cover BOTH vertices.

    ``C·log(1−f_g)`` is the only term bounding ``f_g → 1``. Deleting it must move that answer by a full
    nat of log-odds while barely touching ``f_g → 0``, and symmetrically for the other half. Two claims
    per half, so neither can pass by accident.

    ⭐ Scored in λ, not in ``f_g``: near a vertex ``f_g`` compresses a nat into ~1e-3, and a threshold on
    ``f_g`` would read a large move as a small one.

    ⛔ Stated as a RATIO, not as bit-identity, and that distinction is a measurement rather than
    caution. On the toy's shipped lattice (``n_grid_ss = 256``, ``Δλ = 0.078``) the off-vertex half IS
    bit-identical — but at the 16× finer lattice used here it moves 0.015 nats, so the identity was
    lattice quantisation and not orthogonality (TRAPS: byte-identity-gate: a bit-identity gate has lied in both
    directions). The real claim is that each half is worth **two orders of magnitude more** at its own
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


# ── TRAPS: the-source-does-not-cite-docs — the perturbation that makes TRAPS: one-thing-varied non-vacuous ─────────────────────────────────────────────────


def test_G4_softening_the_reference_moves_lambda_MONOTONICALLY_and_by_a_BOUNDED_step(monkeypatch):
    """⭐⭐ THE PERTURBATION. ``λ* = ½·log(p/C)`` predicts that halving ``C`` buys only ``½·ln 2 = 0.35``
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
    # ⭐ four halvings — a 16x softening — buy less than the ONE decade of precision G1 measured.
    assert (seen[-1] - seen[0]) < 4.0 * _NATS_PER_DECADE, seen


# ── TRAPS: no-enumeration-without-a-census — the CONTROL: on λ, the same claim is honoured at par ────────────────────────────────────────


def test_G5_the_SAME_claim_delivered_on_lambda_is_honoured_AT_PAR():
    """⭐⭐⭐ THE CONTROL THAT MAKES TRAPS: no-magic-numbers A FALSIFICATION (TRAPS: a-gate-that-already-passed).

    ``f_g = 1 − 1e-4`` is ``λ = +9.21``, which the ``L = 10`` grid represents exactly. Delivered on the λ
    axis the claim is honoured; delivered on ``log f_g`` — the SAME claim, the same precision — it is
    not. So the vertex is reachable and the lattice is not the limit: **the coordinate is the price.**

    ⛔ This is what forbids "the grid cannot represent a vertex" as a diagnosis."""
    target = 1.0 - 1e-4
    lam_t = float(_lam(np.array([target]))[0])
    on_lambda = float(_solve(lam_imp_mode=np.full(2, lam_t), lam_imp_prec=np.full(2, 1e3))[0])
    on_log_fg = float(_solve(**_msg_up(1e3))[0])
    assert on_lambda > 0.999, on_lambda
    assert on_log_fg < 0.995, on_log_fg
    assert float(_lam(np.array([on_lambda]))[0] - _lam(np.array([on_log_fg]))[0]) > 2.0


# ── G6 — the root cause: the certified-RNA channel is not read ───────────────────────────────────────


def test_G6_psi_is_BLIND_to_the_certified_RNA_channel():
    """⛔⛔ THE ROOT CAUSE, AS AN OBSERVABLE — and the one a fix MUST move (TRAPS: name-the-observable-per-site).

    ``f_g = 1`` is a claim about ``ρ_r = 0``. No gDNA-side message can establish it at any precision (TRAPS: no-magic-numbers
    prices exactly that attempt). The channel that CAN is the CERTIFIED SPLICE — ``boundary_spliced``, the
    fragments that crossed this boundary contiguously having spliced elsewhere. A spliced fragment cannot be
    gDNA, so this channel has no gDNA term at all; at a silent gene it is exactly zero, which is the
    strongest available evidence for the vertex.

    ⛔ RNA IS RNA. There are two components. Spliced-vs-unspliced is a property of the FRAGMENT (spliced
    ⇒ certified RNA, unspliced ⇒ must be deconvolved), never a second species of RNA.

    ⚠ ``boundary_spliced`` is BOUNDARY-only and is structurally ZERO on every toy geometry (a toy's exons ARE its
    regions, so no spliced molecule crosses an interior boundary contiguously). On the real panel it carries
    85.3 %% of BOUNDARY unspliced mass off capture. So this gate uses ``mass_spliced`` as a pure SWEEP HANDLE
    to prove ψ ignores the channel — it is not a claim about the channel's size.

    ψ does not read it: sweeping ``mass_spliced`` from 0 to 1e5 at fixed strand counts leaves ``f_g``
    **bit-identical**. ⚠ That exclusion is deliberate and its stated reason is sound as far as it goes —
    the spliced mass must not be double-counted into the *unspliced* likelihood. But "do not add it to
    that likelihood" is not "discard the information that it is zero", and this gate is what would fire
    if the second were ever addressed.

    ⭐ If this test starts FAILING, that is the fix landing — re-point it at the new channel, do not
    widen it.

    ⛔⛔ **BUT READ ``test_certified_rna_licence.py`` FIRST — the obvious fix is REFUTED, and one clause
    above is wrong.** "At a silent gene it is exactly zero, which is the strongest available evidence for
    the vertex" does not hold: the chance a crossing RNA fragment shows a *visible* splice, ``q``, is a
    free nuisance, and with it free ``S = 0`` is perfectly explained at ``f_g = 1`` too (take ``q → 0``).
    ``test_certified_rna_licence.py``'s zero-spliced-count gate pins that. The information in this channel lives entirely in ``S > 0`` and points AWAY
    from this vertex, so **nothing here rescues the vertex population.** And the coefficient that would
    make this gate fail — ``(½ + S)`` on the RNA arm — is measured WORSE than the uninformative reference
    on 12 of the 36 ladder conditions (gates ``TRAPS: pure-and-length-censored``/``TRAPS: divide-by-a-probability``). What survives is the algebra, not the
    coefficient: reference + term is exactly Beta(½, ½+S), gate ``TRAPS: opposite-tilts-must-not-pool``."""
    quiet = float(_solve(mass_spliced=np.zeros(2), **_msg_up(10.0))[0])
    loud = float(_solve(mass_spliced=np.full(2, 1e5), **_msg_up(10.0))[0])
    assert loud == pytest.approx(quiet, abs=0.0), (quiet, loud)
    # ⭐ and the object it would matter for is nowhere near the vertex without it.
    assert quiet < 0.95, quiet
