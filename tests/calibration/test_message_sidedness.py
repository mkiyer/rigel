"""STAGE 4d — PER-MESSAGE SIDEDNESS on ψ's RNA channel, falsified clause by clause.

A one-sided RNA claim ("the destination holds AT LEAST this much RNA") penalises only the
contradiction side — `_rna_residual`'s clamp — and until now that was a PROCESS-GLOBAL toggle
(`ONE_SIDED_RNA`, built for the §1.1 arm). Stage 4 needs it PER SLOT: a complete flank's transfer is
a two-sided estimate, an incomplete flank's is the one-sided bound, and both kinds coexist in one
solve. ``rna_one_sided`` rides the message: ``None`` is byte-identical to today, a per-slot bool
selects the clamp slot by slot, and all-True must be bit-identical to the global flag.
"""

from __future__ import annotations

import numpy as np

import rigel.calibration.simplex_logodds as SL
from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

_L = 10.0


def _solve(u_pos, u_neg, allow_pos, allow_neg, rna_mode=None, rna_prec=None, one_sided=None):
    d = _solve_regions_logodds_all(
        np.asarray(u_pos, np.float64),
        np.asarray(u_neg, np.float64),
        np.asarray(allow_pos, bool),
        np.asarray(allow_neg, bool),
        np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64),
        np.zeros(len(u_pos), np.float64),
        kappa=0.95,
        od_g=0.0,
        od_r=0.0,
        n_grid=60,
        n_grid_ss=256,
        n_tilt=None,
        L=_L,
        rna_imp_mode=rna_mode,
        rna_imp_prec=rna_prec,
        rna_one_sided=one_sided,
    )
    return np.asarray(d.gdna_frac, np.float64)


#: two identical single-strand slots holding plenty of RNA, plus two identical AMBIG slots whose
#: strong sense excess reads as RNA+ at κ = 0.95 — so the "at least 30 % RNA+" bound is genuinely
#: SATISFIED on both solver paths, and the per-slot mask is the ONLY thing that can make paired
#: slots answer differently.
U_POS = [80.0, 80.0, 76.0, 76.0]
U_NEG = [4.0, 4.0, 4.0, 4.0]
ALLOW_P = [True, True, True, True]
ALLOW_N = [False, False, True, True]
#: an RNA+ claim well BELOW what the slots' own data support — satisfied as a bound, wrong as an
#: equality — at a precision strong enough to drag a two-sided slot visibly.
MODE = np.log(0.30)
PREC = 50.0


def _msg(prec_mask):
    n = len(U_POS)
    mode = (np.full(n, MODE), np.full(n, MODE))
    prec = (np.where(prec_mask, PREC, 0.0), np.zeros(n))
    return mode, prec


def test_none_is_byte_identical_to_all_false():
    """``rna_one_sided=None`` and an all-False mask are the SAME solve, bit for bit — the field's
    default costs nothing anywhere."""
    mode, prec = _msg(np.ones(4, bool))
    a = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N, mode, prec, one_sided=None)
    b = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N, mode, prec, one_sided=np.zeros(4, bool))
    np.testing.assert_array_equal(a, b)


def test_per_slot_sidedness_on_both_solver_paths():
    """Paired identical slots, one two-sided and one one-sided: the satisfied bound must leave the
    one-sided slot essentially where the no-message solve puts it, while the two-sided twin is
    dragged toward the claim — on the single-strand path AND the AMBIG cube."""
    mode, prec = _msg(np.ones(4, bool))
    sided = np.array([False, True, False, True])
    got = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N, mode, prec, one_sided=sided)
    base = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N)
    for two, one in ((0, 1), (2, 3)):
        pulled = abs(got[two] - base[two])
        held = abs(got[one] - base[one])
        assert pulled > 0.05, f"slots {two}/{one}: the two-sided claim must actually drag"
        assert held < pulled / 5.0, (
            f"slots {two}/{one}: the satisfied one-sided bound must not drag "
            f"(held {held:.4f} vs pulled {pulled:.4f})"
        )


def test_all_true_is_bit_identical_to_the_global_flag():
    """The per-slot mask GENERALIZES the global toggle: all-True with the flag down must equal
    None with the flag up, bit for bit — one mechanism, two spellings, no drift between them."""
    mode, prec = _msg(np.ones(4, bool))
    a = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N, mode, prec, one_sided=np.ones(4, bool))
    assert not SL.ONE_SIDED_RNA[0]
    SL.ONE_SIDED_RNA[0] = True
    try:
        b = _solve(U_POS, U_NEG, ALLOW_P, ALLOW_N, mode, prec, one_sided=None)
    finally:
        SL.ONE_SIDED_RNA[0] = False
    np.testing.assert_array_equal(a, b)


def test_the_contradiction_side_still_penalises():
    """One-sidedness is not a mute: a claim ABOVE what the slot holds ('at least 98 % RNA' at a slot
    with real gDNA evidence) contradicts the data and must drag essentially like the two-sided form —
    NOT exactly, because the two-sided form also penalises the few grid points above the mode, so the
    one-sided posterior keeps slightly MORE of its high-RNA tail (a strict inequality, asserted)."""
    u_pos, u_neg = [10.0, 10.0], [1.0, 1.0]
    ap, an = [True, True], [False, False]
    n = 2
    mode = (np.full(n, np.log(0.98)), np.full(n, np.log(0.98)))
    prec = (np.full(n, PREC), np.zeros(n))
    base = _solve(u_pos, u_neg, ap, an)
    two = _solve(u_pos, u_neg, ap, an, mode, prec, one_sided=np.zeros(n, bool))
    one = _solve(u_pos, u_neg, ap, an, mode, prec, one_sided=np.ones(n, bool))
    assert abs(two[0] - base[0]) > 0.05, "the contradicting claim must drag the two-sided slot"
    np.testing.assert_allclose(one, two, rtol=0.05)
    assert (one <= two).all(), "one-sided keeps more of the high-RNA side, never less"
