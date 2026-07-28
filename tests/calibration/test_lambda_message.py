"""The ψ single-λ + θ message interface (`simplex_logodds`, the M6 rank-1 fix substrate).

The composition is ONE DOF ``λ = logit f_g``; the two-message per-component combine (a gDNA message on
``log f_g`` AND an RNA message on ``log f_R``) counts it TWICE (rank-1, ~2× over-confident,
`message_variance_derivation.md` §4). These pin the replacement: ψ accepts a SINGLE Gaussian on the grid
variable ``λ`` directly (both node classes) + a SEPARATE Gaussian on the tilt ``θ`` (AMBIG only, the genuinely
distinct strand-balance DOF). The three-stream relay that feeds these lives in `bp_solver` (next).
"""

from __future__ import annotations

import numpy as np
from scipy.special import expit

from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

_BASE = dict(kappa=0.5, od_g=0.0, od_r=0.0, n_grid=128, L=10.0, n_tilt=64)


def _nodes():
    # 3 single-strand (+) nodes + 1 AMBIG node, all with data, UNSTRANDED (κ=½ ⇒ the strand is flat, so the
    # λ-message alone drives f_g — a clean isolation of the message channel).
    u_pos = np.array([100.0, 100.0, 100.0, 100.0])
    u_neg = np.array([0.0, 0.0, 0.0, 100.0])
    ap = np.array([True, True, True, True])
    an = np.array([False, False, False, True])
    return u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(4)


def test_lambda_message_pulls_f_g_to_sigma_lambda():
    """A high-precision λ-message pulls each node's f_g to ``σ(λ_target)`` — the direct-on-λ constraint."""
    u_pos, u_neg, ap, an, mu, ms = _nodes()
    lam_target = np.array([-2.0, 0.0, 2.0, 1.0])
    dc = _solve_nodes_logodds_all(
        u_pos,
        u_neg,
        ap,
        an,
        mu,
        ms,
        lam_imp_mode=lam_target,
        lam_imp_prec=np.full(4, 50.0),
        **_BASE,
    )
    for i in range(4):
        assert abs(float(dc.gdna_frac[i]) - float(expit(lam_target[i]))) < 0.02, (
            i,
            dc.gdna_frac[i],
            expit(lam_target[i]),
        )
    assert not np.any(np.isnan(dc.gdna_frac))


def test_lambda_message_zero_precision_is_a_noop():
    """A zero-precision λ-message must not move the solve (no message)."""
    u_pos, u_neg, ap, an, mu, ms = _nodes()
    dc0 = _solve_nodes_logodds_all(u_pos, u_neg, ap, an, mu, ms, **_BASE)
    dc1 = _solve_nodes_logodds_all(
        u_pos,
        u_neg,
        ap,
        an,
        mu,
        ms,
        lam_imp_mode=np.full(4, 3.0),
        lam_imp_prec=np.zeros(4),
        **_BASE,
    )
    assert np.allclose(dc0.gdna_frac, dc1.gdna_frac)


def test_theta_message_tilts_ambig_node():
    """A θ-message (AMBIG only) shifts the strand tilt: θ>0 ⇒ f_pos ≫ f_neg, at fixed f_g (via the λ-message)."""
    u_pos, u_neg, ap, an, mu, ms = _nodes()
    lam = np.array([-2.0, 0.0, 2.0, 1.0])
    common = dict(lam_imp_mode=lam, lam_imp_prec=np.full(4, 50.0), **_BASE)
    dc_flat = _solve_nodes_logodds_all(u_pos, u_neg, ap, an, mu, ms, **common)
    dc_tilt = _solve_nodes_logodds_all(
        u_pos,
        u_neg,
        ap,
        an,
        mu,
        ms,
        theta_imp_mode=np.array([0.0, 0.0, 0.0, 1.2]),
        theta_imp_prec=np.array([0.0, 0.0, 0.0, 50.0]),
        **common,
    )
    ex = 3  # the AMBIG node
    assert dc_tilt.rna_pos_frac[ex] > dc_flat.rna_pos_frac[ex]
    assert dc_tilt.rna_neg_frac[ex] < dc_flat.rna_neg_frac[ex]
    # f_g (the λ axis) is unchanged by the θ message — the tilt is an ORTHOGONAL DOF.
    assert abs(float(dc_tilt.gdna_frac[ex]) - float(dc_flat.gdna_frac[ex])) < 0.02
    assert not np.any(np.isnan(dc_tilt.rna_pos_frac))
