"""Closed-form facts of the enrichment-frame primitives (phase E1), pinned against ALGEBRAIC GROUND TRUTH
rather than against solver behaviour.

Every function in :mod:`rigel.calibration.enrichment_frame` is a pure identity, so each test here is a
closed-form check that cannot drift with the solver — the same contract as ``test_message_frames.py``. The two
load-bearing identities the whole framework rests on (``enrichment_ratio_generalization.md``):

* the composed total-density identity ``ρ_tot = M·(k+1)/(k·E_g + E_r)`` (§1 + §4), which makes
  ``total_density ∘ f_g_from_k ∘ (mass identity)`` mutually exact;
* the r₂ ``(k+1)`` cancellation (§5b, ``junction_enrichment_scaling.md``), verified there to 5.8e−16 over
  20 000 draws — reproduced here as an exact algebraic equality.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import polygamma

from rigel.calibration.enrichment_frame import (
    composition_logvar,
    graft_frame_logvar,
    graft_premise_logvar,
    mismatch_deflate,
    mismatch_gap,
    peel_continue_share,
    peel_rna_logvar,
    conservation_rescale,
    peel_share_logvar,
    residual_level,
    transfer_logvar,
)

# Real-ish effective lengths from the FL models (gDNA ~N(300,60), RNA ~N(200,50)):
# contained regions grow with L; boundary crossings saturate at fl_mean.
EG_REG, ER_REG = 2701.0, 2801.0  # a long contained region (L = 3000)
EG_BND, ER_BND = 300.0, 200.0  # a boundary crossing (E_g = fl_mean_gdna, E_r = fl_mean_rna)


# ---------------------------------------------------------------------------
# §1 / §4 — the composed total-density identity, mutually exact.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# §5b — the r₂ (k+1) cancellation, as an exact algebraic equality.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# §2 — the bounding lemma + the f_g=1 degeneracy (the two estimators coincide there).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# The composition-as-variance gate (decision 1) — constant-free, self-excluding.
# ---------------------------------------------------------------------------


def test_composition_logvar_is_pure_counting_at_the_structural_corner():
    """At ``f_g = 1`` with ``Var(f_g) = 0`` (a structural gDNA node) the composition term vanishes and
    ``Var(log ρ_tot) = 1/n`` — the honest counting precision, no composition penalty."""
    v = float(composition_logvar(1.0, EG_BND, ER_BND, var_fg=0.0, n=25.0))
    assert v == pytest.approx(1.0 / 25.0, rel=1e-14)


def test_composition_logvar_self_excludes_short_regions_without_a_threshold():
    """The composition coefficient ``[(1/E_g − 1/E_r)/B]²`` grows as ``E_g`` collapses on a short region, so a
    short region carries a LARGE composition variance and down-weights itself as an enrichment reference —
    exactly §2's exclusion, with NO ``L ≲ fl_mean`` cut. Ordering: long region ≪ boundary ≪ short region."""
    # isolate the composition term (n→inf ⇒ counting term 0), at the worst-case f_g = ½
    long_reg = float(composition_logvar(0.5, EG_REG, ER_REG, var_fg=1.0, n=np.inf))
    boundary = float(composition_logvar(0.5, EG_BND, ER_BND, var_fg=1.0, n=np.inf))
    short_reg = float(composition_logvar(0.5, 24.4, 101.4, var_fg=1.0, n=np.inf))
    assert long_reg < boundary < short_reg
    # the boundary coefficient is (1/300 − 1/200)/B with B = ½(1/300+1/200) ⇒ (−1/600)/(1/240) = −0.4 ⇒ 0.16
    assert boundary == pytest.approx(0.16, rel=1e-6)


# ---------------------------------------------------------------------------
# §5 — the structural refusals (the one genuine bool).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# The UNIFIED SOLVER theorem (unified_solver_design.md §2): reframe + density-mode ÷ M_dst subsumes the shift,
# is enrichment-invariant, and handles a PARTIAL (set-mismatched) source correctly where the shift cannot.
# ---------------------------------------------------------------------------


def _shift_logfrac(rho_c_src, E_c_src_frame, imputed_masses):
    """The OLD shift factor on log f_c: log( M_c / Σ M_c' ), M_c imputed with the dst's E_c. Reference to
    compare the unified density mode against."""
    M_c = rho_c_src * E_c_src_frame
    return np.log(M_c / np.sum(imputed_masses))


# ═══════════════════════════════════════════════════════════════════════════════════════════════════════════
# The pass-0 message-VARIANCE laws (message_variance_derivation.md M1-M5) — closed-form + MC cross-check.
# The MC ground-truth harness is scratchpad/message_variance_mc.py (independently re-derived + adversarially
# verified, workflow wf_c952640d, <1% in-regime); these pin the arithmetic and one MC per non-trivial law.
# ═══════════════════════════════════════════════════════════════════════════════════════════════════════════

_MC = np.random.default_rng(20260724)


def _lognormal_mc(mean, var_log, size):
    s = np.sqrt(max(float(var_log), 0.0))
    return float(mean) * np.exp(_MC.normal(-0.5 * s * s, s, size=size))


def _beta_mc(mean, var, size):
    m = float(mean)
    v = min(float(var), m * (1.0 - m) * 0.98)
    c = m * (1.0 - m) / v - 1.0
    return _MC.beta(m * c, (1.0 - m) * c, size=size)


# ── M1 transport seed ──


# ── M2 graft (SUM, share-weighted) ──


# ── M3 peel (DIFFERENCE, u-weighted) ──


def test_peel_rna_is_u_weighted_difference():
    """M3: Var(log ρ_ν) = u²·(Var(log ρ_R)+σ²_t) + (u−1)²·v_μ. u=1 (all continues) ⇒ just Var(log T)."""
    assert float(peel_rna_logvar(0.005, 0.01, 0.002, 1.0)) == pytest.approx(0.015, rel=1e-14)
    # u=3: 9·(0.005+0.01) + 4·0.002
    assert float(peel_rna_logvar(0.005, 0.01, 0.002, 3.0)) == pytest.approx(
        9.0 * 0.015 + 4.0 * 0.002, rel=1e-14
    )


def test_peel_rna_matches_mc_in_regime():
    """The u-weighted DIFFERENCE equals the empirical Var(log ρ_ν) at low u (ε≲0.15, the valid regime);
    σ²_transfer is load-bearing."""
    rho_R_x, var_log_rhoR, r, var_log_r, rho_mu, n_s = 40.0, 1.0 / 5000, 200.0, 0.004, 0.10, 1500.0
    Rx = _lognormal_mc(rho_R_x, var_log_rhoR, 400_000)
    rr = _lognormal_mc(r, var_log_r, 400_000)
    rm = _lognormal_mc(rho_mu, 1.0 / n_s, 400_000)
    nu = Rx / rr - rm
    keep = nu > 0
    emp = float(np.var(np.log(nu[keep])))
    T0 = rho_R_x / r
    u = T0 / (T0 - rho_mu)
    pred = float(peel_rna_logvar(var_log_rhoR, var_log_r, 1.0 / n_s, u))
    assert pred == pytest.approx(emp, rel=0.12)


# ── M5 transfer variance (direction-dependent) ──


def test_transfer_logvar_cancels_on_graft_and_adds_on_peel():
    """M5: σ²_transfer = 0 on the GRAFT (r common-mode, cancels), = Var(log ρ_tot^dst)+Var(log ρ_tot^src) on
    the PEEL. Vectorized over a mixed direction mask."""
    dst = np.array([0.1, 0.1, 0.1])
    src = np.array([0.3, 0.3, 0.3])
    graft = np.array([True, False, True])
    out = transfer_logvar(dst, src, graft)
    assert list(out) == [0.0, 0.4, 0.0]


def test_transfer_logvar_feeds_off_composition_logvar():
    """The two inputs are exactly composition_logvar at the two endpoints — one law, reusing the existing
    total-density variance (no new derivation)."""
    vd = float(composition_logvar(0.6, EG_BND, ER_BND, var_fg=0.01, n=100.0))
    vs = float(composition_logvar(0.9, EG_REG, ER_REG, var_fg=0.005, n=50.0))
    assert float(transfer_logvar(vd, vs, graft=False)) == pytest.approx(vd + vs, rel=1e-14)


# ── M8 the graft's frame-mislift variance ──


def test_graft_frame_logvar_is_zero_without_a_frame_change():
    """M8's defining limit: no frame step ⇒ no mislift. This is what makes the term inert off-capture, where
    the shipped graft is measured to be EXACT (required correction log c = +0.009/−0.008/+0.054)."""
    assert float(graft_frame_logvar(1.0)) == 0.0
    assert list(graft_frame_logvar(np.array([1.0, 1.0]))) == [0.0, 0.0]


def test_graft_frame_logvar_is_the_squared_log_step_and_direction_free():
    """``(log r)²`` — the method-of-moments second moment of a single observation of the un-cancelled frame
    step. Symmetric in r ↔ 1/r: a depletion mislifts exactly as badly as an equal enrichment."""
    assert float(graft_frame_logvar(np.e)) == pytest.approx(1.0, rel=1e-14)
    assert float(graft_frame_logvar(6.1)) == pytest.approx(np.log(6.1) ** 2, rel=1e-14)
    assert float(graft_frame_logvar(4.0)) == pytest.approx(float(graft_frame_logvar(0.25)), rel=1e-14)


def test_graft_frame_logvar_guards_a_degenerate_ratio():
    """A node with no frame (r ≤ 0 — no mass, no ρ_tot) must give 0, not a nan: the relay passes such a
    message through at r = 1, so there is no mislift to charge."""
    out = graft_frame_logvar(np.array([0.0, -1.0, 2.0]))
    assert np.all(np.isfinite(out))
    assert list(out[:2]) == [0.0, 0.0]


# ── M4 ÷M_dst precision conversion ──


# ── M7 the cross-cliff COMPOSITION-MISMATCH variance (DerSimonian–Laird) ──


def test_mismatch_gap_is_the_log_ratio_and_flags_only_one_sided_absence():
    """M7: G = log(ρ^msg/ρ^own). ``contradicted`` marks exactly one side absent (an assertion of ``f_c = 0``
    against a node that has the component, or vice versa); BOTH absent is not a contradiction, it is silence."""
    g, c = mismatch_gap(np.array([2.0, 0.0, 5.0, 0.0]), np.array([0.5, 3.0, 0.0, 0.0]))
    assert g[0] == pytest.approx(np.log(4.0), rel=1e-14)
    assert list(c) == [False, True, True, False]
    assert g[3] == 0.0  # both silent ⇒ no gap, nothing to price


def test_mismatch_deflate_is_the_closed_form_and_never_strengthens():
    """M7: p_eff = 1/max(v_msg, G²−v_own) — and a deflation can only ever REDUCE a precision."""
    p = np.array([25.0, 25.0, 25.0])
    gap = np.array([0.0, 0.5, 2.7])
    v_own = np.full(3, 0.56)
    out = mismatch_deflate(p, gap, np.zeros(3, bool), v_own)
    want = 1.0 / np.maximum(1.0 / p, gap * gap - v_own)
    assert out == pytest.approx(want, rel=1e-12)
    assert np.all(out <= p + 1e-12)


def test_mismatch_deflate_is_inert_without_own_evidence():
    """M7's safety property: τ_own = 0 ⇒ v_own = ∞ ⇒ b̂² = 0 ⇒ the message passes BIT-IDENTICALLY. This is the
    AMBIG / unstranded regime — where cross-node messages are the only information — so the term must not
    touch it, and a CONTRADICTED claim is not damped there either (there is no evidence to contradict it)."""
    p = np.array([25.0, 3.0, 0.0])
    out = mismatch_deflate(p, np.array([0.0, 9.9, 4.0]), np.array([False, True, True]), np.full(3, np.inf))
    assert out.tolist() == p.tolist()


def test_mismatch_deflate_kills_a_contradicted_claim_where_there_is_evidence():
    """A message asserting a component is ABSENT at a node whose own evidence says otherwise is the b̂² → ∞
    limit ⇒ precision 0 — expressed as a mask so the numerical zero-test ``_EPS`` never sets the answer."""
    out = mismatch_deflate(np.array([50.0, 50.0]), np.zeros(2), np.array([True, False]), np.array([0.1, 0.1]))
    assert out[0] == 0.0 and out[1] == pytest.approx(50.0, rel=1e-12)


def test_mismatch_deflate_charges_the_full_gap_when_the_node_is_composition_certain():
    """v_own = 0 (a structural pure-gDNA anchor: composition CERTAIN) ⇒ the whole G² is charged, so no message
    can talk it off its composition."""
    out = mismatch_deflate(np.array([100.0]), np.array([1.5]), np.zeros(1, bool), np.array([0.0]))
    assert out[0] == pytest.approx(1.0 / 2.25, rel=1e-12)


def test_mismatch_deflate_pin_safety_invariant():
    """THE governing-principle invariant, as an exact inequality: a message out-weighs the destination's own
    belief iff it agrees to within √2·σ_own. From p_eff = 1/max(v_msg,G²−v_own), p_eff > 1/v_own ⟺ G² < 2·v_own
    — independent of the message's own precision, so no amount of source depth can buy a pin."""
    v_own = 0.32
    p_own = 1.0 / v_own
    for scale, expect_wins in ((0.9, True), (1.1, False)):
        gap = scale * np.sqrt(2.0 * v_own)
        # even an ARBITRARILY confident message obeys the threshold
        out = mismatch_deflate(np.array([1e6]), np.array([gap]), np.zeros(1, bool), np.array([v_own]))
        assert bool(out[0] > p_own) is expect_wins


def test_mismatch_deflate_recovers_the_true_bias_squared_by_mc():
    """M7c: DerSimonian–Laird is a method-of-moments estimator of the between-source variance, so for a REAL
    mismatch b̂² → b² — the load-bearing claim, with no tuned constant. (`scripts/debug/message_variance_mc.py`
    M7c is the full sweep; this pins the identity in-repo.)"""
    rng = np.random.default_rng(20260725)
    v_msg, v_own, b, n = 0.04, 0.30, 2.6, 200_000
    gap = b + rng.normal(0.0, np.sqrt(v_msg + v_own), n)
    p_eff = mismatch_deflate(np.full(n, 1.0 / v_msg), gap, np.zeros(n, bool), np.full(n, v_own))
    b2_hat = float(np.mean(1.0 / p_eff)) - v_msg
    assert b2_hat == pytest.approx(b * b, rel=0.02)


def test_mismatch_deflate_is_finite_over_every_degenerate_input():
    """No nan, no ∞, no 1/0 — over zero precisions, zero gaps, ∞ and 0 own variances, contradicted or not."""
    p = np.array([0.0, 1e-12, 25.0, 25.0, 25.0, 0.0])
    gap = np.array([0.0, 3.0, 0.0, -4.0, 1e3, 1e3])
    contra = np.array([True, False, True, False, True, True])
    v_own = np.array([np.inf, 0.0, 0.5, np.inf, 0.0, 0.0])
    out = mismatch_deflate(p, gap, contra, v_own)
    assert np.all(np.isfinite(out)) and np.all(out >= 0.0) and np.all(out <= p + 1e-12)


# ── M10 the peel as a composition (a share), not a subtraction ──


def test_peel_continue_share_is_enrichment_free():
    """M10's defining property: the continuing SHARE is invariant under a common capture factor, because
    capture multiplies the continuing and the splicing channels alike. This is the whole reason the peel
    becomes a scaling — a scaling commutes with the reframe, a subtraction does not."""
    nu, mu = 3.0, 7.0
    base = float(peel_continue_share(nu, mu))
    for e in (1e-3, 0.5, 1.0, 40.0, 1e4):
        assert float(peel_continue_share(nu * e, mu * e)) == pytest.approx(base, rel=1e-14)
    assert base == pytest.approx(0.3, rel=1e-14)


def test_peel_continue_share_structural_limits():
    """No spliced flux ⇒ nothing is peeled (w = 1). No RNA at the seam at all ⇒ nothing to apportion (w = 1),
    the caller's own gates decide. And w never leaves [0, 1] — which is what retires the old peel's
    zero-truncation (a fully-consumed subtraction emitted ρ_ν = 0 at a LIVE precision)."""
    assert float(peel_continue_share(5.0, 0.0)) == 1.0
    assert float(peel_continue_share(0.0, 0.0)) == 1.0
    assert float(peel_continue_share(0.0, 5.0)) == 0.0
    out = peel_continue_share(np.array([1.0, 0.0, 4.0, 0.0]), np.array([1.0, 3.0, 0.0, 0.0]))
    assert np.all((out >= 0.0) & (out <= 1.0))


def test_peel_share_logvar_is_convex_unlike_the_subtraction():
    """M10 vs M3. The share's delta-method weights are w_μ² ≤ 1 (convex — the mirror of M2's graft SUM),
    where the subtraction carried u² ≥ 1 and AMPLIFIED. At the same operating point the share must therefore
    cost strictly less than the difference it replaces."""
    v_nu, v_mu = 0.02, 1.0 / 4000.0
    for w in (0.5, 0.25, 0.1):
        w_mu = 1.0 - w
        share = float(peel_share_logvar(w_mu, v_nu, v_mu))
        assert share == pytest.approx(w_mu * w_mu * (v_nu + v_mu), rel=1e-14)
        assert share <= v_nu + v_mu  # convex: never worse than its own inputs
        u = 1.0 / w  # the M3 difference at the same continuing fraction
        assert share < float(peel_rna_logvar(v_nu, 0.0, v_mu, u))


def test_peel_share_logvar_vanishes_with_no_spliced():
    """w_μ = 0 (nothing splices away) ⇒ the share contributes NO variance: the message passes through intact."""
    assert float(peel_share_logvar(0.0, 0.5, 0.25)) == 0.0


# ── M11: the LEVEL from the node's own mass + an imputed gDNA density ────────────────────────────────────


def test_residual_level_pure_rna_limit_is_the_count():
    """M11 limit 1 — a gDNA claim of ZERO accounts for none of the crossing, so the whole mass is RNA and the
    only uncertainty left is the node's own Poisson count. This is the limit that makes the level a
    MEASUREMENT in a low-gDNA library, which is exactly where the old no-evidence default silenced the RNA
    channel outright. It is exact only because the arithmetic is done on the FRACTION: the count cancels out
    of ``φ = ρ_g E_g/M``, so ``σ_f → 0`` with ``φ`` and the upper truncation never bites."""
    rho, v_log, v_lin = residual_level(
        mass=1.0e4, n_mass=400.0, rho_g=0.0, E_g=1000.0, E_r=900.0, v_g=1.0e-9
    )
    assert float(rho) == pytest.approx(1.0e4 / 900.0, rel=1e-12)
    assert float(v_log) == pytest.approx(1.0 / 400.0, rel=1e-6)
    assert float(v_lin) == pytest.approx(float(rho) ** 2 * float(v_log), rel=1e-12)


def test_residual_level_gdna_explains_everything_gives_a_tight_linear_zero():
    """M11 limit 2 — a confident gDNA claim that over-explains the crossing drives the RNA density to ~0. The
    LOG-variance is then large (log of a near-zero quantity is unbounded), but the LINEAR statement is tight:
    "below a fraction of a percent of my mass". That asymmetry is the whole reason the consumer fuses levels
    in linear space, and it is what reproduces "intronic unspliced fragments are gDNA until proven otherwise"
    as a measurement rather than as a rule."""
    M, E_r = 1.0e4, 900.0
    rho, v_log, v_lin = residual_level(
        mass=M, n_mass=400.0, rho_g=40.0, E_g=1000.0, E_r=E_r, v_g=1.0e-9
    )
    assert float(rho) / (M / E_r) < 1.0e-6  # essentially no RNA
    assert float(v_log) > np.pi**2 / 6.0  # ...and it knows the log is uninformative
    assert np.sqrt(float(v_lin)) / (M / E_r) < 0.01  # ...while the linear claim is tight


def test_residual_level_ignorance_is_bounded_and_declared():
    """M11 limit 3 — the upper truncation. An imputed gDNA claim carrying ~1 nat of log-variance (routine at
    exon→boundary edges under capture) makes σ_f of order 1. A one-sided positive part would return
    ``E ≈ 0.8σ`` — "most of my mass is RNA" — asserted out of pure ignorance at a confident-looking k ≈ 2.
    Bounded above, the same ignorance degrades to its correct limit: the uniform posterior, ``f_R = ½`` at
    ``k = 3``, wide enough for any real evidence in the fuse to out-weigh it."""
    M, E_r = 1.0e4, 900.0
    rho, v_log, _ = residual_level(
        mass=M, n_mass=1.0e6, rho_g=10.0, E_g=1000.0, E_r=E_r, v_g=100.0
    )
    assert float(rho) / (M / E_r) == pytest.approx(0.5, rel=5e-3)  # an asymptote, approached from below
    assert float(v_log) == pytest.approx(float(polygamma(1, 3.0)), rel=5e-3)


def test_residual_level_needs_the_gdna_claim_to_be_SUPPLIED():
    """M11 — "supplied" is a statement about PRECISION, never about the density's value (the same test the
    λ-emission gate makes). A gDNA claim at infinite log-variance is not a claim, so there is no level at all
    — and, the standing trap, it must not come back as a nan from ``0·inf``."""
    for kwargs in (
        dict(v_g=np.inf, rho_g=1.0),  # no precision on the claim
        dict(v_g=0.01, rho_g=1.0, n_mass=0.0),  # no count
    ):
        base = dict(mass=1.0e4, n_mass=400.0, E_g=1000.0, E_r=900.0)
        rho, v_log, v_lin = residual_level(**{**base, **kwargs})
        assert float(rho) == 0.0
        assert not np.isfinite(float(v_log))
        assert not np.isfinite(float(v_lin))
        assert not np.isnan(float(rho))


def test_residual_level_is_monotone_in_the_gdna_claim():
    """M11 — more imputed gDNA leaves less room for RNA, always. A monotonicity the truncation must not break
    at either bound (the two closed-form tail branches are the places it could)."""
    prev = np.inf
    for rho_g in (0.0, 1.0, 3.0, 5.0, 7.0, 9.0, 9.9, 10.0, 12.0, 40.0, 400.0, 1000.0):
        rho, _, _ = residual_level(
            mass=1.0e4, n_mass=400.0, rho_g=rho_g, E_g=1000.0, E_r=900.0, v_g=0.01
        )
        assert float(rho) <= prev + 1e-12
        prev = float(rho)


# ── M12: the conservation rescale — restore Σ_c ρ_c·E_c = M by moving what is least known ────────────────

_CR_RHO = np.array([[2.0, 3.0, 1.0], [0.4, 12.0, 0.0], [55.0, 0.02, 0.03]])
_CR_EFF = np.array([[2100.0, 1900.0, 1900.0], [110.0, 180.0, 180.0], [9000.0, 8800.0, 8800.0]])
_CR_M = (_CR_RHO * _CR_EFF).sum(axis=1) * np.array([0.7, 1.4, 1.0])  # 0.7x under, 1.4x over, balanced


def _cr_total(rho):
    return (rho * _CR_EFF).sum(axis=1)


def test_conservation_rescale_common_limit_is_the_shipped_pin():
    """M12 limit 1 — when the ONLY error is the shared frame (``w = 0``, i.e. every component's variance is
    the common one) the direction is constant across components, so every component moves by the same factor
    ``k = M/S``. That is exactly what ``_pin_v`` does: the shipped operator is the zero-independent-variance
    limit of this one, not a different operator. Nothing is being replaced."""
    v = np.full_like(_CR_RHO, 0.05)
    out = conservation_rescale(_CR_M, _CR_RHO, _CR_EFF, v, np.full(3, 0.05), _CR_RHO)
    k = _CR_M / _cr_total(_CR_RHO)
    assert out == pytest.approx(_CR_RHO * k[:, None], rel=1e-12)


def test_conservation_rescale_certain_component_does_not_move():
    """M12 limit 2 — a component with NO independent variance is a MEASUREMENT, not an imputation, so it must
    not give any ground: the others absorb the entire residual. This is the limit that matters, because the
    grafted spliced count makes the RNA arm a measurement exactly where pass-0 is confidently wrong."""
    v = np.array([[0.5, 0.0, 0.5], [0.3, 0.0, 0.3], [0.9, 0.0, 0.9]])
    out = conservation_rescale(_CR_M, _CR_RHO, _CR_EFF, v, np.zeros(3), _CR_RHO)
    assert out[:, 1] == pytest.approx(_CR_RHO[:, 1], rel=1e-12)  # untouched
    assert _cr_total(out) == pytest.approx(_CR_M, rel=1e-10)  # ...and the identity still holds


def test_conservation_rescale_restores_the_identity_exactly():
    """M12 — the constraint is an IDENTITY, so it is imposed exactly rather than to first order: the
    magnitude is re-solved along the direction the error model chose. Checked over six orders of magnitude of
    density and of violation, because the linearised form would drift badly at the extremes."""
    rng = np.random.default_rng(20260726)
    for scale in (1e-3, 1.0, 1e3):
        rho = rng.lognormal(0.0, 2.0, size=(64, 3)) * scale
        eff = rng.uniform(50.0, 9000.0, size=(64, 3))
        M = (rho * eff).sum(axis=1) * rng.uniform(0.1, 8.0, size=64)
        v = rng.lognormal(-2.0, 1.5, size=(64, 3))
        out = conservation_rescale(M, rho, eff, v, rng.uniform(0.0, 0.1, 64), rho)
        assert (out * eff).sum(axis=1) == pytest.approx(M, rel=1e-10)


def test_conservation_rescale_keeps_a_partial_claim_partial():
    """M12 — "supplied" is a statement about PRECISION, never about the density's value. An unsupplied
    component (``var = +inf``) contributes the NODE'S OWN density to the mass budget and does not move, so a
    message carrying gDNA only still delivers ``f_g < 1``. Rescaling all components blindly instead regresses
    capture-OFF 3.6x (`scratchpad/derive_2_relay_pin.py`) — this is the load-bearing semantic."""
    rho = np.array([[3.0, 0.0, 0.0]])
    own = np.array([[1.0, 4.0, 2.0]])
    eff = np.array([[2000.0, 1800.0, 1800.0]])
    M = np.array([20000.0])
    v = np.array([[0.2, np.inf, np.inf]])
    out = conservation_rescale(M, rho, eff, v, np.zeros(1), own)
    assert out[0, 1] == 0.0 and out[0, 2] == 0.0  # not supplied ⇒ not moved, still absent
    # the budget used the node's OWN density for the two unsupplied arms, so the gDNA arm absorbs only the
    # residual against them — it is NOT renormalised onto the whole mass:
    assert out[0, 0] * eff[0, 0] == pytest.approx(M[0] - (own[0, 1] + own[0, 2]) * eff[0, 1], rel=1e-10)
    assert not np.isnan(out).any()


def test_conservation_rescale_is_idempotent_when_it_already_balances():
    """M12 — a claim that already accounts for exactly the observed mass is left alone (``μ = 0``), and
    applying the operator twice changes nothing further."""
    rho = np.array([[2.0, 3.0, 1.0]])
    eff = np.array([[2100.0, 1900.0, 1900.0]])
    M = _f_total = (rho * eff).sum(axis=1)
    v = np.array([[0.4, 0.01, 0.4]])
    once = conservation_rescale(M, rho, eff, v, np.zeros(1), rho)
    assert once == pytest.approx(rho, rel=1e-10)
    assert conservation_rescale(M, once, eff, v, np.zeros(1), rho) == pytest.approx(once, rel=1e-12)


def test_conservation_rescale_is_monotone_in_the_violation():
    """M12 — a bigger violation produces a bigger correction, in the same direction, for every component.
    The two closed-form branches (bracket expansion and the dead/fallback corner) are where this could break."""
    rho = np.array([[2.0, 3.0, 1.0]])
    eff = np.array([[2100.0, 1900.0, 1900.0]])
    v = np.array([[0.4, 0.02, 0.4]])
    S = (rho * eff).sum()
    prev = None
    for f in (0.25, 0.5, 0.9, 1.0, 1.2, 2.0, 6.0):
        out = conservation_rescale(np.array([S * f]), rho, eff, v, np.zeros(1), rho)[0]
        if prev is not None:
            assert np.all(out >= prev - 1e-12)
        prev = out


def test_conservation_rescale_degenerate_inputs_are_safe():
    """M12 — no nan, no inf, and no silent identity violation at any degenerate input: zero mass, a zero
    claim, a single supplied component, and a model that admits no uncertainty at all (in which case there is
    nothing to apportion and the common factor is the honest fallback)."""
    eff = np.array([[2100.0, 1900.0, 1900.0]])
    cases = [
        (np.array([0.0]), np.array([[2.0, 3.0, 1.0]]), np.array([[0.1, 0.1, 0.1]])),  # no mass
        (np.array([5000.0]), np.zeros((1, 3)), np.array([[0.1, 0.1, 0.1]])),  # no claim
        (np.array([5000.0]), np.array([[2.0, 0.0, 0.0]]), np.array([[0.1, np.inf, np.inf]])),  # one arm
        (np.array([5000.0]), np.array([[2.0, 3.0, 1.0]]), np.zeros((1, 3))),  # no uncertainty anywhere
    ]
    for M, rho, v in cases:
        out = conservation_rescale(M, rho, eff, v, np.zeros(1), rho)
        assert np.isfinite(out).all()
    # the no-uncertainty corner must still balance, via the common-factor fallback
    M, rho, v = cases[3]
    out = conservation_rescale(M, rho, eff, v, np.zeros(1), rho)
    assert (out * eff).sum() == pytest.approx(M[0], rel=1e-10)


# ── P1d: graft_premise_logvar — the two-seam method-of-moments premise variance ──────────────────────
# ⚠ The SOLVER uses only the pooled (second) return value; the per-edge array is a diagnostic. A variance
# estimated from a single pair is a χ²₁ (CV = √2), so per-edge is mostly noise — see the docstring.


def test_graft_premise_logvar_agreeing_seams_charge_nothing():
    """Two seams that agree carry no detectable premise error — the truncation is the method's own."""
    per, pooled = graft_premise_logvar(
        np.array([3.0, 7.0]), np.array([3.0, 7.0]), np.array([0.01, 0.01]), np.array([0.01, 0.01])
    )
    assert list(per) == [0.0, 0.0]
    assert pooled == 0.0


def test_graft_premise_logvar_is_the_mom_second_moment_halved():
    """``max(0, d² − noise)/2`` exactly, with no coefficient to choose."""
    fa, fb = np.array([np.e**2]), np.array([1.0])  # d = 2 exactly
    va = vb = np.array([0.25])
    per, _ = graft_premise_logvar(fa, fb, va, vb)
    assert float(per[0]) == pytest.approx((4.0 - 0.5) / 2.0, rel=1e-12)


def test_graft_premise_logvar_subtracts_noise_and_is_direction_free():
    """A gap fully explained by the seams' own noise leaves nothing; a↔b swap cannot change the answer."""
    fa, fb = np.array([np.e]), np.array([1.0])  # d = 1
    per, _ = graft_premise_logvar(fa, fb, np.array([0.6]), np.array([0.6]))  # noise 1.2 > d² = 1
    assert float(per[0]) == 0.0
    p1, _ = graft_premise_logvar(fa, fb, np.array([0.1]), np.array([0.2]))
    p2, _ = graft_premise_logvar(fb, fa, np.array([0.2]), np.array([0.1]))
    assert float(p1[0]) == pytest.approx(float(p2[0]), rel=1e-14)


def test_graft_premise_logvar_needs_two_live_seams():
    """One live seam ⇒ no second study ⇒ per-edge 0, and the pooled fit ignores that edge entirely."""
    per, pooled = graft_premise_logvar(
        np.array([5.0, np.e**2]), np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0])
    )
    assert float(per[0]) == 0.0
    assert float(per[1]) == pytest.approx(2.0, rel=1e-12)
    assert pooled == pytest.approx(2.0, rel=1e-12)  # fitted on the ONE edge that has a pair


def test_graft_premise_logvar_pooled_is_the_population_second_moment():
    """The pooled fit is the same estimator over the population — NOT the mean of the per-edge values,
    which would truncate each edge separately and bias the fit upward."""
    fa = np.exp(np.array([2.0, 0.0, 0.0]))  # d² = 4, 0, 0 against a per-seam noise of 0.5+0.5
    fb = np.ones(3)
    v = np.full(3, 0.5)
    per, pooled = graft_premise_logvar(fa, fb, v, v)
    d2 = np.array([4.0, 0.0, 0.0])
    assert list(per) == pytest.approx([1.5, 0.0, 0.0], rel=1e-12)
    assert pooled == pytest.approx(max(0.0, d2.mean() - 1.0) / 2.0, rel=1e-12)
    assert pooled < float(np.mean(per))  # the per-edge truncation really does bias upward


def test_graft_premise_logvar_infinite_noise_is_ignored_not_propagated():
    """A seam with no count (var = inf) must not nan the estimate — it contributes no subtraction."""
    per, pooled = graft_premise_logvar(
        np.array([np.e**2]), np.array([1.0]), np.array([np.inf]), np.array([0.0])
    )
    assert np.isfinite(per).all() and np.isfinite(pooled)
    assert float(per[0]) == pytest.approx(2.0, rel=1e-12)


def test_graft_premise_logvar_never_negative_and_finite():
    rng = np.random.default_rng(11)
    fa, fb = rng.lognormal(0.0, 1.5, 400), rng.lognormal(0.0, 1.5, 400)
    va, vb = rng.gamma(2.0, 0.3, 400), rng.gamma(2.0, 0.3, 400)
    per, pooled = graft_premise_logvar(fa, fb, va, vb)
    assert np.isfinite(per).all() and (per >= 0.0).all()
    assert np.isfinite(pooled) and pooled >= 0.0


def test_graft_premise_logvar_pooled_is_the_load_bearing_return():
    """The solver applies the POOLED scalar to every edge. It must therefore exist and be finite even when
    most edges have no pair at all — the terminal-exon / exon↔exon-boundary case the owner flagged."""
    fa = np.array([np.e, 0.0, 0.0, 0.0, 0.0])  # one pair, four seams with no partner
    fb = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
    v = np.zeros(5)
    per, pooled = graft_premise_logvar(fa, fb, v, v)
    assert list(per[1:]) == [0.0, 0.0, 0.0, 0.0]
    assert pooled == pytest.approx(0.5, rel=1e-12)  # d² = 1, no noise, halved
    assert np.isfinite(pooled)
