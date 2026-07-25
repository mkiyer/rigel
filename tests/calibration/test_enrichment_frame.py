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

from rigel.calibration.enrichment_frame import (
    boundary_unspliced_from_k,
    composition_logvar,
    density_mode_logfrac,
    enrichment_ratio,
    f_g_from_k,
    gdna_fallback_admissible,
    graft_frame_logvar,
    graft_rna_logvar,
    k_from_belief,
    message_precision,
    mismatch_deflate,
    mismatch_gap,
    peel_rna_logvar,
    reframe_density,
    total_density,
    transfer_logvar,
    transport_seed_logvar,
)

# Real-ish effective lengths from the FL models (gDNA ~N(300,60), RNA ~N(200,50)):
# contained regions grow with L; boundary crossings saturate at fl_mean.
EG_REG, ER_REG = 2701.0, 2801.0  # a long contained region (L = 3000)
EG_BND, ER_BND = 300.0, 200.0  # a boundary crossing (E_g = fl_mean_gdna, E_r = fl_mean_rna)


# ---------------------------------------------------------------------------
# §1 / §4 — the composed total-density identity, mutually exact.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 100.0])
@pytest.mark.parametrize("mass", [1.0, 37.0, 1000.0])
def test_total_density_equals_the_composed_kplus1_identity(k, mass):
    """``total_density(M, f_g_from_k(k), E_g, E_r) == M·(k+1)/(k·E_g + E_r)`` — the identity §4 is built on.

    ``f_g`` re-formed from ``k`` in this node's frame, fed back through ``total_density``, must reproduce the
    closed form exactly (to floating point). This is what lets the solver resolve a boundary's total density
    from a transported ``k`` without circularity."""
    fg = float(f_g_from_k(k, EG_BND, ER_BND))
    rho = float(total_density(mass, fg, EG_BND, ER_BND))
    closed = mass * (k + 1.0) / (k * EG_BND + ER_BND)
    assert rho == pytest.approx(closed, rel=1e-14)


@pytest.mark.parametrize("f_g", [0.02, 0.25, 0.5, 0.75, 0.98])
def test_k_and_f_g_are_inverses(f_g):
    """``f_g_from_k(k_from_belief(f_g)) == f_g`` in the SAME frame — transport ``k``, re-form ``f_g``,
    recover the original when the frame does not change."""
    k = float(k_from_belief(f_g, EG_BND, ER_BND))
    fg_back = float(f_g_from_k(k, EG_BND, ER_BND))
    assert fg_back == pytest.approx(f_g, rel=1e-13)


@pytest.mark.parametrize("k", [0.05, 0.5, 2.0, 20.0])
@pytest.mark.parametrize("mass", [5.0, 250.0])
def test_boundary_unspliced_satisfies_the_mass_identity_exactly(k, mass):
    """``boundary_unspliced_from_k`` must satisfy ``ρ_g·E_g + ρ_R·E_r == M`` and ``ρ_g/ρ_R == k`` exactly —
    the observation constraint that closes the step-wise solve (§4 step 1)."""
    rho_g, rho_R = boundary_unspliced_from_k(k, mass, EG_BND, ER_BND)
    rho_g, rho_R = float(rho_g), float(rho_R)
    assert rho_g * EG_BND + rho_R * ER_BND == pytest.approx(mass, rel=1e-13)
    assert rho_g / rho_R == pytest.approx(k, rel=1e-13)


def test_boundary_unspliced_infinite_k_is_pure_gdna_at_the_gdna_frame():
    """The 100 %-gDNA fallback (§4 step 2): ``k = +inf`` ⇒ ``ρ_R = 0`` and ``ρ_g = M/E_g`` — the pure-gDNA
    density in the boundary's own frame, resolving the FRAME without a division-by-zero or a nan."""
    rho_g, rho_R = boundary_unspliced_from_k(np.inf, 300.0, EG_BND, ER_BND)
    assert float(rho_R) == 0.0
    assert float(rho_g) == pytest.approx(300.0 / EG_BND, rel=1e-14)


# ---------------------------------------------------------------------------
# §5b — the r₂ (k+1) cancellation, as an exact algebraic equality.
# ---------------------------------------------------------------------------


def test_r2_intron_to_boundary_kplus1_cancels_over_random_draws():
    """r₂ built from two ``total_density`` calls must equal the closed form where ``(k+1)`` has cancelled
    (``junction_enrichment_scaling.md`` §5b): ``r₂ = (M_B/M_I)·(k·E_gI + E_rI)/(k·E_gB + E_rB)``.

    The boundary borrows the intron's ``k`` (§3.1), so both totals are formed at the SAME ``k`` — the
    intron contained, the boundary crossing. Verified to machine precision over random ``k``, ``M_B``, ``M_I``."""
    rng = np.random.default_rng(0)
    k = rng.lognormal(0.0, 2.0, 20000)
    M_B, M_I = rng.gamma(3.0, 50.0, 20000), rng.gamma(3.0, 50.0, 20000)

    fg_I = f_g_from_k(k, EG_REG, ER_REG)  # intron in its contained frame
    fg_B = f_g_from_k(k, EG_BND, ER_BND)  # boundary in its crossing frame (same k)
    rho_tot_I = total_density(M_I, fg_I, EG_REG, ER_REG)
    rho_tot_B = total_density(M_B, fg_B, EG_BND, ER_BND)

    r2, _ = enrichment_ratio(rho_tot_B, rho_tot_I)
    closed = (M_B / M_I) * (k * EG_REG + ER_REG) / (k * EG_BND + ER_BND)
    assert np.max(np.abs(r2 / closed - 1.0)) < 1e-13


def test_enrichment_ratio_orientation_and_variance_are_additive():
    """``r = ρ_tot(dst)/ρ_tot(src)`` (single canonical orientation) and ``Var(log r) = Var_dst + Var_src``
    (ratio of independent quantities). A transported density is then ``ρ_src·r``."""
    r, vr = enrichment_ratio(6.0, 2.0, var_log_dst=0.1, var_log_src=0.3)
    assert float(r) == pytest.approx(3.0, rel=1e-14)
    assert float(vr) == pytest.approx(0.4, rel=1e-14)
    # transport recovers the destination total exactly
    assert 2.0 * float(r) == pytest.approx(6.0, rel=1e-14)


# ---------------------------------------------------------------------------
# §2 — the bounding lemma + the f_g=1 degeneracy (the two estimators coincide there).
# ---------------------------------------------------------------------------


def test_total_density_equals_M_over_Eg_at_pure_gdna():
    """At ``f_g = 1`` the blend degenerates to ``M/E_g`` — this is WHY family A (intergenic/seam, structurally
    ``f_g ≡ 1``) is bit-identical between the gDNA-frame and blended estimators (E0 finding 2)."""
    for M, Eg, Er in ((123.0, EG_REG, ER_REG), (7.0, EG_BND, ER_BND)):
        assert float(total_density(M, 1.0, Eg, Er)) == pytest.approx(M / Eg, rel=1e-14)
        assert float(total_density(M, 0.0, Eg, Er)) == pytest.approx(M / Er, rel=1e-14)


@pytest.mark.parametrize("Eg,Er", [(EG_REG, ER_REG), (EG_BND, ER_BND), (24.4, 101.4)])
def test_bounding_lemma_holds_for_every_composition(Eg, Er):
    """``M/max(E_g,E_r) ≤ ρ_tot ≤ M/min(E_g,E_r)`` for ANY ``f_g ∈ [0,1]`` (§2) — the convex-combination
    property that makes a totally wrong composition cost only the effective-length ratio."""
    M = 500.0
    lo, hi = M / max(Eg, Er), M / min(Eg, Er)
    for fg in np.linspace(0.0, 1.0, 51):
        rho = float(total_density(M, fg, Eg, Er))
        assert lo - 1e-9 <= rho <= hi + 1e-9


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


def test_gdna_fallback_admissible_only_on_a_clean_junction():
    """The fallback is admissible ONLY on a genuine junction with none of the RNA-rich exceptions; each §5
    condition refuses it independently."""
    # a clean exon↔intron junction: admissible
    assert bool(
        gdna_fallback_admissible(
            is_ambig=False, mature_crosses=False, has_junction=True, retained_intron=False
        )
    )
    # each exception, alone, refuses
    for kw in (
        dict(is_ambig=True, mature_crosses=False, has_junction=True, retained_intron=False),
        dict(is_ambig=False, mature_crosses=True, has_junction=True, retained_intron=False),
        dict(is_ambig=False, mature_crosses=False, has_junction=False, retained_intron=False),
        dict(is_ambig=False, mature_crosses=False, has_junction=True, retained_intron=True),
    ):
        assert not bool(gdna_fallback_admissible(**kw))


# ---------------------------------------------------------------------------
# The UNIFIED SOLVER theorem (unified_solver_design.md §2): reframe + density-mode ÷ M_dst subsumes the shift,
# is enrichment-invariant, and handles a PARTIAL (set-mismatched) source correctly where the shift cannot.
# ---------------------------------------------------------------------------


def _shift_logfrac(rho_c_src, E_c_src_frame, imputed_masses):
    """The OLD shift factor on log f_c: log( M_c / Σ M_c' ), M_c imputed with the dst's E_c. Reference to
    compare the unified density mode against."""
    M_c = rho_c_src * E_c_src_frame
    return np.log(M_c / np.sum(imputed_masses))


@pytest.mark.parametrize("r", [0.01, 1.0, 7.5, 1000.0])
def test_unified_equals_shift_for_a_complete_matched_message(r):
    """§2: for a COMPLETE matched source (the dst's observed mass == the reframed imputed total), the density
    mode ÷ M_dst reproduces the shift EXACTLY, for ANY enrichment ratio r — the reframe cancels."""
    # source component densities (its own frame), dst per-component eff-lengths (crossing frame)
    rho = {"g": 0.6, "p": 0.25, "n": 0.15}
    E = {"g": 300.0, "p": 200.0, "n": 200.0}
    # reframe every component by the SAME r, impute at the dst -> the dst's observed mass is that total
    imp = {c: reframe_density(rho[c], r, 1.0) * E[c] for c in rho}  # rho_tot_dst=r, rho_tot_src=1 -> ratio r
    M_dst = float(sum(imp.values()))  # matched-complete: observed == imputed
    for c in rho:
        rho_reframed = float(reframe_density(rho[c], r, 1.0))
        uni = float(density_mode_logfrac(rho_reframed, E[c], M_dst))
        shift = float(_shift_logfrac(rho[c], E[c], np.array(list({k: rho[k] * E[k] for k in rho}.values()))))
        assert uni == pytest.approx(shift, rel=1e-12), f"component {c}, r={r}"


def test_unified_is_enrichment_invariant():
    """§2: the resulting fraction does not depend on r (the enrichment cancels via ÷ M_dst)."""
    rho = {"g": 0.6, "p": 0.4}
    E = {"g": 300.0, "p": 200.0}
    fracs = []
    for r in (0.02, 1.0, 500.0):
        imp = {c: reframe_density(rho[c], r, 1.0) * E[c] for c in rho}
        M_dst = float(sum(imp.values()))
        f_g = np.exp(float(density_mode_logfrac(reframe_density(rho["g"], r, 1.0), E["g"], M_dst)))
        fracs.append(f_g)
    assert max(fracs) - min(fracs) < 1e-12, "f_g must be invariant to the enrichment ratio"


def test_unified_handles_a_partial_source_where_the_shift_fails():
    """§2: a PARTIAL source (gDNA only — a seam) must give f_g < 1 via ÷ M_dst (the dst's own mass supplies the
    RNA it lacks). The shift's imputed-total normalizer would assert f_g = 1 and be wrong."""
    rho_g_src, E_g = 0.5, 300.0
    # the dst (an exon) actually carries gDNA + RNA, so its observed mass exceeds the gDNA-only imputed mass
    M_dst_gdna = rho_g_src * E_g  # what a gDNA-only impute would predict
    M_dst_observed = 4.0 * M_dst_gdna  # the exon is mostly RNA -> observed is 4x the gDNA mass
    f_g_unified = np.exp(float(density_mode_logfrac(rho_g_src, E_g, M_dst_observed)))
    f_g_shift = np.exp(float(_shift_logfrac(rho_g_src, E_g, np.array([M_dst_gdna]))))  # imputed total = gDNA only
    assert f_g_unified == pytest.approx(0.25, rel=1e-9), "div by observed mass gives the true gDNA fraction"
    assert f_g_shift == pytest.approx(1.0, rel=1e-9), "the shift wrongly asserts pure gDNA on a partial source"


def test_reframe_is_a_pure_scale_and_pass_through_when_source_frame_is_degenerate():
    """The reframe multiplies by ρ_tot(dst)/ρ_tot(src); a zero-mass (frameless) source gives ratio 1
    (pass-through), never a division blow-up."""
    assert float(reframe_density(2.0, 10.0, 5.0)) == pytest.approx(4.0, rel=1e-12)
    assert float(reframe_density(2.0, 10.0, 0.0)) == pytest.approx(2.0, rel=1e-12)  # degenerate src -> pass


def test_enrichment_frame_functions_are_vectorized():
    """Every primitive is an array function (the solver drives them over the whole chain) — shapes broadcast
    and elementwise results match the scalar path."""
    fg = np.array([0.2, 0.8, 1.0])
    eg = np.array([EG_BND, EG_REG, EG_BND])
    er = np.array([ER_BND, ER_REG, ER_BND])
    m = np.array([10.0, 20.0, 30.0])
    rho = total_density(m, fg, eg, er)
    assert rho.shape == (3,)
    for i in range(3):
        assert rho[i] == pytest.approx(float(total_density(m[i], fg[i], eg[i], er[i])), rel=1e-14)


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


def test_transport_seed_is_composition_plus_count():
    """M1: Var(log ρ_c) = Var(log f_c) + 1/n. Measured spliced (var_log_f=0) ⇒ pure 1/n; n=0 ⇒ ∞ (no message)."""
    assert float(transport_seed_logvar(0.02, 500.0)) == pytest.approx(0.02 + 1.0 / 500.0, rel=1e-14)
    assert float(transport_seed_logvar(0.0, 800.0)) == pytest.approx(1.0 / 800.0, rel=1e-14)  # spliced
    assert np.isinf(float(transport_seed_logvar(0.02, 0.0)))  # no count ⇒ no message
    assert np.isinf(float(transport_seed_logvar(np.inf, 500.0)))  # τ_λ=0 ⇒ no composition ⇒ ∞


def test_transport_seed_matches_mc():
    """The seed equals the empirical Var(log(f_g·M/E_g)) — composition (Beta) ⊕ Poisson count (lognormal)."""
    f_g, var_fg, n, E_g = 0.35, 0.003, 900.0, 120.0
    fg = _beta_mc(f_g, var_fg, 400_000)
    M = _lognormal_mc(700.0, 1.0 / n, 400_000)
    emp = float(np.var(np.log(fg * M / E_g)))
    pred = float(transport_seed_logvar(var_fg / f_g**2, n))  # Var(log f_g) = var_fg/f_g² (LOG Jacobian)
    assert pred == pytest.approx(emp, rel=0.06)


# ── M2 graft (SUM, share-weighted) ──


def test_graft_rna_is_share_weighted_sum():
    """M2: Var(log ρ_R) = w_ν²·v_ν + w_μ²·v_μ (convex weights). Degenerate w_μ=0 ⇒ v_ν; w_μ=1 ⇒ v_μ."""
    assert float(graft_rna_logvar(0.01, 0.002, 0.7, 0.3)) == pytest.approx(
        0.49 * 0.01 + 0.09 * 0.002, rel=1e-14
    )
    assert float(graft_rna_logvar(0.01, 0.002, 1.0, 0.0)) == pytest.approx(0.01, rel=1e-14)
    assert float(graft_rna_logvar(0.01, 0.002, 0.0, 1.0)) == pytest.approx(0.002, rel=1e-14)


def test_graft_rna_matches_mc():
    """The share-weighted SUM equals the empirical Var(log(ρ_ν+ρ_μ)) — the item-E rule falls out of the delta
    method, and beats the naive (unweighted) v_ν+v_μ at intermediate w_μ."""
    f_g, var_fg, n, n_s = 0.40, 0.004, 800.0, 600.0
    E_r, E_spl, M_b, S_b = 200.0, 100.0, 900.0, 1500.0
    fg = _beta_mc(f_g, var_fg, 400_000)
    Mb = _lognormal_mc(M_b, 1.0 / n, 400_000)
    Sb = _lognormal_mc(S_b, 1.0 / n_s, 400_000)
    rho_nu = (1.0 - fg) * Mb / E_r
    rho_mu = Sb / E_spl
    emp = float(np.var(np.log(rho_nu + rho_mu)))
    rn0, rm0 = (1.0 - f_g) * M_b / E_r, S_b / E_spl
    rR = rn0 + rm0
    v_nu = transport_seed_logvar(var_fg / (1.0 - f_g) ** 2, n)  # Var(log f_R)=var_fg/(1−f_g)²
    v_mu = transport_seed_logvar(0.0, n_s)
    pred = float(graft_rna_logvar(v_nu, v_mu, rn0 / rR, rm0 / rR))
    assert pred == pytest.approx(emp, rel=0.08)


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


def test_message_precision_is_reciprocal_no_jacobian():
    """M4: p_c = 1/Var(log ρ_c) — NO destination Jacobian, NO 1/n_dst. ∞/0/negative variance ⇒ 0 (no message)."""
    v = np.array([0.02, np.inf, 0.0, -1.0, 0.25])
    p = message_precision(v)
    assert p[0] == pytest.approx(1.0 / 0.02, rel=1e-14)
    assert p[1] == 0.0 and p[2] == 0.0 and p[3] == 0.0  # no message, never a nan/∞
    assert p[4] == pytest.approx(4.0, rel=1e-14)
    assert not np.any(np.isnan(p)) and np.all(np.isfinite(p))


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
