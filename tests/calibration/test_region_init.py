"""Unit tests for the pass-0 per-region INITIALIZATION (`calibration.region_init`).

One test per information source of — MEASURED (Poisson
precision), INTRON FACTORY, STRAND DECONVOLUTION, UNSOLVED default (100% gDNA, ZERO precision) — plus the
pure precision arithmetic (`own_composition_logvar`, `own_precision`). These pin the self-solve that seeds
the unified pass-0 relay.
"""

from __future__ import annotations


import numpy as np
import pytest

from rigel.calibration.region_chain import REGION
from rigel.calibration.messages.variance import composition_logvar, count_logvar
from rigel.calibration.region_geometry import g1_locked, init_beliefs
from _synthetic import make_chain_parts
from rigel.calibration.region_init import (
    build_region_init,
    own_composition_logvar,
    own_precision,
    strand_evidence,
)
from rigel.calibration.density_deconv import (
    GdnaBackground,
    density_factor_precision,
    density_lambda_factor,
)
from rigel.calibration.simplex_logodds import _logodds_grid
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_EXON_NEG,
    TS_NONE,
)


# ── helpers ────────────────────────────────────────────────────────────────────────────────────────────────


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _scenario(kappa=0.9):
    """One ref, 3 regions: intergenic (pure gDNA) | exon+ (single-strand) | AMBIG.

    ⭐ The chain is ``N E N E N`` — 5 slots, the regions at 0/2/4 — and the two boundaries carry no counts, so
    this exercises the region axis alone. The predecessor built ``B R B R B R B`` with four boundary slots
    including two data-free terminals; those do not exist.
    """
    parts = make_chain_parts(
        [TS_NONE, BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG],
        region_size_bp=[1500.0, 900.0, 1200.0],
        # intergenic: symmetric gDNA; exon+: sense-tilted RNA (few antisense); AMBIG: symmetric.
        region_pos=[120.0, 200.0, 90.0],
        region_neg=[110.0, 9.0, 92.0],
        gdna_fl=_delta_pmf(200),
        rna_fl=_delta_pmf(100),
    )
    belief = init_beliefs(
        parts.chain,
        parts.geometry,
        parts.statics,
        rna_sense_frac=kappa,
        gdna_strand_overdispersion=0.2,
        rna_strand_overdispersion=0.1,
        n_grid=60,
        n_grid_ss=256,
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _init(kappa=0.9, n_gdna_obs=230.0):
    chain, statics, geometry, belief, _ = _scenario(kappa)
    ni = build_region_init(
        chain,
        statics,
        geometry,
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        n_tilt=None,
        n_grid_ss=256,
        belief=belief,
    )
    # ⭐ ONE count per slot: it is both the density numerator and the Poisson n.
    return ni, np.asarray(geometry.unspliced_count, float).sum(axis=1)


# ── the pure precision arithmetic ────────────────────────────────────────────────────────────────────────


def test_own_composition_logvar_three_states():
    """Locked → 0 (certain); τ>0 → the Jacobian ((1−f_g)²/τ, f_g²/τ); τ=0 → ∞ (unseen), never a nan."""
    fg = np.array([0.2, 0.2, 0.2])
    tau = np.array(
        [0.0, 4.0, 0.0]
    )  # region0 locked (τ irrelevant); region1 unlocked τ>0; region2 unlocked τ=0
    lock = np.array([True, False, False])
    v_fg, v_fr = own_composition_logvar(fg, tau, lock)
    assert v_fg[0] == 0.0 and v_fr[0] == 0.0  # locked ⇒ certain
    assert np.isclose(v_fg[1], (1 - 0.2) ** 2 / 4.0) and np.isclose(v_fr[1], 0.2**2 / 4.0)
    assert np.isinf(v_fg[2]) and np.isinf(v_fr[2])  # no evidence ⇒ unseen
    assert not np.any(np.isnan(v_fg)) and not np.any(np.isnan(v_fr))


def test_own_precision_monotone_and_zeros():
    """``p = 1/(v + trigamma(n+½))``: monotone in the count; 0 only at ``v=∞`` or a dead component.

    ⭐⭐ **RE-POINTED 2026-08-05 — the line this used to carry, ``assert p[2] == 0.0  # no count``, WAS
    THE DEFECT** (TRAPS: a-zero-count-is-a-measurement). It asserted that an object with zero counts must emit nothing, which
    is true of ``1/n`` and false of the world: a zero count over a known opportunity is a measurement,
    and at a structurally pure-gDNA object it is the strongest one in the library. Measured cost of the
    old behaviour: all 1,298 intergenic regions silent at ``g00``, and 34–38 % phantom gDNA on a library
    with none. ⛔ The two zeros that REMAIN are the real ones — ignorance (``v = ∞``) and impossibility
    (a dead component) — and they are asserted below precisely so this does not become "everything now
    speaks"."""
    n = np.array([10.0, 400.0, 0.0, 50.0])
    v = np.array([0.1, 0.1, 0.1, np.inf])
    live = np.array([True, True, True, True])
    p = own_precision(n, v, live)
    assert p[1] > p[0] > 0.0  # more count ⇒ sharper
    assert 0.0 < p[2] < p[0]  # ⭐ zero counts SPEAK, and speak most quietly of the three
    assert p[3] == 0.0  # no evidence (v=∞), and no nan
    assert (
        own_precision(n, np.full(4, 0.1), np.array([False, True, True, True]))[0] == 0.0
    )  # dead component
    assert not np.any(np.isnan(p))


# ── source 3: strand deconvolution evidence ────────────────────────────────────────────────────────────────


def test_strand_evidence_deadband_kills_unstranded():
    """The derived deadband makes I_strand IDENTICALLY 0 on unstranded data (κ=½), >0 stranded; a gDNA-free
    library (N_gdna=0 ⇒ σ²_d→∞) gates it to 0 even when stranded (migrated from `_compile_strand_evidence`)."""
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    reg = np.array([True, True])
    unl = np.array([False, False])
    base = dict(od_g=0.03, od_r=0.03, n_gdna_obs=1e4, n_rna_obs=1e4, is_region=reg, locked=unl)
    tau_unstr, _ = strand_evidence(u, u, fg, kappa=0.5, **base)
    tau_str, _ = strand_evidence(u, u, fg, kappa=0.99, **base)
    assert np.all(tau_unstr == 0.0)
    assert np.all(tau_str > 0.0)
    tau_nog, _ = strand_evidence(
        u,
        u,
        fg,
        kappa=0.99,
        od_g=0.03,
        od_r=0.03,
        n_gdna_obs=0.0,
        n_rna_obs=1e4,
        is_region=reg,
        locked=unl,
    )
    assert np.all(tau_nog == 0.0)


def test_strand_evidence_struct_lock_regions_only():
    """I_struct (struct_lock) is composition-certainty for LOCKED REGION regions only — never a boundary boundary."""
    z = np.zeros(4)
    _, lock = strand_evidence(
        z,
        z,
        np.full(4, 0.5),
        kappa=0.5,
        od_g=0.03,
        od_r=0.03,
        n_gdna_obs=1e4,
        n_rna_obs=1e4,
        is_region=np.array([True, True, False, False]),
        locked=np.array([True, False, True, False]),
    )
    assert list(lock) == [True, False, False, False]


def test_strand_deconv_single_strand_solves_and_is_precise():
    """A single-strand (TRAPS: one-thing-varied) exon self-solves f_g from the tilt: it carries a live gDNA + sense-RNA own belief
    at finite precision, and NO antisense (the − axis is structurally dead). The strand λ-term (c·a²) applies
    to a single-strand region — the tilt is locked, so the strand PINS f_g (approach E)."""
    ni, _ = _init(kappa=0.95)
    ex = 2  # exon+ region slot
    assert not ni.struct_lock[ex]
    assert ni.tau_lam[ex] > 0.0  # stranded evidence fires (single-strand ⇒ strand pins f_g)
    assert ni.rho_g[ex] > 0.0 and ni.prec_g[ex] > 0.0
    assert ni.rho_pos[ex] > 0.0 and ni.prec_pos[ex] > 0.0
    assert ni.rho_neg[ex] == 0.0 and ni.prec_neg[ex] == 0.0  # − strand dead


def test_ambig_stranded_strand_gives_zero_fg_precision():
    """APPROACH E (the rank-1 fix): the strand Beta-Binomial is rank-1 (informs only p), so for an AMBIG (2-DOF)
    region — where the tilt is a free nuisance — the strand CANCELS out of f_g (the Schur-marginal λ-precision)
    and contributes ZERO. The un-gated `strand_evidence` still computes a POSITIVE single-strand λ-term for the
    AMBIG region (the phantom source), but `build_region_init` GATES it to single-strand regions, so the AMBIG region's
    τ_λ from the strand is 0. Only a density (gDNA) prior can pin an AMBIG region's f_g."""
    chain, statics, geometry, belief, ra = _scenario(kappa=0.9)
    am = 4  # AMBIG region slot (both strands live); no intron prior in _init ⇒ no density evidence
    is_reg = np.asarray(chain.kind) == REGION
    locked = ~(
        (np.asarray(statics.free_pos, bool) | np.asarray(statics.free_neg, bool))
        & (np.asarray(geometry.unspliced_count, float).sum(axis=1) > 0.0)
    )
    i_strand, _ = strand_evidence(
        np.asarray(geometry.unspliced_count, float)[:, 0],
        np.asarray(geometry.unspliced_count, float)[:, 1],
        np.full(chain.n_slots, 0.5),
        kappa=0.9,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=230.0,
        n_rna_obs=85.0,
        is_region=is_reg,
        locked=locked,
    )
    assert (
        i_strand[am] > 0.0
    )  # the un-gated single-strand strand term is positive (the phantom the fix removes)
    ni, _ = _init(kappa=0.9)
    assert not ni.struct_lock[am]
    assert (
        ni.tau_lam[am] == 0.0
    )  # ...but the assembled τ_λ gates the strand term to 0 for the AMBIG region
    assert ni.prec_g[am] == 0.0


# ── source 1: MEASURED (structural pure gDNA) → Poisson precision ────────────────────────────────────────────


def test_measured_intergenic_is_poisson_precision():
    """An intergenic region is structurally pure gDNA (f_g=1, struct_lock, composition variance 0), so its own
    gDNA precision is EXACTLY the Poisson count n — the anchor the whole prior-free pass leans on."""
    ni, n_region = _init(kappa=0.9)
    ig = 0  # intergenic region slot
    assert ni.struct_lock[ig]
    assert ni.f_g[ig] == 1.0
    v_fg, _ = own_composition_logvar(
        ni.f_g[ig : ig + 1], ni.tau_lam[ig : ig + 1], ni.struct_lock[ig : ig + 1]
    )
    assert v_fg[0] == 0.0  # composition CERTAIN
    assert np.isclose(ni.prec_g[ig], n_region[ig])  # precision == the raw count (Poisson)
    assert ni.prec_g[ig] > 0.0


# ── source 4: UNSOLVED → 100% gDNA at ZERO precision ─────────────────────────────────────────────────────────


def test_unsolved_ambig_unstranded_is_zero_precision():
    """An AMBIG region on unstranded data has NO intrinsic gDNA/RNA signal (I_strand=0 by the deadband) and no
    structural lock ⇒ every own precision is 0 (the honest 'no information' default), with no nan."""
    ni, _ = _init(kappa=0.5)  # unstranded
    am = 4  # AMBIG region slot
    assert not ni.struct_lock[am]
    assert ni.tau_lam[am] == 0.0
    assert ni.prec_g[am] == 0.0 and ni.prec_pos[am] == 0.0 and ni.prec_neg[am] == 0.0
    for arr in (ni.rho_g, ni.rho_pos, ni.rho_neg, ni.prec_g, ni.prec_pos, ni.prec_neg):
        assert not np.any(np.isnan(arr))


def test_build_region_init_all_finite_precisions():
    """No init produces a nan/negative precision (the 0·∞ hazard the arithmetic guards against)."""
    for kappa in (0.5, 0.9, 0.99):
        ni, _ = _init(kappa=kappa)
        for arr in (ni.prec_g, ni.prec_pos, ni.prec_neg):
            assert np.all(np.isfinite(arr)) and np.all(arr >= 0.0)


# ── source 2: density-deconvolution factor precision (I_density) ──────────────────────────────────────────────────────────


def test_density_factor_precision_flat_carries_no_evidence():
    """A FLAT λ-factor row carries τ = 0 — NOT the solve grid's own width; None ⇒ None (factory off)."""
    lam, _ = _logodds_grid(60, 10.0)
    assert density_factor_precision(None, lam) is None
    assert np.all(density_factor_precision(np.zeros((4, lam.shape[0])), lam) == 0.0)


def test_density_factor_precision_tracks_curvature_and_count():
    """A sharper factor carries more evidence, and (via NegBinom Var(g)=μ+μ²/α_eff) a high-count intron peels
    more sharply than a low-count one — the self-limiting precision (migrated from `_lambda_factor_precision`)."""
    lam, fg = _logodds_grid(60, 10.0)
    sharp = -0.5 * ((lam - 1.0) ** 2) / 0.05
    diffuse = -0.5 * ((lam - 1.0) ** 2) / 5.0
    assert (
        density_factor_precision(np.stack([sharp]), lam)[0]
        > density_factor_precision(np.stack([diffuse]), lam)[0]
        > 0.0
    )

    bg = GdnaBackground(
        log_mu_bg=float(np.log(0.01)),
        alpha=np.inf,
        sg=1.0e5,
        n0=0.0,
        n_regions=500,
        informative=True,
    )
    eff = np.array([1.0e3, 1.0e5])
    factor = density_lambda_factor(bg, count=0.02 * eff, eff_g=eff, fg_grid=fg)
    tau = density_factor_precision(factor, lam)
    assert tau[1] > tau[0] > 0.0


def test_density_factor_precision_flows_into_region_init():
    """End-to-end: passing an intron λ-factor lifts the intron's own gDNA precision above its strand-only value
    (the factory learning, registered as τ so the intron can propagate)."""
    chain, statics, geometry, belief, region_arrays = _scenario(
        kappa=0.5
    )  # unstranded ⇒ strand τ=0
    # a sharp λ-factor on the AMBIG region (id 5) — stand in for a confident intron peel
    lam, _ = _logodds_grid(60, 10.0)
    prior = np.zeros((chain.n_slots, lam.shape[0]))
    prior[4] = -0.5 * ((lam - 2.0) ** 2) / 0.05
    common = dict(
        kappa=0.5,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=230.0,
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        n_tilt=None,
        n_grid_ss=256,
        belief=belief,
    )
    ni_off = build_region_init(chain, statics, geometry, **common)
    ni_on = build_region_init(chain, statics, geometry, intron_prior=prior, **common)
    assert ni_off.tau_lam[4] == 0.0 and ni_off.prec_g[4] == 0.0  # unstranded, no factory ⇒ silent
    assert ni_on.tau_lam[4] > 0.0 and ni_on.prec_g[4] > 0.0  # factory ⇒ the region can now speak


# ── the scope of struct_lock: STRUCTURALLY PURE gDNA, never merely EMPTY ─────────────────────────────────


#: ⛔⛔⛔ THE DEFECT THESE THREE GATE IS PROVEN AND ITS FIX IS PANEL-NEGATIVE ALONE — so they are STRICT
#: xfails, not deletions, and they must go GREEN rather than be widened.
#:
#: ``struct_lock = locked & is_region`` with ``locked = ~solvable`` declares composition CERTAINTY at every
#: zero-count REGION: measured **19,709** slots on the gDNA ladder against **1,312** that are actually
#: ``g1_locked``, so **18,397** empty exons and introns claim a certainty they have not earned. That
#: contradicts ``strand_evidence``'s own docstring ("scoped to true intergenic REGION regions") and bypasses
#: `region_geometry.g1_locked`, the designated ONE HOME for the predicate.
#:
#: ⚠ It was INERT until 2026-08-06: ``own_precision``'s ``n > 0`` gate silenced every zero-count slot, so the
#: certainty could not leave them. Removing that gate un-masked it.
#:
#: ⛔ Scoping it to ``g1_locked ∧ REGION`` was PROTOTYPED AND MEASURED on the ladder
#: (`ladder_arm_ab.py --arm zc_struct_lock_g1`): the stranded × capture-ON row it was aimed at moved only
#: **−1.2 %**, ``g98`` went **+0.4 %** (worse), and the zero-gDNA control went **+3,207 %** (2,103 →
#: 69,532 fragments) — the mis-scoped mask is load-bearing for the zero-gDNA win. TRAPS: a-cancelling-defect-pair: it is half
#: of a cancelling pair, and the other half is the ``intergenic|exon`` boundary claiming its whole
#: RNA-contaminated crossing mass as gDNA. Price them TOGETHER or not at all.
_STRUCT_LOCK_XFAIL = pytest.mark.xfail(
    strict=True,
    reason="struct_lock is ~solvable & REGION, not g1_locked & REGION — proven, and the scoping fix is "
    "panel-negative alone (zero-gDNA control +3,207 %). Must go green with the boundary-composition arm.",
)


def _empty_exon_scenario(kappa=0.9):
    """The same three-region chain, with the AMBIG EXON's counts set to ZERO.

    ⭐ It must go through `build_region_init`, not through `strand_evidence` directly. The pre-existing gate
    `test_strand_evidence_struct_lock_regions_only` hands ``locked`` in as an argument, so it re-derives the
    caller's own input and cannot see what the caller computes — TRAPS: a-test-that-redefines, and the same hole TRAPS: perturb-every-gate's P3
    perturbation found in ``own_precision``.
    """
    parts = make_chain_parts(
        [TS_NONE, BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG],
        region_size_bp=[1500.0, 900.0, 1200.0],
        region_pos=[120.0, 200.0, 0.0],  # ⭐ the AMBIG exon is EMPTY
        region_neg=[110.0, 9.0, 0.0],
        gdna_fl=_delta_pmf(200),
        rna_fl=_delta_pmf(100),
    )
    belief = init_beliefs(
        parts.chain,
        parts.geometry,
        parts.statics,
        rna_sense_frac=kappa,
        gdna_strand_overdispersion=0.2,
        rna_strand_overdispersion=0.1,
        n_grid=60,
        n_grid_ss=256,
    )
    ni = build_region_init(
        parts.chain,
        parts.statics,
        parts.geometry,
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=230.0,
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        n_tilt=None,
        n_grid_ss=256,
        belief=belief,
    )
    return ni, parts


@_STRUCT_LOCK_XFAIL
def test_an_EMPTY_exon_is_not_composition_certain():
    """⛔⛔ ``struct_lock`` MEANS "STRUCTURALLY PURE gDNA". AN EMPTY EXON IS NOT THAT — IT IS UNMEASURED.

    ``strand_evidence`` is handed ``locked = ~solvable`` with
    ``solvable = (free_pos | free_neg) & (n_region > 0)``, so ``struct_lock = locked & is_region`` was true at
    **every zero-count REGION** — an exon and an intron included. ``own_composition_logvar`` then returns
    ``Var(log f_g) = 0`` there: composition CERTAIN, on a slot whose ``f_g`` is a default belief and whose
    evidence is nothing at all.

    ⭐ The predicate the docstring names is `region_geometry.g1_locked` — neither RNA strand admissible — and
    it is the ONE HOME for it (TRAPS: a-test-that-redefines). An AMBIG exon frees both strands, so it can never be TRAPS: no-magic-numbers
    however few fragments land in it.
    """
    ni, parts = _empty_exon_scenario()
    am = 4  # the AMBIG exon REGION slot
    assert np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)[am] == 0.0  # it IS empty
    assert not g1_locked(parts.statics.free_pos, parts.statics.free_neg)[
        am
    ]  # both strands admissible
    assert not ni.struct_lock[am], "an EMPTY exon was declared composition-CERTAIN"


@_STRUCT_LOCK_XFAIL
def test_struct_lock_is_exactly_g1_locked_on_the_region_axis():
    """The mask must equal ``g1_locked ∧ REGION``, computed from a DIFFERENT expression than the solver's.

    ⭐ ``g1_locked`` is production code and is imported, not restated — so this cannot drift the way TRAPS: a-test-that-redefines'
    two-homes predicate did. The `~np.asarray(...)` form here is the *consumer* side; the point of the gate
    is that the two agree on a fixture where ``~solvable`` and ``g1_locked`` DISAGREE, which is what the
    empty exon supplies."""
    ni, parts = _empty_exon_scenario()
    want = g1_locked(parts.statics.free_pos, parts.statics.free_neg) & (
        np.asarray(parts.chain.kind) == REGION
    )
    assert list(np.asarray(ni.struct_lock, bool)) == list(want)


#: ⛔⛔⛔ THE SECOND, LARGER DEFECT — AND THE PERTURBATION OF THE THREE ABOVE IS WHAT FOUND IT.
#: Scoping ``struct_lock`` to G1 was expected to restore the composition half of ``Var(log ρ_tot)`` at an
#: empty exon. It does not, and the reason is a CORNER DEGENERACY that has nothing to do with
#: ``struct_lock``: an evidence-free slot's default belief is ``f_g = 1`` exactly (the unsolved 100 %-gDNA
#: init), and ``region_sweep`` caps ``Var(f_g)`` at ``f_g(1−f_g)`` — **which is 0 at the corner.** So
#: ``logvar_tot`` at every evidence-free slot is ``count_logvar(n)`` and nothing else, and with
#: ``count_logvar`` finite there is now NO damping there at all. The ``1/n = ∞`` that the 39 % win removed
#: had been the only thing damping those hops.
#:
#: ⭐ The derived repair needs no tuned constant: with ``τ_λ = 0`` there is no evidence, so ``Var(f_g)`` is
#: the variance of ψ's OWN reference — ``Beta(½, ½)`` (`simplex_logodds._JEFFREYS_REF`) — which is exactly
#: **1/8**, and is nonzero wherever the point estimate sits. ``f_g(1−f_g)`` is a BERNOULLI variance and is
#: the right scale only away from the corner; at the corner it asserts a composition known exactly on a slot
#: with no data. ⭐ And the coefficient it multiplies, ``[(1/E_g − 1/E_r)/B]²``, diverges as ``E_g``
#: collapses — so the restored damping is largest exactly at the SHORT capture-depleted slots the stranded ×
#: capture-ON regression lives on.
_CORNER_VAR_XFAIL = pytest.mark.xfail(
    strict=True,
    reason="Var(f_g) is capped at f_g(1-f_g), which is 0 at the f_g=1 default of an evidence-free slot, "
    "so logvar_tot carries no composition term there and the hop is undamped.",
)


@_CORNER_VAR_XFAIL
def test_an_evidence_free_slot_is_not_certain_of_its_composition():
    """⭐⭐ THE OBSERVABLE: ``Var(f_g)`` at a slot with ``τ_λ = 0`` must be the reference prior's variance,
    not zero.

    ⛔ Gated as a strict inequality against ``0``, and separately against the STRUCTURALLY certain twin in
    the same fixture, so it cannot be satisfied by making everything uncertain: the intergenic region is
    genuinely TRAPS: no-magic-numbers and must still be exactly certain.

    ⚠ This reads `sweep.solve_chain`'s expression, restated here because it is a LOCAL there. That is a
    second home for one predicate (TRAPS: a-test-that-redefines) and is the reason the repair must move it into a named
    function beside `own_composition_logvar` when it lands."""
    ni, parts = _empty_exon_scenario()
    ig, am = 0, 4  # the intergenic REGION (truly TRAPS: no-magic-numbers) and the EMPTY AMBIG exon
    fg = np.asarray(ni.f_g, float)
    tau = np.asarray(ni.tau_lam, float)
    assert tau[am] == 0.0 and fg[am] == 1.0  # evidence-free, and parked at the corner
    fgfr = fg * (1.0 - fg)
    var_fg = np.where(
        np.asarray(ni.struct_lock, bool),
        0.0,
        np.where(tau > 1e-9, np.minimum(fgfr * fgfr / np.maximum(tau, 1e-9), fgfr), fgfr),
    )
    assert var_fg[ig] == 0.0  # structural certainty is real and must survive
    assert var_fg[am] > 0.0, "an evidence-free slot declares Var(f_g) = 0 at the f_g = 1 corner"


@_CORNER_VAR_XFAIL
def test_an_evidence_free_slot_pays_a_composition_transfer_variance():
    """⭐⭐⭐ AND THE CONSEQUENCE, in the currency the relay actually spends: ``σ²_transfer``.

    ``region_sweep`` hands that ``Var(f_g)`` to `enrichment_frame.composition_logvar`, whose output IS
    ``σ²_transfer`` (via `transfer_logvar`). With the composition term at 0 the evidence-free slot pays only
    ``count_logvar(n)``, so every message it touches crosses at ``1/(1/p + trigamma(½))`` — which is why a
    zero-mass slot went from a relay BARRIER to a CONDUIT when ``1/n`` was replaced.

    ⛔ The intergenic twin must still pay counting alone: it is composition-CERTAIN by structure, and a
    repair that damps it too would be widening the gate rather than fixing the defect."""
    ni, parts = _empty_exon_scenario()
    ig, am = 0, 4
    E_g = np.asarray(parts.geometry.eff_gdna, float)
    E_r = np.asarray(parts.geometry.eff_rna, float)
    n = np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)
    fg = np.asarray(ni.f_g, float)
    fgfr = fg * (1.0 - fg)
    var_fg = np.where(np.asarray(ni.struct_lock, bool), 0.0, fgfr)  # tau = 0 at both slots
    lv = composition_logvar(fg, E_g, E_r, var_fg, n)
    assert lv[ig] == count_logvar(np.array([n[ig]]))[0]  # certain ⇒ counting only, unchanged
    assert lv[am] > count_logvar(np.array([n[am]]))[0], (
        "an evidence-free slot pays no composition transfer variance, so its hops are undamped"
    )
