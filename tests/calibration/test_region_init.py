"""Unit tests for the pass-0 per-slot INITIALIZATION (`calibration.region_init`).

One test per information source — MEASURED (the structural lock), INTRON FACTORY, STRAND DECONVOLUTION,
UNSOLVED default (100 % gDNA, ZERO evidence) — on the three things the self-solve publishes: the
composition mode ``f_*`` and the own evidence ``tau_lam``. These pin the self-solve
that seeds the sweep.
"""

from __future__ import annotations


import numpy as np

from rigel.calibration.region_geometry import g1_locked, init_beliefs
from _synthetic import make_chain_parts
from rigel.calibration.region_init import (
    build_region_init,
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


# ── source 3: strand deconvolution evidence ────────────────────────────────────────────────────────────────


def test_strand_evidence_deadband_kills_unstranded():
    """The derived deadband makes I_strand IDENTICALLY 0 on unstranded data (κ=½), >0 stranded; a gDNA-free
    library (N_gdna=0 ⇒ σ²_d→∞) gates it to 0 even when stranded (migrated from `_compile_strand_evidence`)."""
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    base = dict(od_g=0.03, od_r=0.03, n_gdna_obs=1e4, n_rna_obs=1e4)
    tau_unstr = strand_evidence(u, u, fg, kappa=0.5, **base)
    tau_str = strand_evidence(u, u, fg, kappa=0.99, **base)
    assert np.all(tau_unstr == 0.0)
    assert np.all(tau_str > 0.0)
    tau_nog = strand_evidence(
        u, u, fg, kappa=0.99, od_g=0.03, od_r=0.03, n_gdna_obs=0.0, n_rna_obs=1e4
    )
    assert np.all(tau_nog == 0.0)


def test_a_single_strand_slot_solves_and_is_precise():
    """A single-strand (TRAPS: one-thing-varied) exon self-solves f_g from the tilt: it carries a live gDNA + sense-RNA own belief
    with own evidence, and NO antisense (the − axis is structurally dead). The strand λ-term (c·a²) applies
    to a single-strand region — the tilt is locked, so the strand PINS f_g (approach E)."""
    ni, _ = _init(kappa=0.95)
    ex = 2  # exon+ region slot
    assert ni.tau_lam[ex] > 0.0  # stranded evidence fires (single-strand ⇒ strand pins f_g)
    assert 0.0 < ni.f_g[ex] < 1.0 and ni.f_pos[ex] > 0.0
    assert ni.f_neg[ex] == 0.0  # − strand dead


def test_ambig_stranded_strand_gives_zero_fg_precision():
    """APPROACH E (the rank-1 fix): the strand Beta-Binomial is rank-1 (informs only p), so for an AMBIG (2-DOF)
    region — where the tilt is a free nuisance — the strand CANCELS out of f_g (the Schur-marginal λ-precision)
    and contributes ZERO. The un-gated `strand_evidence` still computes a POSITIVE single-strand λ-term for the
    AMBIG region (the phantom source), but `build_region_init` GATES it to single-strand regions, so the AMBIG region's
    τ_λ from the strand is 0. Only a density (gDNA) prior can pin an AMBIG region's f_g."""
    chain, statics, geometry, belief, ra = _scenario(kappa=0.9)
    am = 4  # AMBIG region slot (both strands live); no intron prior in _init ⇒ no density evidence
    i_strand = strand_evidence(
        np.asarray(geometry.unspliced_count, float)[:, 0],
        np.asarray(geometry.unspliced_count, float)[:, 1],
        np.full(chain.n_slots, 0.5),
        kappa=0.9,
        od_g=0.2,
        od_r=0.1,
        n_gdna_obs=230.0,
        n_rna_obs=85.0,
    )
    assert (
        i_strand[am] > 0.0
    )  # the un-gated single-strand strand term is positive (the phantom the fix removes)
    ni, _ = _init(kappa=0.9)
    assert (
        ni.tau_lam[am] == 0.0
    )  # ...but the assembled τ_λ gates the strand term to 0 for the AMBIG region


# ── source 1: MEASURED (structural pure gDNA) → composition CERTAIN ──────────────────────────────────────────


def test_measured_intergenic_is_structurally_certain():
    """An intergenic region is structurally pure gDNA: neither strand admissible (`g1_locked`, the one
    predicate) and the self-solve keeps its signature-binary ``f_g = 1`` — the anchor the whole
    prior-free pass leans on, certain without any strand or factory evidence."""
    chain, statics, _geometry, _belief, _ra = _scenario(kappa=0.9)
    ni, _ = _init(kappa=0.9)
    ig = 0  # intergenic region slot
    assert g1_locked(statics.free_pos, statics.free_neg)[ig]
    assert ni.f_g[ig] == 1.0 and ni.f_pos[ig] == 0.0 and ni.f_neg[ig] == 0.0


# ── source 4: UNSOLVED → 100% gDNA at ZERO evidence ──────────────────────────────────────────────────────────


def test_unsolved_ambig_unstranded_has_zero_own_evidence():
    """An AMBIG region on unstranded data has NO intrinsic gDNA/RNA signal (I_strand=0 by the deadband) and no
    structural lock ⇒ its own evidence is 0 (the honest 'no information' default), with no nan anywhere."""
    ni, _ = _init(kappa=0.5)  # unstranded
    am = 4  # AMBIG region slot
    assert ni.tau_lam[am] == 0.0
    for arr in (ni.f_g, ni.f_pos, ni.f_neg, ni.tau_lam):
        assert not np.any(np.isnan(arr))


def test_build_region_init_all_finite():
    """No init produces a nan/negative evidence or an off-simplex composition."""
    for kappa in (0.5, 0.9, 0.99):
        ni, _ = _init(kappa=kappa)
        assert np.all(np.isfinite(ni.tau_lam)) and np.all(ni.tau_lam >= 0.0)
        for arr in (ni.f_g, ni.f_pos, ni.f_neg):
            assert np.all(np.isfinite(arr)) and np.all(arr >= 0.0) and np.all(arr <= 1.0)


# ── source 2: density-deconvolution factor precision (I_density) ──────────────────────────────────────────────────────────


def test_density_factor_precision_flat_carries_no_evidence():
    """A FLAT λ-factor row carries τ = 0 — NOT the solve grid's own width; None ⇒ None (factory off)."""
    lam, _ = _logodds_grid(60, 10.0)
    assert density_factor_precision(None, lam) is None
    assert np.all(density_factor_precision(np.zeros((4, lam.shape[0])), lam) == 0.0)


def test_density_factor_precision_tracks_curvature_and_count():
    """A sharper factor carries more evidence, and (via NegBinom Var(g)=μ+μ²/α_eff) a high-count intron deconvolves
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
        size=1.0e5 + 0.5,  # the Gamma-posterior shape: Σg + ½
        n_regions=500,
        informative=True,
    )
    eff = np.array([1.0e3, 1.0e5])
    factor = density_lambda_factor(bg, count=0.02 * eff, eff_g=eff, fg_grid=fg)
    tau = density_factor_precision(factor, lam)
    assert tau[1] > tau[0] > 0.0


def test_density_factor_precision_flows_into_region_init():
    """End-to-end: passing an intron λ-factor lifts the intron's own evidence above its strand-only value
    (the factory learning, registered as τ so the intron can propagate)."""
    chain, statics, geometry, belief, region_arrays = _scenario(
        kappa=0.5
    )  # unstranded ⇒ strand τ=0
    # a sharp λ-factor on the AMBIG region (id 5) — stand in for a confident intron deconvolve
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
    assert ni_off.tau_lam[4] == 0.0  # unstranded, no factory ⇒ silent
    assert ni_on.tau_lam[4] > 0.0  # factory ⇒ the region can now speak
