"""Unit tests for the pass-0 per-node INITIALIZATION (`calibration.node_init`).

One test per information source of `docs/CARRY_FORWARD.md` — MEASURED (Poisson
precision), INTRON FACTORY, STRAND DECONVOLUTION, UNSOLVED default (100% gDNA, ZERO precision) — plus the
pure precision arithmetic (`own_composition_logvar`, `own_precision`). These pin the self-solve that seeds
the unified pass-0 relay.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.node_chain import NODE, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.node_init import (
    build_node_init,
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
    TS_AMBIG,
    TS_NONE,
    TS_POS,
)


# ── helpers ────────────────────────────────────────────────────────────────────────────────────────────────


def _cview(n_pos, n_neg, spl_sense=None):
    n_pos = np.asarray(n_pos, float)
    n_neg = np.asarray(n_neg, float)
    n = n_pos.shape[0]
    spl = np.zeros(n) if spl_sense is None else np.asarray(spl_sense, float)
    return SimpleNamespace(
        n_unspliced_pos=n_pos,
        n_unspliced_neg=n_neg,
        n_spliced_sense=spl,
        n_spliced_antisense=np.zeros(n),
        mass_unspliced=n_pos + n_neg,
        mass_spliced=spl,
        n_unspliced=n_pos + n_neg,
    )


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _scenario(kappa=0.9):
    """One ref, 3 regions: intergenic (pure gDNA) | exon+ (single-strand) | AMBIG. Region node ids 1, 3, 5."""
    rro = np.array([0, 3])
    rbo = np.array([0, 4])
    chain = build_node_chain(rro, rbo)
    sig = np.array([TS_NONE, BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG], dtype=np.int64)
    sc = np.array([TS_NONE, TS_POS, TS_AMBIG], dtype=np.int8)
    region_arrays = SimpleNamespace(
        strand_class=sc,
        signature=sig,
        region_size_bp=np.array([1500.0, 900.0, 1200.0]),
        n_regions=3,
    )
    # intergenic: symmetric gDNA; exon+: sense-tilted RNA (few antisense); AMBIG: symmetric.
    contained = _cview([120.0, 200.0, 90.0], [110.0, 9.0, 92.0])
    substrate = SimpleNamespace(contained=contained)
    bsub = SimpleNamespace(
        left_region=np.array([-1, 0, 1, 2]),
        right_region=np.array([0, 1, 2, -1]),
        left=_cview([0.0] * 4, [0.0] * 4),
        right=_cview([0.0] * 4, [0.0] * 4),
        junction_strand=np.zeros(4, dtype=np.int8),
    )
    gdna_fl, rna_fl = _delta_pmf(200), _delta_pmf(100)
    geometry = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    statics = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain,
        substrate,
        bsub,
        region_arrays,
        rna_sense_frac=kappa,
        gdna_strand_overdispersion=0.2,
        rna_strand_overdispersion=0.1,
        n_grid=60,
        n_grid_ss=256,
        statics=statics,
    )
    return chain, statics, geometry, belief, region_arrays


def _init(kappa=0.9, n_gdna_obs=230.0):
    chain, statics, geometry, belief, _ = _scenario(kappa)
    ni = build_node_init(
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
    n_node = np.where(
        np.asarray(chain.kind) == NODE,
        geometry.n_unspl_left,
        geometry.n_unspl_left + geometry.n_unspl_right,
    )
    return ni, n_node


# ── the pure precision arithmetic ────────────────────────────────────────────────────────────────────────


def test_own_composition_logvar_three_states():
    """Locked → 0 (certain); τ>0 → the Jacobian ((1−f_g)²/τ, f_g²/τ); τ=0 → ∞ (unseen), never a nan."""
    fg = np.array([0.2, 0.2, 0.2])
    tau = np.array(
        [0.0, 4.0, 0.0]
    )  # node0 locked (τ irrelevant); node1 unlocked τ>0; node2 unlocked τ=0
    lock = np.array([True, False, False])
    v_fg, v_fr = own_composition_logvar(fg, tau, lock)
    assert v_fg[0] == 0.0 and v_fr[0] == 0.0  # locked ⇒ certain
    assert np.isclose(v_fg[1], (1 - 0.2) ** 2 / 4.0) and np.isclose(v_fr[1], 0.2**2 / 4.0)
    assert np.isinf(v_fg[2]) and np.isinf(v_fr[2])  # no evidence ⇒ unseen
    assert not np.any(np.isnan(v_fg)) and not np.any(np.isnan(v_fr))


def test_own_precision_monotone_and_zeros():
    """p = n/(n·v+1): monotone in the count; 0 at n=0 or v=∞ or dead component; no 0·∞ nan."""
    n = np.array([10.0, 400.0, 0.0, 50.0])
    v = np.array([0.1, 0.1, 0.1, np.inf])
    live = np.array([True, True, True, True])
    p = own_precision(n, v, live)
    assert p[1] > p[0] > 0.0  # more count ⇒ sharper
    assert p[2] == 0.0  # no count
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
    """I_struct (struct_lock) is composition-certainty for LOCKED NODE nodes only — never a boundary seam."""
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
    """A single-strand (G2) exon self-solves f_g from the tilt: it carries a live gDNA + sense-RNA own belief
    at finite precision, and NO antisense (the − axis is structurally dead). The strand λ-term (c·a²) applies
    to a single-strand node — the tilt is locked, so the strand PINS f_g (approach E)."""
    ni, _ = _init(kappa=0.95)
    ex = 3  # exon+ region node
    assert not ni.struct_lock[ex]
    assert ni.tau_lam[ex] > 0.0  # stranded evidence fires (single-strand ⇒ strand pins f_g)
    assert ni.rho_g[ex] > 0.0 and ni.prec_g[ex] > 0.0
    assert ni.rho_pos[ex] > 0.0 and ni.prec_pos[ex] > 0.0
    assert ni.rho_neg[ex] == 0.0 and ni.prec_neg[ex] == 0.0  # − strand dead


def test_ambig_stranded_strand_gives_zero_fg_precision():
    """APPROACH E (the rank-1 fix): the strand Beta-Binomial is rank-1 (informs only p), so for an AMBIG (2-DOF)
    node — where the tilt is a free nuisance — the strand CANCELS out of f_g (the Schur-marginal λ-precision)
    and contributes ZERO. The un-gated `strand_evidence` still computes a POSITIVE single-strand λ-term for the
    AMBIG node (the phantom source), but `build_node_init` GATES it to single-strand nodes, so the AMBIG node's
    τ_λ from the strand is 0. Only a density (gDNA) prior can pin an AMBIG node's f_g."""
    chain, statics, geometry, belief, ra = _scenario(kappa=0.9)
    am = 5  # AMBIG region node (both strands live), no intron prior in _init ⇒ no density evidence
    is_reg = np.asarray(chain.kind) == NODE
    locked = ~(
        (np.asarray(statics.free_pos, bool) | np.asarray(statics.free_neg, bool))
        & (np.asarray(statics.mass_unspliced, float) > 0.0)
    )
    i_strand, _ = strand_evidence(
        statics.u_pos,
        statics.u_neg,
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
    )  # ...but the assembled τ_λ gates the strand term to 0 for the AMBIG node
    assert ni.prec_g[am] == 0.0


# ── source 1: MEASURED (structural pure gDNA) → Poisson precision ────────────────────────────────────────────


def test_measured_intergenic_is_poisson_precision():
    """An intergenic node is structurally pure gDNA (f_g=1, struct_lock, composition variance 0), so its own
    gDNA precision is EXACTLY the Poisson count n — the anchor the whole prior-free pass leans on."""
    ni, n_node = _init(kappa=0.9)
    ig = 1  # intergenic region node
    assert ni.struct_lock[ig]
    assert ni.f_g[ig] == 1.0
    v_fg, _ = own_composition_logvar(
        ni.f_g[ig : ig + 1], ni.tau_lam[ig : ig + 1], ni.struct_lock[ig : ig + 1]
    )
    assert v_fg[0] == 0.0  # composition CERTAIN
    assert np.isclose(ni.prec_g[ig], n_node[ig])  # precision == the raw count (Poisson)
    assert ni.prec_g[ig] > 0.0


# ── source 4: UNSOLVED → 100% gDNA at ZERO precision ─────────────────────────────────────────────────────────


def test_unsolved_ambig_unstranded_is_zero_precision():
    """An AMBIG node on unstranded data has NO intrinsic gDNA/RNA signal (I_strand=0 by the deadband) and no
    structural lock ⇒ every own precision is 0 (the honest 'no information' default), with no nan."""
    ni, _ = _init(kappa=0.5)  # unstranded
    am = 5  # AMBIG region node
    assert not ni.struct_lock[am]
    assert ni.tau_lam[am] == 0.0
    assert ni.prec_g[am] == 0.0 and ni.prec_pos[am] == 0.0 and ni.prec_neg[am] == 0.0
    for arr in (ni.rho_g, ni.rho_pos, ni.rho_neg, ni.prec_g, ni.prec_pos, ni.prec_neg):
        assert not np.any(np.isnan(arr))


def test_build_node_init_all_finite_precisions():
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


def test_density_factor_precision_flows_into_node_init():
    """End-to-end: passing an intron λ-factor lifts the intron's own gDNA precision above its strand-only value
    (the factory learning, registered as τ so the intron can propagate)."""
    chain, statics, geometry, belief, region_arrays = _scenario(
        kappa=0.5
    )  # unstranded ⇒ strand τ=0
    # a sharp λ-factor on the AMBIG node (id 5) — stand in for a confident intron peel
    lam, _ = _logodds_grid(60, 10.0)
    prior = np.zeros((chain.n_slots, lam.shape[0]))
    prior[5] = -0.5 * ((lam - 2.0) ** 2) / 0.05
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
    ni_off = build_node_init(chain, statics, geometry, **common)
    ni_on = build_node_init(chain, statics, geometry, intron_prior=prior, **common)
    assert ni_off.tau_lam[5] == 0.0 and ni_off.prec_g[5] == 0.0  # unstranded, no factory ⇒ silent
    assert ni_on.tau_lam[5] > 0.0 and ni_on.prec_g[5] > 0.0  # factory ⇒ the node can now speak
