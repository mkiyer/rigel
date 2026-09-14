"""The pass-0 per-slot initialization, and the signature classification its chain slots carry.

One test per information source — MEASURED (the structural lock), INTRON FACTORY, STRAND DECONVOLUTION,
UNSOLVED default (100 % gDNA, ZERO evidence) — on the three things the self-solve publishes: the
composition mode ``f_*`` and the own evidence ``tau_lam``. These pin the self-solve that seeds the
sweep. The last block gates the ``mrna_active_*`` classification ``build_region_statics`` carries onto
the chain alongside the ``free_*`` (nascent-active / RNA-crossing) masks; that classification is
signature-only, so the counts are irrelevant to it, and it is asserted across every region type and all
four boundary types.
"""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.density_deconv import (
    GdnaBackground,
    density_factor_precision,
    density_lambda_factor,
)
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import build_region_statics, g1_locked, init_beliefs
from rigel.calibration.region_init import (
    build_region_init,
    strand_discriminability,
    strand_evidence,
)
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    TS_NONE,
    transcript_strand_class,
)
from rigel.calibration.simplex_logodds import _logodds_grid

from _synthetic import make_chain_parts


# ── helpers ────────────────────────────────────────────────────────────────────────────────────────────────


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _scenario(kappa=0.9):
    """One ref, 3 regions: intergenic (pure gDNA) | exon+ (single-strand) | AMBIG.

    The chain is ``N E N E N`` — 5 slots, the regions at 0/2/4 — and the two boundaries carry no
    counts, so this exercises the region axis alone. There are no data-free terminal slots: a
    reference with k regions owns k − 1 boundaries.
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
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _init(kappa=0.9):
    chain, statics, geometry, belief, _ = _scenario(kappa)
    ni = build_region_init(
        statics,
        geometry,
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        belief=belief,
    )
    # ONE count per slot: it is both the density numerator and the Poisson n.
    return ni, np.asarray(geometry.unspliced_count, float).sum(axis=1)


# ── source 3: strand deconvolution evidence ────────────────────────────────────────────────────────────────


def test_strand_evidence_is_zero_unstranded_and_positive_stranded():
    """I_strand is identically 0 on unstranded data (κ = ½) and positive on stranded data."""
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    assert np.all(strand_evidence(u, u, fg, kappa=0.5, od_r=0.03, n_rna_obs=1e4) == 0.0)
    assert np.all(strand_evidence(u, u, fg, kappa=0.99, od_r=0.03, n_rna_obs=1e4) > 0.0)


def test_a_gdna_free_stranded_library_keeps_its_strand_channel():
    """gDNA's strand mean is ½ by symmetry and needs no observation, so the channel's liveness is a
    function of the RNA strand fit alone: neither `strand_discriminability` nor `strand_evidence` takes a
    gDNA count or the gDNA overdispersion (the invariant, asserted structurally: a ``1/N_gdna`` term had
    switched every gDNA-free library's channel off, the modal real case), and a stranded library is live
    at any gDNA content, zero included."""
    for fn in (strand_discriminability, strand_evidence):
        names = set(inspect.signature(fn).parameters)
        assert not names & {"n_gdna_obs", "od_g"}, f"{fn.__name__} reads a gDNA quantity: {names}"
    assert strand_discriminability(0.99, 1e4) > 0.0
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    assert np.all(strand_evidence(u, u, fg, kappa=0.99, od_r=0.03, n_rna_obs=1e4) > 0.0)


def test_the_strand_channel_is_a_protocol_decision_not_a_sampling_band():
    """The ladder's unstranded zero control fits κ̂ = 0.500298 from 3,090,137 spliced fragments, 1.05σ
    from ½. An unbiased estimate of (κ−½)² floored at zero — the form this replaced — reads that LIVE (a
    coin toss lost on 32 % of unstranded libraries; 21,484 false gDNA fragments on this one), while the
    Bayes factor of a free κ against κ = ½ exactly reads it dead by 6.7 nats and the stranded rows live
    by 10⁶ nats. At κ = ½ exactly the channel is dead at every N."""
    kappa, n = 0.500298, 3_090_137.0
    assert strand_discriminability(kappa, n) == 0.0
    assert 4.0 * ((kappa - 0.5) ** 2 - 0.25 / n) > 0.0, (
        "the replaced form's verdict on the same fit"
    )
    assert strand_discriminability(0.009883, 3_093_241.0) == pytest.approx(
        4.0 * (0.009883 - 0.5) ** 2
    )
    for n in (10.0, 1e3, 1e6):
        assert strand_discriminability(0.5, n) == 0.0


def test_the_occam_penalty_grows_with_the_spliced_sample():
    """ln BF₁₀ ≈ ½·[z² − ln(2N/π)]: the excursion a sampling fluctuation must clear to read as a protocol
    grows with the sample — 3σ is live at N = 100 and dead at N = 10⁶, 4σ live at 10⁶ and dead at 10⁸ —
    so no multiple of σ is the gate; the free parameter's Occam penalty is."""

    def at(z, n):
        return strand_discriminability(0.5 + z * np.sqrt(0.25 / n), n)

    assert at(3.0, 100.0) > 0.0 and at(3.0, 1e6) == 0.0
    assert at(4.0, 1e6) > 0.0 and at(4.0, 1e8) == 0.0


def test_no_spliced_observations_is_no_protocol_decision():
    """With no spliced fragment the two hypotheses have equal marginal likelihood, so the channel is dead
    whatever κ a caller claims (`calibrate` raises before this on a real library)."""
    assert strand_discriminability(0.99, 0.0) == 0.0
    assert strand_discriminability(0.5, 0.0) == 0.0


def test_a_single_strand_slot_solves_and_is_precise():
    """A single-strand exon self-solves f_g from the tilt: it carries a live gDNA + sense-RNA own
    belief with own evidence, and no antisense, because the − axis is structurally dead. The strand
    λ-term applies to a single-strand region: the tilt is locked, so the strand pins f_g."""
    ni, _ = _init(kappa=0.95)
    ex = 2  # exon+ region slot
    assert ni.tau_lam[ex] > 0.0  # stranded evidence fires (single-strand ⇒ strand pins f_g)
    assert 0.0 < ni.f_g[ex] < 1.0 and ni.f_pos[ex] > 0.0
    assert ni.f_neg[ex] == 0.0  # − strand dead


def test_ambig_stranded_strand_gives_zero_fg_precision():
    """The strand Beta-Binomial is rank-1 — it informs only p — so at an AMBIG region, where the
    tilt is a free nuisance, the strand cancels out of f_g in the Schur-marginal λ-precision and
    contributes zero. The un-gated `strand_evidence` still computes a positive single-strand λ-term
    there, which is the phantom source, but `build_region_init` gates it to single-strand regions,
    so the AMBIG region's τ_λ from the strand is 0. Only a density prior can pin its f_g."""
    chain, statics, geometry, belief, ra = _scenario(kappa=0.9)
    am = 4  # AMBIG region slot (both strands live); no intron prior in _init ⇒ no density evidence
    i_strand = strand_evidence(
        np.asarray(geometry.unspliced_count, float)[:, 0],
        np.asarray(geometry.unspliced_count, float)[:, 1],
        np.full(chain.n_slots, 0.5),
        kappa=0.9,
        od_r=0.1,
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
    """An AMBIG region on unstranded data has no intrinsic gDNA/RNA signal (I_strand=0 by the
    protocol decision) and no structural lock, so its own evidence is 0 — the honest "no information"
    default — with no nan anywhere."""
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
    """A flat λ-factor row carries τ = 0, not the solve grid's own width; None ⇒ None (factory off)."""
    lam, _ = _logodds_grid(60, 10.0)
    assert density_factor_precision(None, lam) is None
    assert np.all(density_factor_precision(np.zeros((4, lam.shape[0])), lam) == 0.0)


def test_density_factor_precision_tracks_curvature_and_count():
    """A sharper factor carries more evidence, and (via NegBinom Var(g)=μ+μ²/α_eff) a high-count
    intron deconvolves more sharply than a low-count one — the self-limiting precision."""
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
        n_rna_obs=85.0,
        n_grid=60,
        logodds_window=10.0,
        belief=belief,
    )
    ni_off = build_region_init(statics, geometry, **common)
    ni_on = build_region_init(statics, geometry, intron_prior=prior, **common)
    assert ni_off.tau_lam[4] == 0.0  # unstranded, no factory ⇒ silent
    assert ni_on.tau_lam[4] > 0.0  # factory ⇒ the region can now speak


# ── the signature-only classification the chain carries: free_* beside mrna_active_* ─────────


def _substrate(n_regions, n_boundaries):
    """A zero-count substrate of the right shape — the classifier under test is signature-only, so the
    counts are irrelevant and only ``boundary_spliced`` is read (for ``spliced_count``)."""

    def view(n):
        return SimpleNamespace(count=np.zeros((n, 2)))

    return SimpleNamespace(
        region_contained=view(n_regions),
        boundary_unspliced=view(n_boundaries),
        boundary_spliced=view(n_boundaries),
    )


def _build_statics(region_sigs):
    """Build a single-reference chain over ``region_sigs`` (genomic order) and return
    ``(chain, statics)``.

    ``k`` regions own ``k − 1`` interior boundaries and there are no terminal slots, so there is no
    "boundary 0 is the reference-start sink" case to classify. A boundary always has a region on
    both sides, which is what lets `build_region_statics` read its flanks straight off the chain's
    adjacency instead of a ``left_region``/``right_region`` array with ``-1`` holes.
    """
    sig = np.asarray(region_sigs, dtype=np.int64)
    n_reg = sig.shape[0]
    region_arrays = SimpleNamespace(
        strand_class=transcript_strand_class(sig),
        signature=sig,
        region_size_bp=np.full(n_reg, 1000.0),
    )
    n_boundary = max(n_reg - 1, 0)
    chain = build_region_chain(np.array([0, n_reg]), np.array([0, n_boundary]))
    return chain, build_region_statics(chain, region_arrays)


def test_classifier_covers_region_and_boundary_types():
    # N0 exon+ | N1 exon+ | N2 intron+ | N3 intergenic | N4 ambig-exon | N5 ambig-exon
    # boundaries (5): E0(N0|N1) E1(N1|N2) E2(N2|N3) E3(N3|N4) E4(N4|N5) — no terminals
    sigs = [
        BIT_EXON_POS,
        BIT_EXON_POS,
        BIT_INTRON_POS,
        0,
        BIT_EXON_POS | BIT_EXON_NEG,
        BIT_EXON_POS | BIT_EXON_NEG,
    ]
    chain, st = _build_statics(sigs)

    # masks are bool and full-length; the whole-chain invariant mrna_active ⇒ free (nascent) holds.
    for m in (st.free_pos, st.free_neg, st.mrna_active_pos, st.mrna_active_neg):
        assert m.dtype == bool and m.shape[0] == st.n_slots
    assert np.all(~st.mrna_active_pos | st.free_pos)  # mature ⇒ nascent-active (+)
    assert np.all(~st.mrna_active_neg | st.free_neg)  # (−)

    kind, ref = np.asarray(chain.kind), np.asarray(chain.obj_idx)
    reg = np.where(kind == REGION)[0]  # N0..N5 (genomic order)
    bnd = np.where(kind == BOUNDARY)[0]  # E0..E4
    np.testing.assert_array_equal(ref[reg], np.arange(6))  # confirm genomic ordering
    np.testing.assert_array_equal(ref[bnd], np.arange(5))

    def state(i):
        return (
            bool(st.free_pos[i]),
            bool(st.free_neg[i]),
            bool(st.mrna_active_pos[i]),
            bool(st.mrna_active_neg[i]),
        )

    # --- regions (free_pos, free_neg, mrna_pos, mrna_neg) ---
    assert state(reg[0]) == (True, False, True, False)  # exon+  : mature-capable +
    assert state(reg[2]) == (
        True,
        False,
        False,
        False,
    )  # intron+: NASCENT-ONLY + (free but not mature)
    assert state(reg[3]) == (False, False, False, False)  # intergenic: gDNA sink
    assert state(reg[4]) == (True, True, True, True)  # ambig-exon: mature-capable both strands

    # --- boundaries: the four types. There is no reference-start terminal slot, so E0 is the first
    # real boundary, N0|N1.
    assert state(bnd[0]) == (True, False, True, False)  # exon↔exon+   : MATURE-CAPABLE
    assert state(bnd[1]) == (True, False, False, False)  # exon↔intron+ : NASCENT-ONLY
    assert state(bnd[2]) == (False, False, False, False)  # intron↔intergenic: SINK (no + crossing)
    assert state(bnd[3]) == (False, False, False, False)  # intergenic↔ambig-exon : SINK
    assert state(bnd[4]) == (True, True, True, True)  # ambig↔ambig  : AMBIG, mature both strands
