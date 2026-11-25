"""The message-free self-solve that seeds the sweep, and the signature classification its chain slots carry.

One test per information source (measured: the structural lock; strand deconvolution, whose evidence the strand
decision gates; and the unsolved default, 100 % gDNA at zero evidence) on what the self-solve publishes: the composition
mode ``f_*`` and the own evidence ``tau_lam``. The last block tests the ``free_*`` (RNA-crossing) masks
``build_region_statics`` carries onto the chain; that classification is signature-only, so counts are irrelevant to it,
and it is asserted across every region type and all four boundary types. The strand decision itself is
``test_strand_balance.py``'s.
"""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np

from rigel.calibration.diagnostics import SweepDiagnostics
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.region_geometry import build_region_statics, g1_locked, init_beliefs
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
)
from rigel.calibration.strand_balance import StrandBalance, strand_discriminability
from rigel.calibration.sweep import solve_chain
from rigel._solve_impl import testing as R

from _synthetic import make_chain_parts


# ── helpers ────────────────────────────────────────────────────────────────────────────────────────────────


def strand_evidence(u_pos, u_neg, fg_loc, *, kappa, n_rna_obs):
    """The kernel's I_strand (`_solve_impl.testing.strand_evidence`): the reference-free strand composition
    evidence at ``fg_loc``, the protocol's discriminability its library-level factor."""
    disc = strand_discriminability(kappa, n_rna_obs)
    f64 = lambda a: np.ascontiguousarray(a, np.float64)  # noqa: E731
    return R.strand_evidence(f64(u_pos), f64(u_neg), f64(fg_loc), float(kappa), float(disc))


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _scenario(kappa=0.9):
    """One ref, 3 regions: intergenic (pure gDNA) | exon+ (single-strand) | AMBIG. The chain is ``N E N E N`` (5 slots,
    regions at 0/2/4), and the two boundaries carry no counts, so this exercises the region axis alone.
    """
    parts = make_chain_parts(
        [0, BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG],
        region_size_bp=[1500.0, 900.0, 1200.0],
        # intergenic: symmetric gDNA; exon+: sense-tilted RNA (few antisense); AMBIG: symmetric.
        region_pos=[120.0, 200.0, 90.0],
        region_neg=[110.0, 9.0, 92.0],
        gdna_fl=_delta_pmf(200),
        rna_fl=_delta_pmf(100),
    )
    belief = init_beliefs(
        parts.geometry,
        parts.statics,
        rna_sense_frac=kappa,
        n_grid=60,
        logodds_window=10.0,
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _init(kappa=0.9):
    """The message-free self-solve of the scenario, read off the sweep under the silent policy: no message is sent, so a
    slot's final solve is its self-solve (the same inputs to ψ), and the capture carries the own evidence ``tau_lam``
    and the self-solve's ``fg_loc`` (`diagnostics.SweepDiagnostics`)."""
    chain, statics, geometry, belief, ra = _scenario(kappa)
    cap = SweepDiagnostics()
    out = solve_chain(
        chain,
        statics,
        geometry,
        belief,
        ra,
        strand=StrandBalance(kappa, 85, False),
        n_grid=60,
        logodds_window=10.0,
        policy=SilentPolicy(),
        diagnostics=cap,
    )
    ni = SimpleNamespace(
        f_g=np.asarray(out.f_g),
        f_pos=np.asarray(out.f_pos),
        f_neg=np.asarray(out.f_neg),
        tau_lam=np.asarray(cap.tau_lam),
        fg_loc=np.asarray(cap.fg_loc),
    )
    # One count per slot: it is both the density numerator and the Poisson n.
    return ni, np.asarray(geometry.unspliced_count, float).sum(axis=1)


# ── source 3: strand deconvolution evidence ────────────────────────────────────────────────────────────────


def test_strand_evidence_is_zero_unstranded_and_positive_stranded():
    """I_strand is identically 0 on unstranded data (κ = ½) and positive on stranded data."""
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    assert np.all(strand_evidence(u, u, fg, kappa=0.5, n_rna_obs=1e4) == 0.0)
    assert np.all(strand_evidence(u, u, fg, kappa=0.99, n_rna_obs=1e4) > 0.0)


def test_a_gdna_free_stranded_library_keeps_its_strand_channel():
    """gDNA's strand mean is ½ by symmetry and needs no observation, so the channel's liveness depends on the RNA strand
    fit alone: neither `strand_discriminability` nor `strand_evidence` takes a gDNA quantity (asserted structurally),
    and a stranded library is live at any gDNA content, zero included.
    """
    names = set(inspect.signature(strand_discriminability).parameters)
    assert "n_gdna_obs" not in names, f"strand_discriminability reads a gDNA quantity: {names}"
    doc = R.strand_evidence.__doc__  # the kernel's binding, its arguments in its signature line
    assert "kappa" in doc and "gdna" not in doc.lower().split("disc")[0], doc
    assert strand_discriminability(0.99, 1e4) > 0.0
    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    assert np.all(strand_evidence(u, u, fg, kappa=0.99, n_rna_obs=1e4) > 0.0)


def test_a_single_strand_slot_solves_and_is_precise():
    """A single-strand exon self-solves f_g from the tilt: it carries a live gDNA + sense-RNA own belief with own
    evidence, and no antisense, because the − axis is structurally dead. The tilt is locked, so the strand λ-term pins
    f_g.
    """
    ni, _ = _init(kappa=0.95)
    ex = 2  # exon+ region slot
    assert ni.tau_lam[ex] > 0.0  # stranded evidence fires (single-strand ⇒ strand pins f_g)
    assert 0.0 < ni.f_g[ex] < 1.0 and ni.f_pos[ex] > 0.0
    assert ni.f_neg[ex] == 0.0  # − strand dead


def test_ambig_stranded_strand_gives_zero_fg_precision():
    """The strand Beta-Binomial is rank-1 (it informs only p), so at an AMBIG region, where the tilt is a free nuisance,
    the strand cancels out of f_g in the Schur-marginal λ-precision. The ungated `strand_evidence` still computes a
    positive single-strand λ-term there, but the kernel gates it to single-strand regions, so the AMBIG region's strand
    τ_λ is 0. Only a density prior can pin its f_g.
    """
    chain, statics, geometry, belief, ra = _scenario(kappa=0.9)
    am = 4  # AMBIG region slot (both strands live); no intron prior in _init ⇒ no density evidence
    i_strand = strand_evidence(
        np.asarray(geometry.unspliced_count, float)[:, 0],
        np.asarray(geometry.unspliced_count, float)[:, 1],
        np.full(chain.n_slots, 0.5),
        kappa=0.9,
        n_rna_obs=85.0,
    )
    assert i_strand[am] > 0.0  # the ungated single-strand strand term is positive
    ni, _ = _init(kappa=0.9)
    assert (
        ni.tau_lam[am] == 0.0
    )  # but the assembled τ_λ gates the strand term to 0 for the AMBIG region


# ── source 1: measured (structural pure gDNA) → composition certain ──────────────────────────────────────────


def test_measured_intergenic_is_structurally_certain():
    """An intergenic region is structurally pure gDNA: neither strand admissible (`g1_locked`) and the self-solve keeps
    its signature-binary ``f_g = 1``, the anchor the prior-free pass leans on, certain without strand or factory
    evidence."""
    chain, statics, _geometry, _belief, _ra = _scenario(kappa=0.9)
    ni, _ = _init(kappa=0.9)
    ig = 0  # intergenic region slot
    assert g1_locked(statics.free_pos, statics.free_neg)[ig]
    assert ni.f_g[ig] == 1.0 and ni.f_pos[ig] == 0.0 and ni.f_neg[ig] == 0.0


# ── source 4: unsolved → 100% gDNA at zero evidence ──────────────────────────────────────────────────────────


def test_unsolved_ambig_unstranded_has_zero_own_evidence():
    """An AMBIG region on unstranded data has no intrinsic gDNA/RNA signal (I_strand=0 by the protocol decision) and no
    structural lock, so its own evidence is 0, the "no information" default, with no nan anywhere."""
    ni, _ = _init(kappa=0.5)  # unstranded
    am = 4  # AMBIG region slot
    assert ni.tau_lam[am] == 0.0
    for arr in (ni.f_g, ni.f_pos, ni.f_neg, ni.tau_lam):
        assert not np.any(np.isnan(arr))


def test_the_self_solve_is_all_finite():
    """No self-solve produces a nan/negative evidence or an off-simplex composition."""
    for kappa in (0.5, 0.9, 0.99):
        ni, _ = _init(kappa=kappa)
        assert np.all(np.isfinite(ni.tau_lam)) and np.all(ni.tau_lam >= 0.0)
        for arr in (ni.f_g, ni.f_pos, ni.f_neg):
            assert np.all(np.isfinite(arr)) and np.all(arr >= 0.0) and np.all(arr <= 1.0)


# ── the signature-only classification the chain carries: free_* ───────────────────────────────


def _build_statics(region_sigs):
    """Build a single-reference chain over ``region_sigs`` (genomic order) and return ``(chain, statics)``. ``k``
    regions own ``k − 1`` interior boundaries and there are no terminal slots, so a boundary always has a region on both
    sides and `build_region_statics` reads its flanks straight off the chain's adjacency.
    """
    sig = np.asarray(region_sigs, dtype=np.int64)
    n_reg = sig.shape[0]
    region_arrays = SimpleNamespace(
        signature=sig,
        region_size_bp=np.full(n_reg, 1000.0),
    )
    n_boundary = max(n_reg - 1, 0)
    chain = build_region_chain(np.array([0, n_reg]), np.array([0, n_boundary]))
    return chain, build_region_statics(chain, region_arrays, np.zeros(n_boundary, np.uint16))


def test_classifier_covers_region_and_boundary_types():
    # N0 exon+ | N1 exon+ | N2 intron+ | N3 intergenic | N4 ambig-exon | N5 ambig-exon
    # boundaries (5): E0(N0|N1) E1(N1|N2) E2(N2|N3) E3(N3|N4) E4(N4|N5), no terminals
    sigs = [
        BIT_EXON_POS,
        BIT_EXON_POS,
        BIT_INTRON_POS,
        0,
        BIT_EXON_POS | BIT_EXON_NEG,
        BIT_EXON_POS | BIT_EXON_NEG,
    ]
    chain, st = _build_statics(sigs)

    # masks are bool and full-length
    for m in (st.free_pos, st.free_neg):
        assert m.dtype == bool and m.shape[0] == st.n_slots

    kind, ref = np.asarray(chain.kind), np.asarray(chain.obj_idx)
    reg = np.where(kind == REGION)[0]  # N0..N5 (genomic order)
    bnd = np.where(kind == BOUNDARY)[0]  # E0..E4
    np.testing.assert_array_equal(ref[reg], np.arange(6))  # confirm genomic ordering
    np.testing.assert_array_equal(ref[bnd], np.arange(5))

    def state(i):
        return (bool(st.free_pos[i]), bool(st.free_neg[i]))

    # --- regions (free_pos, free_neg) ---
    assert state(reg[0]) == (True, False)  # exon+
    assert state(reg[2]) == (True, False)  # intron+
    assert state(reg[3]) == (False, False)  # intergenic: gDNA sink
    assert state(reg[4]) == (True, True)  # ambig-exon: both strands

    # --- boundaries: the four types. E0, N0|N1, is the first boundary.
    assert state(bnd[0]) == (True, False)  # exon↔exon+
    assert state(bnd[1]) == (True, False)  # exon↔intron+
    assert state(bnd[2]) == (False, False)  # intron↔intergenic: sink (no + crossing)
    assert state(bnd[3]) == (False, False)  # intergenic↔ambig-exon: sink
    assert state(bnd[4]) == (True, True)  # ambig↔ambig: AMBIG
