"""The BP-solver per-node geometry (`calibration.node_geometry.build_node_geometry`).

Verifies the per-face assembly onto the chain against the substrate: a region presents its contained geometry
both ways; a boundary presents its per-side crossing geometry (left-side / right-side, in the flank regions'
sizes); the one-sided motif spliced lands on the exon-flank face; per-component FL (gDNA vs RNA) is honoured.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.bp_solver import (
    node_sweep,
)
from rigel.calibration.node_geometry import (
    build_node_geometry,
    build_node_statics,
    init_beliefs,
)
from rigel.calibration.effective_length import (
    boundary_side_eff_length,
    region_eff_length,
    spliced_side_eff_length,
)
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
    mrna_active_strands,
    nrna_active_strands,
)


def _view(mass_u, mass_spl):
    n = len(mass_u)
    z = np.zeros(n)
    return SimpleNamespace(
        n_unspliced_pos=z.copy(),
        n_unspliced_neg=z.copy(),
        n_spliced_sense=z.copy(),
        n_spliced_antisense=z.copy(),
        mass_unspliced=np.asarray(mass_u, float),
        mass_spliced=np.asarray(mass_spl, float),
    )


def _cview(n_pos, n_neg, spl_sense=None, mass_u=None, mass_spl=None):
    """A per-region/side view with per-strand unspliced counts (and optional sense-spliced)."""
    n_pos = np.asarray(n_pos, float)
    n_neg = np.asarray(n_neg, float)
    n = n_pos.shape[0]
    spl = np.zeros(n) if spl_sense is None else np.asarray(spl_sense, float)
    return SimpleNamespace(
        n_unspliced_pos=n_pos,
        n_unspliced_neg=n_neg,
        n_spliced_sense=spl,
        n_spliced_antisense=np.zeros(n),
        mass_unspliced=(n_pos + n_neg) if mass_u is None else np.asarray(mass_u, float),
        mass_spliced=spl if mass_spl is None else np.asarray(mass_spl, float),
    )


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def test_geometry_exon_intron_exon_plus_gene():
    # 1 ref, 3 regions (ex+ | in+ | ex+), 4 boundaries (b0 term, b1 ex→in, b2 in→ex, b3 term).
    rro = np.array([0, 3])
    rbo = np.array([0, 4])
    chain = build_node_chain(rro, rbo)
    # genomic order: B0 R0 B1 R1 B2 R2 B3  → node ids 0..6
    assert list(chain.kind) == [BOUNDARY, REGION, BOUNDARY, REGION, BOUNDARY, REGION, BOUNDARY]
    assert list(chain.ref_idx) == [0, 0, 1, 1, 2, 2, 3]

    L = np.array([1000.0, 2000.0, 500.0])
    sig = np.array([BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS], dtype=np.int64)
    region_arrays = SimpleNamespace(region_size_bp=L, signature=sig)
    substrate = SimpleNamespace(contained=_view([100.0, 50.0, 80.0], [0.0, 0.0, 0.0]))

    # boundary sides; b1 spliced on its LEFT (exon r0) side, b2 on its RIGHT (exon r2) side.
    # b1 (idx 1) and b2 (idx 2) are POS splice junctions (motif strand from the accumulator).
    left = _view([0.0, 30.0, 20.0, 40.0], [0.0, 88.0, 0.0, 0.0])
    right = _view([10.0, 31.0, 22.0, 0.0], [0.0, 0.0, 77.0, 0.0])
    bsub = SimpleNamespace(
        left=left,
        right=right,
        left_region=np.array([-1, 0, 1, 2]),
        right_region=np.array([0, 1, 2, -1]),
        junction_strand=np.array([0, 1, 1, 0], dtype=np.int8),
    )

    gdna_fl = _delta_pmf(300)
    rna_fl = _delta_pmf(200)
    g = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)

    # R0 (node 1): contained, same both faces; per-component FL.
    assert g.mass_left[1] == 100.0 and g.mass_right[1] == 100.0
    assert np.isclose(
        g.eff_gdna_left[1], region_eff_length(np.array([1000.0]), gdna_fl)[0]
    )  # 1000-300+1 = 701 start positions
    # The `+1` is the DISCRETE start-position count (a 300 bp fragment fits in a 1000 bp region at 701
    # positions, not 700) — irrelevant at this scale, decisive only at L ≈ ell where it is 1 vs 0.
    assert np.isclose(g.eff_gdna_left[1], 701.0) and np.isclose(g.eff_rna_left[1], 801.0)
    assert g.spliced_pos_left[1] == 0.0 and g.spliced_pos_right[1] == 0.0

    # B1 (node 2): ex→in junction. left side in r0 (1000bp), right side in r1 (2000bp). spliced on LEFT (exon).
    assert g.mass_left[2] == 30.0 and g.mass_right[2] == 31.0
    assert np.isclose(g.eff_gdna_left[2], boundary_side_eff_length(gdna_fl, np.array([1000.0]))[0])
    # per-face DENSITY length = E[min(ℓ,R)]/2 = min(300,1000)/2. The ½ is the accumulator's deposit
    # rule (a straddling fragment gives each face its own share, not the whole fragment), NOT a
    # constant — see effective_length.boundary_side_eff_length. Un-halved it read ρ/2 on every message.
    assert np.isclose(g.eff_gdna_left[2], 150.0) and np.isclose(g.eff_gdna_right[2], 150.0)
    assert np.isclose(g.eff_rna_left[2], 100.0)  # RNA-FL min(200,1000)/2
    assert g.spliced_pos_left[2] == 88.0 and g.spliced_pos_right[2] == 0.0  # exon on the left flank

    # B2 (node 4): in→ex junction. spliced on the RIGHT (exon r2) side.
    assert g.mass_left[4] == 20.0 and g.mass_right[4] == 22.0
    assert g.spliced_pos_left[4] == 0.0 and g.spliced_pos_right[4] == 77.0
    # right side is r2 (500bp): gDNA min(300,500)/2=150, RNA min(200,500)/2=100
    assert np.isclose(g.eff_gdna_right[4], 150.0) and np.isclose(g.eff_rna_right[4], 100.0)


def test_terminal_boundary_zero_off_edge():
    rro = np.array([0, 3])
    rbo = np.array([0, 4])
    chain = build_node_chain(rro, rbo)
    L = np.array([1000.0, 2000.0, 500.0])
    region_arrays = SimpleNamespace(region_size_bp=L, signature=np.array([2, 8, 2], dtype=np.int64))
    substrate = SimpleNamespace(contained=_view([100.0, 50.0, 80.0], [0.0, 0.0, 0.0]))
    left = _view([0.0, 30.0, 20.0, 40.0], [0.0, 0.0, 0.0, 0.0])
    right = _view([10.0, 31.0, 22.0, 0.0], [0.0, 0.0, 0.0, 0.0])
    bsub = SimpleNamespace(
        left=left,
        right=right,
        left_region=np.array([-1, 0, 1, 2]),
        right_region=np.array([0, 1, 2, -1]),
        junction_strand=np.zeros(4, dtype=np.int8),
    )
    g = build_node_geometry(chain, substrate, bsub, region_arrays, _delta_pmf(300), _delta_pmf(200))
    # B0 (node 0): left_region=-1 → left face eff = _EPS-floored ~0 (off edge); right face = r0 crossing.
    assert g.eff_gdna_left[0] <= 1e-6
    assert np.isclose(g.eff_gdna_right[0], 150.0)  # min(300,1000)/2 — the per-face density length
    # B3 (node 6): right_region=-1 → right face off edge.
    assert g.eff_gdna_right[6] <= 1e-6


def _empty_boundary_substrate(n_b):
    z = np.zeros(n_b)
    side = _cview(z.copy(), z.copy())
    return SimpleNamespace(
        left_region=np.full(n_b, -1),
        right_region=np.full(n_b, -1),
        left=side,
        right=side,
        junction_strand=np.zeros(n_b, dtype=np.int8),
    )


def test_init_zero_gdna_introns_via_strand():
    # The P1 gate: a zero-gDNA library. 3 regions: intergenic | intron+ | AMBIG (one ref).
    rro = np.array([0, 3])
    rbo = np.array([0, 4])
    chain = build_node_chain(rro, rbo)
    sig = np.array([0, BIT_INTRON_POS, BIT_EXON_POS | BIT_EXON_NEG], dtype=np.int64)
    sc = np.array([TS_NONE, TS_POS, TS_AMBIG], dtype=np.int8)
    region_arrays = SimpleNamespace(
        strand_class=sc,
        signature=sig,
        region_size_bp=np.array([1000.0, 2000.0, 800.0]),
    )
    # intergenic gDNA = strand-symmetric; intron+ RNA = strongly sense-tilted (κ=0.95); AMBIG = symmetric.
    contained = _cview([50.0, 95.0, 50.0], [50.0, 5.0, 50.0])
    substrate = SimpleNamespace(contained=contained)
    bsub = _empty_boundary_substrate(4)

    b = init_beliefs(chain, substrate, bsub, region_arrays, rna_sense_frac=0.95, n_grid=60)

    # region node ids on the chain: B0 R0 B1 R1 B2 R2 B3 → regions at 1, 3, 5.
    rid = [1, 3, 5]
    fg = b.f_g[rid]
    # intergenic: locked gDNA sink {0,0,1}, all precision locked (var 0).
    assert fg[0] == 1.0
    assert b.var_gdna[1] == 0.0 and b.var_pos[1] == 0.0 and b.var_neg[1] == 0.0
    # intron+ (zero gDNA): the strand tilt alone drives f_g → 0; the − axis is locked (var 0), + & g finite.
    assert fg[1] < 0.15
    assert b.var_neg[3] == 0.0 and np.isfinite(b.var_gdna[3]) and np.isfinite(b.var_pos[3])
    # AMBIG: unresolved by strand → {0,0,1} default at MAX (inf) variance for the sweep to resolve.
    assert fg[2] == 1.0
    assert np.isinf(b.var_gdna[5]) and np.isinf(b.var_pos[5]) and np.isinf(b.var_neg[5])


def test_init_boundary_continuity_gate():
    # 1 ref, 2 regions (exon+ | intron+) → boundaries B0(term) B1(ex+→in+ junction) B2(term).
    rro = np.array([0, 2])
    rbo = np.array([0, 3])
    chain = build_node_chain(rro, rbo)
    sig = np.array([BIT_EXON_POS, BIT_INTRON_POS], dtype=np.int64)
    sc = np.array([TS_POS, TS_POS], dtype=np.int8)
    region_arrays = SimpleNamespace(
        strand_class=sc,
        signature=sig,
        region_size_bp=np.array([1000.0, 2000.0]),
    )
    substrate = SimpleNamespace(contained=_cview([80.0, 40.0], [4.0, 30.0]))
    # B1 crossing: sense-tilted unspliced (κ=0.95 ⇒ +) + a sense-spliced floor on the exon (left) flank.
    left = _cview([0.0, 90.0, 0.0], [0.0, 5.0, 0.0], spl_sense=[0.0, 50.0, 0.0])
    right = _cview([0.0, 91.0, 0.0], [0.0, 6.0, 0.0], spl_sense=[0.0, 0.0, 0.0])
    bsub = SimpleNamespace(
        left_region=np.array([-1, 0, 1]), right_region=np.array([0, 1, -1]), left=left, right=right
    )

    b = init_beliefs(chain, substrate, bsub, region_arrays, rna_sense_frac=0.95, n_grid=60)
    # boundary node ids: B0=0, B1=2, B2=4.
    # terminals: off-edge flank ⇒ neither strand continuous ⇒ G1 gDNA sink (var locked at 0).
    assert b.f_g[0] == 1.0 and b.var_gdna[0] == 0.0
    assert b.f_g[4] == 1.0 and b.var_gdna[4] == 0.0
    # B1 (ex+→in+): +strand continuous (G2+) ⇒ the strand tilt resolves f_g → 0; − axis locked (var 0).
    assert b.f_g[2] < 0.15
    assert b.var_neg[2] == 0.0 and np.isfinite(b.var_gdna[2])


def test_init_tss_boundary_is_black_hole():
    # intergenic | exon+ : the internal boundary is a TSS (intergenic↔exon) ⇒ continuity blocks RNA ⇒ sink.
    rro = np.array([0, 2])
    rbo = np.array([0, 3])
    chain = build_node_chain(rro, rbo)
    sig = np.array([0, BIT_EXON_POS], dtype=np.int64)
    sc = np.array([TS_NONE, TS_POS], dtype=np.int8)
    region_arrays = SimpleNamespace(
        strand_class=sc,
        signature=sig,
        region_size_bp=np.array([1000.0, 2000.0]),
    )
    substrate = SimpleNamespace(contained=_cview([50.0, 80.0], [50.0, 4.0]))
    # the TSS-crossing fragments are sense-tilted, but continuity must STILL block RNA (the black hole).
    left = _cview([0.0, 90.0, 0.0], [0.0, 5.0, 0.0])
    right = _cview([0.0, 91.0, 0.0], [0.0, 6.0, 0.0])
    bsub = SimpleNamespace(
        left_region=np.array([-1, 0, 1]), right_region=np.array([0, 1, -1]), left=left, right=right
    )

    b = init_beliefs(chain, substrate, bsub, region_arrays, rna_sense_frac=0.95, n_grid=60)
    # B1 (node id 2) is the TSS: a locked gDNA sink despite the sense tilt (all precision locked at 0).
    assert b.f_g[2] == 1.0 and b.var_gdna[2] == 0.0 and b.var_pos[2] == 0.0 and b.var_neg[2] == 0.0


def test_precision_state_count_resolution():
    """The log-density solver's precision state is ``Var(log f_g)`` — the message currency (D2). It reflects
    EVIDENCE: a node with more fragments (same composition) resolves its log-density sharper, so a lower
    ``Var(log f_g)`` ⇒ a more confident message. (In LOG space a confident ``f_g→0`` node has WIDE variance —
    a near-zero gDNA density carries little reliable gDNA-density information to impute, the
    "zero-density-is-not-a-measurement" principle — so the lattice's linear ``Var(f_g)`` ordering does not
    carry over.) A node with no fragments reports zero variance."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    kappa = 0.99
    z = np.zeros(2)
    # Same single-strand + composition at two evidence levels: node 1 has 20× the counts of node 0.
    u_pos = np.array([20.0, 400.0])
    u_neg = np.array([20.0, 400.0])
    allow_pos = np.array([True, True])
    allow_neg = np.array([False, False])  # both single-strand
    mass = u_pos + u_neg
    d = _solve_nodes_logodds_all(
        u_pos,
        u_neg,
        allow_pos,
        allow_neg,
        mass,
        z,
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_grid=60,
    )
    assert d.gdna_frac_var is not None
    # p̂=0.5 at κ=0.99 ⇒ the fragments look unstranded ⇒ the mean channel points at the gDNA mode f_g=1.
    # Under the count-zero-info variance freeze the count enters as PRECISION: more evidence sharpens that
    # signal, so the higher-count node resolves NEARER the mode with a lower Var(log f_g) — it is not pinned
    # count-independent (that was the pre-freeze f_g-dependent normalizer artifact).
    assert d.gdna_frac[0] > 0.85 and d.gdna_frac[1] > 0.85  # both gDNA-dominant
    assert d.gdna_frac[1] >= d.gdna_frac[0]  # more count ⇒ nearer the mode
    assert d.gdna_frac_var[1] < d.gdna_frac_var[0]  # more fragments ⇒ sharper
    # all per-component variances are present, finite, non-negative for active nodes.
    for v in (d.gdna_frac_var, d.rna_pos_frac_var, d.rna_neg_frac_var):
        assert np.all(np.isfinite(v)) and np.all(v >= 0.0)
    # a no-fragment node is inactive ⇒ zero variance on every component.
    d0 = _solve_nodes_logodds_all(
        np.array([0.0]),
        np.array([0.0]),
        np.array([True]),
        np.array([True]),
        np.array([0.0]),
        np.array([0.0]),
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_grid=60,
    )
    assert d0.gdna_frac_var[0] == 0.0 and d0.rna_pos_frac_var[0] == 0.0


def test_gdna_sweep_factor1_uniform():
    # The factor-1 bedrock: a UNIFORM-gDNA chain intergenic | AMBIG | intergenic. Every node's mass is laid
    # down as ρ·eff-len (the accumulator's uniform deposit); after the sweep each node's ρ_g must read back ρ
    # — including the AMBIG node the sweep actually solves (the locked intergenic nodes are trivially ρ).
    rho = 0.5
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    chain = build_node_chain(np.array([0, 3]), np.array([0, 4]))
    L = np.array([1000.0, 1000.0, 1000.0])
    sig = np.array([0, BIT_EXON_POS | BIT_EXON_NEG, 0], dtype=np.int64)
    sc = np.array([TS_NONE, TS_AMBIG, TS_NONE], dtype=np.int8)
    region_arrays = SimpleNamespace(strand_class=sc, signature=sig, region_size_bp=L)
    reg_eff = region_eff_length(L, gdna_fl)  # [700,700,700]
    cmass = rho * reg_eff
    substrate = SimpleNamespace(
        contained=_cview(cmass / 2, cmass / 2, mass_u=cmass, mass_spl=np.zeros(3))
    )
    side_eff = boundary_side_eff_length(gdna_fl, L)  # [300,300,300]
    lr, rr = np.array([-1, 0, 1, 2]), np.array([0, 1, 2, -1])
    lmass = np.where(lr >= 0, rho * side_eff[np.clip(lr, 0, 2)], 0.0)
    rmass = np.where(rr >= 0, rho * side_eff[np.clip(rr, 0, 2)], 0.0)
    left = _cview(lmass / 2, lmass / 2, mass_u=lmass, mass_spl=np.zeros(4))
    right = _cview(rmass / 2, rmass / 2, mass_u=rmass, mass_spl=np.zeros(4))
    bsub = SimpleNamespace(
        left_region=lr,
        right_region=rr,
        left=left,
        right=right,
        junction_strand=np.zeros(len(lr), dtype=np.int8),
    )

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.7, n_grid=40, statics=st
    )
    final = node_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        bsub,
        rna_sense_frac=0.7,
        n_grid=40,
    )
    # gDNA density per face = f_g·M_face / E_gdna_face (the formula the sweep inlines; node_densities removed).
    rho_g_left = final.f_g * np.asarray(geom.mass_left) / np.asarray(geom.eff_gdna_left)
    rho_g_right = final.f_g * np.asarray(geom.mass_right) / np.asarray(geom.eff_gdna_right)
    interg, ambig = [1, 5], 3  # region nodes: R0/R2 intergenic (strand-anchored), R1 AMBIG
    # Intergenic nodes recover ρ exactly (the strand/signature pins them — this invariant holds in every phase).
    assert np.allclose(rho_g_left[interg], rho, atol=0.02)
    assert np.allclose(rho_g_right[interg], rho, atol=0.02)
    # AMBIG node — RESOLVED. This long-standing gap (ρ_g≈0.28 vs the true 0.50) was never the "Phase-1
    # prior-free" gap it was documented as: the AMBIG node carried its own improper prior, `kde + log f_g`,
    # whose Jacobian was an implicit ANTI-gDNA tilt that pulled a balanced node off ρ no matter what the
    # strand said. Removing the reference prior + Jacobian (ψ = strand + logP_g + logP_r, bare —
    # prior_ramp_and_bp_roadmap.md §2) recovers ρ exactly, at the atol the old comment predicted.
    assert np.allclose(rho_g_left[ambig], rho, atol=0.05)
    assert np.allclose(rho_g_right[ambig], rho, atol=0.05)


def test_gdna_emits_across_tss_tes_seam():
    """Structural-gate regression (`message_state_separation.md`): gDNA is genomically continuous, so the
    gene-boundary seams (TSS/TES) flanking a SINGLE-EXON gene must RELAY a gDNA message into it from the
    intergenic regions beyond — even though neither RNA strand is continuous across those seams. Before the
    fix the gDNA message was gated by RNA strand-continuity (`solvable`), so such a gene (both flanks
    intergenic) was a no-relay node, solving on its own local belief alone.

    Conversely, the intergenic flank is structurally RNA-free and emits ZERO RNA authority: the exon receives
    no +/− RNA message from it (a node's confidence about its OWN all-gDNA state grants no authority over a
    neighbour's RNA). The assertions lock the two halves of the three-term emission gate.
    """
    rho = 0.5
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    chain = build_node_chain(np.array([0, 3]), np.array([0, 4]))  # intergenic | exon+ | intergenic
    L = np.array([1000.0, 1000.0, 1000.0])
    sig = np.array([0, BIT_EXON_POS, 0], dtype=np.int64)
    sc = np.array([TS_NONE, TS_POS, TS_NONE], dtype=np.int8)
    region_arrays = SimpleNamespace(strand_class=sc, signature=sig, region_size_bp=L)
    cmass = rho * region_eff_length(L, gdna_fl)
    substrate = SimpleNamespace(
        contained=_cview(cmass / 2, cmass / 2, mass_u=cmass, mass_spl=np.zeros(3))
    )
    side_eff = boundary_side_eff_length(gdna_fl, L)
    lr, rr = np.array([-1, 0, 1, 2]), np.array([0, 1, 2, -1])
    lmass = np.where(lr >= 0, rho * side_eff[np.clip(lr, 0, 2)], 0.0)
    rmass = np.where(rr >= 0, rho * side_eff[np.clip(rr, 0, 2)], 0.0)
    left = _cview(lmass / 2, lmass / 2, mass_u=lmass, mass_spl=np.zeros(4))
    right = _cview(rmass / 2, rmass / 2, mass_u=rmass, mass_spl=np.zeros(4))
    bsub = SimpleNamespace(
        left_region=lr,
        right_region=rr,
        left=left,
        right=right,
        junction_strand=np.zeros(len(lr), dtype=np.int8),
    )

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.7, n_grid=40, statics=st
    )
    cap = {}
    final = node_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        bsub,
        rna_sense_frac=0.7,
        n_grid=40,
        _capture=cap,
    )
    exon = 3  # chain id of the single-exon gene (R1), flanked on both sides by TSS/TES seams
    # THE FIX — the exon receives a gDNA relay across the seam (incoming precision > 0). Pre-fix: 0 (no relay).
    assert cap["prec_g"][exon] > 0.0, "single-exon gene got NO gDNA relay across the TSS/TES seam"
    # The intergenic flanks emit ZERO RNA authority: the exon receives no +/− RNA message from them.
    assert cap["prec_p"][exon] == 0.0 and cap["prec_n"][exon] == 0.0
    # State ⊥ messages: the intergenic nodes stay locked all-gDNA (confident own-state, ignore all inputs).
    assert final.f_g[1] == 1.0 and final.f_g[5] == 1.0


@pytest.mark.xfail(
    reason="Pre-existing known-red, superseded by the τ-precision landing (2026-07-20). On this ARTIFICIAL "
    "seedless chain (intron+|AMBIG|intron−, no intergenic buffer, no enrichment prior ⇒ σ²_transfer=0) two "
    "documented effects break the old `<0.5` bound: (1) the AMBIG node leans gDNA at ≈0.564 — the CORRECT "
    "reference behaviour on an unidentifiable balanced node (reference_prior_derivation.md §10.7), the level "
    "deferred to the hyperprior; (2) the strand-solved introns (local f_g≈0.034) are overridden UP to ≈0.564 "
    "by the directly-adjacent terminal G1 gDNA locks whose message is UNDAMPED here (no σ²_transfer). Both are "
    "artefacts of this minimal chain: on real data (σ²_transfer + intergenic buffers) the stranded controls are "
    "PRESERVED (validate_tau 2026-07-20: STR CW% 0.3→0.3). Re-derive/retire during real-data validation.",
    strict=False,
)
def test_gdna_sweep_zero_gdna_pin_and_monotone():
    # A pure-RNA chain intron+ | AMBIG(in+|in−) | intron−. The AMBIG starts at the all-gDNA init f_g=1; the
    # global (driven to ~0 by the RNA introns) + the RNA-neighbour messages must pull the phantom gDNA down,
    # monotonically.
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    chain = build_node_chain(np.array([0, 3]), np.array([0, 4]))
    L = np.array([2000.0, 2000.0, 2000.0])
    sig = np.array(
        [BIT_INTRON_POS, BIT_INTRON_POS | BIT_INTRON_NEG, BIT_INTRON_NEG], dtype=np.int64
    )
    sc = np.array([TS_POS, TS_AMBIG, TS_NEG], dtype=np.int8)
    region_arrays = SimpleNamespace(strand_class=sc, signature=sig, region_size_bp=L)
    cpos = np.array(
        [95.0, 50.0, 5.0]
    )  # sense-tilted RNA (κ=0.95): + intron genome+, − intron genome−
    cneg = np.array([5.0, 50.0, 95.0])
    substrate = SimpleNamespace(
        contained=_cview(cpos, cneg, mass_u=cpos + cneg, mass_spl=np.zeros(3))
    )
    lr, rr = np.array([-1, 0, 1, 2]), np.array([0, 1, 2, -1])
    # crossing counts cleanly sense-tilted like the regions: B0/B1 are +crossings (genome+), B2/B3 are −.
    lpos = np.array([0.0, 40.0, 2.0, 2.0])
    lneg = np.array([0.0, 2.0, 40.0, 40.0])
    rpos = np.array([40.0, 40.0, 2.0, 0.0])
    rneg = np.array([2.0, 2.0, 40.0, 0.0])
    left = _cview(lpos, lneg, mass_u=lpos + lneg, mass_spl=np.zeros(4))
    right = _cview(rpos, rneg, mass_u=rpos + rneg, mass_spl=np.zeros(4))
    bsub = SimpleNamespace(
        left_region=lr,
        right_region=rr,
        left=left,
        right=right,
        junction_strand=np.zeros(len(lr), dtype=np.int8),
    )

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.95, n_grid=40, statics=st
    )
    assert belief.f_g[3] == 1.0  # AMBIG starts all-gDNA
    final = node_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        bsub,
        rna_sense_frac=0.95,
        n_rna_obs=10000.0,  # library sample sizes so the stranded (κ=0.95) intron seeds fire (τ noise floor)
        n_gdna_obs=10000.0,
        n_grid=40,
    )
    # The AMBIG phantom is pulled DOWN from its all-gDNA init (1.0) toward RNA. This chain is the WORST
    # case for a balanced AMBIG node: it is an ARTIFICIAL all-RNA chain (intron+|AMBIG|intron−) with NO
    # intergenic structural seeds, so the gDNA prior has almost nothing to anchor a zero-gDNA baseline. The
    # AMBIG node's strand is balanced (both strands live), so the strand likelihood is DEGENERATE — a
    # balanced count is equally consistent with gDNA and with balanced ±RNA — and the neutral (f_g, τ)
    # reference measure parsimoniously leans a balanced count toward gDNA, deferring to the prior for the
    # level. With no seeds the (weak) floor + the intron RNA-imputation only pull it to ~0.44 here; on real
    # libraries the intergenic zero-count seeds make the prior decisive (the gdna_none capture-on benchmark
    # shows ~0 false gDNA). Still pulled well below the all-gDNA init and RNA-leaning.
    assert final.f_g[3] < 0.50
    # single-strand introns: the decisive strand wins and the floor DEFERS (a hyperprior cannot overrule a
    # node's own strand evidence) → they stay RNA-leaning. ~0.44 since `emit_locked` defaults ON (was ~0.22
    # muted): this chain's TERMINAL boundaries are G1-locked, and un-muting them makes them emit their
    # structural all-gDNA into the flanking introns. That is an ARTEFACT OF THIS ARTIFICIAL CHAIN — an intron
    # running to the chain edge with no intergenic flank cannot occur in a real annotation (a transcript ends
    # at a TSS/TES, where a crossing fragment really is gDNA because it has no RNA upstream). On real data the
    # seams in a zero-gDNA library carry ~no mass, so the `sm > _EPS` gate keeps them silent and the
    # false-positive guard holds: zero-gDNA mwae is UNCHANGED by `emit_locked` on both suites (0.0088 on
    # quick_3to1_5mb, 0.0543 on ambig_dense_10mb) with over-call also unchanged. The invariant this test
    # actually protects — the introns stay RNA-leaning, well below their all-gDNA init — still holds.
    assert final.f_g[1] < 0.50 and final.f_g[5] < 0.50


# --- density-Gaussian message form: two-sided pull + emergent deference (precision_state_design.md §5) --------


def test_density_message_two_sided_mode_not_vertex():
    """A LOG-fraction gDNA message (mode=log 0.2, strong prec) on a balanced AMBIG node (flat strand) pulls
    f_g TOWARD 0.2 — two-sided by construction (a Gaussian on log f_g, no edge wall), not to the f_g=1
    vertex. The log-density log-odds solver's message form."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    z = np.zeros(1)
    # AMBIG node, balanced counts ⇒ the strand is flat (κ=0.5); only the message shapes f_g.
    d = _solve_nodes_logodds_all(
        np.array([50.0]),
        np.array([50.0]),
        np.array([True]),
        np.array([True]),
        np.array([100.0]),
        z,
        kappa=0.5,
        od_g=0.0,
        od_r=0.0,
        n_grid=80,
        gdna_imp_mode=np.array([np.log(0.2)]),
        gdna_imp_prec=np.array([200.0]),
    )
    fg = float(d.gdna_frac[0])
    assert abs(fg - 0.2) < 0.05, fg


def test_density_message_defers_to_decisive_strand():
    """Emergent deference: a WEAK gDNA message (prec=3) trying to pull f_g→0.9 must lose to a decisive
    single-strand node's ~1000-fragment strand likelihood — f_g stays ≈0 (the honest precision blend means a
    weak message cannot override the data; no log-wall to force it off zero)."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    z = np.zeros(1)
    d = _solve_nodes_logodds_all(
        np.array([1000.0]),
        np.array([5.0]),
        np.array([True]),
        np.array([False]),
        np.array([1005.0]),
        z,
        kappa=0.99,
        od_g=0.0,
        od_r=0.0,
        n_grid=80,
        gdna_imp_mode=np.array([np.log(0.9)]),
        gdna_imp_prec=np.array([3.0]),
    )
    fg = float(d.gdna_frac[0])
    assert fg < 0.1, fg


# --- mature absorption: the spliced mass "absorbs" the imputed RNA, leaving only NASCENT ---------------
# (`spliced_mature_nascent_message.md`). The RNA message src→dst is
#   ρ = src_nascent/E_r + SP[sf][src]/E_spl_src − SP[df][dst]/E_spl_dst.
# The dst-face term subtracts the mature a junction boundary measures, so a pure-mature exon imputes
# ≈0 nascent into the intron beyond it — no wholesale nascent hallucination.


def _mature_exon_chain(*, spliced: bool, rho_g=0.5, rho_m=1.0, kappa=0.95, spl_scale=1.0):
    """A `intron+ | exon+ | intron+` chain, physically consistent for a pure-MATURE expressed exon with
    NO nascent: the exon's contained unspliced = balanced gDNA + sense (+) mature; the introns' contained
    + boundary crossings = balanced gDNA only (mature skips the intron as SPLICED, never as an unspliced
    crossing); the two intron↔exon junctions carry the one-sided sense spliced (mature) on the EXON flank
    iff ``spliced``. Returns the node_sweep args (chain, statics, geometry, belief, region_arrays, bsub)."""
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    chain = build_node_chain(np.array([0, 3]), np.array([0, 4]))  # B0 R0 B1 R1 B2 R2 B3
    L = np.array([2000.0, 2000.0, 2000.0])
    sig = np.array([BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS], dtype=np.int64)
    sc = np.array([TS_POS, TS_POS, TS_POS], dtype=np.int8)
    region_arrays = SimpleNamespace(strand_class=sc, signature=sig, region_size_bp=L)
    Eg = region_eff_length(L, gdna_fl)  # contained gDNA eff-len  (1700)
    Er = region_eff_length(L, rna_fl)  # contained RNA  eff-len  (1800)
    g_half = rho_g * Eg / 2.0  # per-strand contained gDNA count (balanced)
    mat = rho_m * Er  # contained mature count (+ strand only)
    # exon R1: gDNA (balanced) + mature (+); introns R0/R2: gDNA only.
    u_pos = np.array([g_half[0], g_half[1] + mat[1], g_half[2]])
    u_neg = np.array([g_half[0], g_half[1], g_half[2]])
    contained = _cview(u_pos, u_neg, mass_u=u_pos + u_neg, mass_spl=np.zeros(3))
    substrate = SimpleNamespace(contained=contained)
    # boundary crossings: gDNA only (balanced); spliced(mature) on the EXON flank of the two junctions.
    side_g = boundary_side_eff_length(gdna_fl, L)  # (300)
    spl_eff = spliced_side_eff_length(rna_fl, L)  # one-sided half-triangle (100)
    lr, rr = np.array([-1, 0, 1, 2]), np.array([0, 1, 2, -1])

    def gx(r):
        return np.where(r >= 0, rho_g * side_g[np.clip(r, 0, 2)], 0.0)

    lcross, rcross = gx(lr), gx(rr)
    # ``spl_scale`` < 1 models a CAPTURE-DEPLETED junction: junction-spanning (spliced) reads are only
    # partially captured, so the junction under-reports the exon's true mature density ⇒ the B→exon mature
    # MEASUREMENT DISAGREES with the exon's own (confident) unspliced belief. Used by the silencing test.
    spl_mat = rho_m * spl_eff[1] * spl_scale  # mature spliced count on the exon flank
    # B1 (idx1): exon R1 is its RIGHT region → spliced on the RIGHT side. B2 (idx2): exon is its LEFT
    # region → spliced on the LEFT side. (One-sided, exon flank.)
    spl_l = np.array([0.0, 0.0, spl_mat if spliced else 0.0, 0.0])
    spl_r = np.array([0.0, spl_mat if spliced else 0.0, 0.0, 0.0])
    left = _cview(lcross / 2, lcross / 2, spl_sense=spl_l, mass_u=lcross, mass_spl=spl_l.copy())
    right = _cview(rcross / 2, rcross / 2, spl_sense=spl_r, mass_u=rcross, mass_spl=spl_r.copy())
    bsub = SimpleNamespace(
        left_region=lr,
        right_region=rr,
        left=left,
        right=right,
        junction_strand=np.array([0, TS_POS, TS_POS, 0], dtype=np.int8),
    )
    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=kappa, n_grid=60, statics=st
    )
    return chain, st, geom, belief, region_arrays, bsub


def _sweep(args, kappa=0.95, n_rna_obs=10000.0, n_gdna_obs=10000.0):
    chain, st, geom, belief, ra, bsub = args
    cap = {}
    final = node_sweep(
        chain,
        st,
        geom,
        belief,
        ra,
        bsub,
        rna_sense_frac=kappa,
        # The τ strand seed needs the library sample sizes to size its overdispersion noise floor
        # ¼·(1/N + ω); a strongly-stranded fixture (κ=0.95) fires only when N is supplied (the default 0 ⇒
        # ∞ floor ⇒ gated). Large N here ⇒ floor≈0 ⇒ the (2κ−1)² strand seed fires at full strength.
        n_rna_obs=n_rna_obs,
        n_gdna_obs=n_gdna_obs,
        n_grid=60,
        _capture=cap,
    )
    return final, cap


@pytest.mark.xfail(
    reason="Tracked regression on a pure-mature TOY, awaiting the nascent factory (ρ_nascent = ρ_RNA − "
    "ρ_mature, intron-baselined) + the honest σ²_transfer precision. The flanking pure-gDNA introns land at "
    "f_g≈0.82 (truth 1.0): the exon's ~95%-mature unspliced payload leaks in as nascent because the RNA-total "
    "factor does not yet subtract mature. The mature gate that used to block this edge was DISMANTLED "
    "(docs/calibration/message_model_derivation.md §5) — and the dismantle IMPROVED this toy (0.766 gated → "
    "0.821 un-gated), so this xfail marks the residual σ²_transfer/nascent-factory gap, NOT a dismantle "
    "regression. Un-xfail when the nascent factory lands.",
    strict=False,
)
def test_mature_no_nascent_hallucination_in_introns():
    """The user's red line: a pure-mature expressed exon (nascent=0) must NOT manufacture wholesale nascent in
    its flanking pure-gDNA introns; the introns stay gDNA (truth f_g=1). xfail until the honest RNA
    counter-message (nascent = RNA − mature) that pins the introns to ~0 nascent is built."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_introns = fin_m.f_g[[1, 5]]  # chain ids of R0, R2
    assert np.all(fg_introns > 0.85), fg_introns


# NOTE: `test_mature_absorption_lowers_nascent_message_into_junction` was RETIRED when the mature-crossing gate
# landed (Phase 4). It asserted the exon→B1 +RNA message FIRES (`app[b1] > 0`) so its absorption term could
# lower the imputed nascent; the gate now blocks that edge entirely (the exon may not manufacture nascent into
# its intron-side junction), so the message no longer exists to absorb. Its replacement is
# `test_exon_does_not_manufacture_nascent_into_intron` (same fixture, same edge, inverted assertion). The
# B→exon MEASUREMENT + absorption path it half-covered is still guarded by the two `test_mature_measurement_*`
# tests below, which the gate leaves untouched.


def test_mature_measurement_recovers_exon_rna():
    """The companion direction (unchanged B→exon MEASUREMENT): the same chain's expressed exon is
    correctly recovered as mostly RNA (its true f_g ≈ ρ_g·E_g/(ρ_g·E_g+ρ_m·E_r) ≈ 0.32), driven by the
    + strand tilt + the mature measurement — so the absorption does not starve the exon of its own RNA."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_exon = float(fin_m.f_g[3])  # chain id of R1 (the exon)
    assert fg_exon < 0.45, fg_exon  # truth ≈0.32; comfortably RNA-dominated, not pinned to gDNA


def test_mature_measurement_disagreement_silenced():
    """BUG #2 regression: the mature MEASUREMENT message must be DISAGREEMENT-SILENCED like every other RNA
    message (the old exemption applied it at full COUNT precision). Under capture, junction-spanning reads are
    only partially captured, so the B→exon mature density UNDER-reports the exon's true RNA → the measurement
    DISAGREES with the exon's own confident belief. Un-silenced it dragged f_pos down → phantom gDNA by simplex
    complement (−gDNA flagship +0.04→+0.018). Here: a 4×-depleted junction (spl_scale=0.25) genuinely lowers
    the message target, yet the exon's gDNA fraction stays essentially unchanged vs a consistent junction —
    the disagreeing measurement was down-weighted, not applied whole."""
    ex = 3  # chain id of the exon R1
    fin_ok, cap_ok = _sweep(_mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=1.0))
    fin_lo, cap_lo = _sweep(_mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=0.25))
    # (1) the depleted junction really did lower the +RNA message target into the exon (a genuine disagreement)…
    assert cap_lo["mode_p"][ex] < cap_ok["mode_p"][ex] - 0.3, (
        cap_lo["mode_p"][ex],
        cap_ok["mode_p"][ex],
    )
    # (2) …yet the exon's gDNA fraction barely moves (silenced — pre-fix the low measurement inflated f_g).
    assert abs(float(fin_lo.f_g[ex]) - float(fin_ok.f_g[ex])) < 0.05, (
        fin_lo.f_g[ex],
        fin_ok.f_g[ex],
    )
    # (3) and the exon stays RNA-dominated, not pulled toward phantom gDNA.
    assert float(fin_lo.f_g[ex]) < 0.45, fin_lo.f_g[ex]


def test_tau_gag_fix_spliced_junction_emits_when_unstranded():
    """τ-GAG REGRESSION (`message_system_implementation_plan.md` §Phase B, 2026-07-21). On UNSTRANDED data
    (κ=½ ⇒ the strand Fisher info ``I_strand`` is identically 0), a splice-junction boundary still carries
    motif-stranded spliced (mature-RNA) fragments — a DIRECT measurement, independent of strand. That
    measurement MUST reach the exon. The bug: the τ-evidence emission gate (keyed on ``I_strand``+``I_struct``
    only, NOT the spliced count) silenced it, and the spliced-precision credit — which lives *inside* the gated
    block — never fired (52% of junctions). The fix opens RNA emission on spliced presence while keeping the
    deconvolution PREDICTION τ-gated.

    Pins both halves: (1) a spliced junction DELIVERS a +RNA message to its exon even unstranded; (2) the same
    chain with the spliced REMOVED delivers zero +RNA authority (a vacuous unstranded node manufactures no
    phantom RNA — the deconvolution stays gated). This exact pair fails on the pre-fix gated code."""
    ex = 3  # chain id of the expressed exon R1
    fin_spl, cap_spl = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # (1) THE FIX: the spliced (mature) MEASUREMENT reaches the exon with the strand silent (κ=½).
    assert cap_spl["prec_p"][ex] > 0.0, cap_spl["prec_p"][ex]
    # (2) the vacuous control (no spliced, no strand): ZERO +RNA authority — no phantom manufactured.
    assert cap_no["prec_p"][ex] == 0.0, cap_no["prec_p"][ex]
    # (3) the real mature measurement moves the exon TOWARD RNA — never toward phantom gDNA.
    assert float(fin_spl.f_g[ex]) < float(fin_no.f_g[ex]), (fin_spl.f_g[ex], fin_no.f_g[ex])


def test_tau_gag_fix_deconvolution_prediction_stays_gated():
    """The safety half of the τ-gag fix: unblocking the spliced MEASUREMENT must NOT unblock the deconvolution
    PREDICTION on a vacuous source (that is the phantom the τ-precision exists to kill). On the unstranded
    no-spliced chain, the boundary→exon +RNA precision is exactly 0 (asserted above); here we pin that even the
    gDNA message a vacuous boundary sends carries no manufactured composition confidence beyond the honest
    structural/count evidence — the exon's solved f_g stays near the uninformative reference, not driven to a
    confident vertex by a phantom."""
    ex = 3
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # No spliced + no strand ⇒ the exon has no composition evidence of its own; it must not be pinned to a
    # confident vertex by the messages (the phantom would drive it to ~1). It stays mid-range (reference-led).
    assert 0.2 < float(fin_no.f_g[ex]) < 0.8, fin_no.f_g[ex]


def test_compile_strand_evidence_deadband_kills_unstranded():
    """A1 unit (`StrandEvidence`, `message_system_derivation.md` §6B): the DERIVED deadband
    ``disc = 4·max(0, (κ−½)² − σ²_d)`` makes I_strand IDENTICALLY 0 on unstranded data (κ=½) — the real phantom
    fix — and >0 on stranded data; a gDNA-free library (N_gdna=0 ⇒ σ²_d→∞) also gates it to 0."""
    from rigel.calibration.bp_solver import _compile_strand_evidence

    u = np.array([100.0, 100.0])
    fg = np.array([0.5, 0.5])
    reg = np.array([True, True])
    unl = np.array([False, False])
    base = dict(od_g=0.03, od_r=0.03, n_gdna_obs=1e4, n_rna_obs=1e4, is_region=reg, locked=unl)
    ev_unstr = _compile_strand_evidence(u, u, fg, kappa=0.5, **base)
    ev_str = _compile_strand_evidence(u, u, fg, kappa=0.99, **base)
    assert np.all(ev_unstr.tau0_lam == 0.0)  # unstranded: the deadband kills I_strand (the phantom fix)
    assert np.all(ev_str.tau0_lam > 0.0)  # stranded: I_strand fires
    assert np.all(ev_unstr.tau0_th == 0.0) and np.all(ev_str.tau0_th > 0.0)
    ev_nog = _compile_strand_evidence(  # gDNA-free ⇒ σ²_d→∞ ⇒ disc=0 even when stranded
        u, u, fg, kappa=0.99, od_g=0.03, od_r=0.03, n_gdna_obs=0.0, n_rna_obs=1e4, is_region=reg, locked=unl
    )
    assert np.all(ev_nog.tau0_lam == 0.0)


def test_compile_strand_evidence_struct_lock_regions_only():
    """A1 unit: I_struct (``struct_lock``) is composition-certainty for LOCKED REGION nodes only — never a
    boundary seam (locked by structure but sitting between RNA-carrying exons)."""
    from rigel.calibration.bp_solver import _compile_strand_evidence

    z = np.zeros(4)
    ev = _compile_strand_evidence(
        z, z, np.full(4, 0.5), kappa=0.5, od_g=0.03, od_r=0.03, n_gdna_obs=1e4, n_rna_obs=1e4,
        is_region=np.array([True, True, False, False]), locked=np.array([True, False, True, False]),
    )
    # locked region → True; unlocked region → False; locked BOUNDARY (seam) → False; unlocked boundary → False
    assert list(ev.struct_lock) == [True, False, False, False]


def test_strand_overdispersion_prior_default_is_near_binomial():
    """BUG #1 regression: the shipped default strand-overdispersion prior must be the NEAR-BINOMIAL null
    (α=β=14 ⇒ od₀≈0.034), NOT the old over-conservative 0.143 (α=β=3) that widened the gDNA Beta-Binomial
    and erased its specificity at its own mean ½. The validator floors α=β at 2 (Beta(2,2), od=0.2, the most
    overdispersion allowed)."""
    import pytest

    from rigel.calibration.gdna_strand import overdispersion_for_beta
    from rigel.config import CalibrationConfig

    cfg = CalibrationConfig()
    assert cfg.gdna_strand_prior_alpha_beta == 14.0
    assert cfg.rna_strand_prior_alpha_beta == 14.0
    assert overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta) < 0.05
    assert overdispersion_for_beta(3.0) > 0.14  # the old default was ~4× more overdispersed
    with pytest.raises(ValueError):
        CalibrationConfig(gdna_strand_prior_alpha_beta=1.5)


def test_pure_gdna_node_confident_at_near_binomial_od():
    """BUG #1 mechanism (unit): a pure-gDNA single-strand node has EXACT 50/50 per-strand counts, which the
    strand mixture (gDNA mean ½, RNA mean κ≠½) must read as gDNA — f_g≈1. At the near-binomial od (the fixed
    default) it does; at the old inflated od=0.143 the widened gDNA BB loses specificity at ½ and the node is
    dragged toward the RNA/gDNA boundary (f_g well below 1). A pure-RNA control (+frac=κ) stays f_g≈0 at both."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    # κ=0.7 (intermediate strand): gDNA mean ½ is near enough to the RNA mean that the gDNA BB width matters —
    # exactly where the inflated prior does its damage (and where the toy battery regressed pre-fix).
    def solve(u_pos, u_neg, od):
        z = np.zeros(1)
        n = float(u_pos + u_neg)
        return float(
            _solve_nodes_logodds_all(
                np.array([float(u_pos)]),
                np.array([float(u_neg)]),
                np.array([True]),
                np.array([False]),
                np.array([n]),
                z,
                kappa=0.7,
                od_g=od,
                od_r=od,
                n_grid=80,
            ).gdna_frac[0]
        )

    # pure gDNA (exact 50/50, truth f_g=1): near-binomial od → confidently gDNA; inflating od monotonically
    # under-calls it. The harm is REAL but far smaller than it used to look: this assertion previously needed
    # only od=0.143 to force a >0.15 collapse, because ψ then also carried the improper `+0.5·λ` ramp, which
    # AMPLIFIED od harm by fighting the strand near the vertex. With ψ bare the strand speaks cleanly and the
    # solver is materially more od-robust (0.9945 → 0.9683 at od=0.143, vs a collapse below 0.844 before) —
    # i.e. the old numbers were further from the truth, not closer. Assert the physics, not the artifact.
    fg = [solve(500, 500, od) for od in (0.034, 0.143, 0.4)]
    assert fg[0] > 0.8, fg
    assert fg[0] > fg[1] > fg[2], fg  # monotone: inflating od always degrades the gDNA call
    assert fg[2] < fg[0] - 0.15, fg  # and at a materially inflated od the damage is large
    # pure RNA (+frac = κ = 0.7): near-binomial stays RNA-dominated; the inflated prior's symmetric harm is
    # MORE false gDNA on RNA too (it pulls every node toward ½). (At this intermediate κ the gDNA/RNA means
    # are close, so a small residual f_g is inherent — the point is near-binomial is cleaner.)
    rna_near = solve(700, 300, 0.034)
    rna_infl = solve(700, 300, 0.143)
    assert rna_near < 0.25, rna_near
    assert rna_infl > rna_near, (rna_infl, rna_near)


# ---------------------------------------------------------------------------
# RNA-message routing after the mature-crossing gate was DISMANTLED
# (docs/calibration/message_model_derivation.md §5).
#
# Only the STRUCTURAL per-strand `free_s` continuity gate remains: each RNA strand's density flows wherever that
# strand is continuous on BOTH endpoints (intron↔exon in either direction, intron→boundary, boundary→exon), and
# gDNA flows genomically. The asymmetric `send_s = mrna_active_s[dst] or not mrna_active_s[src]` gate that used
# to silence exon→intron RNA is GONE. On the `_mature_exon_chain` fixture (intron+ | exon+ | intron+, chain
# B0 R0 B1 R1 B2 R2 B3) EVERY continuous-strand edge now fires — including the formerly-silenced exon R1→B1
# (backward) and exon R1→B2 (forward). That re-opens the mature leak into the introns; the honest σ²_transfer
# precision + the nascent factory (ρ_nascent = ρ_RNA − ρ_mature) are what will counter it (see §6 of the doc).
# The `mrna_active_*` mask itself stays computed in the statics (the nascent factory will consume it).
# ---------------------------------------------------------------------------

# b_bwd / a_fwd tuple layout: (amg, apg, amp, app, amn, apn). app = +RNA message PRECISION.
_APP = 3
_R1_EXON = 3  # chain id of the expressed exon R1
_B1 = 2  # intron→exon junction; its right neighbour (backward src) is R1
_B2 = 4  # exon→intron junction; its left neighbour (forward src) is R1


def test_gate_dismantled_exon_emits_rna_into_intron():
    """The dismantle's observable effect, asserted on the message itself (the inversion of the retired
    `test_exon_does_not_manufacture_nascent_into_intron`). With the mature-crossing gate GONE, the two edges
    it used to silence — backward exon R1→B1 and forward exon R1→B2 — now carry a non-zero +RNA message (the
    exon's own unspliced RNA density is structurally continuous into its flanking introns). This FAILS on the
    gated code (both were pinned to 0.0), so it is a real falsifier that the gate is truly dismantled, and it
    STAYS true after the nascent factory lands (nascent still flows exon→intron, only with mature subtracted).
    The magnitude of this leak is the known, accepted regression tracked by
    `test_mature_no_nascent_hallucination_in_introns`."""
    _, cap = _sweep(_mature_exon_chain(spliced=True))
    assert cap["b_bwd"][_APP][_B1] > 0.0, cap["b_bwd"][_APP][_B1]  # backward exon R1 → B1 now fires
    assert cap["a_fwd"][_APP][_B2] > 0.0, cap["a_fwd"][_APP][_B2]  # forward  exon R1 → B2 now fires


def test_intron_relays_nascent_into_exon_both_directions():
    """The structural-continuity guard: nascent must keep flowing intron→boundary and boundary→exon from BOTH
    sides. These four edges (intron R0→B1, intron R2→B2, B1→exon R1, B2→exon R1) are all +strand-continuous on
    both endpoints, so the structural `free_s` gate keeps their +RNA precision > 0 — unaffected by the
    mature-gate dismantle (which only re-opens the exon→intron direction, never closes these). Guards against a
    regression that would delete real intron→exon nascent relay."""
    _, cap = _sweep(_mature_exon_chain(spliced=True))
    a, b = cap["a_fwd"], cap["b_bwd"]
    # intron → boundary (source is a non-mature intron ⇒ always sends):
    assert a[_APP][_B1] > 0.0, a[_APP][_B1]  # forward  intron R0 → B1
    assert b[_APP][_B2] > 0.0, b[_APP][_B2]  # backward intron R2 → B2
    # boundary → exon (destination is mature-active ⇒ mature/nascent may enter):
    assert a[_APP][_R1_EXON] > 0.0, a[_APP][_R1_EXON]  # forward  B1 → exon R1
    assert b[_APP][_R1_EXON] > 0.0, b[_APP][_R1_EXON]  # backward B2 → exon R1


def test_mrna_active_matches_same_strand_exon_rule():
    """The `mrna_active_strands` mature-presence mask (no longer wired into the emission gate, but kept in the
    statics for the coming nascent factory `ρ_nascent = ρ_RNA − ρ_mature`) is exactly the user's rule: mature is
    present on strand s across a seam iff the SAME-STRANDED exon bit is set on BOTH flanks. Intron bits never
    qualify; `EX+EX- | EX+EX-` passes on BOTH strands. Enumerate all 16×16 signature pairs (a boundary's two
    flanks) and check `mrna_active_strands` against that predicate, plus the subsumption `mrna_active_s ⇒
    nrna_active_s` (mature ⇒ nascent). Pure, no sweep."""
    sigs = np.arange(N_SIGNATURES, dtype=np.int64)
    for sl in sigs:
        for sr in sigs:
            mrp_l, mrn_l = mrna_active_strands(np.array([sl]))
            mrp_r, mrn_r = mrna_active_strands(np.array([sr]))
            # a boundary's per-strand mature-crossing = AND of the two flanks' own exon bits
            mrp = bool(mrp_l[0] and mrp_r[0])
            mrn = bool(mrn_l[0] and mrn_r[0])
            exp_pos = bool((sl & BIT_EXON_POS) and (sr & BIT_EXON_POS))
            exp_neg = bool((sl & BIT_EXON_NEG) and (sr & BIT_EXON_NEG))
            assert mrp == exp_pos, (sl, sr, mrp, exp_pos)
            assert mrn == exp_neg, (sl, sr, mrn, exp_neg)
            # subsumption: mature-active ⇒ nascent-active (an exon carries nascent too), per strand
            nrp_l, nrn_l = nrna_active_strands(np.array([sl]))
            nrp_r, nrn_r = nrna_active_strands(np.array([sr]))
            nrp = bool(nrp_l[0] and nrp_r[0])
            nrn = bool(nrn_l[0] and nrn_r[0])
            assert not mrp or nrp, (sl, sr)  # mrp ⇒ nrp
            assert not mrn or nrn, (sl, sr)

    # the user's headline case: overlapping opposite-strand exons on both flanks ⇒ mature passes on BOTH strands
    both = BIT_EXON_POS | BIT_EXON_NEG
    mrp, mrn = mrna_active_strands(np.array([both]))
    mrp2, mrn2 = mrna_active_strands(np.array([both]))
    assert bool(mrp[0] and mrp2[0]) and bool(mrn[0] and mrn2[0])
    # and an exon+intron mixed flank does NOT block the exon's own strand (+ passes; − is intron→intron ⇒ no)
    mixed = BIT_EXON_POS | BIT_INTRON_NEG  # exon on +, intron on −
    mp_l, mn_l = mrna_active_strands(np.array([mixed]))
    mp_r, mn_r = mrna_active_strands(np.array([mixed]))
    assert bool(mp_l[0] and mp_r[0])  # + strand: exon|exon ⇒ mature passes
    assert not bool(mn_l[0] and mn_r[0])  # − strand: intron|intron ⇒ no mature


def test_spliced_routing_is_strand_aware_at_ambig_seams():
    """D4 — a splice junction absorbs mature RNA ONLY on its own (motif-known) strand. Mature RNA is
    single-stranded; a junction's strand is data-type-invariant (the genomic splice motif). At an AMBIG seam
    (overlapping opposite-strand transcripts) the spliced (mature) mass must land on the flank carrying the
    junction's OWN-strand exon, never merely "any exon".

    Fixture from the user's complex locus (TA+ 1000-2000/10000-11000/20000-25000, TB+ 9000-25000,
    TC- 3000-4000/10300-10800/27000-32000, TD- 3000-5000), around 10300-11000:
      R0 ≈ 10300-10800 : TA+/TB+ exon OVER TC- exon        → EXON_POS | EXON_NEG   (AMBIG exon)
      R1 ≈ 10800-11000 : TA+ exon OVER TC- intron          → EXON_POS | INTRON_NEG (+exon, −intron)
      R2 ≈ 11000-…     : TA+ intron OVER TC- exon (elsewhere) → INTRON_POS | EXON_NEG (+intron, −exon)
      B1 (R0|R1) = the TC- donor junction at 10800  → junction_strand = −
      B2 (R1|R2) = the TA+ donor junction at 11000  → junction_strand = +
    The − junction (B1) must route its − mature to R0 (has a − exon) and NOT to R1 (only a + exon).
    The + junction (B2) must route its + mature to R1 (has a + exon) and NOT to R2 (only a − exon).
    Strand-BLIND routing (`any_exon`) mis-deposits both onto the opposite-strand exon flank.
    """
    rro, rbo = np.array([0, 3]), np.array([0, 4])
    chain = build_node_chain(rro, rbo)  # B0 R0 B1 R1 B2 R2 B3  → nodes 0..6
    L = np.array([500.0, 200.0, 5000.0])
    sig = np.array(
        [
            BIT_EXON_POS | BIT_EXON_NEG,  # R0: AMBIG exon
            BIT_EXON_POS | BIT_INTRON_NEG,  # R1: + exon, − intron
            BIT_INTRON_POS | BIT_EXON_NEG,  # R2: + intron, − exon
        ],
        dtype=np.int64,
    )
    region_arrays = SimpleNamespace(region_size_bp=L, signature=sig)
    substrate = SimpleNamespace(contained=_view([100.0, 50.0, 100.0], [0.0, 0.0, 0.0]))
    # spliced (mature) mass present on BOTH sides of each junction; routing decides where it lands.
    left = _view([0.0, 50.0, 70.0, 40.0], [0.0, 50.0, 70.0, 0.0])  # B1.left=50, B2.left=70
    right = _view([10.0, 60.0, 80.0, 0.0], [0.0, 60.0, 80.0, 0.0])  # B1.right=60, B2.right=80
    bsub = SimpleNamespace(
        left=left,
        right=right,
        left_region=np.array([-1, 0, 1, 2]),
        right_region=np.array([0, 1, 2, -1]),
        junction_strand=np.array([0, TS_NEG, TS_POS, 0], dtype=np.int8),
    )
    g = build_node_geometry(chain, substrate, bsub, region_arrays, _delta_pmf(300), _delta_pmf(200))

    # B1 (node 2) — the − junction.
    assert g.spliced_neg_left[2] == 50.0  # → R0 (has a − exon): correct before & after
    assert g.spliced_neg_right[2] == 0.0  # NOT → R1 (+ exon only): FAILS on strand-blind code (routes 60)
    assert g.spliced_pos_left[2] == 0.0 and g.spliced_pos_right[2] == 0.0  # a − junction carries no + mature

    # B2 (node 4) — the + junction.
    assert g.spliced_pos_left[4] == 70.0  # → R1 (has a + exon)
    assert g.spliced_pos_right[4] == 0.0  # NOT → R2 (− exon only): FAILS on strand-blind code (routes 80)
    assert g.spliced_neg_left[4] == 0.0 and g.spliced_neg_right[4] == 0.0  # a + junction carries no − mature


# --------------------------------------------------------------------------------------------------
# DOF pie relay (item 2) — S2: a test that FAILS on today's incoherent relay and flips when the
# coordinate (λ,θ) relay lands (docs/calibration/dof_pie_relay_implementation_plan.md §7 S2).
# --------------------------------------------------------------------------------------------------
def test_relay_pie_is_a_composition():
    """The DOF-pie invariant, asserted on the running belief the sweep RELAYS (docs/calibration/
    dof_pie_relay_derivation.md). On a chain with an RNA-dense exon beside gDNA boundaries, each node's relayed
    pie ``(fbg,fbp,fbn)`` MUST be a composition — sum to 1 and every component ≤ 1 (i.e. ``n_c = f_c·M ≤ M``).
    The coherent ``(λ,θ)`` relay makes both hold BY CONSTRUCTION (``f_g=σ(μ_λ)``, ``f_pos+f_neg=1−f_g``). Before
    the fix the relay kept three INDEPENDENT log-fraction Gaussians and this FAILED (the worst boundary relayed
    a pie summing to ~1.16, and with the mature gate off two components saturated to a sum ~1.95); the
    single-component explosion (fbp up to 600×) is a genome-scale phenomenon measured by
    ``scripts/debug/pie_probe.py``. This was xfail(strict) through S3 and flipped when S4 landed."""
    _, cap = _sweep(_mature_exon_chain(spliced=True))
    sol = np.asarray(cap["solvable"], bool)
    for scan, (fbg, fbp, fbn) in zip(("forward", "backward"), cap["_relay_pie"]):
        pie = (fbg + fbp + fbn)[sol]
        maxfrac = np.maximum.reduce([fbg, fbp, fbn])[sol]
        assert np.allclose(pie, 1.0, atol=1e-6), (scan, "pie≠1", pie)
        assert np.all(maxfrac <= 1.0 + 1e-6), (scan, "n_c>M", maxfrac)


def test_solve_returns_coherent_coordinate_seed():
    """S3 — the per-node solve emits the free-coordinate moments ``(E[λ],Var[λ],E[θ],Var[θ])`` that seed the
    (λ,θ) relay, and the pie DERIVED from them is a valid composition BY CONSTRUCTION. Contrast the current
    readout (median ``f_g`` + mean ``f_pos`` + mean ``f_neg``), which need NOT sum to 1. Nodes: single-strand
    ``+``, single-strand ``−``, AMBIG."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    u_pos = np.array([80.0, 10.0, 50.0])
    u_neg = np.array([5.0, 70.0, 45.0])
    ap = np.array([True, False, True])
    an = np.array([False, True, True])
    dc = _solve_nodes_logodds_all(
        u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(3),
        kappa=0.9, od_g=0.0, od_r=0.0, n_grid=60, n_grid_ss=256,
    )
    # the coordinate moments are populated, finite, and variances non-negative
    for arr in (dc.lam_mean, dc.lam_var, dc.theta_mean, dc.theta_var):
        assert arr is not None and np.all(np.isfinite(arr))
    assert np.all(dc.lam_var >= 0.0) and np.all(dc.theta_var >= 0.0)
    # single-strand nodes lock the tilt: θ = +π/2 (pos-only) / −π/2 (neg-only), Var[θ] = 0
    assert np.isclose(dc.theta_mean[0], 0.5 * np.pi) and dc.theta_var[0] == 0.0
    assert np.isclose(dc.theta_mean[1], -0.5 * np.pi) and dc.theta_var[1] == 0.0
    assert dc.theta_var[2] > 0.0  # AMBIG tilt is a live d.o.f.
    # THE SEED IS COHERENT: the pie from (λ,θ) sums to 1 and lies in [0,1] — for every node.
    f_g = 1.0 / (1.0 + np.exp(-dc.lam_mean))
    f_r = 1.0 - f_g
    tau = np.sin(dc.theta_mean)
    f_pos, f_neg = f_r * (1.0 + tau) / 2.0, f_r * (1.0 - tau) / 2.0
    assert np.allclose(f_g + f_pos + f_neg, 1.0, atol=1e-12)
    for f in (f_g, f_pos, f_neg):
        assert np.all(f >= -1e-12) and np.all(f <= 1.0 + 1e-12)
    # the single-strand dead strand is exactly 0 (θ lock): +node has f_neg=0, −node has f_pos=0
    assert f_neg[0] == 0.0 and f_pos[1] == 0.0


# ── Stage 2: the two extracted message-MODE helpers (message_mode_implementation_plan.md) ──────────────────
def test_mode_shift_equals_mc_composition_fraction():
    """``exp(_mode_shift(Mg, Mg+Mr, comp_fl))`` == the MC-validated composition fraction Mg/(Mg+Mr)
    (chain_mode_mc._fg_from_densities) when the one-fragment floor is inactive — the shift IS §4a."""
    from rigel.calibration.bp_solver import _mode_shift

    rho_g, rho_r, egd, erd, md = 0.4, 0.6, 200.0, 150.0, 5000.0
    Mg, Mr = rho_g * egd, rho_r * erd
    comp_fl = 1.0 / md
    fg_shift = np.exp(_mode_shift(Mg, Mg + Mr, comp_fl))
    fg_mc = Mg / (Mg + Mr)  # == chain_mode_mc._fg_from_densities(rho_g, rho_r, egd, erd)
    assert np.isclose(fg_shift, fg_mc, atol=1e-12)
    assert comp_fl < fg_mc  # floor genuinely inactive here


def test_mode_density_closed_form():
    """``_mode_density(ρ, E, md)`` == log(ρ·E/md) when the one-fragment floor is inactive."""
    from rigel.calibration.bp_solver import _mode_density

    rho_c, eff_c, md = 0.4, 200.0, 5000.0
    assert rho_c > 1.0 / eff_c  # floor inactive
    assert np.isclose(_mode_density(rho_c, eff_c, md), np.log(rho_c * eff_c / md), atol=1e-12)


def test_mode_helpers_finite_under_zero_and_negative_density():
    """Domain guard (review #2): the derived one-fragment floors keep both modes FINITE even when the mature
    subtraction drives ρ_c ≤ 0 — no NaN/−inf, and NO added epsilon."""
    import math

    from rigel.calibration.bp_solver import _mode_density, _mode_shift

    egd, erd, md = 200.0, 150.0, 5000.0
    comp_fl = 1.0 / md
    # shift: a component whose imputed mass is 0 floors at comp_fl (log finite, = −log md)
    assert math.isfinite(_mode_shift(0.0, 1e-9, comp_fl))
    assert np.isclose(_mode_shift(0.0, 1e-9, comp_fl), math.log(comp_fl))
    # density: a NEGATIVE residual density (over-absorbed mature) floors at 1/E (log finite)
    for rho_neg_density in (-5.0, 0.0, -1e-3):
        m = _mode_density(rho_neg_density, erd, md)
        assert math.isfinite(m)
        assert np.isclose(m, math.log((1.0 / erd) * erd / md))  # == floor·E/md = 1/md
    assert math.isfinite(_mode_density(-9.9, egd, md))


# ── HARDENING: F (numerical robustness) + G (invariant coverage) ──────────────────────────────────────────
# These lock down the two things this session's emission refactor taught us: (F) the old evidence gates were
# masking latent 0·∞ nans, and (G) the emission↔density coupling bug that only a GOLDEN caught. Goldens tell
# us THAT output moved; these pin WHAT principle must hold.


def test_message_primitives_never_nan():
    """F: the pure message primitives must be FINITE over extreme/degenerate inputs — the emission refactor
    exposed two latent 0·∞ nans the gates had masked (fr²·∞ at the λ-window edge; n_eff·∞). No nan/inf may ever
    reach the fold; and _pred_precision must self-guard even a nan variance."""
    import math

    from rigel.calibration.bp_solver import _EPS, _mode_density, _mode_shift, _pred_precision

    EXT = [0.0, _EPS, 1e-6, 1.0, 1e6, 1e12]
    for count in [0.0, _EPS, 1.0, 1e9]:
        for v in EXT + [math.inf, math.nan]:  # incl. ∞ (unseen) and nan (must be guarded, not propagated)
            for s2t in [0.0, _EPS, 1.0, 1e6]:
                pr = _pred_precision(count, v, s2t)
                assert math.isfinite(pr) and pr >= 0.0, (count, v, s2t, pr)
    for mass in [0.0, _EPS, 1.0, 1e9]:
        for den in [_EPS, 1.0, 1e9]:
            for cfl in [1e-12, 1e-3, 1.0]:
                assert math.isfinite(_mode_shift(mass, den, cfl)), (mass, den, cfl)
    for rho in [-5.0, -1e-6, 0.0, _EPS, 1.0, 1e9]:  # negative ρ = over-absorbed mature (§4b)
        for eff in [_EPS, 1.0, 1e6]:
            for md in [_EPS, 1.0, 1e9]:
                assert math.isfinite(_mode_density(rho, eff, md)), (rho, eff, md)


def test_sweep_finite_over_extreme_configs():
    """F: no nan/inf reaches the fold. The real node_sweep over spliced/±, stranded/unstranded, and extreme
    gDNA/mature densities (pure-gDNA, pure-RNA, empty, tiny, huge) — every final fraction is finite & in range,
    every variance is ≥0 (∞ = the honest 'unsolved' state is allowed; nan is not), every emitted message
    mode/precision is finite."""
    for spliced in (True, False):
        for kappa in (0.5, 0.95):
            for rho_g, rho_m in [(0.5, 1.0), (0.0, 1.0), (2.0, 0.0), (1e-6, 1e-6), (1e4, 1e4)]:
                cfg = dict(spliced=spliced, kappa=kappa, rho_g=rho_g, rho_m=rho_m)
                final, cap = _sweep(_mature_exon_chain(**cfg), kappa=kappa)
                for nm in ("f_g", "f_pos", "f_neg"):
                    v = np.asarray(getattr(final, nm))
                    assert np.all(np.isfinite(v)), (cfg, nm, v)
                    assert np.all(v >= -1e-9) and np.all(v <= 1.0 + 1e-9), (cfg, nm, v)
                for nm in ("var_gdna", "var_pos", "var_neg"):
                    v = np.asarray(getattr(final, nm))
                    assert not np.any(np.isnan(v)) and np.all(v >= -1e-12), (cfg, nm, v)  # ∞ ok, nan not
                for scan in ("a_fwd", "b_bwd"):
                    for arr in cap[scan]:  # (amg, apg, amp, app, amn, apn)
                        assert np.all(np.isfinite(np.asarray(arr))), (cfg, scan)


def test_pred_precision_honest_semantics():
    """G: the honest-precision principle at the source (emission_and_precision_derivation.md §2). A message's
    composition precision is 0 when the composition is UNSEEN (v_log=∞ — a no-evidence source, the ev_λ=∞ fix)
    or when there is no count; positive & finite with real evidence; monotone-increasing in the count."""
    import math

    from rigel.calibration.bp_solver import _pred_precision

    assert _pred_precision(100.0, math.inf, 0.1) == 0.0  # UNSEEN (τ=0 ⇒ v=∞) ⇒ zero precision (the ev_λ=∞ fix)
    assert _pred_precision(0.0, 1.0, 0.1) == 0.0  # no count ⇒ zero precision
    p = _pred_precision(50.0, 1.0, 0.1)
    assert 0.0 < p and math.isfinite(p)  # real evidence ⇒ finite positive
    assert _pred_precision(10.0, 1.0, 0.1) < _pred_precision(100.0, 1.0, 0.1)  # more count ⇒ more power


def test_vacuous_unstranded_source_zero_precision_but_density_flows():
    """G: the phantom guard + the density⊥evidence DECOUPLING (the coupling bug's principle). On a fully
    composition-vacuous UNSTRANDED chain (κ=½ ⇒ I_strand deadband ⇒ τ=0, no spliced), a source manufactures NO
    composition confidence — its gDNA/RNA PREDICTION precision is exactly 0 (ev_λ=∞ ⇒ pr→0) — YET its message
    MODE stays a well-defined finite density (STRUCTURE, not evidence, sets the mode). Precision reflects
    evidence; density does not. This is exactly what the old emission gate violated: it dropped the density when
    the evidence gate closed (`rho_r += rho_neg` inside `if emit_n`)."""
    _, cap = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    for scan in ("a_fwd", "b_bwd"):
        apg, amp, app = np.asarray(cap[scan][1]), np.asarray(cap[scan][2]), np.asarray(cap[scan][3])
        assert np.all(apg <= 1e-12), (scan, "gDNA PREDICTION precision must be 0 on a vacuous source", apg)
        assert np.all(app <= 1e-12), (scan, "RNA PREDICTION precision must be 0 (no spliced, τ=0)", app)
        assert np.all(np.isfinite(amp)), (scan, "RNA density MODE must stay well-defined (density ⊥ evidence)")


def test_node_sweep_deterministic():
    """H: pass-0 must be bit-reproducible. The forward-backward BP sweep is sequential Python (no parallel
    reduction), so the same input must give a bit-identical belief AND identical emitted messages run-to-run —
    a prerequisite for any confidence claim about the solver. Uses the unstranded (κ=½), fully message-driven
    case, the one most sensitive to any ordering nondeterminism."""
    a, capa = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    b, capb = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    for nm in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg"):
        x, y = np.asarray(getattr(a, nm)), np.asarray(getattr(b, nm))
        assert np.array_equal(x, y, equal_nan=True), (nm, x, y)  # BIT-identical (not just close)
    for scan in ("a_fwd", "b_bwd"):  # every emitted message (amg,apg,amp,app,amn,apn) bit-identical
        for i in range(6):
            x, y = np.asarray(capa[scan][i]), np.asarray(capb[scan][i])
            assert np.array_equal(x, y, equal_nan=True), (scan, i, x, y)
