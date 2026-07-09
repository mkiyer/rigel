"""The BP-solver per-node geometry (`calibration.bp_solver.build_node_geometry`).

Verifies the per-face assembly onto the chain against the substrate: a region presents its contained geometry
both ways; a boundary presents its per-side crossing geometry (left-side / right-side, in the flank regions'
sizes); the one-sided motif spliced lands on the exon-flank face; per-component FL (gDNA vs RNA) is honoured.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.bp_solver import (
    NodeBelief,
    NodeGeometry,
    adjacent_disagreement_variance,
    build_node_geometry,
    build_node_statics,
    node_sweep,
    init_beliefs,
    node_densities,
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
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
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
    )  # 1000-300=700
    assert np.isclose(g.eff_gdna_left[1], 700.0) and np.isclose(g.eff_rna_left[1], 800.0)
    assert g.spliced_pos_left[1] == 0.0 and g.spliced_pos_right[1] == 0.0

    # B1 (node 2): ex→in junction. left side in r0 (1000bp), right side in r1 (2000bp). spliced on LEFT (exon).
    assert g.mass_left[2] == 30.0 and g.mass_right[2] == 31.0
    assert np.isclose(
        g.eff_gdna_left[2], boundary_side_eff_length(gdna_fl, np.array([1000.0]))[0]
    )  # min(300,1000)=300
    assert np.isclose(g.eff_gdna_left[2], 300.0) and np.isclose(g.eff_gdna_right[2], 300.0)
    assert np.isclose(g.eff_rna_left[2], 200.0)  # RNA-FL min(200, 1000)
    assert g.spliced_pos_left[2] == 88.0 and g.spliced_pos_right[2] == 0.0  # exon on the left flank

    # B2 (node 4): in→ex junction. spliced on the RIGHT (exon r2) side.
    assert g.mass_left[4] == 20.0 and g.mass_right[4] == 22.0
    assert g.spliced_pos_left[4] == 0.0 and g.spliced_pos_right[4] == 77.0
    # right side is r2 (500bp): gDNA min(300,500)=300, RNA min(200,500)=200
    assert np.isclose(g.eff_gdna_right[4], 300.0) and np.isclose(g.eff_rna_right[4], 200.0)


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
    assert np.isclose(g.eff_gdna_right[0], 300.0)
    # B3 (node 6): right_region=-1 → right face off edge.
    assert g.eff_gdna_right[6] <= 1e-6


def test_node_densities_formula():
    # one node, two faces with different mass/eff-len; + spliced only on the left face.
    g = NodeGeometry(
        n_nodes=1,
        mass_left=np.array([100.0]),
        mass_right=np.array([200.0]),
        eff_gdna_left=np.array([200.0]),
        eff_gdna_right=np.array([400.0]),
        eff_rna_left=np.array([250.0]),
        eff_rna_right=np.array([500.0]),
        eff_spl_left=np.array([125.0]),  # one-sided spliced half-triangle eff-len (distinct from eff_rna)
        eff_spl_right=np.array([250.0]),
        spliced_pos_left=np.array([10.0]),
        spliced_pos_right=np.array([0.0]),
        spliced_neg_left=np.array([0.0]),
        spliced_neg_right=np.array([0.0]),
    )
    b = NodeBelief(
        f_pos=np.array([0.3]),
        f_neg=np.array([0.2]),
        f_g=np.array([0.5]),
        var_pos=np.array([0.0]),
        var_neg=np.array([0.0]),
        var_gdna=np.array([0.0]),
    )
    d = node_densities(b, g)
    assert np.isclose(d.rho_g_left[0], 0.5 * 100 / 200)  # 0.25
    # nascent rides eff_rna (250); the one-sided spliced rides its half-triangle eff_spl (125): 0.12+0.08=0.20
    assert np.isclose(d.rho_pos_left[0], 0.3 * 100 / 250 + 10 / 125)  # 0.20
    assert np.isclose(d.rho_neg_left[0], 0.2 * 100 / 250)  # 0.08 (no spliced on −)
    assert np.isclose(d.rho_g_right[0], 0.5 * 200 / 400)  # 0.25
    assert np.isclose(d.rho_pos_right[0], 0.3 * 200 / 500)  # no spliced on the right face: 0.12


def test_node_densities_factor1_under_uniform():
    # construct mass = ρ·eff-len per face (the accumulator's uniform-density deposit) ⇒ a pure-gDNA node's
    # ρ_g must read back ρ on BOTH faces despite different eff-lens (the factor-1 bedrock, formula level).
    rho = 0.42
    eff_g_l, eff_g_r = np.array([700.0]), np.array([300.0])
    g = NodeGeometry(
        n_nodes=1,
        mass_left=rho * eff_g_l,
        mass_right=rho * eff_g_r,
        eff_gdna_left=eff_g_l,
        eff_gdna_right=eff_g_r,
        eff_rna_left=np.array([800.0]),
        eff_rna_right=np.array([200.0]),
        eff_spl_left=np.array([400.0]),  # inert here (no spliced mass); present for the dataclass
        eff_spl_right=np.array([100.0]),
        spliced_pos_left=np.array([0.0]),
        spliced_pos_right=np.array([0.0]),
        spliced_neg_left=np.array([0.0]),
        spliced_neg_right=np.array([0.0]),
    )
    b = NodeBelief(
        f_pos=np.array([0.0]),
        f_neg=np.array([0.0]),
        f_g=np.array([1.0]),
        var_pos=np.array([0.0]),
        var_neg=np.array([0.0]),
        var_gdna=np.array([0.0]),
    )
    d = node_densities(b, g)
    assert np.isclose(d.rho_g_left[0], rho) and np.isclose(d.rho_g_right[0], rho)


def _empty_boundary_substrate(n_b):
    z = np.zeros(n_b)
    side = _cview(z.copy(), z.copy())
    return SimpleNamespace(
        left_region=np.full(n_b, -1), right_region=np.full(n_b, -1), left=side, right=side,
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
    allow_neg = np.array([False, False])
    strand_obs = allow_pos ^ allow_neg  # both single-strand
    mass = u_pos + u_neg
    d = _solve_nodes_logodds_all(u_pos, u_neg, allow_pos, allow_neg, strand_obs, mass, z,
                                 kappa=kappa, od_g=0.2, od_r=0.1, n_grid=60)
    assert d.gdna_frac_var is not None
    # same composition ⇒ same f_g; more fragments ⇒ sharper (lower Var(log f_g)).
    assert np.isclose(d.gdna_frac[0], d.gdna_frac[1])
    assert d.gdna_frac_var[1] < d.gdna_frac_var[0]
    # all per-component variances are present, finite, non-negative for active nodes.
    for v in (d.gdna_frac_var, d.rna_pos_frac_var, d.rna_neg_frac_var):
        assert np.all(np.isfinite(v)) and np.all(v >= 0.0)
    # a no-fragment node is inactive ⇒ zero variance on every component.
    d0 = _solve_nodes_logodds_all(np.array([0.0]), np.array([0.0]),
                                  np.array([True]), np.array([True]), np.array([False]),
                                  np.array([0.0]), np.array([0.0]), kappa=kappa, od_g=0.2, od_r=0.1,
                                  n_grid=60)
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
        left_region=lr, right_region=rr, left=left, right=right,
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
        max_passes=8,
        convergence_delta=1e-4,
        disagreement_sigma2=adjacent_disagreement_variance(chain, geom),
    )
    dens = node_densities(final, geom)
    interg, ambig = [1, 5], 3  # region nodes: R0/R2 intergenic (strand-anchored), R1 AMBIG
    # Intergenic nodes recover ρ exactly (the strand/signature pins them — this invariant holds in every phase).
    assert np.allclose(dens.rho_g_left[interg], rho, atol=0.02)
    assert np.allclose(dens.rho_g_right[interg], rho, atol=0.02)
    # AMBIG node — the documented PHASE-1 PRIOR-FREE gap (`_PHASE1_PRIOR_FREE`). A balanced AMBIG node has NO
    # intrinsic strand signal (κ-tilt undetermined) and Phase 1 deliberately makes the global gDNA prior
    # stability-only (precision capped at one pseudo-obs) with ê(z) OFF — so nothing pulls the balanced node
    # back to ρ, and the strand-mixture overdispersion biases it low (ρ_g≈0.28 vs 0.50). This is NOT a
    # regression: the AMBIG node is exactly what Phase 2 (the trained gDNA mixture prior + enrichment transfer)
    # resolves. When Phase 2 lands and flips `_PHASE1_PRIOR_FREE`, this recovers ρ=0.50 — restore atol=0.05 then.
    assert 0.20 < dens.rho_g_left[ambig] < 0.40
    assert 0.20 < dens.rho_g_right[ambig] < 0.40


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
        left_region=lr, right_region=rr, left=left, right=right,
        junction_strand=np.zeros(len(lr), dtype=np.int8),
    )

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.7, n_grid=40, statics=st
    )
    cap = {}
    final = node_sweep(
        chain, st, geom, belief, region_arrays, bsub, rna_sense_frac=0.7, n_grid=40, _capture=cap,
        disagreement_sigma2=adjacent_disagreement_variance(chain, geom),
    )
    exon = 3  # chain id of the single-exon gene (R1), flanked on both sides by TSS/TES seams
    # THE FIX — the exon receives a gDNA relay across the seam (incoming precision > 0). Pre-fix: 0 (no relay).
    assert cap["prec_g"][exon] > 0.0, "single-exon gene got NO gDNA relay across the TSS/TES seam"
    # The intergenic flanks emit ZERO RNA authority: the exon receives no +/− RNA message from them.
    assert cap["prec_p"][exon] == 0.0 and cap["prec_n"][exon] == 0.0
    # State ⊥ messages: the intergenic nodes stay locked all-gDNA (confident own-state, ignore all inputs).
    assert final.f_g[1] == 1.0 and final.f_g[5] == 1.0


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
        left_region=lr, right_region=rr, left=left, right=right,
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
        n_grid=40,
        max_passes=10,
        convergence_delta=1e-4,
        disagreement_sigma2=adjacent_disagreement_variance(chain, geom),
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
    # node's own strand evidence) → they stay firmly RNA-leaning (~0.22 under the belief-free σ²_imp message
    # precision on this ARTIFICIAL no-intergenic-seed worst case; real libraries with intergenic seeds pin
    # them harder). The retired legacy σ²_edge basis reached ~0.13.
    assert final.f_g[1] < 0.25 and final.f_g[5] < 0.25


# --- density-Gaussian message form: two-sided pull + emergent deference (precision_state_design.md §5) --------


def test_density_message_two_sided_mode_not_vertex():
    """A LOG-fraction gDNA message (mode=log 0.2, strong prec) on a balanced AMBIG node (flat strand) pulls
    f_g TOWARD 0.2 — two-sided by construction (a Gaussian on log f_g, no edge wall), not to the f_g=1
    vertex. The log-density log-odds solver's message form."""
    from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all

    z = np.zeros(1)
    # AMBIG node, balanced counts ⇒ the strand is flat (κ=0.5); only the message shapes f_g.
    d = _solve_nodes_logodds_all(
        np.array([50.0]), np.array([50.0]),
        np.array([True]), np.array([True]), np.array([False]),
        np.array([100.0]), z, kappa=0.5, od_g=0.0, od_r=0.0, n_grid=80,
        gdna_imp_mode=np.array([np.log(0.2)]), gdna_imp_prec=np.array([200.0]),
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
        np.array([1000.0]), np.array([5.0]),
        np.array([True]), np.array([False]), np.array([True]),
        np.array([1005.0]), z, kappa=0.99, od_g=0.0, od_r=0.0, n_grid=80,
        gdna_imp_mode=np.array([np.log(0.9)]), gdna_imp_prec=np.array([3.0]),
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
    Eg = region_eff_length(L, gdna_fl)            # contained gDNA eff-len  (1700)
    Er = region_eff_length(L, rna_fl)             # contained RNA  eff-len  (1800)
    g_half = rho_g * Eg / 2.0                       # per-strand contained gDNA count (balanced)
    mat = rho_m * Er                                # contained mature count (+ strand only)
    # exon R1: gDNA (balanced) + mature (+); introns R0/R2: gDNA only.
    u_pos = np.array([g_half[0], g_half[1] + mat[1], g_half[2]])
    u_neg = np.array([g_half[0], g_half[1], g_half[2]])
    contained = _cview(u_pos, u_neg, mass_u=u_pos + u_neg, mass_spl=np.zeros(3))
    substrate = SimpleNamespace(contained=contained)
    # boundary crossings: gDNA only (balanced); spliced(mature) on the EXON flank of the two junctions.
    side_g = boundary_side_eff_length(gdna_fl, L)   # (300)
    spl_eff = spliced_side_eff_length(rna_fl, L)     # one-sided half-triangle (100)
    lr, rr = np.array([-1, 0, 1, 2]), np.array([0, 1, 2, -1])
    def gx(r):
        return np.where(r >= 0, rho_g * side_g[np.clip(r, 0, 2)], 0.0)

    lcross, rcross = gx(lr), gx(rr)
    # ``spl_scale`` < 1 models a CAPTURE-DEPLETED junction: junction-spanning (spliced) reads are only
    # partially captured, so the junction under-reports the exon's true mature density ⇒ the B→exon mature
    # MEASUREMENT DISAGREES with the exon's own (confident) unspliced belief. Used by the silencing test.
    spl_mat = rho_m * spl_eff[1] * spl_scale         # mature spliced count on the exon flank
    # B1 (idx1): exon R1 is its RIGHT region → spliced on the RIGHT side. B2 (idx2): exon is its LEFT
    # region → spliced on the LEFT side. (One-sided, exon flank.)
    spl_l = np.array([0.0, 0.0, spl_mat if spliced else 0.0, 0.0])
    spl_r = np.array([0.0, spl_mat if spliced else 0.0, 0.0, 0.0])
    left = _cview(lcross / 2, lcross / 2, spl_sense=spl_l, mass_u=lcross,
                  mass_spl=spl_l.copy())
    right = _cview(rcross / 2, rcross / 2, spl_sense=spl_r, mass_u=rcross,
                   mass_spl=spl_r.copy())
    bsub = SimpleNamespace(
        left_region=lr, right_region=rr, left=left, right=right,
        junction_strand=np.array([0, TS_POS, TS_POS, 0], dtype=np.int8),
    )
    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(chain, substrate, bsub, region_arrays,
                          rna_sense_frac=kappa, n_grid=60, statics=st)
    return chain, st, geom, belief, region_arrays, bsub


def _sweep(args, kappa=0.95):
    chain, st, geom, belief, ra, bsub = args
    cap = {}
    final = node_sweep(chain, st, geom, belief, ra, bsub, rna_sense_frac=kappa, n_grid=60, _capture=cap,
                       disagreement_sigma2=adjacent_disagreement_variance(chain, geom))
    return final, cap


def test_mature_no_nascent_hallucination_in_introns():
    """Regression guard for the user's red line: a pure-mature expressed exon (nascent=0) must NOT
    manufacture wholesale nascent in its flanking pure-gDNA introns. The introns stay gDNA (truth f_g=1).
    (Here the strand is decisive, so the final f_g is also protected by the strand likelihood — the
    message-level proof that the *absorption* is what removes the leak is the next test.)"""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_introns = fin_m.f_g[[1, 5]]  # chain ids of R0, R2
    assert np.all(fg_introns > 0.85), fg_introns


def test_mature_absorption_lowers_nascent_message_into_junction():
    """The mechanism, asserted directly on the message (not the strand-masked final f_g). The exon→junction
    message is the BACKWARD scan's +RNA message into B1 (B1's right neighbour is the exon R1). With the
    junction spliced present, the message density subtracts the measured mature (``− SP[df]/E_spl_dst``),
    so a pure-mature exon imputes ≈0 nascent; zeroing the spliced removes the subtraction and the exon's
    RNA is imputed wholesale as nascent. The mode is a log-fraction target → lower (more negative) ⇒ less
    nascent. This is "the boundary spliced mass absorbs the incoming RNA; the remainder is nascent.\""""
    _, cap_m = _sweep(_mature_exon_chain(spliced=True))
    _, cap_n = _sweep(_mature_exon_chain(spliced=False))
    b1 = 2  # chain id of B1 (the intron→exon junction; its right neighbour R1 is the exon)
    amp_idx = 2  # b_bwd tuple = (amg, apg, amp, app, amn, apn); amp = +RNA message mode
    mode_mature = cap_m["b_bwd"][amp_idx][b1]
    mode_nascent = cap_n["b_bwd"][amp_idx][b1]
    # both must be real messages (emit_p fired across the +continuous junction)…
    assert cap_m["b_bwd"][amp_idx + 1][b1] > 0.0 and cap_n["b_bwd"][amp_idx + 1][b1] > 0.0
    # …and absorption makes the imputed nascent strictly lower (here ≈0 mature-only exon vs full RNA).
    assert mode_mature < mode_nascent - 0.5, (mode_mature, mode_nascent)


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
    assert cap_lo["mode_p"][ex] < cap_ok["mode_p"][ex] - 0.3, (cap_lo["mode_p"][ex], cap_ok["mode_p"][ex])
    # (2) …yet the exon's gDNA fraction barely moves (silenced — pre-fix the low measurement inflated f_g).
    assert abs(float(fin_lo.f_g[ex]) - float(fin_ok.f_g[ex])) < 0.05, (fin_lo.f_g[ex], fin_ok.f_g[ex])
    # (3) and the exon stays RNA-dominated, not pulled toward phantom gDNA.
    assert float(fin_lo.f_g[ex]) < 0.45, fin_lo.f_g[ex]


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
        return float(_solve_nodes_logodds_all(
            np.array([float(u_pos)]), np.array([float(u_neg)]),
            np.array([True]), np.array([False]), np.array([True]),
            np.array([n]), z, kappa=0.7, od_g=od, od_r=od, n_grid=80).gdna_frac[0])

    # pure gDNA (exact 50/50): near-binomial od → confidently gDNA; inflated od → materially under-called.
    fg_near = solve(500, 500, 0.034)
    fg_infl = solve(500, 500, 0.143)
    assert fg_near > 0.8, fg_near
    assert fg_infl < fg_near - 0.15, (fg_infl, fg_near)  # the inflated prior demonstrably degrades the call
    # pure RNA (+frac = κ = 0.7): near-binomial stays RNA-dominated; the inflated prior's symmetric harm is
    # MORE false gDNA on RNA too (it pulls every node toward ½). (At this intermediate κ the gDNA/RNA means
    # are close, so a small residual f_g is inherent — the point is near-binomial is cleaner.)
    rna_near = solve(700, 300, 0.034)
    rna_infl = solve(700, 300, 0.143)
    assert rna_near < 0.25, rna_near
    assert rna_infl > rna_near, (rna_infl, rna_near)
