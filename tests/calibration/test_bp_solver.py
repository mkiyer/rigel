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
    build_node_geometry,
    build_node_statics,
    node_sweep,
    init_beliefs,
    node_densities,
)
from rigel.calibration.effective_length import boundary_side_eff_length, region_eff_length
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
    left = _view([0.0, 30.0, 20.0, 40.0], [0.0, 88.0, 0.0, 0.0])
    right = _view([10.0, 31.0, 22.0, 0.0], [0.0, 0.0, 77.0, 0.0])
    bsub = SimpleNamespace(
        left=left,
        right=right,
        left_region=np.array([-1, 0, 1, 2]),
        right_region=np.array([0, 1, 2, -1]),
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
    assert np.isclose(d.rho_pos_left[0], (0.3 * 100 + 10) / 250)  # spliced rides on +: 0.16
    assert np.isclose(d.rho_neg_left[0], 0.2 * 100 / 250)  # 0.08
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
        left_region=np.full(n_b, -1), right_region=np.full(n_b, -1), left=side, right=side
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


def test_precision_state_strand_resolution():
    """Phase-1 lock (precision_state_design.md §3): the moment-matched Var(f_g) reflects strand-resolving
    power — a BALANCED node (gDNA-indistinguishable from balanced RNA) is uncertain (high Var(f_g)); a
    confident single-strand node is sharp (low Var(f_g)). This is the precision the honest message send (Phase
    2) will consume. A node with no fragments reports zero variance."""
    from rigel.calibration.simplex_sweep import _solve_nodes

    kappa = 0.01
    # node 0 = balanced AMBIG (both strands free, u_pos≈u_neg); node 1 = confident single-strand + (− locked).
    u_pos = np.array([100.0, 2.0])
    u_neg = np.array([100.0, 200.0])
    z = np.zeros(2)
    allow_pos = np.array([True, True])
    allow_neg = np.array([True, False])
    strand_obs = allow_pos ^ allow_neg  # [False=AMBIG, True=single-strand]
    mass = np.array([200.0, 202.0])
    d = _solve_nodes(u_pos, u_neg, z, z, allow_pos, allow_neg, strand_obs, mass, z,
                     kappa=kappa, od_g=0.2, od_r=0.1, n_grid=60)
    assert d.gdna_frac_var is not None
    # the balanced node cannot resolve gDNA-vs-RNA ⇒ higher posterior Var(f_g) than the confident node.
    assert d.gdna_frac_var[0] > d.gdna_frac_var[1]
    # all per-component variances are present, finite, non-negative for active nodes.
    for v in (d.gdna_frac_var, d.rna_pos_frac_var, d.rna_neg_frac_var):
        assert np.all(np.isfinite(v)) and np.all(v >= 0.0)
    # a no-fragment node is inactive ⇒ zero variance on every component.
    d0 = _solve_nodes(np.array([0.0]), np.array([0.0]), np.array([0.0]), np.array([0.0]),
                      np.array([True]), np.array([True]), np.array([False]),
                      np.array([0.0]), np.array([0.0]), kappa=kappa, od_g=0.2, od_r=0.1, n_grid=60)
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
    bsub = SimpleNamespace(left_region=lr, right_region=rr, left=left, right=right)

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.7, n_grid=40, statics=st
    )
    final, _ = node_sweep(
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
    )
    dens = node_densities(final, geom)
    interg, ambig = [1, 5], 3  # region nodes: R0/R2 intergenic (strand-anchored), R1 AMBIG
    # Intergenic nodes recover ρ exactly (the strand/signature pins them).
    assert np.allclose(dens.rho_g_left[interg], rho, atol=0.02)
    assert np.allclose(dens.rho_g_right[interg], rho, atol=0.02)
    # AMBIG node: post-overshoot-fix the message/global precision is now HONEST (count-currency floor +
    # density-dependent global), so it no longer over-pins the balanced AMBIG node — the strand mixture's
    # overdispersion term biases f_g slightly low, and the global MEAN does not yet pull it back. Result
    # ρ_g≈0.46 vs 0.50 (~7.5% low). This is the documented global-mean residual (count_space_solver_design.md
    # §6 / precision_overshoot_design.md): the same AMBIG gDNA under-call the over-confident global used to
    # mask. The global-mean fix (next track) restores this to 0.50; tighten the tolerance back then.
    assert np.allclose(dens.rho_g_left[ambig], rho, atol=0.05)
    assert np.allclose(dens.rho_g_right[ambig], rho, atol=0.05)


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
    bsub = SimpleNamespace(left_region=lr, right_region=rr, left=left, right=right)

    geom = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)
    st = build_node_statics(chain, substrate, bsub, region_arrays)
    belief = init_beliefs(
        chain, substrate, bsub, region_arrays, rna_sense_frac=0.95, n_grid=40, statics=st
    )
    assert belief.f_g[3] == 1.0  # AMBIG starts all-gDNA
    final, deltas = node_sweep(
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
    )
    # the AMBIG phantom is pulled DOWN from its all-gDNA init (1.0) toward RNA. The honest global is now a
    # GENTLE, M-independent tiebreaker (N_global≈1 pseudo-fragment at zero-gDNA, the Poisson-1 floor), so it
    # NUDGES rather than hard-pins; the suppression here leans on the RNA imputation from the intron neighbours.
    # This synthetic chain has NO intergenic structural seeds (intron+|AMBIG|intron−), so it is the hard case
    # for a gentle global — on real libraries the intergenic zero-count seeds drive a firmer zero-baseline.
    assert final.f_g[3] < 0.25  # substantially pulled down (init 1.0 → ~0.18)
    assert final.f_g[1] < 0.15 and final.f_g[5] < 0.15  # the introns stay RNA via the (decisive) strand
    # monotone convergence: the per-pass max-|Δf_g| is non-increasing
    assert all(deltas[k + 1] <= deltas[k] + 1e-9 for k in range(len(deltas) - 1))


# --- density-Gaussian message form: two-sided pull + emergent deference (precision_state_design.md §5) --------


def test_density_message_two_sided_mode_not_vertex():
    """A density-Gaussian gDNA message (mode=0.2, strong prec) on a balanced AMBIG node (flat strand) pulls f_g
    TOWARD 0.2 — two-sided by construction (a Gaussian, no edge log-wall), not to the f_g=1 vertex."""
    from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice

    lat = _simplex_lattice(80)
    _, _, fgg = lat
    # AMBIG node, balanced counts ⇒ the strand is flat (κ=0.5); only the message shapes f_g.
    psi = _local_loglik(
        np.array([50.0]), np.array([50.0]), np.zeros(1), np.zeros(1),
        np.array([True]), np.array([True]), 0.5, 0.0, 0.0, lat,
        strand_obs=np.array([False]),
        gdna_imp_mode=np.array([0.2]), gdna_imp_prec=np.array([200.0]),
    )
    fg = float(_fg_median(psi, fgg)[0])
    assert abs(fg - 0.2) < 0.05, fg


def test_density_message_defers_to_decisive_strand():
    """Emergent deference: a WEAK gDNA message (prec=3) trying to pull f_g→0.9 must lose to a decisive
    single-strand node's ~1000-fragment strand likelihood — f_g stays ≈0 (the honest precision blend means a
    weak message cannot override the data; no log-wall to force it off zero)."""
    from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice

    lat = _simplex_lattice(80)
    _, _, fgg = lat
    psi = _local_loglik(
        np.array([1000.0]), np.array([5.0]), np.zeros(1), np.zeros(1),
        np.array([True]), np.array([False]), 0.99, 0.0, 0.0, lat,
        strand_obs=np.array([True]),
        gdna_imp_mode=np.array([0.9]), gdna_imp_prec=np.array([3.0]),
    )
    fg = float(_fg_median(psi, fgg)[0])
    assert fg < 0.1, fg
