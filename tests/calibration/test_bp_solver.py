"""The BP-solver per-node geometry (`calibration.bp_solver.build_node_geometry`).

Verifies the per-face assembly onto the chain against the substrate: a region presents its contained geometry
both ways; a boundary presents its per-side crossing geometry (left-side / right-side, in the flank regions'
sizes); the one-sided motif spliced lands on the exon-flank face; per-component FL (gDNA vs RNA) is honoured.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.bp_solver import build_node_geometry
from rigel.calibration.effective_length import boundary_side_eff_length, region_eff_length
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS


def _view(mass_u, mass_spl):
    n = len(mass_u)
    z = np.zeros(n)
    return SimpleNamespace(
        n_unspliced_pos=z.copy(), n_unspliced_neg=z.copy(),
        n_spliced_sense=z.copy(), n_spliced_antisense=z.copy(),
        mass_unspliced=np.asarray(mass_u, float), mass_spliced=np.asarray(mass_spl, float),
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
        left=left, right=right,
        left_region=np.array([-1, 0, 1, 2]), right_region=np.array([0, 1, 2, -1]),
    )

    gdna_fl = _delta_pmf(300)
    rna_fl = _delta_pmf(200)
    g = build_node_geometry(chain, substrate, bsub, region_arrays, gdna_fl, rna_fl)

    # R0 (node 1): contained, same both faces; per-component FL.
    assert g.mass_left[1] == 100.0 and g.mass_right[1] == 100.0
    assert np.isclose(g.eff_gdna_left[1], region_eff_length(np.array([1000.0]), gdna_fl)[0])  # 1000-300=700
    assert np.isclose(g.eff_gdna_left[1], 700.0) and np.isclose(g.eff_rna_left[1], 800.0)
    assert g.spliced_pos_left[1] == 0.0 and g.spliced_pos_right[1] == 0.0

    # B1 (node 2): ex→in junction. left side in r0 (1000bp), right side in r1 (2000bp). spliced on LEFT (exon).
    assert g.mass_left[2] == 30.0 and g.mass_right[2] == 31.0
    assert np.isclose(g.eff_gdna_left[2], boundary_side_eff_length(gdna_fl, np.array([1000.0]))[0])  # min(300,1000)=300
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
    bsub = SimpleNamespace(left=left, right=right,
                           left_region=np.array([-1, 0, 1, 2]), right_region=np.array([0, 1, 2, -1]))
    g = build_node_geometry(chain, substrate, bsub, region_arrays, _delta_pmf(300), _delta_pmf(200))
    # B0 (node 0): left_region=-1 → left face eff = _EPS-floored ~0 (off edge); right face = r0 crossing.
    assert g.eff_gdna_left[0] <= 1e-6
    assert np.isclose(g.eff_gdna_right[0], 300.0)
    # B3 (node 6): right_region=-1 → right face off edge.
    assert g.eff_gdna_right[6] <= 1e-6
