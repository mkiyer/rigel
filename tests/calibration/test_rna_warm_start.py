"""Unit tests for the Phase-A RNA warm-start projection (`_build_rna_warm_start` + `RnaWarmStart`).

The production golden fixtures carry no exon↔intron crossing or splice-junction mass, so those branches
would otherwise run only on zeros. These tests drive the projection on a hand-built 3-region chain with
NONZERO contained / crossing / spliced signal on BOTH strands, checking the three node roles independently:

  * CONTAINED — the region node's own belief × its contained mass / RNA eff-length, per strand.
  * CROSSING  — the seam boundary node's belief × the POOLED two-side mass / averaged support, left-keyed,
                and ZERO on the reference's last region (no seam to its right).
  * SPLICED   — the flanking boundary's one-sided motif-stranded mass / half-triangle support, kept in its
                own arrays (never summed into crossing), SINGLE-STRAND: a donor on region r's RIGHT side, an
                acceptor on region r+1's LEFT side, each tagged with its fixed motif strand.
"""

from __future__ import annotations

import numpy as np
import pytest
from types import SimpleNamespace

from rigel.calibration.calibrate import _build_rna_warm_start
from rigel.calibration.node_chain import build_node_chain
from rigel.calibration.node_geometry import NodeBelief, NodeGeometry
from rigel.calibration.result import RnaWarmStart
from rigel.calibration.signature import TS_NEG, TS_POS


def _chain_3_regions():
    """One reference, 3 regions → node order B0 R0 B1 R1 B2 R2 B3 (ids 0..6)."""
    return build_node_chain(
        np.array([0, 3], dtype=np.int64), np.array([0, 4], dtype=np.int64)
    )


def _geometry_7(**overrides) -> NodeGeometry:
    """A length-7 geometry with per-node arrays; every field positive so divisions are finite. Callers
    override the specific faces they exercise. gDNA eff-lengths are inert here (the builder reads only RNA)."""
    z = lambda: np.zeros(7)  # noqa: E731
    o = lambda: np.ones(7)  # noqa: E731
    fields = dict(
        n_nodes=7,
        mass_left=z(), mass_right=z(),
        eff_gdna_left=o(), eff_gdna_right=o(),
        eff_rna_left=o(), eff_rna_right=o(),
        eff_spl_left=o(), eff_spl_right=o(),
        spliced_pos_left=z(), spliced_pos_right=z(),
        spliced_neg_left=z(), spliced_neg_right=z(),
    )
    fields.update(overrides)
    return NodeGeometry(**fields)


def test_build_rna_warm_start_three_roles():
    chain = _chain_3_regions()
    # nodes:            B0    R0    B1    R1    B2    R2    B3
    f_pos = np.array([0.30, 0.70, 0.50, 0.40, 0.60, 0.05, 0.20])
    f_neg = np.array([0.20, 0.10, 0.10, 0.10, 0.10, 0.05, 0.30])
    f_g = 1.0 - f_pos - f_neg
    belief = NodeBelief(f_pos=f_pos, f_neg=f_neg, f_g=f_g,
                        var_pos=np.zeros(7), var_neg=np.zeros(7), var_gdna=np.zeros(7))

    mass_left = np.zeros(7)
    mass_right = np.zeros(7)
    # regions R0/R1/R2 (ids 1/3/5): contained mass, both faces equal.
    mass_left[[1, 3, 5]] = mass_right[[1, 3, 5]] = [100.0, 200.0, 300.0]
    # seam boundaries B1 (id 2, r0↔r1) and B2 (id 4, r1↔r2): per-side crossing mass.
    mass_left[[2, 4]] = mass_right[[2, 4]] = [10.0, 20.0]

    eff_rna_left = np.ones(7)
    eff_rna_right = np.ones(7)
    eff_rna_left[[1, 3, 5]] = eff_rna_right[[1, 3, 5]] = 50.0  # region contained RNA eff-len
    eff_rna_left[[2, 4]] = eff_rna_right[[2, 4]] = [5.0, 8.0]  # seam RNA eff-len

    # SPLICED (single-strand motif): B1 is a + junction (donor→R0, acceptor→R1); B2 is a − junction.
    spliced_pos_left = np.zeros(7)
    spliced_pos_right = np.zeros(7)
    spliced_neg_left = np.zeros(7)
    spliced_neg_right = np.zeros(7)
    spliced_pos_left[2], spliced_pos_right[2] = 8.0, 6.0   # B1 (+): left→R0 donor, right→R1 acceptor
    spliced_neg_left[4], spliced_neg_right[4] = 4.0, 10.0  # B2 (−): left→R1 donor, right→R2 acceptor
    eff_spl_left = np.ones(7)
    eff_spl_right = np.ones(7)
    eff_spl_left[2], eff_spl_right[2] = 4.0, 3.0
    eff_spl_left[4], eff_spl_right[4] = 2.0, 5.0

    geom = _geometry_7(
        mass_left=mass_left, mass_right=mass_right,
        eff_rna_left=eff_rna_left, eff_rna_right=eff_rna_right,
        eff_spl_left=eff_spl_left, eff_spl_right=eff_spl_right,
        spliced_pos_left=spliced_pos_left, spliced_pos_right=spliced_pos_right,
        spliced_neg_left=spliced_neg_left, spliced_neg_right=spliced_neg_right,
    )
    ws = _build_rna_warm_start(chain, belief, geom, SimpleNamespace(n_regions=3))

    # --- CONTAINED: f_s · M_u / eff_rna, per strand ---
    np.testing.assert_allclose(ws.rho_contained_pos, [0.70 * 2, 0.40 * 4, 0.05 * 6])
    np.testing.assert_allclose(ws.rho_contained_neg, [0.10 * 2, 0.10 * 4, 0.05 * 6])

    # --- CROSSING: seam belief × pooled mass / averaged support, left-keyed; last region is 0 (no seam) ---
    # r0 seam B1: m=20, s=0.5*(5+5)=5 → f_s·4;  r1 seam B2: m=40, s=8 → f_s·5;  r2: no seam.
    np.testing.assert_allclose(ws.rho_crossing_pos, [0.50 * 4, 0.60 * 5, 0.0])
    np.testing.assert_allclose(ws.rho_crossing_neg, [0.10 * 4, 0.10 * 5, 0.0])

    # --- SPLICED: single value + fixed motif strand; donor on r's RIGHT, acceptor on r+1's LEFT ---
    # R0: donor 8/4=2 (+) on right; no acceptor on left.
    # R1: acceptor 6/3=2 (+) on left; donor 4/2=2 (−) on right.
    # R2: acceptor 10/5=2 (−) on left; nothing on right.
    np.testing.assert_allclose(ws.rho_spliced_right, [2.0, 2.0, 0.0])
    np.testing.assert_allclose(ws.rho_spliced_left, [0.0, 2.0, 2.0])
    np.testing.assert_array_equal(ws.spliced_strand_right, np.array([TS_POS, TS_NEG, 0], dtype=np.int8))
    np.testing.assert_array_equal(ws.spliced_strand_left, np.array([0, TS_POS, TS_NEG], dtype=np.int8))

    # crossing and spliced are DISTINCT arrays — a splice boundary never inflates the crossing density.
    assert ws.rho_crossing_pos[0] == pytest.approx(2.0)  # R0 crossing (nascent) ...
    assert ws.rho_spliced_right[0] == pytest.approx(2.0)  # ... independent of its donor spliced (mature)


def test_build_rna_warm_start_zero_signal_is_valid_and_zero():
    """No RNA anywhere → all-zero densities and strand 0; the result still validates (the golden regime)."""
    chain = _chain_3_regions()
    belief = NodeBelief(f_pos=np.zeros(7), f_neg=np.zeros(7), f_g=np.ones(7),
                        var_pos=np.zeros(7), var_neg=np.zeros(7), var_gdna=np.zeros(7))
    ws = _build_rna_warm_start(chain, belief, _geometry_7(), SimpleNamespace(n_regions=3))
    for arr in (ws.rho_contained_pos, ws.rho_crossing_pos, ws.rho_spliced_left, ws.rho_spliced_right):
        np.testing.assert_array_equal(arr, np.zeros(3))
    np.testing.assert_array_equal(ws.spliced_strand_left, np.zeros(3, dtype=np.int8))


def _valid_kwargs(n=3):
    f = np.zeros(n)
    return dict(
        rho_contained_pos=f.copy(), rho_contained_neg=f.copy(),
        rho_crossing_pos=f.copy(), rho_crossing_neg=f.copy(),
        rho_spliced_left=f.copy(), rho_spliced_right=f.copy(),
        spliced_strand_left=np.zeros(n, dtype=np.int8),
        spliced_strand_right=np.zeros(n, dtype=np.int8),
    )


def test_rna_warm_start_rejects_bad_strand_value():
    kw = _valid_kwargs()
    kw["spliced_strand_left"] = np.array([0, 3, TS_POS], dtype=np.int8)  # 3 ∉ {0, POS, NEG}
    with pytest.raises(ValueError, match="0, TS_POS, or TS_NEG"):
        RnaWarmStart(**kw)


def test_rna_warm_start_rejects_non_int8_strand():
    kw = _valid_kwargs()
    kw["spliced_strand_right"] = np.zeros(3, dtype=np.int64)
    with pytest.raises(ValueError, match="int8 array"):
        RnaWarmStart(**kw)


def test_rna_warm_start_rejects_negative_density():
    kw = _valid_kwargs()
    kw["rho_contained_pos"] = np.array([-1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="non-negative"):
        RnaWarmStart(**kw)
