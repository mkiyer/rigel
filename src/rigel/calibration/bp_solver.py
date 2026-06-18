"""The belief-propagation calibration solver (density imputation, bidirectional sweep).

`docs/calibration/bp_sweep_rebuild_plan.md`. Rebuilds the per-node deconvolution as a pairwise, bidirectional
belief-propagation sweep over the unified region↔boundary chain (`node_chain`). This module is built in phases:

* **P1 (here):** the per-node static GEOMETRY on the chain (`build_node_geometry`) — each node's unspliced mass,
  per-component (gDNA-FL / RNA-FL) effective lengths, and the boundary's one-sided motif-stranded spliced mass,
  assembled **per face** (left/right) because a boundary's two sides lie in different-sized regions and the
  message uses the side FACING the destination. Regions present the same (contained) geometry to both faces.
  (init fractions + the directional sweep land next.)

Per-component fragment length (the v2-critique fix): gDNA crossings/contained use the gDNA FL; RNA — nascent
unspliced + the spliced — uses the RNA FL. The eff-len primitives are the single-owner `effective_length.py`
functions (`region_eff_length` = `E[max(0,L−ℓ)]` contained support; `boundary_side_eff_length` =
`E[min(ℓ,L_σ)]` per-side crossing support). The spliced is one-sided (its exon flank) and motif-stranded.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .effective_length import boundary_side_eff_length, region_eff_length
from .node_chain import BOUNDARY, REGION, NodeChain
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS

__all__ = ["NodeGeometry", "build_node_geometry"]

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class NodeGeometry:
    """Per-node, per-face static geometry on the chain (length ``n_nodes``). The ``_left``/``_right`` arrays
    are what a node presents to its left/right neighbour: for a region both equal its contained geometry; for a
    boundary they are its left-side / right-side crossing geometry. ``spliced_*`` is the one-sided
    motif-stranded spliced mass (0 on regions and on the non-exon face)."""

    n_nodes: int
    mass_left: np.ndarray  # float64 — unspliced mass facing the left neighbour
    mass_right: np.ndarray
    eff_gdna_left: np.ndarray  # gDNA-FL eff-len facing left (region contained / boundary left-side crossing)
    eff_gdna_right: np.ndarray
    eff_rna_left: np.ndarray  # RNA-FL eff-len facing left
    eff_rna_right: np.ndarray
    spliced_pos_left: np.ndarray  # + motif spliced mass on the left face (boundary, exon-on-left)
    spliced_pos_right: np.ndarray
    spliced_neg_left: np.ndarray
    spliced_neg_right: np.ndarray


def build_node_geometry(
    chain: NodeChain,
    substrate,
    boundary_substrate,
    region_arrays,
    gdna_fl_pmf: np.ndarray,
    rna_fl_pmf: np.ndarray,
) -> NodeGeometry:
    """Assemble the per-node per-face geometry from the region-/boundary-keyed substrates onto the chain."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY

    # --- per-region quantities (contained) ---
    L = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    reg_eff_g = region_eff_length(L, gdna_fl_pmf)  # E_gdna[max(0,L−ℓ)] contained
    reg_eff_r = region_eff_length(L, rna_fl_pmf)
    side_eff_g = boundary_side_eff_length(gdna_fl_pmf, L)  # per region: E_gdna[min(ℓ,L)] side crossing
    side_eff_r = boundary_side_eff_length(rna_fl_pmf, L)
    reg_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    R = L.shape[0]

    # --- per-boundary quantities (two sides) ---
    bsub = boundary_substrate
    blr = np.asarray(bsub.left_region, dtype=np.int64)
    brr = np.asarray(bsub.right_region, dtype=np.int64)
    bmass_l = np.asarray(bsub.left.mass_unspliced, dtype=np.float64)
    bmass_r = np.asarray(bsub.right.mass_unspliced, dtype=np.float64)
    bspl_l = np.asarray(bsub.left.mass_spliced, dtype=np.float64)  # sense+antisense summed
    bspl_r = np.asarray(bsub.right.mass_spliced, dtype=np.float64)
    B = blr.shape[0]
    # a boundary's per-side crossing eff-len = its flank region's E[min(ℓ,L)] (0 at a −1 terminal flank)
    b_eff_g_l = np.where(blr >= 0, side_eff_g[np.clip(blr, 0, R - 1)], 0.0)
    b_eff_g_r = np.where(brr >= 0, side_eff_g[np.clip(brr, 0, R - 1)], 0.0)
    b_eff_r_l = np.where(blr >= 0, side_eff_r[np.clip(blr, 0, R - 1)], 0.0)
    b_eff_r_r = np.where(brr >= 0, side_eff_r[np.clip(brr, 0, R - 1)], 0.0)
    # spliced face: the exon flank for that strand's junction (clean exon↔intron transition); else 0.
    sig_l = np.where(blr >= 0, sig[np.clip(blr, 0, R - 1)], 0)
    sig_r = np.where(brr >= 0, sig[np.clip(brr, 0, R - 1)], 0)

    def _spliced_faces(exon_bit, intron_bit):
        le = (sig_l & exon_bit) != 0
        li = (sig_l & intron_bit) != 0
        re = (sig_r & exon_bit) != 0
        ri = (sig_r & intron_bit) != 0
        side_l = le & ~li & ri & ~re  # exon on the left flank → spliced on the left side
        side_r = li & ~le & re & ~ri  # exon on the right flank → spliced on the right side
        return np.where(side_l, bspl_l, 0.0), np.where(side_r, bspl_r, 0.0)

    b_spl_pos_l, b_spl_pos_r = _spliced_faces(BIT_EXON_POS, BIT_INTRON_POS)
    b_spl_neg_l, b_spl_neg_r = _spliced_faces(BIT_EXON_NEG, BIT_INTRON_NEG)

    # --- gather onto the chain (region nodes ← region arrays; boundary nodes ← boundary arrays) ---
    ri_ = np.clip(idx, 0, R - 1)  # region index where is_reg
    bi_ = np.clip(idx, 0, B - 1)  # boundary index where is_bnd

    def pick(reg_vals, bnd_vals):
        return np.where(is_reg, reg_vals[ri_], np.where(is_bnd, bnd_vals[bi_], 0.0))

    # region presents the same (contained) geometry both ways; boundary presents its per-side geometry.
    mass_left = pick(reg_mass, bmass_l)
    mass_right = pick(reg_mass, bmass_r)
    eff_gdna_left = pick(reg_eff_g, b_eff_g_l)
    eff_gdna_right = pick(reg_eff_g, b_eff_g_r)
    eff_rna_left = pick(reg_eff_r, b_eff_r_l)
    eff_rna_right = pick(reg_eff_r, b_eff_r_r)
    zeros_R = np.zeros(R)
    spliced_pos_left = pick(zeros_R, b_spl_pos_l)
    spliced_pos_right = pick(zeros_R, b_spl_pos_r)
    spliced_neg_left = pick(zeros_R, b_spl_neg_l)
    spliced_neg_right = pick(zeros_R, b_spl_neg_r)

    return NodeGeometry(
        n_nodes=int(chain.n_nodes),
        mass_left=mass_left, mass_right=mass_right,
        eff_gdna_left=np.maximum(eff_gdna_left, _EPS), eff_gdna_right=np.maximum(eff_gdna_right, _EPS),
        eff_rna_left=np.maximum(eff_rna_left, _EPS), eff_rna_right=np.maximum(eff_rna_right, _EPS),
        spliced_pos_left=spliced_pos_left, spliced_pos_right=spliced_pos_right,
        spliced_neg_left=spliced_neg_left, spliced_neg_right=spliced_neg_right,
    )
