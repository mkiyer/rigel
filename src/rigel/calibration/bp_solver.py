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
from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_POS,
)
from .simplex_sweep import _solve_nodes, deconv_regions_sweep

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "NodeBelief",
    "NodeDensities",
    "node_densities",
    "init_beliefs",
]

_EPS = 1.0e-9
_POS_BITS = BIT_EXON_POS | BIT_INTRON_POS
_NEG_BITS = BIT_EXON_NEG | BIT_INTRON_NEG


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


@dataclass(frozen=True, slots=True)
class NodeBelief:
    """Per-node solved state on the chain: the pie `(f_pos, f_neg, f_g)` over the node's UNSPLICED mass + the
    per-component variances `(var_pos, var_neg, var_gdna)` (precision = 1/var). All length ``n_nodes``."""

    f_pos: np.ndarray
    f_neg: np.ndarray
    f_g: np.ndarray
    var_pos: np.ndarray
    var_neg: np.ndarray
    var_gdna: np.ndarray


@dataclass(frozen=True, slots=True)
class NodeDensities:
    """Per-node, per-face component densities — what a node presents to its left/right neighbour. gDNA uses the
    gDNA-FL eff-len; RNA (nascent unspliced + the boundary's one-sided spliced) uses the RNA-FL eff-len. These
    are the identity-mean MESSAGE values (`ARCHITECTURE §1.2`): a source's facing density IS the destination's
    prior mean. All length ``n_nodes``."""

    rho_g_left: np.ndarray
    rho_pos_left: np.ndarray
    rho_neg_left: np.ndarray
    rho_g_right: np.ndarray
    rho_pos_right: np.ndarray
    rho_neg_right: np.ndarray


def node_densities(belief: NodeBelief, geometry: NodeGeometry) -> NodeDensities:
    """Per-component, per-face densities from the belief pie + the static geometry. gDNA:
    ``ρ_g(σ) = f_g·M_σ / E_gdna(σ)``. RNA: ``ρ_s(σ) = (f_s·M_σ + spliced_s(σ)) / E_rna(σ)`` — the spliced
    (boundary-owned, one-sided) rides on its motif strand. The shipped factor-1-under-uniform construction makes
    a region's contained density and a boundary side's crossing density both equal the true ρ under uniform."""
    g = geometry
    fp, fn, fg = belief.f_pos, belief.f_neg, belief.f_g
    return NodeDensities(
        rho_g_left=fg * g.mass_left / g.eff_gdna_left,
        rho_pos_left=(fp * g.mass_left + g.spliced_pos_left) / g.eff_rna_left,
        rho_neg_left=(fn * g.mass_left + g.spliced_neg_left) / g.eff_rna_left,
        rho_g_right=fg * g.mass_right / g.eff_gdna_right,
        rho_pos_right=(fp * g.mass_right + g.spliced_pos_right) / g.eff_rna_right,
        rho_neg_right=(fn * g.mass_right + g.spliced_neg_right) / g.eff_rna_right,
    )


# ---------------------------------------------------------------------------
# Initialization (plan v3 §6) — the signature-binary G1/G2/G3 belief on the chain.
# ---------------------------------------------------------------------------
#
# Variance convention: ``0`` = a LOCKED axis (the fraction is pinned, e.g. a forbidden strand or an
# intergenic sink — certain) and ``inf`` = NO self-information (an allowed-but-unresolved axis; the sweep's
# precision ``1/var`` → 0 there, so the node is fully governed by neighbour messages + the global prior). The
# hard locking itself is enforced by the per-node ``allow_pos``/``allow_neg`` forbid mask in the solve, not by
# ``var=0``; the variance only records confidence. ``inf`` is the honest "no prior" value (magic-free) and is
# only ever consumed as ``1/var`` (= 0), never in a sum or ratio.


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """Build the six per-node belief arrays for ONE node type (regions OR boundaries) from its
    signature-binary classification + its strand-only solve.

    ``free_pos``/``free_neg`` are the per-node booleans for whether each strand's RNA axis is admissible (a
    region's own ±transcript bits; a boundary's ±strand CONTINUITY across the seam). ``deconv`` is the
    strand-only :class:`NodeDeconv` (no global, no imputation). The signature-binary default is all-gDNA
    ``{0,0,1}`` with ``var_gdna=inf``, ``var_pos=inf`` iff ``free_pos`` else ``0`` (locked), symmetric for neg
    (`ARCHITECTURE §3`). The class overrides:

    * **G1** (neither strand free — intergenic region / no-RNA-crossing boundary): a LOCKED gDNA sink — keep
      ``{0,0,1}`` and pin ``var_gdna=0``.
    * **G2** (exactly one strand free, with data): the STRAND DECONVOLUTION alone resolves the pie; the active
      strand axis shares ``f_g``'s posterior variance (a single-strand node is 1-D: ``f_active = 1 − f_g``).
    * **G3** (both strands free — AMBIG): unresolvable by strand → keep the ``{0,0,1}`` default at MAX
      (``inf``) variance; the sweep resolves it from neighbour messages + the global prior.
    """
    n = free_pos.shape[0]
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g = np.ones(n)  # all-gDNA default (count plays NO role; ARCHITECTURE §3)
    var_g = np.full(n, np.inf)
    var_p = np.where(free_pos, np.inf, 0.0)
    var_n = np.where(free_neg, np.inf, 0.0)

    g1 = ~free_pos & ~free_neg
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2-active: take the strand-only solve (median f_g, mean f±, f_g posterior variance).
    fgv = np.asarray(deconv.gdna_frac_var, dtype=np.float64)
    f_g[g2_active] = np.asarray(deconv.gdna_frac, dtype=np.float64)[g2_active]
    f_pos[g2_active] = np.asarray(deconv.rna_pos_frac, dtype=np.float64)[g2_active]
    f_neg[g2_active] = np.asarray(deconv.rna_neg_frac, dtype=np.float64)[g2_active]
    var_g[g2_active] = fgv[g2_active]
    var_p[g2_active & free_pos] = fgv[g2_active & free_pos]
    var_n[g2_active & free_neg] = fgv[g2_active & free_neg]

    # G1 sink: lock the gDNA axis (the fractions are already the {0,0,1} default).
    var_g[g1] = 0.0
    return f_pos, f_neg, f_g, var_p, var_n, var_g


def _boundary_strand_stats(boundary_substrate, region_arrays):
    """Per-boundary strand-solve sufficient statistics + the CONTINUITY-gated node class.

    Count ``N`` (the BB strand power) = ONE side per strand — ``max(left, right)`` (the two sides are equal for
    a contiguous crossing; ``max`` de-dups the straddle and guards a terminal's off-edge zero side). Mass (the
    pie base) = both sides summed (conserved). The spliced floor is the one-sided motif-stranded sense spliced
    on that strand's exon flank (a clean exon↔intron transition), never summed across sides.

    The allow mask is the **transcript-structure CONTINUITY gate** (`rna_imputation_transcript_structure.md`;
    plan v3 §3): a strand-s unspliced crossing is nascent RNA only where strand s is present on BOTH flanks, so
    ``free_s = has_s(left) & has_s(right)``. This BLOCKS RNA at a TSS/TES (intergenic↔exon → neither strand
    continuous → a gDNA sink) and at a mixed exon↔AMBIG seam (the non-shared strand cannot cross) — distinct
    from the legacy flank-OR ``boundary_nodes.boundary_node_class`` (retired in P4).
    """
    bs = boundary_substrate
    lr = np.asarray(bs.left_region, dtype=np.int64)
    rr = np.asarray(bs.right_region, dtype=np.int64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)  # off-edge terminal flank → no bits
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    has_pos_l = (sig_l & _POS_BITS) != 0
    has_pos_r = (sig_r & _POS_BITS) != 0
    has_neg_l = (sig_l & _NEG_BITS) != 0
    has_neg_r = (sig_r & _NEG_BITS) != 0
    free_pos = has_pos_l & has_pos_r  # +strand continuous across the seam
    free_neg = has_neg_l & has_neg_r

    left, right = bs.left, bs.right
    u_pos = np.maximum(left.n_unspliced_pos, right.n_unspliced_pos).astype(np.float64)
    u_neg = np.maximum(left.n_unspliced_neg, right.n_unspliced_neg).astype(np.float64)
    mass_unspl = np.asarray(left.mass_unspliced, float) + np.asarray(right.mass_unspliced, float)
    mass_spliced = np.asarray(left.mass_spliced, float) + np.asarray(right.mass_spliced, float)
    ls = left.n_spliced_sense.astype(np.float64)
    rs = right.n_spliced_sense.astype(np.float64)

    def _floor(exon_bit, intron_bit):
        le = (sig_l & exon_bit) != 0
        li = (sig_l & intron_bit) != 0
        re = (sig_r & exon_bit) != 0
        ri = (sig_r & intron_bit) != 0
        side_l = le & ~li & ri & ~re  # exon on the left → spliced on the left side
        side_r = li & ~le & re & ~ri
        return np.where(side_l, ls, np.where(side_r, rs, 0.0))

    spliced_pos = _floor(BIT_EXON_POS, BIT_INTRON_POS)
    spliced_neg = _floor(BIT_EXON_NEG, BIT_INTRON_NEG)
    return free_pos, free_neg, u_pos, u_neg, spliced_pos, spliced_neg, mass_unspl, mass_spliced


def init_beliefs(
    chain: NodeChain,
    substrate,
    boundary_substrate,
    region_arrays,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
) -> NodeBelief:
    """The signature-binary G1/G2/G3 initial :class:`NodeBelief` on the unified chain (plan v3 §6).

    Region and boundary nodes are classified + strand-solved by their own type, then gathered onto the chain.
    The strand solve is the bare strand likelihood (+ the Jeffreys reference at single-strand nodes) with NO
    global prior and NO imputation — those enter the sweep (P2/P3). Single-strand introns resolve to
    ``f_g≈0`` from the BB tilt alone (the zero-gDNA gate); intergenic / TSS sinks lock at ``{0,0,1}``; AMBIG
    nodes hold ``{0,0,1}`` at MAX variance for the sweep to resolve.
    """
    kw = dict(
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_grid=n_grid,
    )

    # --- region nodes: the strand-only solve (deconv_regions_sweep, no global/imputation) ---
    ts = np.asarray(region_arrays.strand_class)
    reg_free_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    reg_free_neg = (ts == TS_NEG) | (ts == TS_AMBIG)
    reg_deconv = deconv_regions_sweep(substrate, region_arrays, **kw)
    reg_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    r_fp, r_fn, r_fg, r_vp, r_vn, r_vg = _type_belief(reg_free_pos, reg_free_neg, reg_deconv, reg_mass)

    # --- boundary nodes: continuity-gated strand-only solve via the shared per-node core ---
    b_free_pos, b_free_neg, bu_pos, bu_neg, bspl_p, bspl_n, b_mass, b_mass_spl = _boundary_strand_stats(
        boundary_substrate, region_arrays
    )
    b_strand_obs = b_free_pos ^ b_free_neg
    bnd_deconv = _solve_nodes(
        bu_pos, bu_neg, bspl_p, bspl_n, b_free_pos, b_free_neg, b_strand_obs, b_mass, b_mass_spl,
        kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
        n_grid=n_grid,
    )
    b_fp, b_fn, b_fg, b_vp, b_vn, b_vg = _type_belief(b_free_pos, b_free_neg, bnd_deconv, b_mass)

    # --- gather onto the chain (region nodes ← region arrays; boundary nodes ← boundary arrays) ---
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    R = r_fg.shape[0]
    B = b_fg.shape[0]
    ri_ = np.clip(idx, 0, R - 1)
    bi_ = np.clip(idx, 0, B - 1)

    def pick(rv, bv):
        return np.where(is_reg, rv[ri_], np.where(is_bnd, bv[bi_], 0.0))

    return NodeBelief(
        f_pos=pick(r_fp, b_fp), f_neg=pick(r_fn, b_fn), f_g=pick(r_fg, b_fg),
        var_pos=pick(r_vp, b_vp), var_neg=pick(r_vn, b_vn), var_gdna=pick(r_vg, b_vg),
    )
