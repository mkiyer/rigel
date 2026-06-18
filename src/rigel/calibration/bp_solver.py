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
from .simplex_sweep import (
    _axis_mean,
    _fg_median,
    _fg_var,
    _local_loglik,
    _simplex_lattice,
    _solve_nodes,
)
from .strand_deconv import NodeDeconv
from .variance_model import MonotoneVarMean

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "NodeBelief",
    "NodeDensities",
    "node_densities",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
    "global_gdna_prior",
    "fit_gdna_varmean",
    "fit_rna_varmean",
    "node_sweep",
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9
_POS_BITS = BIT_EXON_POS | BIT_INTRON_POS
_NEG_BITS = BIT_EXON_NEG | BIT_INTRON_NEG
_MAD_TO_SIGMA = 1.4826  # normal-consistency constant: σ ≈ 1.4826·MAD (a mathematical constant, not a knob)


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
    continuous → a gDNA sink) and at a mixed exon↔AMBIG seam (the non-shared strand cannot cross).
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


def _region_strand_stats(substrate, region_arrays):
    """Per-region strand-solve sufficient statistics + node class (the region twin of
    :func:`_boundary_strand_stats`). A region's admissible strand axes are its own ±transcript bits
    (``free_s`` from the strand class); its sided spliced floor is the contained sense spliced on a
    single-strand region (0 on AMBIG — AMBIG spliced is not orientable)."""
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    free_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    free_neg = (ts == TS_NEG) | (ts == TS_AMBIG)
    spl = c.n_spliced_sense.astype(np.float64)
    spliced_pos = np.where(ts == TS_POS, spl, 0.0)
    spliced_neg = np.where(ts == TS_NEG, spl, 0.0)
    mass_u = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_s = np.asarray(c.mass_spliced, dtype=np.float64)
    return free_pos, free_neg, u_pos, u_neg, spliced_pos, spliced_neg, mass_u, mass_s


@dataclass(frozen=True, slots=True)
class NodeStatics:
    """Per-chain-node STATIC strand-solve sufficient statistics + node class (length ``n_nodes``) — the
    region- and boundary-keyed stats gathered onto the chain once. The sweep mutates only the dynamic
    :class:`NodeBelief`; these never change. ``free_pos``/``free_neg`` are the admissible strand axes (a
    region's ±bits / a boundary's ±continuity); ``strand_obs = free_pos ^ free_neg`` (a single consistent
    strand). All ``float64`` except the three bool masks."""

    n_nodes: int
    u_pos: np.ndarray
    u_neg: np.ndarray
    spliced_pos: np.ndarray
    spliced_neg: np.ndarray
    free_pos: np.ndarray  # bool
    free_neg: np.ndarray  # bool
    strand_obs: np.ndarray  # bool
    mass_unspliced: np.ndarray
    mass_spliced: np.ndarray


def build_node_statics(chain: NodeChain, substrate, boundary_substrate, region_arrays) -> NodeStatics:
    """Gather the per-region (contained) and per-boundary (continuity-gated, max-of-sides + spliced floor)
    strand-solve statistics onto the unified chain."""
    r_fp, r_fn, r_up, r_un, r_sp, r_sn, r_mu, r_ms = _region_strand_stats(substrate, region_arrays)
    b_fp, b_fn, b_up, b_un, b_sp, b_sn, b_mu, b_ms = _boundary_strand_stats(
        boundary_substrate, region_arrays
    )
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    R = r_fp.shape[0]
    B = b_fp.shape[0]
    ri_ = np.clip(idx, 0, R - 1)
    bi_ = np.clip(idx, 0, B - 1)

    def pick(rv, bv, fill=0.0):
        return np.where(is_reg, rv[ri_], np.where(is_bnd, bv[bi_], fill))

    free_pos = pick(r_fp, b_fp, False)
    free_neg = pick(r_fn, b_fn, False)
    return NodeStatics(
        n_nodes=int(chain.n_nodes),
        u_pos=pick(r_up, b_up), u_neg=pick(r_un, b_un),
        spliced_pos=pick(r_sp, b_sp), spliced_neg=pick(r_sn, b_sn),
        free_pos=free_pos, free_neg=free_neg, strand_obs=free_pos ^ free_neg,
        mass_unspliced=pick(r_mu, b_mu), mass_spliced=pick(r_ms, b_ms),
    )


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
    statics: "NodeStatics | None" = None,
) -> NodeBelief:
    """The signature-binary G1/G2/G3 initial :class:`NodeBelief` on the unified chain (plan v3 §6).

    All chain nodes are strand-solved by one shared per-node core (the bare strand likelihood + the Jeffreys
    reference at single-strand nodes; NO global prior, NO imputation — those enter the sweep, P2/P3). The
    signature-binary class overrides (:func:`_type_belief`) then set the G1/G2/G3 belief. Single-strand
    introns resolve to ``f_g≈0`` from the BB tilt alone (the zero-gDNA gate); intergenic / TSS sinks lock at
    ``{0,0,1}``; AMBIG nodes hold ``{0,0,1}`` at MAX variance for the sweep. Pass ``statics`` to reuse a
    prebuilt :class:`NodeStatics` (the sweep does)."""
    st = statics if statics is not None else build_node_statics(
        chain, substrate, boundary_substrate, region_arrays
    )
    deconv = _solve_nodes(
        st.u_pos, st.u_neg, st.spliced_pos, st.spliced_neg, st.free_pos, st.free_neg, st.strand_obs,
        st.mass_unspliced, st.mass_spliced,
        kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
        n_grid=n_grid,
    )
    f_pos, f_neg, f_g, var_p, var_n, var_g = _type_belief(
        st.free_pos, st.free_neg, deconv, st.mass_unspliced
    )
    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_p, var_neg=var_n, var_gdna=var_g
    )


def global_gdna_prior(region_f_g, region_mass_unspliced, region_eff_gdna, region_count_observable):
    """The population gDNA baseline ``ρ_global`` + its robust density spread ``σ_global`` (plan v3 §6 / §4 step 3).

    ``ρ_global`` is the gDNA density over the **count-observable** (signature-known non-exonic: intergenic /
    intron) regions' current gDNA mass ``f_g·M_u`` — the cleanest gDNA evidence, restricted to avoid the
    AMBIG/exon RNA the strand cannot remove inflating the baseline (`calibrate.py`'s carried-forward rule).
    ``σ_global`` is the robust between-region MAD spread of the per-region gDNA density over ALL regions with
    data (the population spread an unanchored node faces). Both are recomputed each sweep pass on the running
    gDNA estimate. ``ρ_global → 0`` (and ``σ_global → 0``) in a pure-RNA library — the zero-gDNA pin, no
    phantom gDNA. The per-node application (the density→fraction jacobian ``(E/M)²`` + the 1-pseudo-observation
    precision floor) lives in the sweep, which holds each node's own geometry.
    """
    f_g = np.asarray(region_f_g, dtype=np.float64)
    mass = np.asarray(region_mass_unspliced, dtype=np.float64)
    eff = np.maximum(np.asarray(region_eff_gdna, dtype=np.float64), _EPS)
    obs = np.asarray(region_count_observable, dtype=bool)
    gdna_mass = f_g * mass
    rho_global = (
        float(gdna_mass[obs].sum() / max(float(eff[obs].sum()), _EPS)) if obs.any() else 0.0
    )
    mu = gdna_mass / eff  # per-region gDNA density
    active = mu[mass > 0.0]
    sigma_global = (
        _MAD_TO_SIGMA * float(np.median(np.abs(active - np.median(active)))) if active.size else rho_global
    )
    return rho_global, sigma_global


# ---------------------------------------------------------------------------
# P2 — the gDNA-only directional belief-propagation sweep (plan v3 §3/§4/§5/§7).
# ---------------------------------------------------------------------------

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG


def _edge_varmean(chain, rho_left, rho_right, eff_left, eff_right, live):
    """Poisson-offset var~mean points over directed chain edges where BOTH endpoints are ``live`` (a bool mask
    — all nodes for gDNA; per-strand ``free_s`` continuity for RNA). For each directed edge the dst presents
    its facing density/eff and the src the side facing the dst; point = ``(μ=ρ(dst), raw=(ρ(dst)−ρ(src))²,
    offset=ρ(dst)/E(dst)+ρ(src)/E(src))``. Returns ``(means, raws, offs)`` lists (both sweep directions)."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    live = np.asarray(live, dtype=bool)
    means, raws, offs = [], [], []
    for nbr, d_rho, s_rho, d_e, s_e in (
        (left, rho_left, rho_right, eff_left, eff_right),  # from-LEFT: dst left face, src=left[dst] right face
        (right, rho_right, rho_left, eff_right, eff_left),  # from-RIGHT: dst right face, src=right[dst] left
    ):
        idx = np.where((nbr >= 0) & live)[0]
        if not idx.size:
            continue
        s = nbr[idx]
        keep = live[s]  # the source must also be live on this component/strand
        idx, s = idx[keep], s[keep]
        dr, sr, de, se = d_rho[idx], s_rho[s], d_e[idx], s_e[s]
        ok = (dr > 0.0) & (sr > 0.0)
        means.append(dr[ok])
        raws.append((dr[ok] - sr[ok]) ** 2)
        offs.append(dr[ok] / de[ok] + sr[ok] / se[ok])
    return means, raws, offs


def _fit_offset(means, raws, offs):
    cat = lambda parts: np.concatenate(parts) if parts else np.zeros(0)  # noqa: E731
    return MonotoneVarMean.fit_offset(cat(means), cat(raws), cat(offs))


def fit_gdna_varmean(chain: NodeChain, densities: NodeDensities, geometry: NodeGeometry) -> MonotoneVarMean:
    """Fit the gDNA message reliability ``σ²_bio(μ)`` on the FROZEN snapshot densities (plan v3 §7). gDNA flows
    genomically (every edge is ``live``); the ``(dr>0)&(sr>0)`` filter drops pure-RNA pairs. Fitting on the
    frozen previous-pass densities (not in-place) keeps σ²_bio from collapsing to 0 (plan v3 §4)."""
    live = np.ones(int(chain.n_nodes), dtype=bool)
    m, r, o = _edge_varmean(
        chain, densities.rho_g_left, densities.rho_g_right,
        geometry.eff_gdna_left, geometry.eff_gdna_right, live,
    )
    return _fit_offset(m, r, o)


def fit_rna_varmean(
    chain: NodeChain, densities: NodeDensities, geometry: NodeGeometry, statics: NodeStatics
) -> MonotoneVarMean:
    """Fit the RNA message reliability ``σ²_bio(μ)`` on the FROZEN snapshot, POOLING both strands (the symmetric
    RNA process; plan v3 §7). A strand-s edge is ``live`` only where strand s is continuous on BOTH endpoints
    (``free_s`` — the transcript-structure gate), so the curve sees the genuine same-strand cross-node RNA
    dispersion (INCLUDING the AMBIG nodes' per-strand densities — the honest-AMBIG-dispersion the chain now
    provides, the 2c fix). The per-strand RNA density is spliced-inclusive (``node_densities``)."""
    mp, rp, op = _edge_varmean(
        chain, densities.rho_pos_left, densities.rho_pos_right,
        geometry.eff_rna_left, geometry.eff_rna_right, statics.free_pos,
    )
    mn, rn, on = _edge_varmean(
        chain, densities.rho_neg_left, densities.rho_neg_right,
        geometry.eff_rna_left, geometry.eff_rna_right, statics.free_neg,
    )
    return _fit_offset(mp + mn, rp + rn, op + on)


def _message(rho_src, eff_src, eff_dst, mass_dst, rho_dst_cur, varmean, spliced_dst=0.0):
    """One density message (source → dest) as the dest's fraction prior ``(μ_f, τ_f)`` (plan v3 §3/§5). Mean =
    the IDENTITY density ``ρ_src`` matched at the dest: its known spliced (free info, gDNA: 0) is subtracted so
    only the dest's UNSOLVED unspliced fraction is informed — ``μ_f = (ρ_src·E_dst − spliced_dst)/M_dst``.
    Density precision ``τ_ρ = 1/(σ²_bio(ρ_dst) + ρ_src/E_src)`` (the learned biological dispersion + the
    predictor's Poisson sampling noise); jacobian to fraction precision ``τ_f = τ_ρ·(M_dst/E_dst)²``."""
    mu_f = float(np.clip((rho_src * eff_dst - spliced_dst) / max(mass_dst, _EPS), 0.0, 1.0))
    sigma2_bio = float(varmean.predict(np.array([max(rho_dst_cur, _EPS)]))[0])
    tau_rho = 1.0 / max(sigma2_bio + rho_src / max(eff_src, _EPS), _EPS)
    tau_f = tau_rho * (mass_dst / max(eff_dst, _EPS)) ** 2
    return mu_f, tau_f


def node_sweep(
    chain: NodeChain,
    statics: NodeStatics,
    geometry: NodeGeometry,
    belief: NodeBelief,
    region_arrays,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
    max_passes: int,
    convergence_delta: float,
):
    """The directional (Gauss-Seidel) belief-propagation sweep over the chain — gDNA AND per-strand RNA
    messages (plan v3 §3/§4/§5/§7).

    Each outer pass: (1) freeze the snapshot densities + fit the gDNA and RNA ``σ²_bio`` var~means on them;
    (2) recompute the global gDNA prior on the running region beliefs; (3) sweep **L→R** then **R→L**, each node
    pulled by ONE neighbour's message in the sweep direction (its CURRENT — Gauss-Seidel — density). The per-node
    ψ solve integrates: the strand likelihood, the node-class prior (Jeffreys / global), the gDNA imputation
    prior on ``f_g``, and the per-strand RNA imputation priors on ``f±``. The RNA message on strand s flows only
    where strand s is continuous on both endpoints (``free_s``) — so a strand's variance SINK (a TES / a non-s
    flank, ``free_s=False``) blocks it structurally, no extra gate (the bidirectional sweep + the per-strand
    lock state encode strandedness + transcript structure). gDNA and RNA use the SAME machinery; the boundary's
    spliced is free info that rides in its RNA density (``node_densities``) + the dest-spliced subtraction.
    Only G2/G3 nodes with data are solved; G1 sinks + empty nodes are fixed. Returns ``(NodeBelief, deltas)``."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    order = np.asarray(chain.order)
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    fp, fn = statics.free_pos, statics.free_neg
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()

    # per-face geometry as (left, right) pairs, indexed by face (0=left, 1=right).
    EG = (geometry.eff_gdna_left, geometry.eff_gdna_right)
    ER = (geometry.eff_rna_left, geometry.eff_rna_right)
    MS = (geometry.mass_left, geometry.mass_right)
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # per-node "global" gDNA support: region = its contained support; boundary = the shipped AVERAGED
    # per-side density length ½(E_l+E_r) (plan v3 §2) over the total crossing mass.
    eff_global = np.where(is_reg, EG[0], 0.5 * (EG[0] + EG[1]))
    mass_global = np.where(is_reg, MS[0], MS[0] + MS[1])
    geom2_global = (eff_global / np.maximum(mass_global, _EPS)) ** 2

    # region-keyed arrays for the global prior (scatter the region chain nodes → region index).
    reg_nodes = order[is_reg]
    reg_index = np.asarray(chain.ref_idx, dtype=np.int64)[reg_nodes]
    R = int(reg_index.max()) + 1 if reg_index.size else 0
    region_eff_gdna = np.zeros(R)
    region_eff_gdna[reg_index] = EG[0][reg_nodes]
    region_mass = np.zeros(R)
    region_mass[reg_index] = MS[0][reg_nodes]
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    region_count_obs = (sig & _EXON_BITS) == 0

    lattice = _simplex_lattice(int(n_grid))
    fpg, fng, fgg = lattice
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)

    def _solve(i, mu_g, tau_g, mu_p, tau_p, mu_n, tau_n, mu_glob, tau_glob):
        psi = _local_loglik(
            statics.u_pos[i:i + 1], statics.u_neg[i:i + 1],
            statics.spliced_pos[i:i + 1], statics.spliced_neg[i:i + 1],
            fp[i:i + 1], fn[i:i + 1], kappa, od_g, od_r, lattice,
            strand_obs=statics.strand_obs[i:i + 1],
            global_mu=np.array([mu_glob]), global_tau=np.array([tau_glob]),
            gdna_imp_frac=np.array([mu_g]), gdna_imp_precision=np.array([tau_g]),
            rna_imp_frac=(np.array([mu_p]), np.array([mu_n])),
            rna_imp_precision=(np.array([tau_p]), np.array([tau_n])),
        )
        f_g[i] = float(np.clip(_fg_median(psi, fgg)[0], 0.0, 1.0))
        f_pos[i] = float(np.clip(_axis_mean(psi, fpg)[0], 0.0, 1.0))
        f_neg[i] = float(np.clip(_axis_mean(psi, fng)[0], 0.0, 1.0))
        var_g[i] = float(_fg_var(psi, fgg)[0])

    deltas = []
    for _pass in range(int(max_passes)):
        snap = node_densities(
            NodeBelief(f_pos, f_neg, f_g, belief.var_pos, belief.var_neg, var_g), geometry
        )
        SNG = (snap.rho_g_left, snap.rho_g_right)
        SNP = (snap.rho_pos_left, snap.rho_pos_right)
        SNN = (snap.rho_neg_left, snap.rho_neg_right)
        gdna_vm = fit_gdna_varmean(chain, snap, geometry)
        rna_vm = fit_rna_varmean(chain, snap, geometry, statics)
        region_f_g = np.zeros(R)
        region_f_g[reg_index] = f_g[reg_nodes]
        rho_global, sigma_global = global_gdna_prior(
            region_f_g, region_mass, region_eff_gdna, region_count_obs
        )
        mu_global = np.clip(rho_global * eff_global / np.maximum(mass_global, _EPS), 0.0, 1.0)
        tau_global = np.maximum(1.0 / np.maximum(sigma_global**2 * geom2_global, _EPS), 1.0)

        prev = (f_g.copy(), f_pos.copy(), f_neg.copy())
        # df = the dst's face toward its sweep-direction neighbour; sf = the src's face toward the dst.
        for nbr_arr, df, sf in ((left, 0, 1), (right, 1, 0)):  # L→R (from left) then R→L (from right)
            for i in (order if df == 0 else order[::-1]):
                if not solvable[i]:
                    continue
                src = nbr_arr[i]
                if src < 0 or MS[sf][src] <= _EPS:
                    continue  # reference terminal / empty source → no message this direction
                s_mass = MS[sf][src]
                # gDNA message (genomic; spliced_dst=0).
                rho_g_src = f_g[src] * s_mass / max(EG[sf][src], _EPS)
                mu_g, tau_g = _message(rho_g_src, EG[sf][src], EG[df][i], MS[df][i], SNG[df][i], gdna_vm)
                # per-strand RNA messages, gated by free_s on BOTH endpoints; spliced rides in ρ_s + the
                # dest-spliced subtraction. τ=0 (no message) where the strand is not continuous.
                mu_p = tau_p = mu_n = tau_n = 0.0
                if fp[i] and fp[src]:
                    rho_p_src = (f_pos[src] * s_mass + SP[sf][src]) / max(ER[sf][src], _EPS)
                    mu_p, tau_p = _message(
                        rho_p_src, ER[sf][src], ER[df][i], MS[df][i], SNP[df][i], rna_vm, SP[df][i]
                    )
                if fn[i] and fn[src]:
                    rho_n_src = (f_neg[src] * s_mass + SN[sf][src]) / max(ER[sf][src], _EPS)
                    mu_n, tau_n = _message(
                        rho_n_src, ER[sf][src], ER[df][i], MS[df][i], SNN[df][i], rna_vm, SN[df][i]
                    )
                _solve(i, mu_g, tau_g, mu_p, tau_p, mu_n, tau_n, mu_global[i], tau_global[i])
        delta = max(
            (float(np.max(np.abs(cur - p))) if cur.size else 0.0)
            for cur, p in ((f_g, prev[0]), (f_pos, prev[1]), (f_neg, prev[2]))
        )
        deltas.append(delta)
        if delta < convergence_delta:
            break

    return (
        NodeBelief(f_pos=f_pos, f_neg=f_neg, f_g=f_g,
                   var_pos=belief.var_pos, var_neg=belief.var_neg, var_gdna=var_g),
        deltas,
    )


def chain_region_deconv(chain: NodeChain, belief: NodeBelief, substrate) -> NodeDeconv:
    """Project the chain belief's REGION nodes back to a region-keyed :class:`NodeDeconv` — the transitional
    region projection the existing ``CalibrationResult`` / ``priors`` / ``derive`` consume (the per-node
    first-class schema rewire is P4). gDNA / RNA masses from each region's solved ``f_g`` over its contained
    unspliced (+ the always-RNA contained spliced) mass."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    reg = kind == REGION
    mass_u = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    mass_s = np.asarray(substrate.contained.mass_spliced, dtype=np.float64)
    R = mass_u.shape[0]
    f_g = np.zeros(R)
    f_pos = np.zeros(R)
    f_neg = np.zeros(R)
    f_gv = np.zeros(R)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    f_gv[ri] = belief.var_gdna[reg]
    return NodeDeconv(
        gdna_mass=f_g * mass_u, rna_mass=(1.0 - f_g) * mass_u + mass_s,
        gdna_frac=f_g, gdna_frac_var=f_gv, rna_pos_frac=f_pos, rna_neg_frac=f_neg,
    )


def chain_boundary_side_deconv(chain: NodeChain, belief: NodeBelief, substrate):
    """Project the chain belief's BOUNDARY ``f_g`` onto each region's two SIDE views — the boundary-flux that
    ``priors``' pooled-seam gDNA eff-len + ``derive`` consume.

    Region ``r``'s left/right boundary IS its left/right chain neighbour; that boundary's solved ``f_g`` splits
    ``r``'s side crossing mass (``substrate.left[r]`` / ``substrate.right[r]`` — the boundary flux already
    projected onto the region by the D1 side-attribution) into gDNA / RNA (the RNA spliced-inclusive, matching
    the contained projection). One boundary pie applied to its side mass. Returns ``(left, right)`` region-keyed
    :class:`NodeDeconv`."""
    kind = np.asarray(chain.kind)
    reg_nodes = np.asarray(chain.order)[kind == REGION]
    ridx = np.asarray(chain.ref_idx, dtype=np.int64)[reg_nodes]
    R = int(ridx.max()) + 1 if ridx.size else 0
    fg = np.asarray(belief.f_g, dtype=np.float64)
    left_fg = np.zeros(R)
    right_fg = np.zeros(R)
    left_fg[ridx] = fg[np.asarray(chain.left)[reg_nodes]]  # f_g of r's left-flank boundary node
    right_fg[ridx] = fg[np.asarray(chain.right)[reg_nodes]]

    def _side(side_fg, view):
        m_u = np.asarray(view.mass_unspliced, dtype=np.float64)
        m_s = np.asarray(view.mass_spliced, dtype=np.float64)
        return NodeDeconv(
            gdna_mass=side_fg * m_u, rna_mass=(1.0 - side_fg) * m_u + m_s,
            gdna_frac=side_fg, gdna_frac_var=np.zeros(R),
            rna_pos_frac=np.zeros(R), rna_neg_frac=np.zeros(R),
        )

    return _side(left_fg, substrate.left), _side(right_fg, substrate.right)
