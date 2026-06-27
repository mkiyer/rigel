"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a pairwise, directional
(L→R then R→L) Gauss-Seidel sweep. Each per-node solve (`simplex_sweep`) reconciles three sources of
information: the intrinsic strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the
cross-node imputation messages, and the population gDNA prior. Theory: the count-zero-information principle
in `docs/calibration/CALIBRATION_ARCHITECTURE.md`; the message-precision model (count-currency variance
floor) in `docs/calibration/precision_overshoot_design.md`.

Module layout:
* `build_node_geometry` — the per-node, per-FACE static geometry (unspliced mass, per-component gDNA-/RNA-FL
  effective lengths, the boundary's one-sided motif-stranded spliced mass). Per-face because a boundary's
  two sides lie in different-sized regions and a message uses the side FACING the destination; a region
  presents the same contained geometry both ways. Eff-len primitives are the single-owner
  `effective_length.py` (`region_eff_length` = E[max(0,L−ℓ)] contained; `boundary_side_eff_length` =
  E[min(ℓ,L)] per-side crossing). gDNA uses the gDNA FL; RNA (nascent unspliced + spliced) the RNA FL.
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.
* `node_densities` / `build_node_statics` — the per-pass density snapshot + the static per-node inputs.
* `fit_rna_varmean` — the RNA message-reliability σ²_bio(μ), refit each pass on the FROZEN snapshot (RNA is
  well-resolved → benign; the gDNA σ²_g is the seed-based non-circular firewall from `_gdna_seed_estimate`).
* `_message` / `node_sweep` — one imputation message (as a count-currency precision) + the sweep driver.
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .effective_length import (
    boundary_side_eff_length,
    region_eff_length,
    spliced_side_eff_length,
)
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
    _binom_pseudo,
    _local_loglik,
    _node_marginals,
    _simplex_lattice,
    _solve_nodes,
)
from .strand_deconv import NodeDeconv
from .variance_model import MonotoneMean, MonotoneVarMean

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "NodeBelief",
    "NodeDensities",
    "node_densities",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
    "fit_rna_varmean",
    "fit_gdna_varmean",
    "fit_enrichment_transfer",
    "node_sweep",
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9
_MSG_PSEUDOCOUNT = 1.0  # one pseudo-observation: a density from a finite count can never have zero
#                         sampling variance, so a message precision can never escape to ∞ (var stabilizer).
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
    eff_gdna_left: np.ndarray  # gDNA-FL eff-len facing left (contained / boundary left side)
    eff_gdna_right: np.ndarray
    eff_rna_left: np.ndarray  # RNA-FL eff-len facing left (nascent two-sided crossing / contained)
    eff_rna_right: np.ndarray
    eff_spl_left: np.ndarray  # one-sided spliced-crossing RNA eff-len (half-triangle E[min²/2ℓ]) facing left
    eff_spl_right: np.ndarray
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
    """Assemble the per-node per-face geometry from the region-/boundary-keyed substrates onto the chain.

    The one-sided spliced (mature) mass is routed to its exon flank by the boundary's KNOWN junction strand
    (``boundary_substrate.junction_strand``, observed from the splice motif at deposit — ``TS_POS``/``TS_NEG``/0,
    one motif-stranded junction per boundary). This is correct even at AMBIG / exon↔exon seams, where the
    flanking region signatures cannot orient the junction."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY

    # --- per-region quantities (contained) ---
    L = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    reg_eff_g = region_eff_length(L, gdna_fl_pmf)  # E_gdna[max(0,L−ℓ)] contained
    reg_eff_r = region_eff_length(L, rna_fl_pmf)
    side_eff_g = boundary_side_eff_length(
        gdna_fl_pmf, L
    )  # per region: E_gdna[min(ℓ,L)] side crossing
    side_eff_r = boundary_side_eff_length(rna_fl_pmf, L)
    side_eff_spl = spliced_side_eff_length(rna_fl_pmf, L)  # one-sided spliced-crossing half-triangle
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
    b_eff_spl_l = np.where(blr >= 0, side_eff_spl[np.clip(blr, 0, R - 1)], 0.0)
    b_eff_spl_r = np.where(brr >= 0, side_eff_spl[np.clip(brr, 0, R - 1)], 0.0)
    # spliced face: route a boundary side's spliced (mature) mass to that side's EXON flank on the junction's
    # KNOWN genomic strand (``boundary_substrate.junction_strand``, observed from the motif at deposit — one
    # motif-stranded junction per boundary). Correct at AMBIG / exon↔exon seams the signatures cannot orient.
    # One-sided per exon: mass_left → left-exon donor, mass_right → right-exon acceptor.
    sig_l = np.where(blr >= 0, sig[np.clip(blr, 0, R - 1)], 0)
    sig_r = np.where(brr >= 0, sig[np.clip(brr, 0, R - 1)], 0)
    js = np.asarray(boundary_substrate.junction_strand).reshape(-1)  # TS_POS / TS_NEG / 0 (== Strand)
    any_exon_l = (sig_l & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    any_exon_r = (sig_r & (BIT_EXON_POS | BIT_EXON_NEG)) != 0

    def _spliced_faces(strand_val):
        on = js == strand_val  # the junction is on this strand → its exon flank carries the mature floor
        return np.where(on & any_exon_l, bspl_l, 0.0), np.where(on & any_exon_r, bspl_r, 0.0)

    b_spl_pos_l, b_spl_pos_r = _spliced_faces(TS_POS)
    b_spl_neg_l, b_spl_neg_r = _spliced_faces(TS_NEG)

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
    # one-sided spliced-crossing eff-len: boundaries get the half-triangle of their exon flank; regions
    # carry no spliced mass (spliced_* = 0) so their value is inert — use the per-region half-triangle.
    eff_spl_left = pick(side_eff_spl, b_eff_spl_l)
    eff_spl_right = pick(side_eff_spl, b_eff_spl_r)
    zeros_R = np.zeros(R)
    spliced_pos_left = pick(zeros_R, b_spl_pos_l)
    spliced_pos_right = pick(zeros_R, b_spl_pos_r)
    spliced_neg_left = pick(zeros_R, b_spl_neg_l)
    spliced_neg_right = pick(zeros_R, b_spl_neg_r)

    return NodeGeometry(
        n_nodes=int(chain.n_nodes),
        mass_left=mass_left,
        mass_right=mass_right,
        eff_gdna_left=np.maximum(eff_gdna_left, _EPS),
        eff_gdna_right=np.maximum(eff_gdna_right, _EPS),
        eff_rna_left=np.maximum(eff_rna_left, _EPS),
        eff_rna_right=np.maximum(eff_rna_right, _EPS),
        eff_spl_left=np.maximum(eff_spl_left, _EPS),
        eff_spl_right=np.maximum(eff_spl_right, _EPS),
        spliced_pos_left=spliced_pos_left,
        spliced_pos_right=spliced_pos_right,
        spliced_neg_left=spliced_neg_left,
        spliced_neg_right=spliced_neg_right,
    )


@dataclass(frozen=True, slots=True)
class NodeBelief:
    """Per-node solved state on the chain: the composition pie `(f_pos, f_neg, f_g)` over the node's UNSPLICED
    mass + its per-component posterior VARIANCE `(var_pos, var_neg, var_gdna)` = `Var(f_c)`. All length
    ``n_nodes``.

    The variance is the **precision state** (`docs/calibration/precision_state_design.md`): `Var(f_c)=0` ⇒
    locked/certain (e.g. a forbidden strand), `=∞` ⇒ no information (unsolved). It feeds the honest message
    send — a source's outgoing precision is degraded from its own `Var_own = (M/E)²·Var(f_c)` by the
    communication noise, so an unsure node speaks quietly (Phase 2). The composition is stored as a FRACTION
    (the face-invariant quantity — a boundary has two faces but one composition); density `ρ=f·M_face/E_face`
    (`node_densities`) is the message currency, mass `m=f·M_face` (`NodeDeconv`) the output."""

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
    # The RNA density sums two pieces with DIFFERENT crossing geometry: the nascent unspliced (``f·M/eff_rna``,
    # the two-sided ``E[min(ℓ,R)]`` — contiguous crossing) plus the one-sided spliced/mature floor
    # (``spliced/eff_spl``, the half-triangle ``E[min²/2ℓ]``). Sharing one divisor understates the one-sided
    # spliced ~2× (the mature density bug).
    return NodeDensities(
        rho_g_left=fg * g.mass_left / g.eff_gdna_left,
        rho_pos_left=fp * g.mass_left / g.eff_rna_left + g.spliced_pos_left / g.eff_spl_left,
        rho_neg_left=fn * g.mass_left / g.eff_rna_left + g.spliced_neg_left / g.eff_spl_left,
        rho_g_right=fg * g.mass_right / g.eff_gdna_right,
        rho_pos_right=fp * g.mass_right / g.eff_rna_right + g.spliced_pos_right / g.eff_spl_right,
        rho_neg_right=fn * g.mass_right / g.eff_rna_right + g.spliced_neg_right / g.eff_spl_right,
    )


# ---------------------------------------------------------------------------
# Initialization — the signature-binary G1/G2/G3 belief on the chain.
# ---------------------------------------------------------------------------
#
# A strand axis is hard-LOCKED (a forbidden strand, an intergenic sink) by the per-node ``allow_pos`` /
# ``allow_neg`` forbid mask in the solve. The init ALSO sets the per-component precision state ``var(f_c)``
# (`precision_state_design.md`): ``0`` = locked/certain (a forbidden strand, an intergenic gDNA sink), ``inf`` =
# no information (an admissible-but-unsolved axis — it will listen to messages, and emits none until solved).
# A solved single-strand (G2) node takes the strand-solve posterior variance.


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """Build the per-node composition ``(f_pos, f_neg, f_g)`` for ONE node type (regions OR boundaries) from its
    signature-binary classification + its strand-only solve.

    ``free_pos``/``free_neg`` are the per-node booleans for whether each strand's RNA axis is admissible (a
    region's own ±transcript bits; a boundary's ±strand CONTINUITY across the seam). ``deconv`` is the
    strand-only :class:`NodeDeconv` (no global, no imputation). The signature-binary default is all-gDNA
    ``{0,0,1}`` (`ARCHITECTURE §3`). The class overrides:

    * **G1** (neither strand free — intergenic region / no-RNA-crossing boundary): a LOCKED gDNA sink — keep
      ``{0,0,1}``.
    * **G2** (exactly one strand free, with data): the STRAND DECONVOLUTION alone resolves the pie (a
      single-strand node is 1-D: ``f_active = 1 − f_g``).
    * **G3** (both strands free — AMBIG): unresolvable by strand → keep the ``{0,0,1}`` default at MAX
      (``inf``) variance; the sweep resolves it from neighbour messages + the global prior.

    Returns the six per-node arrays ``(f_pos, f_neg, f_g, var_pos, var_neg, var_gdna)`` — the composition + the
    precision state (`precision_state_design.md`): ``var=0`` locked, ``inf`` no-information, else the strand-solve
    posterior variance.
    """
    n = free_pos.shape[0]
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g = np.ones(n)  # all-gDNA default (count plays NO role; ARCHITECTURE §3)
    # precision state: gDNA unsolved (inf); a strand axis is locked (0) iff forbidden, else unsolved (inf).
    var_g = np.full(n, np.inf)
    var_p = np.where(free_pos, np.inf, 0.0)
    var_n = np.where(free_neg, np.inf, 0.0)

    g1 = ~free_pos & ~free_neg
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2-active: take the strand-only solve (median f_g, mean f±, and the posterior variances). G1 sinks + G3
    # AMBIG keep the {0,0,1} default at MAX variance.
    fgv = np.asarray(deconv.gdna_frac_var, dtype=np.float64)
    fpv = np.asarray(deconv.rna_pos_frac_var, dtype=np.float64)
    fnv = np.asarray(deconv.rna_neg_frac_var, dtype=np.float64)
    f_g[g2_active] = np.asarray(deconv.gdna_frac, dtype=np.float64)[g2_active]
    f_pos[g2_active] = np.asarray(deconv.rna_pos_frac, dtype=np.float64)[g2_active]
    f_neg[g2_active] = np.asarray(deconv.rna_neg_frac, dtype=np.float64)[g2_active]
    var_g[g2_active] = fgv[g2_active]
    var_p[g2_active & free_pos] = fpv[g2_active & free_pos]
    var_n[g2_active & free_neg] = fnv[g2_active & free_neg]

    # G1 sink: lock the gDNA axis (the fractions are already the {0,0,1} default).
    var_g[g1] = 0.0
    return f_pos, f_neg, f_g, var_p, var_n, var_g


def _boundary_strand_stats(boundary_substrate, region_arrays):
    """Per-boundary strand-solve sufficient statistics + the CONTINUITY-gated node class.

    Count ``N`` (the BB strand power) = ONE side per strand — ``max(left, right)`` (the two sides are equal for
    a contiguous crossing; ``max`` de-dups the straddle and guards a terminal's off-edge zero side). Mass (the
    pie base) = both sides summed (conserved). The spliced floor is the one-sided motif-stranded sense spliced
    on that strand's exon flank (a clean exon↔intron transition), never summed across sides.

    The allow mask is the **transcript-structure CONTINUITY gate** (`rna_imputation_transcript_structure.md`):
    a strand-s unspliced crossing is nascent RNA only where strand s is present on BOTH flanks, so
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


def build_node_statics(
    chain: NodeChain, substrate, boundary_substrate, region_arrays
) -> NodeStatics:
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
        u_pos=pick(r_up, b_up),
        u_neg=pick(r_un, b_un),
        spliced_pos=pick(r_sp, b_sp),
        spliced_neg=pick(r_sn, b_sn),
        free_pos=free_pos,
        free_neg=free_neg,
        strand_obs=free_pos ^ free_neg,
        mass_unspliced=pick(r_mu, b_mu),
        mass_spliced=pick(r_ms, b_ms),
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
    """The signature-binary G1/G2/G3 initial :class:`NodeBelief` on the unified chain.

    All chain nodes are strand-solved by one shared per-node core (the bare strand likelihood + the Jeffreys
    reference at single-strand nodes; NO global prior, NO imputation — those enter the sweep, P2/P3). The
    signature-binary class overrides (:func:`_type_belief`) then set the G1/G2/G3 belief. Single-strand
    introns resolve to ``f_g≈0`` from the BB tilt alone (the zero-gDNA gate); intergenic / TSS sinks lock at
    ``{0,0,1}``; AMBIG nodes hold ``{0,0,1}`` at MAX variance for the sweep. Pass ``statics`` to reuse a
    prebuilt :class:`NodeStatics` (the sweep does)."""
    st = (
        statics
        if statics is not None
        else build_node_statics(chain, substrate, boundary_substrate, region_arrays)
    )
    deconv = _solve_nodes(
        st.u_pos,
        st.u_neg,
        st.spliced_pos,
        st.spliced_neg,
        st.free_pos,
        st.free_neg,
        st.strand_obs,
        st.mass_unspliced,
        st.mass_spliced,
        kappa=float(rna_sense_frac),
        od_g=gdna_strand_overdispersion,
        od_r=rna_strand_overdispersion,
        n_grid=n_grid,
    )
    f_pos, f_neg, f_g, var_p, var_n, var_g = _type_belief(
        st.free_pos, st.free_neg, deconv, st.mass_unspliced
    )
    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_p, var_neg=var_n, var_gdna=var_g
    )


# ---------------------------------------------------------------------------
# The directional belief-propagation sweep (gDNA + per-strand RNA messages).
# ---------------------------------------------------------------------------


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
        (
            left,
            rho_left,
            rho_right,
            eff_left,
            eff_right,
        ),  # from-LEFT: dst left face, src=left[dst] right face
        (
            right,
            rho_right,
            rho_left,
            eff_right,
            eff_left,
        ),  # from-RIGHT: dst right face, src=right[dst] left
    ):
        idx = np.where((nbr >= 0) & live)[0]
        if not idx.size:
            continue
        s = nbr[idx]
        keep = live[s]  # the source must also be live on this component/strand
        idx, s = idx[keep], s[keep]
        dr, sr, de, se = d_rho[idx], s_rho[s], d_e[idx], s_e[s]
        ok = (dr > 0.0) & (sr > 0.0)
        # Index σ²_bio by the MEAN of the two node densities — monotone-compatible (the cross-node
        # dispersion scales with the average level) AND uses both endpoints, vs indexing by ρ_dst alone
        # (query is then circular at the crushed dest) or ρ_src alone (breaks monotonicity: low-source /
        # high-dest edges are high-residual at low μ).
        means.append(0.5 * (dr[ok] + sr[ok]))
        raws.append((dr[ok] - sr[ok]) ** 2)
        offs.append(dr[ok] / de[ok] + sr[ok] / se[ok])
    return means, raws, offs


def _fit_offset(means, raws, offs):
    cat = lambda parts: np.concatenate(parts) if parts else np.zeros(0)  # noqa: E731
    return MonotoneVarMean.fit_offset(cat(means), cat(raws), cat(offs))


_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_INTRON_BITS = BIT_INTRON_POS | BIT_INTRON_NEG


def _node_region_type(chain, region_arrays):
    """Per-chain-node region type for REGION nodes (0=intergenic, 1=intron, 2=exon; exon>intron), −1 on
    boundary nodes; plus the per-REGION type array. Vectorized bit test (≡ signature.coarse_type_from_signature)."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    R = sig.shape[0]
    rtype = np.where((sig & _EXON_BITS) != 0, 2, np.where((sig & _INTRON_BITS) != 0, 1, 0))  # per REGION
    ri_ = np.clip(idx, 0, R - 1)
    return np.where(kind == REGION, rtype[ri_], -1), rtype


def _gdna_seed_estimate(chain, statics, geometry, region_arrays, boundary_substrate, f_g_init, kappa):
    """The honest, NON-CIRCULAR population gDNA prior, fit ONCE on gDNA-clean seed nodes (§4.3).

    Returns ``(rho_global: float, sigma2_g: MonotoneVarMean)`` — the exposure-pooled gDNA rate + the gDNA
    between-node spread var~mean. Inputs are belief-independent (structural ``M/E`` + the strand-ONLY init
    ``f_g``), so this is computed once before the sweep and never refit (breaks the per-pass circularity).

    Seeds (per-node gDNA density + weight):
      * **structural** (always; the only path for UNSTRANDED data) — intergenic & intron regions and
        intergenic-exon & intron-exon boundary crossings (exon-facing side). Density = observed ``M/E``
        (gDNA-clean by structure under the **nascent-sparse** assumption). Weight 1.
      * **strand** (single-strand nodes not already structural — mainly exons; reach the capture-enriched
        range). Density = strand-deconv ``f_g_init·M/E``. Weight ``(2κ−1)²`` (the strand discriminability) so
        the seeds fade to 0 as κ→½ — stranded data gets the extra coverage, unstranded falls back to the
        structural seeds automatically (no hard threshold)."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    node_rtype, rtype = _node_region_type(chain, region_arrays)
    R = rtype.shape[0]
    EGl, EGr = geometry.eff_gdna_left, geometry.eff_gdna_right
    Ml, Mr = geometry.mass_left, geometry.mass_right

    # boundary flank types → clean (intron/intergenic)-exon boundary + which side is the exon.
    blr = np.asarray(boundary_substrate.left_region, dtype=np.int64)
    brr = np.asarray(boundary_substrate.right_region, dtype=np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(idx, 0, Bn - 1)
    lt = np.where((blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, R - 1)], -1)  # left flank type
    rt = np.where((brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, R - 1)], -1)  # right flank type
    left_clean = (lt == 0) | (lt == 1)
    right_clean = (rt == 0) | (rt == 1)
    clean_exon_bnd = is_bnd & ((left_clean & (rt == 2)) | (right_clean & (lt == 2)))
    exon_on_right = clean_exon_bnd & (rt == 2)

    # representative (mass, gDNA-eff): region = contained; clean-exon boundary = exon-facing side.
    mass = np.where(is_reg, Ml, np.where(exon_on_right, Mr, Ml))
    eff = np.maximum(np.where(is_reg, EGl, np.where(exon_on_right, EGr, EGl)), _EPS)
    rho_obs = mass / eff  # observed total density on the representative face (M/E)

    # Structural (M/E) seeds: intergenic regions (pure gDNA) + clean (intron/intergenic)-exon boundary
    # crossings (the gDNA-clean, capture-OBSERVABLE enriched-edge signal). Intron REGION interiors are NOT
    # structural — their contained mass carries the nascent RNA (so M/E is not gDNA-clean) and they are
    # DEPLETED under capture; a single-strand intron is instead strand-deconvolved below, and the gDNA-clean
    # intron signal the user intends is the intron-exon BOUNDARY crossing (already structural).
    struct_seed = (is_reg & (node_rtype == 0)) | clean_exon_bnd
    strand_seed = np.asarray(statics.strand_obs, bool) & ~struct_seed & (mass > 0.0)
    w_strand = float((2.0 * kappa - 1.0) ** 2)  # strand discriminability ∈[0,1]; →0 at κ→½
    dens = np.where(struct_seed, rho_obs, np.where(strand_seed, np.asarray(f_g_init, float) * rho_obs, 0.0))
    seed_w = np.where(struct_seed, 1.0, np.where(strand_seed, w_strand, 0.0))
    # Structural seeds are kept even at ZERO count (an intergenic region with 0 gDNA fragments over a large
    # exposure E is the STRONGEST evidence that gDNA is scarce — it drives ρ_global→0). Strand seeds need
    # counts to deconvolve, so they alone require mass>0.
    is_seed = struct_seed | strand_seed

    # exposure-pooled gDNA rate with ONE pseudocount of TOTAL gDNA (a=1; the global, not per-node — the
    # Poisson–Gamma posterior Gamma(a+G, E_tot)): G = Σ(w·gcount) (gcount = dens·E = mass for structural,
    # f_g·mass for strand), E_tot = Σ(w·E). ρ_global=(a+G)/E_tot (→1/E_tot≈0⁺ when G=0: gDNA scarce, never
    # zero); var_mean=(a+G)/E_tot² is the rate-estimate variance (1/CV²=a+G floors at the 1 pseudocount).
    sw = seed_w * is_seed
    G = float(np.sum(sw * dens * eff))
    E_tot = max(float(np.sum(sw * eff)), _EPS)
    rho_global = (1.0 + G) / E_tot
    var_mean = (1.0 + G) / (E_tot * E_tot)
    return rho_global, _fit_seed_varmean(chain, dens, eff, is_seed, seed_w), var_mean


def _fit_seed_varmean(chain, dens, eff, is_seed, seed_w):
    """σ²_g(μ) on adjacent SEED-SEED edges from the per-node seed gDNA density (the gDNA between-node spread).
    Edge weight = ``min`` of the two endpoints' seed weights (the weaker endpoint limits reliability)."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    is_seed = np.asarray(is_seed, bool)
    means, raws, offs, wts = [], [], [], []
    for nbr in (left, right):
        idx = np.where((nbr >= 0) & is_seed)[0]
        if not idx.size:
            continue
        s = nbr[idx]
        keep = is_seed[s]
        idx, s = idx[keep], s[keep]
        dr, sr, de, se = dens[idx], dens[s], eff[idx], eff[s]
        ok = (dr > 0.0) & (sr > 0.0)
        means.append(0.5 * (dr[ok] + sr[ok]))
        raws.append((dr[ok] - sr[ok]) ** 2)
        offs.append(dr[ok] / de[ok] + sr[ok] / se[ok])
        wts.append(np.minimum(seed_w[idx][ok], seed_w[s][ok]))
    cat = lambda p: np.concatenate(p) if p else np.zeros(0)  # noqa: E731
    return MonotoneVarMean.fit_offset(cat(means), cat(raws), cat(offs), cat(wts))


def _loglog_r2(z, y):
    """log-log fit R² of ``y`` on ``z`` (both > 0) — the z→ρ_g signal statistic for the significance gate."""
    z = np.asarray(z, float)
    y = np.asarray(y, float)
    ok = (z > 0.0) & (y > 0.0) & np.isfinite(z) & np.isfinite(y)
    if int(ok.sum()) < 4:
        return 0.0
    lx, ly = np.log(z[ok]), np.log(y[ok])
    X = np.vstack([np.ones_like(lx), lx]).T
    cf, *_ = np.linalg.lstsq(X, ly, rcond=None)
    sstot = float(np.sum((ly - ly.mean()) ** 2))
    return 1.0 - float(np.sum((ly - X @ cf) ** 2)) / max(sstot, _EPS)


def fit_enrichment_transfer(
    chain, statics, geometry, region_arrays, boundary_substrate, f_g_init, kappa, rho_global,
    *, n_perm: int = 200, sig_quantile: float = 0.95,
):
    """The per-REGION enrichment transfer ``ê(z) = E[ρ_g | z]`` + its conditional reliability ``σ²_bio(level)``
    (``enrichment_aware_calibration.md`` §J, Phase 2). Fit ONCE pre-loop on the strand-ONLY init ``f_g`` (the
    validated non-circular basis — teacher = single-strand exons, AMBIG withheld), so it carries no belief
    feedback.

    ``z(region)`` = the **region-facing** crossing density of the region's flanking **gDNA-clean**
    (intron/intergenic↔exon) boundaries, both edges averaged — a belief-independent enrichment predictor that
    never touches the node's own ``f_g`` (the decoupling that makes the de-enriched residual non-degenerate).
    ``ê`` is fit on the SINGLE-STRAND EXON regions (the only nodes whose enriched gDNA density is observed
    without a global); the exon edge→interior gradient is exon-specific, so seeds are NOT mixed in here.

    The fit RESPONSE ρ_g is a reliability-weighted BLEND of a STRAND-derived estimate (``f_g_init·T``,
    weight ``(2κ−1)²``) and a SPLICED-derived estimate (from the MOTIF-stranded spliced/pure-mature signal —
    strand-immune, weight ``1−(2κ−1)²``; Part III). κ→1 ⇒ strand; κ→½ ⇒ spliced (the unstranded path); a node
    uses the spliced term only where spliced is present at its flanking clean boundary.

    A **permutation significance gate** (real log-log R² vs the shuffled-``z`` null) collapses ``ê`` to the flat
    ``ρ_global`` when there is no ``z→ρ_g`` signal (off-capture / uniform gDNA) — capture is DETECTED, not
    flagged, and the layer reduces exactly to the genome-wide global. Returns
    ``(ehat: MonotoneMean, z: (n_nodes,) float [NaN where unusable], sigma2_level: MonotoneVarMean | None,
    significant: bool)``."""
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    node_rtype, rtype = _node_region_type(chain, region_arrays)
    R = rtype.shape[0]
    EGl, EGr = geometry.eff_gdna_left, geometry.eff_gdna_right
    Ml, Mr = geometry.mass_left, geometry.mass_right

    # gDNA-clean exon boundaries (intron/intergenic↔exon) — the same cleanliness test as _gdna_seed_estimate.
    blr = np.asarray(boundary_substrate.left_region, dtype=np.int64)
    brr = np.asarray(boundary_substrate.right_region, dtype=np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(idx, 0, Bn - 1)
    lt = np.where((blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, R - 1)], -1)
    rt = np.where((brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, R - 1)], -1)
    clean_exon_bnd = is_bnd & ((((lt == 0) | (lt == 1)) & (rt == 2)) | (((rt == 0) | (rt == 1)) & (lt == 2)))

    # z(region) = mean of the region-FACING crossing densities of its flanking clean boundaries: region i is
    # the RIGHT region of its left boundary and the LEFT region of its right boundary (so an exon reads its own
    # elevated edge, an intron/intergenic region its own depleted side — the (z, ρ_g) pair stays monotone).
    d_left = Ml / np.maximum(EGl, _EPS)   # crossing into the boundary's LEFT region
    d_right = Mr / np.maximum(EGr, _EPS)  # into its RIGHT region
    z = np.full(int(chain.n_nodes), np.nan)
    for i in np.where(is_reg)[0]:
        vals = []
        lb = int(left[i])
        if lb >= 0 and clean_exon_bnd[lb]:
            vals.append(d_right[lb])
        rb = int(right[i])
        if rb >= 0 and clean_exon_bnd[rb]:
            vals.append(d_left[rb])
        if vals:
            z[i] = float(np.mean(vals))

    flat = MonotoneMean.constant(rho_global)  # the off-capture / no-signal fallback (= genome-wide global)
    T = Ml / np.maximum(EGl, _EPS)
    strand_obs = np.asarray(statics.strand_obs, bool)
    ss_exon = is_reg & (node_rtype == 2) & strand_obs & (Ml > 0.0) & np.isfinite(z)
    fi = np.where(ss_exon)[0]
    if fi.size == 0:
        return flat, z, None, False

    # ---- the ê-fit RESPONSE ρ_g: a reliability-weighted BLEND of two estimators (Part III) ----
    # STRAND-derived (ρ_g=f_g_init·T) is gDNA-accurate when κ→1 but collapses to ~0.5·T as κ→½. SPLICED-derived
    # uses the MOTIF-stranded spliced (pure-mature) signal at the flanking clean boundary — strand-immune, so it
    # is the unstranded fallback: ρ_mature=M_spliced/E_rna_crossing(B); ρ_g=clip(M_unspliced−ρ_mature·E_rna_contained,0)/E_gdna.
    # Blend w_str=(2κ−1)² (strand discriminability, →0 at κ=½), w_spl=1−w_str; a node uses the spliced term only
    # where spliced is present at its flanking clean boundary (avail). κ→1 ⇒ strand; κ→½ ⇒ spliced; in between a
    # smooth hand-off. Replaces the old hard strand-contrast gate.
    fg = np.asarray(f_g_init, dtype=np.float64)
    fp = np.asarray(statics.free_pos, bool)  # single-strand: POS iff free_pos (else NEG) — the motif strand
    SPl, SPr = geometry.spliced_pos_left, geometry.spliced_pos_right
    SNl, SNr = geometry.spliced_neg_left, geometry.spliced_neg_right
    ERl = geometry.eff_rna_left  # contained RNA eff (the mature-projection target)
    ESPl, ESPr = geometry.eff_spl_left, geometry.eff_spl_right  # one-sided spliced HALF-TRIANGLE eff-len
    rho_strand = fg[fi] * T[fi]
    rho_spliced = np.zeros(fi.size)
    avail = np.zeros(fi.size, dtype=bool)
    for j, i in enumerate(fi):
        sl, sr = (SPl, SPr) if fp[i] else (SNl, SNr)  # spliced on this region's motif strand
        m_spl, e_spl = 0.0, _EPS
        lb = int(left[i])
        if lb >= 0 and clean_exon_bnd[lb] and sr[lb] > m_spl:  # R is lb's RIGHT region → faces lb's right side
            m_spl, e_spl = float(sr[lb]), float(ESPr[lb])
        rb = int(right[i])
        if rb >= 0 and clean_exon_bnd[rb] and sl[rb] > m_spl:  # R is rb's LEFT region → faces rb's left side
            m_spl, e_spl = float(sl[rb]), float(ESPl[rb])
        if m_spl > 0.0:
            # ρ_mature = spliced mass / its one-sided HALF-TRIANGLE eff-len (eff_spl), consistent with
            # node_densities + the message — was the FULL crossing eff_rna here (a ~2× under-projection of mature
            # ⇒ residual over-attributed to gDNA). Project that mature density onto the exon's contained RNA eff.
            rho_mature = m_spl / max(e_spl, _EPS)
            rho_spliced[j] = max(Ml[i] - rho_mature * ERl[i], 0.0) / max(EGl[i], _EPS)
            avail[j] = True
    w_str = float((2.0 * kappa - 1.0) ** 2)
    w_spl = 1.0 - w_str
    denom = w_str + w_spl * avail.astype(np.float64)
    keep = denom > _EPS  # drop unstranded-AND-no-spliced nodes (no reliable response → fall to the global)
    rho_g = (w_str * rho_strand + w_spl * avail * rho_spliced) / np.maximum(denom, _EPS)
    eff = np.asarray(EGl, dtype=np.float64)[fi]
    w_fit = denom * eff  # per-node total reliability × E (E stabilizes + nets; 0 → excluded)
    fi, zf, rho_g, eff, w_fit = fi[keep], z[fi][keep], rho_g[keep], eff[keep], w_fit[keep]
    if fi.size == 0:
        return flat, z, None, False

    ehat = MonotoneMean.fit(zf, rho_g, weight=w_fit, recal_weight=eff)
    # A DEGENERATE fit (too few points for the spline basis ⇒ MonotoneMean falls back to a flat constant,
    # marked lam=inf) has no enrichment structure AND makes the permutation test unreliable (M-9). Treat it as
    # not-significant ⇒ flat ρ_global. This ties the min-data guard to the spline's OWN basis floor (no new
    # magic number) — small synthetic scenarios stay byte-identical to the genome-wide global.
    if not np.isfinite(ehat.lam):
        return flat, z, None, False
    # significance: real log-log R² vs the permutation null (shuffle z among the fit nodes; seeded → deterministic)
    rng = np.random.default_rng(0)
    real_r2 = _loglog_r2(zf, rho_g)
    null = np.array([_loglog_r2(rng.permutation(zf), rho_g) for _ in range(int(n_perm))])
    if real_r2 <= float(np.quantile(null, sig_quantile)):
        return flat, z, None, False  # no enrichment signal → flat ρ_global (off-capture collapse)

    pred = ehat.predict(zf)
    raw = (rho_g - pred) ** 2                    # residual around the LEARNED mean (not the identity)
    offset = rho_g / np.maximum(eff, _EPS)       # the node's own density Poisson sampling floor
    sigma2_level = MonotoneVarMean.fit_offset(pred, raw, offset, w_fit)  # σ²_bio as a function of ê(z) level
    return ehat, z, sigma2_level, True


def _global_logprior(
    fgg, mass_global, eff_global, rho_global, sigma2_g, var_mean,
    *, ehat=None, z=None, sigma2_level=None, apply_mask=None,
):
    """Precompute the count-space global as a per-node BINOMIAL pseudo-count on f_g, ``(n_nodes, P)`` (§4):

        glob(f_g) = α·log f_g + (N − α)·log(1 − f_g),   α = N·μ,   μ = clip(ρ·E/M, 0, 1),   N = ρ² / σ²_prior

    The genome-wide baseline uses ``ρ = ρ_global`` and ``σ²_prior = var_mean + σ²_g(ρ_global)``. Where the
    ENRICHMENT-AWARE transfer is significant (``apply_mask`` = exon REGION nodes with a usable ``z``;
    :func:`fit_enrichment_transfer`), the per-node ``ρ`` is REPLACED (not added — §L/M-7) by the node's
    enrichment-predicted gDNA density ``ρ̂ = ê(z)``, with the conditional reliability
    ``σ²_prior = σ²_bio(ρ̂) + (ρ̂ + 1/E)/E`` (the level-conditional between-node spread + the node's own density
    Poisson floor). Off-capture / not-significant ⇒ ``apply_mask`` is empty ⇒ BYTE-IDENTICAL to the genome-wide
    baseline (the per-node ``N`` array is then a constant ``n_glob``).

    ``N`` is the **M-INDEPENDENT** population confidence (``= 1/CV²`` of the rate) — a hyperprior is naturally
    imprecise, so it can NEVER overrule a node's own (strand) evidence; the MEAN ``μ`` keeps the honest
    density→fraction ``E/M`` conversion. ``σ²_prior`` = the rate-estimate variance ``var_mean = (1+G)/E_tot²``
    (which carries the 1-pseudocount floor, so ``N → 1`` at uniform zero-gDNA — a gentle zero-baseline, never
    ``0/0``) PLUS the between-node spread ``σ²_g`` (large under capture ⇒ ``N → ρ²/σ²_g ≈`` small, imprecise;
    small under uniform present gDNA ⇒ ``N`` large, confident). Two-sided (mode μ); applied to ALL solvable
    nodes (the strand_obs fork is dissolved)."""
    s2_between = max(float(sigma2_g.predict(np.array([max(rho_global, _EPS)]))[0]), 0.0)
    n_glob = rho_global * rho_global / max(var_mean + s2_between, _EPS)  # scalar, M-independent
    mu = np.clip(rho_global * eff_global / np.maximum(mass_global, _EPS), 0.0, 1.0)
    n_node = np.full(mu.shape, n_glob, dtype=np.float64)  # per-node count (= n_glob everywhere by default)
    if ehat is not None and apply_mask is not None and sigma2_level is not None and bool(np.any(apply_mask)):
        m = np.asarray(apply_mask, bool)
        rho_hat = np.maximum(ehat.predict(np.asarray(z)[m]), 0.0)  # enrichment-predicted gDNA density
        mu[m] = np.clip(rho_hat * eff_global[m] / np.maximum(mass_global[m], _EPS), 0.0, 1.0)
        s2_bio = np.maximum(sigma2_level.predict(rho_hat), 0.0)    # level-conditional between-node spread
        eff_m = np.maximum(eff_global[m], _EPS)
        pois = (rho_hat + _MSG_PSEUDOCOUNT / eff_m) / eff_m        # node density Poisson floor (matches _message)
        n_node[m] = rho_hat * rho_hat / np.maximum(s2_bio + pois, _EPS)
    return _binom_pseudo(fgg[None, :], mu[:, None], n_node[:, None])


def fit_rna_varmean(
    chain: NodeChain, densities: NodeDensities, geometry: NodeGeometry, statics: NodeStatics
) -> MonotoneVarMean:
    """Fit the RNA message reliability ``σ²_bio(μ)`` on the frozen snapshot, POOLING both strands (the symmetric
    RNA process). A strand-s edge is ``live`` only where strand s is continuous on BOTH endpoints (``free_s`` —
    the transcript-structure gate), so the curve sees the genuine same-strand cross-node RNA dispersion
    (INCLUDING the AMBIG nodes' per-strand densities). Refit per pass: unlike gDNA (seed-based, phantom risk →
    a non-circular firewall), RNA is well-resolved by strand+spliced so the swept densities are accurate and the
    fit adapts; a non-circular init-based fit measured WORSE (the AMBIG RNA dispersion only exists after the
    sweep — the gDNA/RNA observability asymmetry, Phase 3a). The per-strand density is spliced-inclusive."""
    mp, rp, op = _edge_varmean(
        chain,
        densities.rho_pos_left,
        densities.rho_pos_right,
        geometry.eff_rna_left,
        geometry.eff_rna_right,
        statics.free_pos,
    )
    mn, rn, on = _edge_varmean(
        chain,
        densities.rho_neg_left,
        densities.rho_neg_right,
        geometry.eff_rna_left,
        geometry.eff_rna_right,
        statics.free_neg,
    )
    return _fit_offset(mp + mn, rp + rn, op + on)


def fit_gdna_varmean(
    chain: NodeChain, densities: NodeDensities, geometry: NodeGeometry, statics: NodeStatics
) -> MonotoneVarMean:
    """Fit the gDNA message-reliability ``σ²_g(μ)`` on the frozen snapshot over **single-strand** edges
    (``statics.strand_obs = free_pos ^ free_neg`` on both endpoints). The twin of :func:`fit_rna_varmean`:
    refit once per OUTER iteration on the converged belief (the nested-loop design), NOT entangled in the inner
    sweep — so neither component is privileged (the firewall asymmetry is gone).

    Why single-strand, not all nodes: fitting on ALL nodes is circular — in zero-gDNA + nascent the AMBIG/exon
    gDNA is the deconv's own (nascent-confounded) guess, so refitting learns the phantom and amplifies it
    (`forward_backward_plan.md` §5). Single-strand nodes are definitively strand-solvable, so their gDNA is
    trustworthy; crucially this INCLUDES single-strand exons — in real data most genes are unexpressed, so many
    exons carry only gDNA, the best signal for the hybrid-capture gDNA-on-exon enrichment (which the intron-exon
    boundary crossings underestimate). AMBIG nodes (not strand-separable) are excluded. PROVISIONAL — revisit
    after the spliced-imputation ~2× fix and FB."""
    live = np.asarray(statics.strand_obs, dtype=bool)  # single-strand: definitively strand-solvable gDNA
    m, r, o = _edge_varmean(
        chain,
        densities.rho_g_left,
        densities.rho_g_right,
        geometry.eff_gdna_left,
        geometry.eff_gdna_right,
        live,
    )
    return _fit_offset(m, r, o)


def node_sweep(
    chain: NodeChain,
    statics: NodeStatics,
    geometry: NodeGeometry,
    belief: NodeBelief,
    region_arrays,
    boundary_substrate,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
    max_passes: int = 1,
    convergence_delta: float = 1e-3,
    max_outer: int = 1,
    outer_convergence_delta: float = 1e-3,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a **FORWARD-BACKWARD** inner solve inside a var~mean OUTER fixed point.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). Entangling the var~mean refit into the sweep made a moving precision
    target (a measured limit cycle); the nested design separates them — the INNER FB solve runs at FIXED
    var~mean, the OUTER loop refits BOTH reliability curves on the solved belief (symmetric — neither gDNA nor
    RNA privileged).

    BEFORE the loop (outer-iteration-0 bootstrap): the non-circular population gDNA prior on gDNA-clean seeds
    (:func:`_gdna_seed_estimate`) — exposure-pooled ``ρ_global`` (anchored, NOT refit) + the seed gDNA spread
    ``σ²_g`` — the enrichment transfer ``ê(z)`` (anchored), and the bootstrap RNA ``σ²_bio`` on the init belief.

    Each OUTER iteration runs ONE FB pass: (A) a batched message-free LOCAL solve → per-component (fraction,
    precision); (B) a FORWARD scan L→R accumulating the left-context message α and the forward belief (local ⊗
    incoming message — NOT the reverse, true tree BP — so a thin node passes the upstream α on: the relay);
    (C) a BACKWARD scan R→L → β; (D) combine α⊗β and one batched FINAL lattice solve. Then refit ``rna_vm`` +
    ``gdna_vm`` on the solved belief and recompute the global NB log-prior. Stop once the belief stops moving
    between outer iterations (``outer_convergence_delta``, capped at ``max_outer``).

    The per-node ψ solve integrates: the strand likelihood, the Jeffreys reference (single-strand nodes), the
    count-space global NB prior (ALL solvable nodes — the fork is dissolved), and the gDNA + per-strand RNA
    imputation messages (the spliced/MATURE floor rides the RNA message via the half-triangle ``ESP`` eff-len).
    The emission gate is a three-term Boolean over the components (gDNA / +RNA / −RNA): each RNA strand flows
    only where that strand is continuous on both endpoints (``free_s``); gDNA is genomically continuous and
    strand-agnostic, so it flows wherever there is unspliced facing mass — including across TSS/TES and other G1
    seams (a locked all-gDNA node is a confident emitter). Only G2/G3 nodes with data are written; G1 sinks +
    empty nodes keep their init; a node with no neighbour message is solved on its own strand+global evidence
    (the intrinsic solve folds into the final batched solve at prec=0). ``max_passes``/``convergence_delta`` are retained for API
    compatibility but unused — FB is single-pass (exact on the chain given the local beliefs). Returns
    ``(NodeBelief, outer_deltas)`` (the per-outer-iteration max-Δf_g — the convergence diagnostic)."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    fp, fn = statics.free_pos, statics.free_neg
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    # precision state (Phase 1: computed + carried; consumed by the honest message send in Phase 2).
    var_pos = np.asarray(belief.var_pos, dtype=np.float64).copy()
    var_neg = np.asarray(belief.var_neg, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()

    # per-face geometry as (left, right) pairs, indexed by face (0=left, 1=right).
    EG = (geometry.eff_gdna_left, geometry.eff_gdna_right)
    ER = (geometry.eff_rna_left, geometry.eff_rna_right)
    ESP = (geometry.eff_spl_left, geometry.eff_spl_right)  # one-sided spliced-crossing half-triangle eff-len
    MS = (geometry.mass_left, geometry.mass_right)
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # per-node "global" gDNA support: region = its contained support; boundary = the AVERAGED per-side density
    # length ½(E_l+E_r) over the total crossing mass.
    eff_global = np.where(is_reg, EG[0], 0.5 * (EG[0] + EG[1]))
    mass_global = np.where(is_reg, MS[0], MS[0] + MS[1])

    lattice = _simplex_lattice(int(n_grid))
    fpg, fng, fgg = lattice
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    # Two distinct gates, both from the region SIGNATURE (never the counts — count-zero-info):
    #   * SOLVE gate (`solvable`): a node deconvolves its own gDNA-vs-RNA split iff it admits ≥1 RNA strand and
    #     has unspliced mass. A G1 node — no admissible RNA strand: an intergenic region, or a gene-boundary seam
    #     (TSS/TES, opposite-strand exon↔exon) — is a LOCKED all-gDNA node ({0,0,1}); it is not solved, keeping
    #     its init (RNA cannot cross a gene boundary, so its unspliced mass is purely gDNA).
    #   * EMISSION gate (per component, in `_scan`): which MESSAGES a node sends. A three-term Boolean over the
    #     components gDNA / +RNA / −RNA, structural and symmetric — defined at the top of the scan loop.
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)

    # BOOTSTRAP (outer-iteration-0 var~mean). The population gDNA prior on gDNA-clean seeds (§4.3): the
    # exposure-pooled rate ρ_global (the ANCHORED global baseline — never refit) + var_mean (its variance) +
    # the seed gDNA spread σ²_g (the bootstrap ``gdna_vm`` message reliability — REFIT each outer iteration on
    # the converged belief by fit_gdna_varmean, symmetric with RNA). Strand seeds use the strand-ONLY init f_g.
    rho_global, gdna_vm, var_mean = _gdna_seed_estimate(
        chain, statics, geometry, region_arrays, boundary_substrate, f_g, kappa
    )
    # The ENRICHMENT-AWARE transfer ê(z), fit once on the strand-only init f_g (anchored, like ρ_global). Under
    # capture it learns the edge→interior gDNA-density transfer; off-capture the significance gate fails and ê
    # collapses to the flat ρ_global. It REPLACES the genome-wide global mean for the exon REGION nodes with a
    # usable z (the AMBIG-exon targets); everything else keeps ρ_global (§J/P3).
    ehat, z_enrich, sigma2_level, enrich_sig = fit_enrichment_transfer(
        chain, statics, geometry, region_arrays, boundary_substrate, f_g, kappa, rho_global
    )
    node_rtype, _rtype = _node_region_type(chain, region_arrays)
    enrich_apply = enrich_sig & is_reg & (node_rtype == 2) & np.isfinite(z_enrich)
    # Bootstrap RNA reliability on the INIT belief (crude — the AMBIG RNA is still parked in gDNA, so this
    # curve is refined by the outer loop once the sweep resolves it). Symmetric bootstrap with the gDNA seed fit.
    rna_vm = fit_rna_varmean(
        chain, node_densities(NodeBelief(f_pos, f_neg, f_g, var_pos, var_neg, var_g), geometry),
        geometry, statics,
    )
    # The count-space global log-prior on f_g (M-INDEPENDENT strength so it can never overrule a node's own
    # strand evidence; ALL solvable nodes; enrichment-aware on the exon override set). Recomputed each outer
    # iteration with the refit gdna_vm (ρ_global / ê stay anchored), so define a helper.
    def _global_lp(gvm):
        return _global_logprior(
            fgg, mass_global, eff_global, rho_global, gvm, var_mean,
            ehat=ehat, z=z_enrich, sigma2_level=sigma2_level, apply_mask=enrich_apply,
        )

    global_lp = _global_lp(gdna_vm)

    # Genomic node order for the forward/backward scans (within each ref path; left/right break at −1). The
    # scans are sequential per node, so iterate as a Python list of ints (faster than numpy scalar indexing).
    order_list = [int(x) for x in np.asarray(chain.order)]
    n_nodes = f_g.shape[0]
    PSEUDO = _MSG_PSEUDOCOUNT

    outer_deltas = []
    prev_outer = None
    for _outer in range(int(max_outer)):
        # FORWARD-BACKWARD inner solve (one pass — on the chain it is exact given the local beliefs + var~mean,
        # so no inner iteration: the outer loop iterates the var~mean fixed point). (A) message-free LOCAL
        # beliefs (one batched lattice solve) → per-component (fraction, precision); (B) FORWARD scan L→R
        # accumulating the left-context message α + the forward belief; (C) BACKWARD scan R→L → β; (D) combine
        # α⊗β and one batched FINAL solve. Each node sees the WHOLE-chain context (relay) in one pass — vs
        # Jacobi's one-hop — because the forward belief at a thin node is dominated by, and passes on, the
        # upstream α. ``true`` tree BP: the forward belief excludes β (no double-counting), unlike Jacobi's
        # full-marginal source.

        # (A) LOCAL message-free beliefs.
        psi_loc = _local_loglik(
            statics.u_pos, statics.u_neg, statics.spliced_pos, statics.spliced_neg,
            fp, fn, kappa, od_g, od_r, lattice,
            strand_obs=statics.strand_obs, global_logprior=global_lp,
        )
        fg_loc, fp_loc, fn_loc, vg_loc, vp_loc, vn_loc = _node_marginals(psi_loc, fpg, fng, fgg)
        pg_loc = 1.0 / np.maximum(vg_loc, _EPS)  # local precision (var floored: a sharp belief ⇒ large finite)
        pp_loc = 1.0 / np.maximum(vp_loc, _EPS)
        pn_loc = 1.0 / np.maximum(vn_loc, _EPS)

        # per-edge σ²_bio — the ONLY varmean.predict, batched on the LOCAL densities (the curve query point need
        # not track the evolving forward belief). Density per node per face incl. the RNA one-sided spliced floor.
        def _face_dens(fc, eff, mass, spl=None, esp=None):
            d = fc * mass / np.maximum(eff, _EPS)
            return d if spl is None else d + spl / np.maximum(esp, _EPS)

        dG = (_face_dens(fg_loc, EG[0], MS[0]), _face_dens(fg_loc, EG[1], MS[1]))
        dP = (_face_dens(fp_loc, ER[0], MS[0], SP[0], ESP[0]),
              _face_dens(fp_loc, ER[1], MS[1], SP[1], ESP[1]))
        dN = (_face_dens(fn_loc, ER[0], MS[0], SN[0], ESP[0]),
              _face_dens(fn_loc, ER[1], MS[1], SN[1], ESP[1]))

        def _edge_sigma2(nbr, sf, df, dens, vm):
            valid = nbr >= 0
            src = np.where(valid, nbr, 0)
            return np.asarray(vm.predict(np.maximum(0.5 * (dens[sf][src] + dens[df]), _EPS)), float)

        # forward edge: src=left, src face=1 (right) → dst face=0 (left). backward: src=right, face 0 → dst 1.
        s2gf, s2gb = _edge_sigma2(left, 1, 0, dG, gdna_vm), _edge_sigma2(right, 0, 1, dG, gdna_vm)
        s2pf, s2pb = _edge_sigma2(left, 1, 0, dP, rna_vm), _edge_sigma2(right, 0, 1, dP, rna_vm)
        s2nf, s2nb = _edge_sigma2(left, 1, 0, dN, rna_vm), _edge_sigma2(right, 0, 1, dN, rna_vm)

        def _scan(seq, nbr, sf, df, s2g, s2p, s2n):
            """Sequential scan: project the running belief from each node's ``nbr`` (src face ``sf`` → dst face
            ``df``) into the dst, then fold it into the dst's running belief (local ⊗ incoming message; NOT the
            reverse message — true tree BP). Returns the per-node message (mode, prec) per component. O(m)."""
            fbg, fbp, fbn = fg_loc.copy(), fp_loc.copy(), fn_loc.copy()  # running belief (starts at local)
            vbg, vbp, vbn = vg_loc.copy(), vp_loc.copy(), vn_loc.copy()
            amg, apg = np.zeros(n_nodes), np.zeros(n_nodes)  # gDNA message (mode, prec)
            amp, app = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-pos
            amn, apn = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-neg
            EGs, EGd, ERs, ERd, ESPs = EG[sf], EG[df], ER[sf], ER[df], ESP[sf]
            MSs, MSd, SPs, SPd, SNs, SNd = MS[sf], MS[df], SP[sf], SP[df], SN[sf], SN[df]
            for i in seq:
                lsrc = nbr[i]
                if lsrc < 0:
                    continue
                md = MSd[i] if MSd[i] > _EPS else _EPS
                egd = EGd[i] if EGd[i] > _EPS else _EPS
                erd = ERd[i] if ERd[i] > _EPS else _EPS
                sm = MSs[lsrc]  # source facing unspliced mass
                # STRUCTURAL emission gate — one Boolean per component; src AND dst must admit it. gDNA is
                # genomically universal (admitted everywhere) ⇒ it gates on facing mass alone; each RNA strand
                # transmits only where THAT strand is continuous across the edge (free_pos / free_neg on both).
                emit_g = sm > _EPS
                emit_p = fp[lsrc] and fp[i] and (sm > _EPS or SPs[lsrc] > _EPS)
                emit_n = fn[lsrc] and fn[i] and (sm > _EPS or SNs[lsrc] > _EPS)
                # gDNA — a G1 seam (intergenic / TSS / TES) is a locked, confident all-gDNA emitter.
                if emit_g:
                    eg = EGs[lsrc] if EGs[lsrc] > _EPS else _EPS
                    rho = fbg[lsrc] * sm / eg
                    vo = (sm / eg) ** 2 * vbg[lsrc]
                    pois = (rho + PSEUDO / eg) / eg
                    pr = (md / egd) ** 2 / max(vo + s2g[i] + pois, _EPS)
                    # mode is a FRACTION ∈ [0,1]: the source density expressed as the dst fraction. Unclipped it
                    # can run to 1e4+ when the source is denser than the dst can hold (cliff / low-count dst),
                    # and `−½·prec·(f−mode)²` over the lattice then amplifies the pull by `mode`× (a vertex slam).
                    mo = min(max(rho * egd / md, 0.0), 1.0)
                    amg[i], apg[i] = mo, pr
                    pt = pg_loc[i] + pr
                    fbg[i], vbg[i] = (pg_loc[i] * fg_loc[i] + pr * mo) / pt, 1.0 / pt
                # RNA-pos — nascent + one-sided spliced/MATURE floor (strand-continuous edge only).
                if emit_p:
                    er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                    esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                    rho = fbp[lsrc] * sm / er + SPs[lsrc] / esp
                    vo = (sm / er) ** 2 * vbp[lsrc] if sm > _EPS else 0.0
                    pois = (rho + PSEUDO / er) / er
                    pr = (md / erd) ** 2 / max(vo + s2p[i] + pois, _EPS)
                    mo = min(max((rho * erd - SPd[i]) / md, 0.0), 1.0)  # fraction ∈[0,1] (also kills spliced over-sub <0)
                    amp[i], app[i] = mo, pr
                    pt = pp_loc[i] + pr
                    fbp[i], vbp[i] = (pp_loc[i] * fp_loc[i] + pr * mo) / pt, 1.0 / pt
                # RNA-neg — symmetric.
                if emit_n:
                    er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                    esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                    rho = fbn[lsrc] * sm / er + SNs[lsrc] / esp
                    vo = (sm / er) ** 2 * vbn[lsrc] if sm > _EPS else 0.0
                    pois = (rho + PSEUDO / er) / er
                    pr = (md / erd) ** 2 / max(vo + s2n[i] + pois, _EPS)
                    mo = min(max((rho * erd - SNd[i]) / md, 0.0), 1.0)  # fraction ∈[0,1]
                    amn[i], apn[i] = mo, pr
                    pt = pn_loc[i] + pr
                    fbn[i], vbn[i] = (pn_loc[i] * fn_loc[i] + pr * mo) / pt, 1.0 / pt
            return amg, apg, amp, app, amn, apn

        a = _scan(order_list, left, 1, 0, s2gf, s2pf, s2nf)        # forward (α: left context)
        b = _scan(order_list[::-1], right, 0, 1, s2gb, s2pb, s2nb)  # backward (β: right context)

        # (D) combine α⊗β (precision-weighted product) per component → one batched FINAL solve.
        def _comb(am_a, ap_a, am_b, ap_b):
            pc = ap_a + ap_b
            return np.where(pc > _EPS, (ap_a * am_a + ap_b * am_b) / np.maximum(pc, _EPS), 0.0), pc

        mode_g, prec_g = _comb(a[0], a[1], b[0], b[1])
        mode_p, prec_p = _comb(a[2], a[3], b[2], b[3])
        mode_n, prec_n = _comb(a[4], a[5], b[4], b[5])
        psi = _local_loglik(
            statics.u_pos, statics.u_neg, statics.spliced_pos, statics.spliced_neg,
            fp, fn, kappa, od_g, od_r, lattice,
            strand_obs=statics.strand_obs, global_logprior=global_lp,
            gdna_imp_mode=mode_g, gdna_imp_prec=prec_g,
            rna_imp_mode=(mode_p, mode_n), rna_imp_prec=(prec_p, prec_n),
        )
        # write back only SOLVABLE nodes (G1 sinks / empty keep their signature-binary init).
        mg_, mp_, mn_, vg_, vp_, vn_ = _node_marginals(psi, fpg, fng, fgg)
        f_g = np.where(solvable, np.clip(mg_, 0.0, 1.0), f_g)
        f_pos = np.where(solvable, np.clip(mp_, 0.0, 1.0), f_pos)
        f_neg = np.where(solvable, np.clip(mn_, 0.0, 1.0), f_neg)
        var_g = np.where(solvable, vg_, var_g)
        var_pos = np.where(solvable, vp_, var_pos)
        var_neg = np.where(solvable, vn_, var_neg)

        if _capture is not None:  # inert diagnostic hook (overwritten each outer iter → ends on the last)
            # strand-ONLY local belief (no global prior, no messages) — to split the LOCAL error into the
            # strand likelihood vs the global gDNA prior contribution.
            fg_strand = _node_marginals(
                _local_loglik(
                    statics.u_pos, statics.u_neg, statics.spliced_pos, statics.spliced_neg,
                    fp, fn, kappa, od_g, od_r, lattice, strand_obs=statics.strand_obs, global_logprior=None,
                ), fpg, fng, fgg,
            )[0]
            _capture.update(
                fg_loc=fg_loc, fg_strand=fg_strand, fp_loc=fp_loc, fn_loc=fn_loc,
                vg_loc=vg_loc, vp_loc=vp_loc, vn_loc=vn_loc,
                a_fwd=a, b_bwd=b, mode_g=mode_g, prec_g=prec_g, mode_p=mode_p, prec_p=prec_p,
                mode_n=mode_n, prec_n=prec_n, f_g=f_g.copy(), f_pos=f_pos.copy(), f_neg=f_neg.copy(),
                var_g=var_g.copy(), solvable=solvable, rho_global=rho_global, gdna_vm=gdna_vm, rna_vm=rna_vm,
                # boundary-emission geometry: gDNA emits iff facing unspliced mass>0 (strand-agnostic);
                # RNA iff free_s on both endpoints & (unspliced or spliced facing mass). Capture the faces.
                mass_l=MS[0], mass_r=MS[1], spl_l=SP[0] + SN[0], spl_r=SP[1] + SN[1],
                free_pos=np.asarray(fp, bool), free_neg=np.asarray(fn, bool),
                # enrichment transfer ê(z): the fitted object + the predictor z + the apply (AMBIG-exon) mask +
                # the global geometry (μ = clip(ρ·eff_global/mass_global, 0, 1) is the implied prior fraction).
                ehat=ehat, z_enrich=z_enrich, enrich_apply=enrich_apply, enrich_sig=enrich_sig,
                eff_global=eff_global, mass_global=mass_global,
            )

        # OUTER LOOP — refit BOTH var~mean reliability curves on the inner-converged belief (SYMMETRIC: neither
        # gDNA nor RNA privileged) + recompute the gDNA global log-prior with the new gdna_vm (ρ_global / ê stay
        # anchored). This is the var~mean fixed point, separated from the inner sweep so the precision target
        # stops moving.
        snap = node_densities(NodeBelief(f_pos, f_neg, f_g, var_pos, var_neg, var_g), geometry)
        rna_vm = fit_rna_varmean(chain, snap, geometry, statics)
        gdna_vm = fit_gdna_varmean(chain, snap, geometry, statics)
        global_lp = _global_lp(gdna_vm)

        # OUTER convergence: did the inner-converged belief stop moving between outer iterations? (The refit
        # changes the var~mean; if that no longer moves the belief, the fixed point is reached.)
        cur_outer = (f_g.copy(), f_pos.copy(), f_neg.copy())
        if prev_outer is not None:
            od = max(
                (float(np.max(np.abs(c - p))) if c.size else 0.0)
                for c, p in zip(cur_outer, prev_outer)
            )
            outer_deltas.append(od)
            if od < outer_convergence_delta:
                break
        prev_outer = cur_outer

    return (
        NodeBelief(
            f_pos=f_pos, f_neg=f_neg, f_g=f_g,
            var_pos=var_pos, var_neg=var_neg, var_gdna=var_g,
        ),
        outer_deltas,
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
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return NodeDeconv(
        gdna_mass=f_g * mass_u,
        rna_mass=(1.0 - f_g) * mass_u + mass_s,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
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
            gdna_mass=side_fg * m_u,
            rna_mass=(1.0 - side_fg) * m_u + m_s,
            gdna_frac=side_fg,
            rna_pos_frac=np.zeros(R),
            rna_neg_frac=np.zeros(R),
        )

    return _side(left_fg, substrate.left), _side(right_fg, substrate.right)
