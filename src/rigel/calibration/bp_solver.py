"""The belief-propagation calibration solver: per-node gDNA-vs-RNA deconvolution by a bidirectional sweep.

Deconvolves each node's unspliced fragment mass into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the unified region↔boundary chain (`node_chain`), by a single forward-backward
(L→R then R→L) belief-propagation pass (exact on the chain, a forest of linear paths). Each per-node solve
(`simplex_logodds`, the log-density log-odds solver) reconciles three sources of information: the intrinsic
strand likelihood (the Beta-Binomial tilt — the only signal a count carries), the cross-node imputation
messages, and the population gDNA prior. Theory: the count-zero-information principle in
`docs/calibration/CALIBRATION_ARCHITECTURE.md`; the message-precision model (per-edge disagreement-aware
σ²_edge) in `docs/calibration/dispersion_aware_message_precision.md`.

Module layout:
* `build_node_geometry` — the per-node, per-FACE static geometry (unspliced mass, per-component gDNA-/RNA-FL
  effective lengths, the boundary's one-sided motif-stranded spliced mass). Per-face because a boundary's
  two sides lie in different-sized regions and a message uses the side FACING the destination; a region
  presents the same contained geometry both ways. Eff-len primitives are the single-owner
  `effective_length.py` (`region_eff_length` = E[max(0,L−ℓ)] contained; `boundary_side_eff_length` =
  E[min(ℓ,L)] per-side crossing). gDNA uses the gDNA FL; RNA (nascent unspliced + spliced) the RNA FL.
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.
* `node_densities` / `build_node_statics` — the per-node density snapshot + the static per-node inputs.
* `_gdna_seed_estimate` — the ANCHORED global gDNA prior: the population baseline
  ρ_global + seed between-node spread σ²_g (the seed-based non-circular firewall) + the enrichment transfer
  ê(z). All fit ONCE before the solve (NOT a per-message reliability — the message precision is per-edge).
* `node_sweep` — the single forward-backward sweep. Message precision is the per-edge **disagreement-aware**
  σ²_edge built inside `_scan` (`dispersion_aware_message_precision.md`): a message is trusted in inverse
  proportion to its surprise vs the destination's message-free local belief — NOT a fitted σ²_bio(μ) curve,
  so there is no var~mean fixed point and no outer loop.
* `chain_region_deconv` / `chain_boundary_side_deconv` — project the converged belief back to the per-region
  / per-boundary-side masses the `CalibrationResult` consumes.
"""

from __future__ import annotations

import math
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
from .simplex_logodds import (
    _logodds_grid,
    _solve_nodes_logodds_all,
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
    "node_sweep",
    "chain_region_deconv",
    "chain_boundary_side_deconv",
]

_EPS = 1.0e-9
_MSG_PSEUDOCOUNT = 1.0  # one pseudo-observation: a density from a finite count can never have zero
#                         sampling variance, so a message precision can never escape to ∞ (var stabilizer).
# ── The pass-1 global gDNA prior is STABILITY-ONLY: it gives low/zero-gDNA nodes a finite baseline so the
#    solve stays sane, capped at one pseudo-observation (`_GLOBAL_STAB_PREC`) so it can NEVER drive a
#    solution. A node resolves from its STRAND likelihood + the messages (strong strand pins the tilt;
#    unstranded relies on the messages) and, in PASS 2, the trained Phase-2 gDNA-density KDE — the only
#    trained prior (it supersedes the earlier enrichment-transfer ê(z), now removed). ──
_GLOBAL_STAB_PREC = 1.0  # one pseudo-observation — the global can never override a node's own data.
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
        # No per-node spliced FLOOR. Spliced (mature) handling is OWNED by the message system (the B→exon
        # MEASUREMENT message + the exon→B absorption); a second floor here would double-count it AND inflate a
        # boundary's UNSPLICED f_pos with mature → phantom nascent into introns (matrix-confirmed: removing it
        # is ≥ keeping it in every κ × capture × ±gDNA regime). Zeroed for ALL nodes (regions carry ~0 contained
        # spliced anyway). The solver's spliced-lower-bound term thus operates on zeros (inert; dead-code
        # cleanup is a follow-up).
        spliced_pos=np.zeros_like(pick(r_sp, b_sp)),
        spliced_neg=np.zeros_like(pick(r_sn, b_sn)),
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
    logodds_window: float = 10.0,
    statics: "NodeStatics | None" = None,
) -> NodeBelief:
    """The signature-binary G1/G2/G3 initial :class:`NodeBelief` on the unified chain.

    All chain nodes are strand-solved by the log-density 1-D/2-D log-odds solver (:mod:`simplex_logodds`,
    ``O(m·K)``; the bare strand likelihood + the Jeffreys reference at single-strand nodes; NO global prior,
    NO imputation — those enter the sweep, P2/P3). The signature-binary class overrides (:func:`_type_belief`)
    then set the G1/G2/G3 belief. Single-strand introns resolve to ``f_g≈0`` from the BB tilt alone (the
    zero-gDNA gate); intergenic / TSS sinks lock at ``{0,0,1}``; AMBIG nodes hold ``{0,0,1}`` at MAX variance
    for the sweep. Pass ``statics`` to reuse a prebuilt :class:`NodeStatics` (the sweep does)."""
    st = (
        statics
        if statics is not None
        else build_node_statics(chain, substrate, boundary_substrate, region_arrays)
    )
    deconv = _solve_nodes_logodds_all(
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
        L=logodds_window,
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
    # var_mean = the global rate-estimate variance, log-density form (D4, the delta method = this design's own
    # 1/count principle): Var(log ρ̂) ≈ Var(ρ̂)/ρ̂² = [(1+G)/E_tot²]/[(1+G)/E_tot]² = 1/(1+G) — the inverse
    # pooled effective count (NOT the opaque trigamma; transparent + consistent with pois_log=1/(count+1)).
    # Gentle N→~0.5 at zero-gDNA (G=0).
    var_mean = 1.0 / (1.0 + G)
    return rho_global, _fit_seed_varmean(chain, dens, eff, is_seed, seed_w), var_mean


def _fit_seed_varmean(chain, dens, eff, is_seed, seed_w):
    """σ²_g(μ) on adjacent SEED-SEED edges from the per-node seed gDNA density (the gDNA between-node spread),
    in the log-density form (twin of :func:`_edge_varmean`). Edge weight = ``min`` of the two endpoints' seed
    weights (the weaker endpoint limits reliability)."""
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
        dr, sr, de, se = dr[ok], sr[ok], de[ok], se[ok]
        means.append(0.5 * (dr + sr))
        raws.append((np.log(dr) - np.log(sr)) ** 2)
        offs.append(1.0 / (dr * de + 1.0) + 1.0 / (sr * se + 1.0))
        wts.append(np.minimum(seed_w[idx][ok], seed_w[s][ok]))
    cat = lambda p: np.concatenate(p) if p else np.zeros(0)  # noqa: E731
    return MonotoneVarMean.fit_offset(cat(means), cat(raws), cat(offs), cat(wts))


def _floor_estimate(chain, geometry, region_arrays, f_g_init, kappa):
    """The DEPLETED gDNA FLOOR — the population background gDNA density from intergenic + intron REGIONS
    (NOT boundaries, NOT exons; the user's directive). Globally these regions are pure gDNA and DEPLETED
    (off-target under capture); introns may carry sparse nascent — handled below. The CONFIDENT floor (a
    coherent depleted population ⇒ tight spread) is applied to the depleted REGION nodes by
    :func:`_global_logprior`: it pins a floor-level intron to ``f_g≈1`` and gives an intron with density
    EXCESS over the floor its nascent (``target = ρ_floor·E/M``) — the nascent-from-self-evidence
    principle, no gate. Enriched exons are excluded (they need ``ê(z)``).

    The per-region gDNA density is **strand-deconvolved where the strand is informative**: a region's gDNA
    fraction is ``(1−w)·1 + w·f_g_init`` with ``w = (2κ−1)²`` (the strand discriminability) — so a
    stranded RNA-rich intron is discounted out of the floor (``f_g_init→0``), while unstranded data falls
    back to the raw ``M/E_g`` (``w→0``: the nascent-sparse assumption — introns ARE the gDNA floor).
    Intergenic regions are locked ``f_g_init=1`` ⇒ their full ``M/E_g`` always counts.

    Returns ``(rho_floor, s2_floor, var_mean_floor, floor_mask)``: the exposure-pooled floor rate, the
    between-region log-density SPREAD (the floor tightness — biological/CNV excess over the per-region
    Poisson floor), the rate-estimate variance ``1/(1+G)``, and the per-chain-node depleted-REGION mask.
    """
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    node_rtype, _ = _node_region_type(chain, region_arrays)
    floor_mask = is_reg & ((node_rtype == 0) | (node_rtype == 1))  # intergenic + intron REGIONS
    EGl = np.maximum(np.asarray(geometry.eff_gdna_left, np.float64), _EPS)  # region: contained gDNA eff-len
    Ml = np.asarray(geometry.mass_left, np.float64)                         # region: contained mass
    # strand-deconvolved gDNA density: discount known-RNA introns where the strand is informative.
    w_str = float((2.0 * kappa - 1.0) ** 2)
    gdna_frac = (1.0 - w_str) + w_str * np.asarray(f_g_init, np.float64)  # ∈ (0,1]; =1 unstranded/intergenic
    dens_g = (Ml / EGl) * gdna_frac
    eff = EGl[floor_mask]
    g_dens = dens_g[floor_mask]
    # Exposure-pooled floor rate (Poisson–Gamma, 1 pseudocount on the TOTAL gDNA): zero-gDNA depleted
    # regions are KEPT (0 gDNA over a large E is the strongest evidence gDNA is scarce → ρ_floor→0⁺).
    G = float(np.sum(g_dens * eff))
    E_tot = max(float(np.sum(eff)), _EPS)
    rho_floor = (1.0 + G) / E_tot
    var_mean_floor = 1.0 / (1.0 + G)
    # Between-region SPREAD of log gDNA-density over the POPULATED floor regions (eff-weighted population
    # variance minus the per-region log-Poisson floor → the excess/biological spread, ≥0). Tight for a
    # coherent depleted population (a confident floor); naturally widens on real data (GC/mappability).
    pop = g_dens > 0.0
    if int(np.sum(pop)) >= 2:
        lr = np.log(g_dens[pop])
        w = eff[pop]
        sw = float(np.sum(w))
        mu = float(np.sum(w * lr) / sw)
        s2_raw = float(np.sum(w * (lr - mu) ** 2) / sw)
        pois = float(np.mean(1.0 / (g_dens[pop] * eff[pop] + 1.0)))
        s2_floor = max(s2_raw - pois, 0.0)
    else:
        s2_floor = 0.0
    return rho_floor, s2_floor, var_mean_floor, floor_mask


def _global_logprior(
    fgg, mass_global, eff_global, rho_global, sigma2_g, var_mean,
    *, floor_mask=None, rho_floor=None, s2_floor_total=None,
):
    """Precompute the count-space global as a per-node BINOMIAL pseudo-count on f_g, ``(n_nodes, P)`` (§4):

        glob(f_g) = α·log f_g + (N − α)·log(1 − f_g),   α = N·μ,   μ = clip(ρ·E/M, 0, 1),   N = ρ² / σ²_prior

    The genome-wide baseline uses ``ρ = ρ_global`` and ``σ²_prior = var_mean + σ²_g(ρ_global)``; the
    depleted intergenic/intron floor override (``floor_mask``/``rho_floor``) pins that population to its
    confident floor.

    ``N`` is the **M-INDEPENDENT** population confidence (``= 1/CV²`` of the rate) — a hyperprior is naturally
    imprecise, so it can NEVER overrule a node's own (strand) evidence; the MEAN ``μ`` keeps the honest
    density→fraction ``E/M`` conversion. ``σ²_prior`` = the rate-estimate variance ``var_mean = (1+G)/E_tot²``
    (which carries the 1-pseudocount floor, so ``N → 1`` at uniform zero-gDNA — a gentle zero-baseline, never
    ``0/0``) PLUS the between-node spread ``σ²_g`` (large under capture ⇒ ``N → ρ²/σ²_g ≈`` small, imprecise;
    small under uniform present gDNA ⇒ ``N`` large, confident). Two-sided (mode μ); applied to ALL solvable
    nodes (the strand_obs fork is dissolved)."""
    s2_between = max(float(sigma2_g.predict(np.array([max(rho_global, _EPS)]))[0]), 0.0)
    s2_flat = var_mean + s2_between
    # ── LOG-DENSITY global: a Gaussian on log f_g (D-plan §1.4). var_mean = 1/(1+G) (D4); the
    #    M-independent precision is N_log = 1/Var(log ρ) directly — NO ρ² factor (s2_flat is ALREADY a
    #    log-variance, not a density variance). target = log(implied gDNA fraction). ──
    eff = np.maximum(eff_global, _EPS)
    mass = np.maximum(mass_global, _EPS)
    n_glob = 1.0 / max(s2_flat, _EPS)
    target = np.log(np.clip(rho_global * eff / mass, _EPS, 1.0))
    n_node = np.full(target.shape, n_glob, dtype=np.float64)
    # DEPLETED-REGION FLOOR (`_floor_estimate`): intergenic + intron REGION nodes are a coherent depleted
    # population, so override them with the CONFIDENT floor (ρ_floor at the tight floor spread). target =
    # ρ_floor·E/M ⇒ a floor-level intron pins to f_g≈1; an intron with density EXCESS over the floor gets
    # nascent (low f_g) — the nascent-from-self-evidence principle, no gate. Exons/boundaries keep the
    # genome-wide baseline. On real data the floor and a node's own strand evidence agree (same physical
    # density); they conflict only in the documented all-RNA-intron case the floor assumption excludes.
    if floor_mask is not None and rho_floor is not None and s2_floor_total is not None:
        fm = np.asarray(floor_mask, bool)
        target[fm] = np.log(np.clip(float(rho_floor) * eff[fm] / mass[fm], _EPS, 1.0))
        n_node[fm] = 1.0 / max(float(s2_floor_total), _EPS)
    # STABILITY-ONLY: cap the WHOLE global (flat + floor) at one pseudo-observation so it cannot drive or
    # drag any node — the single-strand solve is carried by strand + messages, the global only keeps
    # low/zero-gDNA nodes finite. (The target VALUES — floor for depleted, ρ_global elsewhere — stay; only
    # their weight is capped.)
    n_node = np.minimum(n_node, _GLOBAL_STAB_PREC)
    log_fg = np.log(np.maximum(np.asarray(fgg, np.float64), _EPS))  # (K,)
    return -0.5 * n_node[:, None] * (log_fg[None, :] - target[:, None]) ** 2


def _kde_logprior(fgg, mass_global, eff_global, gdna_prior):
    """The GENERATIVE two-density prior term ``(n_nodes, K)`` on the f_g solve grid (design:
    ``ambig_boundary_spliced_deconvolution.md``; derived by the density-prior-integration workflow).

    The node's total density ``d = M/E`` splits into gDNA density ``ρ_g = f_g·d`` and RNA density
    ``ρ_r = (1−f_g)·d``. Two independent density priors:

      * gDNA: the empirical population KDE ``P(log ρ_g)`` — evaluated with **real Gaussian tails**
        (:meth:`GdnaDensityPrior.logpdf_kernel`, NOT the clamped interpolation, whose constant tail lets a
        high-density node drift to ``f_g≈0.5`` → false-positive gDNA).
      * RNA: a scale-free **Jeffreys** prior ``p(ρ_r) ∝ 1/ρ_r`` (RNA spans >10⁴× — no informative scale).
        Its Jacobian into the f_g coordinate is exactly ``1/(1−f_g)`` ⇒ the ``−log(1−f_g)`` term below.

    The Jeffreys term is the crux: without it gDNA is priored but RNA is free, so the cheapest explanation of
    any flat-strand node is "dump mass into free RNA, park gDNA at the tall depleted KDE mode" — the cliff and
    the false-positives. With it, lowering f_g raises ρ_r (penalised), so gDNA is the residual after a
    typical-magnitude RNA. Both terms are on the f_g grid only (the Jeffreys is node-independent, O(K)).
    NO tuned constants (the Jeffreys exponent is 1, the KDE coordinate is native-log)."""
    eff = np.maximum(np.asarray(eff_global, np.float64), _EPS)
    mass = np.maximum(np.asarray(mass_global, np.float64), _EPS)
    log_me = np.log(mass) - np.log(eff)                                # (m,) = log(M/E)
    fg = np.minimum(np.maximum(np.asarray(fgg, np.float64), _EPS), 1.0 - _EPS)  # (K,)
    log_rho = np.log(fg)[None, :] + log_me[:, None]                    # (m,K) = log ρ_g at each grid point
    kde_term = gdna_prior.logpdf_kernel(log_rho.ravel()).reshape(log_rho.shape)  # real quadratic tails
    jeffreys = -np.log1p(-fg)                                          # (K,) RNA Jeffreys 1/(1−f_g)
    return kde_term + jeffreys[None, :]


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
    logodds_window: float = 10.0,
    n_tilt: int | None = None,
    gdna_prior=None,
    _capture: dict | None = None,
):
    """The belief-propagation sweep over the chain — gDNA AND per-strand RNA messages, in COUNT space (no
    ``(M/E)²`` Jacobian) — as a single **FORWARD-BACKWARD** pass.

    The chain is a forest of linear paths, so BP is exact in one forward + one backward pass (vs Gauss-Seidel /
    Jacobi which propagate one hop per pass). The message precision is the per-edge **disagreement-aware**
    ``σ²_edge`` (`dispersion_aware_message_precision.md`; built inside :func:`_scan`), NOT a fitted ``σ²_bio(μ)``
    var~mean curve — so there is no precision to refit and no outer fixed-point loop. The global prior is
    ANCHORED (every input fit once before the solve), so the single FB pass is exact.

    BEFORE the pass: the non-circular population gDNA prior on gDNA-clean seeds (:func:`_gdna_seed_estimate`) —
    exposure-pooled ``ρ_global`` + the seed gDNA spread ``σ²_g`` (the GLOBAL prior's between-node spread, NOT a
    message reliability) — and the enrichment transfer ``ê(z)``; both anchored.

    The pass: (A) a batched message-free LOCAL solve → per-component (fraction, precision) — also the
    disagreement anchor ``lf*_loc`` / ``v*_loc``; (B) a FORWARD scan L→R accumulating the left-context message α
    and the forward belief (local ⊗ incoming message — NOT the reverse, true tree BP — so a thin node passes the
    upstream α on: the relay); (C) a BACKWARD scan R→L → β; (D) combine α⊗β and one batched FINAL solve.

    The per-node ψ solve integrates: the strand likelihood, the Jeffreys reference (single-strand nodes), the
    count-space global NB prior (ALL solvable nodes — the fork is dissolved), and the gDNA + per-strand RNA
    imputation messages (the spliced/MATURE floor is a node-local term in the per-node solver). The emission gate
    is a three-term Boolean over the components (gDNA / +RNA / −RNA): each RNA strand flows only where that strand
    is continuous on both endpoints (``free_s``); gDNA is genomically continuous and strand-agnostic, so it flows
    wherever there is unspliced facing mass — including across TSS/TES and other G1 seams (a locked all-gDNA node
    is a confident emitter). Only G2/G3 nodes with data are written; G1 sinks + empty nodes keep their init; a
    node with no neighbour message is solved on its own strand+global evidence (the intrinsic solve folds into the
    final batched solve at prec=0). ``max_passes``/``convergence_delta`` are retained for API compatibility but
    unused (FB is single-pass). Returns the converged :class:`NodeBelief`."""
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
    ESP = (geometry.eff_spl_left, geometry.eff_spl_right)  # one-sided spliced half-triangle eff-len
    MS = (geometry.mass_left, geometry.mass_right)
    SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)
    SN = (geometry.spliced_neg_left, geometry.spliced_neg_right)
    # per-node "global" gDNA support: region = its contained support; boundary = the AVERAGED per-side density
    # length ½(E_l+E_r) over the total crossing mass.
    eff_global = np.where(is_reg, EG[0], 0.5 * (EG[0] + EG[1]))
    mass_global = np.where(is_reg, MS[0], MS[0] + MS[1])

    # The per-node solve is the log-density 1-D/2-D log-odds solver (simplex_logodds, O(m·K),
    # genome-scale-tractable). The "solve grid" is the f_g axis the global NB prior is evaluated on (the
    # log-odds σ(λ) lattice).
    _lam_lo, _fg_lo = _logodds_grid(int(n_grid), float(logodds_window))
    solve_grid = _fg_lo
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    _zero_spl = np.zeros_like(np.asarray(statics.mass_unspliced, dtype=np.float64))

    def _local_solve(g_arr, gm=None, gp=None, rm=None, rp=None):
        """The per-node local/final solve (log-density log-odds backend). Returns
        ``(f_g median, f_pos mean, f_neg mean, var_g, var_pos, var_neg)``. Phase A calls it message-free;
        phase D passes the FB messages (fraction-space)."""
        dc = _solve_nodes_logodds_all(
            statics.u_pos, statics.u_neg, statics.spliced_pos, statics.spliced_neg, fp, fn,
            statics.strand_obs, statics.mass_unspliced, _zero_spl, kappa=kappa, od_g=od_g, od_r=od_r,
            n_grid=int(n_grid), L=float(logodds_window), n_tilt=n_tilt, global_logprior=g_arr,
            gdna_imp_mode=gm, gdna_imp_prec=gp, rna_imp_mode=rm, rna_imp_prec=rp,
        )
        return (dc.gdna_frac, dc.rna_pos_frac, dc.rna_neg_frac,
                dc.gdna_frac_var, dc.rna_pos_frac_var, dc.rna_neg_frac_var)
    # Two distinct gates, both from the region SIGNATURE (never the counts — count-zero-info):
    #   * SOLVE gate (`solvable`): a node deconvolves its own gDNA-vs-RNA split iff it admits ≥1 RNA strand and
    #     has unspliced mass. A G1 node — no admissible RNA strand: an intergenic region, or a gene-boundary seam
    #     (TSS/TES, opposite-strand exon↔exon) — is a LOCKED all-gDNA node ({0,0,1}); it is not solved, keeping
    #     its init (RNA cannot cross a gene boundary, so its unspliced mass is purely gDNA).
    #   * EMISSION gate (per component, in `_scan`): which MESSAGES a node sends. A three-term Boolean over the
    #     components gDNA / +RNA / −RNA, structural and symmetric — defined at the top of the scan loop.
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)

    # The NON-CIRCULAR population gDNA prior on gDNA-clean seeds (§4.3), ANCHORED — computed ONCE, never refit:
    # the exposure-pooled rate ρ_global (the global baseline) + var_mean (its variance) + the seed gDNA spread
    # σ²_g (``gdna_vm`` — the GLOBAL prior's between-node spread, NOT a message reliability). Strand seeds use the
    # strand-ONLY init f_g. (Message reliability is now per-edge in `_scan` — there is no var~mean message fit.)
    rho_global, gdna_vm, var_mean = _gdna_seed_estimate(
        chain, statics, geometry, region_arrays, boundary_substrate, f_g, kappa
    )
    # The DEPLETED gDNA FLOOR from intergenic + intron REGIONS (the user's empirical-prior directive): a
    # confident floor that pins depleted intron nodes (the nascent-hallucination fix), letting an intron
    # with density excess over the floor carry nascent. Exons/boundaries keep the genome-wide global.
    rho_floor, s2_floor, var_mean_floor, floor_mask = _floor_estimate(
        chain, geometry, region_arrays, f_g, kappa
    )
    # The global gDNA prior on f_g, an (n_nodes, K) log-term on the solve grid, applied to ALL nodes. It is
    # ALWAYS the weak, exposure-pooled stability floor (`_global_logprior`: M-INDEPENDENT, capped at one
    # pseudo-observation, so it never overrules a node's own strand evidence; the depleted-floor override).
    # PASS 2 ADDS the trained Phase-2 mixture ON TOP (`_kde_logprior`). This layering is the design's
    # "floor + mixture" (§3.4): the pooled floor is the always-present "gDNA is scarce" baseline — its
    # downward pull is what resolves a balanced AMBIG node (whose two-strand RNA is strand-degenerate with
    # gDNA) toward ~0 in a gDNA-poor context — while the KDE supplies the enriched-mode structure the floor
    # lacks. Replacing the floor with the KDE (an earlier attempt) lost the downward anchor, because the
    # per-node KDE's RNA-residual mode sits ABOVE the exposure-pooled ρ_global (dissection: gDNA=0 AMBIG
    # went 0.02→0.38). ANCHORED — every input is fit once, so the prior is CONSTANT within a pass.
    # PASS 1 (gdna_prior=None): the weak stability floor only (strand + messages carry the single-strand solve;
    # the floor keeps low/zero-gDNA nodes finite and anchors the depleted population that trains the KDE).
    global_lp = _global_logprior(
        solve_grid, mass_global, eff_global, rho_global, gdna_vm, var_mean,
        floor_mask=floor_mask, rho_floor=rho_floor, s2_floor_total=var_mean_floor + s2_floor,
    )
    # PASS 2 (gdna_prior set): ADD the generative two-density prior — the empirical gDNA-density KDE (real
    # tails) × the Jeffreys RNA prior 1/(1−f_g) (`_kde_logprior`). This is the density-prior INTEGRATION: the
    # strand landscape (ψ, added in the solve) × the population landscape. The Jeffreys term removes the
    # gDNA-priored/RNA-free asymmetry that caused the boundary cliff + the false-positive gDNA; the real KDE
    # tails stop a high-density node drifting to f_g≈0.5. Derivation + validation:
    # ambig_boundary_spliced_deconvolution.md. Applied to ALL solvable nodes (self-scaling: a confident strand
    # dominates ψ, an AMBIG/thin node leans on the population).
    if gdna_prior is not None:
        global_lp = global_lp + _kde_logprior(solve_grid, mass_global, eff_global, gdna_prior)

    # Genomic node order for the forward/backward scans (within each ref path; left/right break at −1). The
    # scans are sequential per node, so iterate as a Python list of ints (faster than numpy scalar indexing).
    order_list = [int(x) for x in np.asarray(chain.order)]
    n_nodes = f_g.shape[0]

    # FORWARD-BACKWARD solve — ONE exact pass on the chain (a forest of linear paths) given the message-free
    # local beliefs + the ANCHORED global prior. (A) message-free LOCAL beliefs (one batched solve) →
    # per-component (fraction, precision); (B) FORWARD scan L→R accumulating the left-context message α + the
    # forward belief; (C) BACKWARD scan R→L → β; (D) combine α⊗β and one batched FINAL solve. Each node sees the
    # WHOLE-chain context (relay) in one pass — a thin node's forward belief is dominated by, and passes on, the
    # upstream α. ``true`` tree BP: the forward belief excludes β (no double-counting). The message precision is
    # per-edge (disagreement-aware, in `_scan`) so there is no var~mean fixed point ⇒ no outer loop.

    # (A) LOCAL message-free beliefs (backend-dispatched).
    fg_loc, fp_loc, fn_loc, vg_loc, vp_loc, vn_loc = _local_solve(global_lp)
    pg_loc = 1.0 / np.maximum(vg_loc, _EPS)  # local precision (var floored: a sharp belief ⇒ large finite)
    pp_loc = 1.0 / np.maximum(vp_loc, _EPS)
    pn_loc = 1.0 / np.maximum(vn_loc, _EPS)


    def _scan(seq, nbr, sf, df):
        """Sequential scan: project the running belief from each node's ``nbr`` (src face ``sf`` → dst face
        ``df``) into the dst, then fold it into the dst's running belief (local ⊗ incoming message; NOT the
        reverse message — true tree BP). Returns the per-node message (mode, prec) per component. O(m).

        DISAGREEMENT-AWARE precision (``dispersion_aware_message_precision.md``): a message's reliability is
        the per-edge process variance ``σ²_edge = max(resid² − (base_var + var_loc[dst]), 0)``, where
        ``resid = mo − logf_loc[dst]`` is the message's surprise vs the dst's MESSAGE-FREE local belief (the
        non-circular anchor — §4.2) and ``base_var = var[src] + pois`` is the sampling + source-belief
        log-variance. ``prec = 1/(base_var + σ²_edge)``: an agreeing edge (``resid≈0``) keeps full precision
        (the smooth relay), a seam edge (large ``resid``) self-silences (``prec ≈ 1/resid²``). Replaces the
        retired ``σ²_bio(μ)`` var~mean curve. Same formula for all three components (gDNA / ±RNA)."""
        fbg, fbp, fbn = fg_loc.copy(), fp_loc.copy(), fn_loc.copy()  # running belief (starts at local)
        vbg, vbp, vbn = vg_loc.copy(), vp_loc.copy(), vn_loc.copy()
        amg, apg = np.zeros(n_nodes), np.zeros(n_nodes)  # gDNA message (mode, prec)
        amp, app = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-pos
        amn, apn = np.zeros(n_nodes), np.zeros(n_nodes)  # RNA-neg
        EGs, EGd, ERs, ERd = EG[sf], EG[df], ER[sf], ER[df]
        MSs, MSd, SPs, SNs = MS[sf], MS[df], SP[sf], SN[sf]
        ESPs = ESP[sf]  # source-face spliced eff-len (for the mature-RNA MEASUREMENT message)
        SPd, SNd, ESPd = SP[df], SN[df], ESP[df]  # DEST-face spliced — the mature ABSORBED at a junction
        # (subtracted from an exon→boundary message so only NASCENT crosses into the intron side).
        # The running belief combines in LOG-fraction space (the message is a Gaussian on log f_c).
        # Precompute the local log-fractions (constant across the scan) for the combine.
        lfg_loc = np.log(np.maximum(fg_loc, _EPS))
        lfp_loc = np.log(np.maximum(fp_loc, _EPS))
        lfn_loc = np.log(np.maximum(fn_loc, _EPS))
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
                n_src = fbg[lsrc] * sm                   # source gDNA COUNT (deconvolved)
                rho = n_src / eg                         # source gDNA DENSITY ρ_g_src — the MESSAGE currency
                # DENSITY message (NO fractions in the wire): the content is the source gDNA density
                # ρ_g_src; the RECEIVER re-expresses it in its OWN log-f_g solve frame via its gDNA
                # density base M_dst/E_gdna_dst (= md/egd), flooring ρ at the dst min-observable density
                # 1/egd. (Fractions are not comparable across nodes — only densities are.)
                mo = math.log(max(rho, 1.0 / egd) / (md / egd))
                # pois = 1/(gDNA COUNT): the source's log-density sampling variance. NO +1 floor — a
                # NON-DETECTION (count→0, a depleted seam) has ~∞ log-density variance → base_var huge →
                # ~0 precision → sends NO message. "Zero density is not a measurement." (user's diagnosis)
                pois = 1.0 / max(n_src, _EPS)
                base_var = vbg[lsrc] + pois              # sampling + source-belief log-var
                s2_edge = max((mo - lfg_loc[i]) ** 2 - (base_var + vg_loc[i]), 0.0)  # per-edge surprise
                pr = 1.0 / max(base_var + s2_edge, _EPS)
                amg[i], apg[i] = mo, pr
                pt = pg_loc[i] + pr
                fbg[i] = math.exp((pg_loc[i] * lfg_loc[i] + pr * mo) / pt)
                vbg[i] = 1.0 / pt
            # RNA-pos — the imputed RNA density is the NASCENT (contiguous, unspliced) field, the only RNA
            # that crosses an edge contiguously. MATURE does NOT flow as an imputation (it skips introns as
            # SPLICED); it enters/leaves only at a junction boundary, via three density terms that sum to the
            # message (the one-sided spliced geometry makes the gating automatic — `spliced_efflen_not_2x_…`):
            #   ρ = (src nascent  ``fbp·sm/E_r``)
            #     + (src-face mature  ``SPs/E_spl_src``)   — ADDED only when the SOURCE is a junction boundary
            #         facing the dst exon (B→exon): a MEASUREMENT of the exon's within-exon mature.
            #     − (dst-face mature  ``SPd/E_spl_dst``)   — ABSORBED only when the DEST is a junction boundary
            #         facing the src exon (exon→B): the exon's mature crosses as spliced, so it must NOT be
            #         imputed into the boundary's (gDNA+nascent) crossing — subtracting it leaves pure nascent.
            # At most one of SPs/SPd is non-zero on an edge (src and dst are never both boundaries), and BOTH
            # are zero on intron↔boundary edges (spliced lives on the exon flank), so introns carry pure
            # nascent with NO gate. Precision: ALL RNA messages (the mature MEASUREMENT and the nascent
            # IMPUTATION alike) use the disagreement-aware ``σ²_edge``. This is self-targeting and exactly
            # realizes the "a confident strand wins, a strand-blind exon snaps" intent: when the message
            # AGREES with the dst's message-free local belief, ``s2_edge = 0`` ⇒ full count precision
            # (unchanged), so a strand-blind exon (uncertain ``f_pos``, large ``vp_loc`` absorbs the residual)
            # still snaps to the mature it observes; when the message DISAGREES with a CONFIDENT local belief
            # (a strong-strand exon whose own deep unspliced reads already pin ``f_pos``), ``s2_edge`` down-
            # weights it. The old exemption (count precision, never silenced) let a capture-biased mature
            # density (junction-spanning reads are only partially captured ⇒ ``n_mat/E_spl`` under-estimates
            # the exon's RNA) OVERRIDE a correct strand-confident exon → phantom gDNA by simplex complement
            # (−gDNA flagship +0.04→+0.018; toy no-gDNA +0.20→+0.005; expressed exons unchanged).
            if emit_p:
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                n_nasc = fbp[lsrc] * sm              # source nascent RNA count (unspliced)
                n_mat = SPs[lsrc]  # source-face mature (>0 only B→exon): MEASURES the exon's mature
                rho_mat_dst = SPd[i] / (ESPd[i] if ESPd[i] > _EPS else _EPS)  # dst-face mature absorbed (exon→B)
                rho = n_nasc / er + n_mat / esp - rho_mat_dst  # NASCENT density (+ MEASUREMENT into an exon)
                mo = math.log(max(rho, 1.0 / erd) / (md / erd))   # → dst log-f_pos frame (floored at min-observable)
                pois = 1.0 / max(n_nasc + n_mat, _EPS)
                base_var = vbp[lsrc] + pois
                s2_edge = max((mo - lfp_loc[i]) ** 2 - (base_var + vp_loc[i]), 0.0)
                pr = 1.0 / max(base_var + s2_edge, _EPS)          # disagreement-aware (MEASUREMENT or IMPUTATION)
                amp[i], app[i] = mo, pr
                pt = pp_loc[i] + pr
                fbp[i] = math.exp((pp_loc[i] * lfp_loc[i] + pr * mo) / pt)
                vbp[i] = 1.0 / pt
            # RNA-neg — symmetric (mature on the −strand junction motif; same 3-term nascent message).
            if emit_n:
                er = ERs[lsrc] if ERs[lsrc] > _EPS else _EPS
                esp = ESPs[lsrc] if ESPs[lsrc] > _EPS else _EPS
                n_nasc = fbn[lsrc] * sm
                n_mat = SNs[lsrc]
                rho_mat_dst = SNd[i] / (ESPd[i] if ESPd[i] > _EPS else _EPS)  # dst-face mature absorbed (exon→B)
                rho = n_nasc / er + n_mat / esp - rho_mat_dst
                mo = math.log(max(rho, 1.0 / erd) / (md / erd))   # → dst log-f_neg frame
                pois = 1.0 / max(n_nasc + n_mat, _EPS)
                base_var = vbn[lsrc] + pois
                s2_edge = max((mo - lfn_loc[i]) ** 2 - (base_var + vn_loc[i]), 0.0)
                pr = 1.0 / max(base_var + s2_edge, _EPS)          # disagreement-aware (MEASUREMENT or IMPUTATION)
                amn[i], apn[i] = mo, pr
                pt = pn_loc[i] + pr
                fbn[i] = math.exp((pn_loc[i] * lfn_loc[i] + pr * mo) / pt)
                vbn[i] = 1.0 / pt
        return amg, apg, amp, app, amn, apn

    a = _scan(order_list, left, 1, 0)         # forward (α: left context)
    b = _scan(order_list[::-1], right, 0, 1)  # backward (β: right context)

    # (D) combine α⊗β (precision-weighted product) per component → one batched FINAL solve.
    def _comb(am_a, ap_a, am_b, ap_b):
        pc = ap_a + ap_b
        return np.where(pc > _EPS, (ap_a * am_a + ap_b * am_b) / np.maximum(pc, _EPS), 0.0), pc

    mode_g, prec_g = _comb(a[0], a[1], b[0], b[1])
    mode_p, prec_p = _comb(a[2], a[3], b[2], b[3])
    mode_n, prec_n = _comb(a[4], a[5], b[4], b[5])
    # (D) FINAL solve with the FB messages (backend-dispatched).
    mg_, mp_, mn_, vg_, vp_, vn_ = _local_solve(
        global_lp, mode_g, prec_g, (mode_p, mode_n), (prec_p, prec_n)
    )
    # write back only SOLVABLE nodes (G1 sinks / empty keep their signature-binary init).
    f_g = np.where(solvable, np.clip(mg_, 0.0, 1.0), f_g)
    f_pos = np.where(solvable, np.clip(mp_, 0.0, 1.0), f_pos)
    f_neg = np.where(solvable, np.clip(mn_, 0.0, 1.0), f_neg)
    var_g = np.where(solvable, vg_, var_g)
    var_pos = np.where(solvable, vp_, var_pos)
    var_neg = np.where(solvable, vn_, var_neg)

    if _capture is not None:  # inert diagnostic hook
        # strand-ONLY local belief (no global prior, no messages) — to split the LOCAL error into the
        # strand likelihood vs the global gDNA prior contribution. Same log-density solver, global=None.
        fg_strand = _solve_nodes_logodds_all(
            statics.u_pos, statics.u_neg, statics.spliced_pos, statics.spliced_neg, fp, fn,
            statics.strand_obs, statics.mass_unspliced, _zero_spl, kappa=kappa, od_g=od_g, od_r=od_r,
            n_grid=int(n_grid), L=float(logodds_window), n_tilt=n_tilt, global_logprior=None,
        ).gdna_frac
        _capture.update(
            fg_loc=fg_loc, fg_strand=fg_strand, fp_loc=fp_loc, fn_loc=fn_loc,
            vg_loc=vg_loc, vp_loc=vp_loc, vn_loc=vn_loc,
            a_fwd=a, b_bwd=b, mode_g=mode_g, prec_g=prec_g, mode_p=mode_p, prec_p=prec_p,
            mode_n=mode_n, prec_n=prec_n, f_g=f_g.copy(), f_pos=f_pos.copy(), f_neg=f_neg.copy(),
            var_g=var_g.copy(), solvable=solvable, rho_global=rho_global, gdna_vm=gdna_vm,
            # boundary-emission geometry: gDNA emits iff facing unspliced mass>0 (strand-agnostic);
            # RNA iff free_s on both endpoints & (unspliced or spliced facing mass). Capture the faces.
            mass_l=MS[0], mass_r=MS[1], spl_l=SP[0] + SN[0], spl_r=SP[1] + SN[1],
            free_pos=np.asarray(fp, bool), free_neg=np.asarray(fn, bool),
            # global geometry (μ = clip(ρ·eff_global/mass_global, 0, 1) is the implied prior fraction).
            eff_global=eff_global, mass_global=mass_global,
            # per-face geometry for message dissection (logodds diagnostics)
            eff_gdna_l=EG[0], eff_gdna_r=EG[1], eff_rna_l=ER[0], eff_rna_r=ER[1],
            rho_floor=rho_floor, floor_mask=floor_mask,
        )

    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g,
        var_pos=var_pos, var_neg=var_neg, var_gdna=var_g,
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
