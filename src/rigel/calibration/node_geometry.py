"""rigel.calibration.node_geometry — per-node chain geometry, beliefs, densities, statics, and init.

The lower layer beneath the belief-propagation solver (`bp_solver`): everything that describes *what a node
is* before any message passing. Pure functions of the accumulator payload + the region/boundary substrates +
the FL pmfs — no sweep state, no global prior.

* `NodeGeometry` / `build_node_geometry` — the per-node, per-FACE static geometry (unspliced mass,
  per-component gDNA-/RNA-FL effective lengths, the boundary's one-sided motif-stranded spliced mass). Per-face
  because a boundary's two sides lie in different-sized regions and a message uses the side FACING the
  destination; a region presents the same contained geometry both ways.
* `NodeBelief` — the per-node pie ``(f_pos, f_neg, f_g)``; the per-component message densities
  ``rho = f·M/E`` are computed inline in ``bp_solver.node_sweep``.
* `NodeStatics` / `build_node_statics` — the static per-node solver inputs (per-strand counts, masks, masses).
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.
* `_node_region_type` — per-chain-node coarse region type (intergenic/intron/exon), shared with the prior
  subsystem and `gdna_density_prior`.

Layering: imports only `node_chain`, `signature`, `effective_length`, `simplex_logodds`, `strand_deconv`
(all lower layers) — never `bp_solver` or `gdna_density_prior`, so it sits cleanly below both (dissolving the
former `gdna_density_prior` → `bp_solver` back-import).
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
    TS_AMBIG,
    TS_NEG,
    TS_POS,
    coarse_type_array,
    mrna_active_strands,
    nrna_active_strands,
)
from .simplex_logodds import _solve_nodes_logodds_all

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "node_global_geometry",
    "node_total_density",
    "NodeBelief",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
]

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
    # The matching integer UNSPLICED COUNT (flux), same faces as the mass above. MASS is the density
    # numerator; COUNT is what a Poisson VARIANCE needs (`Var(log ρ)=1/n`, not `1/mass`) — see the
    # `spliced_n_*` note below, which states the same rule for the spliced channel. The unspliced mass is
    # additionally NOT conserved per node (the accumulator splits one fragment across several nodes), so
    # `1/mass` is wrong in two independent ways; the integer flux is the honest ``n``.
    n_unspl_left: np.ndarray
    n_unspl_right: np.ndarray
    eff_gdna_left: np.ndarray  # gDNA-FL eff-len facing left (contained / boundary left side)
    eff_gdna_right: np.ndarray
    eff_rna_left: np.ndarray  # RNA-FL eff-len facing left (nascent two-sided crossing / contained)
    eff_rna_right: np.ndarray
    eff_spl_left: (
        np.ndarray
    )  # one-sided spliced-crossing RNA eff-len (half-triangle E[min²/2ℓ]) facing left
    eff_spl_right: np.ndarray
    spliced_pos_left: np.ndarray  # + motif spliced MASS on the left face (boundary, exon-on-left)
    spliced_pos_right: np.ndarray
    spliced_neg_left: np.ndarray
    spliced_neg_right: np.ndarray
    # The matching integer spliced COUNT (flux), same faces, same motif/exon routing as the mass above.
    # MASS is the correct numerator for a DENSITY (`spliced/eff_spl`); the COUNT is what a Poisson VARIANCE
    # needs (`Var(log rho_m) = 1/n`, not `1/mass` — mass sums fractional per-fragment shares, so `1/mass`
    # over-states the variance; the Kish n_eff = (sum w)^2/sum w^2 >= mass). Precedent for carrying both:
    # `_boundary_strand_stats` already takes the integer flux for the BB strand power and the float mass for
    # the pie base. Currently UNUSED by the solver — plumbed for the priority-#3 mature-measurement channel
    # (docs/calibration/boundary_spliced_channel_design.md §4.1); inert until that lands.
    spliced_n_pos_left: np.ndarray
    spliced_n_pos_right: np.ndarray
    spliced_n_neg_left: np.ndarray
    spliced_n_neg_right: np.ndarray


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
    side_eff_spl = spliced_side_eff_length(
        rna_fl_pmf, L
    )  # one-sided spliced-crossing half-triangle
    reg_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    R = L.shape[0]

    # --- per-boundary quantities (two sides) ---
    bsub = boundary_substrate
    blr = np.asarray(bsub.left_region, dtype=np.int64)
    brr = np.asarray(bsub.right_region, dtype=np.int64)
    bmass_l = np.asarray(bsub.left.mass_unspliced, dtype=np.float64)
    bmass_r = np.asarray(bsub.right.mass_unspliced, dtype=np.float64)
    # matching integer unspliced flux per side (the Poisson n for the message precision)
    bn_unspl_l = np.asarray(bsub.left.n_unspliced, dtype=np.float64)
    bn_unspl_r = np.asarray(bsub.right.n_unspliced, dtype=np.float64)
    bspl_l = np.asarray(bsub.left.mass_spliced, dtype=np.float64)  # sense+antisense summed
    bspl_r = np.asarray(bsub.right.mass_spliced, dtype=np.float64)
    # the matching integer flux (same channels, summed the same way) — the Poisson count for the variance
    bspn_l = np.asarray(bsub.left.n_spliced_sense, np.float64) + np.asarray(
        bsub.left.n_spliced_antisense, np.float64
    )
    bspn_r = np.asarray(bsub.right.n_spliced_sense, np.float64) + np.asarray(
        bsub.right.n_spliced_antisense, np.float64
    )
    B = blr.shape[0]
    # a boundary's per-side crossing eff-len = its flank region's E[min(ℓ,L)] (0 at a −1 terminal flank)
    b_eff_g_l = np.where(blr >= 0, side_eff_g[np.clip(blr, 0, R - 1)], 0.0)
    b_eff_g_r = np.where(brr >= 0, side_eff_g[np.clip(brr, 0, R - 1)], 0.0)
    b_eff_r_l = np.where(blr >= 0, side_eff_r[np.clip(blr, 0, R - 1)], 0.0)
    b_eff_r_r = np.where(brr >= 0, side_eff_r[np.clip(brr, 0, R - 1)], 0.0)
    # spliced face: route a boundary side's spliced (mature) mass to that side's EXON flank on the junction's
    # KNOWN genomic strand (``boundary_substrate.junction_strand``, observed from the motif at deposit — one
    # motif-stranded junction per boundary). Correct at AMBIG / exon↔exon seams the signatures cannot orient.
    # One-sided per exon: mass_left → left-exon donor, mass_right → right-exon acceptor.
    sig_l = np.where(blr >= 0, sig[np.clip(blr, 0, R - 1)], 0)
    sig_r = np.where(brr >= 0, sig[np.clip(brr, 0, R - 1)], 0)
    js = np.asarray(boundary_substrate.junction_strand).reshape(
        -1
    )  # TS_POS / TS_NEG / 0 (== Strand)
    # STRAND-AWARE routing: mature RNA is single-stranded, and a splice junction is single-stranded with a
    # KNOWN, data-type-invariant strand (the genomic splice motif). A junction absorbs mature ONLY on its own
    # strand, so its spliced mass routes to a flank carrying THAT STRAND's exon bit — never merely "any exon".
    # At an AMBIG seam (overlapping opposite-strand transcripts) a + junction beside a −-only exon must NOT
    # deposit + mature there, and vice versa. Per strand: TS_POS → BIT_EXON_POS flanks, TS_NEG → BIT_EXON_NEG.
    exon_pos_l = (sig_l & BIT_EXON_POS) != 0
    exon_pos_r = (sig_r & BIT_EXON_POS) != 0
    exon_neg_l = (sig_l & BIT_EXON_NEG) != 0
    exon_neg_r = (sig_r & BIT_EXON_NEG) != 0

    def _spliced_faces(strand_val, exon_l, exon_r, vl, vr):
        on = js == strand_val  # the junction is on this strand → its SAME-STRAND exon flank carries the mature
        return np.where(on & exon_l, vl, 0.0), np.where(on & exon_r, vr, 0.0)

    b_spl_pos_l, b_spl_pos_r = _spliced_faces(TS_POS, exon_pos_l, exon_pos_r, bspl_l, bspl_r)
    b_spl_neg_l, b_spl_neg_r = _spliced_faces(TS_NEG, exon_neg_l, exon_neg_r, bspl_l, bspl_r)

    # ── the ONE-SIDED SPLICED eff-length: which divisor, per face (A1/A2, message_layer_derivation.md) ──
    # A spliced fragment credits only its exon flank. Enumerating the accumulator's deposit rule
    # ((slice_len/ℓ)/n_cross) over the fragment's ``a`` bases on this side, with flank length ``R``:
    #     a ≤ R  → the near slice lies inside the flank      → END slice,      n_cross=1 → a/ℓ
    #     a > R  → it overruns into the NEXT genomic region  → INTERIOR slice, n_cross=2 → R/(2ℓ)
    # so the per-unit-density deposit is
    #     Σ_a = R²/(2ℓ) + (ℓ−R)·R/(2ℓ) = min(ℓ,R)/2     when the exon CONTINUES past the flank's far edge
    #     Σ_a = R²/(2ℓ)                                  when the exon TERMINATES there (no fragment can overrun)
    # Brute-force verified to 4 d.p. against both closed forms. The half-triangle
    # (:func:`spliced_side_eff_length`) is therefore correct ONLY for a terminal exon; for every INTERNAL one it
    # is low by exactly ``ℓ/R`` (measured 2.0× at R=100/ℓ=200, 8.0× at R=25, 12.0× at R=25/ℓ=300), and the
    # correct divisor is ``E[min(ℓ,R)]/2`` — which is :func:`boundary_side_eff_length`, i.e. ``side_eff_r``,
    # already computed above on the RNA FL.
    #
    # The selector is STRUCTURAL and carries no constant: can the fragment's coverage extend more than ``R``
    # bases past the junction — i.e. does the MATURE TRANSCRIPT continue past the flank region's far edge?
    #
    # ⚠ Verified against the accumulator reference (`tests/native/_accumulator_reference.py`, the crossing
    # path): a slice is INTERIOR (``n_cross=2``) iff another slice follows it — ``crosses_right = i < n-1`` —
    # which is a property of the FRAGMENT's slice list, not of region geometry. For a spliced fragment the
    # next block lands in a different region whether the far boundary is an exon↔exon seam OR an exon↔intron
    # splice junction. So an INTERNAL exon continues in BOTH cases; only a genuine transcript END terminates.
    # (A first draft of this used ``mrna_active_s`` alone — exon on both sides — which wrongly called every
    # exon↔intron junction "terminating" and so left almost every internal exon face on the half-triangle.)
    #
    #     continues_s(far boundary)  =  mrna_active_s(far)          exon runs contiguously past it
    #                                OR  far is a splice junction on strand s   the transcript splices onward
    #
    # taken on the JUNCTION's own strand, since that is the strand whose mature is deposited on this face.
    lb_of_reg = np.full(R, -1, dtype=np.int64)  # region → the boundary on its LEFT
    rb_of_reg = np.full(R, -1, dtype=np.int64)  # region → the boundary on its RIGHT
    _m = brr >= 0
    lb_of_reg[brr[_m]] = np.nonzero(_m)[0]  # boundary b is the LEFT boundary of region brr[b]
    _m = blr >= 0
    rb_of_reg[blr[_m]] = np.nonzero(_m)[0]  # ... and the RIGHT boundary of region blr[b]

    def _mrna_active(bit_pos, bit_neg):
        """Per-BOUNDARY 'mature crosses contiguously', per strand: exon bit on BOTH flank regions."""
        lp = np.where(blr >= 0, bit_pos[np.clip(blr, 0, R - 1)], False)
        rp = np.where(brr >= 0, bit_pos[np.clip(brr, 0, R - 1)], False)
        ln = np.where(blr >= 0, bit_neg[np.clip(blr, 0, R - 1)], False)
        rn = np.where(brr >= 0, bit_neg[np.clip(brr, 0, R - 1)], False)
        return lp & rp, ln & rn

    _ex_p_reg = (sig & BIT_EXON_POS) != 0
    _ex_n_reg = (sig & BIT_EXON_NEG) != 0
    _mact_p, _mact_n = _mrna_active(_ex_p_reg, _ex_n_reg)

    def _continues(flank_reg, far_bnd_of_reg):
        """Does the MATURE TRANSCRIPT continue past the flank's far edge, on this junction's strand?"""
        ok = flank_reg >= 0
        fr = np.clip(flank_reg, 0, R - 1)
        fb = np.where(ok, far_bnd_of_reg[fr], -1)
        okb = ok & (fb >= 0)
        fbc = np.clip(fb, 0, B - 1)
        contig = np.where(js == TS_POS, _mact_p[fbc], np.where(js == TS_NEG, _mact_n[fbc], False))
        splices_on = (js != 0) & (js[fbc] == js)  # the far boundary is a junction on the SAME strand
        return okb & (contig | splices_on)

    def _eff_spl_face(flank_reg, far_bnd_of_reg):
        ok = flank_reg >= 0
        fr = np.clip(flank_reg, 0, R - 1)
        cont = _continues(flank_reg, far_bnd_of_reg)
        return np.where(ok, np.where(cont, side_eff_r[fr], side_eff_spl[fr]), 0.0)

    b_eff_spl_l = _eff_spl_face(blr, lb_of_reg)  # left face: its flank's FAR boundary is that region's left
    b_eff_spl_r = _eff_spl_face(brr, rb_of_reg)
    # counts: the SAME gate, so `spliced_n_*` is nonzero exactly where `spliced_*` is
    b_spn_pos_l, b_spn_pos_r = _spliced_faces(TS_POS, exon_pos_l, exon_pos_r, bspn_l, bspn_r)
    b_spn_neg_l, b_spn_neg_r = _spliced_faces(TS_NEG, exon_neg_l, exon_neg_r, bspn_l, bspn_r)

    # --- gather onto the chain (region nodes ← region arrays; boundary nodes ← boundary arrays) ---
    ri_ = np.clip(idx, 0, R - 1)  # region index where is_reg
    bi_ = np.clip(idx, 0, B - 1)  # boundary index where is_bnd

    def pick(reg_vals, bnd_vals):
        return np.where(is_reg, reg_vals[ri_], np.where(is_bnd, bnd_vals[bi_], 0.0))

    # region presents the same (contained) geometry both ways; boundary presents its per-side geometry.
    mass_left = pick(reg_mass, bmass_l)
    mass_right = pick(reg_mass, bmass_r)
    reg_n_unspl = np.asarray(substrate.contained.n_unspliced, dtype=np.float64)
    n_unspl_left = pick(reg_n_unspl, bn_unspl_l)
    n_unspl_right = pick(reg_n_unspl, bn_unspl_r)
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
    spliced_n_pos_left = pick(zeros_R, b_spn_pos_l)
    spliced_n_pos_right = pick(zeros_R, b_spn_pos_r)
    spliced_n_neg_left = pick(zeros_R, b_spn_neg_l)
    spliced_n_neg_right = pick(zeros_R, b_spn_neg_r)

    return NodeGeometry(
        n_nodes=int(chain.n_nodes),
        mass_left=mass_left,
        mass_right=mass_right,
        n_unspl_left=n_unspl_left,
        n_unspl_right=n_unspl_right,
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
        spliced_n_pos_left=spliced_n_pos_left,
        spliced_n_pos_right=spliced_n_pos_right,
        spliced_n_neg_left=spliced_n_neg_left,
        spliced_n_neg_right=spliced_n_neg_right,
    )


def node_global_geometry(chain: NodeChain, geometry: NodeGeometry):
    """Per-node 'global' gDNA support ``(mass, eff)``: a REGION uses its contained mass over its contained
    gDNA eff-length; a BOUNDARY uses its both-side crossing mass over the SUMMED per-side density length
    ``E_l + E_r``. This is the basis the enrichment NPMLE (`DensityNPMLE`) is fit on and projected
    onto — shared by :func:`bp_solver.node_sweep` and ``calibrate`` so the fit and the projection use one
    definition.

    The boundary sum is ``E_l + E_r`` (not the old ``½(E_l+E_r)``) because ``eff_gdna_*`` is now the true
    per-face DENSITY length ``E[min(ℓ,R)]/2`` (`effective_length.boundary_side_eff_length`). The old ½ here
    was silently cancelling the ½ that was *missing* from the face length — which is why this frame read the
    correct ρ while every per-face MESSAGE read ρ/2. Both frames are now the same one."""
    is_reg = np.asarray(chain.kind) == REGION
    egl = np.asarray(geometry.eff_gdna_left, dtype=np.float64)
    egr = np.asarray(geometry.eff_gdna_right, dtype=np.float64)
    msl = np.asarray(geometry.mass_left, dtype=np.float64)
    msr = np.asarray(geometry.mass_right, dtype=np.float64)
    mass = np.where(is_reg, msl, msl + msr)
    eff = np.where(is_reg, egl, egl + egr)
    return mass, eff


def node_total_density(chain: NodeChain, geometry: NodeGeometry, f_g):
    """The LAZY, composition-aware node total density (`unified_solver_design.md` §2, owner 2026-07-23):
    the SUM of component densities, each in its OWN FL frame, from the current belief ``f_g``::

        ρ_unspliced = f_g · (M/E_g^gDNA)  +  (1−f_g) · (M/E_r^RNA)      gDNA-FL for gDNA, RNA-FL for RNA
        ρ_spliced   = spliced_mass / E_spl^RNA                          one-sided, boundary only

    Returns ``(rho_unspliced, rho_with_spliced)`` per node. ``rho_with_spliced`` adds ρ_spliced and is the
    total density used to form the enrichment ratio toward the exon/acceptor side (mature-bearing); the
    mature-free ``rho_unspliced`` is used toward the intron side (§6). This is NEVER a pure-gDNA precompute —
    ``f_g`` is the best current composition (self-solve + messages + measured spliced); gDNA-FL alone
    (``f_g = 1``) is only the fallback where composition is genuinely unknown, and the bounding lemma (§2)
    bounds *that* fallback, not this. Mass/eff are the node-level (both-face-pooled) quantities of
    :func:`node_global_geometry`; the RNA eff-length is the RNA-FL twin, summed the same way."""
    is_reg = np.asarray(chain.kind) == REGION
    mass, eff_g = node_global_geometry(chain, geometry)
    erl = np.asarray(geometry.eff_rna_left, dtype=np.float64)
    err_ = np.asarray(geometry.eff_rna_right, dtype=np.float64)
    eff_r = np.where(is_reg, erl, erl + err_)
    fg = np.clip(np.asarray(f_g, dtype=np.float64), 0.0, 1.0)
    rho_unspl = mass * (fg / np.maximum(eff_g, _EPS) + (1.0 - fg) / np.maximum(eff_r, _EPS))
    # one-sided spliced (mature) DENSITY: spliced mass lands on ONE face, so divide by THAT face's E_spl (a
    # summed-eff divisor would under-state it ~2×). Sum the per-face densities (only the acceptor face is nonzero).
    spl_l = np.asarray(geometry.spliced_pos_left, np.float64) + np.asarray(geometry.spliced_neg_left, np.float64)
    spl_r = np.asarray(geometry.spliced_pos_right, np.float64) + np.asarray(geometry.spliced_neg_right, np.float64)
    espl_l = np.maximum(np.asarray(geometry.eff_spl_left, np.float64), _EPS)
    espl_r = np.maximum(np.asarray(geometry.eff_spl_right, np.float64), _EPS)
    rho_spliced = np.where(spl_l > _EPS, spl_l / espl_l, 0.0) + np.where(spl_r > _EPS, spl_r / espl_r, 0.0)
    return rho_unspl, rho_unspl + rho_spliced


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
    is the message currency (computed inline in `bp_solver.node_sweep`), mass `m=f·M_face` (`NodeDeconv`) the
    output."""

    f_pos: np.ndarray
    f_neg: np.ndarray
    f_g: np.ndarray
    var_pos: np.ndarray
    var_neg: np.ndarray
    var_gdna: np.ndarray


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
    ``free_s = nrna_active_s(left) & nrna_active_s(right)``. This BLOCKS RNA at a TSS/TES (intergenic↔exon →
    neither strand continuous → a gDNA sink) and at a mixed exon↔AMBIG seam (the non-shared strand cannot
    cross). ``mrna_active_s = exon_s(left) & exon_s(right)`` is the tighter **mature**-crossing gate
    (contiguous exon on both flanks) that selects the per-node prior (`node_prior_design.md` §3).
    """
    bs = boundary_substrate
    lr = np.asarray(bs.left_region, dtype=np.int64)
    rr = np.asarray(bs.right_region, dtype=np.int64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)  # off-edge terminal flank → no bits
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    free_pos = nrp_l & nrp_r  # +strand nascent continuous across the seam
    free_neg = nrn_l & nrn_r
    mrp_l, mrn_l = mrna_active_strands(sig_l)
    mrp_r, mrn_r = mrna_active_strands(sig_r)
    mrna_pos = mrp_l & mrp_r  # +strand mature (contiguous exon) crosses the seam
    mrna_neg = mrn_l & mrn_r

    left, right = bs.left, bs.right
    u_pos = np.maximum(left.n_unspliced_pos, right.n_unspliced_pos).astype(np.float64)
    u_neg = np.maximum(left.n_unspliced_neg, right.n_unspliced_neg).astype(np.float64)
    mass_unspl = np.asarray(left.mass_unspliced, float) + np.asarray(right.mass_unspliced, float)
    mass_spliced = np.asarray(left.mass_spliced, float) + np.asarray(right.mass_spliced, float)
    return free_pos, free_neg, mrna_pos, mrna_neg, u_pos, u_neg, mass_unspl, mass_spliced


def _region_strand_stats(substrate, region_arrays):
    """Per-region strand-solve sufficient statistics + node class (the region twin of
    :func:`_boundary_strand_stats`). A region's nascent-active axes are its own ±transcript bits (``free_s``
    from the strand class); its ``mrna_active_s`` axes are its own ±exon bits (the prior selector)."""
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    free_pos = (ts == TS_POS) | (ts == TS_AMBIG)
    free_neg = (ts == TS_NEG) | (ts == TS_AMBIG)
    mrna_pos, mrna_neg = mrna_active_strands(np.asarray(region_arrays.signature))  # own ±exon bits
    mass_u = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_s = np.asarray(c.mass_spliced, dtype=np.float64)
    return free_pos, free_neg, mrna_pos, mrna_neg, u_pos, u_neg, mass_u, mass_s


@dataclass(frozen=True, slots=True)
class NodeStatics:
    """Per-chain-node STATIC strand-solve sufficient statistics + node class (length ``n_nodes``) — the
    region- and boundary-keyed stats gathered onto the chain once. The sweep mutates only the dynamic
    :class:`NodeBelief`; these never change. ``free_pos``/``free_neg`` are the nascent-RNA-active axes (a
    region's ±transcript bits / a boundary's ±continuity — the RNA-crossing gate); ``mrna_active_pos``/
    ``mrna_active_neg`` are the mature-RNA-active axes (a region's ±exon bits / a boundary's ±contiguous
    exon) that select the per-node solver prior (`node_prior_design.md` §3). All ``float64`` except the four
    bool masks."""

    n_nodes: int
    u_pos: np.ndarray
    u_neg: np.ndarray
    free_pos: np.ndarray  # bool — nascent-RNA-active (transcript continuity); the RNA-crossing gate
    free_neg: np.ndarray  # bool
    mrna_active_pos: (
        np.ndarray
    )  # bool — mature-RNA-active (contiguous exon); selects the node prior
    mrna_active_neg: np.ndarray  # bool
    mass_unspliced: np.ndarray
    mass_spliced: np.ndarray


def build_node_statics(
    chain: NodeChain, substrate, boundary_substrate, region_arrays
) -> NodeStatics:
    """Gather the per-region (contained) and per-boundary (continuity-gated, max-of-sides) strand-solve
    statistics onto the unified chain."""
    r_fp, r_fn, r_mrp, r_mrn, r_up, r_un, r_mu, r_ms = _region_strand_stats(
        substrate, region_arrays
    )
    b_fp, b_fn, b_mrp, b_mrn, b_up, b_un, b_mu, b_ms = _boundary_strand_stats(
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
        # No per-node spliced FLOOR: spliced (mature) handling is OWNED by the message system (the B→exon
        # MEASUREMENT source + the exon→B absorption in _scan). A node-local floor would double-count it AND
        # inflate a boundary's UNSPLICED f_pos with mature → phantom nascent into introns (matrix-confirmed:
        # removing it is ≥ keeping it in every κ × capture × ±gDNA regime).
        free_pos=free_pos,
        free_neg=free_neg,
        mrna_active_pos=pick(r_mrp, b_mrp, False),
        mrna_active_neg=pick(r_mrn, b_mrn, False),
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
    n_grid_ss: int | None = None,
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
        st.free_pos,
        st.free_neg,
        st.mass_unspliced,
        st.mass_spliced,
        kappa=float(rna_sense_frac),
        od_g=gdna_strand_overdispersion,
        od_r=rna_strand_overdispersion,
        n_grid=n_grid,
        n_grid_ss=n_grid_ss,
        L=logodds_window,
    )
    f_pos, f_neg, f_g, var_p, var_n, var_g = _type_belief(
        st.free_pos, st.free_neg, deconv, st.mass_unspliced
    )
    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_p, var_neg=var_n, var_gdna=var_g
    )


# ---------------------------------------------------------------------------
# Region-type helper (the sweep itself lives in bp_solver.node_sweep).
# ---------------------------------------------------------------------------


def _node_region_type(chain, region_arrays):
    """Per-chain-node region type for REGION nodes (0=intergenic, 1=intron, 2=exon; exon>intron), −1 on
    boundary nodes; plus the per-REGION type array. Single source of truth: :func:`signature.coarse_type_array`."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)  # per REGION
    ri_ = np.clip(idx, 0, rtype.shape[0] - 1)
    return np.where(kind == REGION, rtype[ri_], -1), rtype
