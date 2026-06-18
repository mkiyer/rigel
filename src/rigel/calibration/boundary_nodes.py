"""Boundary nodes of the bipartite R↔B↔R graph (CALIBRATION_PLAN_v6; bipartite_graph_implementation_plan.md).

A boundary node is one genomic seam, solved by the SAME per-node simplex core as a region
(``simplex_sweep._solve_nodes``) with boundary-specific sufficient statistics built from ``BoundarySubstrate``:

* **count `N` (strand BB power) = ONE side** — `max(flux_left, flux_right)` per strand; for unspliced the two
  sides are exactly equal (each contiguous crossing credits both once), so `max` = the distinct crossing count
  and `max` also guards terminals (the off-edge side is 0). NOT the sum (which doubles). §4.1.
* **mass (the pie split base) = `left + right`** — the conserved, partitioned fragment mass. §4.1.
* **the one-sided motif-stranded SPLICED floor** — per strand, the sense spliced on THAT strand's exon flank
  (`_junction_side`); a boundary is one junction = one strand, so its spliced is on exactly one side, never
  summed across sides (§4.16). It enters as the existing sided-spliced lower-bound term in `_local_loglik`.
* **node-class (`allow_pos`/`allow_neg`/`strand_obs`)** from the two flanking region strand classes — flank-OR
  for the allow mask (a {POS,NEG} seam admits both, AMBIG-like); the consistent-sense rule (the same as
  `strand_deconv._side_strand_orientation`) for `strand_obs`.

Solved in ISOLATION in Phase 2 (no chain imputation, not yet wired into calibrate — that is Phase 3).
"""

from __future__ import annotations

import numpy as np

from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
)
from .simplex_sweep import _solve_nodes

__all__ = ["boundary_node_class", "deconv_boundaries_sweep"]

_EPS = 1.0e-9


def boundary_node_class(ts_left, ts_right):
    """Per-boundary ``(allow_pos, allow_neg, strand_obs)`` from the two flanking region strand classes.

    ``allow_s`` = either flank carries strand s (flank-OR; a ``{POS,NEG}`` seam admits BOTH strands, AMBIG-like).
    ``strand_obs`` = the two flanks define a single consistent transcript sense (``{POS,POS}``/``{POS,NONE}``;
    ``TS_NONE`` is a wildcard) — the same rule as ``strand_deconv._side_strand_orientation``. A ``-1`` terminal
    flank is passed as ``TS_NONE`` by the caller (a wildcard; the present flank determines the sense).
    """
    tl = np.asarray(ts_left)
    tr = np.asarray(ts_right)
    has_pos = ((tl == TS_POS) | (tl == TS_AMBIG)) | ((tr == TS_POS) | (tr == TS_AMBIG))
    has_neg = ((tl == TS_NEG) | (tl == TS_AMBIG)) | ((tr == TS_NEG) | (tr == TS_AMBIG))
    either_ambig = (tl == TS_AMBIG) | (tr == TS_AMBIG)
    cons_pos = (
        ~either_ambig
        & ((tl == TS_POS) | (tl == TS_NONE)) & ((tr == TS_POS) | (tr == TS_NONE))
        & ((tl == TS_POS) | (tr == TS_POS))
    )
    cons_neg = (
        ~either_ambig
        & ((tl == TS_NEG) | (tl == TS_NONE)) & ((tr == TS_NEG) | (tr == TS_NONE))
        & ((tl == TS_NEG) | (tr == TS_NEG))
    )
    return has_pos, has_neg, (cons_pos | cons_neg)


def _spliced_floor(sig_l, sig_r, left_sense, right_sense, exon_bit, intron_bit):
    """The one-sided spliced floor count for ONE strand (§4.16): the sense spliced on that strand's exon flank
    (a clean exon↔intron transition, the vectorized `splice_junction._junction_side` rule), else 0. Reads ONE
    side per strand — never sums the two sides (a boundary is one junction = one strand)."""
    le = (sig_l & exon_bit) != 0
    li = (sig_l & intron_bit) != 0
    re = (sig_r & exon_bit) != 0
    ri = (sig_r & intron_bit) != 0
    side_l = le & ~li & ri & ~re  # exon on the left  → spliced on the left side
    side_r = li & ~le & re & ~ri  # exon on the right → spliced on the right side
    return np.where(side_l, left_sense, np.where(side_r, right_sense, 0.0))


def deconv_boundaries_sweep(
    boundary_substrate,
    region_arrays,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    n_grid=20,
    rho_global=0.0,
    boundary_support=None,
    global_tau=None,
    gdna_imp_frac=None,
    gdna_imp_precision=None,
):
    """Solve the boundary nodes by the shared per-node core (`simplex_sweep._solve_nodes`).

    ``boundary_support`` = the per-boundary averaged side density length ``S_b`` (the global-prior density
    divisor); ``None`` ⇒ no global prior (toy tests). The boundary count, mass, spliced floor, and node-class
    are built per the module docstring. Returns a per-boundary ``NodeDeconv``.
    """
    bs = boundary_substrate
    lr = np.asarray(bs.left_region)
    rr = np.asarray(bs.right_region)
    sc = np.asarray(region_arrays.strand_class)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    # Flank strand classes / signatures; a -1 terminal flank → TS_NONE / signature 0 (wildcard; ~no mass).
    ts_l = np.where(lr >= 0, sc[np.clip(lr, 0, None)], TS_NONE)
    ts_r = np.where(rr >= 0, sc[np.clip(rr, 0, None)], TS_NONE)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)

    left, right = bs.left, bs.right
    # COUNT (BB power) = ONE side (max; sides equal for unspliced, max guards terminals — §4.1).
    u_pos = np.maximum(left.n_unspliced_pos, right.n_unspliced_pos).astype(np.float64)
    u_neg = np.maximum(left.n_unspliced_neg, right.n_unspliced_neg).astype(np.float64)
    # MASS (the pie split base) = both sides summed (conserved — §4.1).
    mass_unspl = (
        np.asarray(left.mass_unspliced, dtype=np.float64)
        + np.asarray(right.mass_unspliced, dtype=np.float64)
    )
    mass_spliced = (
        np.asarray(left.mass_spliced, dtype=np.float64)
        + np.asarray(right.mass_spliced, dtype=np.float64)
    )
    # SPLICED floor (§4.16): per strand, the sense spliced on THAT strand's exon flank; never summed.
    ls = left.n_spliced_sense.astype(np.float64)
    rs = right.n_spliced_sense.astype(np.float64)
    spliced_pos = _spliced_floor(sig_l, sig_r, ls, rs, BIT_EXON_POS, BIT_INTRON_POS)
    spliced_neg = _spliced_floor(sig_l, sig_r, ls, rs, BIT_EXON_NEG, BIT_INTRON_NEG)
    allow_pos, allow_neg, strand_obs = boundary_node_class(ts_l, ts_r)

    if boundary_support is not None:
        s_b = np.maximum(np.asarray(boundary_support, dtype=np.float64), _EPS)
        mu_global = np.clip(
            np.where(mass_unspl > 0.0, float(rho_global) * s_b / np.maximum(mass_unspl, _EPS), 0.0),
            0.0, 1.0,
        )
        gtau = 1.0 if global_tau is None else np.asarray(global_tau, dtype=np.float64)
    else:
        mu_global = None
        gtau = 0.0

    return _solve_nodes(
        u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs,
        mass_unspl, mass_spliced,
        kappa=float(rna_sense_frac), od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
        n_grid=n_grid, mu_global=mu_global, global_tau=gtau,
        gdna_imp_frac=gdna_imp_frac, gdna_imp_precision=gdna_imp_precision,
    )
