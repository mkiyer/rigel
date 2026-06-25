"""Per-node deconvolution result type + the boundary-side seeds for the gDNA strand-overdispersion fit.

After the bipartite belief-propagation rebuild (:mod:`rigel.calibration.bp_solver`), the per-node gDNA/RNA
deconvolution lives in the chain sweep. This module retains only the two pieces that still feed it:

* :class:`NodeDeconv` — the per-node deconvolution result (a region or a boundary side), the schema
  ``bp_solver``'s ``chain_region_deconv`` / ``chain_boundary_side_deconv`` and ``simplex_sweep._solve_nodes``
  return, and that ``priors`` / ``derive`` consume.
* :func:`boundary_side_seeds` — the exon–intron / exon–intergenic boundary-side ``(sense, total, gDNA
  weight)`` seeds for the gDNA strand-overdispersion fit (:mod:`rigel.calibration.gdna_strand`), complementing
  the contained-region seeds (needed under hybrid capture, which depletes off-target intergenic / intronic
  gDNA). The weight is the count-clue gDNA fraction (≈ 1 for a count-observable side — gDNA by signature).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .run_fill import same_ref_left_right
from .signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS


@dataclass(frozen=True, slots=True)
class NodeDeconv:
    """Per-node deconvolution result (regions or boundary sides; kept separate)."""

    gdna_mass: np.ndarray  # float64[K] — the consumed output (calibrate/derive read ONLY the *_mass fields)
    rna_mass: np.ndarray  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac: np.ndarray  # float64[K] — the node's gDNA composition (face-invariant; mass = frac·M_face)
    # per-strand RNA fractions of the UNSPLICED mass (posterior means; f_pos+f_neg+gdna_frac = 1), populated
    # by the simplex sweep for the per-strand RNA imputation (the bipartite R↔B↔R chain).
    rna_pos_frac: "np.ndarray | None" = None  # float64[K] — f_pos
    rna_neg_frac: "np.ndarray | None" = None  # float64[K] — f_neg
    # per-component posterior VARIANCES of the fractions (moment-matched from the lattice). The precision
    # state for the honest message send (`precision_state_design.md` §1: Var_own = (M/E)²·Var(f_c)); set by
    # the per-node solve, consumed when a node emits a message. None on the chain region/boundary projections
    # (precision is a chain-node property, not needed by the downstream EM prior).
    gdna_frac_var: "np.ndarray | None" = None  # float64[K] — Var(f_g)
    rna_pos_frac_var: "np.ndarray | None" = None  # float64[K] — Var(f_pos)
    rna_neg_frac_var: "np.ndarray | None" = None  # float64[K] — Var(f_neg)


@dataclass(frozen=True, slots=True)
class _SideQuantities:
    """Per-region boundary-side quantities shared by the seed builder (and the retired deconvolution)."""

    sense: np.ndarray
    antisense: np.ndarray
    n_side: np.ndarray  # pos + neg (count evidence)
    mass: np.ndarray
    mass_spliced: np.ndarray
    count_gdna_frac: np.ndarray  # count module's gDNA fraction for this side (raw density ratio)
    count_gdna_frac_var: np.ndarray  # Poisson floor σ_g² = g²/n_side (cap g(1−g))
    strand_observable: np.ndarray
    count_observable: np.ndarray


def _left_right_neighbors(ts, ref_id, boundary_count_observable):
    """Per-region neighbour strand class + same-ref + boundary count-observability for both sides.

    Region ``r``'s LEFT side is the right side of boundary ``(r−1, r)`` (neighbour ``r−1``); its
    RIGHT side is the left side of boundary ``(r, r+1)`` (neighbour ``r+1``). ``boundary_count_observable[r]``
    describes boundary ``(r, r+1)``. Returns
    ``(left_same, ts_prev, left_observable, right_same, ts_next, right_observable)``.
    """
    r = ts.shape[0]
    left_same, right_same = same_ref_left_right(ref_id)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_observable = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        ts_prev[1:] = ts[:-1]
        left_observable[1:] = boundary_count_observable[:-1]
        ts_next[:-1] = ts[1:]
    return left_same, ts_prev, left_observable, right_same, ts_next, boundary_count_observable


def _side_strand_orientation(view, same, ts_self, ts_other):
    """Per-side sense/antisense split + strand-observability (the ``TS_NONE``-wildcard rule).

    A side is **strand-observable** iff its two flanks define a single consistent transcript sense:
    ``{POS,POS}``/``{NEG,NEG}``, or a gene-edge ``{POS,NONE}``/``{NEG,NONE}`` (intergenic ``TS_NONE`` is a
    wildcard, oriented by the gene side); an opposite-strand ``{POS,NEG}`` or a ``TS_AMBIG`` flank leaves
    the sense undefined. Returns ``(sense, antisense, n_side, strand_observable)`` (all length R).
    """
    pos = view.n_unspliced_pos.astype(np.float64)
    neg = view.n_unspliced_neg.astype(np.float64)
    either_ambig = (ts_self == TS_AMBIG) | (ts_other == TS_AMBIG)
    self_pos_or_none = (ts_self == TS_POS) | (ts_self == TS_NONE)
    other_pos_or_none = (ts_other == TS_POS) | (ts_other == TS_NONE)
    self_neg_or_none = (ts_self == TS_NEG) | (ts_self == TS_NONE)
    other_neg_or_none = (ts_other == TS_NEG) | (ts_other == TS_NONE)
    cons_pos = (
        same & ~either_ambig & self_pos_or_none & other_pos_or_none
        & ((ts_self == TS_POS) | (ts_other == TS_POS))
    )
    cons_neg = (
        same & ~either_ambig & self_neg_or_none & other_neg_or_none
        & ((ts_self == TS_NEG) | (ts_other == TS_NEG))
    )
    strand_observable = cons_pos | cons_neg
    sense = np.where(cons_neg, neg, pos)
    n_side = pos + neg
    return sense, n_side - sense, n_side, strand_observable


def _compute_side(
    view, same, ts_self, ts_other, side_count_observable, eff, region_density
) -> _SideQuantities:
    """Per-side sense split, count-prior density, and strand-observability.

    A side is **strand-observable** iff its boundary's two regions define a single, consistent
    transcript sense (so 'sense' is defined). An **intergenic (TS_NONE) region is a strand wildcard**
    — it carries no transcript, so it is compatible with either strand: a gene-edge boundary
    (stranded exon ↔ intergenic) is therefore oriented by the gene side, ``{POS,NONE}→POS``,
    ``{NEG,NONE}→NEG`` (same as ``{POS,POS}`` / ``{NEG,NEG}``). Only a genuine conflict — opposite
    strands ``{POS,NEG}`` or a ``TS_AMBIG`` flank (overlapping +/− transcripts) — leaves the sense
    undefined. The **count** fraction is the raw crossing density ratio: a count-observable side (no
    shared exon) reads ``count_gdna_frac → 1`` from its own crossing mass; otherwise it borrows the swept
    region density. (No strand cleaning — the count module is raw.)

    The **count variance** is the Poisson floor (a boundary side is a *direct* measurement, not an imputed
    region): the crossing count ``n_side`` is Poisson, so the fraction's relative variance is ``1/n_side``,
    giving ``σ_g² = g_count²·(1/n_side)`` capped at the Bernoulli maximum ``g_count(1−g_count)``.
    """
    mass = view.mass_unspliced
    sense, antisense, n_side, strand_observable = _side_strand_orientation(
        view, same, ts_self, ts_other
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        # Raw crossing density (no strand cleaning); count-observable sides use their own crossing
        # density, others borrow the swept region density.
        own = np.where((mass > 0.0) & (eff > 0.0), mass / np.maximum(eff, 1e-12), 0.0)
        density = np.where(side_count_observable, own, region_density)
        count_gdna_frac = np.clip(
            np.where(mass > 0.0, density * eff / np.maximum(mass, 1e-12), 0.0), 0.0, 1.0
        )
        # Poisson floor: crossing count n_side ~ Poisson ⇒ relative variance 1/n_side;
        # σ_g² = g²·(1/n_side), capped at the Bernoulli maximum g(1−g) (sides with no crossings → 0).
        v_rel = np.where(n_side > 0.0, 1.0 / np.maximum(n_side, 1e-12), 0.0)
        count_gdna_frac_var = np.minimum(
            count_gdna_frac**2 * v_rel, count_gdna_frac * (1.0 - count_gdna_frac)
        )
    return _SideQuantities(
        sense=sense,
        antisense=antisense,
        n_side=n_side,
        mass=mass,
        mass_spliced=view.mass_spliced,
        count_gdna_frac=count_gdna_frac,
        count_gdna_frac_var=count_gdna_frac_var,
        strand_observable=strand_observable,
        count_observable=np.asarray(side_count_observable, dtype=bool),
    )


def boundary_side_seeds(substrate, region_arrays, node_density, boundary_side_eff_len):
    """``(sense, total, gdna_weight)`` seed arrays from count- & strand-observable boundary sides.

    The exon–intron / exon–intergenic seam seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the contained-region seeds (needed under hybrid capture,
    which depletes off-target intergenic/intronic gDNA). The weight is the count-clue gDNA fraction
    ``count_gdna_frac`` (≈ 1 for a count-observable side, gDNA by signature) — the same quantity the
    contained-region seeds use (``node_density.count_gdna_frac``).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )
    sides = (
        _compute_side(substrate.left, l_same, ts, ts_prev, l_obs, eff, region_density),
        _compute_side(substrate.right, r_same, ts, ts_next, r_obs, eff, region_density),
    )
    senses, totals, weights = [], [], []
    for sq in sides:
        seed = sq.count_observable & sq.strand_observable & (sq.n_side > 0.0)
        senses.append(sq.sense[seed])
        totals.append(sq.n_side[seed])
        weights.append(sq.count_gdna_frac[seed])
    return np.concatenate(senses), np.concatenate(totals), np.concatenate(weights)


__all__ = ["NodeDeconv", "boundary_side_seeds"]
