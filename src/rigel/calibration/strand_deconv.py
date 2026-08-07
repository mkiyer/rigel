"""The contiguous-edge seeds for the gDNA strand-overdispersion fit.  LAYER 4 — strand.

After the bipartite belief-propagation rebuild (:mod:`rigel.calibration.sweep`), the per-node gDNA/RNA
deconvolution lives in the chain sweep. What is left here is edge strand geometry:

* :func:`edge_seeds` — the exon–intron / exon–intergenic ``(sense, total, gDNA weight)`` seeds for the
  gDNA strand-overdispersion fit (:mod:`rigel.calibration.gdna_strand`), complementing the contained-node
  seeds (needed under hybrid capture, which depletes off-target intergenic / intronic gDNA).

⭐ **ONE SEED PER LINE, NOT TWO PER BOUNDARY (S5.f).** ``boundary_side_seeds`` emitted a seed for each
of a boundary's two faces — region ``r``'s right side and region ``r+1``'s left side — because the old
accumulator split one crossing fragment's mass across the two flanks. They are the same physical
crossing: pooling both into one method-of-moments estimator doubles its apparent sample size and pairs
every observation with a perfectly correlated twin, which is a dispersion estimate reading its own
duplication. A contiguous edge is a 0-bp line with one count, so there is one seed and the whole
``_SideQuantities`` / ``_compute_side`` / ``_left_right_neighbors`` layer dissolves with the faces.

⭐ **And the seed WEIGHT is now exactly 1, provably.** It was
``clip(density · eff / mass)`` where a count-observable side read its own crossing density
``density = mass / eff`` — algebraically 1 whenever the seed mask admits it, and the borrowed-density
branch could never reach a seed. With ``count`` and ``mass`` the same number that identity is exact
rather than approximate, so the ratio, the effective length it divided by, and the
``boundary_side_eff_len`` argument that carried it all go. A seed's weight is "this line is gDNA by
signature", and that is a 1.
"""

from __future__ import annotations


import numpy as np

# ⭐ :class:`NodeDeconv` MOVED TO `node_chain` (layer 0) on 2026-08-07: three layers were reaching UP for
# it, which `test_layering.py` now forbids. It is re-exported here only so existing importers keep
# working; new code should take it from `node_chain`, where the tool's central datum belongs.
from .node_chain import NodeDeconv  # noqa: F401
from .region_arrays import edge_node_indices
from .signature import TS_NEG, TS_NONE, TS_POS


def edge_strand_orientation(ts_lo: np.ndarray, ts_hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(orient_neg, strand_observable)`` per contiguous edge, from its two flanks' strand classes.

    A line is **strand-observable** iff its two flanks define a single consistent transcript sense:
    ``{POS,POS}`` / ``{NEG,NEG}``, or a gene-edge ``{POS,NONE}`` / ``{NEG,NONE}``. An **intergenic
    (``TS_NONE``) flank is a strand WILDCARD** — it carries no transcript, so it is compatible with
    either strand and the line is oriented by its gene flank. Only a genuine conflict — opposite
    strands ``{POS,NEG}``, or a ``TS_AMBIG`` flank (overlapping ± transcripts) — leaves the sense
    undefined, and such a line cannot seed the fit at all.

    ⭐ **The rule is SYMMETRIC in the two flanks**, which is precisely why the two faces of the old
    boundary always agreed about orientation and differed only in which face's count they read — the
    duplication `edge_seeds` removes.

    ``orient_neg`` selects the NEG genome-strand column as "sense"; where it is ``False`` (and the
    line is observable) the POS column is sense.

    ⚠ **There is no explicit AMBIG guard, and there must not be one.** The predecessor carried an
    ``~either_ambig`` term on both clauses. ``TS_AMBIG`` is a fourth distinct value, so
    ``(ts == TS_POS) | (ts == TS_NONE)`` is already ``False`` on it and the term was **dead** — a
    redundant clause reading as load-bearing, which is worse than absent because it invites the belief
    that the rule lives there. Perturbing it away changed no test, which is how it was found; the rule
    is now pinned by ``test_an_AMBIG_flank_cannot_seed`` instead.
    """
    lo_pos_or_none = (ts_lo == TS_POS) | (ts_lo == TS_NONE)
    hi_pos_or_none = (ts_hi == TS_POS) | (ts_hi == TS_NONE)
    lo_neg_or_none = (ts_lo == TS_NEG) | (ts_lo == TS_NONE)
    hi_neg_or_none = (ts_hi == TS_NEG) | (ts_hi == TS_NONE)
    cons_pos = lo_pos_or_none & hi_pos_or_none & ((ts_lo == TS_POS) | (ts_hi == TS_POS))
    cons_neg = lo_neg_or_none & hi_neg_or_none & ((ts_lo == TS_NEG) | (ts_hi == TS_NEG))
    return cons_neg, cons_pos | cons_neg


def edge_seeds(substrate, region_arrays, node_density) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(sense, total, gdna_weight)`` — ONE seed per count- and strand-observable contiguous edge.

    The exon–intron / exon–intergenic line seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the contained-node seeds (needed under hybrid capture, which
    depletes off-target intergenic / intronic gDNA).

    ⚠ **``gdna_weight`` is identically 1 on every seed**, and that is a derivation rather than a
    simplification — see the module docstring. It is returned as an array because the pooled estimator
    takes a per-seed weight and the contained-node seeds genuinely vary.
    """
    ts = np.asarray(region_arrays.strand_class)
    lo, hi = edge_node_indices(np.asarray(region_arrays.ref_id))
    count = np.asarray(substrate.edge_unspliced.count, dtype=np.float64)
    pos, neg = count[:, 0], count[:, 1]
    total = pos + neg

    orient_neg, strand_observable = edge_strand_orientation(ts[lo], ts[hi])
    count_observable = np.asarray(node_density.edge_count_observable, dtype=bool)
    seed = count_observable & strand_observable & (total > 0.0)

    sense = np.where(orient_neg, neg, pos)[seed]
    return sense, total[seed], np.ones(sense.shape[0], dtype=np.float64)


__all__ = ["edge_seeds", "edge_strand_orientation"]
