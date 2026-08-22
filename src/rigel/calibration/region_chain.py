"""rigel.calibration.region_chain — the region<->boundary chain the belief-propagation sweep traverses.

       Gate: ``tests/calibration/test_region_chain.py``

The calibration graph is a **linear bipartite chain of REGION and BOUNDARY slots, interleaved in genomic
order**. A reference with ``k`` regions owns exactly ``k − 1`` interior boundaries, so its slot sequence is::

    N0  E0  N1  E1  ...  E(k-2)  N(k-1)          2k − 1 slots

⭐ **It starts and ends with a REGION, and there are no terminal slots.** That is the whole shape change
from the predecessor, which ran ``B R B R … R B`` with ``k + 1`` boundary slots per reference — the two
outermost carrying no data and existing only so every region had an object on each side. A contiguous
boundary is the boundary BETWEEN two adjacent regions; there is no such boundary before the first or after the last, so
the object does not exist rather than existing empty. **An boundary therefore always has a region on both
sides**, an invariant the old shape could not state.

Boundary endpoints are **implicit**: boundary ``i`` lies between region ``i`` and region ``i + 1``. Nothing stores
them, and this module is the one place that arithmetic lives.

A slot is addressed by ``(kind, obj_idx)``: ``kind`` is :data:`REGION` or :data:`BOUNDARY`, and ``obj_idx``
indexes the region axis or the contiguous-boundary axis respectively. That keeps every per-object statistic in
its own payload-shaped array — the chain only sequences and links them.

⚠ **SpliceJunction boundaries are NOT chain slots.** The graph is a DAG but not a polytree: every sj boundary
closes an undirected loop, so a sj must be a FACTOR on its endpoint regions and never a message
channel (— never break a cycle by dropping a sj boundary,
that re-isolates the exon the boundary exists for).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["BOUNDARY", "REGION", "RegionChain", "RegionDeconv", "build_region_chain"]

REGION = 0
BOUNDARY = 1


@dataclass(frozen=True, slots=True)
class RegionChain:
    """The genomic-ordered region∪boundary chain and its adjacency. All arrays have length ``n_slots``.

    Slot ids are assigned in genomic visiting order, so ``order`` would be ``arange`` and is not stored.
    ``left``/``right`` give each slot its single adjacent slot of the OTHER kind, ``-1`` at a reference
    terminal — which is a propagation sink, and is now always a REGION.
    """

    kind: np.ndarray  # int8[n_slots] — REGION or BOUNDARY
    obj_idx: (
        np.ndarray
    )  # int64[n_slots] — index into the region axis, or the contiguous-boundary axis
    left: np.ndarray  # int64[n_slots] — adjacent slot id, -1 at a reference start
    right: np.ndarray  # int64[n_slots] — adjacent slot id, -1 at a reference end
    n_regions_total: int
    n_boundaries_total: int

    @property
    def n_slots(self) -> int:
        return int(self.kind.shape[0])

    @property
    def is_region(self) -> np.ndarray:
        return self.kind == REGION

    @property
    def is_boundary(self) -> np.ndarray:
        return self.kind == BOUNDARY


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# ⭐⭐ ONE SLOT'S DECONVOLUTION RESULT — vocabulary, and it lives here because THREE LAYERS need it.
# It was defined in the STRAND family (layer 4) and imported by `region_geometry` and `simplex_logodds`
# (layer 3) and `sweep` (layer 6), so three layers reached UPWARD for a type. `module_census.py` is what
# made that visible, and the repair is the one a layering violation always asks for: the TYPE belongs at
# the bottom, not the code that happened to define it first. ⚠ It is not a strand concept — the pie
# `(f_pos, f_neg, f_g)` is the tool's central datum, and a slot is what carries it.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class RegionDeconv:
    """Per-region deconvolution result. TWO disjoint uses, hence the optional halves:

    * the per-region SOLVE (`simplex_logodds._solve_regions_logodds_all`) returns the **composition** —
      ``*_frac`` + ``*_frac_var`` — and no mass (a region's mass is a per-FACE quantity; the solve is
      face-invariant, so a single ``*_mass`` here would be meaningless);
    * the chain PROJECTION (`sweep.chain_region_deconv` / `chain_boundary_deconv`) returns the
      **mass** the downstream `CalibrationResult` consumes, and no precision.
    """

    gdna_frac: (
        np.ndarray
    )  # float64[K] — the region's gDNA composition (face-invariant; mass = frac·M_face)
    # per-strand RNA fractions of the UNSPLICED mass (posterior means; f_pos+f_neg+gdna_frac = 1), populated
    # by the simplex sweep for the per-strand RNA imputation (the bipartite R↔B↔R chain).
    rna_pos_frac: "np.ndarray | None" = None  # float64[K] — f_pos
    rna_neg_frac: "np.ndarray | None" = None  # float64[K] — f_neg
    # per-component posterior variances in LOG-FRACTION space — `Var(log f_c)`, NOT `Var(f_c)`. They are
    # grid moments of `log f_c` over the λ lattice (`simplex_logodds._solve_regions_logodds`), because the
    # message currency is a log-density and the send precision `1/(Var(log f_c) + 1/n + σ²_transfer)` is
    # log-space throughout. ⚠ They are therefore NOT bounded by ¼ and routinely exceed it — a consumer that
    # needs the LINEAR `Var(f_c)` must convert (delta method: `Var(f_c) ≈ f_c²·Var(log f_c)`, as
    # `sweep.solve_chain` does when it builds `_var_fg` for `composition_logvar`). Set by the per-region
    # solve, consumed when a region emits a message. None on the chain region/boundary projections (precision
    # is a chain-region property, not needed by the downstream EM prior).
    # the PROJECTION's consumed output (calibrate/derive read ONLY these); None on the per-region solve.
    gdna_mass: "np.ndarray | None" = None  # float64[K]
    rna_mass: "np.ndarray | None" = None  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_g)
    rna_pos_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_pos)
    rna_neg_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_neg)


def build_region_chain(
    ref_region_offsets: np.ndarray, ref_boundary_offsets: np.ndarray
) -> RegionChain:
    """Build the chain from the payload's two per-reference CSR offset arrays.

    Reference ``f`` owns regions ``[rno[f], rno[f+1])`` and contiguous boundaries ``[reo[f], reo[f+1])``, with
    ``boundaries == max(regions − 1, 0)``. A reference with no regions contributes nothing at all, which is legal.

    ⚠ Both arrays come from ONE accumulator payload, so a mismatch between them is an accumulator /
    payload inconsistency and **not** a stale index — rebuilding will not fix it, and the error says so.
    """
    region_offsets = np.asarray(ref_region_offsets, dtype=np.int64)
    boundary_offsets = np.asarray(ref_boundary_offsets, dtype=np.int64)
    if region_offsets.shape != boundary_offsets.shape:
        raise ValueError(
            f"ref_region_offsets has shape {region_offsets.shape} and ref_boundary_offsets "
            f"{boundary_offsets.shape}; both are per-reference CSR arrays of length n_refs + 1."
        )
    n_refs = region_offsets.shape[0] - 1
    regions_per_ref = np.diff(region_offsets)
    boundaries_per_ref = np.diff(boundary_offsets)
    expected = np.maximum(regions_per_ref - 1, 0)
    if not np.array_equal(boundaries_per_ref, expected):
        bad = int(np.argmax(boundaries_per_ref != expected))
        raise ValueError(
            f"reference {bad}: the payload reports {int(boundaries_per_ref[bad])} contiguous boundaries for "
            f"{int(regions_per_ref[bad])} regions, but a reference with k regions has exactly k-1 interior "
            f"boundaries (expected {int(expected[bad])}). There are no terminal boundary slots: an boundary is the "
            f"boundary BETWEEN two adjacent regions. Both offset arrays come from ONE accumulator payload, so "
            f"this is an accumulator/payload inconsistency, not a stale index — rebuilding will not fix it."
        )

    n_slots = int(regions_per_ref.sum() + boundaries_per_ref.sum())
    kind = np.empty(n_slots, dtype=np.int8)
    obj_idx = np.empty(n_slots, dtype=np.int64)
    left = np.full(n_slots, -1, dtype=np.int64)
    right = np.full(n_slots, -1, dtype=np.int64)

    slot = 0
    for f in range(n_refs):
        k = int(regions_per_ref[f])
        if k == 0:
            continue
        region_base, boundary_base = int(region_offsets[f]), int(boundary_offsets[f])
        first_slot = slot
        # N0 E0 N1 E1 ... E(k-2) N(k-1): a region, then the boundary to its right, except after the last region
        for i in range(k):
            kind[slot] = REGION
            obj_idx[slot] = region_base + i
            slot += 1
            if i < k - 1:
                kind[slot] = BOUNDARY
                obj_idx[slot] = boundary_base + i
                slot += 1
        # link consecutive slots WITHIN this reference; the two terminals keep -1, so a sweep cannot
        # relay across a reference boundary
        ref_slots = np.arange(first_slot, slot, dtype=np.int64)
        left[ref_slots[1:]] = ref_slots[:-1]
        right[ref_slots[:-1]] = ref_slots[1:]

    return RegionChain(
        kind=kind,
        obj_idx=obj_idx,
        left=left,
        right=right,
        n_regions_total=int(region_offsets[-1]),
        n_boundaries_total=int(boundary_offsets[-1]),
    )
