"""The region/boundary chain the sweep traverses, and its cut into independent locus blocks; layer 0.

The chain interleaves REGION and BOUNDARY slots in genomic order. A reference with ``k`` regions owns exactly
``k - 1`` interior boundaries::

    N0  E0  N1  E1  ...  E(k-2)  N(k-1)          2k - 1 slots

It starts and ends with a REGION, so every boundary has a region on both sides: boundary ``i`` lies between regions
``i`` and ``i + 1``, an implicit relation computed only here. A slot is addressed by ``(kind, obj_idx)``, where
``obj_idx`` indexes the region axis or the contiguous-boundary axis, so every per-object statistic stays in its own
payload-shaped array. Splice junctions are not slots: each closes an undirected loop, so a sj is a factor on its
endpoints, never a message channel.

Gate: ``tests/calibration/test_region_geometry.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

__all__ = [
    "BOUNDARY",
    "REGION",
    "LocusBlock",
    "RegionChain",
    "RegionDeconv",
    "build_region_chain",
    "locus_blocks",
]

REGION = 0
BOUNDARY = 1


@dataclass(frozen=True, slots=True)
class RegionChain:
    """The genomic-ordered region-and-boundary chain and its adjacency, all arrays of length ``n_slots``.

    Slot ids are assigned in genomic visiting order, so ``order`` would be ``arange`` and is not
    stored. ``left`` / ``right`` give each slot its single adjacent slot of the other kind, and
    ``-1`` at a reference terminal — a propagation sink, which is always a REGION.
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


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# One slot's deconvolution result. It lives in the bottom layer because layers 3 and 6 both use it.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class RegionDeconv:
    """Per-slot deconvolution result, with two uses and hence optional fields.

    * The per-slot solve (``simplex_logodds._solve_regions_logodds_all``) returns the fractions and ``gdna_frac_var``
      but no mass, since mass is per face and the solve is face-invariant.
    * The chain projection (``sweep.chain_region_deconv`` / ``chain_boundary_deconv``) returns the masses and
      fractions ``CalibrationResult`` publishes, but no variance.
    """

    gdna_frac: np.ndarray  # float64[K] — the gDNA composition (face-invariant; mass = frac·M_face)
    # per-strand RNA fractions of the unspliced mass, so f_pos + f_neg + gdna_frac = 1
    rna_pos_frac: "np.ndarray | None" = None  # float64[K] — f_pos
    rna_neg_frac: "np.ndarray | None" = None  # float64[K] — f_neg
    # the projection's masses; None on the per-slot solve. ``derive`` and the capture efficiencies read the gDNA one.
    gdna_mass: "np.ndarray | None" = None  # float64[K]
    rna_mass: "np.ndarray | None" = None  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    # Var(log f_g), not Var(f_g), so not bounded by ¼; set by the per-slot solve, read by the landscape's training
    # rule, None on the projections
    gdna_frac_var: "np.ndarray | None" = None  # float64[K]


def build_region_chain(
    ref_region_offsets: np.ndarray, ref_boundary_offsets: np.ndarray
) -> RegionChain:
    """Build the chain from the payload's two per-reference CSR offset arrays.

    Reference ``f`` owns regions ``[rno[f], rno[f+1])`` and contiguous boundaries ``[reo[f], reo[f+1])``, with
    ``boundaries == max(regions - 1, 0)``; a reference with no regions contributes nothing. Both arrays come from one
    accumulator payload, so a mismatch is a payload inconsistency, not a stale index.
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
            f"boundaries (expected {int(expected[bad])}). There are no terminal boundary slots: a boundary is the "
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
        # link consecutive slots within this reference; the two terminals keep -1, so no message crosses references
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


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# The locus: where the chain breaks. Which slots end message passing is the caller's predicate.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


class LocusBlock(NamedTuple):
    """One unit of the sweep's work: the slots ``[start, stop)`` it owns and the slots ``[start, end)`` it reads.

    ``end`` is ``stop + 1`` when the slot at ``stop`` is a terminal linked to the block's last slot: it is read as a
    message source but owned by the next block, which starts on it. A block starts on a terminal or a reference
    start, never inside a locus.
    """

    start: int
    stop: int
    end: int


def locus_blocks(
    chain: RegionChain, terminal: np.ndarray, block_slots: int | None = None
) -> list[LocusBlock]:
    """Cut the chain into blocks (:class:`LocusBlock`) at its message terminals, merging consecutive loci up to
    ``block_slots`` slots per block.

    A locus is the run of slots between two terminals. ``terminal`` flags the slots where message passing ends in
    both directions (in the sweep, a REGION that admits no RNA strand); a reference start is also a cut. The
    two-phase solve is not an iteration, so a node's answer depends only on messages that can reach it, and none
    crosses a terminal: the sweep over the blocks, in any order, is exactly the sweep over the chain.

    ``block_slots`` trades a block's working set against per-block overhead and moves no number; ``None`` makes one
    block. A block never splits a locus, so a longer locus is a block of its own length. Every slot is owned by
    exactly one block, in chain order: consecutive blocks abut and the last stops at ``n_slots``.
    """
    n = int(chain.n_slots)
    left = np.asarray(chain.left, dtype=np.int64)
    term = np.asarray(terminal, dtype=bool)
    if term.shape != (n,):
        raise ValueError(f"terminal has shape {term.shape}; expected ({n},), one flag per slot")
    if n == 0:
        return []
    cut = term | (left < 0)
    cut[0] = True
    cuts = np.flatnonzero(cut).tolist()
    cuts.append(n)
    if block_slots is not None and int(block_slots) < 1:
        raise ValueError(f"block_slots must be >= 1 or None; got {block_slots}")
    limit = n if block_slots is None else int(block_slots)

    def block(start: int, stop: int) -> LocusBlock:
        # the slot at ``stop`` is read as a source iff it is linked to the block's last slot
        reads_next = stop < n and left[stop] >= 0
        return LocusBlock(int(start), int(stop), int(stop + 1 if reads_next else stop))

    blocks: list[LocusBlock] = []
    start = cuts[0]
    for j in range(1, len(cuts)):
        nxt = cuts[j]
        # take the next locus whenever the block is still empty, else only while it fits
        if nxt - start > limit and cuts[j - 1] > start:
            blocks.append(block(start, cuts[j - 1]))
            start = cuts[j - 1]
        if j == len(cuts) - 1:
            blocks.append(block(start, nxt))
    return blocks
