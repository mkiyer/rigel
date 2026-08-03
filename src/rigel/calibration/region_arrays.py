"""rigel.calibration.region_arrays — sorted node geometry + the node↔edge mapping.

Two pieces of pure geometry the calibrator builds on:

* :class:`RegionArrays` — a per-reference-CSR view of the node table,
  sorted by ``(ref_id, start)`` so each reference's rows are contiguous and
  ascending. Carries the structural columns plus the int8 transcript-strand
  class derived from each node's signature (the strand-model input).

* The **node↔contiguous-edge index mapping** — :func:`node_right_edge` and
  :func:`edge_node_indices`.

⭐ **The ``k + 1`` boundary axis is retired (S5.f).** A reference with ``k`` nodes used to own
``k + 1`` boundary slots — the ``k − 1`` interior seams plus two data-free terminals that existed
only so every region had an object on each side. A contiguous edge is the line BETWEEN two adjacent
nodes: there is no such line before the first or after the last, so a reference owns exactly
``k − 1`` of them and **an edge always has a node on both sides**. That kills the ``-1``-terminal
branch, the two-spaces-off-by-one-per-reference arithmetic, and the pair of offset arrays the old
mapping needed — the edge axis is derivable from ``ref_id`` alone.

⚠ The derivation rests on ONE fact: edge ids are assigned per reference in genomic order, in
reference order, which is exactly the order adjacent same-reference node pairs appear in a
``(ref_id, start)``-sorted node table. :func:`~rigel.calibration.node_chain.build_node_chain` lays
out the same numbering by walking the payload's CSR offsets, and
``test_edge_numbering_matches_the_chain_built_from_the_payload_offsets`` pins the two against each
other — a second algorithm, not a second call to the first.

No tunable parameters: this module is index arithmetic only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from .signature import transcript_strand_class


__all__ = [
    "RegionArrays",
    "edge_node_indices",
    "node_right_edge",
]


@dataclass(frozen=True, slots=True)
class RegionArrays:
    """Per-reference-CSR view of the region table.

    Rows are sorted by ``(ref_id, start)`` so each ref's rows are contiguous
    and ascending. ``ref_offsets`` is the per-ref CSR boundary array; the
    slice ``[ref_offsets[r]:ref_offsets[r + 1]]`` is one reference's nodes, sorted
    by ``start`` (and by ``end`` — nodes tile a reference).
    """

    ref_id: np.ndarray  # int32, (R,)
    start: np.ndarray  # int64, (R,)
    end: np.ndarray  # int64, (R,)
    signature: np.ndarray  # uint8, (R,)
    strand_class: np.ndarray  # int8,  (R,) — TS_NONE/TS_POS/TS_NEG/TS_AMBIG
    region_size_bp: np.ndarray  # float64, (R,) — end - start
    ref_offsets: np.ndarray  # int32, (n_refs + 1,)
    n_refs: int

    @property
    def n_regions(self) -> int:
        return int(self.start.shape[0])

    @classmethod
    def from_index(cls, index) -> "RegionArrays":
        """⭐ Build the geometry for **the partition the scanner actually deposits into**.

        This and :func:`~rigel.calibration.splice_graph.build_node_partition_arrays` are the two halves
        of one contract — the calibration geometry must address the payload the scanner produced —
        so they read the same frame (``index.nodes_df``, the v8 splice graph) through one accessor.
        Passing a frame by hand is how the two drift apart, and nothing downstream detects that
        except as a shape error far from its cause.
        """
        return cls.from_frame(index.nodes_df, index.ref_name_to_id)

    @classmethod
    def from_frame(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
    ) -> "RegionArrays":
        """Build from a partition frame directly — for tests that construct a partition by hand.

        Production uses :meth:`from_index`, which cannot address a different partition than the
        scanner did."""
        if "signature" not in region_df.columns:
            raise ValueError(
                "RegionArrays.from_frame: region_df is missing the 'signature' "
                "column. Rebuild the index against the current schema."
            )
        n_refs = len(ref_name_to_id)
        n_regions = len(region_df)

        ref_id = region_df["ref_name"].map(ref_name_to_id).to_numpy()
        if np.any(pd.isna(ref_id)):
            unknown = region_df.loc[pd.isna(ref_id), "ref_name"].unique().tolist()
            raise ValueError(
                f"RegionArrays.from_frame: region_df references {sorted(unknown)} "
                f"which are not in ref_name_to_id. Rebuild the index."
            )
        ref_id = ref_id.astype(np.int32, copy=False)

        order = np.lexsort((region_df["start"].to_numpy(), ref_id))
        ref_id = ref_id[order]
        start = region_df["start"].to_numpy().astype(np.int64, copy=False)[order]
        end = region_df["end"].to_numpy().astype(np.int64, copy=False)[order]
        signature = region_df["signature"].to_numpy().astype(np.uint8, copy=False)[order]
        strand_class = transcript_strand_class(signature)

        counts = np.bincount(ref_id, minlength=n_refs).astype(np.int32, copy=False)
        ref_offsets = np.empty(n_refs + 1, dtype=np.int32)
        ref_offsets[0] = 0
        np.cumsum(counts, out=ref_offsets[1:])
        if int(ref_offsets[-1]) != n_regions:
            raise RuntimeError(  # pragma: no cover — invariant guard
                "RegionArrays.from_frame: ref_offsets sum mismatch."
            )

        return cls(
            ref_id=ref_id,
            start=start,
            end=end,
            signature=signature,
            strand_class=strand_class,
            region_size_bp=(end - start).astype(np.float64, copy=False),
            ref_offsets=ref_offsets,
            n_refs=n_refs,
        )


# NOTE: BoundaryArrays (a per-boundary structural-flags CSR view) and the per-region mature_eligible_{pos,neg}
# columns were deleted in the 2026-07 cleanup — they fed the removed mature/nascent overlay and had no runtime
# reader. The per-boundary annotation flags they mirrored (is_tss/is_tes/is_splice_junction/genomic_sj_strand)
# were dropped from the index schema at the same time (INDEX_FORMAT_VERSION 7); the solver reads junction
# strand from the accumulator splice motif instead.


# ---------------------------------------------------------------------------
# Node ↔ contiguous-edge index mapping
# ---------------------------------------------------------------------------


def node_right_edge(ref_id: np.ndarray) -> np.ndarray:
    """``int64[N]`` — the contiguous edge to the right of each node, ``-1`` at a reference's last node.

    Nodes ``r`` and ``r + 1`` share a line exactly when they are in the same reference, so the edge
    axis is the run of adjacent same-reference pairs, numbered in node order. A reference with one
    node owns no edge; an empty reference contributes nothing.

    ``ref_id`` must be **grouped** (all of a reference's nodes contiguous), which
    :class:`RegionArrays` guarantees by sorting on ``(ref_id, start)``. Ungrouped input would
    manufacture edges that straddle references, so it is refused rather than tolerated.
    """
    ref = np.asarray(ref_id)
    n = int(ref.shape[0])
    out = np.full(n, -1, dtype=np.int64)
    if n < 2:
        return out
    same = ref[:-1] == ref[1:]
    # grouped ⇔ each reference's rows form ONE run ⇔ the number of runs equals the number of
    # distinct references. `np.unique` counts distinct values; `same` counts run breaks.
    if int((~same).sum()) + 1 != int(np.unique(ref).shape[0]):
        raise ValueError(
            "ref_id is not grouped: each reference's nodes must be contiguous. Build the geometry "
            "with RegionArrays.from_index / from_frame, which sorts on (ref_id, start)."
        )
    out[:-1][same] = np.arange(int(same.sum()), dtype=np.int64)
    return out


def edge_node_indices(ref_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(lo_node, hi_node)`` — the two nodes contiguous edge ``e`` lies between, each ``int64[E]``.

    The exact inverse of :func:`node_right_edge`. ``hi_node == lo_node + 1`` always, and both are in
    the same reference by construction — there is no terminal case and no ``-1``.
    """
    right = node_right_edge(ref_id)
    lo = np.flatnonzero(right >= 0).astype(np.int64)
    return lo, lo + 1
