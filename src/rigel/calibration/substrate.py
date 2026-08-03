"""rigel.calibration.substrate — the calibrator-facing view of the accumulator payload.

       Gate: ``tests/calibration/test_substrate.py``

The substrate is the **only** object that knows the payload's encoding. It decodes the fixed point,
widens the integer banks, and hands the calibrator five populations on three axes. Nothing downstream
reads the payload.

THE FIVE POPULATIONS, on three axes off by one from each other per reference::

    nodes             contained        the whole path lies inside the node
                      spanning         one segment covers the node whole
    contiguous edges  unspliced        the mixture being deconvolved
                      spliced          certified RNA -- gDNA cannot be spliced
    junction edges    (one)            pure RNA by construction

⭐ **ONE type, and that is the change.** The predecessor had ``CalibrationSubstrate`` holding three
per-region views (contained / left / right) and ``BoundarySubstrate`` holding the same numbers re-keyed
by boundary. Two classes, one set of numbers, two keyings — and they existed **solely because a boundary
had two sides**. A contiguous edge is a 0-bp line with one set of numbers, so the second class, the
left/right axis, the re-keying identity and ``_make_view`` all dissolve together.

⭐ **THE COLUMNS ARE GENOME STRAND, WITHOUT EXCEPTION.** Sense/antisense is transcript-relative,
derived by a consumer from a junction's own strand, and never stored. The predecessor stored some banks
by genome strand and others by sense, which is how 40–44 % of ``node_spanning`` deposits landed in the
opposite column from the unspliced fragments beside them.

⚠ **THE FIXED POINT IS DECODED HERE, AND ONLY HERE.** ``inv_length_sum`` arrives as
``round(2^32 / placements)`` and leaves as a real number. Doing it at the boundary is what stops every
downstream module needing to know the scale — one place to be wrong instead of many. ``length_sum`` is a
plain sum of integers and is **not** scaled; applying the density decode to it would divide by 4.3
billion and read as zero everywhere.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays

#: The accumulator's fixed-point scale (``rigel::accumulator::kInvLengthScale``,
#: ``_accumulator_reference.INV_LENGTH_SCALE``). Named here because this module is the only decoder.
INV_LENGTH_SCALE = float(1 << 32)

#: A geometry/payload mismatch surfaces as a shape error deep in the solver, which points nowhere near
#: its cause. Name it at the door instead.
PARTITION_MISMATCH_HINT = (
    "Build the geometry with RegionArrays.from_index(index), and re-scan if the payload was cached "
    "against a different index (a cache is keyed by graph_hash, reach_digest and the payload schema)."
)


@dataclass(frozen=True, slots=True)
class PopulationView:
    """One population's three sums, per object and per **genome strand**.

    The three answer different questions and are never interchangeable: ``count`` carries the statistical
    power (a Beta-Binomial needs an integer), ``inv_length_sum`` carries the level — an exact model-free
    density at an edge, and *not* a density at a node — and ``length_sum`` is the second length tilt,
    which is what makes two components with the same mean fragment length separable at all.
    """

    count: np.ndarray  # int64[n, 2] — genome strand: POS then NEG
    inv_length_sum: np.ndarray  # float64[n, 2] — DECODED from the fixed point
    length_sum: np.ndarray  # float64[n, 2] — Sum L, unscaled

    @property
    def n_objects(self) -> int:
        return int(self.count.shape[0])

    @property
    def total_count(self) -> np.ndarray:
        """int64[n] — both strands. Strand-agnostic magnitude."""
        return self.count.sum(axis=1)

    @property
    def total_inv_length_sum(self) -> np.ndarray:
        return self.inv_length_sum.sum(axis=1)

    @property
    def mean_length(self) -> np.ndarray:
        """float64[n] — this object's own mean fragment length, ``length_sum / count``.

        ⚠ **NaN where the count is zero, deliberately.** An object with no fragments has no mean length,
        and that is not 0 — zero would read as "the fragments here are infinitely short" and propagate as
        a confident wrong answer: an object with no opportunity must
        emit nothing, never a floored value.
        """
        count = self.total_count.astype(np.float64)
        out = np.full(count.shape, np.nan, dtype=np.float64)
        np.divide(self.length_sum.sum(axis=1), count, out=out, where=count > 0)
        return out


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Every per-object statistic the calibrator reads, on the payload's own axes."""

    n_nodes: int
    n_edges: int
    n_junctions: int

    strand_class: np.ndarray  # int8[n_nodes] — the node's transcript-strand class
    node_start_count: np.ndarray  # int64[n_nodes] — one per accepted fragment; THE invariant

    node_contained: PopulationView
    node_spanning: PopulationView
    edge_unspliced: PopulationView
    edge_spliced: PopulationView
    junction: PopulationView

    @classmethod
    def from_payload(
        cls, payload: AccumulatorPayload, region_arrays: RegionArrays
    ) -> "CalibrationSubstrate":
        cls._check_alignment(payload, region_arrays)

        def view(count: np.ndarray, inv: np.ndarray, length: np.ndarray) -> PopulationView:
            return PopulationView(
                count=np.asarray(count, dtype=np.int64),
                inv_length_sum=np.asarray(inv, dtype=np.float64) / INV_LENGTH_SCALE,
                length_sum=np.asarray(length, dtype=np.float64),
            )

        return cls(
            n_nodes=payload.n_nodes,
            n_edges=payload.n_edges,
            n_junctions=payload.n_sj,
            strand_class=np.ascontiguousarray(region_arrays.strand_class, dtype=np.int8),
            node_start_count=np.asarray(payload.node_start_count, dtype=np.int64),
            node_contained=view(
                payload.node_contained_count,
                payload.node_contained_inv_length_sum,
                payload.node_contained_length_sum,
            ),
            node_spanning=view(
                payload.node_spanning_count,
                payload.node_spanning_inv_length_sum,
                payload.node_spanning_length_sum,
            ),
            edge_unspliced=view(
                payload.edge_unspliced_count,
                payload.edge_unspliced_inv_length_sum,
                payload.edge_unspliced_length_sum,
            ),
            edge_spliced=view(
                payload.edge_spliced_count,
                payload.edge_spliced_inv_length_sum,
                payload.edge_spliced_length_sum,
            ),
            junction=view(payload.sj_count, payload.sj_inv_length_sum, payload.sj_length_sum),
        )

    @staticmethod
    def _check_alignment(payload: AccumulatorPayload, region_arrays: RegionArrays) -> None:
        """Enforce that the geometry addresses the payload 1:1 — THE one copy of this invariant.

        ⚠ **Matching totals are not sufficient evidence.** The per-reference offsets are checked too,
        because a geometry can have the right object COUNT while slicing it differently — the defect that
        once dropped 476,719 of 476,732 real fragments inside ``deposit()`` while every golden test
        passed.
        """
        if payload is None:
            raise CalibrationSubstrateError(
                "calibration payload is None; BamScanner.set_regions was not called."
            )
        if region_arrays.n_regions != payload.n_nodes:
            raise CalibrationSubstrateError(
                f"node geometry has {region_arrays.n_regions} objects but the payload has "
                f"{payload.n_nodes}. {PARTITION_MISMATCH_HINT}"
            )
        expected = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
        actual = np.asarray(payload.ref_node_offsets, dtype=np.int64)
        if not np.array_equal(expected, actual):
            bad = int(np.argmax(expected != actual)) if expected.shape == actual.shape else -1
            raise CalibrationSubstrateError(
                "node geometry per-reference offsets do not match the payload's ref_node_offsets"
                + (f" (first difference at reference {bad})" if bad >= 0 else "")
                + f". {PARTITION_MISMATCH_HINT}"
            )


__all__ = ["INV_LENGTH_SCALE", "CalibrationSubstrate", "PopulationView"]
