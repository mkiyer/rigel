"""rigel.calibration.substrate — the calibrator-facing view of the accumulator payload.

       Gate: ``tests/calibration/test_substrate.py``

The substrate is the **only** object that knows the payload's encoding. It decodes the fixed point,
widens the integer banks, and hands the calibrator four populations on three axes. Nothing downstream
reads the payload.

THE FOUR POPULATIONS, on three axes off by one from each other per reference — and ⭐ they do NOT all
carry the same channels, because a channel is stored where a named consumer reads it::

    nodes             contained   count  inv_length_sum  length_sum
    contiguous edges  unspliced   count  inv_length_sum  length_sum   the mixture being deconvolved
                      spliced     count                               certified RNA -- gDNA cannot splice
    junction edges    (one)       count  inv_length_sum               pure RNA by construction

⭐ **ONE type, and that is the change.** The predecessor had ``CalibrationSubstrate`` holding three
per-region views (contained / left / right) and ``BoundarySubstrate`` holding the same numbers re-keyed
by boundary. Two classes, one set of numbers, two keyings — and they existed **solely because a boundary
had two sides**. A contiguous edge is a 0-bp line with one set of numbers, so the second class, the
left/right axis, the re-keying identity and ``_make_view`` all dissolve together.

⭐ **THE COLUMNS ARE GENOME STRAND, WITHOUT EXCEPTION.** Sense/antisense is transcript-relative,
derived by a consumer from a junction's own strand, and never stored. The predecessor stored some banks
by genome strand and others by sense, which is how 40–44 % of the SPLICED deposits landed in the
opposite column from the unspliced fragments beside them at the same line.

⭐ **THE COUNTS CARRY BOTH COLUMNS; THE LENGTH MOMENTS CARRY ONE.** Which strand a read aligned to says
nothing about whether the molecule was gDNA or RNA, so the moments are strand-agnostic and every
consumer summed the two columns before using them. The counts keep both because the strand model is a
Beta-Binomial over them, per strand.

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
    """One population's sums. ``count`` is per **genome strand**; the length moments are not.

    They answer different questions and are never interchangeable: ``count`` carries the statistical
    power (a Beta-Binomial needs an integer, per strand), ``inv_length_sum`` carries the level — an exact
    model-free density at an edge, and *not* a density at a node — and ``length_sum`` is the second
    length tilt, which is what makes two components with the same mean fragment length separable at all.

    ⛔ **A population carries only the channels a named consumer reads.** The two the deconvolution
    consumes carry all three; the certified-RNA banks carry fewer, and the absent ones are ``None``
    rather than zeros — see :meth:`_require`.
    """

    #: What this population is called, for the error a missing channel raises. ⚠ A view that cannot say
    #: which population it is turns "no length_sum here" into a traceback nobody can place.
    name: str
    count: np.ndarray  # int64[n, 2] — genome strand: POS then NEG
    #: float64[n] — DECODED from the fixed point. ⛔ ``None`` where the population does not carry it.
    #: ⭐ ONE column, while ``count`` has two: the length moments are strand-AGNOSTIC, and every consumer
    #: summed the two columns before using them.
    inv_length_sum: np.ndarray | None = None
    #: float64[n] — Sum L, unscaled. ⛔ ``None`` where the population does not carry it.
    length_sum: np.ndarray | None = None
    #: float64[n] — ⭐ **THE CONSERVED MASS**, decoded from the fixed point. Sums to ONE per fragment
    #: across the objects it touched, where ``count`` is ``+1`` on each of them. ⛔ ``None`` on the two
    #: node populations, which need no such channel: ``node_contained_count`` is already 1 per contained
    #: fragment, i.e. already the conserved node mass.
    mass: np.ndarray | None = None

    @property
    def n_objects(self) -> int:
        return int(self.count.shape[0])

    @property
    def total_count(self) -> np.ndarray:
        """int64[n] — both strands. Strand-agnostic magnitude."""
        return self.count.sum(axis=1)

    def _require(self, channel: str) -> np.ndarray:
        """⛔ The channel, or an error that names the population and says why it is absent.

        ⚠ **A missing channel is None, never an array of zeros.** Zeros would be a lie in the type: a
        consumer cannot tell "this population does not measure that" from "it measured it and got
        nothing", and the second is an ordinary, meaningful state. The same contract
        ``PopulationView.mean_length`` keeps for a zero count.
        """
        value = getattr(self, channel)
        if value is None:
            raise CalibrationSubstrateError(
                f"population {self.name!r} does not carry {channel!r}. It is stored only where a named "
                f"consumer reads it: the certified-RNA banks carry no length moments, because nothing "
                f"deconvolves a fragment already known to be RNA."
            )
        return value

    @property
    def total_inv_length_sum(self) -> np.ndarray:
        """float64[n]. ⚠ Kept as a named property even though the channel is already strand-summed —
        it is the name every consumer reads, and renaming it would touch them all to say nothing."""
        return self._require("inv_length_sum")

    @property
    def mass_per_crossing(self) -> np.ndarray:
        """float64[n] — ``mass / count``: the mean conserved fragment-mass ONE crossing here carries.

        ⭐ **This is what converts an object-INCIDENCE total into a FRAGMENT count.** It is 1.0 at a line
        whose flanking nodes both exceed every fragment length — a crossing fragment can only cross that
        one line, so its whole 1.0 lands there — and falls toward the node spacing where they do not.
        That gap is the K-inflation, per line.

        ⛔ **1.0 where nothing crossed**, the identity, not 0. There is no mass at such a line to
        rescale, and a 0 would delete whatever mass the deconvolution placed on a line the accumulator
        never saw — the same contract :meth:`mean_length` keeps for a zero count.
        """
        mass = self._require("mass")
        count = self.total_count.astype(np.float64)
        out = np.ones(count.shape, dtype=np.float64)
        np.divide(mass, count, out=out, where=count > 0)
        return out

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
        np.divide(self._require("length_sum"), count, out=out, where=count > 0)
        return out


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Every per-object statistic the calibrator reads, on the payload's own axes."""

    n_nodes: int
    n_edges: int
    n_junctions: int

    strand_class: np.ndarray  # int8[n_nodes] — the node's transcript-strand class
    node_start_count: np.ndarray  # int64[n_nodes] — one per accepted fragment; THE invariant

    #: ⭐⭐ FOUR populations, and they do NOT carry the same channels. A channel is stored where a named
    #: consumer reads it and nowhere else::
    #:
    #:     node_contained   count  inv_length_sum  length_sum
    #:     edge_unspliced   count  inv_length_sum  length_sum  mass
    #:     edge_spliced     count                              mass   certified RNA — not deconvolved
    #:     junction         count  inv_length_sum                     LIVE in second_pass
    #:
    #: ⚠ A fifth, ``node_spanning``, was removed on evidence. ⛔ Its removal means **no spliced fragment
    #: touches the node axis at all** — a spliced fragment can never be *contained*, because both
    #: endpoints of an annotated intron are cuts.
    node_contained: PopulationView
    edge_unspliced: PopulationView
    edge_spliced: PopulationView
    junction: PopulationView

    @classmethod
    def from_payload(
        cls, payload: AccumulatorPayload, region_arrays: RegionArrays
    ) -> "CalibrationSubstrate":
        cls._check_alignment(payload, region_arrays)

        def view(name, count, inv=None, length=None, mass=None) -> PopulationView:
            return PopulationView(
                name=name,
                count=np.asarray(count, dtype=np.int64),
                inv_length_sum=(
                    None if inv is None else np.asarray(inv, dtype=np.float64) / INV_LENGTH_SCALE
                ),
                length_sum=None if length is None else np.asarray(length, dtype=np.float64),
                # ⚠ The mass is a fixed point at the SAME scale, so it decodes the same way. Doing it
                # here keeps this module the one decoder.
                mass=None if mass is None else np.asarray(mass, dtype=np.float64) / INV_LENGTH_SCALE,
            )

        return cls(
            n_nodes=payload.n_nodes,
            n_edges=payload.n_edges,
            n_junctions=payload.n_sj,
            strand_class=np.ascontiguousarray(region_arrays.strand_class, dtype=np.int8),
            node_start_count=np.asarray(payload.node_start_count, dtype=np.int64),
            node_contained=view(
                "node_contained",
                payload.node_contained_count,
                payload.node_contained_inv_length_sum,
                payload.node_contained_length_sum,
            ),
            edge_unspliced=view(
                "edge_unspliced",
                payload.edge_unspliced_count,
                payload.edge_unspliced_inv_length_sum,
                payload.edge_unspliced_length_sum,
                mass=payload.edge_unspliced_mass,
            ),
            edge_spliced=view(
                "edge_spliced", payload.edge_spliced_count, mass=payload.edge_spliced_mass
            ),
            junction=view("junction", payload.sj_count, payload.sj_inv_length_sum),
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
