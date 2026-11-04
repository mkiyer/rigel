"""The calibrator-facing view of the accumulator payload (layer 1).

The substrate is the only object that knows the payload's encoding: it widens the integer banks, folds the strand axis
where a channel is strand-agnostic, and hands the calibrator four populations on three axes. A channel is stored only
where a named consumer reads it::

    regions                contained   count
    contiguous boundaries  unspliced   count  mass   the mixture being deconvolved
                           spliced     count  mass   certified RNA: gDNA cannot splice
    sj boundaries          (one)       count  mass   pure RNA by construction

Count columns are genome strand without exception, and the counts keep both because the strand likelihood is binomial
in them. Sense is transcript-relative and derived by a consumer from an sj's own strand; it is never stored, because
mixing genome-strand and sense banks would put spliced and unspliced deposits at one boundary in opposite columns.

Gate: ``tests/calibration/test_substrate.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays

#: The hint every geometry/payload mismatch raises with, so the error names its cause rather than surfacing as a
#: shape error deep in the solver.
PARTITION_MISMATCH_HINT = (
    "Build the geometry with RegionArrays.from_index(index), and re-scan if the payload was cached "
    "against a different index (a cache is keyed by graph_hash, reach_digest and the payload schema)."
)


@dataclass(frozen=True, slots=True)
class PopulationView:
    """One population's sums. ``count`` is per genome strand; ``mass`` is not.

    A population carries only the channels a named consumer reads; the absent ones are ``None``, not zeros
    (:meth:`_require`).
    """

    #: The population's name, for the error a missing channel raises.
    name: str
    count: np.ndarray  # int64[n, 2] — genome strand: POS then NEG
    #: float64[n]: the conserved mass, strand-agnostic. It sums to one per fragment across the objects the fragment
    #: touched, where ``count`` is +1 on each. ``None`` on the region population, whose count is already one per
    #: contained fragment.
    mass: np.ndarray | None = None

    @property
    def total_count(self) -> np.ndarray:
        """int64[n] — both strands. Strand-agnostic magnitude."""
        return self.count.sum(axis=1)

    def _require(self, channel: str) -> np.ndarray:
        """The channel, or an error that names the population.

        A missing channel is ``None``, never zeros: zeros could not be told apart from "measured and got nothing",
        which is an ordinary state.
        """
        value = getattr(self, channel)
        if value is None:
            raise CalibrationSubstrateError(f"population {self.name!r} does not carry {channel!r}.")
        return value

    @property
    def mass_per_crossing(self) -> np.ndarray:
        """float64[n]: ``mass / count``, the mean conserved fragment mass one crossing here carries.

        It converts an object-incidence total into a fragment count. It is 1.0 at a boundary whose flanking regions
        both exceed every fragment length (a crossing fragment can then cross only that boundary) and falls toward the
        region spacing where they do not. Where nothing crossed it is 1.0, the identity: the only value that neither
        invents nor destroys mass the deconvolution placed there.
        """
        mass = self._require("mass")
        count = self.total_count.astype(np.float64)
        out = np.ones(count.shape, dtype=np.float64)
        np.divide(mass, count, out=out, where=count > 0)
        return out


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Every per-object statistic the calibrator reads, on the payload's own axes."""

    n_regions: int
    n_boundaries: int
    n_sj: int

    #: The four populations of the module docstring's table. No spliced fragment touches the region axis: both
    #: endpoints of an annotated intron are region bounds, so a spliced fragment is never contained.
    region_contained: PopulationView
    boundary_unspliced: PopulationView
    boundary_spliced: PopulationView
    sj: PopulationView

    @classmethod
    def from_payload(
        cls, payload: AccumulatorPayload, region_arrays: RegionArrays
    ) -> "CalibrationSubstrate":
        cls.check_alignment(payload, region_arrays)

        def view(name, count, mass=None) -> PopulationView:
            # the sj mass arrives per strand and is folded to one column: it turns an incidence total into a
            # fragment count, a question with no strand in it
            m = None
            if mass is not None:
                m = np.asarray(mass, dtype=np.float64)
                if m.ndim == 2:
                    m = m.sum(axis=1)
            return PopulationView(
                name=name,
                count=np.asarray(count, dtype=np.int64),
                # the accumulator deposits float64 fractions directly, so there is nothing to decode
                mass=m,
            )

        return cls(
            n_regions=payload.n_regions,
            n_boundaries=payload.n_boundaries,
            n_sj=payload.n_sj,
            region_contained=view("region_contained", payload.region_contained_count),
            boundary_unspliced=view(
                "boundary_unspliced",
                payload.boundary_unspliced_count,
                mass=payload.boundary_unspliced_mass,
            ),
            boundary_spliced=view(
                "boundary_spliced",
                payload.boundary_spliced_count,
                mass=payload.boundary_spliced_mass,
            ),
            sj=view("sj", payload.sj_count, mass=payload.sj_mass),
        )

    @staticmethod
    def check_alignment(payload: AccumulatorPayload, region_arrays: RegionArrays) -> None:
        """Enforce that the geometry addresses the payload 1:1.

        Matching totals are not enough: a geometry can have the right object count but slice it differently across
        references, which would silently discard deposits, so the per-reference offsets are checked too.
        """
        if region_arrays.n_regions != payload.n_regions:
            raise CalibrationSubstrateError(
                f"region geometry has {region_arrays.n_regions} objects but the payload has "
                f"{payload.n_regions}. {PARTITION_MISMATCH_HINT}"
            )
        expected = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
        actual = np.asarray(payload.ref_region_offsets, dtype=np.int64)
        if not np.array_equal(expected, actual):
            bad = int(np.argmax(expected != actual)) if expected.shape == actual.shape else -1
            raise CalibrationSubstrateError(
                "region geometry per-reference offsets do not match the payload's ref_region_offsets"
                + (f" (first difference at reference {bad})" if bad >= 0 else "")
                + f". {PARTITION_MISMATCH_HINT}"
            )


__all__ = ["CalibrationSubstrate", "PopulationView"]
