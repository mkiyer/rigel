"""The calibrator-facing view of the accumulator payload.

Gate: ``tests/calibration/test_substrate.py``.

The substrate is the only object that knows the payload's encoding. It widens the integer banks,
folds the strand axis where the contract says one column, and hands the calibrator four populations
on three axes. Nothing downstream reads the payload.

The four populations sit on axes that are off by one from each other per reference, and they do NOT
all carry the same channels, because a channel is stored only where a named consumer reads it::

    regions                contained   count  inv_opportunity_sum        the 1/(ell-w+1) rule
    contiguous boundaries  unspliced   count  inv_length_sum      mass   the mixture being deconvolved
                           spliced     count                      mass   certified RNA: gDNA cannot splice
    sj boundaries          (one)       count  inv_length_sum      mass   pure RNA by construction

The two reciprocal banks carry two deposit rules with two targets, so they carry two names
(one attribute for both is how the REGION truncation stayed invisible):
``inv_length_sum`` is the boundary/sj crossing rule ``1/(w-1)``, expectation ``rho * P(w>=2) = rho``
on any real library; ``inv_opportunity_sum`` is the region contained rule ``1/(ell-w+1)``,
expectation ``rho * P(w<=ell)``, which is truncated by a per-component pmf functional.

The columns are GENOME strand without exception. Sense/antisense is transcript-relative, is derived
by a consumer from an sj's own strand, and is never stored — storing some banks by genome strand and
others by sense puts spliced and unspliced deposits at the same boundary in opposite columns.

The counts carry both columns; the length moments carry one. Which strand a read aligned to says
nothing about whether the molecule was gDNA or RNA, so the moments are strand-agnostic. The counts
keep both columns because the strand model is a Beta-Binomial over them, per strand. The sj mass
arrives per strand and is folded here, since :attr:`PopulationView.mass` is strand-agnostic by
contract; the per-strand values stay in the payload for the artifact filter that reads them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays

#: A geometry/payload mismatch surfaces as a shape error deep in the solver, which points nowhere near
#: its cause. Name it at the door instead.
PARTITION_MISMATCH_HINT = (
    "Build the geometry with RegionArrays.from_index(index), and re-scan if the payload was cached "
    "against a different index (a cache is keyed by graph_hash, reach_digest and the payload schema)."
)


@dataclass(frozen=True, slots=True)
class PopulationView:
    """One population's sums. ``count`` is per genome strand; the length moments are not.

    They answer different questions and are never interchangeable: ``count`` carries the statistical
    power (a Beta-Binomial needs an integer, per strand) and the reciprocal bank carries the level —
    under two different deposit rules with two different targets, so under two names:
    :attr:`inv_length_sum` is the boundary/sj crossing rule ``1/(w-1)``
    (``E[sum] = rho * P(w>=2) = rho``, an exact model-free density), :attr:`inv_opportunity_sum` the
    region contained rule ``1/(ell-w+1)`` (``E[sum] = rho * P(w<=ell)``, a density SHAPE truncated by
    a per-component pmf functional). A population
    carries exactly one of them.

    A population carries only the channels a named consumer reads; the absent ones are ``None``
    rather than zeros — see :meth:`_require`.
    """

    #: What this population is called, for the error a missing channel raises. A view that cannot say
    #: which population it is turns "no mass here" into a traceback nobody can place.
    name: str
    count: np.ndarray  # int64[n, 2] — genome strand: POS then NEG
    #: float64[n] — the boundary/sj crossing-rule bank ``1/(w-1)``; ``None`` where the population does
    #: not carry it. One column, while ``count`` has two: the length moments are strand-agnostic.
    inv_length_sum: np.ndarray | None = None
    #: float64[n] — the region contained-rule bank ``1/(ell-w+1)``
    #: (payload ``region_contained_inv_opportunity_sum``); ``None`` on every boundary/sj population.
    inv_opportunity_sum: np.ndarray | None = None
    #: float64[n] — the conserved mass, strand-agnostic. It sums to ONE per fragment across the objects
    #: that fragment touched, where ``count`` is ``+1`` on each of them. ``None`` on the two region
    #: populations, which need no such channel: ``region_contained_count`` is already 1 per contained
    #: fragment, i.e. already the conserved region mass.
    mass: np.ndarray | None = None

    @property
    def total_count(self) -> np.ndarray:
        """int64[n] — both strands. Strand-agnostic magnitude."""
        return self.count.sum(axis=1)

    def _require(self, channel: str) -> np.ndarray:
        """The channel, or an error that names the population and says why it is absent.

        A missing channel is None, never an array of zeros. Zeros would be a lie in the type: a
        consumer cannot tell "this population does not measure that" from "it measured it and got
        nothing", and the second is an ordinary, meaningful state. :meth:`mass_per_crossing` keeps
        the same contract for a boundary nothing crossed.
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
    def mass_per_crossing(self) -> np.ndarray:
        """float64[n] — ``mass / count``: the mean conserved fragment-mass one crossing here carries.

        This is what converts an object-INCIDENCE total into a FRAGMENT count. It is 1.0 at a
        boundary whose flanking regions both exceed every fragment length, because a crossing
        fragment can then cross only that one boundary and its whole 1.0 lands there; it falls
        toward the region spacing where they do not, and that gap is the per-boundary inflation of
        the incidence count over the fragment count.

        Where nothing crossed the value is 1.0, the identity, not 0. There is no mass at such a
        boundary to rescale, and a 0 would delete whatever mass the deconvolution placed on a
        boundary the accumulator never saw; the identity is the only value that can neither invent
        nor destroy mass.
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

    strand_class: np.ndarray  # int8[n_regions] — the region's transcript-strand class
    #: int64[n_regions, 2] — the path's FIRST covered base, by genome strand; the column sum equals
    #: qc.deposited, which makes it a ledger. Its opportunity is the region length for every fragment
    #: length, so it is the REGION half of the composition-free total; it is wall-blind only at the
    #: template's DOWNSTREAM end, which the mirror below is side-selected against.
    region_start_count: np.ndarray
    #: int64[n_regions, 2] — the mirror: the path's LAST covered base, summing to qc.deposited too.
    #: Wall-blind only at the template's UPSTREAM end.
    region_end_count: np.ndarray
    #: int64[n_regions, 2] — regions STRICTLY spanned, opportunity ``(w-ell-1)+``, a per-component pmf
    #: functional by design. Its consumers are the ledger invariants: contained <= min(start, end),
    #: and span is identically 0 wherever the region length reaches ``w_max - 1``.
    region_span_count: np.ndarray

    #: Four populations, and they do NOT carry the same channels. A channel is stored where a named
    #: consumer reads it and nowhere else::
    #:
    #:     region_contained     count  inv_opportunity_sum         the 1/(ell-w+1) rule
    #:     boundary_unspliced   count  inv_length_sum       mass   the 1/(w-1) rule
    #:     boundary_spliced     count                       mass   certified RNA, not deconvolved
    #:     sj                   count  inv_length_sum              live in second_pass
    #:
    #: No spliced fragment touches the region axis at all: a spliced fragment can never be *contained*,
    #: because both endpoints of an annotated intron are region bounds.
    region_contained: PopulationView
    boundary_unspliced: PopulationView
    boundary_spliced: PopulationView
    sj: PopulationView

    @classmethod
    def from_payload(
        cls, payload: AccumulatorPayload, region_arrays: RegionArrays
    ) -> "CalibrationSubstrate":
        cls._check_alignment(payload, region_arrays)

        def view(name, count, inv=None, mass=None, inv_opp=None) -> PopulationView:
            # `mass` arrives one-column on two axes and two-column on the sj axis, and is folded to one
            # here. `PopulationView.mass` is strand-agnostic by contract: the mass exists to turn an
            # object-incidence total into a fragment count, a question with no strand in it. The
            # per-strand values are NOT re-exported, because their consumer is artifact filtering, which
            # reads the payload — and a channel with no consumer does not belong in this view.
            m = None
            if mass is not None:
                m = np.asarray(mass, dtype=np.float64)
                if m.ndim == 2:
                    m = m.sum(axis=1)
            return PopulationView(
                name=name,
                count=np.asarray(count, dtype=np.int64),
                inv_length_sum=None if inv is None else np.asarray(inv, dtype=np.float64),
                inv_opportunity_sum=(
                    None if inv_opp is None else np.asarray(inv_opp, dtype=np.float64)
                ),
                # No decode: the accumulator deposits fractions as float64 directly, so there is one
                # numeric convention end to end and nothing here to convert.
                mass=m,
            )

        return cls(
            n_regions=payload.n_regions,
            n_boundaries=payload.n_boundaries,
            n_sj=payload.n_sj,
            strand_class=np.ascontiguousarray(region_arrays.strand_class, dtype=np.int8),
            region_start_count=np.asarray(payload.region_start_count, dtype=np.int64),
            region_end_count=np.asarray(payload.region_end_count, dtype=np.int64),
            region_span_count=np.asarray(payload.region_span_count, dtype=np.int64),
            region_contained=view(
                "region_contained",
                payload.region_contained_count,
                inv_opp=payload.region_contained_inv_opportunity_sum,
            ),
            boundary_unspliced=view(
                "boundary_unspliced",
                payload.boundary_unspliced_count,
                payload.boundary_unspliced_inv_length_sum,
                mass=payload.boundary_unspliced_mass,
            ),
            boundary_spliced=view(
                "boundary_spliced",
                payload.boundary_spliced_count,
                mass=payload.boundary_spliced_mass,
            ),
            sj=view("sj", payload.sj_count, payload.sj_inv_length_sum, mass=payload.sj_mass),
        )

    @staticmethod
    def _check_alignment(payload: AccumulatorPayload, region_arrays: RegionArrays) -> None:
        """Enforce that the geometry addresses the payload 1:1 — the one copy of this invariant.

        Matching totals are not sufficient evidence, so the per-reference offsets are checked too: a
        geometry can have the right object COUNT while slicing it differently across references, and
        that mismatch silently discards deposits rather than raising.
        """
        if payload is None:
            raise CalibrationSubstrateError(
                "calibration payload is None; BamScanner.set_regions was not called."
            )
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
