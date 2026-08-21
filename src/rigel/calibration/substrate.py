"""rigel.calibration.substrate — the calibrator-facing view of the accumulator payload.

       Gate: ``tests/calibration/test_substrate.py``

The substrate is the **only** object that knows the payload's encoding. It decodes the fixed point,
widens the integer banks, and hands the calibrator four populations on three axes. Nothing downstream
reads the payload.

THE FOUR POPULATIONS, on three axes off by one from each other per reference — and ⭐ they do NOT all
carry the same channels, because a channel is stored where a named consumer reads it::

    regions             contained   count  inv_opportunity_sum            the 1/(ell-w+1) rule
    contiguous boundaries  unspliced   count  inv_length_sum  mass        the mixture being deconvolved
                      spliced     count                  mass        certified RNA -- gDNA cannot splice
    sj boundaries    (one)       count  inv_length_sum  mass        pure RNA by construction

⛔ **The two reciprocal banks carry TWO deposit rules with TWO targets, so they carry two names**
(TRAPS: two-masks-one-name — one attribute for both is how the REGION truncation stayed invisible):
``inv_length_sum`` is the BOUNDARY/sj crossing rule ``1/(w−1)``, expectation ``rho * P(w>=2) = rho`` on
any real library; ``inv_opportunity_sum`` is the REGION contained rule ``1/(ell−w+1)``, expectation
``rho * P(w<=ell)`` — truncated by a per-component pmf functional
(TRAPS: a-cancellation-is-conditional-on-its-support).

⭐ **ONE type, and that is the change.** The predecessor had ``CalibrationSubstrate`` holding three
per-region views (contained / left / right) and ``BoundarySubstrate`` holding the same numbers re-keyed
by boundary. Two classes, one set of numbers, two keyings — and they existed **solely because a boundary
had two sides**. A contiguous boundary is a 0-bp boundary with one set of numbers, so the second class, the
left/right axis, the re-keying identity and ``_make_view`` all dissolve together.

⭐ **THE COLUMNS ARE GENOME STRAND, WITHOUT EXCEPTION.** Sense/antisense is transcript-relative,
derived by a consumer from a sj's own strand, and never stored. The predecessor stored some banks
by genome strand and others by sense, which is how 40–44 % of the SPLICED deposits landed in the
opposite column from the unspliced fragments beside them at the same boundary.

⭐ **THE COUNTS CARRY BOTH COLUMNS; THE LENGTH MOMENTS CARRY ONE.** Which strand a read aligned to says
nothing about whether the molecule was gDNA or RNA, so the moments are strand-agnostic and every
consumer summed the two columns before using them. The counts keep both because the strand model is a
Beta-Binomial over them, per strand.

⭐ **AND THE SJ MASS ARRIVES WITH TWO COLUMNS AND LEAVES WITH ONE.** ``sj_mass`` went per-strand
on 2026-08-13 for artifact detection (`SpliceJunction::mass` carries the premise), and this boundary folds
it: :attr:`PopulationView.mass` is strand-agnostic by contract, so nothing downstream of here changed.

⚠ **``length_sum`` and ``mean_length`` are GONE (2026-08-13)** along with the two banks that fed them —
they reached this view and stopped, so the channel had no consumer and its stated justification was
false where it claimed to help (`scan_payload`'s docstring has the retraction).
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
    """One population's sums. ``count`` is per **genome strand**; the length moments are not.

    They answer different questions and are never interchangeable: ``count`` carries the statistical
    power (a Beta-Binomial needs an integer, per strand) and the reciprocal bank carries the level —
    under TWO different deposit rules with TWO different targets, so under TWO names
    (TRAPS: two-masks-one-name): :attr:`inv_length_sum` is the BOUNDARY/sj crossing rule ``1/(w−1)``
    (``E[sum] = rho * P(w>=2) = rho`` — an exact model-free density), :attr:`inv_opportunity_sum` the
    REGION contained rule ``1/(ell−w+1)`` (``E[sum] = rho * P(w<=ell)`` — a density SHAPE truncated by
    a per-component pmf functional, TRAPS: a-cancellation-is-conditional-on-its-support). A population
    carries exactly one of them.

    ⛔ **A population carries only the channels a named consumer reads**, and that rule has teeth: the
    ``length_sum`` channel was deleted in 2026-08-13's schema change precisely because it had none. The
    absent ones are ``None`` rather than zeros — see :meth:`_require`.
    """

    #: What this population is called, for the error a missing channel raises. ⚠ A view that cannot say
    #: which population it is turns "no mass here" into a traceback nobody can place.
    name: str
    count: np.ndarray  # int64[n, 2] — genome strand: POS then NEG
    #: float64[n] — the BOUNDARY/sj crossing-rule bank ``1/(w−1)``. ⛔ ``None`` where the population
    #: does not carry it. ⭐ ONE column, while ``count`` has two: the length moments are strand-AGNOSTIC,
    #: and every consumer summed the two columns before using them.
    inv_length_sum: np.ndarray | None = None
    #: float64[n] — the REGION contained-rule bank ``1/(ell−w+1)``
    #: (payload ``region_contained_inv_opportunity_sum``). ⛔ ``None`` on every boundary/sj population.
    inv_opportunity_sum: np.ndarray | None = None
    #: float64[n] — ⭐ **THE CONSERVED MASS**, strand-agnostic. Sums to ONE per fragment
    #: across the objects it touched, where ``count`` is ``+1`` on each of them. ⛔ ``None`` on the two
    #: region populations, which need no such channel: ``region_contained_count`` is already 1 per contained
    #: fragment, i.e. already the conserved region mass.
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
        :meth:`mass_per_crossing` keeps for a boundary nothing crossed.
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
        """float64[n] — the BOUNDARY/sj crossing-rule bank, already strand-summed. ⚠ A previous
        docstring claimed *"it is the name every consumer reads"* — false: a grep over ``src/`` returns
        this definition and nothing else, and the one live consumer (`region_geometry`) reads the field
        directly. Kept because deleting public surface is its own decision, not a doc fix."""
        return self._require("inv_length_sum")

    @property
    def mass_per_crossing(self) -> np.ndarray:
        """float64[n] — ``mass / count``: the mean conserved fragment-mass ONE crossing here carries.

        ⭐ **This is what converts an object-INCIDENCE total into a FRAGMENT count.** It is 1.0 at a boundary
        whose flanking regions both exceed every fragment length — a crossing fragment can only cross that
        one boundary, so its whole 1.0 lands there — and falls toward the region spacing where they do not.
        That gap is the K-inflation, per boundary.

        ⛔ **1.0 where nothing crossed**, the identity, not 0. There is no mass at such a boundary to
        rescale, and a 0 would delete whatever mass the deconvolution placed on a boundary the accumulator
        never saw — the identity is the only value that cannot invent or destroy mass.
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
    region_start_count: np.ndarray  # int64[n_regions] — one per accepted fragment; THE invariant

    #: ⭐⭐ FOUR populations, and they do NOT carry the same channels. A channel is stored where a named
    #: consumer reads it and nowhere else::
    #:
    #:     region_contained   count  inv_opportunity_sum          the 1/(ell−w+1) rule
    #:     boundary_unspliced   count  inv_length_sum        mass   the 1/(w−1) rule
    #:     boundary_spliced     count                        mass   certified RNA — not deconvolved
    #:     sj         count  inv_length_sum               LIVE in second_pass
    #:
    #: ⚠ A fifth, ``region_spanning``, was removed on evidence. ⛔ Its removal means **no spliced fragment
    #: touches the region axis at all** — a spliced fragment can never be *contained*, because both
    #: endpoints of an annotated intron are region_bounds.
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
            # ⭐⭐ `mass` arrives one-column on two axes and TWO-column on the sj axis, and is
            # folded to one here. `PopulationView.mass` is strand-agnostic by contract — the mass exists
            # to turn an object-incidence total into a fragment count, a question with no strand in it —
            # so folding at the boundary is what keeps every downstream consumer of `sj_mass` unchanged
            # by the strand split. ⛔ The per-strand values are NOT re-exported here: their consumer is
            # artifact filtering, which reads the PAYLOAD, and a channel with no consumer is the defect
            # the two `*_length_sum` banks were deleted for in this same change.
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
                # ⭐ No decode: the accumulator deposits fractions as float64 directly. This module used
                # to be "the one decoder"; under one numeric convention there is nothing to decode.
                mass=m,
            )

        return cls(
            n_regions=payload.n_regions,
            n_boundaries=payload.n_boundaries,
            n_sj=payload.n_sj,
            strand_class=np.ascontiguousarray(region_arrays.strand_class, dtype=np.int8),
            region_start_count=np.asarray(payload.region_start_count, dtype=np.int64),
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
                "boundary_spliced", payload.boundary_spliced_count, mass=payload.boundary_spliced_mass
            ),
            sj=view(
                "sj", payload.sj_count, payload.sj_inv_length_sum, mass=payload.sj_mass
            ),
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
