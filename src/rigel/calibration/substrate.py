"""CalibrationSubstrate — the calibrator-facing view of the accumulator payload.

The substrate is the only object that knows the 4-channel encoding and the
payload topology. It reduces the raw payload into the per-region sufficient
statistics the inference consumes (doc 01 §1, doc 03 §2), exposing **three
parallel per-region views** — *contained*, *left*, *right* — so the E-step
(PR 4) runs the identical code three times.

Channel encoding (``ch = (spliced?2:0) + (strand_pos?0:1)``)::

    ch0 = unspliced & strand+      ch2 = spliced & strand+
    ch1 = unspliced & strand−      ch3 = spliced & strand−

Each view carries integer **counts** (``n_unspliced``, ``n_spliced``,
``k_plus`` = sense among unspliced — the NB/BB likelihood channel) and the
float **mass** magnitude to be deconvolved. For the *contained* view the mass
equals the integer count; for the boundary views the counts come from the
shared integer ``flux`` and the mass from the fractional ``mass_left`` /
``mass_right`` (the D1 side-attribution; see
``00_implementation_plan.md`` §4 D1).

**Strand decodability (NONE vs AMBIG — important).** ``k_plus`` orients the
strand channel to the region's transcript strand (``ts_class``). It is
well-defined only for single-strand regions:

* ``TS_POS`` / ``TS_NEG`` — sense is ``strand+`` / ``strand−`` respectively.
* ``TS_NONE`` (intergenic, no transcript) — gDNA is unstranded, so an
  arbitrary fixed choice (sense = ``strand+``) is **safe**; the strand channel
  contributes ~no net evidence.
* ``TS_AMBIG`` (overlapping transcripts on **both** strands) — every read is
  sense for one transcript and antisense for the other, so **no valid sense
  split exists**. ``k_plus`` for these regions is a non-meaningful placeholder
  and **must not feed strand inference**. AMBIG regions are excluded from the
  strand-balance fit (PR 3) and from the E-step strand log-BF (PR 4); they are
  deconvolved by the count/density channel + the boundary-sweep imputation +
  global fallback (see ``00_implementation_plan.md`` §4 D7). ``ts_class``
  carries this 4-way distinction so downstream can route correctly; the
  strand-agnostic ``n_unspliced`` / ``n_spliced`` / ``mass_*`` remain valid
  for AMBIG.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays, region_boundary_indices
from .signature import TS_NEG


def _reduce_counts(
    counts4: np.ndarray, sense_is_pos: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reduce a 4-channel integer count matrix to (n_unspliced, n_spliced, k_plus)."""
    c = counts4.astype(np.int64, copy=False)
    n_unspliced = c[:, 0] + c[:, 1]
    n_spliced = c[:, 2] + c[:, 3]
    # k_plus = sense-strand count among unspliced. sense = strand+ for TS_POS,
    # strand− for TS_NEG. TS_NONE uses strand+ (arbitrary, safe — gDNA is
    # unstranded). TS_AMBIG also falls here but its value is a non-meaningful
    # placeholder: AMBIG has no valid sense and is excluded from strand
    # inference downstream (keyed off ts_class; see module docstring / D7).
    k_plus = np.where(sense_is_pos, c[:, 0], c[:, 1])
    return n_unspliced, n_spliced, k_plus


def _reduce_mass(mass4: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reduce a 4-channel float mass matrix to (mass_unspliced, mass_spliced)."""
    m = mass4.astype(np.float64, copy=False)
    return m[:, 0] + m[:, 1], m[:, 2] + m[:, 3]


@dataclass(frozen=True, slots=True)
class SubstrateView:
    """One per-region view's sufficient statistics (all length R)."""

    n_unspliced: np.ndarray  # int64[R]  — count channel (NB + BB)
    n_spliced: np.ndarray  # int64[R]
    k_plus: np.ndarray  # int64[R]  — sense among unspliced
    mass_unspliced: np.ndarray  # float64[R] — magnitude to deconvolve
    mass_spliced: np.ndarray  # float64[R]


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Per-region sufficient statistics for the three model views.

    All arrays are length R, aligned to ``RegionArrays`` order. The boundary
    views are already projected onto regions via the D1 side-attribution:
    region ``r`` receives its left boundary's right-side mass and its right
    boundary's left-side mass.
    """

    n_regions: int
    L_eff: np.ndarray  # float64[R] — v1: physical region length (bp)
    ts_class: np.ndarray  # int8[R]    — region transcript-strand class
    contained: SubstrateView
    left: SubstrateView
    right: SubstrateView

    @classmethod
    def from_payload(
        cls,
        payload: AccumulatorPayload,
        region_arrays: RegionArrays,
    ) -> "CalibrationSubstrate":
        cls._check_alignment(payload, region_arrays)
        r = region_arrays.n_regions

        ts_class = np.ascontiguousarray(region_arrays.ts_class, dtype=np.int8)
        # Orientation for the strand channel. TS_NEG → sense is strand−; all
        # others (POS, NONE-arbitrary, AMBIG-placeholder) → strand+. AMBIG's
        # value is excluded downstream via ts_class (see module docstring).
        sense_is_pos = ts_class != TS_NEG

        # Contained view: integer counts; mass equals the count.
        c_u, c_s, c_k = _reduce_counts(payload.region_contained, sense_is_pos)
        contained = SubstrateView(
            n_unspliced=c_u,
            n_spliced=c_s,
            k_plus=c_k,
            mass_unspliced=c_u.astype(np.float64),
            mass_spliced=c_s.astype(np.float64),
        )

        # Boundary projection (D1): region r ← its left boundary's right side
        # + its right boundary's left side. Counts come from the shared integer
        # flux; mass from the fractional side mass. Orient k_plus to region r.
        left_b, right_b = region_boundary_indices(
            payload.ref_region_offsets, payload.ref_boundary_offsets
        )
        left = cls._boundary_view(
            payload.boundary_flux[left_b], payload.boundary_mass_right[left_b], sense_is_pos
        )
        right = cls._boundary_view(
            payload.boundary_flux[right_b], payload.boundary_mass_left[right_b], sense_is_pos
        )

        return cls(
            n_regions=r,
            L_eff=np.ascontiguousarray(region_arrays.region_size_bp, dtype=np.float64),
            ts_class=ts_class,
            contained=contained,
            left=left,
            right=right,
        )

    @staticmethod
    def _boundary_view(
        flux4: np.ndarray, mass4: np.ndarray, sense_is_pos: np.ndarray
    ) -> SubstrateView:
        n_u, n_s, k = _reduce_counts(flux4, sense_is_pos)
        m_u, m_s = _reduce_mass(mass4)
        return SubstrateView(
            n_unspliced=n_u, n_spliced=n_s, k_plus=k, mass_unspliced=m_u, mass_spliced=m_s
        )

    @staticmethod
    def _check_alignment(payload: AccumulatorPayload, region_arrays: RegionArrays) -> None:
        """Enforce that the region geometry addresses the payload 1:1 (doc 04 §4.1)."""
        if region_arrays.n_regions != payload.r_total:
            raise CalibrationSubstrateError(
                f"region geometry has {region_arrays.n_regions} regions but the payload "
                f"has {payload.r_total}; rebuild the index."
            )
        expected = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
        actual = np.asarray(payload.ref_region_offsets, dtype=np.int64)
        if not np.array_equal(expected, actual):
            raise CalibrationSubstrateError(
                "region geometry per-reference offsets do not match the payload "
                "ref_region_offsets; the accumulator region order differs from "
                "index.region_df. Rebuild the index."
            )


__all__ = ["SubstrateView", "CalibrationSubstrate"]
