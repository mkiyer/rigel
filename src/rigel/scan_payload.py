"""rigel.scan_payload — typed dataclass for the BamScanner calibration payload.

The C++ scanner returns ``result["calibration"]`` as a dict of flat ndarrays
(see ``BamScanner::build_result``). ``AccumulatorPayload.from_scan_result``
reshapes those views into the canonical (R_total, 4) / (B_obj_total, 4)
matrices and validates the schema invariants.

Schema (concatenated across all references managed by the scanner):

* ``boundaries`` :: int64[B_pos_total]
    Flat sorted boundary positions, ref-major. The slice for reference ``f``
    is ``boundaries[ref_pos_offsets[f]:ref_pos_offsets[f+1]]``.
* ``ref_pos_offsets`` :: int64[n_refs + 1]
    Offsets into ``boundaries``; ``ref_pos_offsets[n_refs] == B_pos_total``.
* ``ref_region_offsets`` :: int64[n_refs + 1]
    Offsets into ``region_contained``; ``ref_region_offsets[n_refs] == R_total``.
* ``ref_boundary_offsets`` :: int64[n_refs + 1]
    Offsets into the boundary tables; ``ref_boundary_offsets[n_refs] == B_obj_total``.
* ``region_contained`` :: uint32[R_total, 4]
    Per-region per-channel contained fragment counts. Channels are
    ``ch = (spliced ? 2 : 0) + (primary ? 0 : 1)`` where ``primary`` is the
    read's '+' genome strand for unspliced reads and motif-relative SENSE
    (``align_strand == sj_strand``) for spliced reads. So:
    ``ch0 = unspliced & genome+``, ``ch1 = unspliced & genome−``,
    ``ch2 = spliced & sense``, ``ch3 = spliced & antisense``. Unspliced
    channels are genome strand (oriented downstream); spliced channels are
    already transcript-relative (valid even in strand-ambiguous regions).
* ``boundary_mass_left`` :: float32[B_obj_total, 4]
* ``boundary_mass_right`` :: float32[B_obj_total, 4]
* ``boundary_flux_left`` :: uint32[B_obj_total, 4]
* ``boundary_flux_right`` :: uint32[B_obj_total, 4]
    Per-side integer fragment-event counts. A contiguous (unspliced) crossing
    credits both sides of its one boundary; a spliced intron-skip credits one
    side of each flanking boundary (no false exon-intron flux).
* ``n_channels`` :: int (always 4)

Invariant: ``R_total + n_refs == B_obj_total`` because each reference with
``k`` regions contributes ``k`` regions and ``k + 1`` boundaries; for a ref
with zero regions the contribution is ``(0, 0)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


N_CHANNELS = 4


@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    boundaries: np.ndarray
    ref_pos_offsets: np.ndarray
    ref_region_offsets: np.ndarray
    ref_boundary_offsets: np.ndarray
    region_contained: np.ndarray
    boundary_mass_left: np.ndarray
    boundary_mass_right: np.ndarray
    boundary_flux_left: np.ndarray
    boundary_flux_right: np.ndarray
    n_refs: int

    @property
    def n_channels(self) -> int:
        return N_CHANNELS

    @property
    def r_total(self) -> int:
        return int(self.region_contained.shape[0])

    @property
    def b_obj_total(self) -> int:
        return int(self.boundary_mass_left.shape[0])

    @classmethod
    def from_scan_result(cls, scan_result: dict[str, Any]) -> "AccumulatorPayload":
        cal = scan_result.get("calibration")
        if cal is None:
            raise ValueError(
                "scan_result['calibration'] is None; BamScanner.set_regions was not called"
            )

        n_refs = int(cal["n_refs"])
        n_ch = int(cal["n_channels"])
        if n_ch != N_CHANNELS:
            raise ValueError(f"calibration payload has n_channels={n_ch}, expected {N_CHANNELS}")

        boundaries = np.ascontiguousarray(cal["boundaries"], dtype=np.int64)
        ref_pos_offsets = np.ascontiguousarray(cal["ref_pos_offsets"], dtype=np.int64)
        ref_region_offsets = np.ascontiguousarray(cal["ref_region_offsets"], dtype=np.int64)
        ref_boundary_offsets = np.ascontiguousarray(cal["ref_boundary_offsets"], dtype=np.int64)

        if ref_pos_offsets.shape != (n_refs + 1,):
            raise ValueError(f"ref_pos_offsets shape {ref_pos_offsets.shape} != ({n_refs + 1},)")
        if ref_region_offsets.shape != (n_refs + 1,):
            raise ValueError(
                f"ref_region_offsets shape {ref_region_offsets.shape} != ({n_refs + 1},)"
            )
        if ref_boundary_offsets.shape != (n_refs + 1,):
            raise ValueError(
                f"ref_boundary_offsets shape {ref_boundary_offsets.shape} != ({n_refs + 1},)"
            )
        if int(ref_pos_offsets[-1]) != boundaries.shape[0]:
            raise ValueError("ref_pos_offsets[-1] does not match boundaries length")

        r_total = int(ref_region_offsets[-1])
        b_total = int(ref_boundary_offsets[-1])

        region_contained = np.ascontiguousarray(cal["region_contained"], dtype=np.uint32).reshape(
            r_total, N_CHANNELS
        )
        boundary_mass_left = np.ascontiguousarray(
            cal["boundary_mass_left"], dtype=np.float32
        ).reshape(b_total, N_CHANNELS)
        boundary_mass_right = np.ascontiguousarray(
            cal["boundary_mass_right"], dtype=np.float32
        ).reshape(b_total, N_CHANNELS)
        boundary_flux_left = np.ascontiguousarray(
            cal["boundary_flux_left"], dtype=np.uint32
        ).reshape(b_total, N_CHANNELS)
        boundary_flux_right = np.ascontiguousarray(
            cal["boundary_flux_right"], dtype=np.uint32
        ).reshape(b_total, N_CHANNELS)

        # Invariant: each ref with k regions contributes (k, k+1) or (0, 0).
        expected_b_total = r_total + int(np.sum(np.diff(ref_region_offsets) > 0))
        if expected_b_total != b_total:
            raise ValueError(
                f"R_total + nonempty refs ({expected_b_total}) != B_obj_total ({b_total})"
            )

        return cls(
            boundaries=boundaries,
            ref_pos_offsets=ref_pos_offsets,
            ref_region_offsets=ref_region_offsets,
            ref_boundary_offsets=ref_boundary_offsets,
            region_contained=region_contained,
            boundary_mass_left=boundary_mass_left,
            boundary_mass_right=boundary_mass_right,
            boundary_flux_left=boundary_flux_left,
            boundary_flux_right=boundary_flux_right,
            n_refs=n_refs,
        )


__all__ = ["AccumulatorPayload", "N_CHANNELS"]
