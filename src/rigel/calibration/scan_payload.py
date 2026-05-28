"""rigel.calibration.scan_payload \u2014 typed wrapper for the fractional payload.

Phase 3 native cutover replaced the integer 8-mask payload with a
12-channel float32 ``region_counts`` plus 6 FL pools. This module
validates the shape/dtype of the calibration dict returned by
``BamScanner.scan()`` and wraps it in a frozen dataclass.

See ``docs/fineregions/fractional_accumulator_python_cutover.md`` \u00a76.1.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .signature import (
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    channel_index,
)


# Number of fragment-length histogram bins per FL pool. Mirrors
# ``rigel::calibration::CalibrationPayload::kFlBins`` in C++.
FL_HIST_N_BINS = 1024


_ERR_SUFFIX = " Rebuild the index or rerun the scan."


def _check_array(name: str, arr: Any, dtype: type, shape: tuple[int, ...]) -> np.ndarray:
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"calibration payload: {name!r} is not a numpy array.{_ERR_SUFFIX}")
    if arr.dtype != np.dtype(dtype):
        raise ValueError(
            f"calibration payload: {name!r} has dtype {arr.dtype}, "
            f"expected {np.dtype(dtype)}.{_ERR_SUFFIX}"
        )
    if arr.shape != shape:
        raise ValueError(
            f"calibration payload: {name!r} has shape {arr.shape}, expected {shape}.{_ERR_SUFFIX}"
        )
    if not arr.flags.c_contiguous:
        raise ValueError(f"calibration payload: {name!r} is not C-contiguous.{_ERR_SUFFIX}")
    if arr.size and arr.dtype.kind == "f" and not np.all(np.isfinite(arr)):
        raise ValueError(f"calibration payload: {name!r} contains non-finite values.{_ERR_SUFFIX}")
    if arr.size and arr.dtype.kind in "fi" and np.any(arr < 0):
        raise ValueError(f"calibration payload: {name!r} contains negative values.{_ERR_SUFFIX}")
    return arr


@dataclass(frozen=True, slots=True)
class CalibrationScanPayload:
    """Typed view of the fractional C++ calibration scan output."""

    region_counts: np.ndarray  # float32[R, 12]
    channel_mass: np.ndarray  # float64[12]
    signature_mass: np.ndarray  # float64[16]
    fl_pool_mass: np.ndarray  # float64[6, 1024]
    fl_pool_total: np.ndarray  # float64[6]
    # Per-region per-compartment fragment support counts. uint32 is
    # sufficient: per-region per-compartment counts cannot realistically
    # approach 2^32 (bounded by coverage * region_length / fragment_length).
    region_contained_unspliced_support: np.ndarray       # uint32[R]
    region_boundary_left_unspliced_support: np.ndarray   # uint32[R]
    region_boundary_right_unspliced_support: np.ndarray  # uint32[R]
    region_contained_spliced_support: np.ndarray         # uint32[R]
    region_boundary_left_spliced_support: np.ndarray     # uint32[R]
    region_boundary_right_spliced_support: np.ndarray    # uint32[R]
    # Per-region aggregate fragment support, partitioned only by splice
    # class. A fragment is counted at most once per region regardless of
    # how many compartments it touches; a fully-spanning fragment
    # contributes one increment here but two to the per-compartment
    # counters above (left + right). Consumed by the EB exposure /
    # density model.
    region_unspliced_support: np.ndarray  # uint64[R]
    region_spliced_support: np.ndarray    # uint64[R]

    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_excluded_strand_ambig: int
    n_unobserved: int
    n_unannotated_ref: int
    n_fl_unavailable: int

    #: Resolver implicit-splice slack the scanner used. Provenance only;
    #: not a calibration tolerance.
    resolver_splicing_anchor_tolerance: int
    n_regions: int

    @classmethod
    def from_scan_dict(
        cls,
        d: Mapping[str, Any] | None,
        *,
        n_total: int | None = None,
    ) -> "CalibrationScanPayload":
        """Validate a calibration dict produced by ``BamScanner.scan()``."""
        if d is None:
            raise ValueError(
                "calibration payload: scanner returned None \u2014 was "
                f"set_regions() called?{_ERR_SUFFIX}"
            )

        n_regions = int(d["n_regions"])
        if n_regions < 0:
            raise ValueError(f"calibration payload: n_regions ({n_regions}) < 0.{_ERR_SUFFIX}")

        region_counts = _check_array(
            "region_counts", d["region_counts"], np.float32, (n_regions, N_CHANNELS)
        )
        channel_mass = _check_array("channel_mass", d["channel_mass"], np.float64, (N_CHANNELS,))
        signature_mass = _check_array(
            "signature_mass", d["signature_mass"], np.float64, (N_SIGNATURES,)
        )
        fl_pool_mass = _check_array(
            "fl_pool_mass", d["fl_pool_mass"], np.float64, (N_FL_POOLS, FL_HIST_N_BINS)
        )
        fl_pool_total = _check_array("fl_pool_total", d["fl_pool_total"], np.float64, (N_FL_POOLS,))

        def _check_support(name: str) -> np.ndarray:
            return _check_array(name, d[name], np.uint32, (n_regions,))

        region_contained_unspliced_support = _check_support(
            "region_contained_unspliced_support"
        )
        region_boundary_left_unspliced_support = _check_support(
            "region_boundary_left_unspliced_support"
        )
        region_boundary_right_unspliced_support = _check_support(
            "region_boundary_right_unspliced_support"
        )
        region_contained_spliced_support = _check_support(
            "region_contained_spliced_support"
        )
        region_boundary_left_spliced_support = _check_support(
            "region_boundary_left_spliced_support"
        )
        region_boundary_right_spliced_support = _check_support(
            "region_boundary_right_spliced_support"
        )
        region_unspliced_support = _check_array(
            "region_unspliced_support",
            d["region_unspliced_support"],
            np.uint64,
            (n_regions,),
        )
        region_spliced_support = _check_array(
            "region_spliced_support",
            d["region_spliced_support"],
            np.uint64,
            (n_regions,),
        )

        n_observed = int(d["n_observed"])
        n_ex_mm = int(d["n_excluded_multimap"])
        n_ex_ch = int(d["n_excluded_chimera"])
        n_ex_ar = int(d["n_excluded_artifact"])
        n_ex_sa = int(d["n_excluded_strand_ambig"])
        n_unobs = int(d["n_unobserved"])
        n_unann = int(d["n_unannotated_ref"])
        n_flu = int(d["n_fl_unavailable"])
        resolver_K = int(d["resolver_splicing_anchor_tolerance"])

        for label, val in (
            ("n_observed", n_observed),
            ("n_excluded_multimap", n_ex_mm),
            ("n_excluded_chimera", n_ex_ch),
            ("n_excluded_artifact", n_ex_ar),
            ("n_excluded_strand_ambig", n_ex_sa),
            ("n_unobserved", n_unobs),
            ("n_unannotated_ref", n_unann),
            ("n_fl_unavailable", n_flu),
            ("resolver_splicing_anchor_tolerance", resolver_K),
        ):
            if val < 0:
                raise ValueError(f"calibration payload: {label} ({val}) < 0.{_ERR_SUFFIX}")

        # Internal consistency: n_unannotated_ref and n_fl_unavailable are
        # subsets of n_observed.
        if n_unann > n_observed:
            raise ValueError(
                f"calibration payload: n_unannotated_ref ({n_unann}) > "
                f"n_observed ({n_observed}).{_ERR_SUFFIX}"
            )
        if n_flu > n_observed:
            raise ValueError(
                f"calibration payload: n_fl_unavailable ({n_flu}) > "
                f"n_observed ({n_observed}).{_ERR_SUFFIX}"
            )

        # Mass identities. region_counts.sum() should equal
        # channel_mass.sum() == signature_mass.sum() within float32 noise.
        region_total = float(region_counts.sum(dtype=np.float64))
        channel_total = float(channel_mass.sum())
        signature_total = float(signature_mass.sum())
        expected_region_mass = float(n_observed - n_unann)
        # Tolerance: float32 has ~1e-7 relative precision; allow generous
        # slack proportional to the observed mass plus a small absolute floor.
        tol = max(1e-3, 1e-5 * max(region_total, 1.0))
        if abs(region_total - expected_region_mass) > max(1.0, tol):
            raise ValueError(
                f"calibration payload: sum(region_counts) = {region_total} "
                f"!= n_observed - n_unannotated_ref = {expected_region_mass} "
                f"(tol={max(1.0, tol)}).{_ERR_SUFFIX}"
            )
        if abs(channel_total - region_total) > tol:
            raise ValueError(
                f"calibration payload: sum(channel_mass) = {channel_total} "
                f"!= sum(region_counts) = {region_total} (tol={tol}).{_ERR_SUFFIX}"
            )
        if abs(signature_total - region_total) > tol:
            raise ValueError(
                f"calibration payload: sum(signature_mass) = {signature_total} "
                f"!= sum(region_counts) = {region_total} (tol={tol}).{_ERR_SUFFIX}"
            )

        # FL pool internal identity is exact (both built in the same
        # double-precision pass).
        pool_row_sums = fl_pool_mass.sum(axis=1)
        if not np.allclose(pool_row_sums, fl_pool_total, rtol=1e-12, atol=1e-9):
            raise ValueError(
                "calibration payload: fl_pool_mass.sum(axis=1) does not match "
                f"fl_pool_total.{_ERR_SUFFIX}"
            )

        # FL pools cover only unspliced observed fragments with usable FL.
        # Strand-ambiguous fragments are excluded before n_observed is
        # incremented in native code, so they are not subtracted here.
        # Bound generously to absorb float64 rounding.
        fl_total = float(fl_pool_total.sum())
        fl_upper = float(n_observed)
        if fl_total > fl_upper + 1e-3:
            raise ValueError(
                f"calibration payload: sum(fl_pool_total) = {fl_total} > "
                f"n_observed = {fl_upper}.{_ERR_SUFFIX}"
            )

        # Per-region per-compartment support invariants. Each compartment
        # counts each fragment at most once, so no per-cell value can
        # exceed n_observed. (A fully-spanning fragment increments both
        # left and right for the same region; that's still one increment
        # per cell, not two of the same.)
        compartment_supports = (
            ("region_contained_unspliced_support", region_contained_unspliced_support),
            ("region_boundary_left_unspliced_support", region_boundary_left_unspliced_support),
            ("region_boundary_right_unspliced_support", region_boundary_right_unspliced_support),
            ("region_contained_spliced_support", region_contained_spliced_support),
            ("region_boundary_left_spliced_support", region_boundary_left_spliced_support),
            ("region_boundary_right_spliced_support", region_boundary_right_spliced_support),
        )
        for label, arr in compartment_supports:
            if arr.size and int(arr.max()) > n_observed:
                raise ValueError(
                    f"calibration payload: {label} contains a value "
                    f"> n_observed ({n_observed}).{_ERR_SUFFIX}"
                )
        # Aggregate per-region support: distinct fragments touching the
        # region, bounded by n_observed.
        if region_unspliced_support.size and int(region_unspliced_support.max()) > n_observed:
            raise ValueError(
                "calibration payload: region_unspliced_support contains a value "
                f"> n_observed ({n_observed}).{_ERR_SUFFIX}"
            )
        if region_spliced_support.size and int(region_spliced_support.max()) > n_observed:
            raise ValueError(
                "calibration payload: region_spliced_support contains a value "
                f"> n_observed ({n_observed}).{_ERR_SUFFIX}"
            )

        # Per-region per-compartment invariant: positive fractional mass
        # in a (compartment, splice) cell implies positive support in the
        # matching support cell. This is a stricter check than the v6
        # aggregate version: support and mass agree per compartment.
        if n_regions > 0:
            from .signature import (
                CHANNEL_STRAND_NEG,
                CHANNEL_STRAND_POS,
                COMPARTMENT_BOUNDARY_LEFT,
                COMPARTMENT_BOUNDARY_RIGHT,
                COMPARTMENT_CONTAINED,
                SPLICE_SPLICED,
                SPLICE_UNSPLICED,
            )

            strands = (CHANNEL_STRAND_POS, CHANNEL_STRAND_NEG)

            def _cell_mass(compartment: int, splice: int) -> np.ndarray:
                idxs = [channel_index(compartment, splice, s) for s in strands]
                return region_counts[:, idxs].sum(axis=1, dtype=np.float64)

            float32_eps = np.finfo(np.float32).eps
            checks = (
                (
                    "contained_unspliced",
                    _cell_mass(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED),
                    region_contained_unspliced_support,
                ),
                (
                    "boundary_left_unspliced",
                    _cell_mass(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED),
                    region_boundary_left_unspliced_support,
                ),
                (
                    "boundary_right_unspliced",
                    _cell_mass(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED),
                    region_boundary_right_unspliced_support,
                ),
                (
                    "contained_spliced",
                    _cell_mass(COMPARTMENT_CONTAINED, SPLICE_SPLICED),
                    region_contained_spliced_support,
                ),
                (
                    "boundary_left_spliced",
                    _cell_mass(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED),
                    region_boundary_left_spliced_support,
                ),
                (
                    "boundary_right_spliced",
                    _cell_mass(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED),
                    region_boundary_right_spliced_support,
                ),
            )
            for label, mass, support in checks:
                thresh = float32_eps * max(float(mass.sum()), 1.0)
                bad = (support == 0) & (mass > thresh)
                if bool(bad.any()):
                    idx = int(np.argmax(bad))
                    raise ValueError(
                        "calibration payload: region "
                        f"{idx} has {label} fractional mass {float(mass[idx])} > 0 "
                        f"but region_{label}_support == 0. Support and mass "
                        f"must agree per compartment.{_ERR_SUFFIX}"
                    )

        if n_total is not None:
            accounted = n_observed + n_ex_mm + n_ex_ch + n_ex_ar + n_ex_sa + n_unobs
            if accounted != int(n_total):
                raise ValueError(
                    "calibration payload: balance assertion failed: "
                    f"n_observed ({n_observed}) + n_excluded_multimap "
                    f"({n_ex_mm}) + n_excluded_chimera ({n_ex_ch}) + "
                    f"n_excluded_artifact ({n_ex_ar}) + "
                    f"n_excluded_strand_ambig ({n_ex_sa}) + n_unobserved "
                    f"({n_unobs}) = {accounted} != n_total "
                    f"({int(n_total)}).{_ERR_SUFFIX}"
                )

        return cls(
            region_counts=region_counts,
            channel_mass=channel_mass,
            signature_mass=signature_mass,
            fl_pool_mass=fl_pool_mass,
            fl_pool_total=fl_pool_total,
            region_contained_unspliced_support=region_contained_unspliced_support,
            region_boundary_left_unspliced_support=region_boundary_left_unspliced_support,
            region_boundary_right_unspliced_support=region_boundary_right_unspliced_support,
            region_contained_spliced_support=region_contained_spliced_support,
            region_boundary_left_spliced_support=region_boundary_left_spliced_support,
            region_boundary_right_spliced_support=region_boundary_right_spliced_support,
            region_unspliced_support=region_unspliced_support,
            region_spliced_support=region_spliced_support,
            n_observed=n_observed,
            n_excluded_multimap=n_ex_mm,
            n_excluded_chimera=n_ex_ch,
            n_excluded_artifact=n_ex_ar,
            n_excluded_strand_ambig=n_ex_sa,
            n_unobserved=n_unobs,
            n_unannotated_ref=n_unann,
            n_fl_unavailable=n_flu,
            resolver_splicing_anchor_tolerance=resolver_K,
            n_regions=n_regions,
        )


__all__ = [
    "CalibrationScanPayload",
    "FL_HIST_N_BINS",
    "N_CHANNELS",
    "N_FL_POOLS",
    "N_SIGNATURES",
    "channel_index",
]
