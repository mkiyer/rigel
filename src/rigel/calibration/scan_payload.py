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
    if arr.size and not np.all(np.isfinite(arr)):
        raise ValueError(f"calibration payload: {name!r} contains non-finite values.{_ERR_SUFFIX}")
    if arr.size and np.any(arr < 0):
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
