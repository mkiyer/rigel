"""rigel.calibration.scan_payload — typed wrapper for the C++ calibration dict.

The native ``BamScanner.scan()`` returns a ``calibration`` dict (or
``None`` if ``set_regions`` was not called).  This module validates
the shape/dtype of that dict and wraps it in a frozen dataclass.

The validation includes a balance assertion that the per-class
fragment counters (observed + excluded) account for every input
fragment; this is the master regression guard against double-counting
in the C++ observation site.

See ``docs/calibration/m3_implementation_plan.md`` §3.8.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np


#: Number of fragment-length histogram bins per mask category.  Mirrors
#: ``rigel::calibration::CalibrationPayload::kFlBins`` in C++.  Widening
#: this requires touching both sides.
FL_HIST_N_BINS = 1024

#: Number of mask categories (3-bit code: EXON|INTRON|INTERGENIC).  Mirrors
#: ``rigel::calibration::mask::N_STATES`` in C++.
MASK_N_STATES = 8

#: Public mask bit positions (3-bit code).  Mirrors the ``rigel::calibration::mask``
#: enum in C++.  Single source of truth for the calibration package; do not
#: re-declare these as private ``_MASK_*`` constants in callers.
MASK_EXON = 0b001        # bit 0
MASK_INTRON = 0b010      # bit 1
MASK_INTERGENIC = 0b100  # bit 2


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
            f"calibration payload: {name!r} has shape {arr.shape}, "
            f"expected {shape}.{_ERR_SUFFIX}"
        )
    return arr


@dataclass(frozen=True, slots=True)
class CalibrationScanPayload:
    """Typed view of the C++ calibration scan output."""

    global_counts: np.ndarray             # (MASK_N_STATES,)               int64
    per_region_counts: np.ndarray         # (R, MASK_N_STATES)              int64
    fl_hist: np.ndarray                   # (MASK_N_STATES, FL_HIST_N_BINS) int64
    u_left: np.ndarray                    # (R,)                            int64
    u_right: np.ndarray                   # (R,)                            int64
    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_unobserved: int
    n_unannotated_ref: int

    @classmethod
    def from_scan_dict(
        cls,
        d: Mapping[str, Any],
        *,
        n_total: int | None = None,
    ) -> "CalibrationScanPayload":
        """Validate a calibration dict produced by ``BamScanner.scan()``.

        Parameters
        ----------
        d
            The ``result["calibration"]`` dict from the native scanner.
        n_total
            Optional total qname-group count (``stats.n_read_names``).
            If provided, asserts that
            ``n_observed + n_excluded_multimap + n_excluded_chimera +
            n_excluded_artifact + n_unobserved == n_total``.

        Raises
        ------
        ValueError
            On any shape/dtype/balance violation.
        """
        if d is None:
            raise ValueError(
                "calibration payload: scanner returned None — was "
                f"set_regions() called?{_ERR_SUFFIX}"
            )
        per_region_counts = d["per_region_counts"]
        if not isinstance(per_region_counts, np.ndarray) or per_region_counts.ndim != 2:
            raise ValueError(
                f"calibration payload: 'per_region_counts' must be a 2-D ndarray.{_ERR_SUFFIX}"
            )
        n_regions = int(per_region_counts.shape[0])

        global_counts = _check_array("global_counts", d["global_counts"], np.int64, (MASK_N_STATES,))
        per_region_counts = _check_array(
            "per_region_counts", per_region_counts, np.int64, (n_regions, MASK_N_STATES)
        )
        fl_hist = _check_array("fl_hist", d["fl_hist"], np.int64, (MASK_N_STATES, FL_HIST_N_BINS))
        u_left = _check_array("u_left", d["u_left"], np.int64, (n_regions,))
        u_right = _check_array("u_right", d["u_right"], np.int64, (n_regions,))

        n_observed = int(d["n_observed"])
        n_ex_mm = int(d["n_excluded_multimap"])
        n_ex_ch = int(d["n_excluded_chimera"])
        n_ex_ar = int(d["n_excluded_artifact"])
        n_unobs = int(d["n_unobserved"])
        n_unann = int(d["n_unannotated_ref"])

        # Internal consistency: n_unannotated_ref is a subset of n_observed
        # (it counts observed fragments whose every block fell outside any
        # region — mask == 0 after the overlap pass).
        if n_unann > n_observed:
            raise ValueError(
                f"calibration payload: n_unannotated_ref ({n_unann}) > "
                f"n_observed ({n_observed}).{_ERR_SUFFIX}"
            )

        # global_counts must sum to n_observed.
        gc_sum = int(global_counts.sum())
        if gc_sum != n_observed:
            raise ValueError(
                f"calibration payload: sum(global_counts) = {gc_sum} != "
                f"n_observed = {n_observed}.{_ERR_SUFFIX}"
            )

        if n_total is not None:
            accounted = (
                n_observed + n_ex_mm + n_ex_ch + n_ex_ar + n_unobs
            )
            if accounted != int(n_total):
                raise ValueError(
                    "calibration payload: balance assertion failed: "
                    f"n_observed ({n_observed}) + n_excluded_multimap ({n_ex_mm}) "
                    f"+ n_excluded_chimera ({n_ex_ch}) + n_excluded_artifact "
                    f"({n_ex_ar}) + n_unobserved ({n_unobs}) = {accounted} "
                    f"!= n_total ({int(n_total)}).{_ERR_SUFFIX}"
                )

        return cls(
            global_counts=global_counts,
            per_region_counts=per_region_counts,
            fl_hist=fl_hist,
            u_left=u_left,
            u_right=u_right,
            n_observed=n_observed,
            n_excluded_multimap=n_ex_mm,
            n_excluded_chimera=n_ex_ch,
            n_excluded_artifact=n_ex_ar,
            n_unobserved=n_unobs,
            n_unannotated_ref=n_unann,
        )
