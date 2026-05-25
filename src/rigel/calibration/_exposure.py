"""FL-aware geometric exposure primitives for calibration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from ..frag_length_model import FragmentLengthModel
from ._arrays import RegionArrays


__all__ = [
    "_merged_blocks",
    "bp_weighted_mean_exposure_over_blocks",
    "contained_exposure_clipped",
    "fractional_boundary_side_exposure",
    "gdna_eff_len_for_loci",
    "l_eff_contained",
]


def l_eff_contained(
    spans_bp: np.ndarray,
    gdna_fl: FragmentLengthModel,
) -> np.ndarray:
    """FL-PMF-weighted contained effective length, in bp.

    Moved here from ``calibration.density_global`` in v4 Phase 0
    (``docs/fineregions/density_model_impl_plan_v4.md`` §4).
    """
    return gdna_fl.compute_all_transcript_eff_lens(
        np.asarray(spans_bp, dtype=np.int64),
        min_value=0.0,
    )


def _merged_blocks(
    blocks: Sequence[tuple[int, int, int]],
) -> list[tuple[int, int, int]]:
    """Return non-overlapping ``(ref_id, start, end)`` blocks."""
    valid_blocks = [
        (int(ref_id), int(start), int(end))
        for ref_id, start, end in blocks
        if int(end) > int(start)
    ]
    if not valid_blocks:
        return []
    valid_blocks.sort()
    merged: list[tuple[int, int, int]] = []
    cur_ref, cur_start, cur_end = valid_blocks[0]
    for ref_id, start, end in valid_blocks[1:]:
        if ref_id == cur_ref and start <= cur_end:
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_ref, cur_start, cur_end))
            cur_ref, cur_start, cur_end = ref_id, start, end
    merged.append((cur_ref, cur_start, cur_end))
    return merged


def bp_weighted_mean_exposure_over_blocks(
    blocks: Sequence[tuple[int, int, int]],
    region_arrays: RegionArrays,
    exposure,
    *,
    min_weight: float = 1.0e-4,
) -> float:
    """Return bp-weighted mean ``A_r`` over merged genomic blocks.

    ``exposure`` is intentionally structural here: Phase 7 supplies a
    ``RegionExposure`` object with ``mode`` and ``A_r`` attributes. Keeping
    this helper independent of that future class lets Phase 1 retain only the
    stable geometry code.

    The returned weight is the bp-weighted mean of ``exposure.A_r`` over the
    merged blocks, floored at ``min_weight``. It is **not** clipped above 1:
    density-derived ``A_r`` may exceed 1, and clipping would silently
    truncate high-exposure regions.
    """
    if min_weight < 0.0:
        raise ValueError(
            "bp_weighted_mean_exposure_over_blocks: "
            f"min_weight must be >= 0, got {min_weight}."
        )

    merged = _merged_blocks(blocks)
    if not merged:
        return 1.0

    raw_bp = float(sum(end - start for _, start, end in merged))
    if raw_bp <= 0.0:
        return 1.0
    if getattr(exposure, "mode", "uniform") == "uniform":
        return 1.0

    # A_r may exceed 1 when constructed from density evidence.
    weights = np.asarray(getattr(exposure, "A_r"), dtype=np.float64)
    if weights.shape != region_arrays.start.shape:
        raise ValueError(
            "bp_weighted_mean_exposure_over_blocks: exposure.A_r shape "
            f"{weights.shape} != region shape {region_arrays.start.shape}."
        )

    weighted_bp = 0.0
    for ref_id, block_start, block_end in merged:
        if ref_id < 0 or ref_id >= region_arrays.n_refs:
            continue
        lo = int(region_arrays.ref_offsets[ref_id])
        hi = int(region_arrays.ref_offsets[ref_id + 1])
        starts = region_arrays.start[lo:hi]
        ends = region_arrays.end[lo:hi]
        local_lo = int(np.searchsorted(ends, block_start, side="right"))
        local_hi = int(np.searchsorted(starts, block_end, side="left"))
        for local_idx in range(local_lo, local_hi):
            idx = lo + local_idx
            overlap = min(int(region_arrays.end[idx]), block_end) - max(
                int(region_arrays.start[idx]),
                block_start,
            )
            if overlap > 0:
                weighted_bp += float(overlap) * float(weights[idx])

    return float(max(weighted_bp / raw_bp, min_weight))


def _ref_length_for_locus(
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    locus,
) -> int:
    """Return reference length for a Locus-like object."""
    if isinstance(ref_lengths, Mapping):
        if locus.ref in ref_lengths:
            return int(ref_lengths[locus.ref])
        if locus.ref_id in ref_lengths:
            return int(ref_lengths[locus.ref_id])
        raise KeyError(f"Reference length missing for ref={locus.ref!r}, ref_id={locus.ref_id!r}.")

    ref_id = int(locus.ref_id)
    if ref_id < 0 or ref_id >= len(ref_lengths):
        raise KeyError(f"Reference id {ref_id} out of bounds for ref_lengths.")
    return int(ref_lengths[ref_id])


def gdna_eff_len_for_loci(
    loci: tuple | list,
    ref_lengths: Mapping[str | int, int] | Sequence[int],
    fl: FragmentLengthModel,
    *,
    min_value: float = 1.0,
) -> float:
    """FL-marginal gDNA overlap effective length for a ``MultiLocus``."""
    if min_value < 0.0:
        raise ValueError(f"gdna_eff_len_for_loci: min_value must be >= 0, got {min_value}.")
    if not loci:
        return float(min_value)

    pmf = fl.pmf
    positive_ell = np.flatnonzero(pmf[1:] > 0.0) + 1
    if positive_ell.size == 0:
        return float(min_value)

    if len(loci) == 1:
        locus = loci[0]
        span = max(int(locus.end) - int(locus.start), 0)
        ref_len = _ref_length_for_locus(ref_lengths, locus)
        max_ell = int(positive_ell[-1])
        if int(locus.start) - max_ell + 1 >= 0 and int(locus.end) <= ref_len - max_ell + 1:
            probs = pmf[positive_ell]
            positive_mass = float(probs.sum())
            positive_mean = float(np.dot(positive_ell.astype(np.float64), probs))
            eff = float(span) * positive_mass + positive_mean - positive_mass
            return max(eff, float(min_value))

    total = 0.0
    for ell in positive_ell:
        ell_i = int(ell)
        by_ref: dict[int | str, list[tuple[int, int]]] = {}
        for locus in loci:
            ref_len = _ref_length_for_locus(ref_lengths, locus)
            valid_hi = ref_len - ell_i + 1
            if valid_hi <= 0:
                continue
            lo = max(int(locus.start) - ell_i + 1, 0)
            hi = min(int(locus.end), valid_hi)
            if hi <= lo:
                continue
            ref_key = int(locus.ref_id) if hasattr(locus, "ref_id") else locus.ref
            by_ref.setdefault(ref_key, []).append((lo, hi))

        n_valid = 0
        for windows in by_ref.values():
            windows.sort()
            cur_lo, cur_hi = windows[0]
            for lo, hi in windows[1:]:
                if lo <= cur_hi:
                    if hi > cur_hi:
                        cur_hi = hi
                else:
                    n_valid += cur_hi - cur_lo
                    cur_lo, cur_hi = lo, hi
            n_valid += cur_hi - cur_lo

        total += float(pmf[ell_i]) * float(n_valid)

    return max(total, float(min_value))


def contained_exposure_clipped(
    starts: np.ndarray,
    ends: np.ndarray,
    clip_lo: int,
    clip_hi: int,
    fl: FragmentLengthModel,
) -> tuple[np.ndarray, np.ndarray]:
    """FL-PMF-weighted contained-fragment effective length, full and clipped."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if starts.shape != ends.shape:
        raise ValueError(
            f"contained_exposure_clipped: starts.shape {starts.shape} != ends.shape {ends.shape}"
        )
    spans_full = np.maximum(ends - starts, 0)
    clip_starts = np.maximum(starts, np.int64(clip_lo))
    clip_ends = np.minimum(ends, np.int64(clip_hi))
    spans_clip = np.maximum(clip_ends - clip_starts, 0)
    eff_full = fl.compute_all_transcript_eff_lens(spans_full, min_value=0.0)
    eff_clip = fl.compute_all_transcript_eff_lens(spans_clip, min_value=0.0)
    return eff_full, eff_clip


def fractional_boundary_side_exposure(
    lengths_bp: np.ndarray,
    gdna_fl: FragmentLengthModel,
) -> np.ndarray:
    """FL-PMF-weighted per-side boundary exposure for region lengths."""
    lengths = np.maximum(np.asarray(lengths_bp, dtype=np.int64), 0)
    pmf = np.asarray(gdna_fl.pmf, dtype=np.float64)
    if pmf.size == 0:
        return np.zeros(lengths.shape, dtype=np.float64)

    ell_minus_one = np.maximum(np.arange(pmf.size, dtype=np.float64) - 1.0, 0.0)
    positive_mass = pmf.copy()
    positive_mass[0] = 0.0
    cdf = np.cumsum(positive_mass)
    first_moment = np.cumsum(positive_mass * ell_minus_one)
    total_mass = float(cdf[-1])

    threshold = np.minimum(lengths + 1, pmf.size - 1)
    low_moment = first_moment[threshold]
    tail_mass = total_mass - cdf[threshold]
    return 0.5 * (low_moment + lengths.astype(np.float64, copy=False) * tail_mass)
