"""FL-aware geometric exposure primitives for calibration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from ._arrays import RegionArrays
from ..frag_length_model import FragmentLengthModel


MIN_EXPOSURE_FACTOR: float = 1.0e-4
REGION_EXPOSURE_ZERO_WIDTH: np.uint16 = np.uint16(0x1)
REGION_EXPOSURE_UNCOVERED: np.uint16 = np.uint16(0x2)
REGION_EXPOSURE_FLOORED: np.uint16 = np.uint16(0x4)


__all__ = [
    "MIN_EXPOSURE_FACTOR",
    "REGION_EXPOSURE_FLOORED",
    "REGION_EXPOSURE_UNCOVERED",
    "REGION_EXPOSURE_ZERO_WIDTH",
    "RegionWeightedExposure",
    "_merged_blocks",
    "component_bp_weighted_exposure",
    "contained_exposure_clipped",
    "fractional_boundary_side_exposure",
    "gdna_eff_len_for_loci",
    "l_eff_contained",
]


@dataclass(frozen=True, slots=True)
class RegionWeightedExposure:
    """Bp-weighted component exposure over calibration regions."""

    exposure_factor: np.ndarray
    covered_bp: np.ndarray
    weighted_bp: np.ndarray
    flags: np.ndarray


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


def component_bp_weighted_exposure(
    *,
    block_ref_ids: np.ndarray,
    block_starts: np.ndarray,
    block_ends: np.ndarray,
    component_ids: np.ndarray,
    n_components: int,
    region_arrays: RegionArrays,
    omega: np.ndarray,
    strict_coverage: bool = True,
    min_exposure_factor: float = MIN_EXPOSURE_FACTOR,
) -> RegionWeightedExposure:
    """Compute bp-weighted exposure factors for arbitrary component blocks.

    ``omega`` is supplied in region-table order. ``RegionArrays`` stores sorted
    region geometry, so this helper reorders ``omega`` with ``region_arrays.order``
    before sweeping overlaps.
    """
    n_comp = int(n_components)
    if n_comp < 0:
        raise ValueError(f"n_components must be non-negative, got {n_components!r}.")

    min_factor = float(min_exposure_factor)
    if not np.isfinite(min_factor) or min_factor <= 0.0:
        raise ValueError(
            f"min_exposure_factor must be finite and positive, got {min_exposure_factor!r}."
        )

    ref_ids = np.asarray(block_ref_ids, dtype=np.int32)
    starts = np.asarray(block_starts, dtype=np.int64)
    ends = np.asarray(block_ends, dtype=np.int64)
    comp_ids = np.asarray(component_ids, dtype=np.int64)
    if not (ref_ids.shape == starts.shape == ends.shape == comp_ids.shape):
        raise ValueError("block_ref_ids, block_starts, block_ends, and component_ids must align.")

    omega_arr = np.asarray(omega, dtype=np.float64)
    if omega_arr.shape != region_arrays.start.shape:
        raise ValueError(
            f"omega shape {omega_arr.shape} does not match regions {region_arrays.start.shape}."
        )
    if not np.all(np.isfinite(omega_arr)) or np.any(omega_arr <= 0.0):
        raise ValueError("omega must contain only finite positive values.")
    omega_sorted_raw = omega_arr[region_arrays.order]
    omega_floored = omega_sorted_raw < min_factor
    omega_sorted = np.maximum(omega_sorted_raw, min_factor)

    covered_bp = np.zeros(n_comp, dtype=np.float64)
    weighted_bp = np.zeros(n_comp, dtype=np.float64)
    total_bp = np.zeros(n_comp, dtype=np.float64)
    uncovered_bp = np.zeros(n_comp, dtype=np.float64)
    flags = np.zeros(n_comp, dtype=np.uint16)

    if ref_ids.size == 0 or n_comp == 0:
        return RegionWeightedExposure(
            exposure_factor=np.ones(n_comp, dtype=np.float64),
            covered_bp=covered_bp,
            weighted_bp=weighted_bp,
            flags=flags,
        )

    valid_width = ends - starts
    if np.any(valid_width <= 0):
        raise ValueError("component blocks must have end > start.")
    if np.any(ref_ids < 0) or np.any(ref_ids >= int(region_arrays.n_refs)):
        raise ValueError("component block ref ids are out of bounds for RegionArrays.")
    if np.any(comp_ids < 0) or np.any(comp_ids >= n_comp):
        raise ValueError("component_ids are out of bounds for n_components.")

    order = np.lexsort((comp_ids, ends, starts, ref_ids))
    ref_ids = ref_ids[order]
    starts = starts[order]
    ends = ends[order]
    comp_ids = comp_ids[order]

    for ref_id, block_start, block_end, comp_id in zip(ref_ids, starts, ends, comp_ids):
        ref = int(ref_id)
        comp = int(comp_id)
        start = int(block_start)
        end = int(block_end)
        width = float(end - start)
        total_bp[comp] += width

        ref_lo = int(region_arrays.ref_offsets[ref])
        ref_hi = int(region_arrays.ref_offsets[ref + 1])
        if ref_hi <= ref_lo:
            uncovered_bp[comp] += width
            continue

        local_starts = region_arrays.start[ref_lo:ref_hi]
        local_ends = region_arrays.end[ref_lo:ref_hi]
        first = ref_lo + int(np.searchsorted(local_ends, start, side="right"))
        last = ref_lo + int(np.searchsorted(local_starts, end, side="left"))
        if last <= first:
            uncovered_bp[comp] += width
            continue

        region_starts = region_arrays.start[first:last]
        region_ends = region_arrays.end[first:last]
        overlap = np.minimum(region_ends, end) - np.maximum(region_starts, start)
        overlap = np.maximum(overlap, 0)
        covered = float(np.sum(overlap, dtype=np.float64))
        if covered <= 0.0:
            uncovered_bp[comp] += width
            continue

        covered_bp[comp] += covered
        weighted_bp[comp] += float(
            np.dot(overlap.astype(np.float64, copy=False), omega_sorted[first:last])
        )
        if np.any((overlap > 0) & omega_floored[first:last]):
            flags[comp] |= REGION_EXPOSURE_FLOORED
        if covered < width:
            uncovered_bp[comp] += width - covered

    uncovered = uncovered_bp > 1.0e-9
    if np.any(uncovered):
        flags[uncovered] |= REGION_EXPOSURE_UNCOVERED
        if strict_coverage:
            n_bad = int(np.count_nonzero(uncovered))
            worst = float(np.max(uncovered_bp[uncovered]))
            raise ValueError(
                "component_bp_weighted_exposure: region table does not cover "
                f"{n_bad} component(s); max uncovered bp={worst:.0f}."
            )
        covered_bp += uncovered_bp
        weighted_bp += uncovered_bp

    zero = total_bp <= 0.0
    flags[zero] |= REGION_EXPOSURE_ZERO_WIDTH
    factor = np.divide(
        weighted_bp,
        covered_bp,
        out=np.ones(n_comp, dtype=np.float64),
        where=covered_bp > 0.0,
    )
    floored = factor < min_factor
    flags[floored] |= REGION_EXPOSURE_FLOORED
    factor = np.maximum(factor, min_factor)
    if not np.all(np.isfinite(factor)):
        raise ValueError("component exposure factors must be finite after aggregation.")

    return RegionWeightedExposure(
        exposure_factor=factor.astype(np.float64, copy=False),
        covered_bp=covered_bp.astype(np.float64, copy=False),
        weighted_bp=weighted_bp.astype(np.float64, copy=False),
        flags=flags,
    )


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
