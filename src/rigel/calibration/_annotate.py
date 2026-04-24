"""Capture-class annotation for hybrid-capture panels.

Given a user-supplied BED3+ file of target intervals, attribute each
calibration region's mappable bp into two per-bp channels:

* ``E_on``  — mappable bp that overlap the (padded) target intervals.
* ``E_off`` — mappable bp that do not overlap the targets.

By construction ``E_on + E_off = mappable_bp`` for every region.  The
pair is consumed by the composite-Poisson calibration EM (see
``_em.py``) which fits independent ``λ_G_on`` and ``λ_G_off`` rates and
attributes each region's gDNA counts to the two channels via a
sum-of-Poissons superposition.

Design rationale: probes are typically ~120 bp while exonic /
intergenic calibration regions can span tens of kb.  A binary
on/off region label would force the whole region to share one rate,
diluting probe-level enrichment by the region-to-probe length ratio.
Per-bp attribution fixes this without introducing any thresholding.

See ``docs/calibration/capture_class_density_plan.md``.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


# A BED with fewer intervals than this is almost certainly a mistake
# (e.g. per-chromosome list).  Reject with a clear error.
_MIN_INTERVALS = 100

# Target-coverage cap: a BED covering > 50% of the genome is not a
# capture panel and would defeat the purpose of partitioning.
_MAX_COVERAGE_FRAC = 0.5


def _read_bed3(bed_path: Path, pad: int) -> dict[str, np.ndarray]:
    """Read a BED3+ file and return padded intervals grouped by ref."""
    if not bed_path.exists():
        raise FileNotFoundError(f"targets BED not found: {bed_path}")
    if pad < 0:
        raise ValueError(f"targets_pad must be >= 0, got {pad}")

    refs: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    opener = gzip.open if str(bed_path).endswith(".gz") else open
    with opener(bed_path, "rt") as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
                if len(parts) < 3:
                    continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as exc:
                raise ValueError(
                    f"BED {bed_path} line {line_no}: non-integer coordinates"
                ) from exc
            if end <= start:
                continue
            refs.append(parts[0])
            starts.append(max(0, start - pad))
            ends.append(end + pad)

    if len(refs) < _MIN_INTERVALS:
        raise ValueError(
            f"targets BED {bed_path} has {len(refs)} intervals "
            f"(< {_MIN_INTERVALS}); this does not look like a capture "
            f"panel.  Drop --targets for whole-genome / total-RNA data."
        )

    refs_arr = np.asarray(refs, dtype=object)
    starts_arr = np.asarray(starts, dtype=np.int64)
    ends_arr = np.asarray(ends, dtype=np.int64)

    grouped: dict[str, np.ndarray] = {}
    for ref in np.unique(refs_arr):
        mask = refs_arr == ref
        s = starts_arr[mask]
        e = ends_arr[mask]
        order = np.argsort(s)
        grouped[str(ref)] = np.column_stack((s[order], e[order])).astype(np.int64)
    return grouped


def _overlap_sum_sorted(
    region_starts: np.ndarray,
    region_ends: np.ndarray,
    bed_starts: np.ndarray,
    bed_ends: np.ndarray,
) -> np.ndarray:
    """Per-region overlap-bp against a set of BED intervals.

    Merges BED intervals into a disjoint union first so no region is
    double-counted when multiple probes overlap it.
    """
    if bed_starts.size == 0:
        return np.zeros(region_starts.shape[0], dtype=np.int64)

    order = np.argsort(bed_starts)
    bs = bed_starts[order]
    be = bed_ends[order]
    merged_s = [int(bs[0])]
    merged_e = [int(be[0])]
    for s, e in zip(bs[1:], be[1:]):
        s = int(s)
        e = int(e)
        if s <= merged_e[-1]:
            if e > merged_e[-1]:
                merged_e[-1] = e
        else:
            merged_s.append(s)
            merged_e.append(e)
    us = np.asarray(merged_s, dtype=np.int64)
    ue = np.asarray(merged_e, dtype=np.int64)
    out = np.zeros(region_starts.shape[0], dtype=np.int64)
    lo_idx = np.searchsorted(ue, region_starts, side="right")
    hi_idx = np.searchsorted(us, region_ends, side="left")
    for i in range(region_starts.shape[0]):
        rs = int(region_starts[i])
        re = int(region_ends[i])
        total = 0
        for j in range(int(lo_idx[i]), int(hi_idx[i])):
            ov_s = us[j] if us[j] > rs else rs
            ov_e = ue[j] if ue[j] < re else re
            if ov_e > ov_s:
                total += int(ov_e - ov_s)
        out[i] = total
    return out


def annotate_capture_class(
    region_df: pd.DataFrame,
    bed_path: Path | str,
    *,
    pad: int = 150,
    genome_size: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Attribute each region's mappable bp into on/off-target channels.

    Parameters
    ----------
    region_df
        Calibration region table (must contain ``ref``, ``start``,
        ``end`` columns; ``mappable_effective_length`` used if present,
        else ``length``).
    bed_path
        BED3+ file of capture targets.
    pad
        Symmetric padding (bp) applied to each BED interval before
        overlap scoring.
    genome_size
        Optional total genome size (bp) for the coverage-cap sanity
        check.  When ``None``, the cap is derived from ``region_df``.

    Returns
    -------
    (E_on, E_off) : tuple of np.ndarray[float64]
        Per-region mappable-bp attribution arrays with shape
        ``(len(region_df),)`` and ``E_on + E_off = mappable_bp``.  A
        region fully inside a probe has ``E_off == 0``; a region with
        no probe overlap has ``E_on == 0``.
    """
    bed_path = Path(bed_path)
    if not {"ref", "start", "end"}.issubset(region_df.columns):
        raise ValueError(
            "region_df must have columns ref, start, end for capture-class "
            "annotation"
        )

    grouped = _read_bed3(bed_path, pad=pad)

    total_padded_bp = 0
    for arr in grouped.values():
        if arr.size:
            total_padded_bp += int((arr[:, 1] - arr[:, 0]).sum())
    if genome_size is None:
        genome_size = int(
            (region_df["end"].astype(np.int64) - region_df["start"].astype(np.int64)).sum()
        )
    if genome_size > 0 and total_padded_bp > _MAX_COVERAGE_FRAC * genome_size:
        raise ValueError(
            f"targets BED {bed_path} covers {total_padded_bp:,} bp "
            f"(> {100 * _MAX_COVERAGE_FRAC:.0f}% of the reference space, "
            f"{genome_size:,} bp).  This is not a capture panel; drop "
            f"--targets for whole-genome inputs."
        )

    if "mappable_effective_length" in region_df.columns:
        mappable = region_df["mappable_effective_length"].to_numpy(dtype=np.float64)
    else:
        mappable = (
            region_df["end"].astype(np.int64) - region_df["start"].astype(np.int64)
        ).to_numpy(dtype=np.float64)

    refs = region_df["ref"].to_numpy()
    starts = region_df["start"].to_numpy(dtype=np.int64)
    ends = region_df["end"].to_numpy(dtype=np.int64)
    lengths = (ends - starts).astype(np.float64)

    n = len(region_df)
    overlap_bp = np.zeros(n, dtype=np.int64)
    unmatched_refs: set[str] = set()

    ref_order = pd.Series(refs).groupby(refs, sort=False).indices
    for ref, idx in ref_order.items():
        idx_arr = np.asarray(idx, dtype=np.int64)
        if ref not in grouped:
            unmatched_refs.add(str(ref))
            continue
        intervals = grouped[ref]
        overlap_bp[idx_arr] = _overlap_sum_sorted(
            starts[idx_arr], ends[idx_arr], intervals[:, 0], intervals[:, 1],
        )

    # Attribute mappable bp proportionally to the fraction of the
    # region spanned by the (padded) probe union.  Clamp to [0, 1] for
    # numerical safety — overlap_bp ≤ length always holds by
    # construction of the sorted-union overlap sum.
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(lengths > 0, overlap_bp.astype(np.float64) / lengths, 0.0)
    frac = np.clip(frac, 0.0, 1.0)
    E_on = mappable * frac
    E_off = mappable - E_on
    E_on = np.maximum(E_on, 0.0)
    E_off = np.maximum(E_off, 0.0)

    total_map = float(mappable.sum())
    total_on = float(E_on.sum())
    on_regions = int((E_on > 0).sum())
    logger.info(
        "[CAL] targets BED: %s | pad=%d | regions with probe overlap=%d/%d "
        "(%.1f%%) | ΣE_on/ΣE=%.3f",
        bed_path.name, pad, on_regions, n,
        100.0 * on_regions / max(n, 1),
        total_on / max(total_map, 1.0),
    )
    if unmatched_refs:
        sample = sorted(unmatched_refs)[:10]
        logger.warning(
            "[CAL] %d chromosomes in region_df have no matching intervals "
            "in %s (e.g. %s) — those regions default to all-off-target.  "
            "Check that BED and reference use matching chromosome names.",
            len(unmatched_refs), bed_path.name, ", ".join(sample),
        )

    return E_on.astype(np.float64, copy=False), E_off.astype(np.float64, copy=False)
