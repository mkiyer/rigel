"""
rigel.mappability — Per-region effective length from an alignable Zarr.

For every calibration region ``[start, end)`` on reference ``ref``, we
need a single scalar that captures *how much* of that region is
unambiguously coverable by an aligner at a given read length:

    mappable_effective_length(region) = sum_{p in region} f(p)

where ``f(p) = n_frag_ontarget / nh_ontarget`` is the per-base
fractional mappability (1/NH credit on-target) reported by alignable.
This is the correct Poisson exposure for a 1/NH-credit EM and is the
only mappability signal that downstream calibration consumes.

This module is a thin orchestration layer:

* It opens the alignable Zarr store via the documented
  :class:`alignable.api.AlignableStore` API.
* It iterates regions in ``(ref, start)`` order so that adjacent
  queries land on the same Zarr chunk and benefit from the LRU cache.
* It records the alignable provenance (aligner, version, read-length
  bin, store digest) so the index manifest can pin the calibration
  invariants.

The computation is intentionally simple — a single ``sum`` over a
float32 slice per region.  The hot path lives inside alignable
(zarr v3 zstd decompression).  On GRCh38/STAR with ~684k regions the
full sweep takes ~12 minutes and uses ~3GB RAM (cache_size=4).
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MappabilityProvenance:
    """Tag describing the alignable artifact a Rigel index was built from."""

    store_path: str
    aligner: str
    aligner_version: str | None
    read_length: int
    frag_len_mean: int
    frag_len_sd: int
    error_rate: float
    status: str

    def to_dict(self) -> dict:
        return {
            "store_path": self.store_path,
            "aligner": self.aligner,
            "aligner_version": self.aligner_version,
            "read_length": int(self.read_length),
            "frag_len_mean": int(self.frag_len_mean),
            "frag_len_sd": int(self.frag_len_sd),
            "error_rate": float(self.error_rate),
            "status": self.status,
        }


def compute_region_exposures(
    store_path: str | Path,
    region_df: pd.DataFrame,
    *,
    read_length: int = 100,
    cache_size: int = 4,
) -> tuple[np.ndarray, MappabilityProvenance]:
    """Compute per-region mappable effective length from an alignable store.

    Parameters
    ----------
    store_path
        Path to an alignable Zarr store.
    region_df
        Calibration-region table.  Must have columns ``ref``, ``start``,
        ``end`` and a stable row order (the returned array is aligned
        to that row order via ``region_id`` semantics — i.e. ``out[i]``
        is the exposure for the region in row ``i``).
    read_length
        Read-length bin to query.  Must be present in
        ``store.read_length_bins``.
    cache_size
        Number of decompressed Zarr chunks to retain in alignable's
        LRU cache.  Four is enough for a strictly increasing scan;
        increase if the index includes wildly out-of-order chromosome
        groupings.

    Returns
    -------
    exposures : np.ndarray (float32, shape ``(len(region_df),)``)
        ``mappable_effective_length`` per region.  Regions on
        chromosomes the alignable store does not cover get exposure 0.
    provenance : MappabilityProvenance
    """
    try:
        import alignable
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError(
            "Computing mappable effective length requires the 'alignable' "
            "package.  Install it into the rigel environment, or rebuild "
            "the index with --no-mappability."
        ) from exc

    store_path = Path(store_path)
    logger.info(f"[MAP] Opening alignable store: {store_path}")
    store = alignable.open(str(store_path), cache_size=cache_size)

    if read_length not in store.read_length_bins:
        raise ValueError(
            f"read_length={read_length} not in alignable bins "
            f"{store.read_length_bins} (store: {store_path})"
        )
    if store.status != "complete":
        logger.warning(
            f"[MAP] alignable store status is {store.status!r}, expected 'complete'"
        )

    provenance = MappabilityProvenance(
        store_path=str(store_path),
        aligner=str(store.aligner),
        aligner_version=(
            None if store.aligner_version is None else str(store.aligner_version)
        ),
        read_length=int(read_length),
        frag_len_mean=int(store.frag_len_mean),
        frag_len_sd=int(store.frag_len_sd),
        error_rate=float(store.error_rate),
        status=str(store.status),
    )

    n = len(region_df)
    exposures = np.zeros(n, dtype=np.float32)
    if n == 0:
        return exposures, provenance

    refs = region_df["ref"].astype(str).to_numpy()
    starts = region_df["start"].to_numpy(dtype=np.int64)
    ends = region_df["end"].to_numpy(dtype=np.int64)

    chrom_lengths = store.chromosomes
    chrom_set = set(chrom_lengths.keys())

    # Sort by (ref, start) for cache locality, but write back in original order.
    order = np.lexsort((starts, refs))
    n_skipped = 0
    n_clamped = 0
    t0 = time.perf_counter()
    last_log = t0
    last_chrom: str | None = None

    for k, idx in enumerate(order):
        ref = refs[idx]
        if ref not in chrom_set:
            n_skipped += 1
            continue
        s = int(starts[idx])
        e = int(ends[idx])
        if s >= e:
            continue
        clen = chrom_lengths[ref]
        if e > clen or s < 0:
            n_clamped += 1
            s = max(0, s)
            e = min(e, clen)
            if s >= e:
                continue
        fm = store.fractional_mappability(ref, s, e, read_length)
        exposures[idx] = float(fm.sum())

        # Lightweight progress log per chromosome change + every 60s.
        now = time.perf_counter()
        if ref != last_chrom:
            last_chrom = ref
        if now - last_log > 60.0:
            logger.info(
                f"[MAP] processed {k + 1:,} / {n:,} regions "
                f"({(k + 1) / max(now - t0, 1e-9):,.0f}/s, ref={ref})"
            )
            last_log = now

    dt = time.perf_counter() - t0
    cache = store.cache_info()
    logger.info(
        f"[MAP] done: {n:,} regions in {dt:.1f}s "
        f"({n / max(dt, 1e-9):,.0f}/s); skipped={n_skipped}, clamped={n_clamped}; "
        f"cache hits={getattr(cache, 'hits', '?')} misses={getattr(cache, 'misses', '?')}"
    )

    total_span = float((ends - starts).clip(min=0).sum())
    total_eff = float(exposures.sum())
    if total_span > 0:
        logger.info(
            f"[MAP] coverage: {total_eff:,.0f} / {total_span:,.0f} bp "
            f"({100.0 * total_eff / total_span:.1f}%) at read_length={read_length}"
        )

    return exposures, provenance


def uniform_region_exposures(region_df: pd.DataFrame) -> np.ndarray:
    """Return ``mappable_effective_length = length`` for the --no-mappability path.

    Equivalent to assuming ``f(p) = 1`` everywhere — appropriate for
    synthetic genomes, stranded-only benchmarks, and other settings
    where mappability variation is known to be irrelevant.
    """
    n = len(region_df)
    if n == 0:
        return np.zeros(0, dtype=np.float32)
    length = (
        region_df["end"].to_numpy(dtype=np.int64)
        - region_df["start"].to_numpy(dtype=np.int64)
    )
    return length.clip(min=0).astype(np.float32)
