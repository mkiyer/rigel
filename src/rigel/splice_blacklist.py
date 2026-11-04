"""rigel.splice_blacklist — read the splice-artifact blacklist from an alignable store.

The sister tool ``alignable`` tiles the genome with gDNA-like fragments, re-aligns them, and records every spliced
alignment the aligner produces: junctions where it emits an ``N`` operation from plain genomic DNA. Each artifact is
keyed by ``(ref, intron_start, intron_end)``, with the largest left and right anchors (reference-advancing CIGAR bases)
seen at each read length. The C++ scanner (``scan/bam_record.cpp``) rejects a junction as an artifact when either of its
anchors is at most the blacklist's.

This module reads the store (a directory or a packaged ``.zarr.zip``), filters rows by count, and collapses the
per-read-length rows into one conservative anchor envelope per junction for the index. The store's table has one row
per ``(chrom, intron, strand, read_length)``::

    chrom, intron_start, intron_end, strand, read_length, count, max_anchor_left, max_anchor_right

``strand`` is always ``'.'`` (the detection is strand-agnostic), so rows collapse on ``(ref, start, end)``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def aggregate_splice_blacklist(rows: pd.DataFrame, *, min_count: int) -> pd.DataFrame:
    """Aggregate raw alignable splice-blacklist rows into Rigel form.

    Parameters
    ----------
    rows
        One row per ``(chrom, intron, strand, read_length)`` as the alignable store holds them: the
        columns ``chrom``, ``intron_start``, ``intron_end``, ``count``, ``max_anchor_left`` and
        ``max_anchor_right``; any others are ignored.
    min_count
        Only rows with ``count >= min_count`` enter the blacklist;
        ``IndexConfig.splice_blacklist_min_count`` sets it and says why.

    Returns
    -------
    pandas.DataFrame
        One row per unique ``(ref, start, end)`` sj with columns ``ref``, ``start``, ``end``,
        ``max_anchor_left``, ``max_anchor_right`` and ``count``, the anchors aggregated across the
        surviving read-length rows by ``max`` and their counts summed, sorted by ``(ref, start, end)``.
        The resolver reads the anchors; the count is kept so that a threshold can be read off the
        index.
    """
    if min_count < 1:
        raise ValueError(f"min_count must be >= 1, got {min_count}")
    n_raw = len(rows)
    if n_raw == 0:
        logger.info("Splice blacklist: the store's blacklist table is empty")
        return _empty_blacklist_df()

    keep_mask = rows["count"].to_numpy() >= min_count
    n_below = int((~keep_mask).sum())
    kept = rows.loc[keep_mask]
    if kept.empty:
        logger.info(
            f"Splice blacklist: 0 sj retained ({n_raw:,} raw, {n_below:,} below count={min_count})"
        )
        return _empty_blacklist_df()

    kept = kept.rename(columns={"chrom": "ref", "intron_start": "start", "intron_end": "end"})
    agg = (
        kept.groupby(["ref", "start", "end"], sort=False, observed=True)
        .agg(
            max_anchor_left=("max_anchor_left", "max"),
            max_anchor_right=("max_anchor_right", "max"),
            count=("count", "sum"),
        )
        .reset_index()
    )
    agg = agg.sort_values(["ref", "start", "end"], kind="stable").reset_index(drop=True)
    agg["ref"] = agg["ref"].astype(object)
    agg["start"] = agg["start"].astype(np.int32)
    agg["end"] = agg["end"].astype(np.int32)
    agg["max_anchor_left"] = agg["max_anchor_left"].astype(np.int32)
    agg["max_anchor_right"] = agg["max_anchor_right"].astype(np.int32)
    agg["count"] = agg["count"].astype(np.int64)

    logger.info(
        f"Splice blacklist: {n_raw:,} raw rows → {int(keep_mask.sum()):,} kept "
        f"(count>={min_count}, dropped {n_below:,}) → "
        f"{len(agg):,} unique sj"
    )
    return agg


def load_splice_blacklist_from_zarr(
    store_path: str | Path,
    *,
    min_count: int,
) -> pd.DataFrame:
    """Load and aggregate the splice-artifact blacklist from an alignable store.

    Opens the store with :func:`alignable.open` and aggregates its blacklist, a zstd-compressed Feather v2 file at
    ``<store>/mappability.zarr/splice_blacklist.feather``.

    Parameters
    ----------
    store_path
        Path to an alignable output directory or a ``.zarr.zip`` file.
    min_count
        See :func:`aggregate_splice_blacklist`.
    """
    if min_count < 1:
        raise ValueError(f"min_count must be >= 1, got {min_count}")
    if not Path(store_path).exists():
        raise FileNotFoundError(f"No alignable store at {store_path}")

    try:
        import alignable
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError(
            "Reading the splice blacklist from an alignable store requires "
            "the 'alignable' package. Install it into the rigel environment, "
            "or rebuild the index with --no-mappability."
        ) from exc

    logger.info(f"Loading splice blacklist from alignable store: {store_path}")
    store = alignable.open(str(store_path))
    # Arrow → pandas; dict-encoded strings decode automatically.
    return aggregate_splice_blacklist(
        store.splice_blacklist_table().to_pandas(types_mapper=None), min_count=min_count
    )


def _empty_blacklist_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ref": pd.Series([], dtype="object"),
            "start": pd.Series([], dtype=np.int32),
            "end": pd.Series([], dtype=np.int32),
            "max_anchor_left": pd.Series([], dtype=np.int32),
            "max_anchor_right": pd.Series([], dtype=np.int32),
            "count": pd.Series([], dtype=np.int64),
        }
    )
