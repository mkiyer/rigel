"""
rigel.splice_blacklist — Splice-junction artifact blacklist ingestion.

The sister tool ``alignable`` tiles the reference genome with synthetic
gDNA-like fragments, re-aligns them, and records every spliced alignment
the chosen aligner produces.  By construction these are false-positive
splice junctions: coordinates where an aligner (minimap2, STAR, ...) is
known to emit an ``N`` CIGAR op from plain genomic DNA.

Each artifact is keyed by ``(ref, intron_start, intron_end)`` and
characterised by the maximum *left* and *right* anchor (in
reference-advancing CIGAR bases) observed across all false-positive
alignments at each read length.  At BAM-scan time, a fragment's splice
sj is rejected as artifactual when **either** anchor in its CIGAR
is ≤ the blacklist maximum — the sj sits inside the
"plausible-from-gDNA" envelope.

This module handles only blacklist *ingestion* from an alignable store
(unpackaged directory or packaged ``.zarr.zip``): filtering by
per-row count, aggregating per-read-length rows into a single
conservative anchor envelope, and producing a tidy DataFrame for the
Rigel index.  The per-fragment anchor check itself lives in the C++
BAM scanner (``bam_scanner.cpp``).

Alignable blacklist table schema (one row per
``(chrom, intron, strand, read_length)``)::

    chrom, intron_start, intron_end, strand, read_length,
    count, max_anchor_left, max_anchor_right

The ``strand`` column is always ``'.'`` — splice artifacts are strand
agnostic by nature of the detection scheme.  We collapse on
``(ref, start, end)`` only.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


#: Canonical columns in the Rigel blacklist representation.
BLACKLIST_COLUMNS = (
    "ref",
    "start",
    "end",
    "max_anchor_left",
    "max_anchor_right",
)


def aggregate_splice_blacklist(rows: pd.DataFrame, *, min_count: int = 2) -> pd.DataFrame:
    """Aggregate raw alignable splice-blacklist rows into Rigel form.

    Parameters
    ----------
    rows
        One row per ``(chrom, intron, strand, read_length)`` as the alignable store holds them: the
        columns ``chrom``, ``intron_start``, ``intron_end``, ``count``, ``max_anchor_left`` and
        ``max_anchor_right``; any others are ignored.
    min_count
        Only rows with ``count >= min_count`` enter the blacklist. Default ``2`` matches the
        alignable threshold; ``1`` admits singletons, higher values keep only the most reproducible
        artifacts.

    Returns
    -------
    pandas.DataFrame
        One row per unique ``(ref, start, end)`` sj with columns :data:`BLACKLIST_COLUMNS`, the
        anchors aggregated across the surviving read-length rows by ``max``, sorted by
        ``(ref, start, end)``.
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
            max_anchor_left=("max_anchor_left", "max"), max_anchor_right=("max_anchor_right", "max")
        )
        .reset_index()
    )
    agg = agg.sort_values(["ref", "start", "end"], kind="stable").reset_index(drop=True)
    agg["ref"] = agg["ref"].astype(object)
    agg["start"] = agg["start"].astype(np.int32)
    agg["end"] = agg["end"].astype(np.int32)
    agg["max_anchor_left"] = agg["max_anchor_left"].astype(np.int32)
    agg["max_anchor_right"] = agg["max_anchor_right"].astype(np.int32)

    logger.info(
        f"Splice blacklist: {n_raw:,} raw rows → {int(keep_mask.sum()):,} kept "
        f"(count>={min_count}, dropped {n_below:,}) → "
        f"{len(agg):,} unique sj"
    )
    return agg


def load_splice_blacklist_from_zarr(
    store_path: str | Path,
    *,
    min_count: int = 2,
) -> pd.DataFrame:
    """Load and aggregate the splice-artifact blacklist from an alignable store.

    Opens the alignable store (unpackaged directory *or* packaged
    ``.zarr.zip``) via :func:`alignable.open`, pulls the blacklist as a
    zero-copy :class:`pyarrow.Table`, and aggregates it vectorised.
    Per the alignable v0.1+ layout the blacklist is a single
    zstd-compressed Feather v2 file embedded in the Zarr store at
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
    # Arrow → pandas (dict-encoded strings decode automatically). At ~35 M rows × 8 columns this is
    # ~1 GB transient, acceptable for an index build.
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
        }
    )
