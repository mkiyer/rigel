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
junction is rejected as artifactual when **either** anchor in its CIGAR
is ≤ the blacklist maximum — the junction sits inside the
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
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


#: Canonical columns in the Rigel blacklist representation.
BLACKLIST_COLUMNS = (
    "ref", "start", "end", "max_anchor_left", "max_anchor_right",
)


def load_splice_blacklist_from_records(
    records: Iterable[Mapping[str, Any]],
    *,
    min_count: int = 2,
) -> pd.DataFrame:
    """Aggregate raw alignable splice-blacklist records into Rigel form.

    Parameters
    ----------
    records
        Iterable of dict-like rows as returned by
        ``alignable.AlignableStore.splice_blacklist()``.  Each row must
        provide ``chrom`` (str), ``intron_start`` (int), ``intron_end``
        (int), ``count`` (int), ``max_anchor_left`` (int), and
        ``max_anchor_right`` (int).  Other keys are ignored.
    min_count
        Only rows with ``count >= min_count`` enter the blacklist.
        Default ``2`` matches the historical alignable threshold.
        Use ``1`` to admit singletons; higher values keep only the
        most reproducible artifacts.

    Returns
    -------
    pandas.DataFrame
        One row per unique ``(ref, start, end)`` junction with columns
        :data:`BLACKLIST_COLUMNS`.  Anchors are aggregated across
        surviving read-length rows by ``max``.  Sorted by
        ``(ref, start, end)``.
    """
    if min_count < 1:
        raise ValueError(f"min_count must be >= 1, got {min_count}")

    refs: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    a_left: list[int] = []
    a_right: list[int] = []
    n_raw = 0
    n_below = 0
    for row in records:
        n_raw += 1
        c = int(row["count"])
        if c < min_count:
            n_below += 1
            continue
        refs.append(str(row["chrom"]))
        starts.append(int(row["intron_start"]))
        ends.append(int(row["intron_end"]))
        a_left.append(int(row["max_anchor_left"]))
        a_right.append(int(row["max_anchor_right"]))

    if not refs:
        logger.info(
            f"Splice blacklist: 0 junctions retained "
            f"({n_raw:,} raw, {n_below:,} below count={min_count})"
        )
        return _empty_blacklist_df()

    df = pd.DataFrame({
        "ref": refs,
        "start": np.asarray(starts, dtype=np.int32),
        "end": np.asarray(ends, dtype=np.int32),
        "max_anchor_left": np.asarray(a_left, dtype=np.int32),
        "max_anchor_right": np.asarray(a_right, dtype=np.int32),
    })

    agg = (
        df.groupby(["ref", "start", "end"], sort=False, observed=True)
          .agg(
              max_anchor_left=("max_anchor_left", "max"),
              max_anchor_right=("max_anchor_right", "max"),
          )
          .reset_index()
    )
    agg = agg.sort_values(["ref", "start", "end"], kind="stable").reset_index(drop=True)
    agg["start"] = agg["start"].astype(np.int32)
    agg["end"] = agg["end"].astype(np.int32)
    agg["max_anchor_left"] = agg["max_anchor_left"].astype(np.int32)
    agg["max_anchor_right"] = agg["max_anchor_right"].astype(np.int32)

    logger.info(
        f"Splice blacklist: {n_raw:,} raw rows → {len(df):,} kept "
        f"(count>={min_count}, dropped {n_below:,}) → "
        f"{len(agg):,} unique junctions"
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
        See :func:`load_splice_blacklist_from_records`.
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
    table = store.splice_blacklist_table()
    n_raw = table.num_rows

    if n_raw == 0:
        logger.info("Splice blacklist: store has empty blacklist table")
        return _empty_blacklist_df()

    # Arrow → pandas (dict-encoded strings decode automatically).  At
    # ~35 M rows × 8 cols (int32 / dict-string) this is ~1 GB transient
    # — acceptable for index builds and ~1000× faster than the
    # list-of-dicts path.
    df = table.to_pandas(types_mapper=None)
    counts = df["count"].to_numpy()
    keep_mask = counts >= min_count
    n_below = int((~keep_mask).sum())
    df = df.loc[keep_mask]
    if df.empty:
        logger.info(
            f"Splice blacklist: 0 junctions retained "
            f"({n_raw:,} raw, {n_below:,} below count={min_count})"
        )
        return _empty_blacklist_df()

    df = df.rename(columns={"chrom": "ref",
                            "intron_start": "start",
                            "intron_end": "end"})
    agg = (
        df.groupby(["ref", "start", "end"], sort=False, observed=True)
          .agg(max_anchor_left=("max_anchor_left", "max"),
               max_anchor_right=("max_anchor_right", "max"))
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
        f"{len(agg):,} unique junctions"
    )
    return agg


def _empty_blacklist_df() -> pd.DataFrame:
    return pd.DataFrame({
        "ref": pd.Series([], dtype="object"),
        "start": pd.Series([], dtype=np.int32),
        "end": pd.Series([], dtype=np.int32),
        "max_anchor_left": pd.Series([], dtype=np.int32),
        "max_anchor_right": pd.Series([], dtype=np.int32),
    })
