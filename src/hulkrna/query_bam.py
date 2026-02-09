"""
hulkrna.query_bam — BAM scanning and reference index querying (Pass 1).

Parses a coordinate-sorted BAM file, builds ``Fragment`` objects from
read pairs, queries the ``HulkIndex`` for overlapping reference intervals
and exact-match splice junctions, and streams all hits to an intermediate
Feather (Arrow IPC) file in fixed-size chunks so that memory usage stays
bounded regardless of BAM size.

This module is purely a data-generation step: no strand model training,
no insert size learning, no fragment classification, and no count
assignment. All interpretation is deferred to Pass 2.

Intermediate file schema
------------------------
One row per hit (fragment × reference interval match):

============== ====== =================================================
Column         Type   Description
============== ====== =================================================
frag_id        int32  Groups rows belonging to the same fragment
ref            str    Chromosome of the query interval (exon block or SJ)
start          int32  Genomic start (0-based) of the query interval
end            int32  Genomic end (half-open) of the query interval
strand         int8   Genomic strand of the query interval (Strand enum)
interval_type  int8   EXON=0, INTRON=1, INTERGENIC=2, SJ=3
overlap        int32  Overlap bases (exon/intron/intergenic) or 0 (SJ)
t_index        int32  Transcript index in reference
g_index        int32  Gene index in reference
============== ====== =================================================

Reconstructable from hits (not stored):
- ``combined_strand`` — bitwise OR of all ``strand`` values per frag_id
- ``is_spliced`` — any row with ``interval_type == SJ``
- ``is_chimeric`` — multiple ``ref`` values or ambiguous strand pattern
- ``insert_size`` — ``max(end) - min(start)`` per frag_id (same ref)
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa

from .bam import parse_bam_file
from .core import IntervalType, Strand
from .fragment import Fragment
from .index import HulkIndex

logger = logging.getLogger(__name__)

# Column order for the intermediate Feather file
HIT_COLUMNS = [
    "frag_id",
    "ref",
    "start",
    "end",
    "strand",
    "interval_type",
    "overlap",
    "t_index",
    "g_index",
]

# PyArrow schema for streaming writes (Feather v2 = Arrow IPC format)
HIT_SCHEMA = pa.schema([
    ("frag_id", pa.int32()),
    ("ref", pa.string()),
    ("start", pa.int32()),
    ("end", pa.int32()),
    ("strand", pa.int8()),
    ("interval_type", pa.int8()),
    ("overlap", pa.int32()),
    ("t_index", pa.int32()),
    ("g_index", pa.int32()),
])

# Default number of hit rows to buffer before flushing to disk
DEFAULT_CHUNK_SIZE = 500_000


def _query_fragment(
    frag: Fragment,
    frag_id: int,
    index: HulkIndex,
) -> list[tuple]:
    """Query the index for a single fragment and return raw hit tuples.

    For each exon block in the fragment, queries the unified cgranges
    index for overlapping EXON/INTRON/INTERGENIC intervals.

    For each intron (splice junction) in the fragment, performs an
    exact-match lookup against the SJ map in the index.

    Parameters
    ----------
    frag : Fragment
        The consolidated fragment from paired-end reads.
    frag_id : int
        Sequential fragment identifier.
    index : HulkIndex
        The loaded reference index.

    Returns
    -------
    list[tuple]
        Each tuple is ``(frag_id, ref, start, end, strand,
        interval_type, overlap, t_index, g_index)``.
    """
    rows: list[tuple] = []

    # Query each exon block against the unified cgranges index
    for exon in frag.exons:
        for t_idx, g_idx, itype, ovl in index.query_exon(exon):
            rows.append((
                frag_id,
                exon.ref,
                exon.start,
                exon.end,
                exon.strand,
                itype,
                ovl,
                t_idx,
                g_idx,
            ))

    # Exact-match SJ lookup for each intron in the fragment
    for intron in frag.introns:
        key = (intron.ref, intron.start, intron.end, intron.strand)
        match = index.sj_map.get(key)
        if match is not None:
            t_set, g_set = match
            for t_idx in t_set:
                g_idx = int(index.t_to_g_arr[t_idx])
                rows.append((
                    frag_id,
                    intron.ref,
                    intron.start,
                    intron.end,
                    intron.strand,
                    IntervalType.SJ,
                    0,  # exact match, not overlap
                    t_idx,
                    g_idx,
                ))

    return rows


def _flush_buffer(buffer: list[tuple]) -> pa.RecordBatch:
    """Convert a list of hit tuples into a PyArrow RecordBatch.

    Columnar conversion is done in-place to avoid creating intermediate
    DataFrames.
    """
    n = len(buffer)
    frag_ids = np.empty(n, dtype=np.int32)
    refs = [None] * n
    starts = np.empty(n, dtype=np.int32)
    ends = np.empty(n, dtype=np.int32)
    strands = np.empty(n, dtype=np.int8)
    ivtypes = np.empty(n, dtype=np.int8)
    overlaps = np.empty(n, dtype=np.int32)
    t_indices = np.empty(n, dtype=np.int32)
    g_indices = np.empty(n, dtype=np.int32)

    for i, row in enumerate(buffer):
        frag_ids[i] = row[0]
        refs[i] = row[1]
        starts[i] = row[2]
        ends[i] = row[3]
        strands[i] = row[4]
        ivtypes[i] = row[5]
        overlaps[i] = row[6]
        t_indices[i] = row[7]
        g_indices[i] = row[8]

    return pa.record_batch(
        [frag_ids, refs, starts, ends, strands, ivtypes,
         overlaps, t_indices, g_indices],
        schema=HIT_SCHEMA,
    )


def query_bam(
    bam_iter,
    index: HulkIndex,
    output_path: str | Path,
    *,
    skip_duplicates: bool = True,
    multimap_bamfh=None,
    compression: str = "lz4",
    write_tsv: bool = False,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    log_every: int = 1_000_000,
) -> dict:
    """Scan a BAM file, query the index, and stream hits to a Feather file.

    Iterates over read pairs from *bam_iter*, builds ``Fragment`` objects,
    queries the ``HulkIndex`` for overlapping intervals and splice junction
    matches, and writes hit records in chunked RecordBatches to an Arrow
    IPC file (readable via ``pd.read_feather()``).

    Memory usage is bounded by *chunk_size* rows regardless of BAM size.

    Parameters
    ----------
    bam_iter
        An iterator over ``pysam.AlignedSegment`` objects (e.g., from
        ``pysam.AlignmentFile.fetch()``).
    index : HulkIndex
        The loaded reference index.
    output_path : str or Path
        Path for the output Feather (Arrow IPC) file.
    skip_duplicates : bool
        If True (default), discard reads marked as PCR/optical duplicates
        during BAM parsing. Retained duplicates are treated identically to
        other reads.
    multimap_bamfh : pysam.AlignmentFile or None
        Optional BAM handle to write multimapping reads (NH > 1).
    compression : str
        Arrow IPC compression (``'lz4'``, ``'zstd'``, or
        ``'uncompressed'``).
    write_tsv : bool
        If True, write a TSV mirror alongside the Feather file.
    chunk_size : int
        Number of hit rows to buffer before flushing to disk.
    log_every : int
        Log progress every *log_every* fragments at DEBUG level.

    Returns
    -------
    dict
        Stats dict with keys from ``parse_bam_file`` plus
        ``'n_fragments'`` and ``'n_hits'``.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Map compression string to pyarrow codec
    if compression == "uncompressed":
        ipc_compression = None
    else:
        ipc_compression = compression

    options = pa.ipc.IpcWriteOptions(compression=ipc_compression)

    stats: dict = {}
    buffer: list[tuple] = []
    frag_id = 0
    total_hits = 0

    logger.info(f"[START] Scanning BAM and streaming hits to {output_path}")

    with pa.OSFile(str(output_path), "wb") as sink:
        writer = pa.ipc.new_file(sink, HIT_SCHEMA, options=options)

        for r1, r2 in parse_bam_file(
            bam_iter, stats,
            skip_duplicates=skip_duplicates,
            multimap_bamfh=multimap_bamfh,
        ):
            frag = Fragment.from_reads(r1, r2)
            hits = _query_fragment(frag, frag_id, index)
            buffer.extend(hits)
            frag_id += 1

            # Flush buffer when it reaches chunk_size
            if len(buffer) >= chunk_size:
                writer.write_batch(_flush_buffer(buffer))
                total_hits += len(buffer)
                buffer.clear()

            if frag_id % log_every == 0:
                logger.debug(
                    f"  processed {frag_id:,} fragments, "
                    f"{total_hits + len(buffer):,} hits"
                )

        # Flush remaining rows
        if buffer:
            writer.write_batch(_flush_buffer(buffer))
            total_hits += len(buffer)
            buffer.clear()

        writer.close()

    logger.info(f"[DONE] Scanned {frag_id:,} fragments, "
                f"{total_hits:,} total hits → {output_path}")

    stats["n_fragments"] = frag_id
    stats["n_hits"] = total_hits

    # Optional TSV mirror (reads back from the Feather file)
    if write_tsv:
        tsv_path = output_path.with_suffix(".tsv")
        df = pd.read_feather(output_path)
        df.to_csv(tsv_path, sep="\t", index=False)
        logger.debug(f"  TSV mirror: {tsv_path}")

    return stats


def read_hits(path: str | Path) -> pd.DataFrame:
    """Read a hit Feather file back into a DataFrame.

    Convenience wrapper around ``pd.read_feather()`` for the
    intermediate hit file written by ``query_bam()``.

    Parameters
    ----------
    path : str or Path
        Path to the Feather file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns matching ``HIT_COLUMNS``.
    """
    return pd.read_feather(str(path))
