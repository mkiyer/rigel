"""
rigel.index — Build, load, and query the rigel reference index.

The index is constructed from a genome FASTA (with samtools .fai) and a
GENCODE GTF annotation file. It produces four Feather files (with optional
TSV mirrors) in an output directory:

    ref_lengths.feather   — reference names and lengths
    transcripts.feather   — one row per transcript with integer indices
    intervals.feather     — exon/intron/intergenic tiling of the genome
    sj.feather            — annotated splice junctions from transcript introns

The ``TranscriptIndex`` class provides both the ``build()`` method for creating
the index and ``load()`` / query methods for using it during quantification.
"""

import collections
import functools
import logging
import os
from pathlib import Path
from typing import Iterator, Literal

from .native import cgranges as _cgranges_cls
import numpy as np
import pandas as pd
import pysam

from .types import GenomicInterval, IntervalType, AnnotatedInterval, STRAND_POS
from .transcript import Transcript


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Output file names
# ---------------------------------------------------------------------------

REF_LENGTHS_FEATHER = "ref_lengths.feather"
REF_LENGTHS_TSV = "ref_lengths.tsv"

TRANSCRIPTS_FEATHER = "transcripts.feather"
TRANSCRIPTS_TSV = "transcripts.tsv"

INTERVALS_FEATHER = "intervals.feather"
INTERVALS_TSV = "intervals.tsv"

SJ_FEATHER = "sj.feather"
SJ_TSV = "sj.tsv"

REGIONS_FEATHER = "regions.feather"
REGIONS_TSV = "regions.tsv"


# ---------------------------------------------------------------------------
# Build helpers (public for testability)
# ---------------------------------------------------------------------------

def load_reference_lengths(fasta_file: str | Path) -> dict[str, int]:
    """Read reference names and lengths from a FASTA file.

    The FASTA must be indexed (``samtools faidx``). We read the index
    rather than scanning the full file.

    Returns
    -------
    dict[str, int]
        Ordered mapping from reference name to length.
    """
    fasta_file = str(fasta_file)
    with pysam.FastaFile(fasta_file) as fh:
        ref_lengths = {
            fh.references[i]: fh.lengths[i]
            for i in range(fh.nreferences)
        }
    return ref_lengths


def read_transcripts(
    gtf_file: str | Path,
    *,
    gtf_parse_mode: Literal["strict", "warn-skip"] = "strict",
) -> list[Transcript]:
    """Parse a GTF file and return a sorted list of Transcript objects.

    Each transcript is assigned a sequential ``t_index`` and a ``g_index``
    (unique integer per ``g_id``). Transcripts are sorted by
    ``(ref, start, end, strand)`` so that downstream interval generation
    can process them in genomic order.
    """
    transcripts = Transcript.read_gtf(
        str(gtf_file), parse_mode=gtf_parse_mode,
    )

    # Assign integer indices
    g_id_to_index: dict[str, int] = {}
    for ti, t in enumerate(transcripts):
        t.t_index = ti
        if t.g_id not in g_id_to_index:
            g_id_to_index[t.g_id] = len(g_id_to_index)
        t.g_index = g_id_to_index[t.g_id]

    return transcripts


def transcripts_to_dataframe(transcripts: list[Transcript]) -> pd.DataFrame:
    """Convert a list of Transcript objects to a pandas DataFrame."""
    return pd.DataFrame(t.to_dict() for t in transcripts)



# -- Tolerance-based nRNA merging (unified architecture) ----------------------

#: Default merge tolerance (bp) for TSS/TES clustering.
NRNA_MERGE_TOLERANCE: int = 20


def _cluster_coordinates(coords: np.ndarray, tolerance: int) -> np.ndarray:
    """Assign sorted coordinates to clusters within *tolerance* bp.

    Returns an int array of the same length, where each element is the
    cluster id (0-based).  Coordinates must be sorted ascending.
    """
    n = len(coords)
    if n == 0:
        return np.empty(0, dtype=np.intp)
    ids = np.empty(n, dtype=np.intp)
    cid = 0
    anchor = 0  # index of the first element in the current cluster
    for i in range(n):
        if coords[i] - coords[anchor] > tolerance:
            cid += 1
            anchor = i
        ids[i] = cid
    return ids


def create_nrna_transcripts(
    transcripts: list[Transcript],
    tolerance: int = NRNA_MERGE_TOLERANCE,
) -> tuple[list[Transcript], dict[int, tuple]]:
    """Create synthetic nRNA transcripts, detect annotated equivalents.

    This function implements the unified nRNA architecture:

    1. Collect (ref, strand, start, end) for each **multi-exon** transcript.
    2. Cluster nearby TSS and TES within *tolerance* bp independently.
       The merged span uses the **outer envelope** (min start, max end)
       so synthetics are never smaller than any contributing transcript.
    3. Check if an existing **single-exon** transcript already covers
       each merged span (exact match or full containment).
       If so, mark it ``is_nascent_equiv = True`` — no synthetic needed.
    4. For uncovered spans, create a synthetic single-exon ``Transcript``
       flagged ``is_synthetic_nrna = True``.
    5. The synthetics inherit the gene id / name of one contributing
       transcript for output convenience.

    Parameters
    ----------
    transcripts : list[Transcript]
        Sorted transcript list (``t_index`` already assigned).
    tolerance : int
        Maximum distance (bp) for TSS/TES clustering.

    Returns
    -------
    synthetics : list[Transcript]
        Synthetic nRNA transcripts to append to the transcript list.
        Callers must assign ``t_index`` to each returned transcript.
    t_to_span_key : dict[int, tuple]
        Mapping from annotated multi-exon transcript ``t_index`` to its
        merged nRNA span key ``(ref, strand, start, end)``.  Used by
        the caller to set ``nrna_t_index`` after ``t_index`` assignment.
    span_to_syn_idx : dict[tuple, int]
        Mapping from uncovered span key to the index in *synthetics*.
    covered_equiv : dict[tuple, Transcript]
        Mapping from covered span key to the annotated nascent-equiv
        transcript.
    """
    from .types import Interval, Strand

    # -- Step 1: Collect multi-exon transcript spans --------------------------
    multi_exon: list[Transcript] = []
    for t in transcripts:
        if len(t.exons) >= 2:
            multi_exon.append(t)

    if not multi_exon:
        return [], {}, {}, {}

    # -- Step 2: TSS/TES clustering per (ref, strand) ------------------------
    # Group multi-exon transcripts by (ref, strand)
    groups: dict[tuple[str, int], list[Transcript]] = collections.defaultdict(list)
    for t in multi_exon:
        groups[(t.ref, int(t.strand))].append(t)

    # For each transcript, compute merged (representative) start and end.
    # Store as {t_index: (merged_start, merged_end)}.
    t_to_merged_span: dict[int, tuple[int, int]] = {}

    for (ref, strand), grp in groups.items():
        # Cluster starts
        starts = np.array([t.start for t in grp], dtype=np.int64)
        order_s = np.argsort(starts, kind="mergesort")
        sorted_starts = starts[order_s]
        cids_s = _cluster_coordinates(sorted_starts, tolerance)

        # Map cluster_id → min start (outer envelope)
        n_clusters_s = cids_s[-1] + 1 if len(cids_s) else 0
        cid_min_start = np.full(n_clusters_s, np.iinfo(np.int64).max, dtype=np.int64)
        np.minimum.at(cid_min_start, cids_s, sorted_starts)

        # Back-map to original order
        rep_starts = np.empty(len(grp), dtype=np.int64)
        for i, oi in enumerate(order_s):
            rep_starts[oi] = cid_min_start[cids_s[i]]

        # Cluster ends
        ends = np.array([t.end for t in grp], dtype=np.int64)
        order_e = np.argsort(ends, kind="mergesort")
        sorted_ends = ends[order_e]
        cids_e = _cluster_coordinates(sorted_ends, tolerance)

        n_clusters_e = cids_e[-1] + 1 if len(cids_e) else 0
        cid_max_end = np.full(n_clusters_e, np.iinfo(np.int64).min, dtype=np.int64)
        np.maximum.at(cid_max_end, cids_e, sorted_ends)

        rep_ends = np.empty(len(grp), dtype=np.int64)
        for i, oi in enumerate(order_e):
            rep_ends[oi] = cid_max_end[cids_e[i]]

        for j, t in enumerate(grp):
            t_to_merged_span[t.t_index] = (int(rep_starts[j]), int(rep_ends[j]))

    # -- Step 3: Dedup merged spans -------------------------------------------
    # Map merged_key → (first transcript for metadata, set of contributing t_indices)
    merged_spans: dict[tuple, Transcript] = {}  # key → representative transcript
    span_contributor_count: dict[tuple, int] = collections.defaultdict(int)
    t_to_span_key: dict[int, tuple] = {}  # annotated t_index → span key
    for t in multi_exon:
        ms, me = t_to_merged_span[t.t_index]
        key = (t.ref, int(t.strand), ms, me)
        if key not in merged_spans:
            merged_spans[key] = t  # keep first for gene metadata
        span_contributor_count[key] += 1
        t_to_span_key[t.t_index] = key

    # -- Step 4: Detect annotated equivalents ---------------------------------
    # Build single-exon lookup: (ref, strand) → sorted list of (start, end, t)
    se_by_loc: dict[tuple[str, int], list[tuple[int, int, Transcript]]] = (
        collections.defaultdict(list)
    )
    for t in transcripts:
        if len(t.exons) == 1:
            se_by_loc[(t.ref, int(t.strand))].append((t.start, t.end, t))
    for key in se_by_loc:
        se_by_loc[key].sort(key=lambda x: (x[0], x[1]))

    covered: set[tuple] = set()
    covered_equiv: dict[tuple, Transcript] = {}  # span_key → equiv transcript
    for span_key in merged_spans:
        ref, strand, m_start, m_end = span_key
        candidates = se_by_loc.get((ref, strand), [])
        for s_start, s_end, s_tx in candidates:
            if s_start > m_start:
                break  # sorted by start, no more can contain
            if s_start <= m_start and m_end <= s_end:
                s_tx.is_nascent_equiv = True
                s_tx.nrna_n_contributors = span_contributor_count[span_key]
                covered.add(span_key)
                covered_equiv[span_key] = s_tx
                break

    # -- Step 5: Create synthetic transcripts ---------------------------------
    synthetics: list[Transcript] = []
    span_to_syn_idx: dict[tuple, int] = {}  # span_key → index in synthetics list
    for span_key, rep_tx in merged_spans.items():
        if span_key in covered:
            continue
        ref, strand_int, s, e = span_key
        syn = Transcript(
            ref=ref,
            strand=Strand(strand_int),
            exons=[Interval(s, e)],
            length=e - s,
            t_id=f"RIGEL_NRNA_{ref}_{strand_int}_{s}_{e}",
            g_id=rep_tx.g_id,
            g_name=rep_tx.g_name,
            g_type=rep_tx.g_type,
            is_synthetic_nrna=True,
            nrna_n_contributors=span_contributor_count[span_key],
        )
        span_to_syn_idx[span_key] = len(synthetics)
        synthetics.append(syn)

    logger.info(
        f"nRNA: {len(merged_spans)} merged spans "
        f"({len(covered)} annotated equiv, "
        f"{len(synthetics)} synthetic)"
    )
    return synthetics, t_to_span_key, span_to_syn_idx, covered_equiv


def build_splice_junctions(transcripts: list[Transcript]) -> pd.DataFrame:
    """Extract splice junctions (introns) from all transcripts.

    Each intron boundary within a transcript produces one SpliceJunction
    record. The resulting DataFrame is sorted by (ref, start, end, strand).
    """
    rows = [
        AnnotatedInterval(t.ref, start, end, t.strand,
                    IntervalType.SJ, t.t_index)
        for t in transcripts
        for start, end in t.introns()
    ]
    rows.sort(key=lambda sj: (sj.ref, sj.start, sj.end, sj.strand))
    return pd.DataFrame(rows, columns=AnnotatedInterval._fields)


def write_bed12(
    transcripts: list[Transcript],
    bed_path: str | Path,
) -> Path:
    """Write transcripts as BED12 for minimap2 ``-j`` annotation.

    Each transcript becomes one BED12 line with exon blocks.
    Coordinates are already 0-based half-open (BED convention).

    Parameters
    ----------
    transcripts : list of Transcript
        Transcripts to write (must have exons populated).
    bed_path : str or Path
        Output BED12 file path.

    Returns
    -------
    Path
        The written BED12 file path.
    """
    bed_path = Path(bed_path)
    with open(bed_path, "w") as fh:
        for t in transcripts:
            if not t.exons:
                continue
            ref = t.ref
            ref_start = t.exons[0].start
            ref_end = t.exons[-1].end
            name = t.t_id or "."
            score = 0
            strand = t.strand.to_str()
            if strand not in ("+", "-"):
                strand = "+"
            thick_start = ref_start
            thick_end = ref_end
            rgb = "0"
            block_count = len(t.exons)
            block_sizes = ",".join(
                str(e.end - e.start) for e in t.exons
            )
            block_starts = ",".join(
                str(e.start - ref_start) for e in t.exons
            )
            fh.write(
                f"{ref}\t{ref_start}\t{ref_end}\t{name}\t"
                f"{score}\t{strand}\t{thick_start}\t{thick_end}\t"
                f"{rgb}\t{block_count}\t{block_sizes}\t{block_starts}\n"
            )
    return bed_path


def gtf_to_bed12(gtf_path: str | Path, bed_path: str | Path) -> Path:
    """Convert a GTF file to BED12 for minimap2 ``-j`` annotation.

    Convenience wrapper around ``read_transcripts()`` + ``write_bed12()``.

    Parameters
    ----------
    gtf_path : str or Path
        Input GTF file (may be gzipped).
    bed_path : str or Path
        Output BED12 file path.

    Returns
    -------
    Path
        The written BED12 file path.
    """
    transcripts = read_transcripts(gtf_path)
    return write_bed12(transcripts, bed_path)


def _gen_transcript_intervals(t: Transcript) -> Iterator[AnnotatedInterval]:
    """Yield exon and transcript-span intervals for a single transcript.

    Each exon produces an EXON interval.  One TRANSCRIPT interval spans
    the full transcript ``[first_exon.start, last_exon.end)``.
    """
    # Exons
    for e in t.exons:
        yield AnnotatedInterval(t.ref, e.start, e.end, t.strand,
                          IntervalType.EXON, t.t_index)
    # One transcript span (replaces per-gap INTRON intervals)
    if t.exons:
        yield AnnotatedInterval(t.ref, t.exons[0].start, t.exons[-1].end,
                          t.strand, IntervalType.TRANSCRIPT, t.t_index)


def _gen_genomic_intervals(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> Iterator[AnnotatedInterval]:
    """Tile the genome into exon, intron, and intergenic intervals.

    Transcripts must be sorted by ``(ref, start, end)``. Intergenic
    intervals fill the gaps between transcript clusters and extend to
    reference boundaries.
    """
    # Group transcripts by reference
    ref_transcripts: dict[str, list[Transcript]] = collections.defaultdict(list)
    for t in transcripts:
        if t.ref not in ref_lengths:
            raise ValueError(
                f"Transcript {t.t_id} has reference '{t.ref}' "
                f"not found in the FASTA index"
            )
        ref_transcripts[t.ref].append(t)

    # Process each reference
    for ref, ref_length in ref_lengths.items():
        t_list = ref_transcripts.get(ref, [])

        if not t_list:
            # Entire reference is intergenic
            yield AnnotatedInterval(ref, 0, ref_length)
            continue

        end = 0  # tracks the rightmost coordinate reached
        cluster: list[Transcript] = []

        for t in t_list:
            if t.start > end:
                # Intergenic gap before this transcript cluster
                if end < t.start:
                    yield AnnotatedInterval(ref, end, t.start)
                # Emit intervals for the previous cluster
                for tc in cluster:
                    yield from _gen_transcript_intervals(tc)
                cluster = []

            end = max(end, t.end)
            cluster.append(t)

        # Emit final cluster
        for tc in cluster:
            yield from _gen_transcript_intervals(tc)

        # Intergenic gap to end of reference
        if end < ref_length:
            yield AnnotatedInterval(ref, end, ref_length)


def build_genomic_intervals(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> pd.DataFrame:
    """Build the complete interval table covering the genome.

    Returns a DataFrame sorted by (ref, start, end, strand) with columns
    matching the ``AnnotatedInterval`` NamedTuple fields.
    """
    intervals = list(_gen_genomic_intervals(transcripts, ref_lengths))
    intervals.sort(key=lambda iv: (iv.ref, iv.start, iv.end, iv.strand))
    return pd.DataFrame(intervals, columns=AnnotatedInterval._fields)


def build_region_table(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> pd.DataFrame:
    """Build a non-overlapping, reference-complete genomic region partition.

    The partition is constructed by a boundary-sweep algorithm:

    1. Collect all transcript-span and exon start/end coordinates plus
       reference boundaries (0 and ref_length) per reference.
    2. Sort and deduplicate into a boundary array.
    3. Form atomic half-open bins between successive boundaries.
    4. Assign four boolean flags per bin via ``searchsorted``:
       ``exon_pos``, ``exon_neg``, ``tx_pos``, ``tx_neg``.
    5. Assign sequential ``region_id`` across all references.

    Parameters
    ----------
    transcripts : list[Transcript]
        Sorted list of transcripts (must have ``t_index`` assigned).
    ref_lengths : dict[str, int]
        Ordered mapping from reference name to length.

    Returns
    -------
    pd.DataFrame
        Columns: region_id, ref, start, end, length,
        exon_pos, exon_neg, tx_pos, tx_neg.
    """
    # Group transcripts by reference
    ref_transcripts: dict[str, list[Transcript]] = collections.defaultdict(list)
    for t in transcripts:
        ref_transcripts[t.ref].append(t)

    rows: list[dict] = []
    region_id = 0

    for ref, ref_length in ref_lengths.items():
        t_list = ref_transcripts.get(ref, [])

        # -- Collect boundary coordinates ------------------------------------
        boundaries = {0, ref_length}
        for t in t_list:
            boundaries.add(t.start)
            boundaries.add(t.end)
            for e in t.exons:
                boundaries.add(e.start)
                boundaries.add(e.end)

        bounds = np.array(sorted(boundaries), dtype=np.int64)
        n_bins = len(bounds) - 1
        if n_bins <= 0:
            continue

        bin_starts = bounds[:-1]
        bin_ends = bounds[1:]

        # -- Flag arrays (one per bin) ---------------------------------------
        exon_pos = np.zeros(n_bins, dtype=bool)
        exon_neg = np.zeros(n_bins, dtype=bool)
        tx_pos = np.zeros(n_bins, dtype=bool)
        tx_neg = np.zeros(n_bins, dtype=bool)

        for t in t_list:
            is_pos = int(t.strand) == STRAND_POS
            # Transcript span -> tx flag
            lo = int(np.searchsorted(bounds, t.start, side="left"))
            hi = int(np.searchsorted(bounds, t.end, side="left"))
            if is_pos:
                tx_pos[lo:hi] = True
            else:
                tx_neg[lo:hi] = True

            # Exons -> exon flag (also implies tx flag)
            for e in t.exons:
                e_lo = int(np.searchsorted(bounds, e.start, side="left"))
                e_hi = int(np.searchsorted(bounds, e.end, side="left"))
                if is_pos:
                    exon_pos[e_lo:e_hi] = True
                    tx_pos[e_lo:e_hi] = True
                else:
                    exon_neg[e_lo:e_hi] = True
                    tx_neg[e_lo:e_hi] = True

        # -- Merge adjacent bins with identical flags --------------------------
        # Encode the four booleans as a single byte for fast comparison.
        flag_key = (exon_pos.astype(np.uint8)
                    | (exon_neg.astype(np.uint8) << 1)
                    | (tx_pos.astype(np.uint8) << 2)
                    | (tx_neg.astype(np.uint8) << 3))

        # Detect where flag signature changes between adjacent bins.
        changed = np.empty(n_bins, dtype=bool)
        changed[0] = True
        if n_bins > 1:
            changed[1:] = flag_key[1:] != flag_key[:-1]
        group_starts_idx = np.where(changed)[0]

        for gi in range(len(group_starts_idx)):
            first = int(group_starts_idx[gi])
            last = int(group_starts_idx[gi + 1] - 1) if gi + 1 < len(group_starts_idx) else n_bins - 1
            s = int(bin_starts[first])
            e = int(bin_ends[last])
            if s >= e:
                continue
            rows.append({
                "region_id": region_id,
                "ref": ref,
                "start": s,
                "end": e,
                "length": e - s,
                "exon_pos": bool(exon_pos[first]),
                "exon_neg": bool(exon_neg[first]),
                "tx_pos": bool(tx_pos[first]),
                "tx_neg": bool(tx_neg[first]),
            })
            region_id += 1

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# TranscriptIndex — unified index class
# ---------------------------------------------------------------------------

class TranscriptIndex:
    """In-memory reference index for fast transcript/gene overlap queries.

    Use ``TranscriptIndex.build()`` to create the on-disk index from a FASTA +
    GTF, then ``TranscriptIndex.load()`` to read it back for quantification.

    The gene table is derived on-the-fly from the transcript table so
    that only transcripts, intervals, and splice junctions are stored
    on disk.
    """

    def __init__(self):
        self.index_dir: str | None = None
        self.t_df: pd.DataFrame | None = None
        self.g_df: pd.DataFrame | None = None
        self.t_to_g_arr: np.ndarray | None = None
        self.t_to_strand_arr: np.ndarray | None = None
        self.g_to_strand_arr: np.ndarray | None = None

        # Unified cgranges index (collapsed EXON + TRANSCRIPT + INTERGENIC)
        # Label = index into _iv_type and _iv_t_set lookup lists.
        self.cr: _cgranges_cls | None = None
        self._iv_type: list[int] | None = None
        self._iv_t_set: list[frozenset[int]] | None = None

        # Splice-junction exact-match map  (ref, start, end, strand) → frozenset[int]
        self.sj_map: dict | None = None

        # Per-transcript exon intervals for coverage-weight model.
        # Maps t_index → (n_exons, 2) int32 array of [start, end) intervals
        # sorted by genomic start position.
        self._t_exon_intervals: dict[int, np.ndarray] | None = None

    # -- properties -----------------------------------------------------------

    @property
    def num_transcripts(self) -> int:
        return len(self.t_df)

    @property
    def num_genes(self) -> int:
        return len(self.g_df)

    # -- lazy-loaded calibration regions --------------------------------------

    @functools.cached_property
    def region_df(self) -> pd.DataFrame | None:
        """Calibration-region partition table, loaded on first access."""
        if self.index_dir is None:
            return None
        regions_path = os.path.join(self.index_dir, REGIONS_FEATHER)
        if not os.path.exists(regions_path):
            return None
        df = pd.read_feather(regions_path)
        logger.debug(f"Loaded {len(df)} calibration regions (lazy)")
        return df

    @functools.cached_property
    def region_cr(self):
        """cgranges index over calibration regions, built on first access."""
        df = self.region_df
        if df is None or len(df) == 0:
            return None
        _r_refs = df["ref"].values.tolist()
        _r_starts = df["start"].values.tolist()
        _r_ends = df["end"].values.tolist()
        _r_ids = df["region_id"].values.tolist()
        region_cr = _cgranges_cls()
        for i in range(len(_r_ids)):
            region_cr.add(_r_refs[i], _r_starts[i], _r_ends[i], _r_ids[i])
        region_cr.index()
        return region_cr

    # -- build (static) -------------------------------------------------------

    @staticmethod
    def build(
        fasta_file: str | Path,
        gtf_file: str | Path,
        output_dir: str | Path,
        *,
        feather_compression: str = "lz4",
        write_tsv: bool = True,
        gtf_parse_mode: Literal["strict", "warn-skip"] = "strict",
        nrna_tolerance: int = NRNA_MERGE_TOLERANCE,
    ) -> None:
        """Build the rigel reference index and write to disk.

        Parameters
        ----------
        fasta_file : path
            Genome FASTA (must be indexed with ``samtools faidx``).
        gtf_file : path
            Gene annotation in GTF format (GENCODE recommended).
        output_dir : path
            Directory to write index files into (created if needed).
        feather_compression : str
            Compression for Feather files (``'lz4'``, ``'zstd'``, or
            ``'uncompressed'``).
        write_tsv : bool
            If True, write human-readable TSV mirrors alongside Feather files.
        gtf_parse_mode : {"strict", "warn-skip"}
            GTF parsing behavior. ``"strict"`` (default) fails fast on
            malformed lines; ``"warn-skip"`` logs warnings and skips.
        nrna_tolerance : int
            Max distance (bp) for clustering transcript start/end sites
            when building synthetic nascent RNA transcripts.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        feather_kwargs = {"compression": feather_compression}

        # -- Reference lengths ------------------------------------------------
        logger.info(f"[START] Loading references from {fasta_file}")
        ref_lengths = load_reference_lengths(fasta_file)
        logger.info(f"[DONE] Read {len(ref_lengths)} references")

        df = pd.DataFrame(
            list(ref_lengths.items()), columns=["ref", "length"]
        )
        df.to_feather(output_dir / REF_LENGTHS_FEATHER, **feather_kwargs)
        if write_tsv:
            df.to_csv(output_dir / REF_LENGTHS_TSV, sep="\t", index=False)

        # -- Transcripts ------------------------------------------------------
        logger.info(f"[START] Reading transcripts from {gtf_file}")
        transcripts = read_transcripts(
            gtf_file, gtf_parse_mode=gtf_parse_mode,
        )
        logger.info(f"[DONE] Read {len(transcripts)} transcripts")

        # -- Synthetic nRNA transcripts ---------------------------------------
        logger.info("[START] Creating synthetic nRNA transcripts")
        synthetics, t_to_span_key, span_to_syn_idx, covered_equiv = (
            create_nrna_transcripts(transcripts, tolerance=nrna_tolerance)
        )
        # Build g_id → g_index mapping from annotated transcripts
        g_id_to_gindex: dict[str, int] = {}
        for t in transcripts:
            if t.g_id not in g_id_to_gindex:
                g_id_to_gindex[t.g_id] = t.g_index

        # Assign t_index and g_index to each synthetic; build span→t_index
        span_to_nrna_t_index: dict[tuple, int] = {}
        next_t_index = len(transcripts)
        for syn in synthetics:
            syn.t_index = next_t_index
            next_t_index += 1
            syn.g_index = g_id_to_gindex.get(syn.g_id, -1)
            span_key = (syn.ref, int(syn.strand), syn.start, syn.end)
            span_to_nrna_t_index[span_key] = syn.t_index

        # For covered spans, the nRNA entity is the nascent-equiv transcript
        for span_key, equiv_tx in covered_equiv.items():
            span_to_nrna_t_index[span_key] = equiv_tx.t_index
            equiv_tx.nrna_t_index = equiv_tx.t_index

        # Set nrna_t_index on every multi-exon annotated transcript
        for me_tidx, span_key in t_to_span_key.items():
            transcripts[me_tidx].nrna_t_index = span_to_nrna_t_index.get(
                span_key, -1
            )

        transcripts.extend(synthetics)
        logger.info(
            f"[DONE] {len(synthetics)} synthetics added, "
            f"{len(transcripts)} total transcripts"
        )

        # Partition: annotated-only list for region table and splice junctions
        annotated_transcripts = [
            t for t in transcripts if not t.is_synthetic_nrna
        ]

        t_df = transcripts_to_dataframe(transcripts)
        t_df.to_feather(output_dir / TRANSCRIPTS_FEATHER, **feather_kwargs)
        if write_tsv:
            t_df.to_csv(output_dir / TRANSCRIPTS_TSV, sep="\t", index=False)

        # -- Splice junctions -------------------------------------------------
        logger.info("[START] Building splice junctions")
        sj_df = build_splice_junctions(annotated_transcripts)
        logger.info(f"[DONE] Found {len(sj_df)} splice junctions")

        sj_df.to_feather(output_dir / SJ_FEATHER, **feather_kwargs)
        if write_tsv:
            sj_df.to_csv(output_dir / SJ_TSV, sep="\t", index=False)

        # -- Genomic intervals ------------------------------------------------
        logger.info("[START] Building genomic intervals")
        iv_df = build_genomic_intervals(transcripts, ref_lengths)
        logger.info(f"[DONE] Found {len(iv_df)} genomic intervals")

        iv_df.to_feather(output_dir / INTERVALS_FEATHER, **feather_kwargs)
        if write_tsv:
            iv_df.to_csv(output_dir / INTERVALS_TSV, sep="\t", index=False)

        # -- Calibration regions ----------------------------------------------
        logger.info("[START] Building calibration regions")
        region_df = build_region_table(annotated_transcripts, ref_lengths)
        logger.info(f"[DONE] Found {len(region_df)} calibration regions")

        region_df.to_feather(output_dir / REGIONS_FEATHER, **feather_kwargs)
        if write_tsv:
            region_df.to_csv(output_dir / REGIONS_TSV, sep="\t", index=False)

        logger.info(f"Index written to {output_dir}")

    # -- load (classmethod) ---------------------------------------------------

    @staticmethod
    def _build_gene_table(t_df: pd.DataFrame) -> pd.DataFrame:
        """Derive a gene-level summary table from the transcript table."""
        g_df = (
            t_df
            .groupby("g_index", sort=True)
            .agg(
                ref=("ref", "first"),
                start=("start", "min"),
                end=("end", "max"),
                strand=("strand", "first"),
                g_id=("g_id", "first"),
                g_name=("g_name", "first"),
                g_type=("g_type", "first"),
                num_transcripts=("t_index", "count"),
            )
            .reset_index()
        )
        return g_df

    @classmethod
    def load(cls, index_dir: str | Path) -> "TranscriptIndex":
        """Load an index from the Feather files in *index_dir*.

        Builds cgranges interval trees and splice-junction lookup maps
        for fast overlap queries during quantification.
        """
        index_dir = str(index_dir)
        self = cls()
        self.index_dir = index_dir

        # -- transcripts ------------------------------------------------------
        self.t_df = pd.read_feather(
            os.path.join(index_dir, TRANSCRIPTS_FEATHER)
        )
        if "t_index" not in self.t_df.columns:
            raise ValueError(
                f"Invalid index in {index_dir}: missing 't_index' column "
                f"in {TRANSCRIPTS_FEATHER}"
            )
        if not (self.t_df.index == self.t_df["t_index"]).all():
            raise ValueError(
                f"Invalid index in {index_dir}: row index does not match "
                f"'t_index' column in {TRANSCRIPTS_FEATHER}; rebuild index"
            )

        # Backward compat: older indices lack synthetic nRNA columns
        if "is_synthetic_nrna" not in self.t_df.columns:
            self.t_df["is_synthetic_nrna"] = False
        if "is_nascent_equiv" not in self.t_df.columns:
            self.t_df["is_nascent_equiv"] = False

        # -- gene table (derived) ---------------------------------------------
        self.g_df = cls._build_gene_table(self.t_df)
        if "g_index" not in self.g_df.columns:
            raise ValueError(
                f"Invalid derived gene table in {index_dir}: missing "
                f"'g_index' column"
            )
        if not (self.g_df.index == self.g_df["g_index"]).all():
            raise ValueError(
                f"Invalid index in {index_dir}: derived gene table row index "
                f"does not match 'g_index'; rebuild index"
            )

        # -- fast numpy lookup arrays -----------------------------------------
        self.t_to_g_arr = self.t_df["g_index"].values
        self.t_to_strand_arr = self.t_df["strand"].values
        self.g_to_strand_arr = self.g_df["strand"].values

        # -- interval index (unified cgranges) ---------------------------------
        logger.debug("Reading intervals")
        iv_df = pd.read_feather(
            os.path.join(index_dir, INTERVALS_FEATHER)
        )

        # -- collapsed interval index -----------------------------------------
        # Group rows by (ref, start, end, interval_type) and merge
        # transcript indices into frozensets.  Each unique boundary is
        # stored once in cgranges, with a label indexing into _iv_type
        # and _iv_t_set lookup lists.
        #
        # Vectorised: encode refs as integers, sort, detect boundaries via
        # diff, then build frozensets in a single pass.
        logger.debug("Building collapsed interval index")
        _iv_refs = iv_df["ref"].values
        _iv_starts = iv_df["start"].values
        _iv_ends = iv_df["end"].values
        _iv_itypes = iv_df["interval_type"].values
        _iv_tidxs = iv_df["t_index"].values

        # Map string refs → integer codes for fast comparison
        _unique_refs, _ref_codes = np.unique(_iv_refs, return_inverse=True)

        # Sort by (ref_code, start, end, itype) to group identical keys
        _sort_order = np.lexsort((_iv_itypes, _iv_ends, _iv_starts, _ref_codes))
        _rc = _ref_codes[_sort_order]
        _ss = _iv_starts[_sort_order]
        _ee = _iv_ends[_sort_order]
        _it = _iv_itypes[_sort_order]
        _ti = _iv_tidxs[_sort_order]
        _rr = _iv_refs[_sort_order]

        # Detect group boundaries where any key column changes
        n = len(_sort_order)
        _changed = np.empty(n, dtype=bool)
        if n > 0:
            _changed[0] = True
        if n > 1:
            _changed[1:] = (
                (np.diff(_rc) != 0)
                | (np.diff(_ss) != 0)
                | (np.diff(_ee) != 0)
                | (np.diff(_it) != 0)
            )
        _group_starts = np.where(_changed)[0]
        _group_ends = np.append(_group_starts[1:], n)

        # Pre-convert numpy arrays to Python lists once (avoids per-element
        # int()/str() overhead in the 1.2 M-iteration loop below).
        _rr_list = _rr.tolist()
        _ss_list = _ss.tolist()
        _ee_list = _ee.tolist()
        _it_list = _it.tolist()
        _ti_list = _ti.tolist()
        _gs_list = _group_starts.tolist()
        _ge_list = _group_ends.tolist()

        cr = _cgranges_cls()
        iv_type: list[int] = []
        iv_t_set: list[frozenset[int]] = []
        _collapse_keys: list[tuple] = []  # keep for FragmentResolver later
        for label in range(len(_gs_list)):
            s = _gs_list[label]
            e = _ge_list[label]
            ref = _rr_list[s]
            start = _ss_list[s]
            end = _ee_list[s]
            itype = _it_list[s]
            tset = frozenset(t for t in _ti_list[s:e] if t >= 0)
            cr.add(ref, start, end, label)
            iv_type.append(itype)
            iv_t_set.append(tset)
            _collapse_keys.append((ref, start, end, itype))
        cr.index()
        self.cr = cr
        self._iv_type = iv_type
        self._iv_t_set = iv_t_set
        logger.debug(f"Interval index: {len(iv_df)} rows → {len(iv_type)} collapsed")

        # -- per-transcript exon intervals for coverage-weight model ----------
        # Vectorised: extract raw arrays, sort by (t_index, start), then
        # split at group boundaries — avoids pandas groupby overhead.
        _iv_itype_arr = iv_df["interval_type"].values
        _iv_tidx_arr = iv_df["t_index"].values
        _iv_start_arr = iv_df["start"].values
        _iv_end_arr = iv_df["end"].values

        exon_mask = _iv_itype_arr == int(IntervalType.EXON)
        _ex_tidx = _iv_tidx_arr[exon_mask]
        _ex_start = _iv_start_arr[exon_mask]
        _ex_end = _iv_end_arr[exon_mask]

        # Drop negative t_index (intergenic sentinel)
        valid = _ex_tidx >= 0
        _ex_tidx = _ex_tidx[valid]
        _ex_start = _ex_start[valid]
        _ex_end = _ex_end[valid]

        # Sort by (t_index, start) so exons within each transcript are ordered
        order = np.lexsort((_ex_start, _ex_tidx))
        _ex_tidx = _ex_tidx[order]
        _ex_start = _ex_start[order]
        _ex_end = _ex_end[order]

        # Find group boundaries via np.unique
        unique_tidx, group_starts, group_counts = np.unique(
            _ex_tidx, return_index=True, return_counts=True
        )
        # Pre-stack start/end into a single int32 array for fast slicing
        _ex_coords = np.column_stack((_ex_start, _ex_end)).astype(np.int32)

        t_exon_intervals: dict[int, np.ndarray] = {}
        for i in range(len(unique_tidx)):
            s = group_starts[i]
            t_exon_intervals[int(unique_tidx[i])] = _ex_coords[s : s + group_counts[i]]
        self._t_exon_intervals = t_exon_intervals
        logger.debug(f"Cached exon intervals for {len(t_exon_intervals)} transcripts")

        # -- splice junction indexes ------------------------------------------
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(
            os.path.join(index_dir, SJ_FEATHER)
        )

        # Extract raw numpy arrays once
        _sj_refs = sj_df["ref"].values
        _sj_starts = sj_df["start"].values
        _sj_ends = sj_df["end"].values
        _sj_strands = sj_df["strand"].values
        _sj_tidxs = sj_df["t_index"].values

        # Exact-match map: vectorised sort + boundary detection
        logger.debug("Building splice junction exact-match map")
        _sj_uref, _sj_rc = np.unique(_sj_refs, return_inverse=True)
        _sj_order = np.lexsort((_sj_strands, _sj_ends, _sj_starts, _sj_rc))
        _sj_rc_s = _sj_rc[_sj_order]
        _sj_st_s = _sj_starts[_sj_order]
        _sj_en_s = _sj_ends[_sj_order]
        _sj_sd_s = _sj_strands[_sj_order]
        _sj_ti_s = _sj_tidxs[_sj_order]
        _sj_rf_s = _sj_refs[_sj_order]

        _sj_n = len(_sj_order)
        _sj_changed = np.empty(_sj_n, dtype=bool)
        if _sj_n > 0:
            _sj_changed[0] = True
        if _sj_n > 1:
            _sj_changed[1:] = (
                (np.diff(_sj_rc_s) != 0)
                | (np.diff(_sj_st_s) != 0)
                | (np.diff(_sj_en_s) != 0)
                | (np.diff(_sj_sd_s) != 0)
            )
        _sj_gstarts = np.where(_sj_changed)[0]
        _sj_gends = np.append(_sj_gstarts[1:], _sj_n)

        # Pre-convert to Python lists for fast iteration
        _sj_rf_list = _sj_rf_s.tolist()
        _sj_st_list = _sj_st_s.tolist()
        _sj_en_list = _sj_en_s.tolist()
        _sj_sd_list = _sj_sd_s.tolist()
        _sj_ti_list = _sj_ti_s.tolist()
        _sj_gs_list = _sj_gstarts.tolist()
        _sj_ge_list = _sj_gends.tolist()

        sj_map: dict[tuple, frozenset[int]] = {}
        for i in range(len(_sj_gs_list)):
            s = _sj_gs_list[i]
            e = _sj_ge_list[i]
            key = (_sj_rf_list[s], _sj_st_list[s],
                   _sj_en_list[s], _sj_sd_list[s])
            sj_map[key] = frozenset(_sj_ti_list[s:e])
        self.sj_map = sj_map

        logger.debug(f"Splice junctions: {len(self.sj_map)} unique, "
                     f"{len(sj_df)} total")

        # -- C++ FragmentResolver (native fragment resolution) ------------------
        from .native import FragmentResolver
        ctx = FragmentResolver()

        # 1. Overlap index from collapsed data
        len(_collapse_keys)
        cr_refs = [k[0] for k in _collapse_keys]
        cr_starts = [k[1] for k in _collapse_keys]
        cr_ends = [k[2] for k in _collapse_keys]
        # CSR for transcript sets
        tset_flat: list[int] = []
        tset_offsets: list[int] = [0]
        for ts in iv_t_set:
            tset_flat.extend(sorted(ts))
            tset_offsets.append(len(tset_flat))
        ctx.build_overlap_index(
            cr_refs,
            cr_starts,
            cr_ends,
            iv_type,
            tset_flat,
            tset_offsets,
        )

        # 2. SJ exact-match map
        sj_refs_l: list[str] = []
        sj_starts_l: list[int] = []
        sj_ends_l: list[int] = []
        sj_strands_l: list[int] = []
        sj_t_flat: list[int] = []
        sj_t_offsets: list[int] = [0]
        for (ref, start, end, strand), tset in self.sj_map.items():
            sj_refs_l.append(ref)
            sj_starts_l.append(start)
            sj_ends_l.append(end)
            sj_strands_l.append(strand)
            sj_t_flat.extend(sorted(tset))
            sj_t_offsets.append(len(sj_t_flat))
        ctx.build_sj_map(
            sj_refs_l,
            sj_starts_l,
            sj_ends_l,
            sj_strands_l,
            sj_t_flat,
            sj_t_offsets,
        )

        # 3. Per-transcript exon CSR for transcript-space FL computation
        n_t = len(self.t_to_g_arr)
        exon_offsets = [0] * (n_t + 1)
        for t_idx, ivs in (self._t_exon_intervals or {}).items():
            exon_offsets[t_idx + 1] = len(ivs)
        for i in range(n_t):
            exon_offsets[i + 1] += exon_offsets[i]
        total_exons = exon_offsets[n_t]
        exon_starts_flat = [0] * total_exons
        exon_ends_flat = [0] * total_exons
        exon_cumsum_flat = [0] * total_exons
        t_lengths = [0] * n_t
        for t_idx, ivs in (self._t_exon_intervals or {}).items():
            off = exon_offsets[t_idx]
            cum = 0
            for j in range(len(ivs)):
                s, e = int(ivs[j, 0]), int(ivs[j, 1])
                exon_starts_flat[off + j] = s
                exon_ends_flat[off + j] = e
                exon_cumsum_flat[off + j] = cum
                cum += e - s
            t_lengths[t_idx] = cum
        ctx.build_exon_index(
            exon_offsets,
            exon_starts_flat,
            exon_ends_flat,
            exon_cumsum_flat,
            t_lengths,
        )

        # 5. Metadata
        ctx.set_metadata(
            self.t_to_g_arr.astype(np.int32).tolist(),
            len(self.t_to_g_arr),
        )

        self.resolver = ctx
        logger.debug("Built native FragmentResolver for C++ resolution")

        return self

    # -- query methods --------------------------------------------------------

    def query(
        self, exon: GenomicInterval,
    ) -> list[tuple[int, int, int, frozenset[int]]]:
        """Query the collapsed interval index with an aligned exon block.

        Returns one entry per unique ``(ref, start, end, type)`` boundary
        that overlaps the query.  The transcript set is already
        pre-collapsed into a ``frozenset[int]``.

        Parameters
        ----------
        exon : GenomicInterval
            An aligned exon block from a ``Fragment``.

        Returns
        -------
        list[tuple[int, int, int, frozenset[int]]]
            List of ``(hit_start, hit_end, interval_type, t_set)``
            tuples.  ``t_set`` is empty for INTERGENIC intervals.
        """
        hits: list[tuple[int, int, int, frozenset[int]]] = []
        for h_start, h_end, label in self.cr.overlap(exon.ref, exon.start, exon.end):
            hits.append((
                h_start,
                h_end,
                self._iv_type[label],
                self._iv_t_set[label],
            ))
        return hits

    def get_exon_intervals(self, t_idx: int) -> np.ndarray | None:
        """Return sorted exon ``[start, end)`` intervals for a transcript.

        Parameters
        ----------
        t_idx : int
            Global transcript index.

        Returns
        -------
        np.ndarray or None
            ``(n_exons, 2)`` int32 array sorted by genomic start, or
            ``None`` if the transcript has no cached exon intervals.
        """
        if self._t_exon_intervals is None:
            return None
        return self._t_exon_intervals.get(t_idx)

    def build_exon_csr(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Build CSR arrays for per-transcript exon positions.

        Converts the per-transcript exon dict into flat CSR arrays
        suitable for direct C++ consumption, eliminating the need for
        a Python dict → C++ dict-unpacking round-trip.

        Returns
        -------
        offsets : np.ndarray
            int32[num_transcripts + 1] — CSR offsets.
        starts : np.ndarray
            int32[total_exons] — exon start positions.
        ends : np.ndarray
            int32[total_exons] — exon end positions.
        cumsum_before : np.ndarray
            int32[total_exons] — cumulative exon length before each exon.
        """
        n_t = self.num_transcripts
        offsets = np.zeros(n_t + 1, dtype=np.int32)
        empty = np.empty(0, dtype=np.int32)

        if self._t_exon_intervals is None or len(self._t_exon_intervals) == 0:
            return offsets, empty, empty, empty

        # Pass 1: count exons per transcript
        for t_idx, ivs in self._t_exon_intervals.items():
            offsets[t_idx + 1] = len(ivs)
        np.cumsum(offsets, out=offsets)
        total = int(offsets[n_t])

        starts = np.empty(total, dtype=np.int32)
        ends = np.empty(total, dtype=np.int32)
        cumsum_before = np.empty(total, dtype=np.int32)

        # Pass 2: fill exon data
        for t_idx, ivs in self._t_exon_intervals.items():
            off = offsets[t_idx]
            n = len(ivs)
            s = ivs[:, 0]
            e = ivs[:, 1]
            starts[off : off + n] = s
            ends[off : off + n] = e
            lengths = e - s
            cumsum_before[off] = 0
            if n > 1:
                np.cumsum(lengths[:-1], out=cumsum_before[off + 1 : off + n])

        return offsets, starts, ends, cumsum_before

