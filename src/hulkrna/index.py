"""
hulkrna.index — Build, load, and query the hulkrna reference index.

The index is constructed from a genome FASTA (with samtools .fai) and a
GENCODE GTF annotation file. It produces four Feather files (with optional
TSV mirrors) in an output directory:

    ref_lengths.feather   — chromosome names and lengths
    transcripts.feather   — one row per transcript with integer indices
    intervals.feather     — exon/intron/intergenic tiling of the genome
    sj.feather            — annotated splice junctions from transcript introns

The ``HulkIndex`` class provides both the ``build()`` method for creating
the index and ``load()`` / query methods for using it during counting.
"""

import collections
import logging
import os
from pathlib import Path
from typing import Iterator

import cgranges
import numpy as np
import pandas as pd
import pysam

from .types import GenomicInterval, Interval, IntervalType, RefInterval, Strand
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


# ---------------------------------------------------------------------------
# Build helpers (public for testability)
# ---------------------------------------------------------------------------

def load_reference_lengths(fasta_file: str | Path) -> dict[str, int]:
    """Read chromosome names and lengths from a FASTA file.

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


def read_transcripts(gtf_file: str | Path) -> list[Transcript]:
    """Parse a GTF file and return a sorted list of Transcript objects.

    Each transcript is assigned a sequential ``t_index`` and a ``g_index``
    (unique integer per ``g_id``). Transcripts are sorted by
    ``(ref, start, end, strand)`` so that downstream interval generation
    can process them in genomic order.
    """
    transcripts = Transcript.read_gtf(str(gtf_file))

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


def build_splice_junctions(transcripts: list[Transcript]) -> pd.DataFrame:
    """Extract splice junctions (introns) from all transcripts.

    Each intron boundary within a transcript produces one SpliceJunction
    record. The resulting DataFrame is sorted by (ref, start, end, strand).
    """
    rows = [
        RefInterval(t.ref, start, end, t.strand,
                    IntervalType.SJ, t.t_index, t.g_index)
        for t in transcripts
        for start, end in t.introns()
    ]
    rows.sort(key=lambda sj: (sj.ref, sj.start, sj.end, sj.strand))
    return pd.DataFrame(rows, columns=RefInterval._fields)


def _gen_transcript_intervals(t: Transcript) -> Iterator[RefInterval]:
    """Yield exon and intron intervals for a single transcript."""
    # Exons
    for e in t.exons:
        yield RefInterval(t.ref, e.start, e.end, t.strand,
                          IntervalType.EXON, t.t_index, t.g_index)
    # Introns (gaps between consecutive exons)
    for e1, e2 in zip(t.exons[:-1], t.exons[1:]):
        yield RefInterval(t.ref, e1.end, e2.start, t.strand,
                          IntervalType.INTRON, t.t_index, t.g_index)


def _gen_genomic_intervals(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> Iterator[RefInterval]:
    """Tile the genome into exon, intron, and intergenic intervals.

    Transcripts must be sorted by ``(ref, start, end)``. Intergenic
    intervals fill the gaps between transcript clusters and extend to
    chromosome boundaries.
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
            # Entire chromosome is intergenic
            yield RefInterval(ref, 0, ref_length)
            continue

        end = 0  # tracks the rightmost coordinate reached
        cluster: list[Transcript] = []

        for t in t_list:
            if t.start > end:
                # Intergenic gap before this transcript cluster
                if end < t.start:
                    yield RefInterval(ref, end, t.start)
                # Emit intervals for the previous cluster
                for tc in cluster:
                    yield from _gen_transcript_intervals(tc)
                cluster = []

            end = max(end, t.end)
            cluster.append(t)

        # Emit final cluster
        for tc in cluster:
            yield from _gen_transcript_intervals(tc)

        # Intergenic gap to end of chromosome
        if end < ref_length:
            yield RefInterval(ref, end, ref_length)


def build_genomic_intervals(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> pd.DataFrame:
    """Build the complete interval table covering the genome.

    Returns a DataFrame sorted by (ref, start, end, strand) with columns
    matching the ``Interval`` NamedTuple fields.
    """
    intervals = list(_gen_genomic_intervals(transcripts, ref_lengths))
    intervals.sort(key=lambda iv: (iv.ref, iv.start, iv.end, iv.strand))
    return pd.DataFrame(intervals, columns=RefInterval._fields)


# ---------------------------------------------------------------------------
# HulkIndex — unified index class
# ---------------------------------------------------------------------------

class HulkIndex:
    """In-memory reference index for fast transcript/gene overlap queries.

    Use ``HulkIndex.build()`` to create the on-disk index from a FASTA +
    GTF, then ``HulkIndex.load()`` to read it back for counting.

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

        # Unified cgranges index (EXON + INTRON + INTERGENIC intervals)
        # Label = row index into the _iv_* lookup arrays.
        self.cr: cgranges.cgranges | None = None
        self._iv_t_index: np.ndarray | None = None
        self._iv_g_index: np.ndarray | None = None
        self._iv_type: np.ndarray | None = None

        # Splice-junction exact-match map  (ref, start, end, strand) → (t_set, g_set)
        self.sj_map: dict | None = None

        # Splice-junction cgranges index (for gap containment queries / insert size)
        # Label = row index into the _sj_* lookup arrays.
        self.sj_cr: cgranges.cgranges | None = None
        self._sj_t_index: np.ndarray | None = None
        self._sj_g_index: np.ndarray | None = None
        self._sj_strand: np.ndarray | None = None

    # -- properties -----------------------------------------------------------

    @property
    def num_transcripts(self) -> int:
        return len(self.t_df)

    @property
    def num_genes(self) -> int:
        return len(self.g_df)

    # -- build (static) -------------------------------------------------------

    @staticmethod
    def build(
        fasta_file: str | Path,
        gtf_file: str | Path,
        output_dir: str | Path,
        *,
        feather_compression: str = "lz4",
        write_tsv: bool = True,
    ) -> None:
        """Build the hulkrna reference index and write to disk.

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
        transcripts = read_transcripts(gtf_file)
        logger.info(f"[DONE] Read {len(transcripts)} transcripts")

        t_df = transcripts_to_dataframe(transcripts)
        t_df.to_feather(output_dir / TRANSCRIPTS_FEATHER, **feather_kwargs)
        if write_tsv:
            t_df.to_csv(output_dir / TRANSCRIPTS_TSV, sep="\t", index=False)

        # -- Splice junctions -------------------------------------------------
        logger.info("[START] Building splice junctions")
        sj_df = build_splice_junctions(transcripts)
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
    def load(cls, index_dir: str | Path) -> "HulkIndex":
        """Load an index from the Feather files in *index_dir*.

        Builds cgranges interval trees and splice-junction lookup maps
        for fast overlap queries during counting.
        """
        index_dir = str(index_dir)
        self = cls()
        self.index_dir = index_dir

        # -- transcripts ------------------------------------------------------
        self.t_df = pd.read_feather(
            os.path.join(index_dir, TRANSCRIPTS_FEATHER)
        )
        assert (self.t_df.index == self.t_df["t_index"]).all()

        # -- gene table (derived) ---------------------------------------------
        self.g_df = cls._build_gene_table(self.t_df)
        assert (self.g_df.index == self.g_df["g_index"]).all()

        # -- fast numpy lookup arrays -----------------------------------------
        self.t_to_g_arr = self.t_df["g_index"].values
        self.t_to_strand_arr = self.t_df["strand"].values
        self.g_to_strand_arr = self.g_df["strand"].values

        # -- interval index (unified cgranges) ---------------------------------
        logger.debug("Reading intervals")
        iv_df = pd.read_feather(
            os.path.join(index_dir, INTERVALS_FEATHER)
        )
        logger.debug("Building unified interval index")
        cr = cgranges.cgranges()
        for label, row in enumerate(iv_df.itertuples(index=False)):
            cr.add(row.ref, row.start, row.end, label)
        cr.index()
        self.cr = cr
        self._iv_t_index = iv_df["t_index"].values
        self._iv_g_index = iv_df["g_index"].values
        self._iv_type = iv_df["interval_type"].values
        logger.debug(f"Interval index size: {len(iv_df)}")

        # -- splice junction indexes ------------------------------------------
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(
            os.path.join(index_dir, SJ_FEATHER)
        )

        # Exact-match map for counting (ref, start, end, strand) → sets
        logger.debug("Building splice junction exact-match map")
        sj_map: dict[tuple, tuple[set, set]] = {}
        for ref, start, end, strand, _ivtype, t_index, g_index in sj_df[
            list(RefInterval._fields)
        ].itertuples(index=False):
            key = (ref, start, end, strand)
            if key not in sj_map:
                sj_map[key] = (set(), set())
            sj_map[key][0].add(t_index)
            sj_map[key][1].add(g_index)
        self.sj_map = {
            key: (frozenset(val[0]), frozenset(val[1]))
            for key, val in sj_map.items()
        }

        # cgranges tree for containment / gap queries (insert size)
        logger.debug("Building splice junction cgranges index")
        sj_cr = cgranges.cgranges()
        for label, row in enumerate(sj_df.itertuples(index=False)):
            sj_cr.add(row.ref, row.start, row.end, label)
        sj_cr.index()
        self.sj_cr = sj_cr
        self._sj_t_index = sj_df["t_index"].values
        self._sj_g_index = sj_df["g_index"].values
        self._sj_strand = sj_df["strand"].values
        logger.debug(f"Splice junctions: {len(self.sj_map)} unique, "
                     f"{len(sj_df)} total")
        return self

    # -- query methods --------------------------------------------------------

    def query_exon(self, exon: GenomicInterval) -> list[tuple[int, int, int]]:
        """Query the unified interval index with an aligned exon block.

        Intersects the single cgranges tree and returns every reference
        interval (EXON, INTRON, or INTERGENIC) that overlaps the query.

        Because the index tiles the entire genome, a valid query is
        guaranteed to return at least one hit.

        Parameters
        ----------
        exon : GenomicInterval
            An aligned exon block from a ``Fragment``.

        Returns
        -------
        list[tuple[int, int, int]]
            List of (t_index, g_index, interval_type) tuples.
            t_index and g_index are -1 for intergenic intervals.
        """
        hits: list[tuple[int, int, int]] = []
        for _h_start, _h_end, label in self.cr.overlap(exon.ref, exon.start, exon.end):
            hits.append((
                int(self._iv_t_index[label]),
                int(self._iv_g_index[label]),
                int(self._iv_type[label]),
            ))
        return hits

    def query_gap_sjs(
        self, ref: str, start: int, end: int
    ) -> list[tuple[int, int, int, int, int]]:
        """Find annotated splice junctions fully contained in a gap region.

        Used for insert-size computation: given the gap between paired-end
        reads, find all annotated introns (splice junctions) that are
        fully contained within that gap so their lengths can be subtracted
        from the genomic distance.

        Parameters
        ----------
        ref : str
            Chromosome / reference name.
        start : int
            Left edge of the gap (0-based, inclusive).
        end : int
            Right edge of the gap (0-based, exclusive).

        Returns
        -------
        list[tuple[int, int, int, int, int]]
            List of (t_index, g_index, strand, sj_start, sj_end) tuples
            where sj_start >= start and sj_end <= end.
        """
        hits: list[tuple[int, int, int, int, int]] = []
        for h_start, h_end, label in self.sj_cr.overlap(ref, start, end):
            if h_start >= start and h_end <= end:
                hits.append((
                    int(self._sj_t_index[label]),
                    int(self._sj_g_index[label]),
                    int(self._sj_strand[label]),
                    h_start,
                    h_end,
                ))
        return hits
