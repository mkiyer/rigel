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
from typing import Iterator, Literal

import cgranges
import numpy as np
import pandas as pd
import pysam

from .types import GenomicInterval, IntervalType, RefInterval
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


def build_splice_junctions(transcripts: list[Transcript]) -> pd.DataFrame:
    """Extract splice junctions (introns) from all transcripts.

    Each intron boundary within a transcript produces one SpliceJunction
    record. The resulting DataFrame is sorted by (ref, start, end, strand).
    """
    rows = [
        RefInterval(t.ref, start, end, t.strand,
                    IntervalType.SJ, t.t_index)
        for t in transcripts
        for start, end in t.introns()
    ]
    rows.sort(key=lambda sj: (sj.ref, sj.start, sj.end, sj.strand))
    return pd.DataFrame(rows, columns=RefInterval._fields)


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
            chrom = t.ref
            chrom_start = t.exons[0].start
            chrom_end = t.exons[-1].end
            name = t.t_id or "."
            score = 0
            strand = t.strand.to_str()
            if strand not in ("+", "-"):
                strand = "+"
            thick_start = chrom_start
            thick_end = chrom_end
            rgb = "0"
            block_count = len(t.exons)
            block_sizes = ",".join(
                str(e.end - e.start) for e in t.exons
            )
            block_starts = ",".join(
                str(e.start - chrom_start) for e in t.exons
            )
            fh.write(
                f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t"
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


def _gen_transcript_intervals(t: Transcript) -> Iterator[RefInterval]:
    """Yield exon and transcript-span intervals for a single transcript.

    Each exon produces an EXON interval.  One TRANSCRIPT interval spans
    the full transcript ``[first_exon.start, last_exon.end)``.
    """
    # Exons
    for e in t.exons:
        yield RefInterval(t.ref, e.start, e.end, t.strand,
                          IntervalType.EXON, t.t_index)
    # One transcript span (replaces per-gap INTRON intervals)
    if t.exons:
        yield RefInterval(t.ref, t.exons[0].start, t.exons[-1].end,
                          t.strand, IntervalType.TRANSCRIPT, t.t_index)


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
    matching the ``RefInterval`` NamedTuple fields.
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

        # Unified cgranges index (collapsed EXON + TRANSCRIPT + INTERGENIC)
        # Label = index into _iv_type and _iv_t_set lookup lists.
        self.cr: cgranges.cgranges | None = None
        self._iv_type: list[int] | None = None
        self._iv_t_set: list[frozenset[int]] | None = None

        # Splice-junction exact-match map  (ref, start, end, strand) → frozenset[int]
        self.sj_map: dict | None = None

        # Splice-junction cgranges index (for gap containment queries / fragment length)
        # Label = row index into the _sj_* lookup arrays.
        self.sj_cr: cgranges.cgranges | None = None
        self._sj_t_index: np.ndarray | None = None
        self._sj_strand: np.ndarray | None = None

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
        gtf_parse_mode : {"strict", "warn-skip"}
            GTF parsing behavior. ``"strict"`` (default) fails fast on
            malformed lines; ``"warn-skip"`` logs warnings and skips.
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
        logger.debug("Building collapsed interval index")
        _collapse: dict[tuple, tuple[int, set[int]]] = {}
        for row in iv_df.itertuples(index=False):
            key = (row.ref, row.start, row.end, int(row.interval_type))
            if key not in _collapse:
                _collapse[key] = (int(row.interval_type), set())
            t_idx = int(row.t_index)
            if t_idx >= 0:
                _collapse[key][1].add(t_idx)

        cr = cgranges.cgranges()
        iv_type: list[int] = []
        iv_t_set: list[frozenset[int]] = []
        for label, ((ref, start, end, _itype), (itype_val, tset)) in enumerate(
            _collapse.items()
        ):
            cr.add(ref, start, end, label)
            iv_type.append(itype_val)
            iv_t_set.append(frozenset(tset))
        cr.index()
        self.cr = cr
        self._iv_type = iv_type
        self._iv_t_set = iv_t_set
        logger.debug(f"Interval index: {len(iv_df)} rows → {len(iv_type)} collapsed")

        # -- per-transcript exon intervals for coverage-weight model ----------
        exon_mask = iv_df["interval_type"].values == int(IntervalType.EXON)
        exon_rows = iv_df.loc[exon_mask, ["t_index", "start", "end"]]
        t_exon_intervals: dict[int, np.ndarray] = {}
        for t_idx, grp in exon_rows.groupby("t_index"):
            if int(t_idx) < 0:
                continue
            coords = grp[["start", "end"]].values
            coords = coords[coords[:, 0].argsort()]
            t_exon_intervals[int(t_idx)] = coords.astype(np.int32)
        self._t_exon_intervals = t_exon_intervals
        logger.debug(f"Cached exon intervals for {len(t_exon_intervals)} transcripts")

        # -- splice junction indexes ------------------------------------------
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(
            os.path.join(index_dir, SJ_FEATHER)
        )

        # Exact-match map for counting (ref, start, end, strand) → frozenset[int]
        logger.debug("Building splice junction exact-match map")
        sj_map: dict[tuple, set[int]] = {}
        for row in sj_df.itertuples(index=False):
            key = (row.ref, row.start, row.end, row.strand)
            if key not in sj_map:
                sj_map[key] = set()
            sj_map[key].add(int(row.t_index))
        self.sj_map = {
            key: frozenset(val) for key, val in sj_map.items()
        }

        # cgranges tree for containment / gap queries (fragment length)
        logger.debug("Building splice junction cgranges index")
        sj_cr = cgranges.cgranges()
        for label, row in enumerate(sj_df.itertuples(index=False)):
            sj_cr.add(row.ref, row.start, row.end, label)
        sj_cr.index()
        self.sj_cr = sj_cr
        self._sj_t_index = sj_df["t_index"].values
        self._sj_strand = sj_df["strand"].values
        logger.debug(f"Splice junctions: {len(self.sj_map)} unique, "
                     f"{len(sj_df)} total")
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

    def query_gap_sjs(
        self, ref: str, start: int, end: int
    ) -> list[tuple[int, int, int, int]]:
        """Find annotated splice junctions fully contained in a gap region.

        Used for fragment-length computation: given the gap between paired-end
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
        list[tuple[int, int, int, int]]
            List of (t_index, strand, sj_start, sj_end) tuples
            where sj_start >= start and sj_end <= end.
        """
        hits: list[tuple[int, int, int, int]] = []
        for h_start, h_end, label in self.sj_cr.overlap(ref, start, end):
            if h_start >= start and h_end <= end:
                hits.append((
                    int(self._sj_t_index[label]),
                    int(self._sj_strand[label]),
                    h_start,
                    h_end,
                ))
        return hits

    # -- gene-level exon/intron geometry ------------------------------------

    def exon_lengths_per_gene(self) -> list[list[int]]:
        """Merged exon lengths per gene for eCDF effective length.

        For each gene, merges overlapping exonic intervals across all
        transcripts to produce the true exonic footprint, then returns
        the list of individual (non-overlapping) exon lengths.

        Returns
        -------
        list[list[int]]
            ``result[g_index]`` is a list of merged exon lengths for
            gene *g_index*.
        """
        n_genes = self.num_genes
        # Collect per-gene exon intervals from interval table
        gene_exon_ivs: list[list[tuple[int, int]]] = [[] for _ in range(n_genes)]

        # Read exon intervals from the persisted interval file
        if self.index_dir is not None:
            iv_df = pd.read_feather(
                os.path.join(self.index_dir, INTERVALS_FEATHER)
            )
            exon_mask = iv_df["interval_type"].values == int(IntervalType.EXON)
            t_indices = iv_df["t_index"].values[exon_mask]
            starts = iv_df["start"].values[exon_mask]
            ends = iv_df["end"].values[exon_mask]
            valid = t_indices >= 0
            for t_idx, s, e in zip(
                t_indices[valid], starts[valid], ends[valid]
            ):
                g_idx = int(self.t_to_g_arr[int(t_idx)])
                if 0 <= g_idx < n_genes:
                    gene_exon_ivs[g_idx].append((int(s), int(e)))
        else:
            # Fallback: use transcript table with uniform exon length
            for g_idx in range(n_genes):
                mask = self.t_to_g_arr == g_idx
                if mask.any():
                    lengths = self.t_df.loc[mask, "length"].values
                    max_len = int(lengths.max())
                    gene_exon_ivs[g_idx].append((0, max_len))

        # Merge overlapping intervals per gene
        result: list[list[int]] = []
        for g_idx in range(n_genes):
            ivs = gene_exon_ivs[g_idx]
            if not ivs:
                result.append([])
                continue
            # Sort by start
            ivs.sort()
            merged: list[tuple[int, int]] = [ivs[0]]
            for s, e in ivs[1:]:
                if s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            result.append([e - s for s, e in merged])
        return result

    def intron_span_per_gene(self) -> np.ndarray:
        """Intronic span per gene (gene_span minus merged exon length).

        Returns
        -------
        np.ndarray
            float64[num_genes] — intronic span for each gene.
        """
        gene_spans = (
            self.g_df["end"].values - self.g_df["start"].values
        ).astype(np.float64)
        exon_lens = self.exon_lengths_per_gene()
        total_exonic = np.array(
            [float(sum(lens)) for lens in exon_lens],
            dtype=np.float64,
        )
        return np.maximum(gene_spans - total_exonic, 0.0)
