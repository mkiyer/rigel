"""
rigel.index — Build, load, and query the rigel reference index.

The index is constructed from a genome FASTA (with samtools .fai) and a
GENCODE GTF annotation file. It produces feather files (with optional
TSV mirrors) in an output directory:

    ref_lengths.feather        — reference names and lengths
    transcripts.feather        — one row per transcript with integer indices
    intervals.feather          — exon/intron/intergenic tiling of the genome
    regions.feather            — calibration region partition (INTERGENIC/INTRON/EXON)
    sj.feather                 — annotated splice junctions from transcript introns
    splice_blacklist.feather   — (optional) splice-artifact junctions derived
                                 from the alignable Zarr store

The ``TranscriptIndex`` class provides both the ``build()`` method for creating
the index and ``load()`` / query methods for using it during quantification.
"""

import collections
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Literal

from .native import cgranges as _cgranges_cls
import numpy as np
import pandas as pd
import pysam

from .types import GenomicInterval, IntervalType, AnnotatedInterval
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

NODES_FEATHER = "nodes.feather"
NODES_TSV = "nodes.tsv"

EDGES_FEATHER = "edges.feather"
EDGES_TSV = "edges.tsv"

SJ_FEATHER = "sj.feather"
SJ_TSV = "sj.tsv"

SJ_BLACKLIST_FEATHER = "splice_blacklist.feather"
SJ_BLACKLIST_TSV = "splice_blacklist.tsv"

MANIFEST_JSON = "manifest.json"

#: On-disk index format version. Bumped whenever the schema or the
#: meaning of any persisted column changes. Loaders should refuse
#: indexes whose ``format_version`` they do not understand.
#:
#: Version history:
#:   2 — (legacy) intervals.feather + sj.feather + transcripts.feather.
#:   3 — adds regions.feather (calibration partition); ref_lengths.feather
#:        is now mandatory at load time and used to validate the region
#:        partition.
#:   4 — regions.feather stored fine-region signatures plus derived coarse
#:        bridge columns for calibration consumers.
#:   5 — calibration-v6: regions.feather is the minimal merged-signature
#:        partition [region_id, ref_name, start, end, length, signature];
#:        the derived coarse class is recomputed on load from `signature`.
#:   6 — three-channel calibration: regions.feather gained per-strand
#:        `mature_eligible_{pos,neg}` and boundaries.feather carried
#:        per-boundary annotation flags (is_tss/is_tes/is_splice_junction/
#:        genomic_sj_strand).
#:   7 — those v6 precompute columns removed (the message-precision collapse
#:        retired the mature/nascent overlay that consumed them; the solver
#:        reads junction strand from the accumulator motif instead).
#:        regions.feather is again the minimal partition [region_id, ref_name,
#:        start, end, length, signature]; boundaries.feather is [boundary_id,
#:        ref_name, position].
#:   8 — the SPLICE GRAPH replaces the region/boundary partition. nodes.feather +
#:        edges.feather are the only partition artifacts and both are MANDATORY;
#:        regions.feather / boundaries.feather are gone. The scanner is fed the node
#:        cut array by `calibration.splice_graph.build_node_partition_arrays`.
#:        Adjacent nodes may share a signature — that is the point: the merge it
#:        replaces deleted the cut at 53.4 % of real human transcript termini, so the
#:        partition could not see them at all.
INDEX_FORMAT_VERSION = 8


def _rigel_version() -> str:
    try:
        from . import __version__  # type: ignore[attr-defined]

        return str(__version__)
    except Exception:  # pragma: no cover
        return "unknown"


#: Read size for streaming a content digest. An I/O buffer, not a model parameter — it changes how fast
#: a digest is computed and never what the digest is.
_DIGEST_CHUNK_BYTES = 1 << 20


def _sha256_of_file(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_DIGEST_CHUNK_BYTES):
            h.update(chunk)
    return h.hexdigest()


def _sha256_of_directory(path: Path) -> str:
    """Digest a directory source — an unpackaged ``.zarr`` store — over its sorted relative paths.

    Traversal is ``sorted(rglob)``, so the digest is a property of the tree's contents and layout and
    never of the filesystem's iteration order.
    """
    import hashlib

    h = hashlib.sha256()
    for entry in sorted(p for p in path.rglob("*") if p.is_file()):
        h.update(str(entry.relative_to(path)).encode())
        h.update(_sha256_of_file(entry).encode())
    return h.hexdigest()


def source_record(path: str | Path) -> dict:
    """Provenance for ONE build input: where it was, how big it was, and what it contained.

    ⚠ **This is not the stale-cache trap.** That one forbids storing a hash of an artifact
    *beside* that artifact, because the two drift and the stale hash then verifies clean. This hashes an
    input the index **cannot recompute from itself** — the genome and the annotation are not in the index
    — so it is provenance, not a cache key. ``partition_hash`` and ``graph_hash`` stay computed on demand.

    ⭐ ``sha256`` rather than the ``blake2b`` used for the two partition keys, deliberately: those are
    internal cache keys nobody types, while this one is meant to be checked against ``shasum -a 256`` on
    the command line by whoever is trying to find the source again.
    """
    resolved = Path(path).resolve()
    if resolved.is_dir():
        files = [p for p in resolved.rglob("*") if p.is_file()]
        return {
            "path": str(resolved),
            "bytes": sum(p.stat().st_size for p in files),
            "sha256": _sha256_of_directory(resolved),
        }
    return {
        "path": str(resolved),
        "bytes": resolved.stat().st_size,
        "sha256": _sha256_of_file(resolved),
    }


def load_manifest(index_dir: str | Path) -> dict | None:
    """Load ``manifest.json`` from an index directory.

    Returns ``None`` for legacy indexes (no manifest file present).
    """
    path = Path(index_dir) / MANIFEST_JSON
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh)


def load_ref_lengths(path: str | Path) -> dict[str, int]:
    """Read ``ref_lengths.feather`` into an insertion-ordered dict.

    Iteration order matches the on-disk row order, which is the canonical ``ref_id`` assignment —
    the reference-id space the resolver, the scanner and the graph all share.
    """
    df = pd.read_feather(str(path))
    return {str(r): int(L) for r, L in zip(df["ref"], df["length"], strict=True)}


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
        ref_lengths = {fh.references[i]: fh.lengths[i] for i in range(fh.nreferences)}
    return ref_lengths


def _duplicate_transcript_groups(transcripts: list[Transcript]) -> dict[tuple, list[int]]:
    """Group transcripts by identical exon structure.

    Two transcripts are duplicates when they share the same reference sequence,
    strand, and bit-identical sorted ``(start, end)`` exon tuple. Such transcripts
    are mathematically unidentifiable from any read data. Transcripts that share an
    intron chain but differ in 5'/3' UTR boundaries are *not* duplicates.

    Returns ``{(ref, strand, exon_tuple): [transcript list-index, ...]}`` containing
    only groups with two or more members. Exon-less transcripts are ignored.
    """
    from collections import defaultdict

    sig_to_idx: dict[tuple, list[int]] = defaultdict(list)
    for i, t in enumerate(transcripts):
        if not t.exons:
            continue
        exon_tuple = tuple(sorted((iv.start, iv.end) for iv in t.exons))
        sig_to_idx[(t.ref, int(t.strand), exon_tuple)].append(i)

    return {k: v for k, v in sig_to_idx.items() if len(v) > 1}


def _check_no_duplicate_transcripts(transcripts: list[Transcript]) -> None:
    """Raise ``ValueError`` if two transcripts share an identical exon structure.

    Duplicate transcripts (see :func:`_duplicate_transcript_groups`) are
    mathematically unidentifiable from any read data and must be resolved before
    indexing — either by collapsing them in the GTF upstream, or by passing
    ``--collapse-duplicate-transcripts`` to let rigel drop them automatically.
    """
    duplicates = _duplicate_transcript_groups(transcripts)
    if not duplicates:
        return

    examples = []
    for (ref, strand, _exons), idxs in list(duplicates.items())[:5]:
        tids = [transcripts[i].t_id or f"<t_index={transcripts[i].t_index}>" for i in idxs]
        examples.append(f"  ({ref}, strand={strand}): {tids}")
    n_groups = len(duplicates)
    n_tx = sum(len(v) for v in duplicates.values())
    raise ValueError(
        f"GTF contains {n_groups} group(s) totalling {n_tx} transcript(s) "
        f"with identical exon coordinates. Such transcripts are "
        f"mathematically unidentifiable and must be collapsed before indexing.\n"
        f"Pass --collapse-duplicate-transcripts to have rigel drop them "
        f"automatically (keeping the lexicographically-smallest transcript ID per "
        f"group), or collapse them in the GTF upstream.\n"
        f"First {min(5, n_groups)} duplicate group(s):\n" + "\n".join(examples)
    )


def _collapse_duplicate_transcripts(transcripts: list[Transcript]) -> list[Transcript]:
    """Drop duplicate transcripts, keeping the lexicographically-smallest ID per group.

    Duplicate groups are identified exactly as in :func:`_check_no_duplicate_transcripts`
    (same ref, strand, and sorted exon ``(start, end)`` tuple). Because identical-
    coordinate transcripts are mathematically unidentifiable in quantification,
    collapsing each group to a single representative is loss-free for abundance
    estimation. Returns the filtered list (original order preserved); logs a summary of
    what was dropped. When there are no duplicates the input list is returned unchanged.
    """
    duplicates = _duplicate_transcript_groups(transcripts)
    if not duplicates:
        return transcripts

    drop_idx: set[int] = set()
    examples: list[str] = []
    for (ref, strand, _exons), idxs in duplicates.items():
        keep = min(idxs, key=lambda i: transcripts[i].t_id or "")
        dropped = [i for i in idxs if i != keep]
        drop_idx.update(dropped)
        if len(examples) < 5:
            examples.append(
                f"  ({ref}, strand={strand}): kept {transcripts[keep].t_id!r}, "
                f"dropped {[transcripts[i].t_id for i in dropped]}"
            )

    logger.warning(
        "Collapsed %d duplicate transcript(s) across %d group(s) with identical exon "
        "coordinates (kept the lexicographically-smallest transcript ID per group).\n"
        "First %d group(s):\n%s",
        len(drop_idx),
        len(duplicates),
        min(5, len(duplicates)),
        "\n".join(examples),
    )
    return [t for i, t in enumerate(transcripts) if i not in drop_idx]


def read_transcripts(
    gtf_file: str | Path,
    *,
    gtf_parse_mode: Literal["strict", "warn-skip"] = "strict",
    collapse_duplicate_transcripts: bool = False,
) -> list[Transcript]:
    """Parse a GTF file and return a sorted list of Transcript objects.

    Each transcript is assigned a sequential ``t_index`` and a ``g_index``
    (unique integer per ``g_id``). Transcripts are sorted by
    ``(ref, start, end, strand)`` so that downstream interval generation
    can process them in genomic order.

    Parameters
    ----------
    collapse_duplicate_transcripts : bool
        Controls handling of transcripts with identical exon structure
        (same ``(ref, strand, sorted exon coordinates)``). When ``False``
        (default) such duplicates raise ``ValueError``. When ``True`` they
        are collapsed to a single representative — the lexicographically-
        smallest transcript ID per group — via
        :func:`_collapse_duplicate_transcripts`.

    Raises
    ------
    ValueError
        If ``collapse_duplicate_transcripts`` is ``False`` and the GTF
        contains two or more transcripts with identical
        ``(ref, strand, sorted exon coordinates)``. See
        :func:`_check_no_duplicate_transcripts`.
    """
    transcripts = Transcript.read_gtf(
        str(gtf_file),
        parse_mode=gtf_parse_mode,
    )

    if collapse_duplicate_transcripts:
        transcripts = _collapse_duplicate_transcripts(transcripts)
    else:
        _check_no_duplicate_transcripts(transcripts)

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
) -> tuple[list[Transcript], dict[int, tuple], dict[tuple, int], dict[tuple, "Transcript"]]:
    """Create synthetic nRNA transcripts, detect annotated equivalents.

    This function implements the unified nRNA architecture:

    1. Collect (ref, strand, start, end) for each **multi-exon** transcript.
    2. Cluster nearby TSS and TES within *tolerance* bp independently.
       The merged span uses the **outer envelope** (min start, max end)
       so synthetics are never smaller than any contributing transcript.
    3. Check if an existing **single-exon** transcript already covers
       each merged span (exact match or full containment).
       If so, mark it ``is_nrna = True`` — no synthetic needed.
    4. For uncovered spans, create a synthetic single-exon ``Transcript``
       flagged ``is_nrna = True, is_synthetic = True``.
    5. Synthetic transcripts are **gene-neutral** by design: they receive a
       sentinel ``g_id`` equal to their own ``t_id`` and ``g_type =
       "nascent_rna"``.  Rigel is transcript-centric, and a merged nRNA span
       may be contributed to by transcripts from more than one annotated
       gene — inheriting one contributor's gene metadata would bake an
       arbitrary projection into the core index.  Gene-oriented summaries
       are derived downstream from annotated contributors only.

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
        Callers must assign ``t_index`` and ``g_index`` to each returned
        transcript.
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
    # Track one representative transcript per span only to preserve stable
    # iteration order across hash-table rebuilds.  Synthetics are gene-neutral
    # (see docstring) so we do NOT inherit any gene metadata from the rep.
    merged_spans: dict[tuple, Transcript] = {}  # key → first-seen transcript
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
    se_by_loc: dict[tuple[str, int], list[tuple[int, int, Transcript]]] = collections.defaultdict(
        list
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
                s_tx.is_nrna = True
                s_tx.nrna_n_contributors = span_contributor_count[span_key]
                covered.add(span_key)
                covered_equiv[span_key] = s_tx
                break

    # -- Step 5: Create synthetic transcripts ---------------------------------
    synthetics: list[Transcript] = []
    span_to_syn_idx: dict[tuple, int] = {}  # span_key → index in synthetics list
    for span_key in merged_spans:
        if span_key in covered:
            continue
        ref, strand_int, s, e = span_key
        t_id = f"RIGEL_NRNA_{ref}_{strand_int}_{s}_{e}"
        # Synthetic rows are gene-neutral. The sentinel g_id ( = t_id) keeps
        # the g_id column non-null so downstream categorical/string code works
        # unchanged; a dedicated g_index is allocated by the index builder.
        syn = Transcript(
            ref=ref,
            strand=Strand(strand_int),
            exons=[Interval(s, e)],
            length=e - s,
            t_id=t_id,
            g_id=t_id,
            g_name="",
            g_type="nascent_rna",
            is_nrna=True,
            is_synthetic=True,
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
        AnnotatedInterval(t.ref, start, end, t.strand, IntervalType.SJ, t.t_index)
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
            block_sizes = ",".join(str(e.end - e.start) for e in t.exons)
            block_starts = ",".join(str(e.start - ref_start) for e in t.exons)
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

    Synthetic nRNA transcripts are *not* indexed in cgranges -- their
    candidates are derived on-the-fly inside the C++ resolver from each
    contributor mRNA's ``nrna_t_index`` link.  Skipping them removes
    ~17% of the human cgranges interval count, eliminates the
    structural source of the INTRONIC-masking bug, and accelerates
    every overlap query.  Annotated nascent-equiv single-exon
    transcripts (``is_nrna=True, is_synthetic=False``) remain in
    cgranges as normal annotation.
    """
    if t.is_synthetic:
        return
    # Exons
    for e in t.exons:
        yield AnnotatedInterval(t.ref, e.start, e.end, t.strand, IntervalType.EXON, t.t_index)
    # One transcript span (replaces per-gap INTRON intervals)
    if t.exons:
        yield AnnotatedInterval(
            t.ref, t.exons[0].start, t.exons[-1].end, t.strand, IntervalType.TRANSCRIPT, t.t_index
        )


@dataclass(frozen=True, slots=True)
class _IntergenicSpan:
    """A reference span containing no real (non-synthetic) transcript.

    Half-open coordinates ``[start, end)``.
    """

    start: int
    end: int


@dataclass(frozen=True, slots=True)
class _GenicSpan:
    """A reference span covered by one connected (strand-agnostic) transcript cluster.

    ``transcripts`` is the cluster's transcript list in input order; it
    is strand-agnostic and excludes synthetic nRNAs (which are
    constructed downstream of the layout sweep).
    """

    start: int
    end: int
    transcripts: tuple = field(default_factory=tuple)


def _iter_reference_layout(
    ref_length: int,
    transcripts: list[Transcript],
) -> Iterator["_IntergenicSpan | _GenicSpan"]:
    """Yield interleaved ``_IntergenicSpan`` / ``_GenicSpan`` tiling [0, ref_length).

    ``transcripts`` is the per-reference list (already filtered to one
    reference), and **must be sorted by ``(start, end)``**. Synthetic
    transcripts are excluded from the sweep so they cannot coalesce
    real genic spans.

    Invariants:
      - Spans alternate (no two intergenic in a row) and tile the
        reference exactly with no gaps and no overlaps.
      - Each ``_GenicSpan`` contains at least one (non-synthetic)
        transcript with at least one exon.
      - If the reference contains no real transcripts, a single
        ``_IntergenicSpan(0, ref_length)`` is yielded (when
        ``ref_length > 0``).
    """
    real_ts = [t for t in transcripts if not t.is_synthetic]

    if not real_ts:
        if ref_length > 0:
            yield _IntergenicSpan(0, ref_length)
        return

    cursor = 0
    cluster_start = real_ts[0].start
    cluster_end = real_ts[0].end
    cluster: list[Transcript] = [real_ts[0]]

    for t in real_ts[1:]:
        if t.start > cluster_end:
            # Close current cluster.
            if cursor < cluster_start:
                yield _IntergenicSpan(cursor, cluster_start)
            yield _GenicSpan(cluster_start, cluster_end, tuple(cluster))
            cursor = cluster_end
            cluster_start = t.start
            cluster_end = t.end
            cluster = [t]
        else:
            if t.end > cluster_end:
                cluster_end = t.end
            cluster.append(t)

    # Final cluster.
    if cursor < cluster_start:
        yield _IntergenicSpan(cursor, cluster_start)
    yield _GenicSpan(cluster_start, cluster_end, tuple(cluster))
    if cluster_end < ref_length:
        yield _IntergenicSpan(cluster_end, ref_length)


def _emit_genomic_intervals(
    ref: str,
    layout: Iterator["_IntergenicSpan | _GenicSpan"],
) -> Iterator[AnnotatedInterval]:
    """Convert a per-reference layout stream into ``AnnotatedInterval`` rows.

    Each ``_IntergenicSpan`` becomes one INTERGENIC interval. Each
    ``_GenicSpan`` is decomposed via :func:`_gen_transcript_intervals`
    (per-exon EXON intervals plus one TRANSCRIPT span per transcript).
    Synthetic transcripts are filtered out by ``_gen_transcript_intervals``.
    """
    for span in layout:
        if isinstance(span, _IntergenicSpan):
            yield AnnotatedInterval(ref, span.start, span.end)
        else:
            assert isinstance(span, _GenicSpan)
            for tc in span.transcripts:
                yield from _gen_transcript_intervals(tc)


def _group_transcripts_by_ref(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> dict[str, list[Transcript]]:
    """Group transcripts by reference, validating each ref appears in ``ref_lengths``."""
    ref_transcripts: dict[str, list[Transcript]] = collections.defaultdict(list)
    for t in transcripts:
        if t.ref not in ref_lengths:
            raise ValueError(
                f"Transcript {t.t_id} has reference '{t.ref}' not found in the FASTA index"
            )
        ref_transcripts[t.ref].append(t)
    return ref_transcripts


def build_index_artifacts(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build ``intervals.feather`` and the v8 splice graph (``nodes`` + ``edges``).

    The per-reference layout (intergenic / genic spans) feeds ONLY the cgranges-style annotated
    interval table (:func:`_emit_genomic_intervals`). The calibration partition is built
    INDEPENDENTLY, by :func:`rigel.calibration.splice_graph.build_splice_graph`.

    Returns ``(intervals_df, nodes_df, edges_df)``, typed per the on-disk schemas. ``intervals_df``
    carries the columns of :class:`AnnotatedInterval`, sorted by ``(ref, start, end, strand)``.

    Transcripts must be sorted by ``(ref, start, end)``.
    """
    from .calibration.splice_graph import build_splice_graph, validate_graph

    ref_transcripts = _group_transcripts_by_ref(transcripts, ref_lengths)

    intervals: list[AnnotatedInterval] = []

    for ref, ref_length in ref_lengths.items():
        layout_for_intervals = list(
            _iter_reference_layout(ref_length, ref_transcripts.get(ref, []))
        )
        intervals.extend(_emit_genomic_intervals(ref, iter(layout_for_intervals)))

    intervals.sort(key=lambda iv: (iv.ref, iv.start, iv.end, iv.strand))
    iv_df = pd.DataFrame(intervals, columns=AnnotatedInterval._fields)
    nodes_df, edges_df = build_splice_graph(transcripts, ref_lengths)
    validate_graph(nodes_df, edges_df, ref_lengths, transcripts=transcripts)
    return iv_df, nodes_df, edges_df


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
        #: v8 splice graph, or ``None`` on an index built before it existed (plan W1a).
        self.nodes_df = None
        self.edges_df = None
        self.t_df: pd.DataFrame | None = None
        self.g_df: pd.DataFrame | None = None
        self.t_to_g_arr: np.ndarray | None = None
        self.t_to_strand_arr: np.ndarray | None = None
        self.t_to_ref_arr: np.ndarray | None = None
        self.g_to_strand_arr: np.ndarray | None = None

        # Unified cgranges index (collapsed EXON + TRANSCRIPT + INTERGENIC)
        # Label = index into _iv_type and _iv_t_set lookup lists.
        self.cr: _cgranges_cls | None = None
        self._iv_type: list[int] | None = None
        self._iv_t_set: list[frozenset[int]] | None = None

        # Splice-junction exact-match map  (ref, start, end, strand) → frozenset[int]
        self.sj_map: dict | None = None

        # Splice-artifact blacklist size, set at load():
        #   None → index predates the field / never loaded
        #   0    → no blacklist present (artifact detection is OFF)
        #   >0   → number of blacklisted junctions active (detection is ON)
        self.sj_blacklist_size: int | None = None

        # Per-transcript exon intervals for coverage-weight model.
        # Maps t_index → (n_exons, 2) int32 array of [start, end) intervals
        # sorted by genomic start position.
        self._t_exon_intervals: dict[int, np.ndarray] | None = None

        # When True, Python-side structures are kept after C++ projection
        # (for unit tests that inspect _iv_t_set, sj_map, _t_exon_intervals).
        self._retain_test_structures: bool = False

    # -- properties -----------------------------------------------------------

    @property
    def num_transcripts(self) -> int:
        return len(self.t_df)

    @property
    def num_genes(self) -> int:
        """Total rows in ``g_df`` (annotated genes + synthetic nRNA placeholders).

        This is the kernel-facing count used by ``_aggregate_to_gene``; it
        matches ``len(g_df)`` and preserves the ``g_index == row`` invariant.
        User-facing summaries should use :attr:`num_annotated_genes` instead.
        """
        return len(self.g_df)

    @property
    def num_annotated_genes(self) -> int:
        """Number of real (non-synthetic) gene rows in ``g_df``."""
        if "is_synthetic" not in self.g_df.columns:
            return len(self.g_df)
        return int((~self.g_df["is_synthetic"].to_numpy()).sum())

    @property
    def partition_hash(self) -> str:
        """16-hex-char content hash of **the partition the scanner actually sees** — the cache key.

        Hashes exactly the ``(boundary_positions, ref_pos_offsets, region_types)`` triple that
        :meth:`rigel.native.BamScanner.set_regions` receives, plus the reference lengths. Two indexes
        produce the same hash **iff** a scan against them yields an identically-shaped, identically-keyed
        accumulator payload — which is precisely the condition under which a cached payload is reusable.

        ⚠ **Computed on demand, never stored.** A hash written into ``manifest.json`` at build time is a
        derived value that can go stale against the feathers beside it; this one cannot.

        ⚠ **It covers ``nodes.feather`` only, and that is deliberate** — it is the key for a cached
        SCAN, and the scan sees the cut array and nothing else. ``edges.feather`` (the flags and
        reaches) can change without invalidating a payload, and does: the 2026-07-29 flag-filter fix
        rewrote every edge file while leaving every node file byte-identical. Anything that caches an
        *edge*-derived artifact must carry its own provenance; this hash will not catch it.

        Cost at human scale: ~60 ms (8.4 MB of int64 through blake2b plus the groupby that builds
        it) — negligible against the BAM scan it gates.
        """
        import hashlib

        from .calibration.splice_graph import build_node_partition_arrays

        h = hashlib.blake2b(digest_size=8)
        for arr in build_node_partition_arrays(self):
            a = np.ascontiguousarray(arr)
            h.update(str(a.dtype).encode())
            h.update(a.tobytes())
        for name, length in self.ref_lengths.items():
            h.update(f"{name}:{length}|".encode())
        return h.hexdigest()

    @property
    def graph_hash(self) -> str:
        """16-hex-char content hash of **everything the accumulator payload depends on** — nodes AND edges.

        ⚠ **`partition_hash` is not enough for a payload, and that is not an oversight in either of them.**
        That hash keys a cached *scan*, and a scan sees the cut array; this one keys a cached *tally*, whose
        junction axis is meaningless against a different junction CSR. ⭐ The two genuinely differ: the
        2026-07-29 flag fix rewrote every ``edges.feather`` while leaving every ``nodes.feather``
        byte-identical, so a nodes-only key would have verified **clean** against a stale payload and fed
        every downstream comparison the pre-fix junctions.

        So it is ``partition_hash``'s inputs plus the junction CSR — the donor offsets, the acceptor cut
        indices and the annotated strands, i.e. exactly what crosses into ``set_junctions``.

        ⚠ **Computed on demand, never stored.** A hash written beside the data it describes can go stale
        against it; this one cannot.
        """
        import hashlib

        from .calibration.splice_graph import build_junction_edge_arrays

        h = hashlib.blake2b(digest_size=8)
        h.update(self.partition_hash.encode())
        junctions = build_junction_edge_arrays(self)
        # ⚠ `edge_row` is deliberately absent: it is a join key back to `edges_df`, it never crosses the
        # ABI, and hashing it would invalidate a perfectly good payload whenever an unrelated edge row moved.
        for array in (junctions.offsets, junctions.boundary_right, junctions.strand):
            contiguous = np.ascontiguousarray(array)
            h.update(str(contiguous.dtype).encode())
            h.update(contiguous.tobytes())
        return h.hexdigest()

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
        collapse_duplicate_transcripts: bool = False,
        nrna_tolerance: int = NRNA_MERGE_TOLERANCE,
        alignable_zarr_path: str | Path | None = None,
        splice_blacklist_min_count: int = 2,
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
        collapse_duplicate_transcripts : bool
            When ``False`` (default), transcripts with identical exon
            coordinates raise ``ValueError`` (they are unidentifiable in
            quantification). When ``True``, collapse each such group to its
            lexicographically-smallest transcript ID instead of failing.
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
        alignable_zarr_path : path, optional
            Path to an alignable Zarr store built for the same genome+
            aligner.  When provided, the splice-junction artifact
            blacklist is derived from
            ``AlignableStore.splice_blacklist()`` and persisted as
            ``splice_blacklist.feather`` in the index.  When ``None``,
            no blacklist is written.  (Per-region mappability was
            removed in v0.5.0; SRD calibration does not consume it.)
        splice_blacklist_min_count : int
            Minimum per-row count for a (chrom, intron, read_length)
            artifact to enter the blacklist.  Default ``2``.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        feather_kwargs = {"compression": feather_compression}

        # -- Reference lengths ------------------------------------------------
        logger.info(f"[START] Loading references from {fasta_file}")
        ref_lengths = load_reference_lengths(fasta_file)
        logger.info(f"[DONE] Read {len(ref_lengths)} references")

        df = pd.DataFrame(list(ref_lengths.items()), columns=["ref", "length"])
        df.to_feather(output_dir / REF_LENGTHS_FEATHER, **feather_kwargs)
        if write_tsv:
            df.to_csv(output_dir / REF_LENGTHS_TSV, sep="\t", index=False)

        # -- Transcripts ------------------------------------------------------
        logger.info(f"[START] Reading transcripts from {gtf_file}")
        transcripts = read_transcripts(
            gtf_file,
            gtf_parse_mode=gtf_parse_mode,
            collapse_duplicate_transcripts=collapse_duplicate_transcripts,
        )
        logger.info(f"[DONE] Read {len(transcripts)} transcripts")

        # -- Synthetic nRNA transcripts ---------------------------------------
        logger.info("[START] Creating synthetic nRNA transcripts")
        synthetics, t_to_span_key, span_to_syn_idx, covered_equiv = create_nrna_transcripts(
            transcripts, tolerance=nrna_tolerance
        )
        # Synthetic transcripts are gene-neutral (see create_nrna_transcripts
        # docstring). Each synthetic gets its OWN dedicated g_index — we do
        # not map synthetics into any annotated gene, because a merged nRNA
        # span may be contributed to by transcripts from multiple genes and
        # rigel is transcript-centric.  The dedicated g_index keeps the
        # invariant ``t.g_index == row in g_df`` while preventing phantom
        # contributions from polluting real gene-level summaries.
        next_g_index = max((t.g_index for t in transcripts), default=-1) + 1

        # Assign t_index and g_index to each synthetic; build span→t_index
        span_to_nrna_t_index: dict[tuple, int] = {}
        next_t_index = len(transcripts)
        for syn in synthetics:
            syn.t_index = next_t_index
            next_t_index += 1
            syn.g_index = next_g_index
            next_g_index += 1
            span_key = (syn.ref, int(syn.strand), syn.start, syn.end)
            span_to_nrna_t_index[span_key] = syn.t_index

        # For covered spans, the nRNA entity is the nascent-equiv transcript
        for span_key, equiv_tx in covered_equiv.items():
            span_to_nrna_t_index[span_key] = equiv_tx.t_index
            equiv_tx.nrna_t_index = equiv_tx.t_index

        # Set nrna_t_index on every multi-exon annotated transcript
        for me_tidx, span_key in t_to_span_key.items():
            transcripts[me_tidx].nrna_t_index = span_to_nrna_t_index.get(span_key, -1)

        transcripts.extend(synthetics)
        logger.info(
            f"[DONE] {len(synthetics)} synthetics added, {len(transcripts)} total transcripts"
        )

        # Partition: annotated-only list for region table and splice junctions
        annotated_transcripts = [t for t in transcripts if not t.is_synthetic]

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

        # -- Genomic intervals + the v8 splice graph --------------------------
        logger.info("[START] Building genomic intervals + the splice graph")
        iv_df, nodes_df, edges_df = build_index_artifacts(transcripts, ref_lengths)
        logger.info(
            f"[DONE] {len(iv_df)} genomic intervals, {len(nodes_df)} nodes, {len(edges_df)} edges"
        )

        iv_df.to_feather(output_dir / INTERVALS_FEATHER, **feather_kwargs)
        if write_tsv:
            iv_df.to_csv(output_dir / INTERVALS_TSV, sep="\t", index=False)

        nodes_df.to_feather(output_dir / NODES_FEATHER, **feather_kwargs)
        edges_df.to_feather(output_dir / EDGES_FEATHER, **feather_kwargs)
        if write_tsv:
            nodes_df.to_csv(output_dir / NODES_TSV, sep="\t", index=False)
            edges_df.to_csv(output_dir / EDGES_TSV, sep="\t", index=False)

        # -- Splice-junction artifact blacklist (from alignable Zarr) -------
        if alignable_zarr_path is not None:
            from .splice_blacklist import load_splice_blacklist_from_zarr

            bl_df = load_splice_blacklist_from_zarr(
                alignable_zarr_path,
                min_count=splice_blacklist_min_count,
            )
            bl_df.to_feather(output_dir / SJ_BLACKLIST_FEATHER, **feather_kwargs)
            if write_tsv:
                bl_df.to_csv(output_dir / SJ_BLACKLIST_TSV, sep="\t", index=False)

        # -- Manifest --------------------------------------------------------
        # ⛔ Everything needed to REBUILD this index, because the previous manifest recorded neither the
        # sources nor the flags, so a rebuild had to infer both (the index-rebuild
        # entry). Defaults are written out explicitly: a rebuilder must not have to know this rigel
        # version's defaults to reproduce the artifact.
        logger.info("[START] Digesting build sources for the manifest")
        manifest = {
            "format_version": INDEX_FORMAT_VERSION,
            "rigel_version": _rigel_version(),
            "sources": {
                "fasta": source_record(fasta_file),
                "gtf": source_record(gtf_file),
                # `null`, never omitted: absent says "no store was supplied", a missing key would say
                # "this manifest predates the field".
                "alignable_zarr": (
                    None if alignable_zarr_path is None else source_record(alignable_zarr_path)
                ),
            },
            # ⚠ Hand-listed, and `tests/test_index_provenance.py` reads the expected key set off
            # `inspect.signature(build)` — so a NEW build parameter that does not reach here fails that
            # test rather than silently escaping the provenance.
            "build_flags": {
                "collapse_duplicate_transcripts": collapse_duplicate_transcripts,
                "feather_compression": feather_compression,
                "gtf_parse_mode": gtf_parse_mode,
                "nrna_tolerance": nrna_tolerance,
                "splice_blacklist_min_count": splice_blacklist_min_count,
                "write_tsv": write_tsv,
            },
        }
        logger.info("[DONE] Sources digested")
        with open(output_dir / MANIFEST_JSON, "w") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)

        logger.info(f"Index written to {output_dir}")

    # -- load (classmethod) ---------------------------------------------------

    @staticmethod
    def _build_gene_table(t_df: pd.DataFrame) -> pd.DataFrame:
        """Derive a gene-level summary table from the transcript table.

        Every ``g_index`` in ``t_df`` gets one row here, including synthetic
        nRNA gene rows (which are gene-neutral placeholders, one per synthetic
        transcript).  The ``is_synthetic`` column flags those rows so callers
        can filter them out of user-facing gene summaries.  Keeping them in
        the table preserves the invariant ``g_df.index == g_df["g_index"]``
        used by the EM aggregation kernels.
        """
        g_df = (
            t_df.groupby("g_index", sort=True)
            .agg(
                ref=("ref", "first"),
                start=("start", "min"),
                end=("end", "max"),
                strand=("strand", "first"),
                g_id=("g_id", "first"),
                g_name=("g_name", "first"),
                num_transcripts=("t_index", "count"),
                is_synthetic=("is_synthetic", "all"),
            )
            .reset_index()
        )
        # Re-apply categorical encoding (groupby first() returns plain strings).
        for col in ("ref", "g_id", "g_name"):
            if col in g_df.columns:
                g_df[col] = g_df[col].astype("category")
        return g_df

    @classmethod
    def load(
        cls,
        index_dir: str | Path,
        *,
        retain_test_structures: bool = False,
    ) -> "TranscriptIndex":
        """Load an index from the Feather files in *index_dir*.

        Builds cgranges interval trees and splice-junction lookup maps
        for fast overlap queries during quantification.

        Parameters
        ----------
        retain_test_structures : bool
            If True, keep Python-side structures (``_iv_t_set``,
            ``sj_map``) that are normally freed after C++ projection.
            Intended for unit tests that inspect these structures directly.
        """
        index_dir = str(index_dir)
        self = cls()
        self.index_dir = index_dir
        self._retain_test_structures = retain_test_structures

        # -- manifest (version gate) ------------------------------------------
        manifest = load_manifest(index_dir)
        if manifest is None:
            raise RuntimeError(
                f"Index at {index_dir} has no manifest.json. "
                f"Rebuild the index (rigel index --fasta ... --gtf ...). "
                f"Older indexes are not supported."
            )
        on_disk_version = int(manifest.get("format_version", -1))
        if on_disk_version != INDEX_FORMAT_VERSION:
            built_by = manifest.get("rigel_version")
            built_str = f" (built by rigel {built_by})" if built_by else ""
            raise RuntimeError(
                f"Index at {index_dir} has format_version={on_disk_version}"
                f"{built_str}, but this rigel build (v{_rigel_version()}) "
                f"requires {INDEX_FORMAT_VERSION}. "
                f"Rebuild: rigel index --fasta ... --gtf ... -o {index_dir}"
            )

        # -- reference lengths -------------------------------------------------
        ref_lengths_path = os.path.join(index_dir, REF_LENGTHS_FEATHER)
        if not os.path.exists(ref_lengths_path):
            raise RuntimeError(
                f"Index at {index_dir} is missing {REF_LENGTHS_FEATHER}. Rebuild the index."
            )
        self.ref_lengths = load_ref_lengths(ref_lengths_path)
        self.ref_names = list(self.ref_lengths.keys())
        self.ref_name_to_id = {name: i for i, name in enumerate(self.ref_names)}

        # -- transcripts ------------------------------------------------------
        self.t_df = pd.read_feather(os.path.join(index_dir, TRANSCRIPTS_FEATHER))
        if "t_index" not in self.t_df.columns:
            raise ValueError(
                f"Invalid index in {index_dir}: missing 't_index' column in {TRANSCRIPTS_FEATHER}"
            )
        if not (self.t_df.index == self.t_df["t_index"]).all():
            raise ValueError(
                f"Invalid index in {index_dir}: row index does not match "
                f"'t_index' column in {TRANSCRIPTS_FEATHER}; rebuild index"
            )

        # -- compact in-memory representation ---------------------------------
        # String columns → categorical (integer codes + small dictionary).
        for col in ("ref", "t_id", "g_id", "g_name", "g_type"):
            if col in self.t_df.columns:
                self.t_df[col] = self.t_df[col].astype("category")
        # Integer columns → narrowest safe dtype.
        _INT_DOWNCAST = {
            "strand": np.int8,
            "t_index": np.int32,
            "g_index": np.int32,
            "n_exons": np.int16,
            "nrna_t_index": np.int32,
            "nrna_n_contributors": np.int16,
            "length": np.int32,
        }
        for col, dtype in _INT_DOWNCAST.items():
            if col in self.t_df.columns:
                self.t_df[col] = self.t_df[col].astype(dtype)

        # -- gene table (derived) ---------------------------------------------
        self.g_df = cls._build_gene_table(self.t_df)
        if "g_index" not in self.g_df.columns:
            raise ValueError(f"Invalid derived gene table in {index_dir}: missing 'g_index' column")
        if not (self.g_df.index == self.g_df["g_index"]).all():
            raise ValueError(
                f"Invalid index in {index_dir}: derived gene table row index "
                f"does not match 'g_index'; rebuild index"
            )

        # -- fast numpy lookup arrays -----------------------------------------
        self.t_to_g_arr = self.t_df["g_index"].values
        self.t_to_strand_arr = self.t_df["strand"].values
        self.g_to_strand_arr = self.g_df["strand"].values

        # Per-transcript canonical reference id (matches index.ref_name_to_id
        # / BAM tid space, not pandas categorical codes).  Used by the
        # regional-exposure per-unit weight applier and any other code path
        # that needs ref ids without re-mapping categorical codes.
        _ref_cat = self.t_df["ref"].cat
        _cat_to_canonical_ref = np.array(
            [self.ref_name_to_id[str(name)] for name in _ref_cat.categories],
            dtype=np.int32,
        )
        self.t_to_ref_arr = _cat_to_canonical_ref[
            _ref_cat.codes.values.astype(np.int64, copy=False)
        ]

        # -- the splice graph: THE calibration partition ----------------------
        # ⚠ The load-time validation is the GRAPH-INTERNAL half only — I1/I2/I5-I9/I12. I3b, I4, I11
        # and I13 need the transcripts (~3 s to reconstruct at human scale) and run at BUILD. So a
        # `signature` or `flags` column that has drifted from the annotation loads clean; what this
        # catches is a graph that is internally inconsistent or truncated, which is what has
        # historically gone wrong with a stale index.
        from .calibration.splice_graph import load_edges, load_nodes, validate_graph

        logger.debug("Reading splice graph")
        nodes_path = os.path.join(index_dir, NODES_FEATHER)
        edges_path = os.path.join(index_dir, EDGES_FEATHER)
        for name, path in ((NODES_FEATHER, nodes_path), (EDGES_FEATHER, edges_path)):
            if not os.path.exists(path):
                raise RuntimeError(
                    f"Index at {index_dir} is missing {name} (the v8 splice graph). "
                    f"Rebuild the index (rigel index --fasta ... --gtf ... -o {index_dir})."
                )
        self.nodes_df = load_nodes(nodes_path)
        self.edges_df = load_edges(edges_path)
        validate_graph(self.nodes_df, self.edges_df, self.ref_lengths)

        # -- interval index (unified cgranges) ---------------------------------
        logger.debug("Reading intervals")
        iv_df = pd.read_feather(os.path.join(index_dir, INTERVALS_FEATHER))

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

        # Synthetic nRNA transcripts are intentionally excluded from
        # cgranges (see _gen_transcript_intervals) but their single-exon
        # spans must still be recorded for (a) the per-transcript exon
        # CSR consumed by C++ FL/derive-on-demand code and (b) the
        # `_t_exon_intervals` accessor used by tests.  Reconstruct them
        # from t_df (synthetics are single-exon by construction:
        # exon == [start, end)).
        if "is_synthetic" in self.t_df.columns:
            syn_mask = self.t_df["is_synthetic"].to_numpy(dtype=bool)
            if syn_mask.any():
                syn_t_idx = np.where(syn_mask)[0].astype(np.int32)
                syn_starts = self.t_df["start"].to_numpy(dtype=np.int32)[syn_mask]
                syn_ends = self.t_df["end"].to_numpy(dtype=np.int32)[syn_mask]
                for i, t_idx in enumerate(syn_t_idx):
                    t_exon_intervals[int(t_idx)] = np.array(
                        [[syn_starts[i], syn_ends[i]]], dtype=np.int32
                    )

        self._t_exon_intervals = t_exon_intervals
        logger.debug(f"Cached exon intervals for {len(t_exon_intervals)} transcripts")

        # -- splice junction indexes ------------------------------------------
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(os.path.join(index_dir, SJ_FEATHER))

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
            key = (_sj_rf_list[s], _sj_st_list[s], _sj_en_list[s], _sj_sd_list[s])
            sj_map[key] = frozenset(_sj_ti_list[s:e])
        self.sj_map = sj_map

        logger.debug(f"Splice junctions: {len(self.sj_map)} unique, {len(sj_df)} total")

        # -- C++ FragmentResolver (native fragment resolution) ------------------
        from .native import FragmentResolver

        ctx = FragmentResolver()

        # 0. Seed the ref-ID space in canonical ref_names order. Without this the
        #    resolver assigns ref ids by first-seen interval order (gene order),
        #    which scrambles them relative to index.ref_name_to_id and silently
        #    mis-routes the calibration accumulator deposit on multi-reference
        #    genomes (the accumulator partition is laid out in ref_names order;
        #    the deposit indexes it with the resolver's ref id). A single-ref
        #    synthetic index hides this — the two orderings coincide there.
        ctx.set_ref_names([str(name) for name in self.ref_names])

        # 1. Overlap index from collapsed data
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
        # C++ resolver owns the overlap data now; free the Python frozensets.
        if not retain_test_structures:
            self._iv_t_set = None

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
        # C++ resolver owns the SJ data now; free the Python dict.
        if not retain_test_structures:
            self.sj_map = None

        # 2b. Splice-junction artifact blacklist (optional)
        blacklist_path = os.path.join(index_dir, SJ_BLACKLIST_FEATHER)
        if os.path.exists(blacklist_path):
            logger.debug("Loading splice-artifact blacklist")
            bl_df = pd.read_feather(blacklist_path)
            ctx.build_sj_blacklist_map(
                bl_df["ref"].astype(str).tolist(),
                bl_df["start"].astype(np.int32).tolist(),
                bl_df["end"].astype(np.int32).tolist(),
                bl_df["max_anchor_left"].astype(np.int32).tolist(),
                bl_df["max_anchor_right"].astype(np.int32).tolist(),
            )
            self.sj_blacklist_size = int(len(bl_df))
            logger.info(f"Splice artifact blacklist: {len(bl_df):,} junctions active")
        else:
            # No blacklist file → artifact detection is off for this index.
            self.sj_blacklist_size = 0

        # 3. Per-transcript exon CSR for transcript-space FL computation.
        #    build_exon_csr() is the single owner of this CSR (it is also
        #    consumed by scoring), so reuse it here rather than duplicating
        #    the flattening loop — the two paths cannot drift.  NOTE: this
        #    also triggers build_exon_csr()'s free-after-build of
        #    ``_t_exon_intervals`` (unless retain_test_structures); the dict
        #    is no longer needed after this point (get_exon_intervals()
        #    falls back to the cached CSR arrays).
        exon_offsets, exon_starts_flat, exon_ends_flat, exon_cumsum_flat = self.build_exon_csr()
        # Per-transcript spliced length = sum of exon lengths, derived from
        # the same CSR: cumulative segment lengths differenced at the CSR
        # offsets (int64 accumulation, empty groups → 0).
        seg_len = (exon_ends_flat - exon_starts_flat).astype(np.int64)
        seg_cumsum = np.concatenate(([0], np.cumsum(seg_len)))
        t_lengths = seg_cumsum[exon_offsets[1:]] - seg_cumsum[exon_offsets[:-1]]
        ctx.build_exon_index(
            exon_offsets.tolist(),
            exon_starts_flat.tolist(),
            exon_ends_flat.tolist(),
            exon_cumsum_flat.tolist(),
            t_lengths.tolist(),
        )

        # 5. Metadata  (t_to_g_arr is already int32 from downcast above)
        ctx.set_metadata(
            self.t_to_g_arr.tolist(),
            len(self.t_to_g_arr),
        )

        # 6. nRNA mask + parent-index for derive-on-demand resolution.
        #    Synthetic nRNAs are intentionally absent from cgranges
        #    (see _gen_transcript_intervals).  The resolver materialises
        #    each synthetic candidate from real-tx hits via nrna_parent_.
        if "is_synthetic" in self.t_df.columns:
            is_synth = self.t_df["is_synthetic"].to_numpy(dtype=bool)
            ctx.set_nrna_status(is_synth.astype(np.uint8).tolist())
            if "nrna_t_index" in self.t_df.columns:
                nrna_idx = self.t_df["nrna_t_index"].to_numpy(dtype=np.int32)
                parent = np.full(nrna_idx.shape, -1, dtype=np.int32)
                valid = (nrna_idx >= 0) & (nrna_idx < is_synth.size)
                parent_is_synth = np.zeros_like(nrna_idx, dtype=bool)
                parent_is_synth[valid] = is_synth[nrna_idx[valid]]
                parent[parent_is_synth] = nrna_idx[parent_is_synth]
                ctx.set_nrna_parent_index(parent.tolist())

        self.resolver = ctx
        logger.debug("Built native FragmentResolver for C++ resolution")

        return self

    # -- query methods --------------------------------------------------------

    def query(
        self,
        exon: GenomicInterval,
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
        if self._iv_t_set is None:
            raise RuntimeError(
                "query() unavailable: interval sets were freed after C++ "
                "resolver construction. Use the C++ resolver for production "
                "overlap queries."
            )
        hits: list[tuple[int, int, int, frozenset[int]]] = []
        for h_start, h_end, label in self.cr.overlap(exon.ref, exon.start, exon.end):
            hits.append(
                (
                    h_start,
                    h_end,
                    self._iv_type[label],
                    self._iv_t_set[label],
                )
            )
        return hits

    def get_exon_intervals(self, t_idx: int) -> np.ndarray | None:
        """Return sorted exon ``[start, end)`` intervals for a transcript.

        Order-independent: once :meth:`build_exon_csr` has been called it
        frees the source ``_t_exon_intervals`` dict (memory optimization,
        unless ``retain_test_structures``), after which this accessor
        transparently reconstructs the interval from the cached CSR arrays.
        It therefore never returns a misleading ``None`` merely because
        scoring/loading ran first.

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
        if self._t_exon_intervals is not None:
            return self._t_exon_intervals.get(t_idx)
        # Dict was freed after build_exon_csr(); rebuild from the cached CSR.
        cache = getattr(self, "_exon_csr_cache", None)
        if cache is None:
            return None
        offsets, starts, ends, _ = cache
        if t_idx < 0 or t_idx + 1 >= len(offsets):
            return None
        lo = int(offsets[t_idx])
        hi = int(offsets[t_idx + 1])
        if lo == hi:
            return None
        return np.column_stack((starts[lo:hi], ends[lo:hi])).astype(np.int32)

    def build_exon_csr(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Build CSR arrays for per-transcript exon positions.

        Converts the per-transcript exon dict into flat CSR arrays
        suitable for direct C++ consumption, eliminating the need for
        a Python dict → C++ dict-unpacking round-trip.

        The result is cached: subsequent calls return the same arrays
        and the source ``_t_exon_intervals`` dict is freed after the
        first call to reclaim memory (unless ``retain_test_structures``).
        This free-after-build is why :meth:`get_exon_intervals` falls back
        to reconstructing intervals from the cached CSR — callers must not
        rely on ``_t_exon_intervals`` still being populated afterwards.

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
        # Return cached result if available.
        if hasattr(self, "_exon_csr_cache") and self._exon_csr_cache is not None:
            return self._exon_csr_cache

        n_t = self.num_transcripts
        offsets = np.zeros(n_t + 1, dtype=np.int32)
        empty = np.empty(0, dtype=np.int32)

        if self._t_exon_intervals is None or len(self._t_exon_intervals) == 0:
            result = (offsets, empty, empty, empty)
            self._exon_csr_cache = result
            return result

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

        # CSR arrays now hold all exon data; free the per-transcript dict.
        if not self._retain_test_structures:
            self._t_exon_intervals = None
        result = (offsets, starts, ends, cumsum_before)
        self._exon_csr_cache = result

        return result
