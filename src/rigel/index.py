"""Build, load and query the rigel reference index.

The index is built from a genome FASTA (with a samtools .fai) and a GTF annotation. It writes Feather files, with
optional TSV mirrors, to an output directory:

    ref_lengths.feather        — reference names and lengths
    transcripts.feather        — one row per transcript with integer indices
    intervals.feather          — each transcript's exons and span, and the intergenic gaps
    regions.feather            — the splice graph's regions (the calibration partition)
    edges.feather              — the splice graph's boundaries: contiguous and splice-junction
    sj.feather                 — annotated splice junctions from transcript introns
    splice_blacklist.feather   — (optional) splice-artifact sj derived
                                 from the alignable Zarr store
    unobservable.feather       — the unobservable mask: with the alignable store, its zero-
                                 mappability runs at its longest read length; without, the
                                 FASTA's N runs (`rigel.observable`)
    transcript_aliases.feather — dropped → kept transcript IDs of the duplicates
                                 ``--collapse-duplicate-transcripts`` folded
    manifest.json              — the format version, the build sources and flags

:class:`TranscriptIndex` builds the index (``build()``) and loads and queries it for quantification (``load()``).
"""

import collections
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Literal

import numpy as np
import pandas as pd
import pysam

from .config import CONSTANTS, IndexConfig
from .types import AnnotatedInterval, Interval, IntervalType, Strand
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

REGIONS_FEATHER = "regions.feather"
REGIONS_TSV = "regions.tsv"

BOUNDARIES_FEATHER = "edges.feather"
BOUNDARIES_TSV = "boundaries.tsv"

SJ_FEATHER = "sj.feather"
SJ_TSV = "sj.tsv"

SJ_BLACKLIST_FEATHER = "splice_blacklist.feather"
SJ_BLACKLIST_TSV = "splice_blacklist.tsv"

#: The unobservable mask: the runs of the reference no read pair is observed at (`rigel.observable`).
UNOBSERVABLE_FEATHER = "unobservable.feather"
UNOBSERVABLE_TSV = "unobservable.tsv"

#: ``dropped_t_id → kept_t_id`` for every transcript ``--collapse-duplicate-transcripts`` folded into an
#: identical one; written for the user, never read.
TRANSCRIPT_ALIASES_FEATHER = "transcript_aliases.feather"
TRANSCRIPT_ALIASES_TSV = "transcript_aliases.tsv"

MANIFEST_JSON = "manifest.json"

#: On-disk index format version, bumped whenever the schema or the meaning of any persisted column changes;
#: :meth:`TranscriptIndex.load` refuses any other.
INDEX_FORMAT_VERSION = 11


def _rigel_version() -> str:
    from . import __version__

    return str(__version__)


def _sha256_of_file(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(CONSTANTS.resources.digest_chunk_bytes):
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
    """Provenance for one build input: its path, its size and its content's sha256.

    The genome and the annotation are not in the index, so this hash is provenance, not a cache key on the index's
    own content (``partition_hash`` and ``graph_hash`` are computed on demand). It is sha256, unlike those blake2b
    keys, so a user can check it with ``shasum -a 256``.
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


def _blacklist_digest(columns) -> str:
    """A content hash of the splice blacklist exactly as the resolver receives it: the columns `load` passes to
    `build_sj_blacklist_map`, in its argument order and the rows' stored order. Gate:
    `tests/test_splice_blacklist.py`."""
    import hashlib

    refs, *positions = columns
    h = hashlib.blake2b(digest_size=8)
    h.update("\x00".join(refs).encode())
    for values in positions:
        h.update(np.asarray(values, dtype=np.int32).tobytes())
    return h.hexdigest()


def _check_blacklist_references(bl_df: pd.DataFrame, ref_lengths: dict[str, int], store) -> None:
    """A splice blacklist names its junctions by reference. One whose references the FASTA does not hold matches
    nothing, while the index would report splice-artifact detection on, so it refuses the build; references beyond the
    FASTA (decoys, alternates) are named in a warning."""
    if bl_df.empty:
        return
    refs = list(dict.fromkeys(bl_df["ref"].astype(str)))
    beyond = [r for r in refs if r not in ref_lengths]
    some = CONSTANTS.qc.examples
    if len(beyond) == len(refs):
        raise ValueError(
            f"None of the {len(refs):,} references in the splice blacklist of the alignable store {store} is in the "
            f"FASTA ({some(refs)}, against {some(list(ref_lengths))}): the store was built for another genome or "
            "names its references differently."
        )
    if beyond:
        logger.warning(
            "%s of the splice blacklist's %s references are not in the FASTA (%s); their junctions match nothing.",
            f"{len(beyond):,}",
            f"{len(refs):,}",
            some(beyond),
        )


class IndexLoadError(RuntimeError):
    """An index :meth:`TranscriptIndex.load` refuses: no manifest, another format version, a missing file, or a
    transcript table whose ``t_index`` is not its row order or whose ``g_index`` does not number the genes 0..n−1. The
    message names the rebuild."""


def load_manifest(index_dir: str | Path) -> dict | None:
    """Load ``manifest.json`` from an index directory; ``None`` when there is none, which
    :meth:`TranscriptIndex.load` refuses."""
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
# The annotation: references, transcripts and their nascent spans
# ---------------------------------------------------------------------------


def load_reference_lengths(fasta_file: str | Path) -> dict[str, int]:
    """Read reference names and lengths from a FASTA file's index (``samtools faidx``; pysam builds it when it is
    missing) rather than scanning the sequence.

    The lengths are read from the index's own text: pysam reports a length as a C int, clipped at 2^31 - 1. A reference
    longer than that is refused by name, since a coordinate is an int32 from the index to the scan, as it is in BAM.

    Returns
    -------
    dict[str, int]
        Ordered mapping from reference name to length.
    """
    fasta_file = str(fasta_file)
    with pysam.FastaFile(fasta_file) as fh:
        names = list(fh.references)
    with open(fasta_file + ".fai") as fai:
        indexed = dict(line.split("\t")[:2] for line in fai)
    ref_lengths = {name: int(indexed[name]) for name in names}
    longest = int(np.iinfo(np.int32).max)
    too_long = [name for name, length in ref_lengths.items() if length > longest]
    if too_long:
        k = CONSTANTS.qc.report_examples
        listed = "; ".join(f"{name} ({ref_lengths[name]:,} bp)" for name in too_long[:k]) + (
            f"; and {len(too_long) - k:,} more" if len(too_long) > k else ""
        )
        raise ValueError(
            f"{len(too_long):,} reference(s) are longer than {longest:,} bp, the most an int32 coordinate "
            f"addresses (as in BAM): {listed}"
        )
    return ref_lengths


def _duplicate_transcript_groups(transcripts: list[Transcript]) -> dict[tuple, list[int]]:
    """Group transcripts by identical exon structure.

    Two transcripts are duplicates when they share reference, strand and sorted ``(start, end)`` exon tuple; no read
    data can tell them apart. Transcripts that share an intron chain but differ at a UTR end are not duplicates.

    Returns ``{(ref, strand, exon_tuple): [transcript list-index, ...]}`` for the groups of two or more. Exon-less
    transcripts are ignored.
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
    """Raise ``ValueError`` if two transcripts share an identical exon structure (:func:`_duplicate_transcript_groups`).

    They must be resolved before indexing: in the GTF upstream, or with ``--collapse-duplicate-transcripts``.
    """
    duplicates = _duplicate_transcript_groups(transcripts)
    if not duplicates:
        return

    examples = []
    for (ref, strand, _exons), idxs in list(duplicates.items())[: CONSTANTS.qc.report_examples]:
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
        f"First {min(CONSTANTS.qc.report_examples, n_groups)} duplicate group(s):\n"
        + "\n".join(examples)
    )


def _collapse_duplicate_transcripts(transcripts: list[Transcript]) -> list[Transcript]:
    """Drop duplicate transcripts (:func:`_duplicate_transcript_groups`), keeping the lexicographically smallest ID per
    group and recording the dropped IDs as its aliases.

    Identical transcripts are unidentifiable, so the collapse is loss-free. Returns the filtered list in its original
    order (the input itself when there are no duplicates) and logs what was dropped.
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
        transcripts[keep].aliases = tuple(sorted(transcripts[i].t_id for i in dropped))
        if len(examples) < CONSTANTS.qc.report_examples:
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
        min(CONSTANTS.qc.report_examples, len(duplicates)),
        "\n".join(examples),
    )
    return [t for i, t in enumerate(transcripts) if i not in drop_idx]


def _on_the_fasta(
    transcripts: list[Transcript], ref_lengths: dict[str, int], gtf_parse_mode: str
) -> list[Transcript]:
    """The transcripts that lie on the FASTA. One on a reference it lacks, or past a reference's end, is refused by
    name under strict parsing and dropped with a warning under warn-and-skip."""
    off = [t for t in transcripts if t.ref not in ref_lengths or t.end > ref_lengths[t.ref]]
    if not off:
        return transcripts

    def where(t):
        if t.ref not in ref_lengths:
            return f"{t.t_id} on {t.ref}, which the FASTA lacks"
        return f"{t.t_id} ending at {t.ref}:{t.end:,}, past its end at {ref_lengths[t.ref]:,}"

    k = CONSTANTS.qc.report_examples
    names = "; ".join(where(t) for t in off[:k]) + (
        f"; and {len(off) - k:,} more" if len(off) > k else ""
    )
    if gtf_parse_mode == "strict":
        raise ValueError(f"{len(off):,} transcript(s) lie off the FASTA: {names}")
    logger.warning("Skipping %s transcript(s) that lie off the FASTA: %s", f"{len(off):,}", names)
    dropped = {id(t) for t in off}
    return [t for t in transcripts if id(t) not in dropped]


def _one_reference_per_gene(
    transcripts: list[Transcript], ref_lengths: dict[str, int] | None, gtf_parse_mode: str
) -> list[Transcript]:
    """The transcripts of each gene's one reference. A gene's row takes one reference and strand and spans its
    transcripts' extremes, so a gene on two references would span coordinates of both, and the annotation does not
    say which copy is the gene's own. Strict parsing refuses such a gene by name. Warn-and-skip keeps the copy on the
    longest reference, the earlier in the FASTA among equals, and drops the others with a warning: an alternate
    contig, a patch or an unplaced scaffold is a short stand-in for part of a chromosome. The choice is read off the
    FASTA's lengths, so a reader given none refuses in either mode."""
    refs: dict[str, list[str]] = {}
    for t in transcripts:
        on = refs.setdefault(t.g_id, [])
        if t.ref not in on:
            on.append(t.ref)
    split = {g_id: on for g_id, on in refs.items() if len(on) > 1}
    if not split:
        return transcripts
    k = CONSTANTS.qc.report_examples
    shown = list(split.items())[:k]
    more = f"; and {len(split) - k:,} more" if len(split) > k else ""
    if gtf_parse_mode == "strict" or ref_lengths is None:
        names = "; ".join(f"{g_id} on {', '.join(on)}" for g_id, on in shown) + more
        why = (
            ""
            if gtf_parse_mode == "strict"
            else ". Which copy to keep is read off the FASTA's reference lengths, and this reader was given none"
        )
        raise ValueError(
            f"{len(split):,} gene(s) have transcripts on more than one reference: {names}{why}"
        )
    fasta_order = {ref: i for i, ref in enumerate(ref_lengths)}
    kept = {
        g_id: min(on, key=lambda ref: (-ref_lengths[ref], fasta_order[ref]))
        for g_id, on in split.items()
    }
    on_kept = [t for t in transcripts if kept.get(t.g_id, t.ref) == t.ref]
    names = (
        "; ".join(
            f"{g_id} kept on {kept[g_id]}, dropped on {', '.join(ref for ref in on if ref != kept[g_id])}"
            for g_id, on in shown
        )
        + more
    )
    logger.warning(
        "%s gene(s) have transcripts on more than one reference. Each keeps its copy on the longest of them (the "
        "earlier in the FASTA among equals); skipping the %s transcript(s) on the others: %s",
        f"{len(split):,}",
        f"{len(transcripts) - len(on_kept):,}",
        names,
    )
    return on_kept


def read_transcripts(
    gtf_file: str | Path,
    *,
    gtf_parse_mode: Literal["strict", "warn-skip"] = "strict",
    collapse_duplicate_transcripts: bool = False,
    ref_lengths: dict[str, int] | None = None,
) -> list[Transcript]:
    """Parse a GTF file and return a sorted list of Transcript objects.

    Each transcript is assigned a sequential ``t_index`` and a ``g_index``
    (unique integer per ``g_id``). Transcripts are sorted by
    ``(ref, start, end, strand)`` so that downstream interval generation
    can process them in genomic order.

    Parameters
    ----------
    collapse_duplicate_transcripts : bool
        Transcripts with identical exon structure raise ``ValueError`` when ``False`` (the default), and are collapsed
        by :func:`_collapse_duplicate_transcripts` when ``True``.
    ref_lengths : dict, optional
        The FASTA's references and lengths, in the FASTA's order. A transcript on a reference it lacks, or running
        past one's end, follows ``gtf_parse_mode``: strict refuses it by name, warn-and-skip drops it with a warning.
        A gene with transcripts on more than one reference is refused by name under strict parsing; warn-and-skip
        keeps its copy on the longest reference (:func:`_one_reference_per_gene`).

    Raises
    ------
    ValueError
        On duplicate transcripts when ``collapse_duplicate_transcripts`` is ``False``; under strict parsing, on a
        transcript off the FASTA; and on a gene with transcripts on more than one reference, under strict parsing or
        when ``ref_lengths`` is not given.
    """
    transcripts = Transcript.read_gtf(
        str(gtf_file),
        parse_mode=gtf_parse_mode,
    )
    if ref_lengths is not None:
        transcripts = _on_the_fasta(transcripts, ref_lengths, gtf_parse_mode)
    transcripts = _one_reference_per_gene(transcripts, ref_lengths, gtf_parse_mode)

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


# -- Synthetic nRNA spans: start and end sites clustered within a tolerance ---


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


def _clustered_envelope(coords: np.ndarray, tolerance: int, outer: np.ufunc) -> np.ndarray:
    """Each coordinate's cluster envelope, in input order: the coordinates, sorted, cluster within *tolerance* bp
    (:func:`_cluster_coordinates`), and each takes its cluster's outermost value (``np.minimum`` for start sites,
    ``np.maximum`` for end sites)."""
    order = np.argsort(coords, kind="mergesort")
    ids = _cluster_coordinates(coords[order], tolerance)
    fill = np.iinfo(np.int64).max if outer is np.minimum else np.iinfo(np.int64).min
    envelope = np.full(ids[-1] + 1 if len(ids) else 0, fill, dtype=np.int64)
    outer.at(envelope, ids, coords[order])
    out = np.empty(len(coords), dtype=np.int64)
    out[order] = envelope[ids]
    return out


def _merged_span_keys(multi_exon: list[Transcript], tolerance: int) -> dict[int, tuple]:
    """Each multi-exon transcript's merged span ``(ref, strand, start, end)``, by ``t_index``: its start and end sites
    clustered within *tolerance* bp among the transcripts on its reference and strand, independently, the span being
    the clusters' outer envelope, so a synthetic is never smaller than a contributing transcript."""
    groups: dict[tuple[str, int], list[Transcript]] = collections.defaultdict(list)
    for t in multi_exon:
        groups[(t.ref, int(t.strand))].append(t)
    keys: dict[int, tuple] = {}
    for (ref, strand), grp in groups.items():
        starts = np.array([t.start for t in grp], dtype=np.int64)
        ends = np.array([t.end for t in grp], dtype=np.int64)
        spans = zip(
            _clustered_envelope(starts, tolerance, np.minimum).tolist(),
            _clustered_envelope(ends, tolerance, np.maximum).tolist(),
        )
        for t, (start, end) in zip(grp, spans):
            keys[t.t_index] = (ref, strand, start, end)
    return keys


def _covering_transcripts(transcripts: list[Transcript], spans) -> dict[tuple, Transcript]:
    """The annotated single-exon transcript that covers each span one covers (an exact match or full containment): the
    first, by start then end, on the span's reference and strand."""
    by_place: dict[tuple[str, int], list[tuple[int, int, Transcript]]] = collections.defaultdict(
        list
    )
    for t in transcripts:
        if len(t.exons) == 1:
            by_place[(t.ref, int(t.strand))].append((t.start, t.end, t))
    for candidates in by_place.values():
        candidates.sort(key=lambda x: (x[0], x[1]))
    covering: dict[tuple, Transcript] = {}
    for span_key in spans:
        ref, strand, m_start, m_end = span_key
        for s_start, s_end, s_tx in by_place.get((ref, strand), []):
            if s_start > m_start:
                break  # sorted by start, no more can contain
            if m_end <= s_end:
                covering[span_key] = s_tx
                break
    return covering


def _synthetic_span(span_key: tuple, n_contributors: int) -> Transcript:
    """The synthetic single-exon transcript over a merged span. It is gene-neutral: its sentinel ``g_id`` is its own
    ``t_id``, which keeps the ``g_id`` column non-null, and the index builder gives it a ``g_index`` of its own."""
    ref, strand, start, end = span_key
    t_id = f"RIGEL_NRNA_{ref}_{strand}_{start}_{end}"
    return Transcript(
        ref=ref,
        strand=Strand(strand),
        exons=[Interval(start, end)],
        length=end - start,
        t_id=t_id,
        g_id=t_id,
        g_name="",
        g_type="nascent_rna",
        is_single_exon=True,
        is_synthetic=True,
        span_n_contributors=n_contributors,
    )


def create_nrna_transcripts(
    transcripts: list[Transcript],
    tolerance: int = IndexConfig.nascent_tolerance,
) -> tuple[list[Transcript], dict[int, tuple], dict[tuple, Transcript]]:
    """Create synthetic nRNA transcripts and detect their annotated equivalents.

    1. Each multi-exon transcript's start and end sites are clustered within *tolerance* bp, independently, into a
       merged span (:func:`_merged_span_keys`).
    2. A single-exon transcript that already covers a merged span (exact match or full containment) is the span's
       entity, and its ``span_n_contributors`` counts the span's transcripts (:func:`_covering_transcripts`).
    3. Each uncovered span becomes a synthetic single-exon ``Transcript`` (``is_synthetic``, ``is_single_exon``).
       Synthetics are gene-neutral (:func:`_synthetic_span`): a merged span may draw on transcripts of several genes,
       so inheriting one contributor's gene would be arbitrary; gene summaries are derived downstream from annotated
       contributors only.

    *transcripts* is the sorted transcript list, ``t_index`` already assigned. Returns ``(synthetics, t_to_span_key,
    covered_equiv)``: the synthetics to append, in their spans' first-seen order (the caller assigns their ``t_index``
    and ``g_index``); each annotated multi-exon transcript's span key ``(ref, strand, start, end)``, by ``t_index``,
    from which the caller sets ``span_t_index``; and the annotated single-exon transcript covering each covered span.
    """
    multi_exon = [t for t in transcripts if len(t.exons) >= 2]
    if not multi_exon:
        return [], {}, {}

    t_to_span_key = _merged_span_keys(multi_exon, tolerance)
    # each span's contributing transcripts, in the spans' first-seen order
    contributors: dict[tuple, int] = {}
    for t in multi_exon:
        key = t_to_span_key[t.t_index]
        contributors[key] = contributors.get(key, 0) + 1

    covered_equiv = _covering_transcripts(transcripts, contributors)
    for span_key, s_tx in covered_equiv.items():
        s_tx.span_n_contributors += contributors[span_key]
    synthetics = [
        _synthetic_span(span_key, n)
        for span_key, n in contributors.items()
        if span_key not in covered_equiv
    ]
    logger.info(
        f"nRNA: {len(contributors)} merged spans "
        f"({len(covered_equiv)} annotated equiv, "
        f"{len(synthetics)} synthetic)"
    )
    return synthetics, t_to_span_key, covered_equiv


def _add_nascent_spans(transcripts: list[Transcript], tolerance: int) -> None:
    """Append the synthetic spans that model the transcripts' nascent RNA (:func:`create_nrna_transcripts`) and set
    every span entity's ``span_t_index``. Each synthetic gets a ``g_index`` of its own (they are gene-neutral), which
    keeps ``t.g_index == row in g_df`` and keeps synthetics out of annotated genes' summaries."""
    logger.info("[START] Creating synthetic nRNA transcripts")
    synthetics, t_to_span_key, covered_equiv = create_nrna_transcripts(
        transcripts, tolerance=tolerance
    )
    next_g_index = max((t.g_index for t in transcripts), default=-1) + 1
    span_entity_t_index: dict[tuple, int] = {}
    for i, syn in enumerate(synthetics):
        syn.t_index = len(transcripts) + i
        syn.g_index = next_g_index + i
        span_entity_t_index[(syn.ref, int(syn.strand), syn.start, syn.end)] = syn.t_index
    # A covered span's entity is the annotated single-exon transcript that covers it
    for span_key, equiv_tx in covered_equiv.items():
        span_entity_t_index[span_key] = equiv_tx.t_index
        equiv_tx.span_t_index = equiv_tx.t_index
    for me_tidx, span_key in t_to_span_key.items():
        transcripts[me_tidx].span_t_index = span_entity_t_index.get(span_key, -1)
    transcripts.extend(synthetics)
    logger.info(f"[DONE] {len(synthetics)} synthetics added, {len(transcripts)} total transcripts")


def build_splice_sj(transcripts: list[Transcript]) -> pd.DataFrame:
    """Extract splice junctions (introns) from all transcripts.

    Each intron of each transcript produces one row.
    The resulting DataFrame is sorted by (ref, start, end, strand).
    """
    rows = [
        (t.ref, start, end, t.strand, t.t_index) for t in transcripts for start, end in t.introns()
    ]
    rows.sort(key=lambda sj: (sj[0], sj[1], sj[2], sj[3]))
    return pd.DataFrame(rows, columns=["ref", "start", "end", "strand", "t_index"])


def write_bed12(
    transcripts: list[Transcript],
    bed_path: str | Path,
) -> Path:
    """Write transcripts as BED12 for minimap2 ``-j`` annotation.

    Each transcript with exons becomes one BED12 line, one block per exon.
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


def _gen_transcript_intervals(t: Transcript) -> Iterator[AnnotatedInterval]:
    """Yield exon and transcript-span intervals for a single transcript.

    Each exon produces an EXON interval.  One TRANSCRIPT interval spans
    the full transcript ``[first_exon.start, last_exon.end)``.

    Synthetic nRNA transcripts yield nothing: the C++ resolver derives
    their candidates from the transcripts a fragment hits, through each
    one's ``span_t_index``.  An annotated single-exon transcript that
    covers a span (``is_synthetic=False``) yields intervals like any
    other transcript.
    """
    if t.is_synthetic:
        return
    # Exons
    for e in t.exons:
        yield AnnotatedInterval(t.ref, e.start, e.end, t.strand, IntervalType.EXON, t.t_index)
    # One transcript span
    if t.exons:
        yield AnnotatedInterval(
            t.ref, t.exons[0].start, t.exons[-1].end, t.strand, IntervalType.TRANSCRIPT, t.t_index
        )


@dataclass(frozen=True, slots=True)
class IntergenicSpan:
    """A reference span containing no real (non-synthetic) transcript.

    Half-open coordinates ``[start, end)``.
    """

    start: int
    end: int


@dataclass(frozen=True, slots=True)
class GenicSpan:
    """A reference span covered by one connected (strand-agnostic) transcript cluster.

    ``transcripts`` is the cluster's transcript list in input order; it
    is strand-agnostic and excludes synthetic nRNAs (which are
    constructed downstream of the layout sweep).
    """

    start: int
    end: int
    transcripts: tuple = field(default_factory=tuple)


def iter_reference_layout(
    ref_length: int,
    transcripts: list[Transcript],
) -> Iterator["IntergenicSpan | GenicSpan"]:
    """Yield interleaved ``IntergenicSpan`` / ``GenicSpan`` tiling [0, ref_length).

    ``transcripts`` is the per-reference list (already filtered to one
    reference), and **must be sorted by ``(start, end)``**. Synthetic
    transcripts are excluded from the sweep so they cannot coalesce
    real genic spans.

    Invariants:
      - Spans alternate (no two intergenic in a row) and tile the
        reference exactly with no gaps and no overlaps.
      - Each ``GenicSpan`` contains at least one (non-synthetic)
        transcript with at least one exon.
      - If the reference contains no real transcripts, a single
        ``IntergenicSpan(0, ref_length)`` is yielded (when
        ``ref_length > 0``).
    """
    real_ts = [t for t in transcripts if not t.is_synthetic]

    if not real_ts:
        if ref_length > 0:
            yield IntergenicSpan(0, ref_length)
        return

    cursor = 0
    cluster_start = real_ts[0].start
    cluster_end = real_ts[0].end
    cluster: list[Transcript] = [real_ts[0]]

    for t in real_ts[1:]:
        if t.start > cluster_end:
            # Close current cluster.
            if cursor < cluster_start:
                yield IntergenicSpan(cursor, cluster_start)
            yield GenicSpan(cluster_start, cluster_end, tuple(cluster))
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
        yield IntergenicSpan(cursor, cluster_start)
    yield GenicSpan(cluster_start, cluster_end, tuple(cluster))
    if cluster_end < ref_length:
        yield IntergenicSpan(cluster_end, ref_length)


def _emit_genomic_intervals(
    ref: str,
    layout: Iterator["IntergenicSpan | GenicSpan"],
) -> Iterator[AnnotatedInterval]:
    """Convert a per-reference layout stream into ``AnnotatedInterval`` rows.

    Each ``IntergenicSpan`` becomes one INTERGENIC interval. Each
    ``GenicSpan`` is decomposed via :func:`_gen_transcript_intervals`
    (per-exon EXON intervals plus one TRANSCRIPT span per transcript).
    Synthetic transcripts are filtered out by ``_gen_transcript_intervals``.
    """
    for span in layout:
        if isinstance(span, IntergenicSpan):
            yield AnnotatedInterval(ref, span.start, span.end)
        else:
            assert isinstance(span, GenicSpan)
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
    """Build ``intervals.feather`` and the splice graph (``regions`` + ``boundaries``).

    The per-reference layout (intergenic / genic spans) feeds only the annotated interval table
    (:func:`_emit_genomic_intervals`); the calibration partition is built separately, by
    :func:`rigel.calibration.splice_graph.build_splice_graph`.

    Returns ``(intervals_df, regions_df, edges_df)``, typed per the on-disk schemas. ``intervals_df``
    carries the columns of :class:`AnnotatedInterval`, sorted by ``(ref, start, end, strand)``.

    Transcripts must be sorted by ``(ref, start, end)``.
    """
    from .calibration.splice_graph import build_splice_graph, validate_graph

    ref_transcripts = _group_transcripts_by_ref(transcripts, ref_lengths)

    intervals: list[AnnotatedInterval] = []

    for ref, ref_length in ref_lengths.items():
        layout = iter_reference_layout(ref_length, ref_transcripts.get(ref, []))
        intervals.extend(_emit_genomic_intervals(ref, layout))

    intervals.sort(key=lambda iv: (iv.ref, iv.start, iv.end, iv.strand))
    iv_df = pd.DataFrame(intervals, columns=AnnotatedInterval._fields)
    regions_df, edges_df = build_splice_graph(transcripts, ref_lengths)
    validate_graph(regions_df, edges_df, ref_lengths, transcripts=transcripts)
    return iv_df, regions_df, edges_df


# ---------------------------------------------------------------------------
# Build steps
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _TableWriter:
    """Writes each index table as Feather, with its TSV mirror when the build writes them."""

    output_dir: Path
    compression: str
    write_tsv: bool

    def __call__(self, df: pd.DataFrame, feather_name: str, tsv_name: str) -> None:
        df.to_feather(self.output_dir / feather_name, compression=self.compression)
        if self.write_tsv:
            df.to_csv(self.output_dir / tsv_name, sep="\t", index=False)


@dataclass(frozen=True, slots=True)
class _GenomeInputs:
    """What the build reads of the genome: the references, the splice-artifact blacklist (``None`` without an
    alignable store) and the unobservable mask with the manifest's record of what it is (`rigel.observable`)."""

    ref_lengths: dict[str, int]
    blacklist: pd.DataFrame | None
    mask: pd.DataFrame
    mask_record: dict

    @property
    def masked_bases(self) -> int:
        return int((self.mask["end"] - self.mask["start"]).sum())

    @property
    def genome_bases(self) -> int:
        return int(sum(self.ref_lengths.values()))


def _read_genome(
    fasta_file: str | Path,
    alignable_zarr_path: str | Path | None,
    splice_blacklist_min_count: int,
) -> _GenomeInputs:
    """The references, the blacklist and the mask, each read and checked before the first write: a store that cannot
    be read, is incomplete, or was built for other references or another assembly refuses the build and leaves the
    index it would replace intact."""
    logger.info(f"[START] Loading references from {fasta_file}")
    ref_lengths = load_reference_lengths(fasta_file)
    logger.info(f"[DONE] Read {len(ref_lengths)} references")

    blacklist = None
    if alignable_zarr_path is not None:
        from .splice_blacklist import load_splice_blacklist_from_zarr

        blacklist = load_splice_blacklist_from_zarr(
            alignable_zarr_path,
            min_count=splice_blacklist_min_count,
        )
        _check_blacklist_references(blacklist, ref_lengths, alignable_zarr_path)

    from .observable import mask_from_fasta, mask_from_store

    if alignable_zarr_path is None:
        logger.info("[START] Recording the reference's unknown runs (no alignable store)")
        mask = mask_from_fasta(fasta_file, ref_lengths)
        record = {"rule": "the reference's unknown bases (N or n)", "read_length": None}
    else:
        logger.info(f"[START] Recording the unobservable mask from {alignable_zarr_path}")
        mask, record = mask_from_store(alignable_zarr_path, ref_lengths, fasta_file)
    genome = _GenomeInputs(ref_lengths, blacklist, mask, record)
    logger.info(
        f"[DONE] {len(mask):,} unobservable runs: {genome.masked_bases:,} of {genome.genome_bases:,} bp"
    )
    return genome


def _write_annotation_tables(
    write: _TableWriter, transcripts: list[Transcript], ref_lengths: dict[str, int]
) -> None:
    """The transcript table and its alias map, the annotated splice junctions, the genomic intervals and the splice
    graph's regions and boundaries."""
    write(transcripts_to_dataframe(transcripts), TRANSCRIPTS_FEATHER, TRANSCRIPTS_TSV)
    aliases_df = pd.DataFrame(
        [(alias, t.t_id) for t in transcripts for alias in t.aliases],
        columns=["dropped_t_id", "kept_t_id"],
        dtype=object,
    )
    write(aliases_df, TRANSCRIPT_ALIASES_FEATHER, TRANSCRIPT_ALIASES_TSV)

    logger.info("[START] Building splice junctions")
    sj_df = build_splice_sj([t for t in transcripts if not t.is_synthetic])
    logger.info(f"[DONE] Found {len(sj_df)} splice junctions")
    write(sj_df, SJ_FEATHER, SJ_TSV)

    logger.info("[START] Building genomic intervals + the splice graph")
    iv_df, regions_df, edges_df = build_index_artifacts(transcripts, ref_lengths)
    logger.info(
        f"[DONE] {len(iv_df)} genomic intervals, {len(regions_df)} regions, {len(edges_df)} boundaries"
    )
    write(iv_df, INTERVALS_FEATHER, INTERVALS_TSV)
    write(regions_df, REGIONS_FEATHER, REGIONS_TSV)
    write(edges_df, BOUNDARIES_FEATHER, BOUNDARIES_TSV)


def _write_genome_tables(write: _TableWriter, genome: _GenomeInputs) -> None:
    """The splice blacklist, when there is a store (one left by an earlier build is removed first, so a rebuild
    without a store leaves none behind), and the unobservable mask."""
    (write.output_dir / SJ_BLACKLIST_FEATHER).unlink(missing_ok=True)
    (write.output_dir / SJ_BLACKLIST_TSV).unlink(missing_ok=True)
    if genome.blacklist is not None:
        write(genome.blacklist, SJ_BLACKLIST_FEATHER, SJ_BLACKLIST_TSV)
    write(genome.mask, UNOBSERVABLE_FEATHER, UNOBSERVABLE_TSV)


def _write_manifest(
    output_dir: Path,
    sources: dict[str, str | Path | None],
    build_flags: dict,
    genome: _GenomeInputs,
) -> None:
    """``manifest.json``: everything needed to rebuild this index — each source's path, size and digest (``null`` for
    an alignable store not supplied), every flag with its default written out, so a rebuild need not know this
    version's defaults, and what the mask is and how much of the reference it holds."""
    logger.info("[START] Digesting build sources for the manifest")
    manifest = {
        "format_version": INDEX_FORMAT_VERSION,
        "rigel_version": _rigel_version(),
        # `null`, never omitted: null says "no store was supplied", a missing key would say the field is unknown
        "sources": {
            name: None if path is None else source_record(path) for name, path in sources.items()
        },
        "build_flags": build_flags,
        # the runs are `unobservable.feather`
        "observable": {
            **genome.mask_record,
            "unobservable_bases": genome.masked_bases,
            "observable_bases": genome.genome_bases - genome.masked_bases,
        },
    }
    logger.info("[DONE] Sources digested")
    with open(output_dir / MANIFEST_JSON, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# Load steps
# ---------------------------------------------------------------------------

#: The files every index has, checked before any is read. The splice blacklist is read only when the manifest
#: records its store, and the alias table never.
_REQUIRED_FILES = (
    REF_LENGTHS_FEATHER,
    TRANSCRIPTS_FEATHER,
    REGIONS_FEATHER,
    BOUNDARIES_FEATHER,
    INTERVALS_FEATHER,
    SJ_FEATHER,
    UNOBSERVABLE_FEATHER,
)

#: The transcript table's integer columns, each held at its narrowest safe dtype.
_TRANSCRIPT_INT_DTYPES = {
    "strand": np.int8,
    "t_index": np.int32,
    "g_index": np.int32,
    "span_t_index": np.int32,
    "span_n_contributors": np.int16,
    "length": np.int32,
}


def _checked_manifest(index_dir: str) -> dict:
    """The index's manifest, refused when absent or of another format version, and every required file present."""
    manifest = load_manifest(index_dir)
    if manifest is None:
        raise IndexLoadError(
            f"Index at {index_dir} has no manifest.json. "
            f"Rebuild the index (rigel index --fasta ... --gtf ...). "
            f"Older indexes are not supported."
        )
    on_disk_version = int(manifest.get("format_version", -1))
    if on_disk_version != INDEX_FORMAT_VERSION:
        built_by = manifest.get("rigel_version")
        built_str = f" (built by rigel {built_by})" if built_by else ""
        raise IndexLoadError(
            f"Index at {index_dir} has format_version={on_disk_version}"
            f"{built_str}, but this rigel build (v{_rigel_version()}) "
            f"requires {INDEX_FORMAT_VERSION}. "
            f"Rebuild: rigel index --fasta ... --gtf ... -o {index_dir}"
        )
    missing = [
        name for name in _REQUIRED_FILES if not os.path.exists(os.path.join(index_dir, name))
    ]
    if missing:
        raise IndexLoadError(
            f"Index at {index_dir} is missing {', '.join(missing)}. "
            f"Rebuild the index (rigel index --fasta ... --gtf ... -o {index_dir})."
        )
    return manifest


def _load_transcript_table(index_dir: str) -> pd.DataFrame:
    """The transcript table, its rows in ``t_index`` order, held compactly: string columns categorical (integer codes
    and a small dictionary), integer columns at their narrowest safe dtype."""
    t_df = pd.read_feather(os.path.join(index_dir, TRANSCRIPTS_FEATHER))
    if not (t_df.index == t_df["t_index"]).all():
        raise IndexLoadError(
            f"Invalid index in {index_dir}: row index does not match "
            f"'t_index' column in {TRANSCRIPTS_FEATHER}; rebuild index"
        )
    for col in ("ref", "t_id", "g_id", "g_name", "g_type"):
        t_df[col] = t_df[col].astype("category")
    for col, dtype in _TRANSCRIPT_INT_DTYPES.items():
        t_df[col] = t_df[col].astype(dtype)
    return t_df


def _key_runs(*keys: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sort rows by *keys*, the first the most significant, and find the runs of rows that agree on every key:
    returns the sort order and each run's first and one-past-last position in it."""
    order = np.lexsort(keys[::-1])
    n = len(order)
    changed = np.empty(n, dtype=bool)
    if n > 0:
        changed[0] = True
    if n > 1:
        changed[1:] = np.logical_or.reduce([np.diff(key[order]) != 0 for key in keys])
    first = np.where(changed)[0]
    return order, first, np.append(first[1:], n)


def _collapsed_intervals(iv_df: pd.DataFrame) -> tuple[list, list, list, list, list, list]:
    """The resolver's overlap index, as the six lists ``FragmentResolver.build_overlap_index`` takes: one interval per
    unique ``(ref, start, end, interval_type)`` — refs, starts, ends and types — and, in CSR, the sorted transcripts
    of every row it merges, the intergenic sentinel -1 dropped."""
    refs = iv_df["ref"].values
    starts = iv_df["start"].values
    ends = iv_df["end"].values
    itypes = iv_df["interval_type"].values
    _unique_refs, ref_codes = np.unique(refs, return_inverse=True)
    order, first, last = _key_runs(ref_codes, starts, ends, itypes)
    t_sorted = iv_df["t_index"].values[order].tolist()
    t_flat: list[int] = []
    t_offsets: list[int] = [0]
    for lo, hi in zip(first.tolist(), last.tolist()):
        t_flat.extend(sorted({t for t in t_sorted[lo:hi] if t >= 0}))
        t_offsets.append(len(t_flat))
    logger.debug(f"Interval index: {len(iv_df)} rows → {len(first)} collapsed")
    return (
        refs[order][first].tolist(),
        starts[order][first].tolist(),
        ends[order][first].tolist(),
        itypes[order][first].tolist(),
        t_flat,
        t_offsets,
    )


def _exon_intervals(iv_df: pd.DataFrame, t_df: pd.DataFrame) -> dict[int, np.ndarray]:
    """Each transcript's exons, ``(n_exons, 2)`` int32 ``[start, end)`` sorted by start, by ``t_index``. Synthetics are
    absent from the interval table (see :func:`_gen_transcript_intervals`), but the exon CSR the native code reads
    still needs their single exon, ``[start, end)`` from the transcript table."""
    exon = (iv_df["interval_type"].values == int(IntervalType.EXON)) & (
        iv_df["t_index"].values >= 0
    )
    ex_tidx = iv_df["t_index"].values[exon]
    ex_start = iv_df["start"].values[exon]
    ex_end = iv_df["end"].values[exon]
    # Sorted by (t_index, start) so exons within each transcript are ordered, then split at the group boundaries
    order = np.lexsort((ex_start, ex_tidx))
    ex_tidx = ex_tidx[order]
    coords = np.column_stack((ex_start[order], ex_end[order])).astype(np.int32)
    unique_tidx, group_starts, group_counts = np.unique(
        ex_tidx, return_index=True, return_counts=True
    )
    t_exon_intervals: dict[int, np.ndarray] = {
        int(t): coords[s : s + n] for t, s, n in zip(unique_tidx, group_starts, group_counts)
    }
    syn_mask = t_df["is_synthetic"].to_numpy(dtype=bool)
    if syn_mask.any():
        syn_starts = t_df["start"].to_numpy(dtype=np.int32)[syn_mask]
        syn_ends = t_df["end"].to_numpy(dtype=np.int32)[syn_mask]
        for t_idx, start, end in zip(np.where(syn_mask)[0].tolist(), syn_starts, syn_ends):
            t_exon_intervals[t_idx] = np.array([[start, end]], dtype=np.int32)
    logger.debug(f"Cached exon intervals for {len(t_exon_intervals)} transcripts")
    return t_exon_intervals


def _sj_exact_match_map(sj_df: pd.DataFrame) -> dict[tuple, frozenset[int]]:
    """``(ref, start, end, strand) → frozenset(t_index)``: every annotated junction with the transcripts holding it."""
    refs = sj_df["ref"].values
    starts = sj_df["start"].values
    ends = sj_df["end"].values
    strands = sj_df["strand"].values
    _unique_refs, ref_codes = np.unique(refs, return_inverse=True)
    order, first, last = _key_runs(ref_codes, starts, ends, strands)
    keys = list(
        zip(
            refs[order].tolist(),
            starts[order].tolist(),
            ends[order].tolist(),
            strands[order].tolist(),
        )
    )
    t_sorted = sj_df["t_index"].values[order].tolist()
    sj_map = {keys[lo]: frozenset(t_sorted[lo:hi]) for lo, hi in zip(first.tolist(), last.tolist())}
    logger.debug(f"Splice junctions: {len(sj_map)} unique, {len(sj_df)} total")
    return sj_map


def _sj_map_csr(sj_map: dict[tuple, frozenset[int]]) -> tuple[list, list, list, list, list, list]:
    """The junction map as the six lists ``FragmentResolver.build_sj_map`` takes: refs, starts, ends, strands and,
    in CSR, each junction's sorted transcripts."""
    refs, starts, ends, strands, t_flat, t_offsets = [], [], [], [], [], [0]
    for (ref, start, end, strand), tset in sj_map.items():
        refs.append(ref)
        starts.append(start)
        ends.append(end)
        strands.append(strand)
        t_flat.extend(sorted(tset))
        t_offsets.append(len(t_flat))
    return refs, starts, ends, strands, t_flat, t_offsets


def _install_blacklist(resolver, index_dir: str, manifest: dict) -> tuple[int, str]:
    """Install the splice-artifact blacklist, used only when the manifest records the store it was built from;
    returns the active junctions' count (0 is detection off) and the digest of the rows handed to the resolver
    (``"none"`` when none are), which the scan cache's key reads, because the blacklist changes what a scan tallies.
    A blacklist with no store recorded, or a recorded store with no blacklist, leaves artifact detection off with no
    other sign, so the load warns."""
    blacklist_path = os.path.join(index_dir, SJ_BLACKLIST_FEATHER)
    store = manifest.get("sources", {}).get("alignable_zarr")
    has_blacklist = os.path.exists(blacklist_path)
    if store is not None and has_blacklist:
        logger.debug("Loading splice-artifact blacklist")
        bl_df = pd.read_feather(blacklist_path)
        columns = (
            bl_df["ref"].astype(str).tolist(),
            bl_df["start"].astype(np.int32).tolist(),
            bl_df["end"].astype(np.int32).tolist(),
            bl_df["max_anchor_left"].astype(np.int32).tolist(),
            bl_df["max_anchor_right"].astype(np.int32).tolist(),
        )
        resolver.build_sj_blacklist_map(*columns)
        logger.info(f"Splice artifact blacklist: {len(bl_df):,} sj active")
        return int(len(bl_df)), _blacklist_digest(columns)
    if has_blacklist:
        logger.warning(
            f"{SJ_BLACKLIST_FEATHER} is not applied: this index's manifest records no alignable "
            "store it was built from, so splice-artifact detection is off. Rebuild the index "
            "with --alignable-zarr to use a blacklist."
        )
    elif store is not None:
        logger.warning(
            f"This index was built from the alignable store {store['path']}, but its "
            f"{SJ_BLACKLIST_FEATHER} is missing, so splice-artifact detection is off. "
            "Rebuild the index."
        )
    return 0, "none"


def _nrna_parents(t_df: pd.DataFrame) -> np.ndarray:
    """Each transcript's synthetic span, by ``t_index``, -1 where its span entity is no synthetic. Synthetics are absent
    from the interval table (see :func:`_gen_transcript_intervals`), so the resolver derives each synthetic candidate
    from the annotated transcripts a fragment hits, through this index."""
    is_synth = t_df["is_synthetic"].to_numpy(dtype=bool)
    span_idx = t_df["span_t_index"].to_numpy(dtype=np.int32)
    parent = np.full(span_idx.shape, -1, dtype=np.int32)
    valid = (span_idx >= 0) & (span_idx < is_synth.size)
    parent_is_synth = np.zeros_like(span_idx, dtype=bool)
    parent_is_synth[valid] = is_synth[span_idx[valid]]
    parent[parent_is_synth] = span_idx[parent_is_synth]
    return parent


# ---------------------------------------------------------------------------
# TranscriptIndex — unified index class
# ---------------------------------------------------------------------------


class TranscriptIndex:
    """The in-memory reference index.

    ``TranscriptIndex.build()`` writes the on-disk index from a FASTA and a GTF; ``TranscriptIndex.load()`` reads it
    back for quantification. The gene table is derived from the transcript table at load, never stored.
    """

    def __init__(self):
        self.index_dir: str | None = None
        #: the splice graph's region partition; :meth:`load` requires it (older formats are refused).
        self.regions_df = None
        #: the unobservable mask, ``ref, start, end`` runs (`rigel.observable.check_mask`); :meth:`load` requires it.
        self.unobservable_df: pd.DataFrame | None = None
        self.edges_df = None
        self.t_df: pd.DataFrame | None = None
        self.g_df: pd.DataFrame | None = None
        self.t_to_g_arr: np.ndarray | None = None
        self.t_to_strand_arr: np.ndarray | None = None

        # Splice-junction exact-match map  (ref, start, end, strand) → frozenset[int]
        self.sj_map: dict | None = None

        # Splice-artifact blacklist size, set at load():
        #   None → not loaded
        #   0    → no blacklisted sj active (artifact detection is off)
        #   >0   → number of blacklisted sj active (detection is on)
        self.sj_blacklist_size: int | None = None
        # The digest of the blacklist rows handed to the resolver, "none" when none are: set at load(), and read
        # by the scan cache's key, because the blacklist changes what a scan tallies.
        self.sj_blacklist_digest: str | None = None

        # Per-transcript exon intervals for coverage-weight model.
        # Maps t_index → (n_exons, 2) int32 array of [start, end) intervals
        # sorted by genomic start position.
        self._t_exon_intervals: dict[int, np.ndarray] | None = None
        # The same exons as flat CSR arrays, built once by build_exon_csr().
        self._exon_csr_cache: tuple | None = None

        # When True, load() keeps the Python-side copies of what the native resolver owns
        # (sj_map, _t_exon_intervals).
        self._retain_lookup_maps: bool = False

    # -- properties -----------------------------------------------------------

    @property
    def num_transcripts(self) -> int:
        return len(self.t_df)

    @property
    def num_genes(self) -> int:
        """Rows in ``g_df``, annotated genes and synthetic nRNA placeholders, so ``g_index == row``.

        User-facing summaries use :attr:`num_annotated_genes`.
        """
        return len(self.g_df)

    @property
    def num_annotated_genes(self) -> int:
        """Number of real (non-synthetic) gene rows in ``g_df``."""
        return int((~self.g_df["is_synthetic"].to_numpy()).sum())

    @property
    def partition_hash(self) -> str:
        """16-hex-char content hash of the region partition a scan deposits into: the payload's region axis.

        Hashes the ``(boundary_positions, ref_pos_offsets, region_types)`` triple that
        :meth:`rigel.native.BamScanner.set_regions` receives, plus the reference lengths, so two indexes share a
        hash iff they give the accumulator payload the same region axis. Computed on demand, never stored, so it
        cannot go stale against the Feather files.

        It is not a scan's key. The scan reads more of the index than the region bounds: resolution reads which
        transcripts hold a fragment's blocks and junctions, and the junctions and the splice blacklist decide what
        it offers, so two indexes that share this hash can tally one BAM differently. ``edges.feather`` can change
        while every region stays byte-identical too; :attr:`graph_hash` adds the junction axis.
        """
        import hashlib

        from .calibration.splice_graph import build_region_partition_arrays

        h = hashlib.blake2b(digest_size=8)
        for arr in build_region_partition_arrays(self):
            a = np.ascontiguousarray(arr)
            h.update(str(a.dtype).encode())
            h.update(a.tobytes())
        for name, length in self.ref_lengths.items():
            h.update(f"{name}:{length}|".encode())
        return h.hexdigest()

    @property
    def graph_hash(self) -> str:
        """16-hex-char content hash of the accumulator payload's two axes: its regions and its junctions.

        A tally's sj axis is meaningless against a different sj CSR, and the boundaries can change while the regions
        stay byte-identical, so this hashes ``partition_hash``'s inputs plus the sj CSR that crosses into ``set_sj``:
        the donor offsets, the acceptor region-bound indices and the annotated strands. A payload records it, so a
        cached tally can be checked against the index it is read with. It is not all a scan depends on: the
        transcript models resolution reads, the splice blacklist and the scan settings are not in it. Computed on
        demand, never stored.
        """
        import hashlib

        from .calibration.splice_graph import build_sj_arrays

        h = hashlib.blake2b(digest_size=8)
        h.update(self.partition_hash.encode())
        sj = build_sj_arrays(self)
        # `edge_row` is left out: a join key back to `edges_df` that never crosses the ABI, so hashing it would
        # invalidate a good payload whenever an unrelated boundary row moved
        for array in (sj.offsets, sj.boundary_right, sj.strand):
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
        feather_compression: str = IndexConfig.feather_compression,
        write_tsv: bool = IndexConfig.write_tsv,
        gtf_parse_mode: Literal["strict", "warn-skip"] = IndexConfig.gtf_parse_mode,
        collapse_duplicate_transcripts: bool = IndexConfig.collapse_duplicate_transcripts,
        nascent_tolerance: int = IndexConfig.nascent_tolerance,
        alignable_zarr_path: str | Path | None = None,
        splice_blacklist_min_count: int = IndexConfig.splice_blacklist_min_count,
    ) -> None:
        """Build the rigel reference index and write to disk.

        The settings are validated as one :class:`~rigel.config.IndexConfig`, so a value it refuses raises
        ``ValueError`` before anything is written, whoever the caller.

        Parameters
        ----------
        fasta_file : path
            Genome FASTA (must be indexed with ``samtools faidx``).
        gtf_file : path
            Gene annotation in GTF format (GENCODE recommended). One that yields no transcripts
            raises ``ValueError`` before any file is written.
        output_dir : path
            Directory to write index files into (created if needed). An existing
            ``manifest.json`` in it is deleted before the first write.
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
        nascent_tolerance : int
            Max distance (bp) for clustering multi-exon transcripts' start
            and end sites into the synthetic spans that model their nascent RNA.
        alignable_zarr_path : path, optional
            Path to an alignable Zarr store built for the same genome+
            aligner.  When provided, the splice-junction artifact
            blacklist is derived by
            :func:`rigel.splice_blacklist.load_splice_blacklist_from_zarr` and persisted as
            ``splice_blacklist.feather`` in the index.  When ``None``,
            no blacklist is written and one left by an earlier build is removed.
            :meth:`load` applies the blacklist only when the manifest records the store.
        splice_blacklist_min_count : int
            Minimum per-row count for a (chrom, intron, read_length)
            artifact to enter the blacklist.  Default ``2``.
        """
        IndexConfig(
            feather_compression=feather_compression,
            write_tsv=write_tsv,
            gtf_parse_mode=gtf_parse_mode,
            collapse_duplicate_transcripts=collapse_duplicate_transcripts,
            nascent_tolerance=nascent_tolerance,
            splice_blacklist_min_count=splice_blacklist_min_count,
        )
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        genome = _read_genome(fasta_file, alignable_zarr_path, splice_blacklist_min_count)
        logger.info(f"[START] Reading transcripts from {gtf_file}")
        transcripts = read_transcripts(
            gtf_file,
            gtf_parse_mode=gtf_parse_mode,
            collapse_duplicate_transcripts=collapse_duplicate_transcripts,
            ref_lengths=genome.ref_lengths,
        )
        if not transcripts:
            raise ValueError(f"GTF {gtf_file} yields no transcripts; an index needs at least one")
        logger.info(f"[DONE] Read {len(transcripts)} transcripts")

        # The manifest is written last, so removing an old one before the first write leaves an
        # interrupted rebuild with none, which `load` refuses, rather than an old manifest over a
        # mix of old and new files.
        (output_dir / MANIFEST_JSON).unlink(missing_ok=True)
        write = _TableWriter(output_dir, feather_compression, write_tsv)
        write(
            pd.DataFrame(list(genome.ref_lengths.items()), columns=["ref", "length"]),
            REF_LENGTHS_FEATHER,
            REF_LENGTHS_TSV,
        )
        _add_nascent_spans(transcripts, nascent_tolerance)
        _write_annotation_tables(write, transcripts, genome.ref_lengths)
        _write_genome_tables(write, genome)
        _write_manifest(
            output_dir,
            {"fasta": fasta_file, "gtf": gtf_file, "alignable_zarr": alignable_zarr_path},
            # Every parameter of `build` must appear here; a test checks this set against build's signature
            {
                "collapse_duplicate_transcripts": collapse_duplicate_transcripts,
                "feather_compression": feather_compression,
                "gtf_parse_mode": gtf_parse_mode,
                "nascent_tolerance": nascent_tolerance,
                "splice_blacklist_min_count": splice_blacklist_min_count,
                "write_tsv": write_tsv,
            },
            genome,
        )
        logger.info(f"Index written to {output_dir}")

    # -- load (classmethod) ---------------------------------------------------

    @staticmethod
    def _build_gene_table(t_df: pd.DataFrame) -> pd.DataFrame:
        """Derive the gene table from the transcript table.

        Every ``g_index`` gets one row, including each synthetic transcript's gene-neutral placeholder, flagged
        ``is_synthetic`` for callers to filter; keeping them preserves ``g_df.index == g_df["g_index"]``.
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
            g_df[col] = g_df[col].astype("category")
        return g_df

    @classmethod
    def load(
        cls,
        index_dir: str | Path,
        *,
        retain_lookup_maps: bool = False,
    ) -> "TranscriptIndex":
        """Load an index from the Feather files in *index_dir*: the version gate, the references and the mask, the
        transcript and gene tables, the splice graph, and the native resolver built from the intervals and junctions.

        Parameters
        ----------
        retain_lookup_maps : bool
            If True, keep the Python-side copies of what the native resolver owns, ``sj_map`` and
            the per-transcript exon intervals; :meth:`get_exon_intervals` raises without them.
        """
        index_dir = str(index_dir)
        self = cls()
        self.index_dir = index_dir
        self._retain_lookup_maps = retain_lookup_maps
        manifest = _checked_manifest(index_dir)

        self.ref_lengths = load_ref_lengths(os.path.join(index_dir, REF_LENGTHS_FEATHER))
        self.ref_names = list(self.ref_lengths.keys())
        self.ref_name_to_id = {name: i for i, name in enumerate(self.ref_names)}
        from .observable import check_mask

        try:
            self.unobservable_df = check_mask(
                pd.read_feather(os.path.join(index_dir, UNOBSERVABLE_FEATHER)), self.ref_lengths
            )
        except ValueError as exc:
            raise IndexLoadError(
                f"Index at {index_dir}: {exc}. Rebuild the index (rigel index --fasta ... --gtf ... -o {index_dir})."
            ) from exc

        self.t_df = _load_transcript_table(index_dir)
        self.g_df = cls._build_gene_table(self.t_df)
        if not (self.g_df.index == self.g_df["g_index"]).all():
            raise IndexLoadError(
                f"Invalid index in {index_dir}: derived gene table row index "
                f"does not match 'g_index'; rebuild index"
            )
        self.t_to_g_arr = self.t_df["g_index"].values
        self.t_to_strand_arr = self.t_df["strand"].values

        # The splice graph, the calibration partition. Load checks only that the graph is internally consistent and
        # complete; the checks against the transcripts run at build, so a `signature` or `flags` column that has
        # drifted from the annotation loads clean
        from .calibration.splice_graph import load_boundaries, load_regions, validate_graph

        logger.debug("Reading splice graph")
        self.regions_df = load_regions(os.path.join(index_dir, REGIONS_FEATHER))
        self.edges_df = load_boundaries(os.path.join(index_dir, BOUNDARIES_FEATHER))
        validate_graph(self.regions_df, self.edges_df, self.ref_lengths)

        logger.debug("Reading intervals")
        iv_df = pd.read_feather(os.path.join(index_dir, INTERVALS_FEATHER))
        logger.debug("Building collapsed interval index")
        overlap_index = _collapsed_intervals(iv_df)
        self._t_exon_intervals = _exon_intervals(iv_df, self.t_df)
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(os.path.join(index_dir, SJ_FEATHER))
        logger.debug("Building splice junction exact-match map")
        self.sj_map = _sj_exact_match_map(sj_df)

        self.resolver = self._build_resolver(overlap_index, manifest)
        logger.debug("Built native FragmentResolver for C++ resolution")
        return self

    def _build_resolver(self, overlap_index: tuple, manifest: dict):
        """The native ``FragmentResolver``: every index fact fragment resolution reads, installed once, so every use
        of the resolver has them, a scan or not. It owns the junction map afterwards, so the Python copy is freed
        unless the index retains its lookup maps."""
        from .native import FragmentResolver

        resolver = FragmentResolver()
        # The ref-ID space, seeded in ref_names order. The resolver would otherwise number refs in first-seen interval
        # order, and the accumulator deposit, which indexes a partition laid out in ref_names order by the resolver's
        # ref id, would be mis-routed on a multi-reference genome (one reference hides it).
        resolver.set_ref_names([str(name) for name in self.ref_names])
        resolver.build_overlap_index(*overlap_index)
        resolver.build_sj_map(*_sj_map_csr(self.sj_map))
        if not self._retain_lookup_maps:
            self.sj_map = None
        self.sj_blacklist_size, self.sj_blacklist_digest = _install_blacklist(
            resolver, self.index_dir, manifest
        )
        # Each candidate's introns, for its path through the gaps and its length. build_exon_csr() is the one owner of
        # this CSR (scoring reads it too), and it frees ``_t_exon_intervals`` unless retain_lookup_maps.
        exon_offsets, exon_starts_flat, exon_ends_flat, _ = self.build_exon_csr()
        resolver.build_exon_index(
            exon_offsets.tolist(),
            exon_starts_flat.tolist(),
            exon_ends_flat.tolist(),
        )
        # The transcript count, which sizes the resolver's scratch buffers.
        resolver.set_metadata(self.num_transcripts)
        # Each transcript's strand. It sets every fragment's ambiguous-strand flag (its candidates disagree on
        # strand), the strand each gap hypothesis infers from the transcripts that imply it, and the truth label of
        # the all-exonic diagnostic strand model.
        resolver.set_transcript_strands(self.t_to_strand_arr.tolist())
        resolver.set_nrna_parent_index(_nrna_parents(self.t_df).tolist())
        return resolver

    def get_exon_intervals(self, t_idx: int) -> np.ndarray | None:
        """Return sorted exon ``[start, end)`` intervals for a transcript.

        Needs the per-transcript dict, which :meth:`build_exon_csr` frees at load unless the index was
        loaded with ``retain_lookup_maps``.

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
            raise RuntimeError(
                "get_exon_intervals() unavailable: the exon intervals were freed after the exon CSR "
                "was built. Load the index with retain_lookup_maps=True."
            )
        return self._t_exon_intervals.get(t_idx)

    def build_exon_csr(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Build CSR arrays for per-transcript exon positions.

        Flattens the per-transcript exon dict into CSR arrays the native code reads directly. The result is
        cached, and the first call frees ``_t_exon_intervals`` unless the index keeps ``retain_lookup_maps``.

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
        if self._exon_csr_cache is not None:
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
        if not self._retain_lookup_maps:
            self._t_exon_intervals = None
        result = (offsets, starts, ends, cumsum_before)
        self._exon_csr_cache = result

        return result
