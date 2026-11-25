"""Build, load and query the rigel reference index.

The index is built from a genome FASTA (with a samtools .fai) and a GTF annotation. It writes Feather files, with
optional TSV mirrors, to an output directory:

    ref_lengths.feather        — reference names and lengths
    transcripts.feather        — one row per transcript with integer indices
    intervals.feather          — each transcript's exons and span, and the intergenic gaps
    regions.feather            — the splice graph's regions (the calibration partition)
    edges.feather              — the splice graph's boundaries: contiguous and splice-junction
    sj.feather                 — annotated splice junctions from transcript introns
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
from typing import Iterator

import numpy as np
import pandas as pd
import pysam

from .config import CONSTANTS, GtfParseMode, IndexConfig
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

EDGES_FEATHER = "edges.feather"
EDGES_TSV = "edges.tsv"

SJ_FEATHER = "sj.feather"
SJ_TSV = "sj.tsv"

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
INDEX_FORMAT_VERSION = 13


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
        listed = CONSTANTS.qc.examples(
            (f"{name} ({ref_lengths[name]:,} bp)" for name in too_long), sep="; "
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

    names = CONSTANTS.qc.examples((where(t) for t in off), sep="; ")
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
    if gtf_parse_mode == "strict" or ref_lengths is None:
        names = CONSTANTS.qc.examples(
            (f"{g_id} on {', '.join(on)}" for g_id, on in split.items()), sep="; "
        )
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
    names = CONSTANTS.qc.examples(
        (
            f"{g_id} kept on {kept[g_id]}, dropped on {', '.join(ref for ref in on if ref != kept[g_id])}"
            for g_id, on in split.items()
        ),
        sep="; ",
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
    gtf_parse_mode: GtfParseMode = "strict",
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
    gtf_file : str or Path
        The annotation.
    gtf_parse_mode : {"strict", "warn-skip"}
        ``"strict"`` refuses a malformed record; ``"warn-skip"`` logs and skips it.
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


# -- Synthetic nascent spans: one entity per extent within a tolerance at both ends ---


def _linked_components(points: list[tuple[int, int]], tolerance: int) -> list[int]:
    """Each ``(start, end)`` point's component: two points whose starts and ends each lie within *tolerance* bp are
    linked, and a component is a connected set of links (union-find over a start-sorted sweep)."""
    parent = list(range(len(points)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    order = sorted(range(len(points)), key=lambda i: points[i])
    for a, i in enumerate(order):
        for j in order[a + 1 :]:
            if points[j][0] - points[i][0] > tolerance:
                break
            if abs(points[j][1] - points[i][1]) <= tolerance:
                parent[find(i)] = find(j)
    return [find(i) for i in range(len(points))]


def _entity_labels(points: list[tuple[int, int]], tolerance: int) -> list[int]:
    """Each point's entity: its component (:func:`_linked_components`), then the components' envelopes linked again
    until no two envelopes lie within *tolerance* at both ends. An envelope takes its start from one member and its
    end from another, so two components no member of which is near the other's can still have near envelopes."""
    labels = _linked_components(points, tolerance)
    while True:
        ids = sorted(set(labels))
        position = {label: k for k, label in enumerate(ids)}
        envelope = [[np.iinfo(np.int64).max, np.iinfo(np.int64).min] for _ in ids]
        for (start, end), label in zip(points, labels):
            env = envelope[position[label]]
            env[0], env[1] = min(env[0], start), max(env[1], end)
        merged = _linked_components([(a, b) for a, b in envelope], tolerance)
        if len(set(merged)) == len(ids):
            return labels
        labels = [merged[position[label]] for label in labels]


def _merged_spans(
    transcripts: list[Transcript], tolerance: int
) -> tuple[dict[int, tuple], dict[tuple, int], dict[tuple, Transcript]]:
    """The nascent entities. On each reference and strand, every multi-exon transcript's extent and every
    single-exon transcript are points; points within *tolerance* at both ends are one entity, transitively
    (:func:`_entity_labels`). An entity holding a multi-exon transcript is a span, its multi-exon members'
    envelope. A single-exon transcript of the entity whose ends lie within *tolerance* of the span's own is that
    span: a single-exon transcript is its own nascent RNA, so no span equal to one is made; where several are, the one
    nearest the span, then the smallest ``t_id``. One that is near a member but not the span (a chain's far end) is not
    the span's equal and stands for nothing. Returns each multi-exon transcript's span key ``(ref, strand, start,
    end)`` by ``t_index``, each span's contributor count, and each span's single-exon transcript where it has one."""
    places: dict[tuple[str, int], list[Transcript]] = collections.defaultdict(list)
    for t in transcripts:
        places[(t.ref, int(t.strand))].append(t)
    t_to_span_key: dict[int, tuple] = {}
    contributors: dict[tuple, int] = {}
    single_exon: dict[tuple, Transcript] = {}
    for (ref, strand), members in places.items():
        labels = _entity_labels([(t.start, t.end) for t in members], tolerance)
        entities: dict[int, list[Transcript]] = collections.defaultdict(list)
        for t, label in zip(members, labels):
            entities[label].append(t)
        for entity in entities.values():
            multi = [t for t in entity if len(t.exons) >= 2]
            if not multi:
                continue
            key = (ref, strand, min(t.start for t in multi), max(t.end for t in multi))
            for t in multi:
                t_to_span_key[t.t_index] = key
            contributors[key] = len(multi)
            singles = [
                t
                for t in entity
                if len(t.exons) == 1
                and max(abs(t.start - key[2]), abs(t.end - key[3])) <= tolerance
            ]
            if singles:
                single_exon[key] = min(
                    singles, key=lambda t: (max(abs(t.start - key[2]), abs(t.end - key[3])), t.t_id)
                )
    return t_to_span_key, contributors, single_exon


def _containing_transcripts(transcripts: list[Transcript], spans) -> dict[tuple, Transcript]:
    """The single-exon transcript that contains each span whole, beyond the tolerance (the first, by start then
    end, on the span's reference and strand): it carries the span too."""
    by_place: dict[tuple[str, int], list[tuple[int, int, Transcript]]] = collections.defaultdict(
        list
    )
    for t in transcripts:
        if len(t.exons) == 1:
            by_place[(t.ref, int(t.strand))].append((t.start, t.end, t))
    for candidates in by_place.values():
        candidates.sort(key=lambda x: (x[0], x[1]))
    containing: dict[tuple, Transcript] = {}
    for span_key in spans:
        ref, strand, m_start, m_end = span_key
        for s_start, s_end, s_tx in by_place.get((ref, strand), []):
            if s_start > m_start:
                break  # sorted by start, no more can contain
            if m_end <= s_end:
                containing[span_key] = s_tx
                break
    return containing


#: Leads every synthetic nascent span's ``t_id``, followed by ``{ref}_{strand}_{start}_{end}``.
NASCENT_SPAN_PREFIX = "RIGEL_NASCENT_"


def _synthetic_span(span_key: tuple, n_contributors: int) -> Transcript:
    """The synthetic single-exon transcript over a merged span. It is gene-neutral: its sentinel ``g_id`` is its own
    ``t_id``, which keeps the ``g_id`` column non-null, and the index builder gives it a ``g_index`` of its own."""
    ref, strand, start, end = span_key
    t_id = f"{NASCENT_SPAN_PREFIX}{ref}_{strand}_{start}_{end}"
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


def create_nascent_spans(
    transcripts: list[Transcript],
    tolerance: int = IndexConfig.end_tolerance,
) -> tuple[list[Transcript], dict[int, tuple], dict[tuple, Transcript]]:
    """Create the synthetic nascent spans and detect their annotated equivalents.

    1. Each multi-exon transcript's extent and each single-exon transcript, on one reference and strand, are one
       entity with every other within *tolerance* bp at both ends, transitively (:func:`_merged_spans`); a span is
       its multi-exon members' envelope.
    2. A single-exon transcript of a span's entity whose ends lie within *tolerance* of the span's own is the span: a
       single-exon transcript is its own nascent RNA. A span no single-exon transcript is near is carried by one that
       contains it whole, if any
       (:func:`_containing_transcripts`). Either way its ``span_n_contributors`` counts the span's transcripts.
    3. Each span not carried becomes a synthetic single-exon ``Transcript`` (``is_synthetic``, ``is_single_exon``).
       Synthetics are gene-neutral (:func:`_synthetic_span`): a merged span may draw on transcripts of several genes,
       so inheriting one contributor's gene would be arbitrary; gene summaries are derived downstream from annotated
       contributors only.

    *transcripts* is the sorted transcript list, ``t_index`` already assigned. Returns ``(synthetics, t_to_span_key,
    carriers)``: the synthetics to append, in their spans' first-seen order (the caller assigns their ``t_index``
    and ``g_index``); each annotated multi-exon transcript's span key ``(ref, strand, start, end)``, by ``t_index``,
    from which the caller sets ``span_t_index``; and the annotated single-exon transcript carrying each carried span.
    """
    if not any(len(t.exons) >= 2 for t in transcripts):
        return [], {}, {}
    t_to_span_key, contributors, carriers = _merged_spans(transcripts, tolerance)
    carriers.update(
        _containing_transcripts(transcripts, [k for k in contributors if k not in carriers])
    )
    for span_key, s_tx in carriers.items():
        s_tx.span_n_contributors += contributors[span_key]
    synthetics = [
        _synthetic_span(span_key, n)
        for span_key, n in contributors.items()
        if span_key not in carriers
    ]
    logger.info(
        f"Nascent: {len(contributors)} merged spans "
        f"({len(carriers)} annotated equiv, "
        f"{len(synthetics)} synthetic)"
    )
    return synthetics, t_to_span_key, carriers


def _add_nascent_spans(transcripts: list[Transcript], tolerance: int) -> None:
    """Append the synthetic spans that model the transcripts' nascent RNA (:func:`create_nascent_spans`) and set
    every span entity's ``span_t_index``. Each synthetic gets a ``g_index`` of its own (they are gene-neutral), which
    keeps ``t.g_index == row in g_df`` and keeps synthetics out of annotated genes' summaries."""
    logger.info("[START] Creating synthetic nascent spans")
    synthetics, t_to_span_key, carriers = create_nascent_spans(transcripts, tolerance=tolerance)
    next_g_index = max((t.g_index for t in transcripts), default=-1) + 1
    span_entity_t_index: dict[tuple, int] = {}
    for i, syn in enumerate(synthetics):
        syn.t_index = len(transcripts) + i
        syn.g_index = next_g_index + i
        span_entity_t_index[(syn.ref, int(syn.strand), syn.start, syn.end)] = syn.t_index
    # A carried span's entity is the annotated single-exon transcript that carries it: near its ends, or containing it
    for span_key, equiv_tx in carriers.items():
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
    the full transcript ``[first_exon.start, last_exon.end)``. Only annotated
    transcripts reach it (:func:`iter_reference_layout` filters the synthetic
    spans out): the C++ resolver derives a synthetic candidate from the
    transcripts a fragment hits, through each one's ``span_t_index``.
    """
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
    is strand-agnostic and excludes the synthetic spans, which
    :func:`iter_reference_layout` filters out.
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
    The synthetic spans never reach it: :func:`iter_reference_layout` filters them out.
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
    """What the build reads of the genome: the references and the unobservable mask with the manifest's record of
    what it is (`rigel.observable`)."""

    ref_lengths: dict[str, int]
    mask: pd.DataFrame
    mask_record: dict

    @property
    def masked_bases(self) -> int:
        return int((self.mask["end"] - self.mask["start"]).sum())

    @property
    def genome_bases(self) -> int:
        return int(sum(self.ref_lengths.values()))


def _read_genome(fasta_file: str | Path, alignable_zarr_path: str | Path | None) -> _GenomeInputs:
    """The references and the mask, each read and checked before the first write: a store that cannot be read, is
    incomplete, or was built for other references or another assembly refuses the build and leaves the index it would
    replace intact."""
    logger.info(f"[START] Loading references from {fasta_file}")
    ref_lengths = load_reference_lengths(fasta_file)
    logger.info(f"[DONE] Read {len(ref_lengths)} references")

    from .observable import mask_from_fasta, mask_from_store

    if alignable_zarr_path is None:
        logger.info("[START] Recording the reference's unknown runs (no alignable store)")
        mask = mask_from_fasta(fasta_file, ref_lengths)
        record = {"rule": "the reference's unknown bases (N or n)", "read_length": None}
    else:
        logger.info(f"[START] Recording the unobservable mask from {alignable_zarr_path}")
        mask, record = mask_from_store(alignable_zarr_path, ref_lengths, fasta_file)
    genome = _GenomeInputs(ref_lengths, mask, record)
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
    write(edges_df, EDGES_FEATHER, EDGES_TSV)


def _write_genome_tables(write: _TableWriter, genome: _GenomeInputs) -> None:
    """The unobservable mask."""
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

#: The files every index has, checked before any is read. The alias table is never read.
_REQUIRED_FILES = (
    REF_LENGTHS_FEATHER,
    TRANSCRIPTS_FEATHER,
    REGIONS_FEATHER,
    EDGES_FEATHER,
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


def _ref_codes(refs) -> np.ndarray:
    """Each row's reference as its rank among the distinct names in sorted order, the codes ``np.unique`` gives, by
    hashing rather than sorting the strings."""
    return pd.factorize(pd.Series(refs), sort=True)[0].astype(np.int64)


def _run_members(
    order: np.ndarray,
    first: np.ndarray,
    last: np.ndarray,
    values: np.ndarray,
    *,
    drop_negative: bool,
):
    """``(flat, offsets)``: each run's distinct *values* in ascending order, the negative ones dropped where asked, as a
    CSR over the runs. The runs are ``_key_runs``' over the rows *order* sorts."""
    run = np.repeat(np.arange(first.shape[0], dtype=np.int64), last - first)
    v = np.asarray(values)[order]
    if drop_negative:
        kept = v >= 0
        run, v = run[kept], v[kept]
    by = np.lexsort((v, run))
    run, v = run[by], v[by]
    new = np.ones(v.shape[0], dtype=bool)
    new[1:] = (run[1:] != run[:-1]) | (v[1:] != v[:-1])
    run, v = run[new], v[new]
    offsets = np.zeros(first.shape[0] + 1, dtype=np.int64)
    np.cumsum(np.bincount(run, minlength=first.shape[0]), out=offsets[1:])
    return v, offsets


def _collapsed_intervals(iv_df: pd.DataFrame) -> tuple[list, list, list, list, list, list]:
    """The resolver's overlap index, as the six lists ``FragmentResolver.build_overlap_index`` takes: one interval per
    unique ``(ref, start, end, interval_type)`` — refs, starts, ends and types — and, in CSR, the sorted transcripts
    of every row it merges, the intergenic sentinel -1 dropped."""
    refs = iv_df["ref"].values
    starts = iv_df["start"].values
    ends = iv_df["end"].values
    itypes = iv_df["interval_type"].values
    order, first, last = _key_runs(_ref_codes(refs), starts, ends, itypes)
    t_flat, t_offsets = _run_members(
        order, first, last, iv_df["t_index"].values, drop_negative=True
    )
    logger.debug(f"Interval index: {len(iv_df)} rows → {len(first)} collapsed")
    return (
        refs[order][first].tolist(),
        starts[order][first].tolist(),
        ends[order][first].tolist(),
        itypes[order][first].tolist(),
        t_flat.tolist(),
        t_offsets.tolist(),
    )


def _exon_csr(
    iv_df: pd.DataFrame, t_df: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(offsets, starts, ends, cumsum_before)``: every transcript's exons as CSR over ``t_index``, int32, sorted by
    start, with each exon's cumulative exon length before it in its transcript. Synthetics are absent from the interval
    table (:func:`iter_reference_layout` filters them out), but the native code still needs their single exon,
    ``[start, end)`` from the transcript table; a synthetic's own span replaces anything the table holds for it."""
    n_t = int(t_df.shape[0])
    syn = t_df["is_synthetic"].to_numpy(dtype=bool)
    tidx = iv_df["t_index"].values
    exon = (iv_df["interval_type"].values == int(IntervalType.EXON)) & (tidx >= 0)
    exon[exon] = ~syn[tidx[exon]]
    syn_t = np.flatnonzero(syn)
    t = np.concatenate([tidx[exon].astype(np.int64), syn_t])
    start = np.concatenate(
        [iv_df["start"].values[exon].astype(np.int32), t_df["start"].to_numpy(dtype=np.int32)[syn]]
    )
    end = np.concatenate(
        [iv_df["end"].values[exon].astype(np.int32), t_df["end"].to_numpy(dtype=np.int32)[syn]]
    )
    order = np.lexsort((start, t))
    t, start, end = t[order], start[order], end[order]
    offsets = np.zeros(n_t + 1, dtype=np.int32)
    np.cumsum(np.bincount(t, minlength=n_t), out=offsets[1:])
    length = (end - start).astype(np.int64)
    within = np.cumsum(length) - length
    cumsum_before = (within - within[offsets[:-1]][t] if t.size else within).astype(np.int32)
    return offsets, start, end, cumsum_before


def _sj_csr(sj_df: pd.DataFrame) -> tuple[list, list, list, list, list, list]:
    """Every annotated junction with the transcripts holding it, as the six lists ``FragmentResolver.build_sj_map``
    takes: refs, starts, ends, strands and, in CSR, each junction's sorted transcripts; the junctions in
    ``(ref, start, end, strand)`` order."""
    refs = sj_df["ref"].values
    starts = sj_df["start"].values
    ends = sj_df["end"].values
    strands = sj_df["strand"].values
    order, first, last = _key_runs(_ref_codes(refs), starts, ends, strands)
    t_flat, t_offsets = _run_members(
        order, first, last, sj_df["t_index"].values, drop_negative=False
    )
    logger.debug(f"Splice junctions: {len(first)} unique, {len(sj_df)} total")
    return (
        refs[order][first].tolist(),
        starts[order][first].tolist(),
        ends[order][first].tolist(),
        strands[order][first].tolist(),
        t_flat.tolist(),
        t_offsets.tolist(),
    )


def _nascent_parents(t_df: pd.DataFrame) -> np.ndarray:
    """Each transcript's synthetic span, by ``t_index``, -1 where its span entity is no synthetic. Synthetics are absent
    from the interval table (:func:`iter_reference_layout` filters them out), so the resolver derives each synthetic
    candidate from the annotated transcripts a fragment hits, through this index."""
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
        #: every transcript's exons as flat CSR arrays (:func:`_exon_csr`), built once by :meth:`load`.
        self._exon_csr_cache: tuple | None = None

    # -- properties -----------------------------------------------------------

    @property
    def num_transcripts(self) -> int:
        return len(self.t_df)

    @property
    def num_genes(self) -> int:
        """Rows in ``g_df``, annotated genes and synthetic nascent spans, so ``g_index == row``.

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
        transcripts hold a fragment's blocks and junctions, and the junctions decide what it offers, so two indexes
        that share this hash can tally one BAM differently. ``edges.feather`` can change
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
        the offsets over each low end, the high ends' region-bound indices and the annotated strands. A payload
        records it, so a cached tally can be checked against the index it is read with. It is not all a scan depends
        on: the transcript models resolution reads and the scan settings are not in it. Computed on demand, never
        stored.
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
        gtf_parse_mode: GtfParseMode = IndexConfig.gtf_parse_mode,
        collapse_duplicate_transcripts: bool = IndexConfig.collapse_duplicate_transcripts,
        end_tolerance: int = IndexConfig.end_tolerance,
        alignable_zarr_path: str | Path | None = None,
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
        end_tolerance : int
            Nascent entities whose start and end sites each lie within this many bp are one
            (``IndexConfig.end_tolerance``): the synthetic spans that model the transcripts' nascent RNA are merged
            so, and none is made equal to a single-exon transcript.
        alignable_zarr_path : path, optional
            Path to an alignable Zarr store built for the same genome and aligner, whose zero set is the
            unobservable mask (`rigel.observable`).  When ``None``, the mask is the FASTA's unknown runs.
        """
        IndexConfig(
            feather_compression=feather_compression,
            write_tsv=write_tsv,
            gtf_parse_mode=gtf_parse_mode,
            collapse_duplicate_transcripts=collapse_duplicate_transcripts,
            end_tolerance=end_tolerance,
        )
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        genome = _read_genome(fasta_file, alignable_zarr_path)
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
        _add_nascent_spans(transcripts, end_tolerance)
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
                "end_tolerance": end_tolerance,
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
    def load(cls, index_dir: str | Path) -> "TranscriptIndex":
        """Load an index from the Feather files in *index_dir*: the version gate, the references and the mask, the
        transcript and gene tables, the splice graph, and the native resolver built from the intervals and junctions.
        """
        index_dir = str(index_dir)
        self = cls()
        self.index_dir = index_dir
        _checked_manifest(index_dir)

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
        from .calibration.splice_graph import load_edges, load_regions, validate_graph

        logger.debug("Reading splice graph")
        self.regions_df = load_regions(os.path.join(index_dir, REGIONS_FEATHER))
        self.edges_df = load_edges(os.path.join(index_dir, EDGES_FEATHER))
        validate_graph(self.regions_df, self.edges_df, self.ref_lengths)

        logger.debug("Reading intervals")
        iv_df = pd.read_feather(os.path.join(index_dir, INTERVALS_FEATHER))
        logger.debug("Building collapsed interval index")
        overlap_index = _collapsed_intervals(iv_df)
        self._exon_csr_cache = _exon_csr(iv_df, self.t_df)
        logger.debug("Reading splice junctions")
        sj_df = pd.read_feather(os.path.join(index_dir, SJ_FEATHER))
        logger.debug("Building splice junction exact-match map")
        sj_csr = _sj_csr(sj_df)

        self.resolver = self._build_resolver(overlap_index, sj_csr)
        logger.debug("Built native FragmentResolver for C++ resolution")
        return self

    def _build_resolver(self, overlap_index: tuple, sj_csr: tuple):
        """The native ``FragmentResolver``: every index fact fragment resolution reads, installed once, so every use
        of the resolver has them, a scan or not. It owns the junction map."""
        from .native import FragmentResolver

        resolver = FragmentResolver()
        # The ref-ID space, seeded in ref_names order. The resolver would otherwise number refs in first-seen interval
        # order, and the accumulator deposit, which indexes a partition laid out in ref_names order by the resolver's
        # ref id, would be mis-routed on a multi-reference genome (one reference hides it).
        resolver.set_ref_names([str(name) for name in self.ref_names])
        resolver.build_overlap_index(*overlap_index)
        resolver.build_sj_map(*sj_csr)
        # Each candidate's introns, for its path through the gaps and its length: the exon CSR, which scoring reads too
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
        resolver.set_nascent_parent_index(_nascent_parents(self.t_df).tolist())
        return resolver

    def get_exon_intervals(self, t_idx: int) -> np.ndarray:
        """Transcript ``t_idx``'s exons, ``(n_exons, 2)`` int32 ``[start, end)`` sorted by start: its slice of the exon
        CSR (:meth:`build_exon_csr`). A ``t_idx`` outside the transcript table raises ``IndexError``."""
        if not 0 <= t_idx < self.num_transcripts:
            raise IndexError(
                f"transcript index {t_idx} is outside the index's {self.num_transcripts} transcripts"
            )
        offsets, starts, ends, _ = self.build_exon_csr()
        lo, hi = int(offsets[t_idx]), int(offsets[t_idx + 1])
        return np.column_stack((starts[lo:hi], ends[lo:hi]))

    def build_exon_csr(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """The per-transcript exon positions as the CSR arrays the native code reads directly, built once by
        :meth:`load` (:func:`_exon_csr`).

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
        return self._exon_csr_cache
