"""The simulator's truth tables, read off the read names (``rigel.sim.read_name``). A helper module, not a script.

There are two transcript tables, both written by :func:`_write_transcript_table` from :func:`_transcript_metadata`:

- the molar table, one per nascent setting (:func:`write_molar_truth`): the molecules each transcript was simulated
  at, under the simulator's input names, so the file is also a valid abundance input;
- one per condition (:func:`write_post_capture_truth`): the same molecules (``pre_capture_*``) beside the fragments
  the condition realised (``observed_*_fragments``), with a fragment-length table and a summary JSON.
"""

from __future__ import annotations

import gzip
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Iterable

import pysam

from rigel.sim.read_name import Origin, OriginKind, parse_origin

logger = logging.getLogger(__name__)

#: The molar table's columns. ``total_rna`` is the mature plus the nascent molecules.
_MOLAR_COLUMNS = (
    "transcript_id",
    "gene_id",
    "gene_name",
    "ref",
    "strand",
    "annotated_abundance",
    "nascent_abundance",
    "total_rna",
    "n_exons",
    "spliced_length",
    "genomic_span",
)

#: A condition's table's columns: the molecules before capture, then the fragments realised.
_POST_CAPTURE_COLUMNS = (
    "transcript_id",
    "gene_id",
    "gene_name",
    "ref",
    "strand",
    "pre_capture_annotated_abundance",
    "pre_capture_nascent_abundance",
    "pre_capture_total_rna",
    "observed_annotated_fragments",
    "observed_nascent_fragments",
    "observed_total_rna_fragments",
    "n_exons",
    "spliced_length",
    "genomic_span",
)

#: The fragment-length table's origin pools, in the order it lists them.
_LENGTH_POOLS = ("annotated", "nascent", "rna", "gdna", "all")


# ── the read names ───────────────────────────────────────────────────────────────────────────────


def _bam_qnames(bam_path: Path) -> set[str]:
    """One query name per fragment of a name-sorted or unsorted BAM."""
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        return {read.query_name for read in bam}


def _fastq_qnames(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 0:
                yield line


def _origins(*, bam_path: Path | None, fastq_path: Path | None) -> Iterable[Origin]:
    """Every fragment's origin, from the BAM when it exists, else from the FASTQ."""
    if bam_path is not None and bam_path.exists():
        qnames: Iterable[str] = _bam_qnames(bam_path)
    elif fastq_path is not None and fastq_path.exists():
        qnames = _fastq_qnames(fastq_path)
    else:
        raise FileNotFoundError("post-capture truth requires an existing BAM or FASTQ source")
    return (parse_origin(qname) for qname in qnames)


# ── the tables ───────────────────────────────────────────────────────────────────────────────────


def _length_stats(length_counts: Counter[int]) -> dict[str, float | int | None]:
    n = int(sum(length_counts.values()))
    if n == 0:
        return {"n": 0, "mean": None, "std": None, "min": None, "max": None}
    lengths = sorted(length_counts)
    mean = float(sum(length * count for length, count in length_counts.items())) / float(n)
    var = sum(((length - mean) ** 2) * count for length, count in length_counts.items()) / float(n)
    return {"n": n, "mean": mean, "std": var**0.5, "min": int(lengths[0]), "max": int(lengths[-1])}


def _transcript_metadata(transcript) -> dict[str, object]:
    """A transcript's identity and structure columns and the molecules it was simulated at (mature, nascent and
    their sum), as both truth tables write them."""
    length = transcript.length or transcript.compute_length()
    # 0 for a transcript with no exon, whose start and end are both 0
    genomic_span = transcript.end - transcript.start
    pre_capture_annotated = float(transcript.abundance or 0.0)
    pre_capture_nascent = float(transcript.nascent_abundance)
    return {
        "transcript_id": transcript.t_id,
        "gene_id": transcript.g_id,
        "gene_name": transcript.g_name,
        "ref": transcript.ref,
        "strand": transcript.strand.to_str(),
        "n_exons": len(transcript.exons),
        "spliced_length": int(length),
        "genomic_span": int(genomic_span),
        "pre_capture_annotated_abundance": pre_capture_annotated,
        "pre_capture_nascent_abundance": pre_capture_nascent,
        "pre_capture_total_rna": pre_capture_annotated + pre_capture_nascent,
    }


def _write_transcript_table(path: Path, columns: tuple[str, ...], rows: Iterable[dict]) -> None:
    """Write a header of ``columns``, then one tab-separated line per row, each value as ``str`` gives it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in columns) + "\n")


def _molar_row(transcript) -> dict[str, object]:
    meta = _transcript_metadata(transcript)
    return {
        **meta,
        "annotated_abundance": f"{meta['pre_capture_annotated_abundance']:.4f}",
        "nascent_abundance": f"{meta['pre_capture_nascent_abundance']:.4f}",
        "total_rna": f"{meta['pre_capture_total_rna']:.4f}",
    }


def write_molar_truth(transcripts, path: Path) -> None:
    """Write the molar truth table: the molecules each transcript was simulated at, mature and nascent, before
    capture and sampling, to four decimals. They depend only on the nascent setting, so the grid writes one table
    per setting. The columns are the ``sim`` abundance format (``_sim_abundance``), so the table can be given back
    as ``abundance.file``. Gates: ``tests/test_sim_truth.py``."""
    _write_transcript_table(path, _MOLAR_COLUMNS, (_molar_row(t) for t in transcripts))
    logger.info("Wrote truth abundances to %s", path)


def summarize_post_capture_truth(
    transcripts,
    *,
    bam_path: Path | None = None,
    fastq_path: Path | None = None,
) -> dict[str, object]:
    """Count the fragments each transcript realised, by origin, from the read names of a BAM or FASTQ.

    Under capture these realised counts are the abundance truth, since capture changes which templates and starts
    are drawn; the pre-capture molecules are kept beside them as provenance.
    """
    by_transcript: dict[OriginKind, Counter[str]] = {"annotated": Counter(), "nascent": Counter()}
    origin_counts: Counter[OriginKind] = Counter()
    length_counts: dict[str, Counter[int]] = {pool: Counter() for pool in _LENGTH_POOLS}

    for origin in _origins(bam_path=bam_path, fastq_path=fastq_path):
        origin_counts[origin.kind] += 1
        if origin.kind in by_transcript and origin.transcript_id is not None:
            by_transcript[origin.kind][origin.transcript_id] += 1
        if origin.start is None or origin.end is None or origin.end <= origin.start:
            continue
        length = int(origin.end) - int(origin.start)
        length_counts[origin.kind][length] += 1
        if origin.kind in by_transcript:
            length_counts["rna"][length] += 1
        length_counts["all"][length] += 1

    transcript_rows = []
    for transcript in transcripts:
        meta = _transcript_metadata(transcript)
        annotated = int(by_transcript["annotated"].get(str(meta["transcript_id"]), 0))
        nascent = int(by_transcript["nascent"].get(str(meta["transcript_id"]), 0))
        transcript_rows.append(
            {
                **meta,
                "observed_annotated_fragments": annotated,
                "observed_nascent_fragments": nascent,
                "observed_total_rna_fragments": annotated + nascent,
            }
        )

    return {
        "transcript_rows": transcript_rows,
        "origin_counts": {kind: int(origin_counts.get(kind, 0)) for kind in ("annotated", "nascent", "gdna")},
        "fragment_length_counts": length_counts,
        "fragment_lengths": {pool: _length_stats(counts) for pool, counts in length_counts.items()},
    }


def _write_fragment_length_truth(summary: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    length_counts = summary["fragment_length_counts"]  # type: ignore[index]
    with open(path, "w") as handle:
        handle.write("kind\tfragment_length\tcount\tfraction\n")
        for pool in _LENGTH_POOLS:
            counts = length_counts[pool]
            total = float(sum(counts.values()))
            for length in sorted(counts):
                count = int(counts[length])
                fraction = count / total if total > 0.0 else 0.0
                handle.write(f"{pool}\t{int(length)}\t{count}\t{fraction:.12g}\n")


def write_post_capture_truth(
    transcripts,
    abundance_path: Path,
    fl_path: Path,
    summary_path: Path,
    *,
    bam_path: Path | None = None,
    fastq_path: Path | None = None,
    condition: str | None = None,
    pre_capture_abundances: str | None = None,
    gdna_strand_overdispersion: float | None = None,
) -> dict[str, object]:
    """Write one condition's truth from the origins its read names carry.

    ``abundance_path`` gets one row per transcript, the molecules it was simulated at (``pre_capture_*``) beside the
    fragments realised (``observed_*_fragments``); ``fl_path`` the realised fragment lengths per origin; and
    ``summary_path`` the origin counts, the length statistics and the three tables' paths
    (``pre_capture_abundances``, ``post_capture_abundances``, ``post_capture_fragment_lengths``).
    """
    summary = summarize_post_capture_truth(transcripts, bam_path=bam_path, fastq_path=fastq_path)
    _write_transcript_table(abundance_path, _POST_CAPTURE_COLUMNS, summary["transcript_rows"])  # type: ignore[arg-type]
    _write_fragment_length_truth(summary, fl_path)

    json_summary = {
        "condition": condition,
        "truth_kind": "post_capture_empirical",
        "gdna_strand_overdispersion_true": gdna_strand_overdispersion,
        "pre_capture_abundances": pre_capture_abundances,
        "post_capture_abundances": str(abundance_path),
        "post_capture_fragment_lengths": str(fl_path),
        "origin_counts": summary["origin_counts"],
        "fragment_lengths": summary["fragment_lengths"],
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as handle:
        json.dump(json_summary, handle, indent=2)
    return json_summary
