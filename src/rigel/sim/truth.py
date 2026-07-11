"""Fragment-origin aggregation + post-capture truth-table writing.

The read-name *parsing* (``Origin`` / ``parse_origin``) lives in :mod:`read_name`; this module
consumes it to count origins and write per-condition truth tables. ``Origin``/``OriginKind``/
``parse_origin`` are re-exported here for the many existing ``from rigel.sim.truth import ...``
call sites.
"""

from __future__ import annotations

import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

import pysam

from .read_name import Origin, OriginKind, parse_origin  # re-exported (see module docstring)

__all__ = [
    "Origin",
    "OriginKind",
    "count_mrna_by_transcript_from_bam",
    "count_mrna_by_transcript_from_fastq",
    "count_origins_from_bam",
    "count_origins_from_fastq",
    "deduplicate_bam_qnames",
    "parse_origin",
    "summarize_post_capture_truth",
    "write_post_capture_truth",
]


def deduplicate_bam_qnames(bam_path: Path) -> list[str]:
    """Return one query name per fragment from a name-sorted or unsorted BAM."""
    seen: set[str] = set()
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for read in bam:
            seen.add(read.query_name)
    return list(seen)


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _iter_fastq_qnames(path: Path) -> Iterable[str]:
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 0:
                yield line


def _iter_origins_from_source(
    *,
    bam_path: Path | None = None,
    fastq_path: Path | None = None,
) -> Iterable[Origin]:
    if bam_path is not None and bam_path.exists():
        for qname in deduplicate_bam_qnames(bam_path):
            yield parse_origin(qname)
        return
    if fastq_path is not None and fastq_path.exists():
        for qname in _iter_fastq_qnames(fastq_path):
            yield parse_origin(qname)
        return
    raise FileNotFoundError("post-capture truth requires an existing BAM or FASTQ source")


def _count_origins(qnames: Iterable[str]) -> Counter[OriginKind]:
    counts: Counter[OriginKind] = Counter()
    for qname in qnames:
        counts[parse_origin(qname).kind] += 1
    return counts


def _count_mrna_by_transcript(qnames: Iterable[str]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for qname in qnames:
        origin = parse_origin(qname)
        if origin.kind == "mrna" and origin.transcript_id is not None:
            counts[origin.transcript_id] += 1
    return counts


def count_origins_from_fastq(path: Path) -> Counter[OriginKind]:
    """Count mRNA, nRNA, and gDNA fragments from FASTQ read names."""
    return _count_origins(_iter_fastq_qnames(path))


def count_origins_from_bam(path: Path) -> Counter[OriginKind]:
    """Count mRNA, nRNA, and gDNA fragments from BAM query names."""
    return _count_origins(deduplicate_bam_qnames(path))


def count_mrna_by_transcript_from_fastq(path: Path) -> Counter[str]:
    """Count mRNA fragments per transcript from FASTQ read names."""
    return _count_mrna_by_transcript(_iter_fastq_qnames(path))


def count_mrna_by_transcript_from_bam(path: Path) -> Counter[str]:
    """Count mRNA fragments per transcript from BAM query names."""
    return _count_mrna_by_transcript(deduplicate_bam_qnames(path))


def _origin_length(origin: Origin) -> int | None:
    if origin.start is None or origin.end is None:
        return None
    length = int(origin.end) - int(origin.start)
    return length if length > 0 else None


def _length_stats(length_counts: Counter[int]) -> dict[str, float | int | None]:
    n = int(sum(length_counts.values()))
    if n == 0:
        return {
            "n": 0,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
        }
    lengths = sorted(length_counts)
    total = float(sum(length * count for length, count in length_counts.items()))
    mean = total / float(n)
    var = sum(((length - mean) ** 2) * count for length, count in length_counts.items()) / float(n)
    return {
        "n": n,
        "mean": mean,
        "std": var ** 0.5,
        "min": int(lengths[0]),
        "max": int(lengths[-1]),
    }


def _transcript_metadata(transcript) -> dict[str, object]:
    length = transcript.length or transcript.compute_length()
    genomic_span = transcript.end - transcript.start if transcript.end and transcript.start else 0
    pre_capture_mrna = float(transcript.abundance or 0.0)
    pre_capture_nrna = float(transcript.nrna_abundance)
    return {
        "transcript_id": transcript.t_id,
        "gene_id": transcript.g_id,
        "gene_name": transcript.g_name,
        "ref": transcript.ref,
        "strand": transcript.strand.to_str(),
        "n_exons": len(transcript.exons),
        "spliced_length": int(length),
        "genomic_span": int(genomic_span),
        "pre_capture_mrna_abundance": pre_capture_mrna,
        "pre_capture_nrna_abundance": pre_capture_nrna,
        "pre_capture_total_rna": pre_capture_mrna + pre_capture_nrna,
    }


def summarize_post_capture_truth(
    transcripts,
    *,
    bam_path: Path | None = None,
    fastq_path: Path | None = None,
) -> dict[str, object]:
    """Summarize empirical post-capture truth from simulated read origins.

    The simulator encodes the selected template and interval in every read name.
    For hybrid-capture benchmarking, these observed fragment origins are the
    relevant abundance truth: capture changes which transcripts and starts are
    sampled, so molecular pre-capture abundances are retained only as provenance.
    """
    mrna_by_tx: Counter[str] = Counter()
    nrna_by_tx: Counter[str] = Counter()
    origin_counts: Counter[OriginKind] = Counter()
    length_counts: dict[str, Counter[int]] = {
        "mrna": Counter(),
        "nrna": Counter(),
        "rna": Counter(),
        "gdna": Counter(),
        "all": Counter(),
    }

    for origin in _iter_origins_from_source(bam_path=bam_path, fastq_path=fastq_path):
        origin_counts[origin.kind] += 1
        if origin.kind == "mrna" and origin.transcript_id is not None:
            mrna_by_tx[origin.transcript_id] += 1
        elif origin.kind == "nrna" and origin.transcript_id is not None:
            nrna_by_tx[origin.transcript_id] += 1

        length = _origin_length(origin)
        if length is not None:
            length_counts[origin.kind][length] += 1
            if origin.kind in {"mrna", "nrna"}:
                length_counts["rna"][length] += 1
            length_counts["all"][length] += 1

    transcript_rows = []
    for transcript in transcripts:
        meta = _transcript_metadata(transcript)
        transcript_id = str(meta["transcript_id"])
        mrna_count = int(mrna_by_tx.get(transcript_id, 0))
        nrna_count = int(nrna_by_tx.get(transcript_id, 0))
        total_count = mrna_count + nrna_count
        transcript_rows.append({
            **meta,
            "post_capture_mrna_fragments": mrna_count,
            "post_capture_nrna_fragments": nrna_count,
            "post_capture_total_rna_fragments": total_count,
            "mrna_abundance": float(mrna_count),
            "nrna_abundance": float(nrna_count),
            "total_rna": float(total_count),
            "observed_mrna_fragments": mrna_count,
            "observed_nrna_fragments": nrna_count,
            "observed_total_rna_fragments": total_count,
        })

    return {
        "transcript_rows": transcript_rows,
        "origin_counts": {
            kind: int(origin_counts.get(kind, 0)) for kind in ("mrna", "nrna", "gdna")
        },
        "fragment_length_counts": length_counts,
        "fragment_lengths": {
            kind: _length_stats(counts) for kind, counts in length_counts.items()
        },
    }


def _write_post_capture_abundances(summary: dict[str, object], path: Path) -> None:
    rows = list(summary["transcript_rows"])  # type: ignore[index]
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "ref",
        "strand",
        "mrna_abundance",
        "nrna_abundance",
        "total_rna",
        "pre_capture_mrna_abundance",
        "pre_capture_nrna_abundance",
        "pre_capture_total_rna",
        "post_capture_mrna_fragments",
        "post_capture_nrna_fragments",
        "post_capture_total_rna_fragments",
        "observed_mrna_fragments",
        "observed_nrna_fragments",
        "observed_total_rna_fragments",
        "n_exons",
        "spliced_length",
        "genomic_span",
    ]
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in columns) + "\n")


def _write_fragment_length_truth(summary: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    length_counts = summary["fragment_length_counts"]  # type: ignore[index]
    with open(path, "w") as handle:
        handle.write("kind\tfragment_length\tcount\tfraction\n")
        for kind in ("mrna", "nrna", "rna", "gdna", "all"):
            counts = length_counts[kind]
            total = float(sum(counts.values()))
            for length in sorted(counts):
                count = int(counts[length])
                fraction = count / total if total > 0.0 else 0.0
                handle.write(f"{kind}\t{int(length)}\t{count}\t{fraction:.12g}\n")


def write_post_capture_truth(
    transcripts,
    abundance_path: Path,
    fl_path: Path,
    summary_path: Path,
    *,
    bam_path: Path | None = None,
    fastq_path: Path | None = None,
    condition: str | None = None,
    molecular_truth: str | None = None,
    gdna_strand_overdispersion: float | None = None,
) -> dict[str, object]:
    """Write condition-local empirical truth tables after simulation.

    ``abundance_path`` is compatible with existing analysis code: the
    ``mrna_abundance`` and ``nrna_abundance`` columns contain observed
    post-capture fragment counts. Molecular pre-capture abundances are retained
    in separate columns for provenance.
    """
    summary = summarize_post_capture_truth(
        transcripts,
        bam_path=bam_path,
        fastq_path=fastq_path,
    )
    _write_post_capture_abundances(summary, abundance_path)
    _write_fragment_length_truth(summary, fl_path)

    json_summary = {
        "condition": condition,
        "truth_kind": "post_capture_empirical",
        "gdna_strand_overdispersion_true": gdna_strand_overdispersion,
        "pre_capture_abundances": molecular_truth,
        "post_capture_abundances": str(abundance_path),
        "post_capture_fragment_lengths": str(fl_path),
        "origin_counts": summary["origin_counts"],
        "fragment_lengths": summary["fragment_lengths"],
        "files": {
            "pre_capture_abundances": molecular_truth,
            "post_capture_abundances": str(abundance_path),
            "post_capture_fragment_lengths": str(fl_path),
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as handle:
        json.dump(json_summary, handle, indent=2)
    return json_summary