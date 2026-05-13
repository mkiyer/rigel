"""Read-name truth parsing and fragment-origin aggregation."""

from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import pysam

OriginKind = Literal["mrna", "nrna", "gdna"]

__all__ = [
    "Origin",
    "OriginKind",
    "count_mrna_by_transcript_from_bam",
    "count_mrna_by_transcript_from_fastq",
    "count_origins_from_bam",
    "count_origins_from_fastq",
    "deduplicate_bam_qnames",
    "parse_origin",
]


@dataclass(frozen=True)
class Origin:
    """Ground-truth origin encoded in a simulator read name."""

    kind: OriginKind
    transcript_id: str | None
    ref: str | None
    start: int | None
    end: int | None
    strand: str | None
    index: int | None


def _normalize_qname(qname: str) -> str:
    qname = qname.strip()
    if qname.startswith("@"):
        qname = qname[1:]
    qname = qname.split()[0]
    if qname.endswith("/1") or qname.endswith("/2"):
        qname = qname[:-2]
    return qname


def _parse_interval(text: str) -> tuple[int, int]:
    start_text, end_text = text.split("-", 1)
    return int(start_text), int(end_text)


def _parse_index(text: str) -> int | None:
    try:
        return int(text)
    except ValueError:
        return None


def parse_origin(qname: str) -> Origin:
    """Parse a simulator FASTQ/BAM read name into a structured origin."""
    qname = _normalize_qname(qname)
    parts = qname.split(":")

    if not parts or not parts[0]:
        raise ValueError("empty simulator read name")

    if parts[0] == "gdna":
        if len(parts) == 4:
            ref = None
            interval_text, strand, index_text = parts[1:]
        elif len(parts) == 5:
            ref = parts[1]
            interval_text, strand, index_text = parts[2:]
        else:
            raise ValueError(f"invalid gDNA read name: {qname!r}")
        start, end = _parse_interval(interval_text)
        return Origin(
            kind="gdna",
            transcript_id=None,
            ref=ref,
            start=start,
            end=end,
            strand=strand,
            index=_parse_index(index_text),
        )

    if len(parts) != 4:
        raise ValueError(f"invalid RNA read name: {qname!r}")

    name, interval_text, strand, index_text = parts
    start, end = _parse_interval(interval_text)
    if name.startswith("nrna_"):
        return Origin(
            kind="nrna",
            transcript_id=name.removeprefix("nrna_"),
            ref=None,
            start=start,
            end=end,
            strand=strand,
            index=_parse_index(index_text),
        )
    return Origin(
        kind="mrna",
        transcript_id=name,
        ref=None,
        start=start,
        end=end,
        strand=strand,
        index=_parse_index(index_text),
    )


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