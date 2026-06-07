"""Simulator read-name encoding: the ground-truth origin parser.

The simulator encodes each fragment's ground truth in its read name
(``{t_id}:{start}-{end}:{strand}:{index}`` for RNA, ``gdna:[ref:]{start}-{end}:{strand}:{index}``
for gDNA). This module owns the *parsing* side — a small, dependency-free decoder consumed by the
benchmark/analysis tooling. Truth-table *writing* (and the counting helpers that aggregate parsed
origins) live in :mod:`truth`, which imports :func:`parse_origin` from here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

OriginKind = Literal["mrna", "nrna", "gdna"]

__all__ = ["Origin", "OriginKind", "parse_origin"]


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
