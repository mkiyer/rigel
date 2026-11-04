"""Simulator read names: each fragment's ground-truth origin, built by :func:`format_rna_name` and
:func:`format_gdna_name` and read back by :func:`parse_origin`. The three forms are::

    mature RNA    {t_id}:{start}-{end}:{strand}:{index}
    nascent RNA   nrna_{t_id}:{start}-{end}:{strand}:{index}
    gDNA          gdna:{ref}:{start}-{end}:{strand}:{index}

``start``-``end`` is the fragment's 0-based half-open interval: on the transcript's own 5'-to-3' sequence for RNA
(a nascent fragment's template is a single-exon nRNA entity), on reference ``ref`` for gDNA. ``strand`` is ``f`` or
``r``: the transcript's strand for RNA (``r`` for minus, ``f`` otherwise), the fragment's genomic strand for gDNA.
``index`` counts the writer's fragments: within an RNA row, and across one writer's gDNA chunks. A FASTQ header's
leading ``@`` and ``/1`` or ``/2`` suffix are stripped by the parser.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

OriginKind = Literal["mrna", "nrna", "gdna"]

#: The strand alphabet: ``f`` forward (+), ``r`` reverse (−).
STRAND_FORWARD = "f"
STRAND_REVERSE = "r"
STRANDS = (STRAND_FORWARD, STRAND_REVERSE)

#: Glued to the transcript id of a nascent-RNA name; a mature-RNA name has no prefix.
NASCENT_PREFIX = "nrna_"
#: A gDNA name's leading field of its own, followed by ``:{ref}``.
GDNA_PREFIX = "gdna"

__all__ = [
    "GDNA_PREFIX",
    "NASCENT_PREFIX",
    "STRANDS",
    "STRAND_FORWARD",
    "STRAND_REVERSE",
    "Origin",
    "OriginKind",
    "format_gdna_name",
    "format_rna_name",
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


def format_rna_name(
    t_id: str, start: int, end: int, strand: str, index: int, *, nascent: bool
) -> str:
    """The read name of an RNA fragment over ``[start, end)`` of transcript ``t_id``'s sequence; ``nascent`` marks a
    fragment of an nRNA entity's nascent row."""
    lead = f"{NASCENT_PREFIX}{t_id}" if nascent else t_id
    return f"{lead}:{start}-{end}:{strand}:{index}"


def format_gdna_name(ref: str, start: int, end: int, strand: str, index: int) -> str:
    """The read name of a gDNA fragment over ``[start, end)`` of reference ``ref``."""
    return f"{GDNA_PREFIX}:{ref}:{start}-{end}:{strand}:{index}"


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
    """Parse a simulator FASTQ/BAM read name into a structured origin.

    The name is split from the right, since only the leading field (a transcript id, or ``gdna:{ref}``) may
    contain ``:``. A strand outside ``f`` / ``r`` raises ``ValueError``: such a name is not the simulator's.
    """
    qname = _normalize_qname(qname)
    parts = qname.rsplit(":", 3)
    if len(parts) != 4 or not parts[0]:
        raise ValueError(f"invalid simulator read name: {qname!r}")
    name, interval_text, strand, index_text = parts
    if strand not in STRANDS:
        raise ValueError(
            f"invalid strand {strand!r} in simulator read name {qname!r}: the alphabet is "
            f"{STRAND_FORWARD!r} / {STRAND_REVERSE!r}"
        )
    start, end = _parse_interval(interval_text)
    index = _parse_index(index_text)

    prefix, _, ref = name.partition(":")
    if prefix == GDNA_PREFIX:
        if not ref:
            raise ValueError(f"gDNA read name without a reference: {qname!r}")
        return Origin(
            kind="gdna",
            transcript_id=None,
            ref=ref,
            start=start,
            end=end,
            strand=strand,
            index=index,
        )
    if name.startswith(NASCENT_PREFIX):
        return Origin(
            kind="nrna",
            transcript_id=name.removeprefix(NASCENT_PREFIX),
            ref=None,
            start=start,
            end=end,
            strand=strand,
            index=index,
        )
    return Origin(
        kind="mrna",
        transcript_id=name,
        ref=None,
        start=start,
        end=end,
        strand=strand,
        index=index,
    )
