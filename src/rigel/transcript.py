"""
rigel.transcript — Transcript dataclass and GTF-to-Transcript parsing.

The Transcript is the atomic unit of transcription. Each transcript belongs
to a gene (identified by g_id / g_index) and contains an ordered list of
Exon intervals in 0-based half-open coordinates.
"""

import logging
from dataclasses import dataclass, field
from typing import Literal

from .config import CONSTANTS
from .types import Interval, Strand
from .gtf import GTFRecord


def _malformed_exon(parse_mode: str, message: str) -> None:
    """An exon record the annotation cannot place: strict raises, warn-and-skip logs and drops that exon."""
    if parse_mode == "strict":
        raise ValueError(f"Malformed GTF: {message}")
    logging.warning("Skipping GTF %s", message)


def _disjoint_exons(t: "Transcript", parse_mode: str) -> list[Interval]:
    """A transcript's exons in order, disjoint. Exons with no base between them (RefSeq writes them at genome indels)
    are one block and merge; an exon that overlaps or repeats the one before it is malformed, so strict refuses it and
    warn-and-skip drops it."""
    exons: list[Interval] = []
    for e in sorted(t.exons):
        if exons and e.start < exons[-1].end:
            _malformed_exon(
                parse_mode,
                f"exon {t.ref}:{e.start + 1}-{e.end} of transcript {t.t_id} overlaps its exon "
                f"{exons[-1].start + 1}-{exons[-1].end}",
            )
            continue
        if exons and e.start == exons[-1].end:
            exons[-1] = Interval(exons[-1].start, e.end)
            continue
        exons.append(e)
    return exons


def _first_attr(val):
    """Return the first element when attrs are lists; otherwise the value itself."""
    if isinstance(val, list):
        return val[0] if val else None
    return val


@dataclass(slots=True)
class Transcript:
    ref: str | None = None
    strand: Strand = Strand.NONE
    exons: list[Interval] = field(default_factory=list)
    length: int | None = None
    t_id: str | None = None
    g_id: str | None = None
    g_name: str | None = None
    g_type: str | None = None
    t_index: int = -1
    g_index: int = -1
    is_basic: bool = False
    is_mane: bool = False
    is_ccds: bool = False
    is_nrna: bool = False
    is_synthetic: bool = False
    nrna_t_index: int = -1
    """t_index of the associated nRNA entity (-1 = none).

    For multi-exon transcripts: points to the nRNA transcript
    (synthetic or annotated single-exon) that covers this transcript's
    merged TSS/TES span.
    For single-exon transcripts: points to itself if it is the
    nascent-equiv covering a merged multi-exon span, else -1.
    For synthetics: -1.
    """
    nrna_n_contributors: int = 0
    """Number of annotated multi-exon transcripts merged into this nRNA entity.

    Only meaningful for synthetic nRNA transcripts and nascent equivalents.
    Zero for regular annotated transcripts.
    """
    abundance: float | None = None
    nrna_abundance: float = 0.0
    """Nascent RNA (pre-mRNA) molecular abundance.

    Independent of ``abundance`` (mature mRNA molecular abundance).
    Read fractions are proportional to::

        mRNA  : abundance      × spliced_effective_length
        nRNA  : nrna_abundance × genomic_span_effective_length

    Setting ``nrna_abundance = k × abundance`` means there are *k* times
    as many nascent RNA molecules as mature mRNA molecules.
    """
    aliases: tuple[str, ...] = ()
    """The IDs of the identical transcripts ``--collapse-duplicate-transcripts`` folded into this one."""

    @classmethod
    def from_gtf(cls, feature) -> "Transcript":
        """Create a Transcript object from a GTF feature"""
        t = cls()
        if not hasattr(feature, "seqname"):
            raise AttributeError("GTF feature missing seqname")
        t.ref = feature.seqname
        t.strand = Strand.from_str(getattr(feature, "strand", "."))
        t.abundance = getattr(feature, "score", None)

        attrs = getattr(feature, "attrs", None) or {}
        t.t_id = _first_attr(attrs.get("transcript_id"))
        t.g_id = _first_attr(attrs.get("gene_id"))
        t.g_name = _first_attr(attrs.get("gene_name"))
        # GENCODE writes gene_type, Ensembl gene_biotype
        t.g_type = _first_attr(attrs.get("gene_type")) or _first_attr(attrs.get("gene_biotype"))
        # gencode specific tags
        tags = getattr(feature, "tags", set()) or set()
        t.is_basic = "basic" in tags
        t.is_mane = ("MANE_Select" in tags) or ("MANE_Plus_Clinical" in tags)
        t.is_ccds = "CCDS" in tags
        return t

    @property
    def start(self) -> int:
        """Start coordinate of the first exon (0 if empty)."""
        return self.exons[0].start if self.exons else 0

    @property
    def end(self) -> int:
        """End coordinate of the last exon (0 if empty)."""
        return self.exons[-1].end if self.exons else 0

    def compute_length(self):
        if self.length is not None:
            return self.length
        self.length = sum((e.end - e.start) for e in self.exons)
        return self.length

    def introns(self):
        for i in range(1, len(self.exons)):
            yield self.exons[i - 1].end, self.exons[i].start

    def to_dict(self):
        return {
            "ref": self.ref,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "length": self.length,
            "t_id": self.t_id,
            "g_id": self.g_id,
            "t_index": self.t_index,
            "g_index": self.g_index,
            "g_name": self.g_name,
            "g_type": self.g_type,
            "is_basic": self.is_basic,
            "is_mane": self.is_mane,
            "is_ccds": self.is_ccds,
            "is_nrna": self.is_nrna,
            "is_synthetic": self.is_synthetic,
            "nrna_t_index": self.nrna_t_index,
            "nrna_n_contributors": self.nrna_n_contributors,
        }

    @staticmethod
    def read_gtf(
        gtf_file: str,
        *,
        parse_mode: Literal["strict", "warn-skip"] = "strict",
    ) -> list["Transcript"]:
        """Read GTF and construct list of Transcript objects."""
        transcripts: dict[str, Transcript] = {}
        logging.debug("[Transcript] Reading GTF file: %s", gtf_file)
        n_exon_features = 0
        for f in GTFRecord.parse_file(gtf_file, parse_mode=parse_mode):
            if f.feature != "exon":
                continue
            t_id = _first_attr(f.attrs.get("transcript_id"))
            where = f"exon {f.seqname}:{f.start + 1}-{f.end}({f.strand})"
            if t_id is None:
                _malformed_exon(parse_mode, f"{where} has no transcript_id")
                continue
            if t_id not in transcripts:
                t = Transcript.from_gtf(f)
                transcripts[t_id] = t
            else:
                t = transcripts[t_id]
                # a transcript's exons share its reference and strand; a stray one would join two chromosomes
                if f.seqname != t.ref or Strand.from_str(f.strand) != t.strand:
                    _malformed_exon(
                        parse_mode,
                        f"{where} is off transcript {t_id}'s {t.ref}({t.strand.to_str()})",
                    )
                    continue
                # and its gene: two genes sharing a transcript ID would merge behind a fake intron
                g_id = _first_attr(f.attrs.get("gene_id"))
                if g_id != t.g_id:
                    _malformed_exon(
                        parse_mode,
                        f"{where} is in gene {g_id}, but transcript {t_id} is in {t.g_id}",
                    )
                    continue
            t.exons.append(Interval(f.start, f.end))
            n_exon_features += 1
            if n_exon_features % CONSTANTS.resources.gtf_log_interval == 0:
                logging.debug("[Transcript] Read %d GTF exon features", n_exon_features)
        logging.debug("[Transcript] Done reading GTF: %d exon features", n_exon_features)

        logging.debug("[Transcript] Processing transcripts")
        for t in transcripts.values():
            t.exons = _disjoint_exons(t, parse_mode)
            t.compute_length()
            if len(t.exons) == 1:
                t.is_nrna = True
        sorted_transcripts = sorted(
            transcripts.values(),
            key=lambda t: (t.ref, t.start, t.end, t.strand),
        )
        logging.debug("[Transcript] Read %d transcripts", len(sorted_transcripts))
        return sorted_transcripts
