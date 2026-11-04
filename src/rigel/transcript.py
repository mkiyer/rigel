"""rigel.transcript — the Transcript dataclass, read from a GTF.

A transcript belongs to one gene (``g_id`` / ``g_index``) and holds its exons in order, 0-based half-open.
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
    is_single_exon: bool = False
    """The transcript has one exon, as every synthetic span does."""
    is_synthetic: bool = False
    span_t_index: int = -1
    """t_index of the transcript's span entity, -1 for none. A multi-exon transcript's is the entity over its merged
    TSS/TES span (a synthetic span, or the annotated single-exon transcript covering it); a single-exon transcript's
    is itself if it covers a merged span; a synthetic's is -1."""
    span_n_contributors: int = 0
    """The number of annotated multi-exon transcripts whose merged spans this entity stands for; zero except on
    synthetic spans and the single-exon transcripts that cover a span."""
    abundance: float | None = None
    nrna_abundance: float = 0.0
    """A simulator input, read only by ``rigel.sim``: the abundance of RNA that has not spliced, sampled contiguously
    over the genomic span (pooled onto the span's single-exon entity, ``sim.whole_genome.assign_nrna_to_entities``).
    Independent of ``abundance``; read fractions are proportional to::

        spliced template : abundance      × spliced_effective_length
        genomic span     : nrna_abundance × genomic_span_effective_length
    """
    aliases: tuple[str, ...] = ()
    """The IDs of the identical transcripts ``--collapse-duplicate-transcripts`` folded into this one."""

    @classmethod
    def from_gtf(cls, feature) -> "Transcript":
        """A Transcript from a GTF feature's reference, strand, IDs and tags."""
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
        # GENCODE tags
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
            "is_single_exon": self.is_single_exon,
            "is_synthetic": self.is_synthetic,
            "span_t_index": self.span_t_index,
            "span_n_contributors": self.span_n_contributors,
        }

    @staticmethod
    def read_gtf(
        gtf_file: str,
        *,
        parse_mode: Literal["strict", "warn-skip"] = "strict",
    ) -> list["Transcript"]:
        """Read a GTF's exon records into transcripts, sorted by ``(ref, start, end, strand)``."""
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
                t.is_single_exon = True
        sorted_transcripts = sorted(
            transcripts.values(),
            key=lambda t: (t.ref, t.start, t.end, t.strand),
        )
        logging.debug("[Transcript] Read %d transcripts", len(sorted_transcripts))
        return sorted_transcripts
