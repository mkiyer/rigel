"""
hulkrna.transcript — Transcript dataclass and GTF-to-Transcript parsing.

The Transcript is the atomic unit of transcription. Each transcript belongs
to a gene (identified by g_id / g_index) and contains an ordered list of
Exon intervals in 0-based half-open coordinates.
"""

import collections
import logging
from dataclasses import dataclass, field

from .types import Interval, Strand
from .gtf import GTF


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
    abundance: float | None = None

    @classmethod
    def from_gtf(cls, feature) -> 'Transcript':
        '''Create a Transcript object from a GTF feature'''
        t = cls()
        if not hasattr(feature, 'seqname'):
            raise AttributeError("GTF feature missing seqname")
        t.ref = feature.seqname
        t.strand = Strand.from_str(getattr(feature, 'strand', '.'))
        t.abundance = getattr(feature, 'score', None)

        attrs = getattr(feature, 'attrs', None) or {}
        t.t_id = _first_attr(attrs.get('transcript_id'))
        t.g_id = _first_attr(attrs.get('gene_id'))
        t.g_name = _first_attr(attrs.get('gene_name'))
        t.g_type = _first_attr(attrs.get('gene_type'))
        # gencode specific tags
        tags = getattr(feature, 'tags', set()) or set()
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
            yield self.exons[i-1].end, self.exons[i].start

    def to_dict(self):
        return {
            'ref': self.ref,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'length': self.length,
            't_id': self.t_id,
            'g_id': self.g_id,
            't_index': self.t_index,
            'g_index': self.g_index,
            'g_name': self.g_name,
            'g_type': self.g_type,
            'is_basic': self.is_basic,
            'is_mane': self.is_mane,
            'is_ccds': self.is_ccds,
            'abundance': self.abundance
        }
    
    @staticmethod
    def read_gtf(gtf_file: str) -> list['Transcript']:
        '''
        read GTF and construct list of Transcript objects
        '''
        transcripts = collections.OrderedDict()
        # read gtf exons into a transcript dictionary
        logging.debug(f'[Transcript] Reading GTF file: {gtf_file}')
        num_lines = 0
        for f in GTF.parse_file(gtf_file):
            if f.feature != 'exon':
                continue
            t_id = f.attrs['transcript_id']
            # lookup transcript, create if new
            if t_id not in transcripts:
                t = Transcript.from_gtf(f)
                transcripts[t_id] = t
            else:
                # get existing transcript
                t = transcripts[t_id]
            # update transcript exons
            t.exons.append(Interval(f.start, f.end))
            num_lines += 1
            if num_lines % 100000 == 0:
                logging.debug(f'[Transcript] Read {num_lines} GTF features')
        logging.debug(f'[Transcript] Done reading GTF and found {num_lines} features')

        # sort transcript exons
        logging.debug(f'[Transcript] Processing transcripts')
        for t in transcripts.values():
            t.exons.sort()
            t.compute_length()
        # sort transcripts by ref, start, end, strand
        transcripts = sorted(transcripts.values(), key=lambda t: (t.ref, t.start, t.end, t.strand))
        logging.debug(f'[Transcript] Read {len(transcripts)} transcripts')
        return transcripts
