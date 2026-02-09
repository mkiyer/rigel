"""
hulkrna.transcript — Transcript dataclass and GTF-to-Transcript parsing.

The Transcript is the atomic unit of transcription. Each transcript belongs
to a gene (identified by g_id / g_index) and contains an ordered list of
Exon intervals in 0-based half-open coordinates.
"""

import collections
import logging
from dataclasses import dataclass, field
from typing import Generator, List, Optional

import pysam

from .core import Interval, Strand
from .gtf import GTF


def _first_attr(val):
    """Return the first element when attrs are lists; otherwise the value itself."""
    if isinstance(val, list):
        return val[0] if val else None
    return val


@dataclass(slots=True)
class Transcript:
    ref: Optional[str] = None
    strand: Strand = Strand.NONE
    exons: List[Interval] = field(default_factory=list)
    length: Optional[int] = None
    t_id: Optional[str] = None
    g_id: Optional[str] = None
    g_name: Optional[str] = None
    g_type: Optional[str] = None
    t_index: int = -1
    g_index: int = -1
    is_basic: bool = False
    is_mane: bool = False
    is_ccds: bool = False
    seq: Optional[str] = None
    abundance: Optional[float] = None

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

    def to_gtf_lines(self, source: str = "transcriptor") -> Generator[str, None, None]:
        """
        Yields GTF-formatted string lines for the exons in this transcript.
        Coordinates are converted to 1-based inclusive.
        """
        # format score
        score_str = '.' if self.abundance is None else str(self.abundance)
        # collect attributes
        attrs = []
        if self.g_id: attrs.append(f'gene_id "{self.g_id}";')
        if self.t_id: attrs.append(f'transcript_id "{self.t_id}";')
        if self.g_name: attrs.append(f'gene_name "{self.g_name}";')
        if self.g_type: attrs.append(f'gene_type "{self.g_type}";')
        # add tags
        if self.is_basic: attrs.append('tag "basic";')
        if self.is_mane: attrs.append('tag "MANE_Select";')
        if self.is_ccds: attrs.append('tag "CCDS";')
        # build attributes string
        attr_str = " ".join(attrs)
        # yield exon lines
        for e in self.exons:
            fields = [
                str(self.ref) if self.ref else ".",
                source,
                "exon",
                str(e.start + 1),  # Convert 0-based to 1-based
                str(e.end),
                score_str, # score
                self.strand.to_str(), # strand
                ".", # phase
                attr_str
            ]
            yield "\t".join(fields)

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

    def fetch_seq(self, fasta_fh: pysam.FastaFile):
        """
        fetch the transcript sequence using a pysam FastaFile object
        regardless of transcript strand, the returned sequence will be on 
        the forward strand. this is so that exon intervals correspond 
        directly to the sequence coordinates
        """
        exon_seqs = []
        # length = 0
        for e in self.exons:
            s = fasta_fh.fetch(self.ref, e.start, e.end)
            exon_seqs.append(s)
            # length += (e.end - e.start)
        self.seq = ''.join(exon_seqs)
        # self.length = length

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
            'seq': self.seq,
            'abundance': self.abundance
        }
    
    @staticmethod
    def read_gtf(gtf_file: str) -> List['Transcript']:
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
