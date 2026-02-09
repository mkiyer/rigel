import collections
import logging
import gzip
import itertools
import bisect

# local imports
from hulkrna.gtf import GTF
    
# exon as namedtuple
Exon = collections.namedtuple('Exon', ['start', 'end'])

class Transcript:
    def __init__(self):
        self.ref = None
        self.strand = '.'
        self.exons = []
        self.cds = []
        self.t_id = None
        self.gene_id = None
        self.gene_name = None
        self.gene_type = None
        self.is_ccds = False
        self.is_mane = False
        self.is_basic = False
        self.length = 0

    @staticmethod
    def from_gtf(feature) -> 'Transcript':
        '''Create a Transcript object from a GTF feature'''
        self = Transcript()
        self.ref = feature.seqid
        self.strand = feature.strand
        self.t_id = feature.attrs['transcript_id']
        self.gene_id = feature.attrs['gene_id']
        self.gene_name = feature.attrs.get('gene_name', None)
        self.gene_type = feature.attrs.get('gene_type', None)
        # gencode specific tags
        self.is_basic = "basic" in feature.tags
        self.is_mane = ("MANE_Select" in feature.tags) or ("MANE_Plus_Clinical" in feature.tags)
        self.is_ccds = "CCDS" in feature.tags
        return self

    def __str__(self):
        return (f'{self.__class__.__name__}(t_id={self.t_id}, '
                f'gene_id={self.gene_id}, chrom={self.chrom}, '
                f'strand={self.strand}, exons={self.exons})')

    @property
    def start(self):
        return self.exons[0].start if self.exons else 0
    @property
    def end(self):
        return self.exons[-1].end if self.exons else 0

    def compute_length(self):
        return sum((e.end - e.start) for e in self.exons)

    def introns(self):
        for i in range(1, len(self.exons)):
            yield self.exons[i-1].end, self.exons[i].start

    def find_exon_index(self, pos, exon_start=0):
        return bisect.bisect_right(self.exon_cumsum, pos, lo=exon_start)


class Gene:
    def __init__(self):
        self.gene_id = None
        self.gene_name = None
        self.gene_type = None
        self.ref = None
        self.start = None
        self.end = None
        self.strand = '.'
        self.length = 0
        self.num_exons = 0
        self.num_introns = 0
        self.t_ids = set()

    @staticmethod
    def from_transcript(t: Transcript) -> 'Gene':
        '''Create a Gene object from a Transcript object'''
        self = Gene()
        self.gene_id = t.gene_id
        self.gene_name = t.gene_name
        self.gene_type = t.gene_type
        return self
    
    @staticmethod
    def from_gtf(feature) -> 'Gene':
        '''Create a Gene object from a GTF feature'''
        self = Gene()
        self.gene_id = feature.attrs['gene_id']
        self.gene_name = feature.attrs.get('gene_name', None)
        self.gene_type = feature.attrs.get('gene_type', None)
        return self
    
    def __str__(self):
        return (f'{self.__class__.__name__}(gene_id={self.gene_id}, '
                f'gene_name={self.gene_name}, gene_type={self.gene_type}, '
                f't_ids={self.t_ids}')



# open gtf file based on extension
def _gtf_open(filename):
    mode = 'rt'
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
    
# improve readability by moving parsing code out of main loop
def _gtf_parse_features(filename, feature_types=None):
    if feature_types is None:
        feature_types = {'exon', 'CDS'}
    with _gtf_open(filename) as f:
        for feature in GTF.parse(f):
            if feature.feature in feature_types:
                yield feature


def cluster_intervals(intervals):
    # sort based on start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    # initialize the first cluster
    clusters = []
    cluster_start, cluster_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start >= cluster_end:
            # finalize cluster
            clusters.append((cluster_start, cluster_end))
            # start new cluster
            cluster_start, cluster_end = (start, end)
        else:
            cluster_end = max(cluster_end, end)
    # finalize final cluster
    clusters.append((cluster_start, cluster_end))
    return clusters


class GencodeDB:
    def __init__(self):
        self.transcripts = collections.OrderedDict()
        self.genes = collections.OrderedDict()

    @staticmethod
    def from_gencode_gtf(gtf_file):
        self = GencodeDB()
        logging.info(f'Reading GTF file: {gtf_file}')
        num_lines = 0
        # read gtf exons into a transcript dictionary
        for f in _gtf_parse_features(gtf_file, {'exon', 'CDS'}):
            t_id = f.attrs['transcript_id']
            # create transcript if not seen before
            if t_id not in self.transcripts:
                t = Transcript.from_gtf(f)
                self.transcripts[t_id] = t
            else:
                t = self.transcripts[t_id]  
            # update transcript
            if f.feature == 'exon':
                t.exons.append(Exon(f.start, f.end))
            elif f.feature == 'CDS':
                t.cds.extend((f.start, f.end))
            # log
            num_lines += 1
            if num_lines % 100000 == 0:
                logging.debug(f'Read {num_lines} features')
        logging.info(f'[DONE] Read {num_lines} features')
        self.compute_transcript_attributes()
        self.compute_gene_attributes()
        return self

    def compute_transcript_attributes(self):
        # compute transcript attributes
        logging.info(f'[START] Computing transcript attributes')
        for t in self.transcripts.values():
            t.length = t.compute_length()
            t.exons.sort()
            t.exon_sizes = [(e.end - e.start) for e in t.exons]
            t.exon_cumsum = list(itertools.accumulate(t.exon_sizes))
            t.exon_starts = [(t.exon_cumsum[i] - t.exon_sizes[i]) for i in range(len(t.exons))]
            if len(t.cds) == 0:
                t.cds = Exon(t.start, t.end)
            else:
                t.cds.sort()
                t.cds = Exon(t.cds[0], t.cds[-1])
            cds_start_index = bisect.bisect_right(t.exons, t.cds[0], key=lambda e: e.end)
            cds_end_index = bisect.bisect_left(t.exons, t.cds[1], key=lambda e: e.end)
            t.cds_start = t.exon_starts[cds_start_index] + (t.cds[0] - t.exons[cds_start_index].start)
            t.cds_end = t.exon_starts[cds_end_index] + (t.cds[1] - t.exons[cds_end_index].start)
            assert t.cds_start < t.cds_end
        logging.info(f'[DONE] Computing transcript attributes')

    def compute_gene_attributes(self):
        # erase existing attributes
        self.genes = collections.OrderedDict()
        # create genes from transcripts
        for t in self.transcripts.values():
            t_id = t.t_id
            gene_id = t.gene_id
            # create gene if not seen before
            if gene_id not in self.genes:
                g = Gene.from_transcript(t)
                self.genes[gene_id] = g
            else:
                g = self.genes[gene_id]
            g.t_ids.add(t_id)
        # compute gene attributes
        logging.info(f'[START] Computing gene attributes')
        for gene_id in self.genes.keys():
            g = self.genes[gene_id]
            refs = set()
            strands = set()
            exons = set()
            introns = set()
            for t_id in g.t_ids:
                t = self.transcripts[t_id]
                refs.add(t.ref)
                strands.add(t.strand)
                exons.update(t.exons)
                introns.update(t.introns())
            assert len(refs) == 1
            assert len(strands) == 1
            # calculate gene attributes
            exons = cluster_intervals(exons)
            g.ref = refs.pop()
            g.strand = strands.pop()
            g.start = exons[0][0]
            g.end = exons[-1][1]
            g.length = sum((e[1] - e[0]) for e in exons)
            g.num_exons = len(exons)
            g.num_introns = len(introns)        
        logging.info(f'[DONE] Computing gene attributes')