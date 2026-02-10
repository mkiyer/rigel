import multiprocessing as mp
import logging
import re
from enum import IntEnum

import numpy as np
import pysam
import zarr
from numcodecs import blosc, Blosc


class PileupBuffer:    
    DEFAULT_SIZE = 10000000
    DEFAULT_DTYPE = 'u4'

    def __init__(self, start=0, end=0, size=0, nchannels=1, dtype=None):
        self.start = start
        self.end = max(start, end)
        buflen = size if size > 0 else PileupBuffer.DEFAULT_SIZE
        buflen = max(buflen, end - start)
        dtype = dtype if dtype is not None else PileupBuffer.DEFAULT_DTYPE
        nchannels = max(nchannels, 1)
        self.buf = np.zeros((buflen, nchannels), dtype=dtype)

    def isempty(self):
        return self.start == self.end

    def get(self):
        j = self.end - self.start
        return (self.start, self.end, self.buf[:j])
    
    def reset(self, start=0):
        self.buf[:] = 0
        self.start = start
        self.end = start

    def advance(self, start, end):
        if start > end:
            raise ValueError(f'start={start} must be <= end={end}')
        if start < self.start:
            raise ValueError(f'start={start} must be >= current buffer start={self.start}. is the input sorted?')
        # check if not enough room within buffer
        buflen = self.buf.shape[0]
        if (end - self.start) > buflen:
            # flush buffer if buffer not empty and start position advances
            if (self.end > self.start) and (start > self.start):
                flushpos = min(start, self.end)
                i = flushpos - self.start           
                yield (self.start, flushpos, self.buf[:i])
                # reset buffer
                self.buf[:-i] = self.buf[i:]
                self.buf[-i:] = 0
            # check if buffer large enough
            if (end - start) > buflen:
                # resize buffer
                newlen = end - start
                nchannels = self.buf.shape[1]
                newbuf = np.zeros((newlen, nchannels), dtype=self.buf.dtype)
                newbuf[:buflen] = self.buf
                self.buf = newbuf
            # reset start/end
            self.start = start
            # dont advance beyond current end where data has been written
            self.end = max(self.end, start)

    def add(self, start, end, value=1, channel=0):
        if start > end:
            raise ValueError(f'start={start} must be <= end={end}')
        bufend = self.start + self.buf.shape[0]
        if start < self.start or end > bufend:
            raise ValueError(f'start={start} end={end} must lie within current buffer [{self.start}:{bufend}]')
        i = start - self.start
        j = end - self.start
        self.buf[i:j, channel] += value
        self.end = max(self.end, end)



STRAND_NONE = 0
STRAND_POS = 1
STRAND_NEG = 2
STRAND_INT_TO_STR = ['.', '+', '-']
STRAND_STR_TO_INT = {'.': 0, '+': 1, '-': 2}

_interval_string_re = re.compile(r'^(?P<ref>\w+)(?:\[(?P<strand>[+-.])\])?(?::(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?)?')

def _parse_interval_string(interval_string):
    """ref[strand]:start-end

    `strand` can be '+', '-', or '.'
    """
    m = _interval_string_re.match(interval_string)
    if m is None:
        ref, start, end, strand = None, None, None, STRAND_NONE
        # if interval string is just a single character,
        # try to extract something useful
        if len(interval_string) == 1:
            strand = STRAND_STR_TO_INT[interval_string[0]]
    else:
        ref = m.group('ref')
        start = m.group('start')
        end = m.group('end')
        strand = m.group('strand')
        if start is not None:
            start = int(start.replace(",", ""))
        if end is not None:
            end = int(end.replace(",", ""))
        if strand is None:
            strand = STRAND_NONE
        else:
            strand = STRAND_STR_TO_INT[strand]
    return ref, start, end, strand


def _parse_interval_tuple(interval):
    if interval is None:
        return None, None, None, STRAND_NONE
    if len(interval) == 1:
        return interval[0], None, None, STRAND_NONE
    elif len(interval) == 2:
        return interval[0], interval[1], interval[1] + 1, STRAND_NONE
    elif len(interval) == 3:
        return interval[0], interval[1], interval[2], STRAND_NONE
    else:
        if isinstance(interval[3], str):
            strand = STRAND_STR_TO_INT(interval[3])
        else:
            strand = interval[3]
        return interval[0], interval[1], interval[2], strand


def parse_interval(interval):
    """parses a genomic interval specifier in either the string
    format `<ref[strand]:start-end>` or tuple (ref, start, end, strand)

    :param interval: interval specifier
    :type interval: str or tuple

    >>> parse_interval("chr1:100-200")
    ("chr1", 100, 200, ".")
    >>> parse_interval("chr1:1,000-20,000")
    ("chr1", 1000, 20000)
    >>> parse_interval("chr1[-]:1,000-20,000")
    ("chr1", 1000, 20000, "-")
    >>> parse_interval("chrX")
    ("chrX", None, None, ".")
    >>> parse_interval("chrX:5")
    ("chrX", 5, 6, ".")
    >>> parse_interval(("chr2", 100, 300))
    ("chr2", 100, 300, ".")
    """
    if isinstance(interval, str):
        interval = _parse_interval_string(interval)
    interval = _parse_interval_tuple(interval)
    return interval


class GenomeTrack:
    VERSION = 'v0.0.1a'

    ATTRS = ['num_unique_mapped_r1', 
             'num_unique_mapped_r2', 
             'num_spliced_r1',
             'num_spliced_r2',
             'num_bp',
             'num_spliced_bp']
    NCHANNELS = 4

    STRANDED_NO = 0
    STRANDED_YES = 1
    STRANDED_REV = 2
    STRANDED_STR_TO_INT = {'no': STRANDED_NO, 'yes': STRANDED_YES, 'reverse': STRANDED_REV}

    # stranded, readnum, is_reverse
    UNSPLICED_CHANNEL_DICT = {
        (STRANDED_NO, 1, False): 0,
        (STRANDED_NO, 1, True): 1,
        (STRANDED_NO, 2, False): 0,
        (STRANDED_NO, 2, True): 1,
        (STRANDED_YES, 1, False): 0,
        (STRANDED_YES, 1, True): 1,
        (STRANDED_YES, 2, False): 1,
        (STRANDED_YES, 2, True): 0,
        (STRANDED_REV, 1, False): 1,
        (STRANDED_REV, 1, True): 0,
        (STRANDED_REV, 2, False): 0,
        (STRANDED_REV, 2, True): 1,
    }

    # storage options
    CHUNK_SIZE = 1<<20
    DTYPE = 'u4'

    # blosc compressor options
    BLOSC_CNAME = 'zstd'
    BLOSC_CLEVEL = 5
    BLOSC_SHUFFLE = Blosc.BITSHUFFLE
    BLOSC_COMPRESSOR = Blosc(cname=BLOSC_CNAME, 
                             clevel=BLOSC_CLEVEL, 
                             shuffle=BLOSC_SHUFFLE)
    # dimension separator
    DIM_SEP = '/'

    def __init__(self, filename):
        self._z = zarr.open(filename, mode='r')

    def get_refs(self):
        return ((x, self._z[x].shape[0]) for x in self._z)
    
    def fetch(self, interval):
        ref, start, end, strand = parse_interval(interval)
        if ref not in self._z:
            raise KeyError(f'ref={ref} not found')
        length = self._z[ref].shape[0]
        if start is not None and end is not None:
            if start < 0 or start >= end or end > length:
                raise ValueError(f'invalid query ref={ref} start={start} end={end} ({ref} length={length})')
            return self._z[ref][start:end]
        return self._z[ref][:]
    

def _get_unspliced_channel(stranded_code, readnum, is_reverse):
    channel = GenomeTrack.UNSPLICED_CHANNEL_DICT[(stranded_code, readnum, is_reverse)]
    return channel + 2    


def bam_se_coverage_worker(args):
    # unpack args
    bam_file, tmpdir, ref, length, stranded = args

    # TODO: zarr documentation suggests we turn off blosc threads when using 
    # multiprocessing?
    # blosc.use_threads = False

    # tmp storage
    store_path = f'{tmpdir}/{ref}.zarr'
    store = zarr.DirectoryStore(store_path, dimension_separator=GenomeTrack.DIM_SEP)

    stranded_code = GenomeTrack.STRANDED_STR_TO_INT[stranded]
    nchannels = GenomeTrack.NCHANNELS
    chunk_length = min(GenomeTrack.CHUNK_SIZE, length)
    z = zarr.zeros(path=ref,
                   store=store, 
                   overwrite=True,
                   shape=(length, nchannels),
                   chunks=(chunk_length, nchannels),
                   dtype=GenomeTrack.DTYPE,
                   write_empty_chunks=False,
                   compressor=GenomeTrack.BLOSC_COMPRESSOR)
    # pileup buffer
    buf = PileupBuffer(nchannels=nchannels)
    # open bam file
    bamfh = pysam.AlignmentFile(bam_file)
    # iterate through bam file
    num_unique_mapped_se = 0
    num_unique_mapped_r1 = 0
    num_unique_mapped_r2 = 0
    num_bp = 0
    num_spliced_se = 0
    num_spliced_r1 = 0
    num_spliced_r2 = 0
    num_spliced_bp = 0
    logging.debug(f'reference={ref} length={length}')
    for read in bamfh.fetch(ref):
        # filter reads
        if read.is_unmapped or read.is_qcfail or read.is_secondary or read.is_supplementary:
            continue
        # only use unique mapping reads
        if read.get_tag('NH') > 1:
            continue

        for start, end, arr in buf.advance(read.reference_start, read.reference_end):
            logging.debug(f'[worker] flushing ref={ref} start={start} end={end}')
            z[start:end] = arr

        spliced = False
        if read.has_tag('XS'):
            spliced = True
            strand = read.get_tag('XS')
            channel = 1 if strand == '-' else 0

        if not read.is_paired:
            num_unique_mapped_se += 1
            if spliced:
                num_spliced_se += 1
            else:
                channel = _get_unspliced_channel(stranded_code, 1, read.is_reverse)
        elif read.is_read1:
            num_unique_mapped_r1 += 1
            if spliced: 
                num_spliced_r1 += 1
            else:  
                channel = _get_unspliced_channel(stranded_code, 1, read.is_reverse)
        elif read.is_read2:
            num_unique_mapped_r2 += 1
            if spliced: 
                num_spliced_r2 += 1
            else:
                channel = _get_unspliced_channel(stranded_code, 1, read.is_reverse)
        else:
            pass

        for start, end in read.get_blocks():
            bp = end - start
            num_bp += bp
            if spliced:
                num_spliced_bp += bp
            buf.add(start, end, value=1, channel=channel)
    bamfh.close()
    # flush and write whatever is left in buffer
    start, end, arr = buf.get()
    if end > start:
        logging.debug(f'[worker] flushing ref={ref} start={start} end={end}')
        z[start:end] = arr
    logging.info(f'[worker] completed ref={ref}')
    z.attrs['num_unique_mapped_se'] = num_unique_mapped_se
    z.attrs['num_unique_mapped_r1'] = num_unique_mapped_r1
    z.attrs['num_unique_mapped_r2'] = num_unique_mapped_r2
    z.attrs['num_spliced_se'] = num_spliced_se
    z.attrs['num_spliced_r1'] = num_spliced_r1
    z.attrs['num_spliced_r2'] = num_spliced_r2
    z.attrs['num_bp'] = num_bp
    z.attrs['num_spliced_bp'] = num_spliced_bp
    store.close()
    return((ref, store_path))


def bam_se_coverage_mp(bam_file, output_file, tmpdir, stranded=None, num_processes=1):
    # transcript strand for strand-specific data
    if stranded is None:
        stranded = 'no'
    stranded_code = GenomeTrack.STRANDED_STR_TO_INT[stranded]

    # read references
    bamfh = pysam.AlignmentFile(bam_file)
    refs = list(zip(map(str, bamfh.references), map(int, bamfh.lengths)))
    bamfh.close()

    # reduce store
    tmpstore = zarr.TempStore(suffix='out', dir=tmpdir, 
                              dimension_separator=GenomeTrack.DIM_SEP)
    root = zarr.group(store=tmpstore, overwrite=True)
    root.attrs['version'] = GenomeTrack.VERSION
    root.attrs['stranded'] = stranded
    root.attrs['stranded_code'] = stranded_code

    # create a parallel process pool
    with mp.Pool(num_processes) as pool:
        logging.info(f'created process pool with n={num_processes} workers')
        job_args = []
        for ref, length in refs:
            job_args.append((bam_file, tmpdir, ref, length, stranded))
        # use imap_unordered to parallelize the function with multiple arguments
        results = pool.imap_unordered(bam_se_coverage_worker, job_args)
        for ref, source_store_path in results:
            logging.info(f'ref={ref} copying to store')
            source_store = zarr.DirectoryStore(source_store_path, dimension_separator=GenomeTrack.DIM_SEP)
            #zarr.copy_store(source_store, tmpstore, source_path=ref, dest_path=ref, log=sys.stdout)
            zarr.copy_store(source_store, tmpstore, source_path=ref, dest_path=ref)
            logging.info(f'ref={ref} finished')

    # aggregate attributes
    logging.info('[root] aggregating attributes')
    for k in GenomeTrack.ATTRS:
        v = sum(root[ref].attrs[k] for ref in root)
        root.attrs[k] = v
        logging.debug(f'attr={k} value={v}')
    logging.info('[root] done aggregating attributes')

    # output as a zipstore
    logging.info(f'[root] writing output {output_file} as zipstore')
    zipstore = zarr.ZipStore(output_file, mode='w', dimension_separator=GenomeTrack.DIM_SEP)
    zarr.copy_store(tmpstore, zipstore)
    zipstore.close()
    logging.info(f'[root] done writing {output_file}')