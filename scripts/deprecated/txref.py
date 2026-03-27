import logging
import os
from enum import IntEnum
from typing import NamedTuple

import pandas as pd
import cgranges

from hulkrna.transcript import Strand


class SJ(NamedTuple):
    ref: str
    start: int
    end: int
    strand: int
    t_index: int
    g_index: int

class IntervalType(IntEnum):
    EXON = 0
    INTRON = 1
    INTERGENIC = 2

class Interval(NamedTuple):
    '''Genomic interval with metadata'''
    ref: str
    start: int
    end: int
    strand: int = Strand.NONE
    interval_type: int = IntervalType.INTERGENIC
    t_index: int = -1
    g_index: int = -1


def _merge_sets_with_relaxation(t_sets, g_sets):
    '''
    Merge sets of transcript and gene indices with a progressive 
    relaxation of constraints

    Progressive relaxation:
    1. Intersection of all sets (most specific)
    2. Intersection of non-empty sets (less specific)
    3. Union of all sets (most sensitive)    
    '''   
    # 1. cumulative intersection of sets (most specific)
    t_inds = frozenset.intersection(*t_sets) if t_sets else frozenset()
    g_inds = frozenset.intersection(*g_sets) if g_sets else frozenset()
    if t_inds or g_inds:
        return t_inds, g_inds

    # 2. intersection of non-empty sets (skip empty sets in intersection)
    t_sets_nonempty = [x for x in t_sets if x]
    g_sets_nonempty = [x for x in g_sets if x]
    t_inds = frozenset.intersection(*t_sets_nonempty) if t_sets_nonempty else frozenset()
    g_inds = frozenset.intersection(*g_sets_nonempty) if g_sets_nonempty else frozenset()
    if t_inds or g_inds:
        return t_inds, g_inds

    # 3. union of all sets (most sensitive)
    t_inds = frozenset.union(*t_sets) if t_sets else frozenset()
    g_inds = frozenset.union(*g_sets) if g_sets else frozenset()
    return t_inds, g_inds


class TxIndex:
    TRANSCRIPTS_FEATHER_FILE = 'transcripts.feather'
    GENES_FEATHER_FILE = 'genes.feather'
    INTERVALS_FEATHER_FILE = 'intervals.feather'
    SJ_FEATHER_FILE = 'sj.feather'

    def __init__(self, index_dir: str):
        # index attributes
        self.index_dir = index_dir
        self.t_df = None
        self.g_df = None
        self.t_to_g_arr = None
        self.cr_indexes = None
        self.sj_map = None

    @property
    def num_transcripts(self):
        return len(self.t_df)
    
    @property
    def num_genes(self):
        # This is much faster than running np.unique() on a large array.
        return len(self.g_df)   

    @staticmethod
    def load(index_dir: str):
        self = TxIndex(index_dir)
        self.t_df = pd.read_feather(os.path.join(index_dir, TxIndex.TRANSCRIPTS_FEATHER_FILE))
        assert(all(self.t_df.index == self.t_df['t_index']))
        self.g_df = pd.read_feather(os.path.join(index_dir, TxIndex.GENES_FEATHER_FILE))
        assert(all(self.g_df.index == self.g_df['g_index']))
        
        # check dataframes
        # TODO: can remove this eventually
        for trow in self.t_df.itertuples(index=False):
            grow = self.g_df.iloc[trow.g_index]
            assert trow.g_id == grow.g_id
            assert trow.t_index in grow.t_index
            assert trow.t_id in grow.t_id

        # conversion from transcript to gene index
        # 1D numpy array, position-aligned for direct lookup
        self.t_to_g_arr = self.t_df['g_index'].values

        # conversion from transcript to strand
        # 1D numpy array, position-aligned for direct lookup
        self.t_to_strand_arr = self.t_df['strand'].values

        # conversion from gene to strand
        # 1D numpy array, position-aligned for direct lookup
        self.g_to_strand_arr = self.g_df['strand'].values

        # read intervals
        logging.debug('Reading intervals')
        df = pd.read_feather(os.path.join(index_dir, TxIndex.INTERVALS_FEATHER_FILE))
        # build interval index
        logging.debug('Building interval index')
        cr_indexes = [cgranges.cgranges() for itype in IntervalType]
        for row in df.itertuples(index=False):
            cr = cr_indexes[row.interval_type]
            cr.add(row.ref, row.start, row.end, row.t_index)
        for cr in cr_indexes:
            cr.index()
        logging.debug(f'Interval index size: {len(df)}')
        self.cr_indexes = cr_indexes

        # build splice junction index
        logging.debug(f'Reading splice junctions')
        df = pd.read_feather(os.path.join(index_dir, TxIndex.SJ_FEATHER_FILE))
        logging.debug(f'Building splice junction index')
        sj_map = {}
        for ref, start, end, strand, t_index, g_index in df[['ref', 'start', 'end', 'strand', 't_index', 'g_index']].itertuples(index=False):
            key = (ref, start, end, strand)
            if key not in sj_map:
                sj_map[key] = (set(), set())
            sj_map[key][0].add(t_index)
            sj_map[key][1].add(g_index)
        # convert sets to frozensets at the end for immutability
        self.sj_map = {
            key: (frozenset(val[0]), frozenset(val[1])) 
            for key, val in sj_map.items()
        }
        logging.debug(f'Splice junctions found: {len(self.sj_map)}')
        return self
    

    def search_intervals(self, ref, intervals, interval_type):
        """
        Progressive search for matches, prioritizing same-strand over antisense:
        1. Intersection of all sets (most specific)
        2. Intersection of non-empty sets (skip empty, less specific)
        3. Union of all sets (most sensitive)
        Searches same strand first; will only use antisense if same-strand yields empty.
        Returns:
            t_inds, g_inds: sets of transcript/gene indices
        """       
        cr = self.cr_indexes[interval_type]
        t_sets = []
        g_sets = []

        for start, end in intervals:
            # search for overlapping intervals
            t_set = frozenset(hit[2] for hit in cr.overlap(ref, start, end))
            t_sets.append(t_set)
            # gather gene indexes corresponding to the transcripts
            g_set = frozenset(self.t_to_g_arr[list(t_set)]) if t_set else frozenset()
            g_sets.append(g_set)           
        t_inds, g_inds = _merge_sets_with_relaxation(t_sets, g_sets)
        return t_inds, g_inds


    def search_splice_junctions(self, ref, sjs):
        """
        Progressive search for splice junction matches:
        1. Intersection of all sets (most specific)
        2. Intersection of non-empty sets (less specific)
        3. Union of all sets (most sensitive)
        Returns:
            t_inds, g_inds: sets of transcript/gene indices
        """
        # gather sets for each interval, may be empty
        tg_set_pairs = [
            self.sj_map.get((ref, start, end, strand), (frozenset(), frozenset())) 
            for start, end, strand in sjs
        ]
        if tg_set_pairs:
            t_sets, g_sets = zip(*tg_set_pairs)
        else:
            # ensure we have at least one empty set
            t_sets, g_sets = [frozenset()], [frozenset()]
        t_inds, g_inds = _merge_sets_with_relaxation(t_sets, g_sets)
        return t_inds, g_inds


    def search_intergenic(self, ref, intervals):
        '''return True if all intervals are intergenic'''
        cr = self.cr_indexes[IntervalType.INTERGENIC]
        for start, end in intervals:
            for _, _, t_index in cr.overlap(ref, start, end):
                assert t_index == -1, "Expected intergenic hit with t_index -1"
                return True
        return False