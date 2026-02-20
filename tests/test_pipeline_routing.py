"""Tests for pipeline routing counters and multimapper shadow behavior."""

from dataclasses import dataclass

import numpy as np

from hulkrna.buffer import (
    FRAG_GENE_AMBIG,
    FRAG_ISOFORM_AMBIG,
    FRAG_MULTIMAPPER,
    FRAG_UNIQUE,
)
from hulkrna.categories import SpliceType
from hulkrna.counter import ReadCounter
from hulkrna.insert_model import InsertSizeModels
from hulkrna.pipeline import _scan_and_build_em_data
from hulkrna.stats import PipelineStats
from hulkrna.strand_model import StrandModels
from hulkrna.types import Strand


@dataclass
class _BF:
    t_inds: np.ndarray
    splice_type: int
    exon_strand: int
    insert_size: int
    exon_bp: np.ndarray
    intron_bp: np.ndarray
    frag_length: int


class _Chunk:
    def __init__(self, bfs, fragment_classes, frag_ids):
        self._bfs = bfs
        self.fragment_classes = np.array(fragment_classes, dtype=np.uint8)
        self.frag_id = np.array(frag_ids, dtype=np.int64)
        self.size = len(bfs)

    def __getitem__(self, idx):
        return self._bfs[idx]


class _Buffer:
    def __init__(self, chunks):
        self._chunks = chunks
        self.total_fragments = sum(ch.size for ch in chunks)

    def iter_chunks(self):
        yield from self._chunks


class _Index:
    def __init__(self, t_to_g, t_to_strand, g_to_strand):
        import pandas as pd
        self.t_to_g_arr = np.array(t_to_g, dtype=np.int64)
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int64)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int64)
        self.num_transcripts = len(t_to_g)
        self.num_genes = len(g_to_strand)
        self.t_df = pd.DataFrame({
            "t_id": [f"t{i}" for i in range(self.num_transcripts)],
            "ref": ["chr1"] * self.num_transcripts,
        })


def _make_env(index):
    strand_models = StrandModels()
    # Provide enough observations so _best_rna_model() doesn't error.
    for _ in range(20):
        strand_models.exonic_spliced.observe(Strand.POS, Strand.POS)
    insert_models = InsertSizeModels(max_size=1000)
    insert_models.observe(200, SpliceType.UNSPLICED)
    counter = ReadCounter(index.num_transcripts, index.num_genes, seed=1)
    stats = PipelineStats()
    return strand_models, insert_models, counter, stats


def test_multimapper_spliced_annot_skips_shadows():
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, insert_models, counter, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            insert_size=200,
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            frag_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            insert_size=200,
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            frag_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[10, 10],
    )
    buffer = _Buffer([chunk])

    em = _scan_and_build_em_data(
        buffer,
        index,
        strand_models,
        insert_models,
        counter,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # Spliced multimappers: only mRNA candidates (no nRNA/gDNA)
    assert np.all(em.t_indices < em.nrna_base_index)


def test_multimapper_unspliced_adds_nrna_shadows():
    """Unspliced multimappers get nRNA candidates in global CSR.

    gDNA candidates are no longer in the global CSR — they are added
    per-locus during locus EM construction.
    """
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, insert_models, counter, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            insert_size=200,
            exon_bp=np.array([90], dtype=np.int16),
            intron_bp=np.array([10], dtype=np.int16),
            frag_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            insert_size=200,
            exon_bp=np.array([95], dtype=np.int16),
            intron_bp=np.array([5], dtype=np.int16),
            frag_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[11, 11],
    )
    buffer = _Buffer([chunk])

    em = _scan_and_build_em_data(
        buffer,
        index,
        strand_models,
        insert_models,
        counter,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # nRNA candidates present (index >= nrna_base_index)
    assert np.any(em.t_indices >= em.nrna_base_index)
    # gDNA candidates NOT in global CSR (added per-locus)
    # gDNA log-likelihood is pre-computed per unit for unspliced
    assert em.gdna_log_liks[0] > -np.inf


def test_route_counters_are_exclusive_per_unit():
    index = _Index(
        t_to_g=[0, 0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS), int(Strand.NEG)],
        g_to_strand=[int(Strand.POS), int(Strand.NEG)],
    )
    strand_models, insert_models, counter, stats = _make_env(index)

    bfs = [
        # Deterministic unique: FRAG_UNIQUE + SPLICED_ANNOT
        _BF(np.array([0], dtype=np.int32), int(SpliceType.SPLICED_ANNOT), int(Strand.POS), 200,
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Unique routed to EM: FRAG_UNIQUE + UNSPLICED
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), 200,
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Isoform ambiguous
        _BF(np.array([0, 1], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), 200,
            np.array([95, 90], dtype=np.int16), np.array([5, 10], dtype=np.int16), 100),
        # Gene ambiguous
        _BF(np.array([0, 2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), 200,
            np.array([90, 85], dtype=np.int16), np.array([10, 15], dtype=np.int16), 100),
        # Multimapper group (2 hits, same frag_id)
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), 200,
            np.array([95], dtype=np.int16), np.array([5], dtype=np.int16), 100),
        _BF(np.array([2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), 200,
            np.array([90], dtype=np.int16), np.array([10], dtype=np.int16), 100),
    ]

    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[
            FRAG_UNIQUE,
            FRAG_UNIQUE,
            FRAG_ISOFORM_AMBIG,
            FRAG_GENE_AMBIG,
            FRAG_MULTIMAPPER,
            FRAG_MULTIMAPPER,
        ],
        frag_ids=[1, 2, 3, 4, 5, 5],
    )
    buffer = _Buffer([chunk])

    em = _scan_and_build_em_data(
        buffer,
        index,
        strand_models,
        insert_models,
        counter,
        stats,
        log_every=1_000_000,
    )

    assert stats.deterministic_unique_units == 1
    assert stats.em_routed_unique_units == 1
    assert stats.em_routed_isoform_ambig_units == 1
    assert stats.em_routed_gene_ambig_units == 1
    assert stats.em_routed_multimapper_units == 1

    # Back-compat aliases map to authoritative route counters
    assert stats.n_counted_truly_unique == 1
    assert stats.n_counted_isoform_ambig == 1
    assert stats.n_counted_gene_ambig == 1
    assert stats.n_counted_multimapper == 1

    total_units = (
        stats.deterministic_unique_units
        + stats.em_routed_unique_units
        + stats.em_routed_isoform_ambig_units
        + stats.em_routed_gene_ambig_units
        + stats.em_routed_multimapper_units
    )
    assert total_units == 5
    assert em.n_units == 4
