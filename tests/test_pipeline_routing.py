"""Tests for pipeline routing counters and multimapper shadow behavior."""

from dataclasses import dataclass

import numpy as np
import pytest

from hulkrna.buffer import (
    FRAG_GENE_AMBIG,
    FRAG_ISOFORM_AMBIG,
    FRAG_MULTIMAPPER,
    FRAG_UNIQUE,
)
from hulkrna.categories import SpliceType
from hulkrna.config import EMConfig
from hulkrna.estimator import AbundanceEstimator
from hulkrna.frag_length_model import FragmentLengthModels
from hulkrna.scoring import ScoringContext
from hulkrna.scan import EmDataBuilder
from hulkrna.stats import PipelineStats
from hulkrna.strand_model import StrandModels
from hulkrna.types import Strand


@dataclass
class _BF:
    t_inds: np.ndarray
    splice_type: int
    exon_strand: int
    frag_lengths: np.ndarray | None
    exon_bp: np.ndarray
    intron_bp: np.ndarray
    read_length: int
    nm: int = 0
    unambig_intron_bp: np.ndarray = None
    genomic_footprint: int = 200
    genomic_start: int = -1

    def __post_init__(self):
        if self.unambig_intron_bp is None:
            self.unambig_intron_bp = np.zeros_like(self.intron_bp)


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
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int8)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int8)
        self.num_transcripts = len(t_to_g)
        self.num_genes = len(g_to_strand)
        self.t_df = pd.DataFrame({
            "t_id": [f"t{i}" for i in range(self.num_transcripts)],
            "ref": ["chr1"] * self.num_transcripts,
            "start": np.zeros(self.num_transcripts, dtype=np.int32),
            "end": np.full(self.num_transcripts, 10000, dtype=np.int32),
            "length": np.full(self.num_transcripts, 1000, dtype=np.int32),
        })

    def get_exon_intervals(self, t_idx):
        return None


def _make_env(index):
    strand_models = StrandModels()
    # Provide enough observations so _best_rna_model() doesn't error.
    for _ in range(20):
        strand_models.exonic_spliced.observe(Strand.POS, Strand.POS)
    frag_length_models = FragmentLengthModels(max_size=1000)
    frag_length_models.observe(200, SpliceType.UNSPLICED)
    counter = AbundanceEstimator(index.num_transcripts, index.num_genes,
                                  em_config=EMConfig(seed=1))
    stats = PipelineStats()
    return strand_models, frag_length_models, counter, stats


def _scan_em_data(
    buffer, index, strand_models, frag_length_models, counter, stats,
    log_every=1_000_000, *, gdna_splice_penalties=None,
    overhang_log_penalty=None, mismatch_log_penalty=None, annotations=None,
):
    """Build EM data from a buffer using ScoringContext + EmDataBuilder."""
    ctx = ScoringContext.from_models(
        strand_models, frag_length_models, index, counter,
        overhang_log_penalty=overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        gdna_splice_penalties=gdna_splice_penalties,
    )
    builder = EmDataBuilder(
        ctx, counter, stats, index, strand_models,
        annotations=annotations,
    )
    return builder.scan(buffer, log_every)


def test_multimapper_spliced_annot_skips_shadows():
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, frag_length_models, counter, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            read_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[10, 10],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
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
    strand_models, frag_length_models, counter, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([90], dtype=np.int16),
            intron_bp=np.array([10], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            exon_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([95], dtype=np.int16),
            intron_bp=np.array([5], dtype=np.int16),
            read_length=100,
        ),
    ]
    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[11, 11],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
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
    strand_models, frag_length_models, counter, stats = _make_env(index)

    bfs = [
        # Deterministic unique: FRAG_UNIQUE + SPLICED_ANNOT
        _BF(np.array([0], dtype=np.int32), int(SpliceType.SPLICED_ANNOT), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Unique routed to EM: FRAG_UNIQUE + UNSPLICED
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16), np.array([0], dtype=np.int16), 100),
        # Isoform ambiguous
        _BF(np.array([0, 1], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200, 200], dtype=np.int32),
            np.array([95, 90], dtype=np.int16), np.array([5, 10], dtype=np.int16), 100),
        # Gene ambiguous
        _BF(np.array([0, 2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200, 200], dtype=np.int32),
            np.array([90, 85], dtype=np.int16), np.array([10, 15], dtype=np.int16), 100),
        # Multimapper group (2 hits, same frag_id)
        _BF(np.array([0], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
            np.array([95], dtype=np.int16), np.array([5], dtype=np.int16), 100),
        _BF(np.array([2], dtype=np.int32), int(SpliceType.UNSPLICED), int(Strand.POS), np.array([200], dtype=np.int32),
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

    em = _scan_em_data(
        buffer,
        index,
        strand_models,
        frag_length_models,
        counter,
        stats,
        log_every=1_000_000,
    )

    assert stats.deterministic_unique_units == 1
    assert stats.em_routed_unique_units == 1
    assert stats.em_routed_isoform_ambig_units == 1
    assert stats.em_routed_gene_ambig_units == 1
    assert stats.em_routed_multimapper_units == 1

    total_units = (
        stats.deterministic_unique_units
        + stats.em_routed_unique_units
        + stats.em_routed_isoform_ambig_units
        + stats.em_routed_gene_ambig_units
        + stats.em_routed_multimapper_units
    )
    assert total_units == 5
    assert em.n_units == 4


def test_nm_penalty_discriminates_multimapper_hits():
    """Multimapper with different NM values → NM penalty produces different log-liks.

    Hit A: maps to t0 with NM=0 (perfect match)
    Hit B: maps to t1 with NM=4 (4 mismatches)
    Both have identical geometry (oh=0), strand, fragment length.
    The NM penalty should make t0 score significantly higher than t1.
    """
    import math
    index = _Index(
        t_to_g=[0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS), int(Strand.POS)],
    )
    strand_models, frag_length_models, counter, stats = _make_env(index)

    # Hit A: maps to t0 with NM=0
    bf_a = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=0,
    )
    # Hit B: maps to t1 with NM=4
    bf_b = _BF(
        t_inds=np.array([1], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=4,
    )

    chunk = _Chunk(
        bfs=[bf_a, bf_b],
        fragment_classes=[FRAG_MULTIMAPPER, FRAG_MULTIMAPPER],
        frag_ids=[1, 1],
    )
    buffer = _Buffer([chunk])

    em = _scan_em_data(
        buffer, index, strand_models, frag_length_models,
        counter, stats, log_every=1_000_000,
    )

    assert em.n_units == 1
    # Both transcripts should be candidates
    mRNA_mask = em.t_indices < counter.nrna_base_index
    mRNA_t = em.t_indices[mRNA_mask]
    mRNA_ll = em.log_liks[mRNA_mask]
    assert len(mRNA_t) == 2
    # t0 (NM=0) should have higher log-lik than t1 (NM=4)
    t0_idx = np.where(mRNA_t == 0)[0][0]
    t1_idx = np.where(mRNA_t == 1)[0][0]
    # Default mismatch_alpha = 0.1 → log(0.1) ≈ -2.3 per mismatch
    # 4 mismatches → Δ ≈ -9.2
    assert mRNA_ll[t0_idx] > mRNA_ll[t1_idx]
    assert mRNA_ll[t0_idx] - mRNA_ll[t1_idx] == pytest.approx(
        -4 * math.log(0.1), abs=0.01,
    )


def test_nm_penalty_zero_when_disabled():
    """With mismatch_alpha=1.0 (no penalty), NM has no effect on log-liks.

    Two identical fragments differing only in NM should score identically.
    """
    index = _Index(
        t_to_g=[0],
        t_to_strand=[int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, frag_length_models, counter_a, stats_a = _make_env(index)
    _, _, counter_b, stats_b = _make_env(index)

    bf_nm0 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=0,
    )
    bf_nm5 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        exon_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=5,
    )

    # Run with mismatch_alpha=1.0 (disabled)
    import math
    mismatch_lp = math.log(1.0)  # = 0.0

    chunk_a = _Chunk(bfs=[bf_nm0], fragment_classes=[FRAG_UNIQUE], frag_ids=[1])
    chunk_b = _Chunk(bfs=[bf_nm5], fragment_classes=[FRAG_UNIQUE], frag_ids=[1])

    em_a = _scan_em_data(
        _Buffer([chunk_a]), index, strand_models, frag_length_models,
        counter_a, stats_a, log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )
    em_b = _scan_em_data(
        _Buffer([chunk_b]), index, strand_models, frag_length_models,
        counter_b, stats_b, log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )

    # Log-liks should be identical since penalty is 0
    assert em_a.log_liks[0] == pytest.approx(em_b.log_liks[0])
