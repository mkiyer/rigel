"""How the pipeline routes a fragment into EM units, and what it counts on the way.

A multimapper's EM unit carries transcript candidates only, spliced or not — synthetic nascent spans
are ordinary transcripts, and the gDNA candidate is appended per locus inside the EM — and an unspliced
one carries a finite per-unit gDNA log-likelihood; a multimapper's gDNA likelihood is the mean over
the hits gDNA can explain, not a sum and not over `NH`; the route counters are exclusive, so
each unit is counted once and the totals mean something; and the `NM` penalty discriminates between
a multimapper's hits when enabled and is exactly inert when not.
"""

import math
from dataclasses import dataclass

import numpy as np
import pytest

from rigel.buffer import (
    FRAG_AMBIG_OPP_STRAND,
    FRAG_AMBIG_SAME_STRAND,
    FRAG_MULTIMAPPER,
    FRAG_UNAMBIG,
)
from rigel.splice import SpliceType
from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.frag_length_model import FragmentLengthModel
from rigel.scoring import FragmentScorer
from rigel.scan import FragmentRouter
from rigel.stats import PipelineStats
from rigel.strand_model import StrandModel, StrandModels
from rigel.types import Strand


@dataclass
class _BF:
    t_inds: np.ndarray
    splice_type: int
    align_strand: int
    frag_lengths: np.ndarray | None
    exon_bp: np.ndarray
    intron_bp: np.ndarray
    read_length: int
    nm: int = 0
    genomic_footprint: int = 200
    genomic_start: int = -1


class _Chunk:
    """Columnar mock chunk compatible with the C++ native scan path.

    Has ``t_offsets`` and other columnar arrays so ``_scan_native``
    fires.  Also supports ``__getitem__`` for multimapper handling.
    """

    def __init__(self, bfs, fragment_classes, frag_ids):
        self._bfs = bfs
        n = len(bfs)

        # Per-fragment columnar arrays
        self.splice_type = np.array([bf.splice_type for bf in bfs], dtype=np.uint8)
        self.align_strand = np.array([bf.align_strand for bf in bfs], dtype=np.uint8)
        self.fragment_classes = np.array(fragment_classes, dtype=np.uint8)
        self.frag_id = np.array(frag_ids, dtype=np.int64)
        self.read_length = np.array([bf.read_length for bf in bfs], dtype=np.uint16)
        self.genomic_footprint = np.array([bf.genomic_footprint for bf in bfs], dtype=np.int32)
        self.genomic_start = np.array([bf.genomic_start for bf in bfs], dtype=np.int32)
        self.nm = np.array([bf.nm for bf in bfs], dtype=np.uint16)
        self.size = n

        # CSR transcript indices (variable-length per fragment)
        offsets = [0]
        flat_t = []
        flat_fl = []
        flat_exon = []
        for bf in bfs:
            n_cand = len(bf.t_inds)
            flat_t.extend(bf.t_inds)
            if bf.frag_lengths is not None:
                flat_fl.extend(bf.frag_lengths)
            else:
                flat_fl.extend([200] * n_cand)
            flat_exon.extend(bf.exon_bp)
            offsets.append(len(flat_t))

        self.t_offsets = np.array(offsets, dtype=np.int32)
        self.t_indices = np.array(flat_t, dtype=np.int32)
        self.frag_lengths = np.array(flat_fl, dtype=np.int32)
        self.exon_bp = np.array(flat_exon, dtype=np.uint16)

    def to_scoring_arrays(self):
        return (
            np.ascontiguousarray(self.t_offsets, dtype=np.int32),
            np.ascontiguousarray(self.t_indices, dtype=np.int32),
            np.ascontiguousarray(self.frag_lengths, dtype=np.int32),
            np.ascontiguousarray(self.exon_bp, dtype=np.uint16),
            np.ascontiguousarray(self.splice_type, dtype=np.uint8),
            np.ascontiguousarray(self.align_strand, dtype=np.uint8),
            np.ascontiguousarray(self.fragment_classes, dtype=np.uint8),
            np.ascontiguousarray(self.frag_id, dtype=np.int64),
            np.ascontiguousarray(self.read_length, dtype=np.uint16),
            np.ascontiguousarray(self.genomic_footprint, dtype=np.int32),
            np.ascontiguousarray(self.genomic_start, dtype=np.int32),
            np.ascontiguousarray(self.nm, dtype=np.uint16),
        )


class _Buffer:
    def __init__(self, chunks):
        self._chunks = list(chunks)
        self.total_fragments = sum(ch.size for ch in self._chunks)

    def iter_chunks_consuming(self):
        while self._chunks:
            yield self._chunks.pop(0)


class _Index:
    def __init__(self, t_to_g, t_to_strand, g_to_strand):
        import pandas as pd

        self.t_to_g_arr = np.array(t_to_g, dtype=np.int64)
        self.t_to_strand_arr = np.array(t_to_strand, dtype=np.int8)
        self.g_to_strand_arr = np.array(g_to_strand, dtype=np.int8)
        self.num_transcripts = len(t_to_g)
        self.num_genes = len(g_to_strand)
        self.t_df = pd.DataFrame(
            {
                "t_id": [f"t{i}" for i in range(self.num_transcripts)],
                "ref": ["chr1"] * self.num_transcripts,
                "start": np.zeros(self.num_transcripts, dtype=np.int32),
                "end": np.full(self.num_transcripts, 10000, dtype=np.int32),
                "length": np.full(self.num_transcripts, 1000, dtype=np.int32),
                "is_nrna": np.zeros(self.num_transcripts, dtype=bool),
                "is_synthetic": np.zeros(self.num_transcripts, dtype=bool),
            }
        )

    def build_exon_csr(self):
        n_t = self.num_transcripts
        offsets = np.zeros(n_t + 1, dtype=np.int32)
        empty = np.empty(0, dtype=np.int32)
        return offsets, empty, empty, empty


def _make_env(index):
    strand_models = StrandModels(
        exonic_spliced=StrandModel.from_labels([int(Strand.POS)] * 20, [int(Strand.POS)] * 20)
    )
    # Only the SIZE is needed here — an int, rather than a whole length-model container threaded
    # through several frames as an object.
    max_frag_size = 1000
    estimator = AbundanceEstimator(index.num_transcripts, em_config=EMConfig(seed=1))
    stats = PipelineStats()
    return strand_models, max_frag_size, estimator, stats


def _scan_em_data(
    buffer,
    index,
    strand_models,
    max_frag_size,
    estimator,
    stats,
    log_every=1_000_000,
    *,
    gdna_splice_penalties=None,
    overhang_log_penalty=None,
    mismatch_log_penalty=None,
    annotations=None,
    fl=None,
):
    """Build EM data from a buffer using FragmentScorer + FragmentRouter."""
    # Routing tests don't exercise FL scoring: pass unfinalized (empty) FL models
    # so the scorer's FL LUT is inert (_log_prob is None), matching the prior use
    # of the container's empty rna/gdna models. A test that must see the length term passes ``fl``.
    empty_fl = fl if fl is not None else FragmentLengthModel(max_size=max_frag_size)
    ctx = FragmentScorer.from_models(
        strand_models,
        empty_fl,
        empty_fl,
        index,
        overhang_log_penalty=overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        gdna_splice_penalties=gdna_splice_penalties,
    )
    builder = FragmentRouter(
        ctx,
        estimator,
        stats,
        index,
        annotations=annotations,
    )
    return builder.scan(buffer, log_every)


def test_multimapper_spliced_annot_skips_shadows():
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, max_frag_size, estimator, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            align_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([100], dtype=np.int16),
            intron_bp=np.array([0], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.SPLICED_ANNOT),
            align_strand=int(Strand.POS),
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
        max_frag_size,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # Spliced multimappers: only transcript candidates
    assert np.all(em.t_indices < index.num_transcripts)


def test_multimapper_unspliced_has_transcript_candidates_only_and_a_gdna_log_lik():
    """An unspliced multimapper's unit carries transcript candidates only in the global CSR.

    Synthetic nRNA transcripts are ordinary transcripts, so nothing is appended for them, and the
    gDNA candidate is appended per locus inside the locus EM. What the unit does carry is a finite
    per-unit gDNA log-likelihood.
    """
    index = _Index(
        t_to_g=[0, 0],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, max_frag_size, estimator, stats = _make_env(index)

    bfs = [
        _BF(
            t_inds=np.array([0], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            align_strand=int(Strand.POS),
            frag_lengths=np.array([200], dtype=np.int32),
            exon_bp=np.array([90], dtype=np.int16),
            intron_bp=np.array([10], dtype=np.int16),
            read_length=100,
        ),
        _BF(
            t_inds=np.array([1], dtype=np.int32),
            splice_type=int(SpliceType.UNSPLICED),
            align_strand=int(Strand.POS),
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
        max_frag_size,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    assert stats.em_routed_multimapper_units == 1
    # All candidates are transcript-level (no nRNA shadows)
    assert np.all(em.t_indices < index.num_transcripts)
    # gDNA candidates NOT in global CSR (added per-locus)
    # gDNA log-likelihood is pre-computed per unit for unspliced
    assert em.gdna_log_liks[0] > -np.inf


def test_multimapper_gdna_likelihood_is_the_mean_over_its_gdna_explainable_hits():
    """Identical gDNA-eligible MM hits should collapse to one scalar.

    The scorer emits ``logsumexp(hit_scores) - log(nh_gdna)`` over the hits gDNA
    can explain. If two hits have identical gDNA scores, the emitted
    scalar must match the one-hit case rather than gaining ``log(2)`` mass.
    """
    index = _Index(
        t_to_g=[0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS)],
        g_to_strand=[int(Strand.POS), int(Strand.POS)],
    )

    def scan_for(bfs, frag_ids):
        strand_models, max_frag_size, estimator, stats = _make_env(index)
        chunk = _Chunk(
            bfs=bfs,
            fragment_classes=[FRAG_MULTIMAPPER] * len(bfs),
            frag_ids=frag_ids,
        )
        return _scan_em_data(
            _Buffer([chunk]),
            index,
            strand_models,
            max_frag_size,
            estimator,
            stats,
        )

    hit0 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
        exon_bp=np.array([90], dtype=np.int16),
        intron_bp=np.array([10], dtype=np.int16),
        read_length=100,
        genomic_footprint=200,
    )
    hit1 = _BF(
        t_inds=np.array([1], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
        exon_bp=np.array([90], dtype=np.int16),
        intron_bp=np.array([10], dtype=np.int16),
        read_length=100,
        genomic_footprint=200,
    )

    one_hit = scan_for([hit0], [21])
    two_hits = scan_for([hit0, hit1], [22, 22])

    assert one_hit.n_units == 1
    assert two_hits.n_units == 1
    assert two_hits.gdna_log_liks[0] == pytest.approx(one_hit.gdna_log_liks[0])


def test_route_counters_are_exclusive_per_unit():
    index = _Index(
        t_to_g=[0, 0, 1],
        t_to_strand=[int(Strand.POS), int(Strand.POS), int(Strand.NEG)],
        g_to_strand=[int(Strand.POS), int(Strand.NEG)],
    )
    strand_models, max_frag_size, estimator, stats = _make_env(index)

    bfs = [
        # Deterministic unique: FRAG_UNAMBIG + SPLICED_ANNOT
        _BF(
            np.array([0], dtype=np.int32),
            int(SpliceType.SPLICED_ANNOT),
            int(Strand.POS),
            np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16),
            np.array([0], dtype=np.int16),
            100,
        ),
        # Unique routed to EM: FRAG_UNAMBIG + UNSPLICED
        _BF(
            np.array([0], dtype=np.int32),
            int(SpliceType.UNSPLICED),
            int(Strand.POS),
            np.array([200], dtype=np.int32),
            np.array([100], dtype=np.int16),
            np.array([0], dtype=np.int16),
            100,
        ),
        # Isoform ambiguous
        _BF(
            np.array([0, 1], dtype=np.int32),
            int(SpliceType.UNSPLICED),
            int(Strand.POS),
            np.array([200, 200], dtype=np.int32),
            np.array([95, 90], dtype=np.int16),
            np.array([5, 10], dtype=np.int16),
            100,
        ),
        # Gene ambiguous
        _BF(
            np.array([0, 2], dtype=np.int32),
            int(SpliceType.UNSPLICED),
            int(Strand.POS),
            np.array([200, 200], dtype=np.int32),
            np.array([90, 85], dtype=np.int16),
            np.array([10, 15], dtype=np.int16),
            100,
        ),
        # Multimapper group (2 hits, same frag_id)
        _BF(
            np.array([0], dtype=np.int32),
            int(SpliceType.UNSPLICED),
            int(Strand.POS),
            np.array([200], dtype=np.int32),
            np.array([95], dtype=np.int16),
            np.array([5], dtype=np.int16),
            100,
        ),
        _BF(
            np.array([2], dtype=np.int32),
            int(SpliceType.UNSPLICED),
            int(Strand.POS),
            np.array([200], dtype=np.int32),
            np.array([90], dtype=np.int16),
            np.array([10], dtype=np.int16),
            100,
        ),
    ]

    chunk = _Chunk(
        bfs=bfs,
        fragment_classes=[
            FRAG_UNAMBIG,
            FRAG_UNAMBIG,
            FRAG_AMBIG_SAME_STRAND,
            FRAG_AMBIG_OPP_STRAND,
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
        max_frag_size,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert stats.deterministic_unambig_units == 1
    assert stats.em_routed_unambig_units == 1
    assert stats.em_routed_ambig_same_strand_units == 1
    assert stats.em_routed_ambig_opp_strand_units == 1
    assert stats.em_routed_multimapper_units == 1

    total_units = (
        stats.deterministic_unambig_units
        + stats.em_routed_unambig_units
        + stats.em_routed_ambig_same_strand_units
        + stats.em_routed_ambig_opp_strand_units
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
    strand_models, max_frag_size, estimator, stats = _make_env(index)

    # Hit A: maps to t0 with NM=0
    bf_a = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.POS),
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
        align_strand=int(Strand.POS),
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
        buffer,
        index,
        strand_models,
        max_frag_size,
        estimator,
        stats,
        log_every=1_000_000,
    )

    assert em.n_units == 1
    # Both transcripts should be candidates
    mRNA_t = em.t_indices
    mRNA_ll = em.log_liks
    assert len(mRNA_t) == 2
    # t0 (NM=0) should have higher log-lik than t1 (NM=4)
    t0_idx = np.where(mRNA_t == 0)[0][0]
    t1_idx = np.where(mRNA_t == 1)[0][0]
    # Default mismatch_alpha = 0.1 → log(0.1) ≈ -2.3 per mismatch
    # 4 mismatches → Δ ≈ -9.2
    assert mRNA_ll[t0_idx] > mRNA_ll[t1_idx]
    assert mRNA_ll[t0_idx] - mRNA_ll[t1_idx] == pytest.approx(
        -4 * math.log(0.1),
        abs=0.01,
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
    strand_models, max_frag_size, estimator_a, stats_a = _make_env(index)
    _, _, estimator_b, stats_b = _make_env(index)

    bf_nm0 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
        exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=0,
    )
    bf_nm5 = _BF(
        t_inds=np.array([0], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
        exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=5,
    )

    # Run with mismatch_alpha=1.0 (disabled)
    import math

    mismatch_lp = math.log(1.0)  # = 0.0

    chunk_a = _Chunk(bfs=[bf_nm0], fragment_classes=[FRAG_UNAMBIG], frag_ids=[1])
    chunk_b = _Chunk(bfs=[bf_nm5], fragment_classes=[FRAG_UNAMBIG], frag_ids=[1])

    em_a = _scan_em_data(
        _Buffer([chunk_a]),
        index,
        strand_models,
        max_frag_size,
        estimator_a,
        stats_a,
        log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )
    em_b = _scan_em_data(
        _Buffer([chunk_b]),
        index,
        strand_models,
        max_frag_size,
        estimator_b,
        stats_b,
        log_every=1_000_000,
        mismatch_log_penalty=mismatch_lp,
    )

    # Log-liks should be identical since penalty is 0
    assert em_a.log_liks[0] == pytest.approx(em_b.log_liks[0])


# ── which fragments may be gDNA: one rule, by alignment ───────────────────────────────────────────
#
# gDNA can explain an alignment it can produce: an unspliced one; an IMPLICIT one — an unspliced pair whose
# mate gap could hold an annotated intron, where the gap is unobserved and only the fragment length tells the
# readings apart; and an ARTIFACT one, whose every junction the blacklist rejected, which is unspliced. A
# multimapper can be gDNA if ANY of its alignments can, and its gDNA term is the mean over those alignments.
# gDNA's term is the unspliced term at the footprint whatever the label, and gDNA then COMPETES exactly as an
# RNA candidate does: it is pruned when its likelihood falls more than -log(pruning_min_posterior) below the
# unit's best RNA candidate. No length rule of its own: a footprint the length law cannot produce is pruned.


def _one_unit(bfs, fragment_classes, frag_ids, fl=None):
    index = _Index(
        t_to_g=[0, 0, 0],
        t_to_strand=[int(Strand.POS)] * 3,
        g_to_strand=[int(Strand.POS)],
    )
    strand_models, max_frag_size, estimator, stats = _make_env(index)
    em = _scan_em_data(
        _Buffer([_Chunk(bfs=bfs, fragment_classes=fragment_classes, frag_ids=frag_ids)]),
        index,
        strand_models,
        max_frag_size,
        estimator,
        stats,
        fl=fl,
    )
    assert em.n_units == 1
    return em, max_frag_size


def _hit(splice_type, footprint, t=0, nm=0):
    return _BF(
        t_inds=np.array([t], dtype=np.int32),
        splice_type=int(splice_type),
        align_strand=int(Strand.POS),
        frag_lengths=np.array([200], dtype=np.int32),
        exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        nm=nm,
        genomic_footprint=footprint,
    )


def _unspliced_term(footprint, nm=0, fl=None):
    em, _ = _one_unit([_hit(SpliceType.UNSPLICED, footprint, nm=nm)], [FRAG_UNAMBIG], [1], fl=fl)
    return float(em.gdna_log_liks[0])


def _live_fl():
    """A length law that varies with length on [50, 500] and has no mass beyond, as the ladder's does."""
    pmf = np.zeros(1001)
    pmf[50:501] = np.linspace(1.0, 2.0, 451)
    return FragmentLengthModel.from_pmf(pmf, 1000)


@pytest.mark.parametrize("footprint", [300, 5000])
def test_an_ARTIFACT_splice_is_unspliced_and_a_gdna_candidate_whatever_its_footprint(footprint):
    em, _ = _one_unit([_hit(SpliceType.SPLICE_ARTIFACT, footprint)], [FRAG_UNAMBIG], [1])
    assert not em.is_spliced[0]
    assert float(em.gdna_log_liks[0]) == _unspliced_term(footprint)


@pytest.mark.parametrize("splice_type", [SpliceType.SPLICED_ANNOT, SpliceType.SPLICED_UNANNOT])
def test_a_sequenced_junction_that_survived_the_blacklist_is_RNA_only(splice_type):
    em, _ = _one_unit([_hit(splice_type, 300)], [FRAG_AMBIG_SAME_STRAND], [1])
    assert em.is_spliced[0]
    assert em.gdna_log_liks[0] == -np.inf


@pytest.mark.parametrize("spliced", [SpliceType.SPLICED_ANNOT, SpliceType.SPLICED_UNANNOT])
def test_a_multimapper_with_ONE_unspliced_alignment_is_a_gdna_candidate(spliced):
    """One spliced alignment no longer removes gDNA from the group: the unspliced alignment is a place
    gDNA could have come from, and the gDNA term is taken over the eligible alignments alone."""
    em, _ = _one_unit(
        [_hit(spliced, 300, t=0), _hit(SpliceType.UNSPLICED, 250, t=1)],
        [FRAG_MULTIMAPPER] * 2,
        [7, 7],
    )
    assert not em.is_spliced[0]
    assert float(em.gdna_log_liks[0]) == _unspliced_term(250)


def test_a_multimapper_of_IMPLICIT_and_ARTIFACT_alignments_is_a_gdna_candidate():
    em, _ = _one_unit(
        [_hit(SpliceType.SPLICED_IMPLICIT, 400, t=0), _hit(SpliceType.SPLICE_ARTIFACT, 400, t=1)],
        [FRAG_MULTIMAPPER] * 2,
        [8, 8],
    )
    assert not em.is_spliced[0]
    assert float(em.gdna_log_liks[0]) == _unspliced_term(400)


@pytest.mark.parametrize(
    "splice_type", [SpliceType.UNSPLICED, SpliceType.SPLICED_IMPLICIT, SpliceType.SPLICE_ARTIFACT]
)
def test_the_gdna_term_is_the_unspliced_term_AT_THE_FOOTPRINT_whatever_the_label(splice_type):
    fl = _live_fl()
    assert _unspliced_term(250, fl=fl) != _unspliced_term(300, fl=fl)  # the law is live
    em, _ = _one_unit([_hit(splice_type, 300)], [FRAG_UNAMBIG], [1], fl=fl)
    assert float(em.gdna_log_liks[0]) == pytest.approx(fl._log_prob[300] + math.log(0.5), abs=1e-5)


def test_a_multimapper_s_gdna_term_excludes_its_spliced_hit_and_keeps_its_reported_type():
    fl = _live_fl()
    em, _ = _one_unit(
        [_hit(SpliceType.SPLICED_ANNOT, 400, t=0), _hit(SpliceType.UNSPLICED, 250, t=1, nm=2)],
        [FRAG_MULTIMAPPER] * 2,
        [7, 7],
        fl=fl,
    )
    assert float(em.gdna_log_liks[0]) == pytest.approx(_unspliced_term(250, nm=2, fl=fl), abs=1e-5)
    assert em.splice_type[0] == SpliceType.SPLICED_ANNOT


def test_a_multimapper_s_gdna_term_is_the_MEAN_of_its_eligible_hits_terms():
    fl = _live_fl()
    a = _unspliced_term(300, nm=0, fl=fl)
    b = _unspliced_term(300, nm=2, fl=fl)
    em, _ = _one_unit(
        [
            _hit(SpliceType.SPLICED_IMPLICIT, 300, t=0),
            _hit(SpliceType.SPLICE_ARTIFACT, 300, t=1, nm=2),
        ],
        [FRAG_MULTIMAPPER] * 2,
        [8, 8],
        fl=fl,
    )
    assert float(em.gdna_log_liks[0]) == pytest.approx(np.logaddexp(a, b) - math.log(2), abs=1e-5)


def _geometric_fl():
    """A law whose mass falls geometrically with length, so a footprint moves gDNA's term by whole nats."""
    pmf = np.zeros(1001)
    lengths = np.arange(50, 1001)
    pmf[50:] = np.exp(-lengths / 20.0)
    return FragmentLengthModel.from_pmf(pmf, 1000)


@pytest.mark.parametrize(
    "splice_type", [SpliceType.UNSPLICED, SpliceType.SPLICED_IMPLICIT, SpliceType.SPLICE_ARTIFACT]
)
def test_gdna_is_PRUNED_exactly_as_an_RNA_candidate_is(splice_type):
    """The unit's RNA candidate reads its own length (200); gDNA reads the footprint. gDNA stays a candidate
    iff its term is within -log(1e-4) of the best RNA candidate's, and that bound — not a length — decides:
    both outcomes occur across the footprints, and each matches the bound computed from the two terms."""
    fl = _geometric_fl()
    delta = -math.log(1e-4)
    kept, pruned = 0, 0
    for footprint in range(200, 700, 7):
        em, _ = _one_unit([_hit(splice_type, footprint)], [FRAG_UNAMBIG], [1], fl=fl)
        best_rna = float(np.max(em.log_liks))
        gdna_term = fl._log_prob[footprint] + math.log(0.5)
        competes = best_rna - gdna_term <= delta
        assert bool(em.is_spliced[0]) == (not competes), footprint
        if competes:
            kept += 1
            assert float(em.gdna_log_liks[0]) == pytest.approx(gdna_term, abs=1e-5)
        else:
            pruned += 1
            assert em.gdna_log_liks[0] == -np.inf
    assert kept and pruned


def test_an_ARTIFACT_whose_footprint_the_law_cannot_produce_is_pruned():
    """Its footprint still spans the rejected intron (the splicing-artifact work edits that); at 5,000 bp the
    law has no mass, so gDNA does not compete."""
    em, _ = _one_unit([_hit(SpliceType.SPLICE_ARTIFACT, 5000)], [FRAG_UNAMBIG], [1], fl=_live_fl())
    assert em.is_spliced[0]
    assert em.gdna_log_liks[0] == -np.inf


def test_a_multimapper_s_gdna_mean_takes_every_eligible_hit_even_one_the_law_cannot_produce():
    fl = _live_fl()
    a = fl._log_prob[1000] + math.log(0.5)  # an IMPLICIT hit at 1,000 bp: the floor
    b = _unspliced_term(250, nm=2, fl=fl)
    em, _ = _one_unit(
        [_hit(SpliceType.SPLICED_IMPLICIT, 1000, t=0), _hit(SpliceType.UNSPLICED, 250, t=1, nm=2)],
        [FRAG_MULTIMAPPER] * 2,
        [7, 7],
        fl=fl,
    )
    assert float(em.gdna_log_liks[0]) == pytest.approx(np.logaddexp(a, b) - math.log(2), abs=1e-5)


def test_a_multimapper_is_RNA_only_when_no_alignment_can_be_gdna_or_its_gdna_is_pruned():
    fl = _live_fl()
    certified, _ = _one_unit(
        [_hit(SpliceType.SPLICED_ANNOT, 300, t=0), _hit(SpliceType.SPLICED_UNANNOT, 300, t=1)],
        [FRAG_MULTIMAPPER] * 2,
        [9, 9],
        fl=fl,
    )
    pruned, _ = _one_unit(
        [_hit(SpliceType.SPLICED_ANNOT, 300, t=0), _hit(SpliceType.SPLICED_IMPLICIT, 1500, t=1)],
        [FRAG_MULTIMAPPER] * 2,
        [9, 9],
        fl=fl,
    )
    for em in (certified, pruned):
        assert em.is_spliced[0]
        assert em.gdna_log_liks[0] == -np.inf


# ---------------------------------------------------------------------------
# The coverage weight is read along the transcript, 5'→3'
# ---------------------------------------------------------------------------

_TX_START, _TX_LEN, _FLEN = 1000, 1000, 200
_ONE_EXON = [(_TX_START, _TX_START + _TX_LEN)]


class _ExonIndex(_Index):
    """Transcript 0 on the minus strand and transcript 1 on the plus strand, both with the same ``exons``
    (``[start, end)`` in genomic order), present in the exon CSR the scorer reads."""

    def __init__(self, exons):
        strands = [int(Strand.NEG), int(Strand.POS)]
        super().__init__(t_to_g=[0, 1], t_to_strand=strands, g_to_strand=strands)
        self.exons = [tuple(e) for e in exons]
        self.t_df["start"] = self.exons[0][0]
        self.t_df["end"] = self.exons[-1][1]
        self.t_df["length"] = sum(e - s for s, e in self.exons)

    def build_exon_csr(self):
        k = len(self.exons)
        starts = np.array([s for s, _ in self.exons] * 2, dtype=np.int32)
        ends = np.array([e for _, e in self.exons] * 2, dtype=np.int32)
        before = np.concatenate([[0], np.cumsum(ends[:k] - starts[:k])[:-1]]).astype(np.int32)
        return (
            np.array([0, k, 2 * k], dtype=np.int32),
            starts,
            ends,
            np.concatenate([before, before]),
        )


def _fragment(t, genomic_start, flen=_FLEN):
    """A ``flen``-bp fragment on transcript ``t`` (sense) whose genomic start is ``genomic_start``."""
    return _BF(
        t_inds=np.array([t], dtype=np.int32),
        splice_type=int(SpliceType.UNSPLICED),
        align_strand=int(Strand.NEG) if t == 0 else int(Strand.POS),
        frag_lengths=np.array([flen], dtype=np.int32),
        exon_bp=np.array([100], dtype=np.int16),
        intron_bp=np.array([0], dtype=np.int16),
        read_length=100,
        genomic_footprint=flen,
        genomic_start=genomic_start,
    )


def _weights(exons, bfs, fragment_classes, frag_ids):
    index = _ExonIndex(exons)
    strand_models, max_frag_size, estimator, stats = _make_env(index)
    em = _scan_em_data(
        _Buffer([_Chunk(bfs=bfs, fragment_classes=fragment_classes, frag_ids=frag_ids)]),
        index,
        strand_models,
        max_frag_size,
        estimator,
        stats,
    )
    return {
        (u, int(em.t_indices[j])): float(em.coverage_weights[j])
        for u in range(em.n_units)
        for j in range(int(em.offsets[u]), int(em.offsets[u + 1]))
    }


def _trapezoid_weight(lo, hi, length):
    """The coverage model the scorer documents, integrated exactly: under uniform fragmentation the coverage
    at transcript position x is min(x, w, L − x) with w = min(f, L/2), and a fragment on [lo, hi) weighs w over
    its mean coverage, floored at 1. The coverage is piecewise linear, so a trapezoid rule over its kinks is
    exact."""
    f = hi - lo
    if f <= 0:
        return 1.0
    w = min(f, length / 2)

    def cov(x):
        return min(x, w, length - x)

    points = sorted({lo, hi, *[k for k in (w, length - w) if lo < k < hi]})
    area = sum((b - a) * (cov(a) + cov(b)) / 2 for a, b in zip(points, points[1:]))
    return w / max(area / f, 1.0)


def _forward_start(genomic_start, exons):
    """The fragment's start in the transcript's forward coordinates, measured as the resolver measures a
    fragment's length (``resolve_context.h``'s ``tx_frag_length``): inside an exon it is the exon's offset; in an
    intron or before the first exon it sits the overhang's length before the next exon's start; past the last
    exon it sits that far past the transcript's end."""
    before = 0
    for s, e in exons:
        if genomic_start < s:
            return before - (s - genomic_start)
        if genomic_start < e:
            return before + (genomic_start - s)
        before += e - s
    return before + (genomic_start - exons[-1][1])


def _expected_weight(exons, genomic_start, flen, strand):
    """The coverage model on the interval the fragment occupies in the transcript's own 5'→3' coordinates, cut to
    the transcript: ``[a, a + f)`` forward from its start ``a``, the mirror ``[L − a − f, L − a)`` on the minus
    strand, whose 5'→3' runs against the genome."""
    length = sum(e - s for s, e in exons)
    a = _forward_start(genomic_start, exons)
    lo, hi = (a, a + flen) if strand == "plus" else (length - a - flen, length - a)
    return _trapezoid_weight(max(lo, 0), min(hi, length), length)


# Offsets of the fragment's genomic start into a 1,000-bp exon: both ramps (the first and last _FLEN bases), their
# edges, the plateau between, and fragments overhanging either end of the transcript.
_OFFSETS = [-150, -1, 0, 1, 50, 199, 200, 400, 600, 601, 750, 800, 801, 950]


@pytest.mark.parametrize("offset", _OFFSETS)
def test_a_minus_strand_fragment_weighs_what_its_plus_strand_mirror_weighs(offset):
    """The trapezoid coverage is read along the transcript 5'→3', and on the minus strand that runs against
    the genome: a fragment whose genomic start sits ``offset`` bp into the exon covers transcript positions
    ``[L − offset − f, L − offset)``, exactly the interval a plus-strand fragment at ``L − offset − f`` covers,
    so the genome's mirror image of a fragment weighs what the fragment weighs. Anchored on the flipped START
    instead, the interval ran off the transcript's 3' end: at an end of the minus-strand transcript the weight
    read 1 (the plateau's) instead of 2."""
    mirror = _TX_LEN - offset - _FLEN
    w = _weights(
        _ONE_EXON,
        [_fragment(0, _TX_START + offset), _fragment(1, _TX_START + mirror)],
        [FRAG_UNAMBIG, FRAG_UNAMBIG],
        [1, 2],
    )
    assert w[(0, 0)] == pytest.approx(w[(1, 1)], rel=1e-12)
    # At either end the fragment fills one ramp: coverage averages f/2 against the plateau's f, weight 2.
    if offset in (0, _TX_LEN - _FLEN):
        assert w[(0, 0)] == pytest.approx(2.0, rel=1e-12)


@pytest.mark.parametrize("offset", _OFFSETS)
def test_a_multimapper_s_hits_are_weighed_along_their_own_transcripts(offset):
    """The multimapper path weighs each hit on its own transcript by the same rule as a single fragment."""
    mirror = _TX_LEN - offset - _FLEN
    w = _weights(
        _ONE_EXON,
        [_fragment(0, _TX_START + offset), _fragment(1, _TX_START + mirror)],
        [FRAG_MULTIMAPPER] * 2,
        [3, 3],
    )
    assert w[(0, 0)] == pytest.approx(w[(0, 1)], rel=1e-12)
    if offset in (0, _TX_LEN - _FLEN):
        assert w[(0, 0)] == pytest.approx(2.0, rel=1e-12)


# (exons, genomic start, fragment length): every offset above; a transcript shorter than two fragments, where a
# fragment overhanging an end weighs differently cut (its part on the transcript) than shifted wholly onto it; and
# two exons, with starts before the first exon, inside each, in the intron — where the start is measured back from
# the next exon's start, as the resolver measures the fragment's length — and past the end.
_SHORT = [(1000, 1300)]
_TWO_EXONS = [(1000, 1300), (2000, 2150)]
_CASES = (
    [(_ONE_EXON, _TX_START + o, _FLEN) for o in _OFFSETS]
    + [(_SHORT, g, 200) for g in (850, 900, 1000, 1050, 1100, 1250)]
    + [(_TWO_EXONS, g, 110) for g in (950, 1000, 1200, 1250, 1300, 1990, 2000, 2040, 2100)]
)


@pytest.mark.parametrize("exons,genomic_start,flen", _CASES)
@pytest.mark.parametrize("strand", ["minus", "plus"])
def test_every_fragment_weighs_what_the_coverage_model_says(exons, genomic_start, flen, strand):
    """Against an independent evaluation of the coverage model on the interval the fragment occupies. The weight
    is stored as a float32, so the bound is float32's."""
    t = 0 if strand == "minus" else 1
    w = _weights(exons, [_fragment(t, genomic_start, flen)], [FRAG_UNAMBIG], [1])
    assert w[(0, t)] == pytest.approx(
        _expected_weight(exons, genomic_start, flen, strand), rel=1e-6
    )
