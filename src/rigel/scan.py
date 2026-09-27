"""rigel.scan — CSR builder / buffer scanner.

Converts a ``FragmentBuffer`` into ``ScoredFragments`` for the EM solver.
The ``FragmentRouter`` class builds per-fragment scoring data, routes
unambiguous fragments to deterministic assignment, and accumulates
ambiguous fragments into CSR-formatted equivalence-class data for
the locus-level EM.
"""

import logging

import numpy as np

from .buffer import FragmentBuffer, FRAG_UNAMBIG, FRAG_CHIMERIC
from .estimator import AbundanceEstimator
from .scored_fragments import ScoredFragments
from .index import TranscriptIndex
from .scoring import FragmentScorer
from .splice import SPLICE_ANNOT
from .stats import PipelineStats
from .annotate import AF_CHIMERIC, winner_flag

logger = logging.getLogger(__name__)


class FragmentRouter:
    """Builds the global CSR (ScoredFragments) from a FragmentBuffer.

    Drives the native StreamingScorer over buffer chunks: it scores each
    fragment against its compatible transcripts in C++, routes
    deterministic-unique fragments to direct assignment, and packs the
    remaining ambiguous fragments into the global CSR ``ScoredFragments``
    for the locus-level EM.

    Parameters
    ----------
    ctx : FragmentScorer
        Pre-computed scoring parameters.
    estimator : AbundanceEstimator
        Accumulates unambig counts during scan.
    stats : PipelineStats
        Pipeline statistics accumulator.
    index : TranscriptIndex
        Reference index.
    annotations : AnnotationTable or None
        If provided, deterministic-unambig and chimeric fragments are
        annotated immediately during the scan.
    """

    def __init__(
        self,
        ctx: FragmentScorer,
        estimator: AbundanceEstimator,
        stats: PipelineStats,
        index: TranscriptIndex,
        annotations=None,
    ):
        self.ctx = ctx
        self.estimator = estimator
        self.stats = stats
        self.index = index
        self.annotations = annotations

    # ------------------------------------------------------------------
    # Main scan driver
    # ------------------------------------------------------------------

    def scan(self, buffer: FragmentBuffer, log_every: int) -> ScoredFragments:
        """Single pass over buffer.  Returns packed ScoredFragments.

        Deterministic-unique (SPLICED_ANNOT + FRAG_UNAMBIG) fragments
        are assigned directly via the C++ scoring path.  All other
        exonic fragments build EM units.  Chimeric fragments are
        recorded in the annotation table (if active) and skipped.

        Chunks are consumed from the buffer in order and scored by the
        native ``StreamingScorer``.
        """
        native_scorer = self.ctx._native_ctx
        estimator = self.estimator
        stats = self.stats
        index = self.index
        t_to_g = self.ctx.t_to_g
        annotations = self.annotations

        from .native import StreamingScorer

        t_strand_arr = np.ascontiguousarray(index.t_to_strand_arr, dtype=np.int8)

        # ---- Streaming scoring: one chunk at a time ----
        scorer = StreamingScorer(
            native_scorer,
            t_strand_arr,
            estimator.unambig_counts,
        )

        n_processed = 0
        n_total = buffer.total_fragments
        for chunk in buffer.iter_chunks_consuming():
            arrays = chunk.to_scoring_arrays()
            scorer.score_chunk(arrays)
            n_processed += chunk.size
            if n_processed % log_every < chunk.size:
                logger.debug(f"  Scan: {n_processed:,} / {n_total:,}")
            del arrays, chunk

        result = scorer.finish()

        (
            offsets,
            t_indices,
            log_liks,
            count_cols,
            coverage_weights,
            is_spliced_raw,
            gdna_log_liks,
            frag_ids,
            frag_classes,
            splice_types,
            det_tids,
            det_fids,
            chim_fids,
            chim_stypes,
            n_det,
            n_em_u,
            n_em_as,
            n_em_ao,
            n_gated,
            _,
            n_mm,
        ) = result

        stats.deterministic_unambig_units += int(n_det)
        stats.em_routed_unambig_units += int(n_em_u)
        stats.em_routed_ambig_same_strand_units += int(n_em_as)
        stats.em_routed_ambig_opp_strand_units += int(n_em_ao)
        stats.n_gated_out += int(n_gated)
        stats.em_routed_multimapper_units += int(n_mm)

        # ---- Chimeric annotations (from C++ accumulators) ----
        if annotations is not None and len(chim_fids) > 0:
            for j in range(len(chim_fids)):
                annotations.add(
                    frag_id=int(chim_fids[j]),
                    best_tid=-1,
                    best_gid=-1,
                    tx_flags=AF_CHIMERIC,
                    posterior=0.0,
                    frag_class=FRAG_CHIMERIC,
                    n_candidates=0,
                    splice_type=int(chim_stypes[j]),
                )

        # Det-unambig annotations
        if annotations is not None and len(det_tids) > 0:
            is_nrna_arr = index.t_df["is_nrna"].values.astype(bool)
            is_synth_arr = index.t_df["is_synthetic"].values.astype(bool)
            for j in range(len(det_tids)):
                tid = int(det_tids[j])
                gid = int(t_to_g[tid])
                flags = winner_flag(bool(is_nrna_arr[tid]), bool(is_synth_arr[tid]))
                annotations.add(
                    frag_id=int(det_fids[j]),
                    best_tid=tid,
                    best_gid=gid,
                    tx_flags=flags,
                    posterior=1.0,
                    frag_class=FRAG_UNAMBIG,
                    n_candidates=1,
                    splice_type=int(SPLICE_ANNOT),
                )

        return ScoredFragments(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=log_liks,
            count_cols=count_cols,
            coverage_weights=coverage_weights,
            is_spliced=np.asarray(is_spliced_raw, dtype=np.int8).astype(bool),
            gdna_log_liks=gdna_log_liks,
            frag_ids=frag_ids,
            frag_class=frag_classes,
            splice_type=splice_types,
            n_units=int(len(offsets) - 1),
            n_candidates=int(len(t_indices)),
        )
