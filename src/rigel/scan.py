"""rigel.scan — CSR builder / buffer scanner.

Converts a ``FragmentBuffer`` into ``ScoredFragments`` for the EM solver.
The ``FragmentRouter`` class builds per-fragment scoring data, routes
unambiguous fragments to deterministic assignment, and accumulates
ambiguous fragments into CSR-formatted equivalence-class data for
the locus-level EM.
"""

import math

import numpy as np

from .buffer import FragmentBuffer, FRAG_UNAMBIG, FRAG_CHIMERIC
from .estimator import AbundanceEstimator, ScoredFragments
from .index import TranscriptIndex
from .scoring import FragmentScorer
from .splice import SPLICE_UNSPLICED, SPLICE_ANNOT
from .stats import PipelineStats
from .annotate import POOL_CODE_MRNA, POOL_CODE_CHIMERIC

# Default splice penalty for splice types not found in the penalty map.
_DEFAULT_SPLICE_PENALTY = 1.0


class FragmentRouter:
    """Builds the global CSR (ScoredFragments) from a FragmentBuffer.

    State that was previously captured via closures is now explicit
    instance attributes.  Every method is a candidate for Cython / C
    acceleration.

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
        strand_models,
        annotations=None,
    ):
        self.ctx = ctx
        self.estimator = estimator
        self.stats = stats
        self.index = index
        self.strand_models = strand_models
        self.annotations = annotations

        # Native scoring context (set by FragmentScorer.from_models)
        self._native_ctx = getattr(ctx, "_native_ctx", None)

    # ------------------------------------------------------------------
    # Main scan driver
    # ------------------------------------------------------------------

    def scan(self, buffer: FragmentBuffer, log_every: int) -> ScoredFragments:
        """Single pass over buffer.  Returns packed ScoredFragments.

        Deterministic-unique (SPLICED_ANNOT + FRAG_UNAMBIG) fragments
        are assigned directly via ``estimator.assign_unambig``.  All other
        exonic fragments build EM units.  Chimeric fragments are
        recorded in the annotation table (if active) and skipped.

        Chunks are streamed from the buffer one at a time so that only
        one spilled chunk is loaded from disk at any moment during
        conversion to contiguous numpy arrays.
        """
        native_scorer = self._native_ctx
        if native_scorer is None:
            raise RuntimeError("NativeFragmentScorer not available; cannot scan")
        return self._scan_native(buffer, log_every)

    def _scan_native(
        self,
        buffer: FragmentBuffer,
        log_every: int,
    ) -> ScoredFragments:
        """Chunk-level C++ scan path (fast).

        Streams chunks from the buffer iterator, converting each to
        contiguous numpy arrays and releasing the raw chunk before
        loading the next.  Spilled chunks are loaded from disk on
        demand, so peak memory during conversion is bounded to one
        extra chunk.
        """
        import logging

        logger = logging.getLogger(__name__)

        native_scorer = self._native_ctx
        estimator = self.estimator
        stats = self.stats
        index = self.index
        t_to_g = self.ctx.t_to_g
        annotations = self.annotations

        from .scoring import LOG_SAFE_FLOOR

        gdna_unspliced_penalty = self.ctx.gdna_splice_penalties.get(
            SPLICE_UNSPLICED,
            _DEFAULT_SPLICE_PENALTY,
        )
        gdna_log_penalty = math.log(max(gdna_unspliced_penalty, LOG_SAFE_FLOOR))

        t_strand_arr = np.ascontiguousarray(index.t_to_strand_arr, dtype=np.int8)

        # ---- Stream chunks: convert one at a time, release raw ----
        # Each iteration loads at most one spilled chunk from disk,
        # converts it to contiguous arrays, then lets the raw chunk
        # be garbage-collected before the next iteration.
        chunk_arrays = []
        chunk_sizes = []
        for chunk in buffer.iter_chunks():
            chunk_arrays.append(chunk.to_scoring_arrays())
            chunk_sizes.append(chunk.size)

        # ---- Fused two-pass C++ scoring ----
        result = native_scorer.fused_score_buffer(
            chunk_arrays,
            t_strand_arr,
            estimator.transcript_unspliced_sense,
            estimator.transcript_unspliced_antisense,
            estimator.transcript_intronic_sense,
            estimator.transcript_intronic_antisense,
            estimator.unambig_counts,
            gdna_log_penalty,
        )

        (
            offsets,
            t_indices,
            log_liks,
            count_cols,
            coverage_weights,
            tx_starts,
            tx_ends,
            locus_t_indices,
            locus_count_cols,
            is_spliced_raw,
            gdna_log_liks,
            genomic_footprints,
            frag_ids,
            frag_classes,
            splice_types,
            det_tids,
            det_fids,
            n_det,
            n_em_u,
            n_em_as,
            n_em_ao,
            n_gated,
            n_chim,
            n_mm,
        ) = result

        stats.deterministic_unambig_units += int(n_det)
        stats.em_routed_unambig_units += int(n_em_u)
        stats.em_routed_ambig_same_strand_units += int(n_em_as)
        stats.em_routed_ambig_opp_strand_units += int(n_em_ao)
        stats.n_gated_out += int(n_gated)
        stats.em_routed_multimapper_units += int(n_mm)

        # ---- Annotations (from chunk_arrays — no raw chunks) ----
        # Tuple layout: [6]=splice_type, [8]=fragment_classes, [9]=frag_id
        n_processed = 0
        for ci in range(len(chunk_arrays)):
            chunk_frag_classes = chunk_arrays[ci][8]

            if annotations is not None:
                chimeric_indices = np.where(chunk_frag_classes == FRAG_CHIMERIC)[0]
                for idx in chimeric_indices:
                    annotations.add(
                        frag_id=int(chunk_arrays[ci][9][idx]),
                        best_tid=-1,
                        best_gid=-1,
                        pool=POOL_CODE_CHIMERIC,
                        posterior=0.0,
                        frag_class=FRAG_CHIMERIC,
                        n_candidates=0,
                        splice_type=int(chunk_arrays[ci][6][idx]),
                    )

            n_processed += chunk_sizes[ci]
            if n_processed % log_every < chunk_sizes[ci]:
                logger.debug(f"  Scan: {n_processed:,} / {buffer.total_fragments:,}")

        # Det-unambig annotations
        if annotations is not None and len(det_tids) > 0:
            for j in range(len(det_tids)):
                tid = int(det_tids[j])
                gid = int(t_to_g[tid])
                annotations.add(
                    frag_id=int(det_fids[j]),
                    best_tid=tid,
                    best_gid=gid,
                    pool=POOL_CODE_MRNA,
                    posterior=1.0,
                    frag_class=FRAG_UNAMBIG,
                    n_candidates=1,
                    splice_type=int(SPLICE_ANNOT),
                )

        # Free chunk_arrays — no longer needed
        del chunk_arrays

        return ScoredFragments(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=log_liks,
            count_cols=count_cols,
            coverage_weights=coverage_weights,
            tx_starts=tx_starts,
            tx_ends=tx_ends,
            locus_t_indices=locus_t_indices,
            locus_count_cols=locus_count_cols,
            is_spliced=np.asarray(is_spliced_raw, dtype=np.int8).astype(bool),
            gdna_log_liks=gdna_log_liks,
            genomic_footprints=genomic_footprints,
            frag_ids=frag_ids,
            frag_class=frag_classes,
            splice_type=splice_types,
            n_units=int(len(offsets) - 1),
            n_candidates=int(len(t_indices)),
            nrna_base_index=self.ctx.nrna_base,
        )
