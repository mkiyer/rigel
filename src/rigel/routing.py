"""rigel.routing — score every buffered fragment and route it.

:class:`FragmentRouter` turns a ``FragmentBuffer`` into ``ScoredFragments``: it drives the native scorer over the
buffer, assigns unambiguous fragments deterministically, and packs the rest into the CSR the locus-level EM reads.
"""

import logging

import numpy as np

from .buffer import FragmentBuffer, FRAG_MULTIMAPPER, FRAG_UNAMBIG
from .calibration.capture_eff_length import place_log_weights
from .estimator import AbundanceEstimator
from .scored_fragments import ScoredFragments
from .index import TranscriptIndex
from .stats import PipelineStats
from .annotate import AF_GDNA_INTERGENIC, AF_UNRESOLVED, winner_flag

logger = logging.getLogger(__name__)


class FragmentRouter:
    """Builds the global CSR (``ScoredFragments``) from a ``FragmentBuffer`` with the native ``StreamingScorer``.

    Parameters
    ----------
    scorer : NativeFragmentScorer
        The run's scorer (``rigel.scoring.build_scorer``).
    estimator : AbundanceEstimator
        Accumulates unambig counts during routing.
    stats : PipelineStats
        Pipeline statistics accumulator.
    index : TranscriptIndex
        Reference index.
    annotations : AnnotationTable or None
        If provided, deterministic-unambig fragments are annotated as they are routed.
    place_efficiencies : tuple or None
        ``(c_region, c_boundary, c_sj)`` (``capture_eff_length.place_efficiencies``): each multimapper's place weighs
        its candidates by its deposit masses at them (``capture_eff_length.place_log_weights``). ``None`` weighs every
        place 1.
    """

    def __init__(
        self,
        scorer,
        estimator: AbundanceEstimator,
        stats: PipelineStats,
        index: TranscriptIndex,
        annotations=None,
        place_efficiencies=None,
    ):
        self.scorer = scorer
        self.estimator = estimator
        self.stats = stats
        self.index = index
        self.annotations = annotations
        self.place_efficiencies = place_efficiencies

    # ------------------------------------------------------------------
    # The routing pass
    # ------------------------------------------------------------------

    def route(self, buffer: FragmentBuffer, log_every: int) -> ScoredFragments:
        """Score and route every buffered fragment in one pass, consuming the buffer's chunks in order.

        A deterministic-unique fragment (spliced-annotated and ``FRAG_UNAMBIG``) is assigned directly in C++; every
        other genic fragment becomes an EM unit.
        """
        native_scorer = self.scorer
        estimator = self.estimator
        stats = self.stats
        index = self.index
        t_to_g = index.t_to_g_arr
        annotations = self.annotations

        from .native import StreamingScorer

        scorer = StreamingScorer(native_scorer, estimator.unambig_counts)

        n_processed = 0
        n_total = buffer.total_fragments
        for chunk in buffer.iter_chunks_consuming():
            arrays = chunk.to_scoring_arrays()
            if self.place_efficiencies is not None and chunk.place_offsets is not None:
                c_region, c_boundary, c_sj = self.place_efficiencies
                arrays["log_weight"] = place_log_weights(
                    chunk.place_offsets,
                    chunk.place_kind,
                    chunk.place_index,
                    chunk.place_mass,
                    intergenic=np.diff(chunk.t_offsets) == 0,
                    c_region=c_region,
                    c_boundary=c_boundary,
                    c_sj=c_sj,
                )
            scorer.score_chunk(arrays)
            n_processed += chunk.size
            if n_processed % log_every < chunk.size:
                logger.debug(f"  Routed: {n_processed:,} / {n_total:,}")
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
            gated_fids,
            gated_classes,
            n_det,
            n_em_u,
            n_em_as,
            n_em_ao,
            n_gated,
            n_mm,
            intergenic_log_liks,
            intergenic_fids,
        ) = result

        stats.deterministic_unambig_units += int(n_det)
        stats.em_routed_unambig_units += int(n_em_u)
        stats.em_routed_ambig_same_strand_units += int(n_em_as)
        stats.em_routed_ambig_opp_strand_units += int(n_em_ao)
        stats.n_gated_out += int(n_gated)
        stats.em_routed_multimapper_units += int(n_mm)
        # A multimapper whose only places left after the gates are intergenic counts as intergenic
        stats.n_intergenic += len(intergenic_fids)

        # deterministic-unambiguous annotations
        if annotations is not None and len(det_tids) > 0:
            is_synth_arr = index.t_df["is_synthetic"].values.astype(bool)
            for j in range(len(det_tids)):
                tid = int(det_tids[j])
                gid = int(t_to_g[tid])
                flags = winner_flag(bool(is_synth_arr[tid]))
                annotations.add(
                    frag_id=int(det_fids[j]),
                    best_tid=tid,
                    best_gid=gid,
                    tx_flags=flags,
                    posterior=1.0,
                    frag_class=FRAG_UNAMBIG,
                    n_candidates=1,
                )

        # A fragment a hard gate dropped reaches no pool: unresolved, with its class
        if annotations is not None and len(gated_fids) > 0:
            n = len(gated_fids)
            annotations.add_batch(
                frag_ids=np.asarray(gated_fids, dtype=np.int64),
                best_tids=np.full(n, -1, dtype=np.int32),
                best_gids=np.full(n, -1, dtype=np.int32),
                tx_flags=np.full(n, AF_UNRESOLVED, dtype=np.uint8),
                posteriors=np.zeros(n, dtype=np.float32),
                frag_classes=np.asarray(gated_classes, dtype=np.int8),
                n_candidates=np.zeros(n, dtype=np.int16),
            )
        if annotations is not None and len(intergenic_fids) > 0:
            n = len(intergenic_fids)
            annotations.add_batch(
                frag_ids=np.asarray(intergenic_fids, dtype=np.int64),
                best_tids=np.full(n, -1, dtype=np.int32),
                best_gids=np.full(n, -1, dtype=np.int32),
                tx_flags=np.full(n, AF_GDNA_INTERGENIC, dtype=np.uint8),
                posteriors=np.ones(n, dtype=np.float32),
                frag_classes=np.full(n, FRAG_MULTIMAPPER, dtype=np.int8),
                n_candidates=np.zeros(n, dtype=np.int16),
            )

        return ScoredFragments(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=log_liks,
            count_cols=count_cols,
            coverage_weights=coverage_weights,
            is_spliced=np.asarray(is_spliced_raw, dtype=np.int8).astype(bool),
            gdna_log_liks=gdna_log_liks,
            intergenic_log_liks=intergenic_log_liks,
            frag_ids=frag_ids,
            frag_class=frag_classes,
            splice_type=splice_types,
            n_units=int(len(offsets) - 1),
            n_candidates=int(len(t_indices)),
        )
