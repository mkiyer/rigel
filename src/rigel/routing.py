"""rigel.routing — score every buffered fragment and route it.

:class:`FragmentRouter` turns a ``FragmentBuffer`` into ``ScoredFragments``: it drives the native scorer over the
buffer, assigns unambiguous fragments deterministically, and packs the rest into the CSR the locus-level EM reads.
"""

import logging
from typing import NamedTuple

import numpy as np

from .buffer import FragmentBuffer, FRAG_MULTIMAPPER, FRAG_UNAMBIG
from .calibration.capture_eff_length import place_log_weights
from .estimator import AbundanceEstimator
from .scored_fragments import ScoredFragments
from .index import TranscriptIndex
from .stats import PipelineStats
from .annotate import AF_GDNA_INTERGENIC, AF_UNRESOLVED, winner_flags

logger = logging.getLogger(__name__)


class _Routed(NamedTuple):
    """What ``StreamingScorer.finish`` returns, in its order: the EM units' CSR and per-unit arrays, the
    deterministic fragments, the gated ones, the routing counts and the intergenic multimappers."""

    offsets: np.ndarray
    t_indices: np.ndarray
    log_liks: np.ndarray
    count_cols: np.ndarray
    coverage_weights: np.ndarray
    is_spliced: np.ndarray
    gdna_log_liks: np.ndarray
    frag_ids: np.ndarray
    frag_classes: np.ndarray
    splice_types: np.ndarray
    #: the deterministic-unique fragments' transcripts and fragment ids
    det_tids: np.ndarray
    det_fids: np.ndarray
    #: the fragments a hard gate dropped, and their classes
    gated_fids: np.ndarray
    gated_classes: np.ndarray
    n_deterministic: int
    n_em_unambig: int
    n_em_ambig_same_strand: int
    n_em_ambig_opp_strand: int
    n_gated: int
    n_em_multimapper: int
    intergenic_log_liks: np.ndarray
    #: the multimappers whose only places left after the gates are intergenic
    intergenic_fids: np.ndarray


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
        If provided, every fragment routing decides on its own — deterministic, gated or intergenic — is annotated.
    place_efficiencies : tuple or None
        ``capture_eff_length.place_efficiencies``' tuple, the efficiencies and each junction's two boundaries: each
        multimapper's place weighs
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

    def route(self, buffer: FragmentBuffer, log_every: int) -> ScoredFragments:
        """Score and route every buffered fragment in one pass, consuming the buffer's chunks in order.

        A deterministic-unique fragment (spliced-annotated and ``FRAG_UNAMBIG``) is assigned directly in C++; every
        other genic fragment becomes an EM unit.
        """
        from .native import StreamingScorer

        scorer = StreamingScorer(self.scorer, self.estimator.unambig_counts)
        n_processed = 0
        n_total = buffer.total_fragments
        for chunk in buffer.iter_chunks_consuming():
            scorer.score_chunk(self._scoring_arrays(chunk))
            n_processed += chunk.size
            if n_processed % log_every < chunk.size:
                logger.debug(f"  Routed: {n_processed:,} / {n_total:,}")
            del chunk

        routed = _Routed(*scorer.finish())
        self._count(routed)
        if self.annotations is not None:
            self._annotate(routed)
        return ScoredFragments(
            offsets=routed.offsets,
            t_indices=routed.t_indices,
            log_liks=routed.log_liks,
            count_cols=routed.count_cols,
            coverage_weights=routed.coverage_weights,
            is_spliced=np.asarray(routed.is_spliced, dtype=np.int8).astype(bool),
            gdna_log_liks=routed.gdna_log_liks,
            intergenic_log_liks=routed.intergenic_log_liks,
            frag_ids=routed.frag_ids,
            frag_class=routed.frag_classes,
            splice_type=routed.splice_types,
            n_units=int(len(routed.offsets) - 1),
            n_candidates=int(len(routed.t_indices)),
        )

    def _scoring_arrays(self, chunk) -> dict:
        """The chunk's arrays as the scorer reads them, with each multimapper place's log capture weight when the run
        weighs places."""
        arrays = chunk.to_scoring_arrays()
        if self.place_efficiencies is not None and chunk.place_offsets is not None:
            (
                c_region,
                c_boundary,
                c_sj,
                c_sj_crossing,
                sj_low_boundary,
                sj_high_boundary,
            ) = self.place_efficiencies
            arrays["log_weight"] = place_log_weights(
                chunk.place_offsets,
                chunk.place_kind,
                chunk.place_index,
                chunk.place_mass,
                intergenic=np.diff(chunk.t_offsets) == 0,
                c_region=c_region,
                c_boundary=c_boundary,
                c_sj=c_sj,
                c_sj_crossing=c_sj_crossing,
                sj_low_boundary=sj_low_boundary,
                sj_high_boundary=sj_high_boundary,
            )
        return arrays

    def _count(self, routed: _Routed) -> None:
        stats = self.stats
        stats.deterministic_unambig_units += int(routed.n_deterministic)
        stats.em_routed_unambig_units += int(routed.n_em_unambig)
        stats.em_routed_ambig_same_strand_units += int(routed.n_em_ambig_same_strand)
        stats.em_routed_ambig_opp_strand_units += int(routed.n_em_ambig_opp_strand)
        stats.n_gated_out += int(routed.n_gated)
        stats.em_routed_multimapper_units += int(routed.n_em_multimapper)
        # A multimapper whose only places left after the gates are intergenic counts as intergenic
        stats.n_intergenic += len(routed.intergenic_fids)

    def _annotate(self, routed: _Routed) -> None:
        """Annotate the fragments routing decides on its own: a deterministic-unique fragment is its transcript's with
        posterior 1; a fragment a hard gate dropped reaches no pool, unresolved with its class; and a multimapper left
        only intergenic places is intergenic gDNA."""
        annotations = self.annotations
        tids = np.asarray(routed.det_tids, dtype=np.int32)
        n = len(tids)
        if n > 0:
            is_synthetic = self.index.t_df["is_synthetic"].values.astype(bool)
            annotations.add_batch(
                frag_ids=np.asarray(routed.det_fids, dtype=np.int64),
                best_tids=tids,
                best_gids=np.asarray(self.index.t_to_g_arr[tids], dtype=np.int32),
                tx_flags=winner_flags(is_synthetic[tids]),
                posteriors=np.ones(n, dtype=np.float32),
                frag_classes=np.full(n, FRAG_UNAMBIG, dtype=np.int8),
                n_candidates=np.ones(n, dtype=np.int32),
            )
        n = len(routed.gated_fids)
        if n > 0:
            annotations.add_batch(
                frag_ids=np.asarray(routed.gated_fids, dtype=np.int64),
                best_tids=np.full(n, -1, dtype=np.int32),
                best_gids=np.full(n, -1, dtype=np.int32),
                tx_flags=np.full(n, AF_UNRESOLVED, dtype=np.uint8),
                posteriors=np.zeros(n, dtype=np.float32),
                frag_classes=np.asarray(routed.gated_classes, dtype=np.int8),
                n_candidates=np.zeros(n, dtype=np.int32),
            )
        n = len(routed.intergenic_fids)
        if n > 0:
            annotations.add_batch(
                frag_ids=np.asarray(routed.intergenic_fids, dtype=np.int64),
                best_tids=np.full(n, -1, dtype=np.int32),
                best_gids=np.full(n, -1, dtype=np.int32),
                tx_flags=np.full(n, AF_GDNA_INTERGENIC, dtype=np.uint8),
                posteriors=np.ones(n, dtype=np.float32),
                frag_classes=np.full(n, FRAG_MULTIMAPPER, dtype=np.int8),
                n_candidates=np.zeros(n, dtype=np.int32),
            )
