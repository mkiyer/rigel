"""
rigel.pipeline — Single-pass BAM quantification pipeline with locus-level EM.

Architecture
------------
The pipeline reads the BAM file **once** via the C++ native BAM scanner
(``_bam_impl``) and processes fragments in two logical stages:

**BAM Scan** (``scan_and_buffer``): Parse all fragments in C++ via htslib,
resolve each against the reference index, train strand and fragment-length
models from unique-mapper fragments, and buffer all resolved fragments into
a memory-efficient columnar buffer (``FragmentBuffer``).

**Quantification** (``quant_from_buffer``): Iterate the buffer once
to assign unambig counts, build CSR EM data (via ``scan.FragmentRouter``),
construct loci (via ``locus.build_loci``), compute Empirical Bayes gDNA
priors, and run per-locus EM with ``2*n_t + 1`` components.

Scoring functions live in ``scoring.py``.  Locus construction and EM
initialization live in ``locus.py``.  The CSR builder lives in
``scan.py``.  This module is a thin orchestrator.
"""

from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass, replace as _replace
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .annotate import (
    AF_GDNA_EM,
    AF_UNRESOLVED,
    winner_flags,
)
from .buffer import FragmentBuffer, _FinalizedChunk
from .config import (
    EMConfig,
    PipelineConfig,
    BamScanConfig,
    FragmentScoringConfig,
    TranscriptGeometry,
)
from .estimator import AbundanceEstimator
from .frag_length_model import FragmentLengthModel, FragmentLengthModels
from .index import TranscriptIndex
from .locus import build_loci, compute_locus_priors_from_partitions
from .native import BamScanner as _NativeBamScanner
from .native import detect_sj_strand_tag as _native_detect_sj_tag
from .locus_partition import partition_and_free
from .scan import FragmentRouter
from .scoring import FragmentScorer
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from .annotate import AnnotationTable
    from .calibration import CalibrationResult
    from .scored_fragments import ScoredFragments

logger = logging.getLogger(__name__)

# Padding and minimum capacity for the annotation table.
_ANNOTATION_TABLE_PADDING = 1024
_ANNOTATION_TABLE_MIN_CAPACITY = 4096

#: Fallback mean fragment length when no observations are available.
_DEFAULT_MEAN_FRAG: float = 200.0

# ---------------------------------------------------------------------------
# Pipeline result
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """Complete output of the pipeline."""

    stats: PipelineStats
    strand_models: StrandModels
    frag_length_models: FragmentLengthModels
    estimator: AbundanceEstimator
    pipeline_config: "PipelineConfig" = None
    calibration: "CalibrationResult" = None


def _sj_tag_to_spec(sj_strand_tag) -> str:
    """Convert BamScanConfig.sj_strand_tag to the string spec for BamScanner."""
    if isinstance(sj_strand_tag, str):
        return sj_strand_tag if sj_strand_tag else "none"
    if isinstance(sj_strand_tag, (list, tuple)):
        return ",".join(sj_strand_tag) if sj_strand_tag else "none"
    return "none"


def _replay_strand_observations(
    strand_dict: dict,
    strand_models: StrandModels,
) -> None:
    """Replay C++ strand observation arrays into Python StrandModels."""
    for prefix, model in [
        ("exonic_spliced", strand_models.exonic_spliced),
        ("exonic", strand_models.exonic),
    ]:
        obs = strand_dict.get(f"{prefix}_obs", [])
        truth = strand_dict.get(f"{prefix}_truth", [])
        if len(obs) > 0:
            model.observe_batch(obs, truth)


def _replay_fraglen_observations(
    fraglen_dict: dict,
    frag_length_models: FragmentLengthModels,
) -> None:
    """Replay C++ fragment-length observation arrays into Python models."""
    lengths = fraglen_dict.get("lengths", [])
    splice_types = fraglen_dict.get("splice_types", [])
    if len(lengths) > 0:
        frag_length_models.observe_batch(
            np.asarray(lengths, dtype=np.intp),
            np.asarray(splice_types, dtype=np.intp),
        )


def _apply_scan_stats(stats: PipelineStats, stats_dict: dict) -> None:
    """Apply C++ scan statistics to PipelineStats.

    Copies all matching keys from the C++ stats dict to the dataclass.
    """
    for key in (
        # BAM-level
        "total",
        "qc_fail",
        "unmapped",
        "secondary",
        "supplementary",
        "duplicate",
        "n_read_names",
        "unique",
        "multimapping",
        "proper_pair",
        "improper_pair",
        "mate_unmapped",
        # Fragment-level
        "n_fragments",
        "n_chimeric",
        "n_chimeric_trans",
        "n_chimeric_cis_strand_same",
        "n_chimeric_cis_strand_diff",
        "n_intergenic_unspliced",
        "n_intergenic_spliced",
        "n_with_exon",
        "n_with_annotated_sj",
        "n_with_unannotated_sj",
        "n_same_strand",
        "n_ambig_strand",
        # Strand model training
        "n_strand_trained",
        "n_strand_skipped_no_sj",
        "n_strand_skipped_ambig_strand",
        "n_strand_skipped_ambiguous",
        # Fragment length model training
        "n_frag_length_unambiguous",
        "n_frag_length_ambiguous",
        # Multimapper
        "n_multimapper_groups",
        "n_multimapper_alignments",
        # Splice-artifact blacklist
        "n_sj_observed",
        "n_sj_blacklisted",
    ):
        setattr(stats, key, stats_dict.get(key, 0))


def scan_and_buffer(
    bam_path: str,
    index: TranscriptIndex,
    scan: "BamScanConfig",
) -> tuple[
    PipelineStats,
    StrandModels,
    FragmentLengthModels,
    FragmentBuffer,
]:
    """Single-pass C++ BAM scan: resolve, train models, buffer — all in one pass.

    The entire BAM parsing, fragment construction, overlap resolution,
    model training and columnar buffering happens in C++ via htslib.

    Parameters
    ----------
    bam_path : str
        Path to the name-sorted / collated BAM file.
    index : TranscriptIndex
        The loaded reference index.
    scan : BamScanConfig
        BAM scanning and buffering configuration.

    Returns
    -------
    tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]
    """
    stats = PipelineStats()
    strand_models = StrandModels()
    frag_length_models = FragmentLengthModels(max_size=scan.max_frag_length)
    buffer = FragmentBuffer(
        t_strand_arr=index.t_to_strand_arr,
        chunk_size=scan.chunk_size,
        max_memory_bytes=scan.max_memory_bytes,
        spill_dir=scan.spill_dir,
    )

    logger.info("[START] Native C++ BAM scan → resolve + train + buffer")

    # Provide gene strand info for exonic fallback strand model training
    resolve_ctx = index.resolver
    resolve_ctx.set_gene_strands(index.g_to_strand_arr.tolist())
    resolve_ctx.set_transcript_strands(index.t_to_strand_arr.tolist())

    # Provide nRNA status so FL training excludes synthetic nRNA candidates
    # (whose fragment lengths represent genomic spans, not real fragment sizes).
    # Only synthetic nRNAs should be excluded; annotated single-exon transcripts
    # have legitimate fragment lengths and must contribute to FL training.
    nrna_arr = index.t_df["is_synthetic"].values.astype("uint8")
    resolve_ctx.set_nrna_status(nrna_arr.tolist())

    # nRNA parent-index wiring (set_nrna_parent_index) is performed by
    # TranscriptIndex.load() at index-load time, so no need to repeat
    # it here.

    # Create the native scanner
    sj_spec = _sj_tag_to_spec(scan.sj_strand_tag)
    scanner = _NativeBamScanner(
        resolve_ctx,
        sj_spec,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
    )

    # Streaming chunk callback — receives zero-copy dict from C++
    def _on_chunk(raw: dict) -> None:
        chunk = _FinalizedChunk.from_raw(raw)
        buffer.inject_chunk(chunk)

    # Execute the full BAM scan in C++ with streaming chunk output
    n_scan = scan.n_scan_threads
    if n_scan == 0:
        n_scan = os.cpu_count() or 1
    result = scanner.scan(
        bam_path,
        chunk_callback=_on_chunk,
        t_strand_arr=index.t_to_strand_arr.tolist(),
        chunk_size=scan.chunk_size,
        n_workers=n_scan,
        n_decomp_threads=scan.n_decomp_threads,
    )

    # Replay stats
    _apply_scan_stats(stats, result["stats"])

    # Replay strand observations into Python models
    _replay_strand_observations(result["strand_observations"], strand_models)

    # Replay fragment-length observations into Python models
    _replay_fraglen_observations(
        result["frag_length_observations"],
        frag_length_models,
    )

    logger.info(
        f"[DONE] Native scan: {stats.n_fragments:,} fragments → "
        f"{stats.n_same_strand:,} same-strand, "
        f"{stats.n_ambig_strand:,} ambig-strand, "
        f"{stats.n_intergenic:,} intergenic, "
        f"{buffer.total_fragments:,} buffered "
        f"({buffer.memory_bytes / 1024**2:.1f} MB in memory, "
        f"{buffer.n_spilled} chunks spilled), "
        f"{stats.n_multimapper_groups:,} multimapper molecules"
    )
    strand_models.log_summary()
    frag_length_models.log_summary()

    return stats, strand_models, frag_length_models, buffer


# ---------------------------------------------------------------------------
# Quantify from buffer (locus-level EM, no global EM)
# ---------------------------------------------------------------------------


def _setup_geometry_and_estimator(
    index: TranscriptIndex,
    frag_length_models: FragmentLengthModels,
    em_config: EMConfig,
) -> tuple["TranscriptGeometry", AbundanceEstimator]:
    """Compute transcript geometry and create the AbundanceEstimator."""
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    if frag_length_models.rna_model.n_observations > 0:
        effective_lengths = frag_length_models.rna_model.compute_all_transcript_eff_lens(
            exonic_lengths.astype(np.int64),
        )
    else:
        effective_lengths = np.maximum(exonic_lengths - _DEFAULT_MEAN_FRAG + 1.0, 1.0)

    transcript_spans = (index.t_df["end"].values - index.t_df["start"].values).astype(np.float64)

    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        transcript_spans=transcript_spans,
    )

    estimator = AbundanceEstimator(
        index.num_transcripts,
        em_config=em_config,
        geometry=geometry,
        is_nrna=index.t_df["is_nrna"].values,
        is_synthetic=index.t_df["is_synthetic"].values,
    )
    return geometry, estimator


def _score_fragments(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    frag_length_models: FragmentLengthModels,
    stats: PipelineStats,
    estimator: AbundanceEstimator,
    scoring: "FragmentScoringConfig",
    log_every: int,
    annotations,
) -> "ScoredFragments":
    """Build FragmentScorer, scan buffer, and return ScoredFragments."""
    ctx = FragmentScorer.from_models(
        strand_models,
        frag_length_models,
        index,
        estimator,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
        pruning_min_posterior=scoring.pruning_min_posterior,
    )
    builder = FragmentRouter(
        ctx,
        estimator,
        stats,
        index,
        strand_models,
        annotations=annotations,
    )
    em_data = builder.scan(buffer, log_every)

    # Scanner accumulators no longer needed; buffer was consumed during scan
    del builder, ctx
    return em_data


def _assign_locus_ids(estimator: AbundanceEstimator, loci: list) -> None:
    """Stamp ``locus_id`` onto every transcript on the estimator.

    Required by the nRNA-fraction prior cascade in the C++ EM.
    """
    for locus in loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.locus_id


def _populate_em_annotations(
    batch_parts,
    out_winner_tid,
    out_winner_post,
    out_n_candidates,
    annotations,
    index,
):
    """Populate AnnotationTable from EM assignment output and partition metadata."""
    if out_winner_tid is None:
        return

    t_to_g = index.t_to_g_arr
    is_nrna_arr = index.t_df["is_nrna"].values.astype(bool)
    is_synth_arr = index.t_df["is_synthetic"].values.astype(bool)

    # Concatenate per-locus metadata in same order as C++ output
    frag_ids = np.concatenate([p.frag_ids for p in batch_parts])
    frag_class = np.concatenate([p.frag_class for p in batch_parts])
    splice_type_arr = np.concatenate([p.splice_type for p in batch_parts])

    n = len(frag_ids)
    best_tid = np.asarray(out_winner_tid, dtype=np.int32)
    posteriors = np.asarray(out_winner_post, dtype=np.float32)
    n_cand = np.asarray(out_n_candidates, dtype=np.int16)

    # Gene index: map valid transcript winners to gene
    valid_t = best_tid >= 0
    best_gid = np.full(n, -1, dtype=np.int32)
    best_gid[valid_t] = t_to_g[best_tid[valid_t]]

    # ZF assignment flags bitfield
    tx_flags = np.full(n, AF_UNRESOLVED, dtype=np.uint8)
    tx_flags[best_tid == -2] = AF_GDNA_EM
    if valid_t.any():
        winner_tids = best_tid[valid_t]
        tx_flags[valid_t] = winner_flags(
            is_nrna_arr[winner_tids],
            is_synth_arr[winner_tids],
        )

    # Clean sentinel: -2 (gDNA) -> -1
    best_tid_clean = best_tid.copy()
    best_tid_clean[best_tid == -2] = -1

    annotations.add_batch(
        frag_ids=frag_ids,
        best_tids=best_tid_clean,
        best_gids=best_gid,
        tx_flags=tx_flags,
        posteriors=posteriors,
        frag_classes=frag_class.view(np.int8),
        n_candidates=n_cand,
        splice_types=splice_type_arr,
    )


def _run_locus_em_partitioned(
    estimator: AbundanceEstimator,
    partitions: dict,
    loci: list,
    index: TranscriptIndex,
    alpha_gdna: np.ndarray,
    alpha_rna: np.ndarray,
    em_config: EMConfig,
    *,
    emit_locus_stats: bool = False,
    annotations: "AnnotationTable | None" = None,
) -> None:
    """Run batch locus EM from partitioned data with incremental freeing."""
    t_to_g = index.t_to_g_arr
    # ``is_synthetic_g`` marks synthetic (gene-neutral) gene rows so they can
    # be excluded from the user-facing ``n_genes`` per locus.
    if "is_synthetic" in index.g_df.columns:
        is_synthetic_g = index.g_df["is_synthetic"].to_numpy()
    else:
        is_synthetic_g = np.zeros(len(index.g_df), dtype=bool)
    n_threads = em_config.n_threads or os.cpu_count() or 1
    emit_assignments = annotations is not None

    def _build_locus_meta(locus, *, rna_total, gdna, alpha_g, alpha_r):
        gene_set = {
            int(t_to_g[int(t_idx)])
            for t_idx in locus.transcript_indices
            if not is_synthetic_g[int(t_to_g[int(t_idx)])]
        }
        return {
            "locus_id": locus.locus_id,
            "locus_span_bp": locus.gdna_span,
            "n_transcripts": len(locus.transcript_indices),
            "n_genes": len(gene_set),
            "n_em_fragments": len(locus.unit_indices),
            # ``rna_total`` = sum of posterior mass assigned to every
            # non-gDNA component in the locus. The C++ EM treats annotated
            # mRNA and synthetic nRNA transcripts identically, so this
            # value is NOT mRNA alone — it is (annotated mRNA + synthetic
            # nRNA), i.e. total RNA. The annotated vs synthetic split is
            # done downstream in ``AbundanceEstimator.get_loci_df`` using
            # per-transcript ``is_synthetic`` flags from the index.
            "rna_total": float(rna_total),
            "gdna": float(gdna),
            "alpha_gdna": float(alpha_g),
            "alpha_rna": float(alpha_r),
        }

    def _call_batch_em(parts, batch_loci, batch_alpha_gdna, batch_alpha_rna):
        """Pack tuples, call C++, record results."""
        partition_tuples = [
            (
                p.offsets,
                p.t_indices,
                p.log_liks,
                p.coverage_weights,
                p.tx_starts,
                p.tx_ends,
                p.count_cols,
                p.is_spliced,
                p.gdna_log_liks,
                p.locus_t_indices,
                p.locus_count_cols,
            )
            for p in parts
        ]
        locus_t_lists = [loc.transcript_indices for loc in batch_loci]

        return estimator.run_batch_locus_em_partitioned(
            partition_tuples,
            locus_t_lists,
            batch_alpha_gdna,
            batch_alpha_rna,
            index,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
            emit_locus_stats=emit_locus_stats,
            emit_assignments=emit_assignments,
        )

    # Classify mega vs normal
    locus_work = {
        loc.locus_id: len(loc.transcript_indices) * partitions[loc.locus_id].n_units for loc in loci
    }
    total_work = sum(locus_work.values())
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [loc for loc in loci if locus_work[loc.locus_id] >= fair_share],
        key=lambda loc: locus_work[loc.locus_id],
        reverse=True,
    )
    mega_ids = {loc.locus_id for loc in mega_loci}

    total_gdna_em = 0.0

    # Phase A: Mega-loci (one at a time, free after each)
    for locus in mega_loci:
        part = partitions.pop(locus.locus_id)
        lid = locus.locus_id
        em_result = _call_batch_em(
            [part],
            [locus],
            np.array([alpha_gdna[lid]], dtype=np.float64),
            np.array([alpha_rna[lid]], dtype=np.float64),
        )
        gdna_em, rna_arr, gdna_arr = em_result[0], em_result[1], em_result[2]
        total_gdna_em += gdna_em
        if annotations is not None and len(em_result) > 3:
            _populate_em_annotations(
                [part],
                em_result[3],
                em_result[4],
                em_result[5],
                annotations,
                index,
            )
        estimator.locus_results.append(
            _build_locus_meta(
                locus,
                rna_total=rna_arr[0],
                gdna=gdna_arr[0],
                alpha_g=alpha_gdna[lid],
                alpha_r=alpha_rna[lid],
            )
        )
        del part
        logger.debug(
            f"[MEGA] Locus {locus.locus_id}: "
            f"{len(locus.transcript_indices)} transcripts, "
            f"{len(locus.unit_indices)} units"
        )

    # Phase B: Normal loci (one batched call)
    normal_loci = [loc for loc in loci if loc.locus_id not in mega_ids]
    if normal_loci:
        # Pop partitions out of the dict so the only references during
        # the batched C++ call live in ``normal_parts``.  After the call
        # completes and annotations are written we drop ``normal_parts``
        # to release per-locus arrays before EM downstream phases run.
        normal_parts = [partitions.pop(loc.locus_id) for loc in normal_loci]
        normal_ag = np.array([alpha_gdna[loc.locus_id] for loc in normal_loci], dtype=np.float64)
        normal_ar = np.array([alpha_rna[loc.locus_id] for loc in normal_loci], dtype=np.float64)
        em_result = _call_batch_em(
            normal_parts,
            normal_loci,
            normal_ag,
            normal_ar,
        )
        gdna_em, rna_arr, gdna_arr = em_result[0], em_result[1], em_result[2]
        total_gdna_em += gdna_em
        if annotations is not None and len(em_result) > 3:
            _populate_em_annotations(
                normal_parts,
                em_result[3],
                em_result[4],
                em_result[5],
                annotations,
                index,
            )
        for i, locus in enumerate(normal_loci):
            estimator.locus_results.append(
                _build_locus_meta(
                    locus,
                    rna_total=rna_arr[i],
                    gdna=gdna_arr[i],
                    alpha_g=normal_ag[i],
                    alpha_r=normal_ar[i],
                )
            )
        # Release per-locus partition arrays before downstream phases.
        del normal_parts, em_result

    del partitions

    estimator._gdna_em_total = total_gdna_em

    n_total_units = sum(len(loc.unit_indices) for loc in loci)
    logger.info(
        f"[DONE] Per-locus EM (partitioned): {len(loci)} loci "
        f"({len(mega_loci)} mega), "
        f"{n_total_units:,} ambiguous fragments, "
        f"gDNA EM={total_gdna_em:.0f}"
    )


def quant_from_buffer(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    frag_length_models: FragmentLengthModels,
    stats: PipelineStats,
    calibration: "CalibrationResult",
    *,
    em_config: EMConfig | None = None,
    scoring: FragmentScoringConfig | None = None,
    log_every: int = 1_000_000,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
    fl_prior_ess: float | None = None,
) -> AbundanceEstimator:
    """Quantify transcripts from buffered fragments via locus-level EM.

    Parameters
    ----------
    buffer : FragmentBuffer
    index : TranscriptIndex
    strand_models : StrandModels
    frag_length_models : FragmentLengthModels
    stats : PipelineStats
    calibration : CalibrationResult
        Required.  Provides the gDNA fragment-length model used for
        scoring and the per-locus gDNA Dirichlet prior parameters.
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    scoring : FragmentScoringConfig or None
        Scoring penalty configuration.  Defaults to ``FragmentScoringConfig()``.
    log_every : int
        Log progress every N fragments (default 1M).
    annotations : AnnotationTable or None
        If provided, record per-fragment assignment annotations.

    Returns
    -------
    AbundanceEstimator
    """
    if calibration is None or calibration.gdna_fl_model is None:
        raise ValueError(
            "quant_from_buffer() requires a calibrated CalibrationResult "
            "with a populated gdna_fl_model; got "
            f"{calibration!r}.  Run the calibration stage before "
            "locus-level quantification."
        )
    if em_config is None:
        em_config = EMConfig()
    if scoring is None:
        scoring = FragmentScoringConfig()

    # Phase 1: Geometry + estimator
    geometry, estimator = _setup_geometry_and_estimator(
        index,
        frag_length_models,
        em_config,
    )

    logger.info(
        f"[START] Quantifying {buffer.total_fragments:,} buffered fragments "
        f"(locus-level EM: mRNA/nRNA/gDNA)"
    )

    # Apply calibrated gDNA fragment-length model for scoring.
    # Copy to avoid aliasing (calibration may return the intergenic
    # model object directly).  Re-finalize with global Dirichlet prior
    # so gDNA and RNA FL models share the same prior baseline.
    cal_gdna = calibration.gdna_fl_model
    gdna_copy = FragmentLengthModel(max_size=cal_gdna.max_size)
    gdna_copy.counts = cal_gdna.counts.copy()
    gdna_copy._total_weight = cal_gdna._total_weight
    gdna_copy.finalize(prior_counts=frag_length_models.global_model.counts, prior_ess=fl_prior_ess)
    frag_length_models.gdna_model = gdna_copy
    logger.info(
        f"[CAL-FL] gDNA FL model: "
        f"mean={calibration.gdna_fl_model.mean:.1f}, "
        f"obs={calibration.gdna_fl_model.n_observations}"
    )

    # Phase 2: Score fragments
    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        frag_length_models,
        stats,
        estimator,
        scoring,
        log_every,
        annotations,
    )

    logger.info(
        f"[DONE] Scan: {stats.deterministic_unambig_units:,} unique, "
        f"{em_data.n_units:,} ambiguous units "
        f"({em_data.n_candidates:,} RNA candidates)"
    )

    # Phase 3: Locus construction, priors, and EM
    if em_data.n_units > 0:
        loci = build_loci(em_data, index)

        if loci:
            max_locus_t = max(len(loc.transcript_indices) for loc in loci)
            max_locus_u = max(len(loc.unit_indices) for loc in loci)
        else:
            max_locus_t = max_locus_u = 0

        logger.info(
            f"[LOCI] {len(loci)} loci from {em_data.n_units:,} units "
            f"(largest: {max_locus_t} transcripts, {max_locus_u} units)"
        )

        _assign_locus_ids(estimator, loci)

        # Phase 4 (NEW): Fused scatter into per-locus tuples
        # Phase 4 (NEW): Array-by-array scatter + incremental free
        partitions = partition_and_free(em_data, loci)
        del em_data

        # SRD v1 Pass 4: per-locus Dirichlet prior from per-fragment
        # gDNA posteriors aggregated over the partitioned data.  See
        # docs/calibration/srd_v1_implementation.md §2.6.
        alpha_gdna, alpha_rna = compute_locus_priors_from_partitions(
            partitions,
            loci,
            pi_pool=float(calibration.pi_pool),
        )

        # Phase 5 (NEW): Streaming locus EM with incremental partition freeing
        _run_locus_em_partitioned(
            estimator,
            partitions,
            loci,
            index,
            alpha_gdna,
            alpha_rna,
            em_config,
            emit_locus_stats=emit_locus_stats,
            annotations=annotations,
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")
        del em_data

    _gdna_em = estimator.gdna_em_count
    stats.n_gdna_em = 0 if (math.isnan(_gdna_em) or math.isinf(_gdna_em)) else int(_gdna_em)

    logger.info(
        f"[gDNA] total={estimator.gdna_em_count:.0f} "
        f"(EM={_gdna_em:.0f}), "
        f"contamination rate={estimator.gdna_contamination_rate:.2%}"
    )

    return estimator


# ---------------------------------------------------------------------------
# Full pipeline orchestration
# ---------------------------------------------------------------------------


def run_pipeline(
    bam_path,
    index: TranscriptIndex,
    config: PipelineConfig | None = None,
) -> PipelineResult:
    """Run the complete quantification pipeline with locus-level EM.

    Opens the BAM **once** via the C++ BAM scanner, scans all fragments,
    trains models, buffers resolved fragments, then quantifies via per-locus EM.

    Parameters
    ----------
    bam_path : str or Path
        Path to the name-sorted / collated BAM file.
    index : TranscriptIndex
        The loaded reference index.
    config : PipelineConfig or None
        Complete pipeline configuration.  If ``None``, uses defaults.

    Returns
    -------
    PipelineResult
    """
    if config is None:
        config = PipelineConfig()

    bam_path = str(bam_path)

    # -- Resolve sj_strand_tag "auto" → concrete tag(s) --
    # Use the C++ implementation (htslib) instead of the pysam-based
    # ``detect_sj_strand_tag`` because pysam's Cython tracing hooks
    # are corrupted when a cProfile profiler has been active during
    # a prior nanobind C++ extension call in the same process.
    # The native function returns a spec string ("XS", "ts", "XS,ts",
    # or "none") which ``_sj_tag_to_spec`` already handles.
    scan = config.scan
    if scan.sj_strand_tag == "auto":
        detected_spec = _native_detect_sj_tag(bam_path)
        scan = _replace(scan, sj_strand_tag=detected_spec)

    # -- Single BAM pass (C++ native scanner) --
    stats, strand_models, frag_length_models, buffer = scan_and_buffer(
        bam_path, index, scan
    )

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
    frag_length_models.build_scoring_models()
    frag_length_models.finalize(prior_ess=config.calibration.fl_prior_ess)

    # -- gDNA calibration (SRD v1: single buffer walk, geometric pool) --
    from .calibration import calibrate_gdna

    cal_cfg = config.calibration
    strand_ci_eps = strand_models.strand_specificity_ci_epsilon(confidence=0.99)
    logger.info(
        "[CAL] Strand trainer: n_spliced_obs=%d  ss_est=%.6f  ε_CI(99%%)=%.4g",
        strand_models.n_observations,
        strand_models.strand_specificity,
        strand_ci_eps,
    )

    calibration = calibrate_gdna(
        buffer,
        index,
        frag_length_models,
        strand_models.strand_specificity,
        read1_sense=bool(strand_models.read1_sense),
        exon_fit_tolerance_bp=cal_cfg.exon_fit_tolerance_bp,
        fl_prior_ess=cal_cfg.fl_prior_ess,
        max_iter=cal_cfg.max_iter,
        tol=cal_cfg.tol,
    )
    cal_summary = calibration.to_summary_dict()
    logger.info(
        f"[CAL] gDNA calibration (SRD v1): quality={cal_summary['gdna_fl_quality']}, "
        f"π_pool={cal_summary['pi_pool']:.4f}, "
        f"n_pool={cal_summary['n_pool']}, "
        f"SS={cal_summary['strand_specificity']:.3f}"
    )

    # -- Annotation table for second BAM pass (opt-in) --
    annotations = None
    if config.annotated_bam_path is not None:
        from .annotate import AnnotationTable

        annotations = AnnotationTable.create(
            capacity=max(
                buffer.total_fragments + _ANNOTATION_TABLE_PADDING,
                _ANNOTATION_TABLE_MIN_CAPACITY,
            )
        )

    try:
        estimator = quant_from_buffer(
            buffer,
            index,
            strand_models,
            frag_length_models,
            stats,
            em_config=config.em,
            scoring=config.scoring,
            log_every=scan.log_every,
            annotations=annotations,
            calibration=calibration,
            emit_locus_stats=config.emit_locus_stats,
            fl_prior_ess=config.calibration.fl_prior_ess,
        )
    finally:
        buffer.cleanup()

    # -- Second BAM pass: write annotated BAM (opt-in) --
    if config.annotated_bam_path is not None and annotations is not None:
        from .annotate import write_annotated_bam

        write_annotated_bam(
            bam_path,
            str(config.annotated_bam_path),
            annotations,
            index,
            skip_duplicates=scan.skip_duplicates,
            include_multimap=scan.include_multimap,
            sj_strand_tag=scan.sj_strand_tag,
            locus_id_per_transcript=estimator.locus_id_per_transcript,
        )

    return PipelineResult(
        stats=stats,
        strand_models=strand_models,
        frag_length_models=frag_length_models,
        estimator=estimator,
        pipeline_config=config,
        calibration=calibration,
    )
