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

import gc
import logging
import math
import os
from dataclasses import dataclass, replace as _replace
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .buffer import FragmentBuffer, _FinalizedChunk
from .config import (
    EMConfig,
    PipelineConfig,
    BamScanConfig,
    FragmentScoringConfig,
    TranscriptGeometry,
)
from .estimator import AbundanceEstimator
from .frag_length_model import FragmentLengthModels
from .index import TranscriptIndex
from .locus import build_loci, compute_gdna_locus_gammas
from .native import BamScanner as _NativeBamScanner
from .native import detect_sj_strand_tag as _native_detect_sj_tag
from .partition import partition_and_free
from .scan import FragmentRouter
from .scoring import FragmentScorer
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from .annotate import AnnotationTable
    from .calibration import GDNACalibration
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
    calibration: "GDNACalibration" = None


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
        ("intergenic", strand_models.intergenic),
    ]:
        obs = strand_dict.get(f"{prefix}_obs", [])
        truth = strand_dict.get(f"{prefix}_truth", [])
        if len(obs) > 0:
            model.observe_batch(
                np.asarray(obs, dtype=np.int8),
                np.asarray(truth, dtype=np.int8),
            )


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

    intergenic_lengths = fraglen_dict.get("intergenic_lengths", [])
    if len(intergenic_lengths) > 0:
        frag_length_models.observe_intergenic_batch(
            np.asarray(intergenic_lengths, dtype=np.intp),
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
        "n_frag_length_intergenic",
        # Multimapper
        "n_multimapper_groups",
        "n_multimapper_alignments",
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
    pd.DataFrame | None,
    pd.DataFrame | None,
]:
    """Single-pass C++ BAM scan: resolve, train models, buffer — all in one pass.

    The entire BAM parsing, fragment construction, overlap resolution,
    model training and columnar buffering happens in C++ via htslib.
    When the index contains a region partition, region-level fragment
    evidence is accumulated during the same scan pass.

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
    tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer,
          pd.DataFrame | None, pd.DataFrame | None]
        The last two elements are *region_counts* and *fl_table* when the
        index contains a region partition, or ``None`` otherwise.
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

    # Provide nRNA status so FL training excludes nRNA candidates
    nrna_arr = index.t_df["is_synthetic_nrna"].values.astype("uint8")
    resolve_ctx.set_nrna_status(nrna_arr.tolist())

    # Build region partition index for gDNA calibration (if available)
    has_regions = index.region_df is not None and len(index.region_df) > 0
    if has_regions:
        rdf = index.region_df
        resolve_ctx.build_region_index(
            rdf["ref"].values.tolist(),
            rdf["start"].values.tolist(),
            rdf["end"].values.tolist(),
            rdf["region_id"].values.tolist(),
        )
        logger.debug(f"Built region index with {len(rdf)} regions for C++ scanner")

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

    # Extract region evidence from the scan result (if accumulated)
    region_counts: pd.DataFrame | None = None
    fl_table: pd.DataFrame | None = None
    if "region_evidence" in result:
        re = result["region_evidence"]
        n_regions = re["n_regions"]
        raw_counts = np.asarray(re["counts"], dtype=np.float64).reshape(n_regions, 4)
        region_counts = pd.DataFrame(
            {
                "region_id": np.arange(n_regions, dtype=np.int32),
                "n_unspliced_pos": raw_counts[:, 0].astype(np.float32),
                "n_unspliced_neg": raw_counts[:, 1].astype(np.float32),
                "n_spliced_pos": raw_counts[:, 2].astype(np.float32),
                "n_spliced_neg": raw_counts[:, 3].astype(np.float32),
            }
        )
        fl_region_ids = np.asarray(re["fl_region_ids"], dtype=np.int32)
        fl_frag_lens = np.asarray(re["fl_frag_lens"], dtype=np.int32)
        fl_table = pd.DataFrame({"region_id": fl_region_ids, "frag_len": fl_frag_lens})
        logger.info(f"[DONE] Region evidence: {n_regions} regions, {len(fl_table)} FL observations")

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

    return stats, strand_models, frag_length_models, buffer, region_counts, fl_table


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
        is_synthetic_nrna=index.t_df["is_synthetic_nrna"].values,
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

    # Free scanner accumulators + buffer: scan is done
    del builder, ctx
    buffer.release()
    return em_data


def _compute_priors(
    estimator: AbundanceEstimator,
    loci: list,
    index: TranscriptIndex,
    calibration: "GDNACalibration" = None,
) -> np.ndarray:
    """Compute calibration locus gammas.

    Returns the per-locus locus_gammas array (float64[n_loci] ∈ [0, 1]).
    """
    # Assign locus_id to every transcript (needed by nrna_frac prior cascade)
    for locus in loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.locus_id

    # Aggregate calibration region posteriors per locus
    locus_gammas = compute_gdna_locus_gammas(loci, index, calibration)

    return locus_gammas


def _run_locus_em_partitioned(
    estimator: AbundanceEstimator,
    partitions: dict,
    loci: list,
    index: TranscriptIndex,
    locus_gammas: np.ndarray,
    em_config: EMConfig,
    *,
    emit_locus_stats: bool = False,
) -> None:
    """Run batch locus EM from partitioned data with incremental freeing."""
    t_to_g = index.t_to_g_arr
    n_threads = em_config.n_threads or os.cpu_count() or 1

    def _build_locus_meta(locus, *, mrna, gdna, locus_gamma):
        gene_set = {int(t_to_g[int(t_idx)]) for t_idx in locus.transcript_indices}
        return {
            "locus_id": locus.locus_id,
            "locus_span_bp": locus.gdna_span,
            "n_transcripts": len(locus.transcript_indices),
            "n_genes": len(gene_set),
            "n_em_fragments": len(locus.unit_indices),
            "mrna": float(mrna),
            "gdna": float(gdna),
            "gdna_init": float(locus_gamma),
        }

    def _call_batch_em(parts, batch_loci, batch_gammas, batch_spans):
        """Pack tuples, call C++, record results."""
        partition_tuples = [
            (
                p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
                p.tx_starts, p.tx_ends, p.count_cols,
                p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
                p.locus_t_indices, p.locus_count_cols,
            )
            for p in parts
        ]
        locus_t_lists = [l.transcript_indices for l in batch_loci]

        return estimator.run_batch_locus_em_partitioned(
            partition_tuples,
            locus_t_lists,
            batch_gammas,
            batch_spans,
            index,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
            emit_locus_stats=emit_locus_stats,
        )

    # Classify mega vs normal
    total_work = sum(
        len(l.transcript_indices) * partitions[l.locus_id].n_units
        for l in loci
    )
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [l for l in loci
         if len(l.transcript_indices) * partitions[l.locus_id].n_units >= fair_share],
        key=lambda l: len(l.transcript_indices) * partitions[l.locus_id].n_units,
        reverse=True,
    )
    mega_ids = {l.locus_id for l in mega_loci}

    total_gdna_em = 0.0

    # Phase A: Mega-loci (one at a time, free after each)
    for locus in mega_loci:
        part = partitions.pop(locus.locus_id)
        gdna_em, mrna_arr, gdna_arr = _call_batch_em(
            [part], [locus],
            np.array([locus_gammas[locus.locus_id]], dtype=np.float64),
            np.array([locus.gdna_span], dtype=np.int64),
        )
        total_gdna_em += gdna_em
        estimator.locus_results.append(
            _build_locus_meta(
                locus,
                mrna=mrna_arr[0],
                gdna=gdna_arr[0],
                locus_gamma=locus_gammas[locus.locus_id],
            )
        )
        del part
        gc.collect()
        logger.debug(
            f"[MEGA] Locus {locus.locus_id}: "
            f"{len(locus.transcript_indices)} transcripts, "
            f"{len(locus.unit_indices)} units"
        )

    # Phase B: Normal loci (one batched call)
    normal_loci = [l for l in loci if l.locus_id not in mega_ids]
    if normal_loci:
        normal_parts = [partitions[l.locus_id] for l in normal_loci]
        normal_gammas = np.array(
            [locus_gammas[l.locus_id] for l in normal_loci], dtype=np.float64
        )
        normal_spans = np.array(
            [l.gdna_span for l in normal_loci], dtype=np.int64
        )
        gdna_em, mrna_arr, gdna_arr = _call_batch_em(
            normal_parts, normal_loci, normal_gammas, normal_spans,
        )
        total_gdna_em += gdna_em
        for i, locus in enumerate(normal_loci):
            estimator.locus_results.append(
                _build_locus_meta(
                    locus,
                    mrna=mrna_arr[i],
                    gdna=gdna_arr[i],
                    locus_gamma=locus_gammas[locus.locus_id],
                )
            )

    del partitions
    gc.collect()

    estimator._gdna_em_total = total_gdna_em

    n_total_units = sum(len(l.unit_indices) for l in loci)
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
    *,
    em_config: EMConfig | None = None,
    scoring: FragmentScoringConfig | None = None,
    log_every: int = 1_000_000,
    annotations: "AnnotationTable | None" = None,
    calibration: "GDNACalibration" = None,
    emit_locus_stats: bool = False,
) -> AbundanceEstimator:
    """Quantify transcripts from buffered fragments via locus-level EM.

    Parameters
    ----------
    buffer : FragmentBuffer
    index : TranscriptIndex
    strand_models : StrandModels
    frag_length_models : FragmentLengthModels
    stats : PipelineStats
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
    # Calibration always provides a gDNA FL model (EM-based when ESS
    # is sufficient, intergenic fallback otherwise).
    frag_length_models.gdna_model = calibration.gdna_fl_model
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

        locus_gammas = _compute_priors(
            estimator,
            loci,
            index,
            calibration=calibration,
        )

        # Phase 4 (NEW): Array-by-array scatter + incremental free
        n_em_units = em_data.n_units
        partitions = partition_and_free(em_data, loci)
        del em_data
        gc.collect()

        # Phase 5 (NEW): Streaming locus EM with incremental partition freeing
        _run_locus_em_partitioned(
            estimator,
            partitions,
            loci,
            index,
            locus_gammas,
            em_config,
            emit_locus_stats=emit_locus_stats,
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")
        del em_data

    # Phase 6: Cleanup
    gc.collect()

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
    stats, strand_models, frag_length_models, buffer, region_counts, fl_table = scan_and_buffer(
        bam_path, index, scan
    )

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
    frag_length_models.build_scoring_models()
    frag_length_models.finalize()

    # -- gDNA calibration (post-model, pre-scoring) --
    if region_counts is None or fl_table is None or index.region_df is None:
        raise RuntimeError(
            "gDNA calibration requires region data. Rebuild index with "
            "'rigel index' to generate regions.feather."
        )
    from .calibration import calibrate_gdna

    cal_cfg = config.calibration
    calibration = calibrate_gdna(
        region_counts,
        fl_table,
        index.region_df,
        strand_models.strand_specificity,
        max_iterations=cal_cfg.max_iterations,
        convergence_tol=cal_cfg.convergence_tol,
        density_percentile=cal_cfg.density_percentile,
        min_gdna_regions=cal_cfg.min_gdna_regions,
        min_fl_ess=cal_cfg.min_fl_ess,
        intergenic_fl_model=frag_length_models.intergenic,
    )
    logger.info(
        f"[CAL] gDNA calibration: π={calibration.mixing_proportion:.3f}, "
        f"κ_strand={calibration.kappa_strand:.1f}, λ_G={calibration.gdna_density_global:.2e}, "
        f"converged={calibration.converged} ({calibration.n_iterations} iters)"
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
