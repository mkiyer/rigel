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
construct loci (via ``locus.build_multi_loci``), compute Empirical Bayes gDNA
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
from .frag_length_model import FragmentLengthModels
from .index import TranscriptIndex
from .locus import build_multi_loci
from .calibration.locus_prior import assemble_priors
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
    from .calibration._regional_exposure import (
        RegionalGdnaExposure,
        RegionalWeightApplicationStats,
    )
    from .calibration.scan_payload import CalibrationScanPayload
    from .scored_fragments import ScoredFragments

logger = logging.getLogger(__name__)

# Padding and minimum capacity for the annotation table.
_ANNOTATION_TABLE_PADDING = 1024
_ANNOTATION_TABLE_MIN_CAPACITY = 4096

#: Fallback mean fragment length when no observations are available.
_DEFAULT_MEAN_FRAG: float = 200.0


# ---------------------------------------------------------------------------
# Regional gDNA exposure: per-unit weight application
# ---------------------------------------------------------------------------


def _apply_unit_gdna_weights(
    em_data: "ScoredFragments",
    exposure: "RegionalGdnaExposure",
    index: TranscriptIndex,
) -> "RegionalWeightApplicationStats":
    """Apply ``log A_r`` to each unit's gDNA log-likelihood in place.

    Units with finite ``gdna_log_liks`` and a valid genomic midpoint get
    ``log A_r`` added at their midpoint.  Units whose candidate transcripts
    span multiple references are skipped (left at ``A=1``) and counted as
    ``n_units_cross_ref_skipped``.  In uniform-exposure mode this is a
    no-op besides bookkeeping.
    """
    from .calibration._regional_exposure import RegionalWeightApplicationStats

    n_units = int(em_data.n_units)
    if exposure.mode == "uniform" or n_units == 0:
        return RegionalWeightApplicationStats(n_units_seen=n_units)

    # Invariant guard (R6): every EM unit must have >=1 candidate.
    # ``np.*.reduceat`` silently misbehaves on empty groups.
    assert (np.diff(em_data.offsets) > 0).all(), "empty EM unit detected"

    finite = np.isfinite(em_data.gdna_log_liks)
    n_no_gdna = int((~finite).sum())

    int64_min = np.iinfo(np.int64).min
    has_mid = em_data.genomic_midpoint != int64_min
    n_missing_mid = int((finite & ~has_mid).sum())

    base_mask = finite & has_mid

    if not base_mask.any():
        return RegionalWeightApplicationStats(
            n_units_seen=n_units,
            n_units_no_gdna=n_no_gdna,
            n_units_missing_midpoint=n_missing_mid,
        )

    if index.t_to_ref_arr is None:
        raise RuntimeError("TranscriptIndex.t_to_ref_arr not populated")
    candidate_refs = index.t_to_ref_arr[em_data.t_indices]
    unit_starts = em_data.offsets[:-1]
    unit_min_ref = np.minimum.reduceat(candidate_refs, unit_starts)
    unit_max_ref = np.maximum.reduceat(candidate_refs, unit_starts)
    same_ref = unit_min_ref == unit_max_ref
    n_cross_ref = int((base_mask & ~same_ref).sum())

    mask = base_mask & same_ref
    n_weighted = int(mask.sum())
    if n_weighted > 0:
        log_w = exposure.log_weights_for_positions(
            unit_min_ref[mask], em_data.genomic_midpoint[mask]
        )
        em_data.gdna_log_liks[mask] += log_w.astype(
            em_data.gdna_log_liks.dtype, copy=False
        )

    return RegionalWeightApplicationStats(
        n_units_seen=n_units,
        n_units_weighted=n_weighted,
        n_units_no_gdna=n_no_gdna,
        n_units_missing_midpoint=n_missing_mid,
        n_units_cross_ref_skipped=n_cross_ref,
    )


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
    calibration_payload: "object" = None  # CalibrationScanPayload | None


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


def _strand_summary_identifiable(
    strand_summary,
    *,
    confidence: float = 0.99,
) -> bool:
    from .calibration.density_global import STRAND_CONTRAST_NUMERICAL_FLOOR

    effective_min = max(
        STRAND_CONTRAST_NUMERICAL_FLOOR,
        strand_summary.signed_strand_contrast_margin(confidence=confidence),
    )
    return abs(strand_summary.signed_strand_contrast) >= effective_min


def _calibration_strand_summary(strand_models: StrandModels):
    """Choose the strand summary used by calibration density correction."""
    from .calibration._orient import StrandSummary

    return StrandSummary.from_model(strand_models.exonic_spliced)


def _warn_if_calibration_strand_unidentifiable(strand_models: StrandModels) -> None:
    """Warn when calibration cannot identify strand from spliced RNA evidence."""
    from .calibration._orient import StrandSummary

    primary = StrandSummary.from_model(strand_models.exonic_spliced)
    if _strand_summary_identifiable(primary):
        return

    logger.warning(
        "[CAL] Spliced strand model is not identifiable at 99%% confidence "
        "(n_spliced_obs=%d, p_r1_sense=%.6f, signed_contrast=%.6f, margin_99=%.6f). "
        "Calibration will use the unstranded count/exposure estimator for channels "
        "without identifiable strand correction. If this was expected to be a stranded "
        "RNA-seq library, inspect splice evidence and contamination/nascent-RNA levels.",
        primary.n_observations,
        primary.p_r1_sense,
        primary.signed_strand_contrast,
        primary.signed_strand_contrast_margin(confidence=0.99),
    )

    diagnostic = StrandSummary.from_model(strand_models.exonic)
    if _strand_summary_identifiable(diagnostic):
        logger.warning(
            "[CAL] Exonic diagnostic strand signal is identifiable "
            "(n_exonic_obs=%d, p_r1_sense=%.6f), but it is not used for calibration. "
            "Unspliced exonic fragments can include gDNA or nascent RNA and are not an "
            "RNA-pure strand-training source.",
            diagnostic.n_observations,
            diagnostic.p_r1_sense,
        )


def scan_and_buffer(
    bam_path: str,
    index: TranscriptIndex,
    scan: "BamScanConfig",
) -> tuple[
    PipelineStats,
    StrandModels,
    FragmentLengthModels,
    FragmentBuffer,
    "object",  # CalibrationScanPayload | None
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
    tuple
        ``(stats, strand_models, frag_length_models, buffer,
        calibration_payload)`` — the last is ``None`` if the index has
        no ``regions.feather`` (legacy indexes pre-v3).
    """
    stats = PipelineStats()
    strand_models = StrandModels()
    frag_length_models = FragmentLengthModels(max_size=scan.max_frag_length)
    buffer = FragmentBuffer(
        t_strand_arr=index.t_to_strand_arr,
        chunk_size=scan.fragments_per_chunk,
        max_memory_bytes=scan.buffer_size_bytes,
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

    # Splicing-anchor tolerance K (bp): one-sided slack used by both the
    # SPLICED_IMPLICIT per-intron whole-containment discriminant and the
    # boundary-flux calibration accumulator. Set on the resolver before
    # the scanner so all _resolve_core calls see the configured K.
    resolve_ctx.set_splicing_anchor_tolerance(int(scan.splicing_anchor_tolerance))

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

    # Calibration region wiring (M3): install the per-genome region
    # partition into the scanner so per-fragment 8-state observations
    # are collected during the scan.  Skipped when the index has no
    # regions.feather (legacy index, INDEX_FORMAT_VERSION < 3).
    region_df = getattr(index, "region_df", None)
    if region_df is not None and len(region_df) > 0:
        _wire_calibration_regions(
            scanner,
            index,
            region_df,
            splicing_anchor_tolerance=int(scan.splicing_anchor_tolerance),
        )

    # Streaming chunk callback — receives zero-copy dict from C++
    def _on_chunk(raw: dict) -> None:
        chunk = _FinalizedChunk.from_raw(raw)
        buffer.inject_chunk(chunk)

    # Execute the full BAM scan in C++ with streaming chunk output
    n_scan, n_bgzf = scan.resolved_scan_threads()
    if n_bgzf != scan.bgzf_threads:
        logger.info(
            "[scan] Capped BGZF decompression threads from %d to %d to fit "
            "within total thread budget %d.",
            scan.bgzf_threads,
            n_bgzf,
            scan.resolved_total_threads(),
        )
    result = scanner.scan(
        bam_path,
        chunk_callback=_on_chunk,
        t_strand_arr=index.t_to_strand_arr.tolist(),
        chunk_size=scan.fragments_per_chunk,
        n_workers=n_scan,
        n_decomp_threads=n_bgzf,
        qname_batch_size=scan.read_name_batch_size,
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

    # Calibration payload (optional — present iff regions were wired)
    cal_dict = result.get("calibration") if isinstance(result, dict) else None
    if cal_dict is None:
        calibration_payload = None
    else:
        from .calibration.scan_payload import CalibrationScanPayload

        # Balance basis: n_read_names is the count of non-empty qname
        # groups, which is exactly the number of times the calibration
        # dispatch site fires (multimappers are noted before the early
        # return; non-multimappers reach the end-of-function policy).
        calibration_payload = CalibrationScanPayload.from_scan_dict(
            cal_dict,
            n_total=stats.n_read_names,
        )

    return stats, strand_models, frag_length_models, buffer, calibration_payload


def _wire_calibration_regions(
    scanner,
    index: TranscriptIndex,
    region_df,
    *,
    splicing_anchor_tolerance: int = 0,
) -> None:
    """Install the index's region partition into a native BamScanner.

    Translates ``ref_name`` to the resolver's internal ref_id (the same
    numbering the BAM scanner uses for ``ExonBlock.ref_id``).  Regions
    on references unknown to the resolver are dropped — the resolver
    only registers refs that have transcripts, so dropped regions are
    on intergenic-only chromosomes whose BAM fragments map with
    ``ref_id == -1`` and would never overlap any region anyway.
    """
    resolver_ref_to_id = index.resolver.get_ref_to_id()
    region_ref_names = region_df["ref_name"].astype(str).to_numpy()
    ref_ids = np.fromiter(
        (resolver_ref_to_id.get(name, -1) for name in region_ref_names),
        dtype=np.int32,
        count=len(region_ref_names),
    )
    keep = ref_ids >= 0
    if not keep.all():
        ref_ids = ref_ids[keep]
        region_df = region_df.loc[keep]

    starts = np.ascontiguousarray(region_df["start"].to_numpy(np.int64))
    ends = np.ascontiguousarray(region_df["end"].to_numpy(np.int64))
    types = region_df["type"].to_numpy(np.uint8)
    strands = np.ascontiguousarray(region_df["strand"].to_numpy(np.uint8))
    # Bit layout: bit 0 = EXON (type 2), bit 1 = INTRON (type 1),
    # bit 2 = INTERGENIC (type 0).
    type_masks = (np.uint8(1) << (np.uint8(2) - types)).astype(np.uint8)

    # set_regions requires sorted (ref_id, start) order.  region_df is
    # already in genomic order per ref but ref iteration order may not
    # match resolver IDs — re-sort to be safe.
    order = np.lexsort((starts, ref_ids))
    if not np.all(order == np.arange(len(order))):
        ref_ids = ref_ids[order]
        starts = starts[order]
        ends = ends[order]
        type_masks = type_masks[order]
        strands = strands[order]

    n_refs = max(int(resolver_ref_to_id[name]) for name in resolver_ref_to_id) + 1
    scanner.set_regions(
        np.ascontiguousarray(ref_ids),
        np.ascontiguousarray(starts),
        np.ascontiguousarray(ends),
        np.ascontiguousarray(type_masks),
        np.ascontiguousarray(strands),
        n_refs,
        int(splicing_anchor_tolerance),
    )


# ---------------------------------------------------------------------------
# Quantify from buffer (locus-level EM, no global EM)
# ---------------------------------------------------------------------------


def _setup_geometry_and_estimator(
    index: TranscriptIndex,
    rna_fl,
    em_config: EMConfig,
) -> tuple["TranscriptGeometry", AbundanceEstimator]:
    """Compute transcript geometry and create the AbundanceEstimator."""
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    if rna_fl.n_observations > 0:
        effective_lengths = rna_fl.compute_all_transcript_eff_lens(
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
    rna_fl,
    gdna_fl,
    stats: PipelineStats,
    estimator: AbundanceEstimator,
    scoring: "FragmentScoringConfig",
    log_every: int,
    annotations,
) -> "ScoredFragments":
    """Build FragmentScorer, scan buffer, and return ScoredFragments."""
    ctx = FragmentScorer.from_models(
        strand_models,
        rna_fl,
        gdna_fl,
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


def _assign_locus_ids(estimator: AbundanceEstimator, multi_loci: list) -> None:
    """Stamp ``multi_locus_id`` onto every transcript on the estimator.

    Required by the nRNA-fraction prior cascade in the C++ EM.
    """
    for locus in multi_loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.multi_locus_id


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
    multi_loci: list,
    index: TranscriptIndex,
    gdna_prior_count: np.ndarray,
    gdna_eff_len: np.ndarray,
    em_config: EMConfig,
    *,
    enable_gdna: np.ndarray | None = None,
    emit_locus_stats: bool = False,
    annotations: "AnnotationTable | None" = None,
    gdna_eff_len_unweighted: np.ndarray | None = None,
    gdna_prior_count_regional: np.ndarray | None = None,
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

    def _build_locus_meta(
        locus,
        *,
        rna_total,
        gdna,
        gdna_prior,
        gdna_leff,
        gdna_leff_unweighted=None,
        gdna_prior_regional=None,
    ):
        gene_set = {
            int(t_to_g[int(t_idx)])
            for t_idx in locus.transcript_indices
            if not is_synthetic_g[int(t_to_g[int(t_idx)])]
        }
        gdna_leff_f = float(gdna_leff)
        if gdna_leff_unweighted is None:
            gdna_leff_unw_f = gdna_leff_f
        else:
            gdna_leff_unw_f = float(gdna_leff_unweighted)
        weight_ratio = gdna_leff_f / gdna_leff_unw_f if gdna_leff_unw_f > 0.0 else 1.0
        gdna_prior_reg_f = (
            float(gdna_prior) if gdna_prior_regional is None else float(gdna_prior_regional)
        )
        return {
            "locus_id": locus.multi_locus_id,
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
            "gdna_prior_count": float(gdna_prior),
            "gdna_eff_len": gdna_leff_f,
            "gdna_eff_len_per_bp": gdna_leff_f / max(float(locus.gdna_span), 1.0),
            "gdna_eff_len_unweighted": gdna_leff_unw_f,
            "gdna_eff_len_weight_ratio": float(weight_ratio),
            "gdna_prior_count_regional": gdna_prior_reg_f,
        }

    def _call_batch_em(
        parts,
        batch_loci,
        batch_gdna_prior_count,
        batch_gdna_eff_len,
        batch_enable_gdna=None,
    ):
        """Pack tuples, call C++, record results."""
        partition_tuples = [
            (
                p.offsets,
                p.t_indices,
                p.log_liks,
                p.coverage_weights,
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
            batch_gdna_prior_count,
            index,
            gdna_eff_len=batch_gdna_eff_len,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
            emit_locus_stats=emit_locus_stats,
            emit_assignments=emit_assignments,
            enable_gdna=batch_enable_gdna,
        )

    # Classify mega vs normal
    locus_work = {
        loc.multi_locus_id: len(loc.transcript_indices) * partitions[loc.multi_locus_id].n_units
        for loc in multi_loci
    }
    total_work = sum(locus_work.values())
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [loc for loc in multi_loci if locus_work[loc.multi_locus_id] >= fair_share],
        key=lambda loc: locus_work[loc.multi_locus_id],
        reverse=True,
    )
    mega_ids = {loc.multi_locus_id for loc in mega_loci}

    total_gdna_em = 0.0

    # Phase A: Mega-loci (one at a time, free after each)
    for locus in mega_loci:
        part = partitions.pop(locus.multi_locus_id)
        lid = locus.multi_locus_id
        em_result = _call_batch_em(
            [part],
            [locus],
            np.array([gdna_prior_count[lid]], dtype=np.float64),
            np.array([gdna_eff_len[lid]], dtype=np.float64),
            batch_enable_gdna=(
                np.array([enable_gdna[lid]], dtype=np.uint8) if enable_gdna is not None else None
            ),
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
                gdna_prior=gdna_prior_count[lid],
                gdna_leff=gdna_eff_len[lid],
                gdna_leff_unweighted=(
                    gdna_eff_len_unweighted[lid]
                    if gdna_eff_len_unweighted is not None
                    else None
                ),
                gdna_prior_regional=(
                    gdna_prior_count_regional[lid]
                    if gdna_prior_count_regional is not None
                    else None
                ),
            )
        )
        del part
        logger.debug(
            f"[MEGA] Locus {locus.multi_locus_id}: "
            f"{len(locus.transcript_indices)} transcripts, "
            f"{len(locus.unit_indices)} units"
        )

    # Phase B: Normal loci (one batched call)
    normal_loci = [loc for loc in multi_loci if loc.multi_locus_id not in mega_ids]
    if normal_loci:
        # Pop partitions out of the dict so the only references during
        # the batched C++ call live in ``normal_parts``.  After the call
        # completes and annotations are written we drop ``normal_parts``
        # to release per-locus arrays before EM downstream phases run.
        normal_parts = [partitions.pop(loc.multi_locus_id) for loc in normal_loci]
        normal_gp = np.array(
            [gdna_prior_count[loc.multi_locus_id] for loc in normal_loci],
            dtype=np.float64,
        )
        normal_gdna_eff_len = np.array(
            [gdna_eff_len[loc.multi_locus_id] for loc in normal_loci],
            dtype=np.float64,
        )
        em_result = _call_batch_em(
            normal_parts,
            normal_loci,
            normal_gp,
            normal_gdna_eff_len,
            batch_enable_gdna=(
                np.array(
                    [enable_gdna[loc.multi_locus_id] for loc in normal_loci],
                    dtype=np.uint8,
                )
                if enable_gdna is not None
                else None
            ),
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
            lid = locus.multi_locus_id
            estimator.locus_results.append(
                _build_locus_meta(
                    locus,
                    rna_total=rna_arr[i],
                    gdna=gdna_arr[i],
                    gdna_prior=normal_gp[i],
                    gdna_leff=normal_gdna_eff_len[i],
                    gdna_leff_unweighted=(
                        gdna_eff_len_unweighted[lid]
                        if gdna_eff_len_unweighted is not None
                        else None
                    ),
                    gdna_prior_regional=(
                        gdna_prior_count_regional[lid]
                        if gdna_prior_count_regional is not None
                        else None
                    ),
                )
            )
        # Release per-locus partition arrays before downstream phases.
        del normal_parts, em_result

    del partitions

    estimator._gdna_em_total = total_gdna_em

    n_total_units = sum(len(loc.unit_indices) for loc in multi_loci)
    logger.info(
        f"[DONE] Per-locus EM (partitioned): {len(multi_loci)} loci "
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
    calibration_payload: "CalibrationScanPayload",
    *,
    em_config: EMConfig | None = None,
    scoring: FragmentScoringConfig | None = None,
    log_every: int = 1_000_000,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> tuple[AbundanceEstimator, "CalibrationResult"]:
    """Quantify transcripts from buffered fragments via locus-level EM.

    Consumes the v6 :class:`CalibrationResult` produced by
    :func:`rigel.calibration.calibrate`.  Per-``MultiLocus`` priors are
    assembled here (after ``build_multi_loci``) via
    :func:`assemble_priors` and back-filled into the calibration result
    with :meth:`CalibrationResult.with_priors`; the populated result is
    returned alongside the estimator so callers (CLI, tests) can read
    ``gdna_prior_count`` and diagnostic dataframes.

    Parameters
    ----------
    buffer : FragmentBuffer
    index : TranscriptIndex
    strand_models : StrandModels
    frag_length_models : FragmentLengthModels
        Scanner-trained accumulator.  Carries raw histograms used for
        index-side geometry; the FL distributions consumed by scoring
        come from ``calibration.fl_models``, NOT this object.
    stats : PipelineStats
    calibration : CalibrationResult
        v6 calibration result.  Must carry a populated
        ``fl_models`` field; per-locus priors may be empty (they are
        filled in here).
    calibration_payload : CalibrationScanPayload
        The C++ scanner's payload.  Required by
        :func:`assemble_priors` for per-region counts and boundary-flux
        counters.
    em_config, scoring, log_every, annotations, emit_locus_stats
        Standard pipeline knobs.

    Returns
    -------
    (AbundanceEstimator, CalibrationResult)
        The fitted estimator and the calibration result with priors
        backfilled (zero-locus priors when no MultiLoci were built).
    """
    if calibration is None:
        raise ValueError(
            "quant_from_buffer() requires a v6 CalibrationResult "
            "(got None).  Run rigel.calibration.calibrate(...) before "
            "locus-level quantification."
        )
    if em_config is None:
        em_config = EMConfig()
    if scoring is None:
        scoring = FragmentScoringConfig()

    # The scorer reads finalized FL distributions directly from
    # CalibrationResult.fl_models.  No backflow mutation of
    # frag_length_models.
    fl_models = calibration.fl_models

    # Phase 1: Geometry + estimator
    geometry, estimator = _setup_geometry_and_estimator(
        index,
        fl_models.rna,
        em_config,
    )

    logger.info(
        f"[START] Quantifying {buffer.total_fragments:,} buffered fragments "
        f"(locus-level EM: mRNA/nRNA/gDNA)"
    )

    logger.info(
        f"[CAL-FL] gDNA FL: mean={fl_models.gdna.mean:.1f}, "
        f"quality={fl_models.gdna_quality} "
        f"(n_pool={fl_models.n_gdna})"
    )

    # Phase 2: Score fragments
    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        fl_models.rna,
        fl_models.gdna,
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
        multi_loci = build_multi_loci(em_data, index)

        if multi_loci:
            max_locus_t = max(len(loc.transcript_indices) for loc in multi_loci)
            max_locus_u = max(len(loc.unit_indices) for loc in multi_loci)
        else:
            max_locus_t = max_locus_u = 0

        logger.info(
            f"[LOCI] {len(multi_loci)} loci from {em_data.n_units:,} units "
            f"(largest: {max_locus_t} transcripts, {max_locus_u} units)"
        )

        _assign_locus_ids(estimator, multi_loci)

        # v6 calibration: assemble the PriorTable for this batch and
        # back-fill it into the calibration result.  Pulls global
        # density posteriors from ``calibration.global_densities`` and
        # the gDNA FL mean from ``fl_models.gdna``.
        # NOTE: must run BEFORE ``partition_and_free`` because the
        # latter nulls out ``em_data`` arrays as it scatters them.
        prior_table = assemble_priors(
            multi_loci,
            em_data,
            index,
            calibration_payload,
            calibration.global_densities,
            gdna_fl=fl_models.gdna,
            splicing_anchor_tolerance=int(
                getattr(calibration_payload, "splicing_anchor_tolerance", 0)
            ),
            regional_exposure=calibration.regional_exposure,
        )
        calibration = calibration.with_priors(prior_table)

        # Regional gDNA exposure: attenuate per-unit gDNA log-likelihoods
        # using the per-region weight ``log A_r`` at the unit midpoint.
        # Must run after ``assemble_priors`` (single owner of the global
        # ``gdna_log_liks`` array) and before ``partition_and_free``.
        if calibration.regional_exposure is not None:
            weight_stats = _apply_unit_gdna_weights(
                em_data, calibration.regional_exposure, index
            )
            calibration = calibration.with_regional_weighting_stats(weight_stats)
        # Log the prior-derived calibration summary now that the
        # PriorTable is back-filled (the earlier "[CAL] v6 quality=..."
        # line ran before assembly and could not include these fields).
        _post_summary = calibration.to_summary_dict()
        logger.info(
            "[CAL] v6 priors  mean_pi_gdna=%.4f  n_multi_loci=%d",
            float(_post_summary["mean_pi_gdna"]),
            int(_post_summary["n_multi_loci"]),
        )
        gdna_prior_count = prior_table.gdna_prior_count
        gdna_eff_len = prior_table.gdna_eff_len
        enable_gdna_arr = prior_table.enable_gdna

        # Phase 4 (NEW): Fused scatter into per-locus tuples
        # Phase 4 (NEW): Array-by-array scatter + incremental free
        partitions = partition_and_free(em_data, multi_loci)
        del em_data

        # Phase 5 (NEW): Streaming locus EM with incremental partition freeing
        _run_locus_em_partitioned(
            estimator,
            partitions,
            multi_loci,
            index,
            gdna_prior_count,
            gdna_eff_len,
            em_config,
            enable_gdna=enable_gdna_arr,
            emit_locus_stats=emit_locus_stats,
            annotations=annotations,
            gdna_eff_len_unweighted=prior_table.gdna_eff_len_unweighted,
            gdna_prior_count_regional=prior_table.gdna_prior_count_regional,
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")
        del em_data

    _gdna_em = estimator.gdna_em_count
    stats.n_gdna_em = 0 if (math.isnan(_gdna_em) or math.isinf(_gdna_em)) else int(_gdna_em)

    # Global gDNA estimate from calibration densities (projected over
    # the full genome, including exonic regions invisible to direct
    # classification).
    region_df = getattr(index, "region_df", None)
    if region_df is not None and len(region_df) > 0 and calibration is not None:
        from rigel.calibration.density_global import estimate_global_gdna_fragments

        gdna_global = estimate_global_gdna_fragments(calibration.global_densities, region_df)
    else:
        gdna_global = float(_gdna_em)
    stats.n_gdna_global = (
        int(gdna_global) if not (math.isnan(gdna_global) or math.isinf(gdna_global)) else 0
    )

    total_frags = float(estimator.unambig_counts.sum() + estimator.em_counts.sum()) + gdna_global
    gdna_rate_global = gdna_global / total_frags if total_frags > 0 else 0.0

    logger.info(
        f"[gDNA] em_locus={_gdna_em:.0f}, "
        f"global_projected={gdna_global:.0f}, "
        f"contamination_rate={gdna_rate_global:.2%}"
    )

    return estimator, calibration


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
    stats, strand_models, frag_length_models, buffer, calibration_payload = scan_and_buffer(
        bam_path, index, scan
    )

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
    # NOTE: v6 calibration (rigel.calibration.calibrate) builds its own
    # finalised FLModels (RNA + gDNA + global) inside CalibrationResult.
    # We do NOT call ``frag_length_models.build_scoring_models()`` or
    # ``.finalize(...)`` here — the scanner-trained accumulator is kept
    # raw and only consulted for index-side geometry.

    # -- v6 calibration (single FL build, global density posteriors) --
    from .calibration import calibrate

    cal_cfg = config.calibration
    strand_summary = _calibration_strand_summary(strand_models)
    _warn_if_calibration_strand_unidentifiable(strand_models)
    strand_ci_eps = strand_models.strand_specificity_ci_epsilon(confidence=0.99)
    logger.info(
        "[CAL] Strand trainer: n_spliced_obs=%d  ss_est=%.6f  ε_CI(99%%)=%.4g",
        strand_models.n_observations,
        strand_models.strand_specificity,
        strand_ci_eps,
    )

    calibration = calibrate(
        index=index,
        payload=calibration_payload,
        scan_trained=frag_length_models,
        fl_prior_ess=cal_cfg.prior_ess,
        pool_quality_good=cal_cfg.pool_quality_good,
        pool_quality_weak=cal_cfg.pool_quality_weak,
        strand_summary=strand_summary,
        regional_exposure_enabled=cal_cfg.regional_exposure_enabled,
    )
    cal_summary = calibration.to_summary_dict()
    logger.info(
        "[CAL] v6 quality=%s",
        cal_summary["fl_models"]["gdna_quality"],
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
        estimator, calibration = quant_from_buffer(
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
            calibration_payload=calibration_payload,
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
        calibration_payload=calibration_payload,
    )
