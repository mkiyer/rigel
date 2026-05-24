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

**Quantification** (``quant_from_buffer``): currently fails fast after
calibration while the Phase 4 fractional density and per-locus prior
pipeline is being developed. The scorer / locus-EM helpers below remain
as the Phase 4 re-entry points.

Scoring functions live in ``scoring.py``.  Locus construction and EM
initialization live in ``locus.py``.  The CSR builder lives in
``scan.py``.  This module is a thin orchestrator.
"""

from __future__ import annotations

import logging
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
from .native import BamScanner as _NativeBamScanner
from .native import detect_sj_strand_tag as _native_detect_sj_tag
from .scan import FragmentRouter
from .scoring import FragmentScorer
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from .annotate import AnnotationTable
    from .calibration import CalibrationResult
    from .calibration.scan_payload import CalibrationScanPayload
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
    from .calibration.strand_summary import STRAND_CONTRAST_NUMERICAL_FLOOR

    effective_min = max(
        STRAND_CONTRAST_NUMERICAL_FLOOR,
        strand_summary.signed_strand_contrast_margin(confidence=confidence),
    )
    return abs(strand_summary.signed_strand_contrast) >= effective_min


def _calibration_strand_summary(strand_models: StrandModels):
    """Choose the strand summary used by calibration density correction."""
    from .calibration.strand_summary import StrandSummary

    return StrandSummary.from_model(strand_models.exonic_spliced)


def _warn_if_calibration_strand_unidentifiable(strand_models: StrandModels) -> None:
    """Warn when calibration cannot identify strand from spliced RNA evidence."""
    from .calibration.strand_summary import StrandSummary

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

    # Splicing-anchor tolerance K (bp): resolver-only one-sided slack used
    # by the SPLICED_IMPLICIT per-intron whole-containment discriminant.
    # The fractional calibration accumulator records raw compartment mass
    # directly and does not consume this tolerance.
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

    # Calibration region wiring: install the per-genome fine-region
    # partition into the scanner so fractional per-region evidence is
    # collected during the scan.
    region_df = getattr(index, "region_df", None)
    if region_df is not None and len(region_df) > 0:
        _wire_calibration_regions(
            scanner,
            index,
            region_df,
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
) -> None:
    """Install the index's region partition into a native BamScanner.

    Translates ``ref_name`` to the index's canonical reference ids. This
    mapping matches ``index.ref_name_to_id`` / BAM tid space and includes
    references that have no transcripts, which the resolver-local map may
    omit.
    """
    ref_name_to_id = index.ref_name_to_id
    region_ref_names = region_df["ref_name"].astype(str).to_numpy()
    ref_ids = np.fromiter(
        (ref_name_to_id.get(name, -1) for name in region_ref_names),
        dtype=np.int32,
        count=len(region_ref_names),
    )
    keep = ref_ids >= 0
    if not keep.all():
        dropped_refs = sorted(set(region_ref_names[~keep].tolist()))
        logger.warning(
            "[scan] Dropping %d calibration regions on %d references absent "
            "from index.ref_name_to_id: %s",
            int((~keep).sum()),
            len(dropped_refs),
            ", ".join(dropped_refs[:10]) + ("..." if len(dropped_refs) > 10 else ""),
        )
        ref_ids = ref_ids[keep]
        region_df = region_df.loc[keep]
    if len(region_df) == 0:
        logger.warning("[scan] No calibration regions remain after reference-id filtering.")
        return

    starts = np.ascontiguousarray(region_df["start"].to_numpy(np.int64))
    ends = np.ascontiguousarray(region_df["end"].to_numpy(np.int64))
    signatures = np.ascontiguousarray(region_df["signature"].to_numpy(np.uint8))
    left_signatures = np.ascontiguousarray(region_df["left_signature"].to_numpy(np.uint8))
    right_signatures = np.ascontiguousarray(region_df["right_signature"].to_numpy(np.uint8))
    boundary_kind_left = np.ascontiguousarray(region_df["boundary_kind_left"].to_numpy(np.uint8))
    boundary_kind_right = np.ascontiguousarray(region_df["boundary_kind_right"].to_numpy(np.uint8))
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
        signatures = signatures[order]
        left_signatures = left_signatures[order]
        right_signatures = right_signatures[order]
        boundary_kind_left = boundary_kind_left[order]
        boundary_kind_right = boundary_kind_right[order]
        type_masks = type_masks[order]
        strands = strands[order]

    n_refs = len(ref_name_to_id)
    scanner.set_regions(
        np.ascontiguousarray(ref_ids),
        np.ascontiguousarray(starts),
        np.ascontiguousarray(ends),
        np.ascontiguousarray(signatures),
        np.ascontiguousarray(left_signatures),
        np.ascontiguousarray(right_signatures),
        np.ascontiguousarray(boundary_kind_left),
        np.ascontiguousarray(boundary_kind_right),
        np.ascontiguousarray(type_masks),
        np.ascontiguousarray(strands),
        n_refs,
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
        effective_lengths_em=None,
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
    gdna_prior_count_em: np.ndarray,
    gdna_eff_len: np.ndarray,
    em_config: EMConfig,
    *,
    enable_gdna: np.ndarray | None = None,
    emit_locus_stats: bool = False,
    annotations: "AnnotationTable | None" = None,
    gdna_eff_len_unweighted: np.ndarray | None = None,
    gdna_prior_count_raw: np.ndarray | None = None,
    gdna_em_exposure_weight: np.ndarray | None = None,
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
        gdna_prior_em,
        gdna_leff,
        gdna_leff_unweighted=None,
        gdna_prior_raw=None,
        gdna_weight=None,
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
        gdna_prior_em_f = float(gdna_prior_em)
        gdna_prior_raw_f = gdna_prior_em_f if gdna_prior_raw is None else float(gdna_prior_raw)
        gdna_weight_f = weight_ratio if gdna_weight is None else float(gdna_weight)
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
            "gdna_prior_count": gdna_prior_raw_f,
            "gdna_prior_count_em": gdna_prior_em_f,
            "gdna_eff_len": gdna_leff_f,
            "gdna_eff_len_per_bp": gdna_leff_f / max(float(locus.gdna_span), 1.0),
            "gdna_eff_len_unweighted": gdna_leff_unw_f,
            "gdna_eff_len_weight_ratio": float(weight_ratio),
            "gdna_em_exposure_weight": gdna_weight_f,
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
            np.array([gdna_prior_count_em[lid]], dtype=np.float64),
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
                gdna_prior_em=gdna_prior_count_em[lid],
                gdna_leff=gdna_eff_len[lid],
                gdna_leff_unweighted=(
                    gdna_eff_len_unweighted[lid] if gdna_eff_len_unweighted is not None else None
                ),
                gdna_prior_raw=(
                    gdna_prior_count_raw[lid] if gdna_prior_count_raw is not None else None
                ),
                gdna_weight=(
                    gdna_em_exposure_weight[lid] if gdna_em_exposure_weight is not None else None
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
            [gdna_prior_count_em[loc.multi_locus_id] for loc in normal_loci],
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
                    gdna_prior_em=normal_gp[i],
                    gdna_leff=normal_gdna_eff_len[i],
                    gdna_leff_unweighted=(
                        gdna_eff_len_unweighted[lid]
                        if gdna_eff_len_unweighted is not None
                        else None
                    ),
                    gdna_prior_raw=(
                        gdna_prior_count_raw[lid] if gdna_prior_count_raw is not None else None
                    ),
                    gdna_weight=(
                        gdna_em_exposure_weight[lid]
                        if gdna_em_exposure_weight is not None
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
    """Fail fast at the Phase 6 locus-EM boundary.

    The scan and calibration stages are complete when this function is
    called. Locus-level scoring, prior assembly, and EM are intentionally
    blocked until the Phase 6 wiring lands.

    Parameters
    ----------
    calibration : CalibrationResult
        Fractional calibration result. Must carry a populated ``fl_models``
        and diagnostics block.
    Other parameters
        Retained for API stability and Phase 4 re-entry.

    Returns
    -------
    Never returns until Phase 6 replaces the fail-fast boundary.
    """
    if calibration is None:
        raise ValueError(
            "quant_from_buffer() requires a v6 CalibrationResult "
            "(got None).  Run rigel.calibration.calibrate(...) before "
            "locus-level quantification."
        )

    raise NotImplementedError("rigel quant: locus EM lands in Phase 6")


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
        rna_lower_confidence=cal_cfg.rna_lower_confidence,
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
