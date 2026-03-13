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

import gc
import logging
import math
import os
from dataclasses import dataclass, replace as _replace
from pathlib import Path

import numpy as np
import pandas as pd

# Padding and minimum capacity for the annotation table.
_ANNOTATION_TABLE_PADDING = 1024
_ANNOTATION_TABLE_MIN_CAPACITY = 4096

from .buffer import FragmentBuffer, _FinalizedChunk
from .splice import SpliceType
from .estimator import (
    AbundanceEstimator,
    compute_global_gdna_density,
)
from .types import ChimeraType, Strand
from .index import TranscriptIndex
from .frag_length_model import FragmentLengthModels
from .stats import PipelineStats
from .strand_model import StrandModels
from .native import BamScanner as _NativeBamScanner
from .native import detect_sj_strand_tag as _native_detect_sj_tag

# --- New modular imports ---
from .scoring import FragmentScorer
from .locus import (
    build_loci,
    compute_nrna_init,
    compute_gdna_density_from_strand,
    compute_gdna_density_hybrid,
    compute_eb_gdna_priors,
)
from .scan import FragmentRouter
from .config import (
    EMConfig,
    PipelineConfig,
    BamScanConfig,
    FragmentScoringConfig,
    TranscriptGeometry,
)

logger = logging.getLogger(__name__)

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


def _sj_tag_to_spec(sj_strand_tag) -> str:
    """Convert BamScanConfig.sj_strand_tag to the string spec for BamScanner."""
    if isinstance(sj_strand_tag, str):
        return sj_strand_tag if sj_strand_tag else "none"
    if isinstance(sj_strand_tag, (list, tuple)):
        return ",".join(sj_strand_tag) if sj_strand_tag else "none"
    return "none"


def _compute_intergenic_density(
    stats: PipelineStats,
    index: TranscriptIndex,
) -> float:
    """Compute background intergenic density (frags / bp).

    Uses reference lengths from the index and the merged
    genic union span to derive intergenic territory, then normalises
    the intergenic fragment count to obtain per-bp density.

    Returns 0.0 if reference lengths are unavailable (e.g. in unit
    tests) or no intergenic fragments exist.
    """
    n_intergenic = stats.n_intergenic
    if n_intergenic == 0:
        return 0.0

    # Load reference lengths from the index directory
    index_dir = getattr(index, "index_dir", None)
    if index_dir is None:
        return 0.0
    ref_path = os.path.join(index_dir, "ref_lengths.feather")
    if not os.path.exists(ref_path):
        return 0.0

    ref_df = pd.read_feather(ref_path)
    total_genome_bp = float(ref_df["length"].sum())
    if total_genome_bp <= 0:
        return 0.0

    # Compute genic union span (merge overlapping transcript intervals)
    genic_bp = 0
    for _, grp in index.t_df.groupby("ref"):
        starts = grp["start"].values
        ends = grp["end"].values
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        # Merge overlapping intervals
        ms = int(starts[0])
        me = int(ends[0])
        for j in range(1, len(starts)):
            if int(starts[j]) <= me:
                me = max(me, int(ends[j]))
            else:
                genic_bp += me - ms
                ms = int(starts[j])
                me = int(ends[j])
        genic_bp += me - ms

    intergenic_bp = max(total_genome_bp - genic_bp, 1.0)
    return n_intergenic / intergenic_bp


def _replay_strand_observations(
    strand_dict: dict,
    strand_models: StrandModels,
) -> None:
    """Replay C++ strand observation arrays into Python StrandModels."""
    # Exonic-spliced model
    obs = strand_dict.get("exonic_spliced_obs", [])
    truth = strand_dict.get("exonic_spliced_truth", [])
    if len(obs) > 0:
        strand_models.exonic_spliced.observe_batch(
            np.asarray(obs, dtype=np.int8),
            np.asarray(truth, dtype=np.int8),
        )

    # Exonic fallback model
    obs = strand_dict.get("exonic_obs", [])
    truth = strand_dict.get("exonic_truth", [])
    if len(obs) > 0:
        strand_models.exonic.observe_batch(
            np.asarray(obs, dtype=np.int8),
            np.asarray(truth, dtype=np.int8),
        )

    # Intergenic model
    obs = strand_dict.get("intergenic_obs", [])
    truth = strand_dict.get("intergenic_truth", [])
    if len(obs) > 0:
        strand_models.intergenic.observe_batch(
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

    unspliced_same = fraglen_dict.get("unspliced_same_strand_lengths", [])
    if len(unspliced_same) > 0:
        frag_length_models.unspliced_same_strand.observe_batch(
            np.asarray(unspliced_same, dtype=np.intp),
        )
    unspliced_opp = fraglen_dict.get("unspliced_opp_strand_lengths", [])
    if len(unspliced_opp) > 0:
        frag_length_models.unspliced_opp_strand.observe_batch(
            np.asarray(unspliced_opp, dtype=np.intp),
        )


def _apply_scan_stats(stats: PipelineStats, stats_dict: dict) -> None:
    """Apply C++ scan statistics to PipelineStats."""
    # BAM-level stats (stored as bam_stats sub-dict fields)
    stats.total = stats_dict.get("total", 0)
    stats.qc_fail = stats_dict.get("qc_fail", 0)
    stats.unmapped = stats_dict.get("unmapped", 0)
    stats.secondary = stats_dict.get("secondary", 0)
    stats.supplementary = stats_dict.get("supplementary", 0)
    stats.duplicate = stats_dict.get("duplicate", 0)
    stats.n_read_names = stats_dict.get("n_read_names", 0)
    stats.unique = stats_dict.get("unique", 0)
    stats.multimapping = stats_dict.get("multimapping", 0)
    stats.proper_pair = stats_dict.get("proper_pair", 0)
    stats.improper_pair = stats_dict.get("improper_pair", 0)
    stats.mate_unmapped = stats_dict.get("mate_unmapped", 0)

    # Fragment-level stats
    stats.n_fragments = stats_dict.get("n_fragments", 0)
    stats.n_chimeric = stats_dict.get("n_chimeric", 0)
    stats.n_chimeric_trans = stats_dict.get("n_chimeric_trans", 0)
    stats.n_chimeric_cis_strand_same = stats_dict.get(
        "n_chimeric_cis_strand_same", 0)
    stats.n_chimeric_cis_strand_diff = stats_dict.get(
        "n_chimeric_cis_strand_diff", 0)
    stats.n_intergenic_unspliced = stats_dict.get(
        "n_intergenic_unspliced", 0)
    stats.n_intergenic_spliced = stats_dict.get("n_intergenic_spliced", 0)
    stats.n_with_exon = stats_dict.get("n_with_exon", 0)
    stats.n_with_annotated_sj = stats_dict.get("n_with_annotated_sj", 0)
    stats.n_with_unannotated_sj = stats_dict.get(
        "n_with_unannotated_sj", 0)
    stats.n_same_strand = stats_dict.get("n_same_strand", 0)
    stats.n_ambig_strand = stats_dict.get("n_ambig_strand", 0)
    stats.n_strand_trained = stats_dict.get("n_strand_trained", 0)
    stats.n_strand_skipped_no_sj = stats_dict.get(
        "n_strand_skipped_no_sj", 0)
    stats.n_strand_skipped_ambig_strand = stats_dict.get(
        "n_strand_skipped_ambig_strand", 0)
    stats.n_strand_skipped_ambiguous = stats_dict.get(
        "n_strand_skipped_ambiguous", 0)
    stats.n_frag_length_unambiguous = stats_dict.get(
        "n_frag_length_unambiguous", 0)
    stats.n_frag_length_ambiguous = stats_dict.get(
        "n_frag_length_ambiguous", 0)
    stats.n_frag_length_intergenic = stats_dict.get(
        "n_frag_length_intergenic", 0)
    stats.n_multimapper_groups = stats_dict.get("n_multimapper_groups", 0)
    stats.n_multimapper_alignments = stats_dict.get(
        "n_multimapper_alignments", 0)


def scan_and_buffer(
    bam_path: str,
    index: TranscriptIndex,
    scan: "BamScanConfig",
) -> tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]:
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
    strand_models.strand_prior_kappa = scan.strand_prior_kappa
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

    # Create the native scanner
    sj_spec = _sj_tag_to_spec(scan.sj_strand_tag)
    scanner = _NativeBamScanner(
        resolve_ctx,
        sj_spec,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
    )

    # Execute the full BAM scan in C++
    n_scan = scan.n_scan_threads
    if n_scan == 0:
        import os
        n_scan = os.cpu_count() or 1
    result = scanner.scan(bam_path, n_workers=n_scan)

    # Replay stats
    _apply_scan_stats(stats, result["stats"])

    # Replay strand observations into Python models
    _replay_strand_observations(result["strand_observations"], strand_models)

    # Replay fragment-length observations into Python models
    _replay_fraglen_observations(
        result["frag_length_observations"], frag_length_models,
    )

    # Finalize the C++ accumulator to raw bytes, then build buffer chunk
    acc_size = result["accumulator_size"]
    if acc_size > 0:
        raw = scanner.finalize_accumulator(index.t_to_strand_arr.tolist())

        chunk = _FinalizedChunk.from_raw(raw)

        # # Sort by frag_id to restore deterministic order after
        # # multi-threaded scanning (frag_id is assigned sequentially
        # # by the single-threaded BAM reader).
        # chunk = chunk.sort_by_frag_id()

        buffer.inject_chunk(chunk)

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
    if frag_length_models.global_model.n_observations > 0:
        mean_frag = frag_length_models.global_model.mean
    else:
        from .estimator import _DEFAULT_MEAN_FRAG
        mean_frag = _DEFAULT_MEAN_FRAG

    if frag_length_models.global_model.n_observations > 0:
        effective_lengths = (
            frag_length_models.global_model.compute_all_transcript_eff_lens(
                exonic_lengths.astype(np.int64),
            )
        )
    else:
        effective_lengths = np.maximum(
            exonic_lengths - mean_frag + 1.0, 1.0,
        )

    transcript_spans = (
        index.t_df["end"].values - index.t_df["start"].values
    ).astype(np.float64)

    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        mean_frag=mean_frag,
        transcript_spans=transcript_spans,
    )

    estimator = AbundanceEstimator(
        index.num_transcripts,
        num_nrna=index.num_nrna,
        t_to_nrna=index.t_to_nrna_arr,
        em_config=em_config,
        geometry=geometry,
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
        strand_models, frag_length_models, index, estimator,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
    )
    builder = FragmentRouter(
        ctx, estimator, stats, index, strand_models,
        annotations=annotations,
    )
    em_data = builder.scan(buffer, log_every)

    # Free scanner accumulators + buffer: scan is done
    del builder, ctx
    buffer.release()
    gc.collect()
    return em_data


def _compute_priors(
    estimator: AbundanceEstimator,
    em_data: "ScoredFragments",
    loci: list,
    index: TranscriptIndex,
    strand_models: StrandModels,
    stats: PipelineStats,
    geometry: "TranscriptGeometry",
    em_config: EMConfig,
) -> list[float]:
    """Compute nRNA init, EB gDNA priors, and nRNA-frac Beta priors.

    Returns the per-locus gdna_inits array.
    """
    nrna_spans = (
        index.nrna_df["end"].values - index.nrna_df["start"].values
    ).astype(np.float64)
    nrna_max_exonic = np.zeros(index.num_nrna, dtype=np.float64)
    np.maximum.at(nrna_max_exonic, index.t_to_nrna_arr, geometry.exonic_lengths)
    nrna_init = compute_nrna_init(
        estimator.transcript_intronic_sense,
        estimator.transcript_intronic_antisense,
        nrna_spans,
        nrna_max_exonic,
        strand_models,
    )
    estimator.nrna_init = nrna_init

    strand_spec = strand_models.strand_specificity
    w_strand = (2.0 * strand_spec - 1.0) ** 2
    logger.info(
        f"[INFO] Strand specificity: {strand_spec:.3f} "
        f"(strand weight={w_strand:.3f}, density weight={1.0 - w_strand:.3f})"
    )

    # Assign locus_id to every transcript (needed by nrna_frac prior cascade)
    for locus in loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.locus_id

    # EB gDNA priors (must precede nrna_frac priors)
    intergenic_density = _compute_intergenic_density(stats, index)
    gdna_inits = compute_eb_gdna_priors(
        loci, em_data, estimator, index, strand_models,
        intergenic_density=intergenic_density,
        kappa_ref=em_config.gdna_kappa_ref,
        kappa_locus=em_config.gdna_kappa_locus,
        mom_min_evidence_ref=em_config.gdna_mom_min_evidence_ref,
        mom_min_evidence_locus=em_config.gdna_mom_min_evidence_locus,
        kappa_min=em_config.gdna_kappa_min,
        kappa_max=em_config.gdna_kappa_max,
        kappa_fallback=em_config.gdna_kappa_fallback,
        kappa_min_obs=em_config.gdna_kappa_min_obs,
    )

    return gdna_inits


def _run_locus_em(
    estimator: AbundanceEstimator,
    em_data: "ScoredFragments",
    loci: list,
    index: TranscriptIndex,
    gdna_inits: list[float],
    em_config: EMConfig,
) -> None:
    """Run batch locus-level EM and record per-locus results."""
    t_refs = index.t_df["ref"].values
    t_to_g = index.t_to_g_arr

    def _build_locus_meta(locus, *, mrna, nrna, gdna, gdna_init):
        """Build one locus_results record."""
        ref_counts: dict[str, int] = {}
        for t_idx in locus.transcript_indices:
            ref = str(t_refs[int(t_idx)])
            ref_counts[ref] = ref_counts.get(ref, 0) + 1
        primary_ref = (
            max(ref_counts, key=ref_counts.get) if ref_counts else ""
        )
        gene_set = {int(t_to_g[int(t_idx)]) for t_idx in locus.transcript_indices}
        return {
            "locus_id": locus.locus_id,
            "ref": primary_ref,
            "n_transcripts": len(locus.transcript_indices),
            "n_genes": len(gene_set),
            "n_em_fragments": len(locus.unit_indices),
            "mrna": float(mrna),
            "nrna": float(nrna),
            "gdna": float(gdna),
            "gdna_init": float(gdna_init),
        }

    (
        total_gdna_em,
        locus_mrna_arr,
        locus_nrna_arr,
        locus_gdna_arr,
    ) = estimator.run_batch_locus_em(
        loci,
        em_data,
        index,
        np.asarray(gdna_inits, dtype=np.float64),
        em_iterations=em_config.iterations,
        em_convergence_delta=em_config.convergence_delta,
        confidence_threshold=em_config.confidence_threshold,
    )

    for i, locus in enumerate(loci):
        estimator.locus_results.append(_build_locus_meta(
            locus,
            mrna=locus_mrna_arr[i],
            nrna=locus_nrna_arr[i],
            gdna=locus_gdna_arr[i],
            gdna_init=gdna_inits[i],
        ))

    estimator._gdna_em_total = total_gdna_em

    logger.info(
        f"[DONE] Per-locus EM: {len(loci)} loci, "
        f"{em_data.n_units:,} ambiguous fragments, "
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
        index, frag_length_models, em_config,
    )

    logger.info(
        f"[START] Quantifying {buffer.total_fragments:,} buffered fragments "
        f"(locus-level EM: mRNA/nRNA/gDNA)"
    )

    # Phase 2: Score fragments
    em_data = _score_fragments(
        buffer, index, strand_models, frag_length_models,
        stats, estimator, scoring, log_every, annotations,
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
            max_locus_t = max(len(l.transcript_indices) for l in loci)
            max_locus_u = max(len(l.unit_indices) for l in loci)
        else:
            max_locus_t = max_locus_u = 0

        logger.info(
            f"[LOCI] {len(loci)} loci from {em_data.n_units:,} units "
            f"(largest: {max_locus_t} transcripts, {max_locus_u} units)"
        )

        gdna_inits = _compute_priors(
            estimator, em_data, loci, index, strand_models,
            stats, geometry, em_config,
        )

        _run_locus_em(
            estimator, em_data, loci, index, gdna_inits, em_config,
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")

    # Phase 4: Cleanup
    del em_data
    gc.collect()

    _gdna_em = estimator.gdna_em_count
    stats.n_gdna_em = 0 if (math.isnan(_gdna_em) or math.isinf(_gdna_em)) else int(_gdna_em)

    logger.info(
        f"[gDNA] total={estimator.gdna_total:.0f} "
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
    stats, strand_models, frag_length_models, buffer = (
        scan_and_buffer(bam_path, index, scan)
    )

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
    frag_length_models.mix_models(
        s_rna=strand_models.strand_specificity,
        p_r1_sense=strand_models.exonic_spliced.p_r1_sense,
    )
    frag_length_models.finalize()

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
        )

    return PipelineResult(
        stats=stats,
        strand_models=strand_models,
        frag_length_models=frag_length_models,
        estimator=estimator,
    )
