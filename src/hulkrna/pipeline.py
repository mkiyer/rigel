"""
hulkrna.pipeline — Single-pass BAM counting pipeline with locus-level EM.

Architecture
------------
The pipeline reads the BAM file **once** and processes counts in
two logical stages:

**BAM Scan** (``scan_and_buffer``): Parse all fragments, resolve each
against the reference index, train strand and fragment-length models
from unique-mapper fragments, and buffer all resolved fragments into
a memory-efficient columnar buffer (``FragmentBuffer``).

**Count Assignment** (``count_from_buffer``): Iterate the buffer once
to assign unique counts, build CSR EM data (via ``scan.EmDataBuilder``),
construct loci (via ``locus.build_loci``), compute Empirical Bayes gDNA
priors, and run per-locus EM with ``2*n_t + 1`` components.

Scoring functions live in ``scoring.py``.  Locus construction and EM
initialization live in ``locus.py``.  The CSR builder lives in
``scan.py``.  This module is a thin orchestrator.
"""

import logging
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bam import parse_bam_file, detect_sj_strand_tag
from .buffer import FragmentBuffer
from .categories import SpliceType
from .estimator import AbundanceEstimator
from .fragment import Fragment
from .types import ChimeraType, Strand
from .index import HulkIndex
from .frag_length_model import FragmentLengthModels
from .resolution import pair_multimapper_reads, resolve_fragment
from .stats import PipelineStats
from .strand_model import StrandModels

# --- New modular imports ---
from .scoring import (
    ScoringContext,
    overhang_alpha_to_log_penalty,
    DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    DEFAULT_OVERHANG_ALPHA,
    DEFAULT_MISMATCH_ALPHA,
    GDNA_SPLICE_PENALTIES,
    SPLICE_UNSPLICED,
    SPLICE_UNANNOT,
    SPLICE_ANNOT,
    score_gdna_standalone,
)
from .locus import (
    build_loci,
    build_locus_em_data,
    compute_nrna_init,
    compute_gdna_rate_from_strand,
    compute_eb_gdna_priors,
    EB_K_LOCUS,
    EB_K_CHROM,
)
from .scan import EmDataBuilder

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Backward-compatible re-exports
# (tests, scripts, and external tools import these names from pipeline)
# ---------------------------------------------------------------------------
_GDNA_SPLICE_PENALTIES = GDNA_SPLICE_PENALTIES
_SPLICE_UNSPLICED = SPLICE_UNSPLICED
_SPLICE_UNANNOT = SPLICE_UNANNOT
_SPLICE_ANNOT = SPLICE_ANNOT
_DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT = DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT

_score_gdna_candidate = score_gdna_standalone
_build_loci = build_loci
_build_locus_em_data = build_locus_em_data
_compute_nrna_init = compute_nrna_init
_compute_gdna_rate_from_strand = compute_gdna_rate_from_strand
_compute_eb_gdna_priors = compute_eb_gdna_priors
_EB_K_LOCUS = EB_K_LOCUS
_EB_K_CHROM = EB_K_CHROM


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


# ---------------------------------------------------------------------------
# BAM Scan: resolve, train models, buffer
# ---------------------------------------------------------------------------

def scan_and_buffer(
    bam_iter,
    index: HulkIndex,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    log_every: int = 1_000_000,
    max_frag_length: int = 1000,
    chunk_size: int = 1_000_000,
    max_memory_bytes: int = 2 * 1024**3,
    spill_dir: Path | None = None,
    sj_strand_tag: str | tuple[str, ...] = "XS",
    overlap_min_frac: float = 0.99,
    min_spliced_observations: int = 10,
) -> tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]:
    """Single-pass BAM scan: resolve fragments, train models, buffer results.

    Reads the BAM once.  For each fragment:

    1. Build a ``Fragment`` from the read pair.
    2. Resolve against the reference index.
    3. Train strand and fragment-length models (unique mappers only).
    4. Buffer all resolved fragments into the columnar buffer.

    Two strand models are trained (unique mappers only):

    - **exonic_spliced** *(pure RNA)*: from SPLICED_ANNOT fragments
      with unique gene and unambiguous strands.
    - **intergenic** *(~100% gDNA)*: from intergenic fragments,
      POS reference convention.

    Parameters
    ----------
    bam_iter
        Iterator over ``pysam.AlignedSegment`` objects.
    index : HulkIndex
        The loaded reference index.
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include multimapping reads (NH > 1) for counting
        (default False).
    log_every : int
        Log progress every *log_every* read-name groups.
    max_frag_length : int
        Maximum fragment length for the histogram model.
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory for in-memory buffer chunks before disk spill
        (default 2 GiB).
    spill_dir : Path or None
        Directory for spilled chunk files.
    sj_strand_tag : str or tuple of str
        BAM tag(s) for splice-junction strand.

    Returns
    -------
    tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]
    """
    stats = PipelineStats()
    bam_stats = stats.as_bam_stats_dict()

    strand_models = StrandModels()
    strand_models.min_spliced_observations = max(1, int(min_spliced_observations))
    frag_length_models = FragmentLengthModels(max_size=max_frag_length)
    buffer = FragmentBuffer(
        t_to_g_arr=index.t_to_g_arr,
        chunk_size=chunk_size,
        max_memory_bytes=max_memory_bytes,
        spill_dir=spill_dir,
    )

    logger.info("[START] Scanning BAM → resolve + train models + buffer")

    frag_id = 0
    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        bam_iter,
        bam_stats,
        skip_duplicates=skip_duplicates,
        include_multimap=include_multimap,
    ):
        secondary_pairs: list[tuple[list, list]] = []
        if sec_r1 or sec_r2:
            secondary_pairs = pair_multimapper_reads(
                sec_r1, sec_r2, index,
                sj_strand_tag=sj_strand_tag,
                overlap_min_frac=overlap_min_frac,
            )

        all_hits = list(hits) + secondary_pairs
        num_hits = max(nh, len(all_hits))
        is_unique_mapper = num_hits == 1
        n_buffered_mm = 0

        for r1_reads, r2_reads in all_hits:
            frag = Fragment.from_reads(
                r1_reads, r2_reads, sj_strand_tag=sj_strand_tag,
            )
            stats.n_fragments += 1

            if not frag.exons:
                continue

            result = resolve_fragment(
                frag, index, overlap_min_frac=overlap_min_frac,
            )

            if result is None:
                # --- Intergenic: deterministic gDNA assignment ---
                if frag.introns:
                    stats.n_intergenic_spliced += 1
                else:
                    stats.n_intergenic_unspliced += 1
                if is_unique_mapper and not frag.introns:
                    # Only unspliced intergenic fragments are reliable
                    # for training.  Spliced intergenic fragments may be
                    # novel unannotated transcripts or alignment artifacts.
                    flen = frag.genomic_footprint
                    if flen > 0:
                        frag_length_models.observe(flen, splice_type=None)
                        stats.n_frag_length_intergenic += 1
                if is_unique_mapper:
                    intergenic_strand = Strand.NONE
                    for eb in frag.exons:
                        intergenic_strand |= Strand(eb.strand)
                    if intergenic_strand in (Strand.POS, Strand.NEG):
                        strand_models.intergenic.observe(
                            intergenic_strand, Strand.POS,
                        )
                continue

            result.num_hits = num_hits
            result.nm = frag.nm

            # --- Handle chimeric fragments ---
            if result.is_chimeric:
                ct = result.chimera_type
                stats.n_chimeric += 1
                if ct == ChimeraType.TRANS:
                    stats.n_chimeric_trans += 1
                elif ct == ChimeraType.CIS_STRAND_SAME:
                    stats.n_chimeric_cis_strand_same += 1
                elif ct == ChimeraType.CIS_STRAND_DIFF:
                    stats.n_chimeric_cis_strand_diff += 1
                buffer.append(result, frag_id=frag_id)
                continue

            # --- Stats: overlap type ---
            stats.n_with_exon += 1
            if result.splice_type == SpliceType.SPLICED_ANNOT:
                stats.n_with_annotated_sj += 1
            elif result.splice_type == SpliceType.SPLICED_UNANNOT:
                stats.n_with_unannotated_sj += 1

            if result.is_unique_gene:
                stats.n_unique_gene += 1
            else:
                stats.n_multi_gene += 1

            if is_unique_mapper:
                # --- Train models (unique mappers only) ---
                if result.is_strand_qualified:
                    strand_models.exonic_spliced.observe(
                        result.exon_strand, result.sj_strand,
                    )
                    stats.n_strand_trained += 1
                elif result.splice_type != SpliceType.SPLICED_ANNOT:
                    stats.n_strand_skipped_no_sj += 1
                elif not result.is_unique_gene:
                    stats.n_strand_skipped_multi_gene += 1
                else:
                    stats.n_strand_skipped_ambiguous += 1

                # Train exonic fallback model from ALL exonic fragments
                # with unique gene assignment and unambiguous strand.
                # Uses annotated gene strand as truth (not SJ strand).
                if (
                    result.is_unique_gene
                    and result.exon_strand in (Strand.POS, Strand.NEG)
                ):
                    t_idx = next(iter(result.t_inds))
                    gene_idx = index.t_to_g_arr[t_idx]
                    gene_strand = Strand(index.g_to_strand_arr[gene_idx])
                    if gene_strand in (Strand.POS, Strand.NEG):
                        strand_models.exonic.observe(
                            result.exon_strand, gene_strand,
                        )

                # Train fragment-length model.  Only use reads where
                # ALL candidate transcripts agree on a single fragment
                # length (i.e. no ambiguity from different SJ corrections).
                fl_vals = result.frag_lengths
                if fl_vals and len(set(fl_vals.values())) == 1:
                    fl = next(iter(fl_vals.values()))
                    if fl > 0:
                        frag_length_models.observe(fl, result.splice_type)
                        stats.n_frag_length_unambiguous += 1
                    else:
                        stats.n_frag_length_ambiguous += 1
                else:
                    stats.n_frag_length_ambiguous += 1

            # --- Buffer ALL resolved fragments with frag_id ---
            buffer.append(result, frag_id=frag_id)

            if not is_unique_mapper:
                n_buffered_mm += 1

        if n_buffered_mm > 0:
            stats.n_multimapper_groups += 1
            stats.n_multimapper_alignments += n_buffered_mm

        frag_id += 1
        if frag_id % log_every == 0:
            logger.debug(
                f"  {frag_id:,} read-name groups, "
                f"{stats.n_fragments:,} fragments, "
                f"{strand_models.n_observations:,} strand obs, "
                f"{buffer.total_fragments:,} buffered"
            )

    buffer.finalize()

    logger.info(
        f"[DONE] Scan: {stats.n_fragments:,} fragments → "
        f"{stats.n_unique_gene:,} unique-gene, "
        f"{stats.n_multi_gene:,} multi-gene, "
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
# Count from buffer (locus-level EM, no global EM)
# ---------------------------------------------------------------------------

def count_from_buffer(
    buffer: FragmentBuffer,
    index: HulkIndex,
    strand_models: StrandModels,
    frag_length_models: FragmentLengthModels,
    stats: PipelineStats,
    *,
    seed: int | None = None,
    em_pseudocount: float = 0.5,
    em_iterations: int = 1000,
    em_convergence_delta: float = 1e-6,
    log_every: int = 1_000_000,
    gdna_splice_penalties: dict | None = None,
    confidence_threshold: float = 0.95,
    overhang_log_penalty: float | None = None,
    mismatch_log_penalty: float | None = None,
    annotations: "AnnotationTable | None" = None,
) -> AbundanceEstimator:
    """Assign counts from buffered fragments via locus-level EM.

    Architecture per locus: 2*n_t + 1 components
    (mRNA + nRNA + ONE gDNA shadow).

    Spliced fragments NEVER compete for gDNA.
    Only unspliced fragments have a gDNA candidate.

    Parameters
    ----------
    buffer : FragmentBuffer
    index : HulkIndex
    strand_models : StrandModels
    frag_length_models : FragmentLengthModels
    stats : PipelineStats
    seed : int or None
    em_pseudocount : float
    em_iterations : int
    em_convergence_delta : float
    log_every : int
    gdna_splice_penalties : dict or None
    confidence_threshold : float
    overhang_log_penalty : float or None
    mismatch_log_penalty : float or None
    annotations : AnnotationTable or None

    Returns
    -------
    AbundanceEstimator
    """
    # --- Compute geometry ---
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    mean_frag = (
        frag_length_models.global_model.mean
        if frag_length_models.global_model.n_observations > 0
        else 200.0
    )

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

    gene_spans = (
        index.g_df["end"].values - index.g_df["start"].values
    ).astype(np.float64)

    transcript_spans = (
        index.t_df["end"].values - index.t_df["start"].values
    ).astype(np.float64)

    intronic_spans = index.intron_span_per_gene()

    counter = AbundanceEstimator(
        index.num_transcripts, index.num_genes,
        seed=seed, alpha=em_pseudocount,
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        gene_spans=gene_spans,
        mean_frag=mean_frag,
        intronic_spans=intronic_spans,
        transcript_spans=transcript_spans,
    )

    logger.info(
        f"[START] Counting {buffer.total_fragments:,} buffered fragments "
        f"(locus-level EM: mRNA/nRNA/gDNA)"
    )

    # --- Build ScoringContext and scan buffer ---
    ctx = ScoringContext.from_models(
        strand_models, frag_length_models, index, counter,
        overhang_log_penalty=overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        gdna_splice_penalties=gdna_splice_penalties,
    )
    builder = EmDataBuilder(
        ctx, counter, stats, index, strand_models,
        annotations=annotations,
    )
    em_data = builder.scan(buffer, log_every)

    # --- Per-transcript nRNA init from intronic sense excess ---
    nrna_init = compute_nrna_init(
        counter.transcript_intronic_sense,
        counter.transcript_intronic_antisense,
        transcript_spans,
        exonic_lengths,
        mean_frag,
        strand_models,
    )
    counter.nrna_init = nrna_init

    ss = strand_models.strand_specificity
    if ss <= 0.6:
        logger.warning(
            f"[WARN] Strand specificity {ss:.3f} ≤ 0.6: "
            "strand-based gDNA/nRNA estimation unreliable"
        )

    logger.info(
        f"[DONE] Scan: {stats.n_counted_truly_unique:,} unique, "
        f"{em_data.n_units:,} ambiguous units "
        f"({em_data.n_candidates:,} RNA candidates)"
    )

    # --- Locus-based EM ---
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

        gdna_inits = compute_eb_gdna_priors(
            loci, em_data, counter, index, strand_models,
        )

        # Per-locus EM + posterior assignment
        total_gdna_em = 0.0
        _unit_ann_list: list | None = [] if annotations is not None else None
        for i, locus in enumerate(loci):
            if len(locus.unit_indices) == 0:
                continue
            locus_em = build_locus_em_data(
                locus, em_data, counter, index, mean_frag,
                gdna_init=gdna_inits[i],
            )
            theta, alpha = counter.run_locus_em(
                locus_em,
                em_iterations=em_iterations,
                em_convergence_delta=em_convergence_delta,
            )
            gdna_assigned = counter.assign_locus_ambiguous(
                locus_em, theta,
                confidence_threshold=confidence_threshold,
                unit_annotations=_unit_ann_list,
            )
            total_gdna_em += gdna_assigned

            counter.gdna_locus_results.append({
                "locus_id": locus.locus_id,
                "gdna_count": gdna_assigned,
            })

            if gdna_assigned > 0 and len(locus.gene_indices) > 0:
                per_gene = gdna_assigned / len(locus.gene_indices)
                for g_idx in locus.gene_indices:
                    counter.gdna_gene_summary[int(g_idx)] += per_gene

            if len(locus.transcript_indices) > 100:
                logger.debug(
                    f"  Locus {locus.locus_id}: "
                    f"{len(locus.transcript_indices)} transcripts, "
                    f"{len(locus.unit_indices)} units, "
                    f"gDNA={gdna_assigned:.0f}"
                )

        counter._gdna_em_total = total_gdna_em

        # --- Map EM unit annotations to frag_ids ---
        if annotations is not None and _unit_ann_list:
            for (
                global_unit_idx, g_tid, g_gid,
                pool_code, posterior, n_cand,
            ) in _unit_ann_list:
                fid = int(em_data.frag_ids[global_unit_idx])
                fc = int(em_data.frag_class[global_unit_idx])
                st = int(em_data.splice_type[global_unit_idx])
                annotations.add(
                    frag_id=fid,
                    best_tid=g_tid,
                    best_gid=g_gid,
                    pool=pool_code,
                    posterior=posterior,
                    frag_class=fc,
                    n_candidates=n_cand,
                    splice_type=st,
                )

        logger.info(
            f"[DONE] Per-locus EM: {len(loci)} loci, "
            f"{em_data.n_units:,} ambiguous fragments, "
            f"gDNA EM={total_gdna_em:.0f}"
        )
    else:
        logger.info("[SKIP] No ambiguous fragments for EM")

    # --- Update stats ---
    _gdna_em = counter.gdna_em_count
    stats.n_gdna_em = 0 if (math.isnan(_gdna_em) or math.isinf(_gdna_em)) else int(_gdna_em)

    logger.info(
        f"[gDNA] total={counter.gdna_total:.0f} "
        f"(EM={_gdna_em:.0f}), "
        f"contamination rate={counter.gdna_contamination_rate:.2%}"
    )

    return counter


# ---------------------------------------------------------------------------
# Backward-compatible wrapper for _scan_and_build_em_data
# ---------------------------------------------------------------------------

def _scan_and_build_em_data(
    buffer,
    index,
    strand_models,
    frag_length_models,
    counter,
    stats,
    log_every,
    gdna_splice_penalties=None,
    overhang_log_penalty=None,
    mismatch_log_penalty=None,
    annotations=None,
):
    """Backward-compatible wrapper around EmDataBuilder.scan().

    New code should use ``ScoringContext`` + ``EmDataBuilder`` directly.
    This wrapper exists so that ``test_pipeline_routing.py`` and
    profiling scripts can continue importing ``_scan_and_build_em_data``
    from ``pipeline`` without changes.
    """
    ctx = ScoringContext.from_models(
        strand_models, frag_length_models, index, counter,
        overhang_log_penalty=overhang_log_penalty,
        mismatch_log_penalty=mismatch_log_penalty,
        gdna_splice_penalties=gdna_splice_penalties,
    )
    builder = EmDataBuilder(
        ctx, counter, stats, index, strand_models,
        annotations=annotations,
    )
    return builder.scan(buffer, log_every)


# ---------------------------------------------------------------------------
# Full pipeline orchestration
# ---------------------------------------------------------------------------

def run_pipeline(
    bam_path,
    index: HulkIndex,
    *,
    skip_duplicates: bool = True,
    include_multimap: bool = False,
    max_frag_length: int = 1000,
    log_every: int = 1_000_000,
    chunk_size: int = 1_000_000,
    max_memory_bytes: int = 2 * 1024**3,
    spill_dir: Path | str | None = None,
    sj_strand_tag: str | tuple[str, ...] = "auto",
    seed: int | None = None,
    em_pseudocount: float = 0.5,
    em_iterations: int = 1000,
    em_convergence_delta: float = 1e-6,
    gdna_splice_penalty_unannot: float = DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    confidence_threshold: float = 0.95,
    overlap_min_frac: float = 0.99,
    overhang_alpha: float | None = None,
    overhang_log_penalty: float | None = None,
    mismatch_alpha: float | None = None,
    mismatch_log_penalty: float | None = None,
    min_spliced_observations: int = 10,
    annotated_bam_path: str | Path | None = None,
) -> PipelineResult:
    """Run the complete counting pipeline with locus-level EM.

    Opens the BAM **once**, scans all fragments, trains models,
    buffers resolved fragments, then counts via per-locus EM.

    Parameters
    ----------
    bam_path : str or Path
        Path to the name-sorted / collated BAM file.
    index : HulkIndex
        The loaded reference index.
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include multimapping reads (default False).
    max_frag_length : int
        Maximum fragment length for histogram models.
    log_every : int
        Log progress every *log_every* read-name groups.
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory before disk spill (default 2 GiB).
    spill_dir : Path, str, or None
        Directory for disk-spilled buffer chunks.
    sj_strand_tag : str, tuple of str, or ``"auto"``
        BAM tag(s) for splice-junction strand.
    seed : int or None
        Random seed for reproducibility.
    em_pseudocount : float
        Dirichlet prior pseudocount (default 0.5).
    em_iterations : int
        Maximum number of EM iterations (default 1000).
    em_convergence_delta : float
        Convergence threshold for EM theta updates (default 1e-6).
    gdna_splice_penalty_unannot : float
        gDNA splice penalty for SPLICED_UNANNOT fragments.
    confidence_threshold : float
        Posterior threshold for high-confidence assignment.
    min_spliced_observations : int
        Minimum spliced observations before trusting pure model.
    annotated_bam_path : str, Path, or None
        If provided, write an annotated BAM with per-fragment assignment
        tags to this path (second BAM pass).

    Returns
    -------
    PipelineResult
    """
    import pysam

    bam_path = str(bam_path)
    if spill_dir is not None:
        spill_dir = Path(spill_dir)

    # -- Resolve sj_strand_tag --
    if sj_strand_tag == "auto":
        detected = detect_sj_strand_tag(bam_path)
        if len(detected) == 1:
            sj_strand_tag = detected[0]
        elif len(detected) > 1:
            sj_strand_tag = detected
        else:
            sj_strand_tag = ()

    # -- Single BAM pass --
    bamfh = pysam.AlignmentFile(bam_path, "rb")
    try:
        stats, strand_models, frag_length_models, buffer = (
            scan_and_buffer(
                bamfh.fetch(until_eof=True),
                index,
                skip_duplicates=skip_duplicates,
                include_multimap=include_multimap,
                log_every=log_every,
                max_frag_length=max_frag_length,
                chunk_size=chunk_size,
                max_memory_bytes=max_memory_bytes,
                spill_dir=spill_dir,
                sj_strand_tag=sj_strand_tag,
                overlap_min_frac=overlap_min_frac,
                min_spliced_observations=min_spliced_observations,
            )
        )
    finally:
        bamfh.close()

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
    frag_length_models.finalize()

    # -- Count assignment via locus-level EM --
    gdna_splice_penalties = dict(GDNA_SPLICE_PENALTIES)
    gdna_splice_penalties[SPLICE_UNANNOT] = gdna_splice_penalty_unannot

    if overhang_log_penalty is None and overhang_alpha is not None:
        overhang_log_penalty = overhang_alpha_to_log_penalty(overhang_alpha)
    if mismatch_log_penalty is None and mismatch_alpha is not None:
        mismatch_log_penalty = overhang_alpha_to_log_penalty(mismatch_alpha)

    # -- Annotation table for second BAM pass (opt-in) --
    annotations = None
    if annotated_bam_path is not None:
        from .annotate import AnnotationTable
        annotations = AnnotationTable.create(
            capacity=max(buffer.total_fragments + 1024, 4096)
        )

    try:
        counter = count_from_buffer(
            buffer,
            index,
            strand_models,
            frag_length_models,
            stats,
            seed=seed,
            em_pseudocount=em_pseudocount,
            em_iterations=em_iterations,
            em_convergence_delta=em_convergence_delta,
            log_every=log_every,
            gdna_splice_penalties=gdna_splice_penalties,
            confidence_threshold=confidence_threshold,
            overhang_log_penalty=overhang_log_penalty,
            mismatch_log_penalty=mismatch_log_penalty,
            annotations=annotations,
        )
    finally:
        buffer.cleanup()

    # -- Second BAM pass: write annotated BAM (opt-in) --
    if annotated_bam_path is not None and annotations is not None:
        from .annotate import write_annotated_bam
        write_annotated_bam(
            bam_path,
            str(annotated_bam_path),
            annotations,
            index,
            skip_duplicates=skip_duplicates,
            include_multimap=include_multimap,
            sj_strand_tag=sj_strand_tag,
            overlap_min_frac=overlap_min_frac,
        )

    return PipelineResult(
        stats=stats,
        strand_models=strand_models,
        frag_length_models=frag_length_models,
        estimator=counter,
    )
