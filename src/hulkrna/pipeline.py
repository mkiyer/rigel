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
from dataclasses import dataclass, replace as _replace
from pathlib import Path

import numpy as np

# Transcript count above which a locus triggers debug logging.
_LARGE_LOCUS_LOG_THRESHOLD = 100
# Padding and minimum capacity for the annotation table.
_ANNOTATION_TABLE_PADDING = 1024
_ANNOTATION_TABLE_MIN_CAPACITY = 4096

import os

from .bam import parse_bam_file, detect_sj_strand_tag
from .buffer import FragmentBuffer, _FinalizedChunk
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
from .scoring import ScoringContext
from .locus import (
    build_loci,
    build_locus_em_data,
    compute_nrna_init,
    compute_gdna_rate_from_strand,
    compute_eb_gdna_priors,
)
from .scan import EmDataBuilder
from .config import (
    EMConfig,
    PipelineConfig,
    ScanConfig,
    ScoringConfig,
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


# ---------------------------------------------------------------------------
# Native BAM scanner availability
# ---------------------------------------------------------------------------

_HAS_NATIVE_BAM_SCANNER = False
try:
    from ._bam_impl import BamScanner as _NativeBamScanner  # type: ignore[import-not-found]
    from ._bam_impl import detect_sj_strand_tag as _native_detect_sj_tag  # type: ignore[import-not-found]
    _HAS_NATIVE_BAM_SCANNER = True
except ImportError:
    pass

# Force Python fallback via env var (useful for testing / debugging).
if os.environ.get("HULKRNA_FORCE_PYTHON_BAM", "").strip().lower() in (
    "1", "true", "yes",
):
    _HAS_NATIVE_BAM_SCANNER = False


def _sj_tag_to_spec(sj_strand_tag) -> str:
    """Convert ScanConfig.sj_strand_tag to the string spec for BamScanner."""
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
    # Exonic-spliced model
    obs = strand_dict.get("exonic_spliced_obs", [])
    truth = strand_dict.get("exonic_spliced_truth", [])
    for o, t in zip(obs, truth):
        strand_models.exonic_spliced.observe(Strand(int(o)), Strand(int(t)))

    # Exonic fallback model
    obs = strand_dict.get("exonic_obs", [])
    truth = strand_dict.get("exonic_truth", [])
    for o, t in zip(obs, truth):
        strand_models.exonic.observe(Strand(int(o)), Strand(int(t)))

    # Intergenic model
    obs = strand_dict.get("intergenic_obs", [])
    truth = strand_dict.get("intergenic_truth", [])
    for o, t in zip(obs, truth):
        strand_models.intergenic.observe(Strand(int(o)), Strand(int(t)))


def _replay_fraglen_observations(
    fraglen_dict: dict,
    frag_length_models: FragmentLengthModels,
) -> None:
    """Replay C++ fragment-length observation arrays into Python models."""
    lengths = fraglen_dict.get("lengths", [])
    splice_types = fraglen_dict.get("splice_types", [])
    for fl, st in zip(lengths, splice_types):
        frag_length_models.observe(int(fl), SpliceType(int(st)))

    intergenic_lengths = fraglen_dict.get("intergenic_lengths", [])
    for fl in intergenic_lengths:
        frag_length_models.observe(int(fl), splice_type=None)


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
    stats.n_unique_gene = stats_dict.get("n_unique_gene", 0)
    stats.n_multi_gene = stats_dict.get("n_multi_gene", 0)
    stats.n_strand_trained = stats_dict.get("n_strand_trained", 0)
    stats.n_strand_skipped_no_sj = stats_dict.get(
        "n_strand_skipped_no_sj", 0)
    stats.n_strand_skipped_multi_gene = stats_dict.get(
        "n_strand_skipped_multi_gene", 0)
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


def scan_and_buffer_native(
    bam_path: str,
    index: HulkIndex,
    scan: "ScanConfig",
) -> tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]:
    """Native C++ BAM scan: resolve, train models, buffer — all in one pass.

    This is the fast path that replaces ``scan_and_buffer`` when the
    ``_bam_impl`` C++ extension is available.  The entire BAM parsing,
    fragment construction, overlap resolution, model training and
    columnar buffering happens in C++ via htslib.

    Parameters
    ----------
    bam_path : str
        Path to the name-sorted / collated BAM file.
    index : HulkIndex
        The loaded reference index.
    scan : ScanConfig
        BAM scanning and buffering configuration.

    Returns
    -------
    tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]
    """
    stats = PipelineStats()
    strand_models = StrandModels()
    strand_models.min_spliced_observations = max(
        1, int(scan.min_spliced_observations),
    )
    frag_length_models = FragmentLengthModels(max_size=scan.max_frag_length)
    buffer = FragmentBuffer(
        t_to_g_arr=index.t_to_g_arr,
        chunk_size=scan.chunk_size,
        max_memory_bytes=scan.max_memory_bytes,
        spill_dir=scan.spill_dir,
    )

    logger.info("[START] Native C++ BAM scan → resolve + train + buffer")

    # Provide gene strand info for exonic fallback strand model training
    resolve_ctx = index._resolve_ctx
    resolve_ctx.set_gene_strands(index.g_to_strand_arr.tolist())

    # Create the native scanner
    sj_spec = _sj_tag_to_spec(scan.sj_strand_tag)
    scanner = _NativeBamScanner(
        resolve_ctx,
        sj_spec,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
    )

    # Execute the full BAM scan in C++
    result = scanner.scan(bam_path)

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
        raw = scanner.finalize_accumulator(index.t_to_g_arr.tolist())
        size = raw["size"]

        chunk = _FinalizedChunk(
            splice_type=np.frombuffer(
                raw["splice_type"], dtype=np.uint8).copy(),
            exon_strand=np.frombuffer(
                raw["exon_strand"], dtype=np.uint8).copy(),
            sj_strand=np.frombuffer(
                raw["sj_strand"], dtype=np.uint8).copy(),
            num_hits=np.frombuffer(
                raw["num_hits"], dtype=np.uint16).copy(),
            merge_criteria=np.frombuffer(
                raw["merge_criteria"], dtype=np.uint8).copy(),
            chimera_type=np.frombuffer(
                raw["chimera_type"], dtype=np.uint8).copy(),
            t_offsets=np.frombuffer(
                raw["t_offsets"], dtype=np.int64).copy(),
            t_indices=np.frombuffer(
                raw["t_indices"], dtype=np.int32).copy(),
            frag_lengths=np.frombuffer(
                raw["frag_lengths"], dtype=np.int32).copy(),
            exon_bp=np.frombuffer(
                raw["exon_bp"], dtype=np.int32).copy(),
            intron_bp=np.frombuffer(
                raw["intron_bp"], dtype=np.int32).copy(),
            unambig_intron_bp=np.frombuffer(
                raw["unambig_intron_bp"], dtype=np.int32).copy(),
            n_genes=np.frombuffer(
                raw["n_genes"], dtype=np.uint8).copy(),
            frag_id=np.frombuffer(
                raw["frag_id"], dtype=np.int64).copy(),
            read_length=np.frombuffer(
                raw["read_length"], dtype=np.uint32).copy(),
            genomic_footprint=np.frombuffer(
                raw["genomic_footprint"], dtype=np.int32).copy(),
            genomic_start=np.frombuffer(
                raw["genomic_start"], dtype=np.int32).copy(),
            nm=np.frombuffer(
                raw["nm"], dtype=np.uint16).copy(),
            size=size,
        )

        buffer._chunks.append(chunk)
        buffer._total_size += size
        buffer._memory_bytes += chunk.memory_bytes

        # Spill if over memory budget
        if buffer.max_memory_bytes > 0:
            while buffer._memory_bytes > buffer.max_memory_bytes:
                if not buffer._spill_oldest():
                    break

    logger.info(
        f"[DONE] Native scan: {stats.n_fragments:,} fragments → "
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
# BAM Scan: resolve, train models, buffer
# ---------------------------------------------------------------------------

def scan_and_buffer(
    bam_iter,
    index: HulkIndex,
    scan: ScanConfig,
) -> tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]:
    """Single-pass BAM scan: resolve fragments, train models, buffer results.

    Parameters
    ----------
    bam_iter
        Iterator over ``pysam.AlignedSegment`` objects.
    index : HulkIndex
        The loaded reference index.
    scan : ScanConfig
        BAM scanning and buffering configuration.

    Returns
    -------
    tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]
    """
    stats = PipelineStats()
    bam_stats = stats.as_bam_stats_dict()

    strand_models = StrandModels()
    strand_models.min_spliced_observations = max(
        1, int(scan.min_spliced_observations)
    )
    frag_length_models = FragmentLengthModels(max_size=scan.max_frag_length)
    buffer = FragmentBuffer(
        t_to_g_arr=index.t_to_g_arr,
        chunk_size=scan.chunk_size,
        max_memory_bytes=scan.max_memory_bytes,
        spill_dir=scan.spill_dir,
    )

    logger.info("[START] Scanning BAM → resolve + train models + buffer")

    frag_id = 0
    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        bam_iter,
        bam_stats,
        skip_duplicates=scan.skip_duplicates,
        include_multimap=scan.include_multimap,
    ):
        secondary_pairs: list[tuple[list, list]] = []
        if sec_r1 or sec_r2:
            secondary_pairs = pair_multimapper_reads(
                sec_r1, sec_r2, index,
                sj_strand_tag=scan.sj_strand_tag,
            )

        all_hits = list(hits) + secondary_pairs
        num_hits = max(nh, len(all_hits))
        is_unique_mapper = num_hits == 1
        n_buffered_mm = 0

        for r1_reads, r2_reads in all_hits:
            frag = Fragment.from_reads(
                r1_reads, r2_reads, sj_strand_tag=scan.sj_strand_tag,
            )
            stats.n_fragments += 1

            if not frag.exons:
                continue

            result = resolve_fragment(frag, index)

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
                    t_idx = (
                        result.first_t_ind
                        if hasattr(result, 'first_t_ind')
                        else next(iter(result.t_inds))
                    )
                    gene_idx = index.t_to_g_arr[t_idx]
                    gene_strand = Strand(index.g_to_strand_arr[gene_idx])
                    if gene_strand in (Strand.POS, Strand.NEG):
                        strand_models.exonic.observe(
                            result.exon_strand, gene_strand,
                        )

                # Train fragment-length model.  Only use reads where
                # ALL candidate transcripts agree on a single fragment
                # length (i.e. no ambiguity from different SJ corrections).
                _ufl = getattr(result, 'unique_frag_length', None)
                if _ufl is not None:
                    # Native C++ ResolvedResult path
                    if _ufl > 0:
                        frag_length_models.observe(_ufl, result.splice_type)
                        stats.n_frag_length_unambiguous += 1
                    else:
                        stats.n_frag_length_ambiguous += 1
                else:
                    # Python ResolvedFragment fallback
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
        if frag_id % scan.log_every == 0:
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
    em_config: EMConfig | None = None,
    scoring: ScoringConfig | None = None,
    log_every: int = 1_000_000,
    annotations: "AnnotationTable | None" = None,
) -> AbundanceEstimator:
    """Assign counts from buffered fragments via locus-level EM.

    Parameters
    ----------
    buffer : FragmentBuffer
    index : HulkIndex
    strand_models : StrandModels
    frag_length_models : FragmentLengthModels
    stats : PipelineStats
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    scoring : ScoringConfig or None
        Scoring penalty configuration.  Defaults to ``ScoringConfig()``.
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
        scoring = ScoringConfig()

    # --- Compute geometry ---
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

    gene_spans = (
        index.g_df["end"].values - index.g_df["start"].values
    ).astype(np.float64)

    transcript_spans = (
        index.t_df["end"].values - index.t_df["start"].values
    ).astype(np.float64)

    intronic_spans = index.intron_span_per_gene()

    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        gene_spans=gene_spans,
        mean_frag=mean_frag,
        intronic_spans=intronic_spans,
        transcript_spans=transcript_spans,
    )

    counter = AbundanceEstimator(
        index.num_transcripts, index.num_genes,
        em_config=em_config,
        geometry=geometry,
    )

    logger.info(
        f"[START] Counting {buffer.total_fragments:,} buffered fragments "
        f"(locus-level EM: mRNA/nRNA/gDNA)"
    )

    # --- Build ScoringContext and scan buffer ---
    ctx = ScoringContext.from_models(
        strand_models, frag_length_models, index, counter,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
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
        geometry.transcript_spans,
        geometry.exonic_lengths,
        geometry.mean_frag,
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
        f"[DONE] Scan: {stats.deterministic_unique_units:,} unique, "
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
                locus, em_data, counter, index, geometry.mean_frag,
                gdna_init=gdna_inits[i],
            )
            theta, alpha = counter.run_locus_em(
                locus_em,
                em_iterations=em_config.iterations,
                em_convergence_delta=em_config.convergence_delta,
            )
            gdna_assigned = counter.assign_locus_ambiguous(
                locus_em, theta,
                confidence_threshold=em_config.confidence_threshold,
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

            if len(locus.transcript_indices) > _LARGE_LOCUS_LOG_THRESHOLD:
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
# Full pipeline orchestration
# ---------------------------------------------------------------------------

def _scan_with_pysam(
    bam_path: str,
    index: HulkIndex,
    scan: "ScanConfig",
) -> tuple[PipelineStats, StrandModels, FragmentLengthModels, FragmentBuffer]:
    """Python/pysam BAM scan fallback."""
    import pysam

    bamfh = pysam.AlignmentFile(bam_path, "rb")
    try:
        return scan_and_buffer(bamfh.fetch(until_eof=True), index, scan)
    finally:
        bamfh.close()


def run_pipeline(
    bam_path,
    index: HulkIndex,
    config: PipelineConfig | None = None,
) -> PipelineResult:
    """Run the complete counting pipeline with locus-level EM.

    Opens the BAM **once**, scans all fragments, trains models,
    buffers resolved fragments, then counts via per-locus EM.

    When the native C++ BAM scanner (``_bam_impl``) is available, the
    entire BAM scan runs in C++ via htslib for dramatically better
    performance.  If the native scanner is not available or fails, the
    function falls back to the Python/pysam path automatically.  Set
    ``HULKRNA_FORCE_PYTHON_BAM=1`` to force the fallback.

    Parameters
    ----------
    bam_path : str or Path
        Path to the name-sorted / collated BAM file.
    index : HulkIndex
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
    scan = config.scan
    if scan.sj_strand_tag == "auto":
        detected = detect_sj_strand_tag(bam_path)
        if len(detected) == 1:
            resolved_tag = detected[0]
        elif len(detected) > 1:
            resolved_tag = detected
        else:
            resolved_tag = ()
        scan = _replace(scan, sj_strand_tag=resolved_tag)

    # -- Single BAM pass (native C++ or Python fallback) --
    if _HAS_NATIVE_BAM_SCANNER:
        try:
            stats, strand_models, frag_length_models, buffer = (
                scan_and_buffer_native(bam_path, index, scan)
            )
        except Exception:
            logger.warning(
                "Native BAM scanner failed; falling back to Python/pysam",
                exc_info=True,
            )
            stats, strand_models, frag_length_models, buffer = (
                _scan_with_pysam(bam_path, index, scan)
            )
    else:
        stats, strand_models, frag_length_models, buffer = (
            _scan_with_pysam(bam_path, index, scan)
        )

    # -- Finalize models: cache derived values for fast scoring --
    strand_models.finalize()
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
        counter = count_from_buffer(
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
        estimator=counter,
    )
