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

**Quantification** (``quant_from_buffer``): score buffered fragments,
construct MultiLoci, assemble the per-locus gDNA/RNA split prior from the
calibration result (``calibration.priors.assemble_priors``), partition the
CSR data, and run locus-level native EM.

Scoring functions live in ``scoring.py``.  Locus construction and EM
initialization live in ``locus.py``.  The CSR builder lives in
``scan.py``.  This module is a thin orchestrator.
"""

from __future__ import annotations

import logging
import time
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
from .frag_length_model import FragmentLengthModel
from .index import TranscriptIndex
from .native import BamScanner as _NativeBamScanner
from .native import detect_sj_strand_tag as _native_detect_sj_tag
from .scan import FragmentRouter
from .scoring import FragmentScorer
from .splice import SpliceType, census_field
from .stats import PipelineStats
from .strand_model import StrandModels

if TYPE_CHECKING:
    from .annotate import AnnotationTable
    from .calibration import CalibrationResult
    from .calibration.fl import FLModels
    from .calibration.region_arrays import RegionArrays
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
    estimator: AbundanceEstimator
    pipeline_config: "PipelineConfig" = None
    calibration: "CalibrationResult" = None
    # Library-wide FL distributions actually built + used for scoring/calibration
    # (global / RNA / gDNA); the empirical views feed the summary QC report.
    fl_models: "FLModels" = None
    # Per-region genome-wide gDNA track (ref/start/end + gDNA mass/density/frac),
    # built from the calibration result; feeds the report's genome track + bedGraph.
    calibration_track: "object" = None
    # Report-facing calibration diagnostics (the fitted gDNA-density KDE); None
    # when the Phase-2 KDE was not fit (e.g. tiny / toy scenarios).
    calibration_diagnostics: "object" = None


def _sj_tag_to_spec(sj_strand_tag) -> str:
    """Convert BamScanConfig.sj_strand_tag to the string spec for BamScanner."""
    if isinstance(sj_strand_tag, str):
        return sj_strand_tag if sj_strand_tag else "none"
    if isinstance(sj_strand_tag, (list, tuple)):
        return ",".join(sj_strand_tag) if sj_strand_tag else "none"
    return "none"


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
        # Multimapper
        "n_multimapper_groups",
        "n_multimapper_alignments",
        # Splice-artifact blacklist
        "n_sj_observed",
        "n_sj_blacklisted",
        # The deposit hold-out that closes the census identity
        "n_deposit_not_offered",
    ):
        setattr(stats, key, stats_dict.get(key, 0))

    # ⭐ The splice census, enumerated from the enum rather than listed by hand — the whole point of
    # deriving the key name on both sides is that no third list can fall out of step with it.
    for stype in SpliceType:
        setattr(stats, census_field(stype), stats_dict.get(census_field(stype), 0))


def _warn_if_calibration_strand_unidentifiable(strand_models: StrandModels) -> None:
    """Warn when calibration cannot identify strand from spliced RNA evidence."""
    from .calibration.strand_summary import StrandSummary

    primary = StrandSummary.from_model(strand_models.exonic_spliced)
    if primary.is_identifiable():
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
    if diagnostic.is_identifiable():
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
    FragmentBuffer,
    "object",  # AccumulatorPayload | None
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
        ``(stats, strand_models, buffer, calibration_payload)`` — the last
        is ``None`` if the index has no ``regions.feather`` (legacy indexes
        pre-v3).

    ⭐ **No fragment-length model comes out of here.** The scanner used to train its own histogram
    during this pass; it measured length by two rules that were neither each other nor the
    accumulator's ``L``, over a population that was never stated. TRAPS: pure-and-length-censored deleted it. Every
    fragment-length distribution the tool uses is now built from the payload by
    :func:`rigel.calibration.fl.build_fl_models`.
    """
    stats = PipelineStats()
    buffer = FragmentBuffer(
        t_strand_arr=index.t_to_strand_arr,
        chunk_size=scan.fragments_per_chunk,
        max_memory_bytes=scan.buffer_size_bytes,
        spill_dir=scan.spill_dir,
    )
    logger.info("[START] Native C++ BAM scan → resolve + train + buffer")

    # Per-transcript strands drive the all-exonic diagnostic strand model.
    # (A parallel per-GENE array was pushed here too; the resolver stored it and
    # never read it, so both the setter and this full-length list build were
    # deleted 2026-07-28.)
    resolve_ctx = index.resolver
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
    # ⭐ `detect_chimera` needs the library's fragment-length limit to tell an ordinary genomic
    # molecule (contiguous, spanning whatever transcripts lie under it) from a real rearrangement.
    resolve_ctx.set_max_fragment_length(int(scan.max_frag_length))

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

    # Calibration region wiring: install the per-genome v8 region partition into the scanner so
    # per-region evidence is collected during the scan.
    #
    # ⚠ This used to be `getattr(index, "region_df", None)` — a SILENT skip. An index that could not
    # supply a partition disabled calibration with no error anywhere, and the whole pipeline then ran
    # as though the library had no gDNA. The two cases are now separated: a MISSING graph is a broken
    # index and raises; an EMPTY one (every reference of length 0) is a degenerate genome with
    # genuinely nothing to deposit, and is skipped exactly as before.
    # A genome whose references are all zero-length has no regions and nothing to deposit; the loader
    # guarantees the graph itself is present.
    if len(index.regions_df) > 0:
        _wire_calibration_regions(
            scanner,
            index,
            scan.max_frag_length,
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

    # Build the strand models from the scanner's observations. Immutable and built
    # once — the spliced 2×2 is the marginal of its per-sj SJ strand table.
    strand_models = StrandModels.from_scan(result["strand_observations"])
    # ONE strand-qualified fragment credits ONE sj. That is what makes the 2×2
    # exactly the table's marginal, so both halves of the RNA strand Beta-Binomial are
    # fitted on one population. The
    # C++ counts the fragments independently of the table it builds, so this is a real
    # cross-check, and it is the invariant to break if the crediting rule ever changes.
    n_credited = strand_models.sj_table.n_observations
    if n_credited != stats.n_strand_trained:
        raise RuntimeError(
            f"SJ strand table credited {n_credited:,} observations but the scanner "
            f"qualified {stats.n_strand_trained:,} fragments; one qualified fragment must "
            "credit exactly one sj."
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

    # Calibration payload — Phase B: build the AccumulatorPayload from
    # the C++ scanner's fractional-accumulator results.
    from .scan_payload import AccumulatorPayload

    if result.get("calibration") is not None:
        # ⚠ The provenance covers regions AND boundaries. The payload is boundary-keyed by construction — its
        # sj axis is meaningless against a different sj CSR — and `partition_hash` covers
        # `regions.feather` only, deliberately.
        calibration_payload = AccumulatorPayload.from_scan_result(
            result, graph_hash=index.graph_hash
        )
    else:
        calibration_payload = None

    return stats, strand_models, buffer, calibration_payload


def _drain_side_buffer(
    payload, index: TranscriptIndex, strand_models, *, seed: int, _lift: dict | None = None
):
    """⭐ **THE SECOND PASS.** Score the held fragments, draw one hypothesis each, re-deposit.

         (where this sits), §5 (the draw), §6 (the drain)

    Pass 1 holds every fragment whose unsequenced gap has more than one surviving explanation — 2–3.5 % of
    a library, and systematically the **long** ones, because a longer gap admits more hypotheses. Nothing
    is lost, but nothing is tallied either until this runs.

    ⭐ **IT RUNS HERE, BETWEEN THE SCAN AND CALIBRATION, AND THAT IS THE STRUCTURAL DECISION** (§2). Every
    input the score needs — the densities at each object, the fragment-length models, the strand model —
    comes from pass 1 alone, so no calibration output is required. Which is what lets **calibration run
    exactly once, on the complete tally**, instead of once before the drain and again after it.

    ⚠ **The models the SCORER uses are pass one's; the models CALIBRATION uses are the drained tally's.**
    That ordering is §7.1's no-iteration rule made concrete: fit once, score once, drain once, stop. The
    confident set is biased short, so feeding the drained anchor back into the score would prefer the
    shorter — that is, the more-spliced — path, and that loop can run away.

    ⭐ **``_lift`` — the out-parameter an origin-split ORACLE needs, and why it is not a second
    implementation.** TRAPS: draining-breaks-the-oracle: the drain conditions on the WHOLE tally, so partitions drained
    independently do not sum to the whole drained. `second_pass.lift_choices` repairs that by replaying
    the whole's already-drawn choices inside each partition — which needs the choices, the undrained
    whole they were drawn on, and the two index-derived arrays `drain` takes. All four exist only inside
    this function, so it publishes them into ``_lift`` rather than letting a caller re-derive them and
    drift (TRAPS: a-test-that-redefines). Same convention as ``calibrate(_debug=)`` / ``solve_chain(_capture=)``, and
    inert in production, where nobody passes it. ⚠ An empty side buffer leaves ``_lift`` UNTOUCHED — the
    early return below is the "nothing was drained" signal on this path too.
    """
    from .calibration.fl import build_fl_models
    from .calibration.gdna_opportunity import gdna_opportunity_from_index
    from .calibration.sj_opportunity import crossing_probability_from_index
    from .calibration.splice_graph import build_sj_arrays, build_region_partition_arrays
    from .second_pass import choose_hypotheses, drain, score_held_fragments

    held = payload.deferred.n_fragments
    if held == 0:
        # ⚠ Not an error and not a no-op to hide: a library with no annotated intron in any mate gap is a
        # real state. Left undrained, so `payload.drain is None` continues to mean "pass one only".
        logger.info("[SP2] nothing held; the side buffer is empty and the drain is skipped")
        return payload

    start = time.perf_counter()
    _region_bounds, _offsets, region_types = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)
    scores = score_held_fragments(
        payload,
        # ⭐ The SAME de-tilted RNA pool the calibrator will read. The scorer weighs a candidate
        # length by `f(L)`, so handing it the tilted pool would make it prefer the longer hypothesis
        # for the same reason the pool is long in the first place — one definition of the RNA length
        # distribution, or the second pass and the calibration disagree about the library.
        fl_models=build_fl_models(
            payload,
            sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
            gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
        ),
        # ⚠ `P(align_strand agrees | RNA)`, and on an R1-antisense (dUTP) library — which real cfRNA is —
        # this is ≈ 0.01, so DISAGREEMENT is the likely case.
        rna_sense_frac=strand_models.p_r1_sense,
        region_types=region_types,
        sj=sj,
    )
    choices = choose_hypotheses(scores, payload, seed=seed)
    drained = drain(payload, choices, region_types=region_types, sj=sj)
    if _lift is not None:
        # ⛔ ``undrained`` is the payload as it entered here, NOT ``drained``: the drained bank is empty by
        # design ("after it nothing is held"), so it cannot supply `lift_choices`' key pool.
        _lift.update(undrained=payload, choices=choices, region_types=region_types, sj=sj)
    report = drained.drain
    logger.info(
        "[SP2] drained %d held fragments in %.1f s: %d deposited, %d dropped "
        "(%d chose the genomic path, %d a spliced one); %d were undecided and drawn uniformly",
        report.offered,
        time.perf_counter() - start,
        report.deposited,
        report.dropped,
        report.chose_genomic,
        report.chose_spliced,
        scores.n_undecided,
    )
    return drained


def _wire_calibration_regions(
    scanner,
    index: TranscriptIndex,
    max_frag_length: int,
) -> None:
    """Install the index's v8 splice graph into a native BamScanner: the region partition, then the sj.

    :func:`~rigel.calibration.splice_graph.build_region_partition_arrays` flattens the per-reference
    partition into the ``(region_bounds, ref_region_bound_offsets, n_refs, region_types, max_length)`` ABI, and
    :func:`~rigel.calibration.splice_graph.build_sj_arrays` builds the sj CSR keyed by the
    flat region_bound index. Both are derived from ``index.regions_df``/``index.edges_df`` and ordered to match
    ``index.ref_names``, which is the resolver's reference-id space.

    ⚠ **Two calls, and both are required.** ``set_regions`` refuses to run twice, which is why the sj
    are separate; and ``scan`` refuses to run if the second call is missing, because a missing sj table
    is invisible — every observed intron would simply read as unannotated, so all 404,168 sj boundaries and
    both spliced banks would come back empty from a scan that looked perfectly well-formed.
    """
    from .calibration.splice_graph import build_sj_arrays, build_region_partition_arrays

    region_bounds, ref_region_bound_offsets, region_types = build_region_partition_arrays(index)
    n_refs = len(index.ref_names)
    scanner.set_regions(
        np.ascontiguousarray(region_bounds, dtype=np.int64),
        np.ascontiguousarray(ref_region_bound_offsets, dtype=np.int64),
        n_refs,
        np.ascontiguousarray(region_types, dtype=np.uint8),
        int(max_frag_length),
    )
    # ⚠ ``edge_row`` is deliberately NOT passed. It is a join key back to ``index.edges_df``, not the
    # sj-boundary id — that IS the CSR slot — and using it to index a sj bank would write 1.04 M
    # rows past the end of a 404,168-entry array.
    sj = build_sj_arrays(index)
    scanner.set_sj(
        np.ascontiguousarray(sj.offsets, dtype=np.int64),
        np.ascontiguousarray(sj.boundary_right, dtype=np.int64),
        np.ascontiguousarray(sj.strand, dtype=np.int8),
    )


# ---------------------------------------------------------------------------
# Quantify from buffer (locus-level EM, no global EM)
# ---------------------------------------------------------------------------


def _setup_geometry_and_estimator(
    index: TranscriptIndex,
    rna_fl,
    em_config: EMConfig,
    calibration: "CalibrationResult | None" = None,
    region_arrays: "RegionArrays | None" = None,
) -> tuple["TranscriptGeometry", AbundanceEstimator]:
    """Compute transcript geometry and create the AbundanceEstimator.

    When ``calibration`` + ``region_arrays`` are given, the EM effective lengths
    (``effective_lengths_em``) are capture-contracted by the per-region gDNA enrichment over each
    transcript's region set (exon regions for mRNA, full span for nRNA) — so mRNA/nRNA compete with the
    gDNA component on equal (contracted) footing under capture. The output ``effective_lengths`` (used for
    TPM) stay the full FL-marginal length. Uniform enrichment (capture-off) ⇒ ``effective_lengths_em`` ==
    ``effective_lengths`` (bit-identical). See ``calibration/capture_eff_length.py``.
    """
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    if rna_fl.n_observations > 0:
        effective_lengths = rna_fl.compute_all_transcript_eff_lens(
            exonic_lengths.astype(np.int64),
        )
    else:
        effective_lengths = np.maximum(exonic_lengths - _DEFAULT_MEAN_FRAG + 1.0, 1.0)

    effective_lengths_em = None
    if calibration is not None and region_arrays is not None:
        from .calibration.capture_eff_length import transcript_capture_eff_lengths

        effective_lengths_em = transcript_capture_eff_lengths(
            calibration, region_arrays, index, effective_lengths
        )

    transcript_spans = (index.t_df["end"].values - index.t_df["start"].values).astype(np.float64)

    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        transcript_spans=transcript_spans,
        effective_lengths_em=effective_lengths_em,
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
    # True per-fragment locus-subproblem id (so in-locus gDNA winners keep their
    # locus in the ZL tag, instead of collapsing to -1 via the winning transcript).
    locus_ids = np.concatenate(
        [np.full(p.frag_ids.shape[0], p.locus_id, dtype=np.int32) for p in batch_parts]
    )

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
        locus_ids=locus_ids,
    )


def _run_locus_em_partitioned(
    estimator: AbundanceEstimator,
    partitions: dict,
    multi_loci: list,
    index: TranscriptIndex,
    gdna_prior_count: np.ndarray,
    rna_prior_count: np.ndarray,
    gdna_eff_len: np.ndarray,
    *,
    em_config: EMConfig,
    rna_prior_weight: np.ndarray | None = None,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> None:
    """Run the per-locus batch EM and record lean per-locus results.

    Packs each ``LocusPartition`` into the C++ 9-tuple, runs the whole batch in
    one ``run_batch_locus_em_partitioned`` call (the solver is OpenMP-parallel
    internally), and appends a per-locus dict to ``estimator.locus_results`` for
    ``get_loci_df``. The calibration prior enters as the two per-locus alpha
    scalars; ``enable_gdna`` is the structural eligibility (any unspliced unit
    carrying a finite gDNA log-lik — the rule the C++ extractor uses).

    ⭐ ``rna_prior_weight`` is the one PER-TRANSCRIPT array here and it is passed **flat**, while every
    other prior array is subscripted by ``ids`` into per-locus order. That asymmetry is the point: the
    RNA prior's TOTAL is a property of the locus, its ALLOCATION is a property of the transcripts, and
    the C++ remaps the latter itself when it builds each sub-problem. ⛔ Subscripting it by ``ids``
    would silently read transcript rows at locus positions.
    """
    if not multi_loci:
        return

    parts = [partitions[loc.multi_locus_id] for loc in multi_loci]
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
    locus_t_lists = [loc.transcript_indices for loc in multi_loci]
    ids = [loc.multi_locus_id for loc in multi_loci]
    gdna_prior = np.ascontiguousarray(gdna_prior_count[ids], dtype=np.float64)
    rna_prior = np.ascontiguousarray(rna_prior_count[ids], dtype=np.float64)
    g_eff = np.ascontiguousarray(gdna_eff_len[ids], dtype=np.float64)
    enable_gdna = np.array(
        [
            bool(p.is_spliced.size)
            and bool(np.any((p.is_spliced == 0) & np.isfinite(p.gdna_log_liks)))
            for p in parts
        ],
        dtype=np.uint8,
    )

    em_result = estimator.run_batch_locus_em_partitioned(
        partition_tuples,
        locus_t_lists,
        gdna_prior,
        index,
        rna_prior_count=rna_prior,
        # ⛔ FLAT and per-TRANSCRIPT — deliberately not `[ids]`. ⚠ This boundary was dropped once while
        # deleting an unrelated caller, leaving the parameter accepted and silently ignored: every
        # allocation, however extreme, produced byte-identical output. A FIRE COUNTER CANNOT SEE THAT —
        # it counts nonzero weights, not weights the solver read. What caught it was an arm that
        # injected a maximally WRONG allocation and still matched `base` exactly
        # (`quant_accuracy.py --arm oracle_alloc_flip`, TRAPS: could-the-arm-have-fired).
        rna_prior_weight=rna_prior_weight,
        gdna_eff_len=g_eff,
        enable_gdna=enable_gdna,
        em_iterations=em_config.iterations,
        em_convergence_delta=em_config.convergence_delta,
        emit_locus_stats=emit_locus_stats,
        emit_assignments=annotations is not None,
    )
    total_gdna_em, rna_arr, gdna_arr = em_result[0], em_result[1], em_result[2]
    if annotations is not None and len(em_result) > 3:
        _populate_em_annotations(
            parts, em_result[3], em_result[4], em_result[5], annotations, index
        )

    t_to_g = index.t_to_g_arr
    if "is_synthetic" in index.g_df.columns:
        is_synth_g = index.g_df["is_synthetic"].to_numpy()
    else:
        is_synth_g = np.zeros(len(index.g_df), dtype=bool)

    for i, loc in enumerate(multi_loci):
        lid = loc.multi_locus_id
        gene_set = {
            int(t_to_g[int(t)])
            for t in loc.transcript_indices
            if not is_synth_g[int(t_to_g[int(t)])]
        }
        estimator.locus_results.append(
            {
                "locus_id": lid,
                "locus_span_bp": loc.gdna_span,
                "n_transcripts": len(loc.transcript_indices),
                "n_genes": len(gene_set),
                "n_em_fragments": len(loc.unit_indices),
                "rna_total": float(rna_arr[i]),
                "gdna": float(gdna_arr[i]),
                "gdna_prior_count": float(gdna_prior_count[lid]),
                "rna_prior_count": float(rna_prior_count[lid]),
                "enable_gdna": int(enable_gdna[i]),
                "gdna_eff_len_em": float(gdna_eff_len[lid]),
            }
        )

    estimator._gdna_em_total = total_gdna_em
    n_units = sum(len(loc.unit_indices) for loc in multi_loci)
    logger.info(
        "[DONE] Per-locus EM: %d loci, %d ambiguous fragments, gDNA EM=%.0f",
        len(multi_loci),
        n_units,
        total_gdna_em,
    )


def quant_from_buffer(
    buffer: FragmentBuffer,
    index: TranscriptIndex,
    strand_models: StrandModels,
    fl_models: "FLModels",
    region_arrays: "RegionArrays",
    stats: PipelineStats,
    calibration: "CalibrationResult",
    calibration_payload: "object",
    *,
    em_config: EMConfig | None = None,
    scoring: FragmentScoringConfig | None = None,
    log_every: int = 1_000_000,
    annotations: "AnnotationTable | None" = None,
    emit_locus_stats: bool = False,
) -> tuple[AbundanceEstimator, "CalibrationResult"]:
    """Quantify buffered fragments: calibration prior → per-locus EM (PR 6).

    Scores the buffer (RNA/gDNA FL models built from the calibrated pmfs),
    builds connected-component loci, turns the per-region ``CalibrationResult``
    into the per-locus gDNA/RNA split prior (``assemble_priors``), partitions
    the global CSR, and runs the per-locus EM. Returns the populated
    ``AbundanceEstimator`` and the (unchanged) ``CalibrationResult``.
    """
    from .calibration.priors import assemble_priors
    from .locus import build_multi_loci
    from .locus_partition import partition_and_free

    em_config = em_config or EMConfig()
    scoring_cfg = scoring or FragmentScoringConfig()

    # Scorer FL models from the calibrated pmfs (PR 4c FLModels → scoring LUTs).
    rna_fl = FragmentLengthModel.from_pmf(fl_models.rna_pmf, fl_models.max_size)
    gdna_fl = FragmentLengthModel.from_pmf(fl_models.gdna_pmf, fl_models.max_size)

    geometry, estimator = _setup_geometry_and_estimator(
        index, rna_fl, em_config, calibration=calibration, region_arrays=region_arrays
    )

    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        rna_fl,
        gdna_fl,
        stats,
        estimator,
        scoring_cfg,
        log_every,
        annotations,
    )

    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)

    if getattr(em_data, "n_units", 0) == 0 or not multi_loci:
        return estimator, calibration

    priors = assemble_priors(calibration, region_arrays, multi_loci)
    partitions = partition_and_free(em_data, multi_loci)
    _run_locus_em_partitioned(
        estimator,
        partitions,
        multi_loci,
        index,
        priors.gdna_prior_count,
        priors.rna_prior_count,
        priors.gdna_eff_len,
        em_config=em_config,
        annotations=annotations,
        emit_locus_stats=emit_locus_stats,
    )
    del partitions
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
    stats, strand_models, buffer, calibration_payload = scan_and_buffer(bam_path, index, scan)

    # -- The second pass: drain the side buffer, BEFORE calibration --
    if calibration_payload is not None:
        calibration_payload = _drain_side_buffer(
            calibration_payload, index, strand_models, seed=config.second_pass_seed
        )

    # -- Calibration (acyclic) --
    # Build the region geometry, verify it boundaries up 1:1 with the accumulator
    # payload, then hand both (plus the trained strand model and the gDNA FL pmf)
    # to the calibrator. Single feed-forward pass: deconvolve each region into
    # gDNA/RNA and derive ρ_0 + per-region exposure. See
    from .calibration import calibrate
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import (
        build_boundary_flags_array,
        build_contiguous_boundary_reach_arrays,
        build_mature_wall_distances,
        build_sj_geometry_arrays,
    )

    _warn_if_calibration_strand_unidentifiable(strand_models)
    strand_ci_eps = strand_models.strand_specificity_ci_epsilon(confidence=0.99)
    logger.info(
        "[CAL] Strand trainer: n_spliced_obs=%d  ss_est=%.6f  ε_CI(99%%)=%.4g",
        strand_models.n_observations,
        strand_models.strand_specificity,
        strand_ci_eps,
    )

    # ⚠ No alignment check here: CalibrationSubstrate.from_payload runs the identical one inside
    # calibrate(), microseconds later, and two copies of an invariant is one too many.
    region_arrays = RegionArrays.from_index(index)
    boundary_flags = build_boundary_flags_array(index)
    # ⭐ The two annotation-only WALL inputs, beside the other index-derived arrays. They are consulted
    # only when `CalibrationConfig.background_abundance == "measured_total"`, and building them
    # unconditionally keeps production and `scan_cache.index_derived_inputs` on ONE code path — the
    # alternative is a flag-shaped branch here that the instruments do not take.
    mature_walls = build_mature_wall_distances(index, region_arrays)
    boundary_reach = build_contiguous_boundary_reach_arrays(index)
    # The sj axis, in the accumulator's own sj slot order: where each sj attaches,
    # its TRANSCRIPT strand, and its exonic reach either side. The calibrator places it as a FACTOR on
    # its two endpoint regions — never as a message channel, since every sj closes an undirected
    # loop and the graph is not a polytree.
    sj = build_sj_geometry_arrays(index)

    # The two COMPONENT fragment-length models the calibrator's effective lengths need, each fitted
    # from a pool that is PURE BY CONSTRUCTION: gDNA from fragments
    # contained in an intergenic or intronic region, RNA from fragments that used an annotated sj
    # with the splice OBSERVED. Both are smooth-EB shrunk toward the unconditional global FL.
    #
    # ⭐ ALL THREE come from the PAYLOAD — one object, one frame, one definition of length. The
    # scanner's spliced histogram is transcript-space and needs a UNIQUE transcript; the
    # accumulator's pool is a structural rule over a larger population and is binned at the same L as
    # everything else. A fragment enters a pool when exactly ONE hypothesis survived, so its `L` is
    # not in doubt however it was arrived at — determinacy, not provenance.
    #
    # ⭐ TRAPS: pure-and-length-censored.1: the ANCHOR moved here too, off the scanner's histogram and onto the accumulator's own
    # `deposited_lengths`. Until then the pools were accumulator-frame and the anchor they were
    # shrunk toward was not.
    #
    # ⭐ And the RNA pool is de-tilted by its own sj opportunity: "used an annotated sj"
    # is a length-dependent selection, so the raw pool is measurably longer than the library.
    from .calibration.fl import build_fl_models
    from .calibration.gdna_opportunity import gdna_opportunity_from_index
    from .calibration.sj_opportunity import crossing_probability_from_index

    fl_models = build_fl_models(
        calibration_payload,
        sj_opportunity=crossing_probability_from_index(
            index, int(calibration_payload.max_length)
        ),
        gdna_opportunity=gdna_opportunity_from_index(index, int(calibration_payload.max_length)),
    )
    gdna_fl_pmf = fl_models.gdna_pmf
    _bins = np.arange(gdna_fl_pmf.size, dtype=np.float64)
    logger.info(
        "[CAL] FL models: n_gdna_pool=%.0f gdna_fl_mean=%.1f rna_fl_mean=%.1f global_fl_mean=%.1f",
        fl_models.n_gdna,
        float(np.dot(_bins, fl_models.gdna_pmf)),
        float(np.dot(_bins, fl_models.rna_pmf)),
        float(np.dot(_bins, fl_models.global_pmf)),
    )

    # NOTE: the buffer is NOT freed here — quant_from_buffer scans it below.
    _calib_diag: dict = {}
    calibration = calibrate(
        payload=calibration_payload,
        region_arrays=region_arrays,
        strand_model=strand_models,
        gdna_fl_pmf=gdna_fl_pmf,
        rna_fl_pmf=fl_models.rna_pmf,
        config=config.calibration,
        sj=sj,
        diagnostics_out=_calib_diag,
        boundary_flags=boundary_flags,
        mature_walls=mature_walls,
        boundary_reach=boundary_reach,
    )
    calibration_diagnostics = _calib_diag.get("calibration")

    # ⭐ FRAGMENTS, not object incidences. The sum still runs over all three axes — gDNA is contained in
    # a region or crosses a boundary, RNA also jumps, and at a donor boundary the sj flux IS the gene's whole
    # mature output — but each axis is converted by its own population's mass-per-crossing first. Adding
    # the raw banks counted one fragment once per boundary it crossed AND once per sj it used, which
    # read `f_gdna` 0.3851 against a truth of 0.5085 on ladder g50 capture_off.
    logger.info(
        "calibration: N=%d E=%d J=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gDNA_fragments=%.0f RNA_fragments=%.0f (sj incidences %.0f)",
        calibration.n_regions,
        calibration.n_boundaries,
        calibration.n_sj,
        calibration.gdna_density_global,
        calibration.rna_sense_frac,
        calibration.library_gdna_fragments,
        calibration.library_rna_fragments,
        float(calibration.count_rna_sj.sum()),
    )

    # -- Quantification: calibration prior → per-locus EM --
    # Annotation table for the optional second BAM pass (opt-in).
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
            fl_models,
            region_arrays,
            stats,
            calibration,
            calibration_payload,
            em_config=config.em,
            scoring=config.scoring,
            log_every=scan.log_every,
            annotations=annotations,
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

    # Genome-wide gDNA track for the QC report (+ bedGraph). Pure persistence of
    # the per-region calibration solution; skipped if calibration did not run.
    calibration_track = None
    if calibration is not None:
        from .calibration.track import build_gdna_track

        calibration_track = build_gdna_track(calibration, region_arrays, index.ref_names)

    return PipelineResult(
        stats=stats,
        strand_models=strand_models,
        estimator=estimator,
        pipeline_config=config,
        calibration=calibration,
        fl_models=fl_models,
        calibration_track=calibration_track,
        calibration_diagnostics=calibration_diagnostics,
    )
