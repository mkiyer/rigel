#!/usr/bin/env python3
"""Profile the partition-native pipeline on CAPAN-1.

Exercises the actual quant_from_buffer() path (which now uses
partition_and_free + batch_locus_em_partitioned) with per-phase
RSS tracking and timing.
"""
import gc
import json
import logging
import resource
import sys
import time
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# RSS helper
# ---------------------------------------------------------------------------

def rss_mb():
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def peak_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

BAM = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_runs/human/ccle_capan_1_pancreas/rigel/annotated.bam"
INDEX = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
OUTDIR = Path("/home/mkiyer/proj/rigel/results/profile_partition_v1")

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = {}
    rss_snapshots = {}

    rss_snapshots["start"] = rss_mb()
    logger.info(f"RSS start: {rss_snapshots['start']:.0f} MB")

    # ── Load index ──────────────────────────────────────────
    from rigel.index import TranscriptIndex

    t0 = time.perf_counter()
    index = TranscriptIndex.load(INDEX)
    t_index = time.perf_counter() - t0
    results["index_load_sec"] = t_index
    rss_snapshots["after_index"] = rss_mb()
    logger.info(f"Index loaded: {index.num_transcripts} transcripts, "
                f"{index.num_genes} genes ({t_index:.1f}s, "
                f"RSS={rss_snapshots['after_index']:.0f} MB)")

    # ── Scan + buffer ───────────────────────────────────────
    from rigel.config import BamScanConfig, CalibrationConfig
    from rigel.native import detect_sj_strand_tag
    from rigel.pipeline import scan_and_buffer

    sj_tag = detect_sj_strand_tag(BAM)
    scan_cfg = BamScanConfig(sj_strand_tag=sj_tag, include_multimap=True)

    t0 = time.perf_counter()
    stats, strand_models, frag_length_models, buffer, region_counts, fl_table = (
        scan_and_buffer(BAM, index, scan_cfg)
    )
    t_scan = time.perf_counter() - t0
    results["scan_and_buffer_sec"] = t_scan
    rss_snapshots["after_scan"] = rss_mb()
    logger.info(f"Scan: {stats.total:,} alignments, {buffer.total_fragments:,} "
                f"buffered ({t_scan:.1f}s, RSS={rss_snapshots['after_scan']:.0f} MB)")

    # ── Finalize models ─────────────────────────────────────
    t0 = time.perf_counter()
    strand_models.finalize()
    frag_length_models.build_scoring_models()
    frag_length_models.finalize()
    t_finalize = time.perf_counter() - t0
    results["finalize_models_sec"] = t_finalize

    # ── Calibration ─────────────────────────────────────────
    from rigel.calibration import calibrate_gdna

    cal_cfg = CalibrationConfig()
    t0 = time.perf_counter()
    calibration = calibrate_gdna(
        region_counts, fl_table, index.region_df,
        strand_models.strand_specificity,
        max_iterations=cal_cfg.max_iterations,
        convergence_tol=cal_cfg.convergence_tol,
        density_percentile=cal_cfg.density_percentile,
        min_gdna_regions=cal_cfg.min_gdna_regions,
        min_fl_ess=cal_cfg.min_fl_ess,
        intergenic_fl_model=frag_length_models.intergenic,
    )
    frag_length_models.gdna_model = calibration.gdna_fl_model
    t_cal = time.perf_counter() - t0
    results["calibration_sec"] = t_cal
    rss_snapshots["after_calibration"] = rss_mb()
    logger.info(f"Calibration: {t_cal:.1f}s, RSS={rss_snapshots['after_calibration']:.0f} MB")

    # ── quant_from_buffer (NEW partitioned path) ────────────
    from rigel.config import EMConfig, FragmentScoringConfig, TranscriptGeometry
    from rigel.estimator import AbundanceEstimator
    from rigel.locus import build_loci, compute_locus_priors
    from rigel.partition import partition_and_free
    from rigel.scan import FragmentRouter
    from rigel.scoring import FragmentScorer

    em_config = EMConfig()
    scoring = FragmentScoringConfig()

    # 3a: Geometry
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    effective_lengths = frag_length_models.rna_model.compute_all_transcript_eff_lens(
        exonic_lengths.astype(np.int64))
    transcript_spans = (index.t_df["end"].values - index.t_df["start"].values).astype(np.float64)
    geometry = TranscriptGeometry(
        effective_lengths=effective_lengths,
        exonic_lengths=exonic_lengths,
        t_to_g=index.t_to_g_arr,
        transcript_spans=transcript_spans,
    )
    estimator = AbundanceEstimator(
        index.num_transcripts, em_config=em_config, geometry=geometry,
        is_nrna=index.t_df["is_nrna"].values,
        is_synthetic=index.t_df["is_synthetic"].values,
    )
    ctx = FragmentScorer.from_models(
        strand_models, frag_length_models, index, estimator,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
    )

    # 3d: Fragment router scan (scoring)
    t0 = time.perf_counter()
    builder = FragmentRouter(ctx, estimator, stats, index, strand_models)
    em_data = builder.scan(buffer, log_every=1_000_000)
    t_route = time.perf_counter() - t0
    results["fragment_router_scan_sec"] = t_route
    rss_snapshots["after_router_scan"] = rss_mb()
    logger.info(f"Router scan: {em_data.n_units:,} units, {em_data.n_candidates:,} "
                f"candidates ({t_route:.1f}s, RSS={rss_snapshots['after_router_scan']:.0f} MB)")

    # Free scanner accumulators + buffer
    del builder, ctx
    buffer.release()
    gc.collect()
    rss_snapshots["after_buffer_release"] = rss_mb()
    logger.info(f"After buffer release: RSS={rss_snapshots['after_buffer_release']:.0f} MB")

    results["em_n_units"] = em_data.n_units
    results["em_n_candidates"] = em_data.n_candidates

    # 3e: Build loci
    t0 = time.perf_counter()
    loci = build_loci(em_data, index)
    t_loci = time.perf_counter() - t0
    results["build_loci_sec"] = t_loci
    results["n_loci"] = len(loci)
    max_t = max(len(l.transcript_indices) for l in loci)
    max_u = max(len(l.unit_indices) for l in loci)
    results["max_locus_transcripts"] = max_t
    results["max_locus_units"] = max_u
    logger.info(f"Build loci: {len(loci)} loci ({t_loci:.1f}s), "
                f"largest: {max_t} transcripts, {max_u:,} units")

    # 3f: Compute priors
    for locus in loci:
        for t_idx in locus.transcript_indices:
            estimator.locus_id_per_transcript[int(t_idx)] = locus.locus_id

    t0 = time.perf_counter()
    alpha_gdna, alpha_rna = compute_locus_priors(loci, index, calibration=calibration)
    t_priors = time.perf_counter() - t0
    results["compute_priors_sec"] = t_priors
    rss_snapshots["after_priors"] = rss_mb()
    logger.info(f"Priors computed: {t_priors:.1f}s")

    # ── NEW: Partition and free ─────────────────────────────
    t0 = time.perf_counter()
    partitions = partition_and_free(em_data, loci)
    t_partition = time.perf_counter() - t0
    results["partition_and_free_sec"] = t_partition
    rss_snapshots["after_partition"] = rss_mb()
    logger.info(f"Partition + free: {t_partition:.1f}s, "
                f"RSS={rss_snapshots['after_partition']:.0f} MB")

    del em_data
    gc.collect()
    rss_snapshots["after_em_data_del"] = rss_mb()
    logger.info(f"After em_data del: RSS={rss_snapshots['after_em_data_del']:.0f} MB")

    # ── NEW: Partitioned EM (mega + normal) ─────────────────
    import os
    n_threads = em_config.n_threads or os.cpu_count() or 1

    # Classify mega vs normal
    total_work = sum(
        len(l.transcript_indices) * partitions[l.locus_id].n_units for l in loci
    )
    fair_share = total_work // n_threads if n_threads > 1 else total_work + 1

    mega_loci = sorted(
        [l for l in loci
         if len(l.transcript_indices) * partitions[l.locus_id].n_units >= fair_share],
        key=lambda l: len(l.transcript_indices) * partitions[l.locus_id].n_units,
        reverse=True,
    )
    mega_ids = {l.locus_id for l in mega_loci}
    normal_loci = [l for l in loci if l.locus_id not in mega_ids]
    results["n_mega_loci"] = len(mega_loci)
    results["n_normal_loci"] = len(normal_loci)

    def _pack_and_call(parts, batch_loci, batch_ag, batch_ar, batch_spans, emit_stats=False):
        partition_tuples = [
            (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
             p.tx_starts, p.tx_ends, p.count_cols,
             p.is_spliced, p.gdna_log_liks, p.genomic_footprints,
             p.locus_t_indices, p.locus_count_cols)
            for p in parts
        ]
        locus_t_lists = [l.transcript_indices for l in batch_loci]
        return estimator.run_batch_locus_em_partitioned(
            partition_tuples, locus_t_lists,
            batch_ag, batch_ar, batch_spans, index,
            em_iterations=em_config.iterations,
            em_convergence_delta=em_config.convergence_delta,
            emit_locus_stats=emit_stats,
        )

    # Phase A: Mega-loci
    t0 = time.perf_counter()
    total_gdna_em = 0.0
    for i, locus in enumerate(mega_loci):
        t_mega_start = time.perf_counter()
        part = partitions.pop(locus.locus_id)
        gdna_em, mrna_arr, gdna_arr = _pack_and_call(
            [part], [locus],
            np.array([alpha_gdna[locus.locus_id]], dtype=np.float64),
            np.array([alpha_rna[locus.locus_id]], dtype=np.float64),
            np.array([locus.gdna_span], dtype=np.int64),
            emit_stats=True,
        )
        total_gdna_em += gdna_em
        del part
        gc.collect()
        t_mega_elapsed = time.perf_counter() - t_mega_start
        rss_now = rss_mb()
        logger.info(f"  Mega-locus {i} (id={locus.locus_id}): "
                    f"{len(locus.transcript_indices)} tx, "
                    f"{len(locus.unit_indices):,} units, "
                    f"{t_mega_elapsed:.1f}s, RSS={rss_now:.0f} MB")
    t_mega = time.perf_counter() - t0
    results["mega_em_sec"] = t_mega
    rss_snapshots["after_mega_em"] = rss_mb()
    logger.info(f"Mega EM: {len(mega_loci)} loci ({t_mega:.1f}s, "
                f"RSS={rss_snapshots['after_mega_em']:.0f} MB)")

    # Phase B: Normal loci
    t0 = time.perf_counter()
    if normal_loci:
        normal_parts = [partitions[l.locus_id] for l in normal_loci]
        normal_ag = np.array(
            [alpha_gdna[l.locus_id] for l in normal_loci], dtype=np.float64)
        normal_ar = np.array(
            [alpha_rna[l.locus_id] for l in normal_loci], dtype=np.float64)
        normal_spans = np.array(
            [l.gdna_span for l in normal_loci], dtype=np.int64)
        gdna_em, mrna_arr, gdna_arr = _pack_and_call(
            normal_parts, normal_loci, normal_ag, normal_ar, normal_spans,
            emit_stats=True,
        )
        total_gdna_em += gdna_em
    t_normal = time.perf_counter() - t0
    results["normal_em_sec"] = t_normal
    rss_snapshots["after_normal_em"] = rss_mb()
    logger.info(f"Normal EM: {len(normal_loci)} loci ({t_normal:.1f}s, "
                f"RSS={rss_snapshots['after_normal_em']:.0f} MB)")

    del partitions
    gc.collect()
    rss_snapshots["after_em_cleanup"] = rss_mb()
    logger.info(f"After EM cleanup: RSS={rss_snapshots['after_em_cleanup']:.0f} MB")

    results["total_locus_em_sec"] = t_mega + t_normal
    results["total_gdna_em"] = total_gdna_em

    # ── Summary ─────────────────────────────────────────────
    total_wall = (results["scan_and_buffer_sec"] + results["finalize_models_sec"]
                  + results["calibration_sec"] + results["fragment_router_scan_sec"]
                  + results["build_loci_sec"] + results["compute_priors_sec"]
                  + results["partition_and_free_sec"] + results["total_locus_em_sec"])
    results["total_wall_sec"] = total_wall
    results["peak_rss_mb"] = peak_rss_mb()
    results["rss_snapshots"] = rss_snapshots

    # Save locus stats
    locus_stats = getattr(estimator, "locus_stats", None)
    if locus_stats:
        results["n_locus_stats"] = len(locus_stats)
        # Summarize top-10 loci by time
        top10 = sorted(locus_stats, key=lambda d: d.get("total_us", 0), reverse=True)[:10]
        results["top10_loci"] = top10
        # Save full locus stats
        with open(OUTDIR / "locus_stats.json", "w") as f:
            json.dump(locus_stats, f, indent=2, default=str)
        logger.info(f"Saved {len(locus_stats)} locus stats")

        # Convergence analysis
        max_iters = em_config.iterations // 3  # SQUAREM budget divisor
        at_max = [s for s in locus_stats if s.get("squarem_iterations", 0) >= max_iters]
        results["loci_at_max_iterations"] = len(at_max)
        results["max_iterations_threshold"] = max_iters
        if at_max:
            at_max_time = sum(s.get("total_us", 0) for s in at_max) / 1e6
            results["time_at_max_iters_sec"] = at_max_time

    # Print summary
    logger.info("")
    logger.info("=" * 60)
    logger.info("PROFILE SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total wall time:        {total_wall:.1f}s")
    logger.info(f"Peak RSS:               {results['peak_rss_mb']:.0f} MB")
    logger.info(f"")
    logger.info(f"scan_and_buffer:        {results['scan_and_buffer_sec']:.1f}s")
    logger.info(f"finalize_models:        {results['finalize_models_sec']:.1f}s")
    logger.info(f"calibration:            {results['calibration_sec']:.1f}s")
    logger.info(f"fragment_router_scan:   {results['fragment_router_scan_sec']:.1f}s")
    logger.info(f"build_loci:             {results['build_loci_sec']:.1f}s")
    logger.info(f"compute_priors:         {results['compute_priors_sec']:.1f}s")
    logger.info(f"partition_and_free:     {results['partition_and_free_sec']:.1f}s")
    logger.info(f"mega_em ({results['n_mega_loci']} loci):        {results['mega_em_sec']:.1f}s")
    logger.info(f"normal_em ({results['n_normal_loci']} loci):    {results['normal_em_sec']:.1f}s")
    logger.info(f"total_locus_em:         {results['total_locus_em_sec']:.1f}s")
    logger.info(f"")
    logger.info(f"EM units:               {results['em_n_units']:,}")
    logger.info(f"EM candidates:          {results['em_n_candidates']:,}")
    logger.info(f"Loci:                   {results['n_loci']}")
    logger.info(f"Mega-loci:              {results['n_mega_loci']}")
    if "loci_at_max_iterations" in results:
        logger.info(f"Loci at max iters:      {results['loci_at_max_iterations']}")
        if "time_at_max_iters_sec" in results:
            logger.info(f"Time at max iters:      {results['time_at_max_iters_sec']:.1f}s")
    logger.info(f"")
    logger.info("RSS snapshots:")
    for k, v in rss_snapshots.items():
        logger.info(f"  {k:30s} {v:,.0f} MB")

    # Save results
    with open(OUTDIR / "profile_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info(f"\nResults saved to {OUTDIR}")


if __name__ == "__main__":
    main()
