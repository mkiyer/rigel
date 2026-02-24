#!/usr/bin/env python
"""Profile hulkrna pipeline stages and hot functions.

Usage:
    PYTHONPATH=src python scripts/profile_pipeline.py \
        --bam /path/to/reads_namesort.bam \
        --index /path/to/hulkrna_index/ \
        [--cprofile]  # Enable cProfile for hot-function analysis
"""
import argparse
import cProfile
import io
import logging
import math
import pstats
import sys
import time
from pathlib import Path

import numpy as np

# Ensure hulkrna is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.bam import detect_sj_strand_tag, parse_bam_file
from hulkrna.buffer import (
    FragmentBuffer,
    FRAG_UNIQUE,
    FRAG_ISOFORM_AMBIG,
    FRAG_MULTIMAPPER,
    FRAG_CHIMERIC,
)
from hulkrna.fragment import Fragment
from hulkrna.index import HulkIndex
from hulkrna.pipeline import (
    scan_and_buffer,
    count_from_buffer,
    _scan_and_build_em_data,
    _build_loci,
    _build_locus_em_data,
    _compute_eb_gdna_priors,
    _compute_nrna_init,
    _GDNA_SPLICE_PENALTIES,
    _SPLICE_UNANNOT,
    _DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    overhang_alpha_to_log_penalty,
    DEFAULT_OVERHANG_ALPHA,
    DEFAULT_MISMATCH_ALPHA,
)
from hulkrna.estimator import AbundanceEstimator
from hulkrna.frag_length_model import FragmentLengthModels
from hulkrna.resolution import pair_multimapper_reads, resolve_fragment
from hulkrna.stats import PipelineStats
from hulkrna.strand_model import StrandModels

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


class Timer:
    """Simple context-manager timer."""
    def __init__(self, label: str):
        self.label = label
        self.elapsed = 0.0
    def __enter__(self):
        self._start = time.perf_counter()
        return self
    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self._start


def profile_stages(bam_path: str, index_dir: str, enable_cprofile: bool = False):
    """Profile each pipeline stage independently."""
    import pysam

    print("=" * 72)
    print("HULKRNA PIPELINE PROFILER")
    print("=" * 72)
    print(f"BAM:   {bam_path}")
    print(f"Index: {index_dir}")
    print()

    # ---------------------------------------------------------------
    # Stage 0: Load index
    # ---------------------------------------------------------------
    with Timer("Load index") as t0:
        index = HulkIndex.load(index_dir)
    print(f"[Stage 0] Load index: {t0.elapsed:.3f}s")
    print(f"  Transcripts: {index.num_transcripts:,}")
    print(f"  Genes:       {index.num_genes:,}")
    print()

    # ---------------------------------------------------------------
    # Stage 1: Detect SJ strand tag
    # ---------------------------------------------------------------
    with Timer("Detect SJ strand tag") as t1:
        detected = detect_sj_strand_tag(bam_path)
        if len(detected) == 1:
            sj_strand_tag = detected[0]
        elif len(detected) > 1:
            sj_strand_tag = detected
        else:
            sj_strand_tag = ()
    print(f"[Stage 1] Detect SJ strand tag: {t1.elapsed:.3f}s")
    print(f"  Tag(s): {sj_strand_tag}")
    print()

    # ---------------------------------------------------------------
    # Stage 2: scan_and_buffer (the big BAM pass)
    # ---------------------------------------------------------------
    bamfh = pysam.AlignmentFile(bam_path, "rb")
    try:
        with Timer("scan_and_buffer") as t2:
            stats, strand_models, frag_length_models, buffer = scan_and_buffer(
                bamfh.fetch(until_eof=True),
                index,
                skip_duplicates=True,
                include_multimap=False,
                log_every=10_000,
                max_frag_length=1000,
                sj_strand_tag=sj_strand_tag,
            )
    finally:
        bamfh.close()

    print(f"[Stage 2] scan_and_buffer: {t2.elapsed:.3f}s")
    print(f"  Fragments:       {stats.n_fragments:,}")
    print(f"  Buffered:        {buffer.total_fragments:,}")
    print(f"  Unique gene:     {stats.n_unique_gene:,}")
    print(f"  Multi gene:      {stats.n_multi_gene:,}")
    print(f"  Intergenic:      {stats.n_intergenic:,}")
    print(f"  Buffer memory:   {buffer.memory_bytes / 1024**2:.1f} MB")
    print(f"  Spilled chunks:  {buffer.n_spilled}")
    print()

    # ---------------------------------------------------------------
    # Stage 2b: Sub-profile scan_and_buffer internals
    # ---------------------------------------------------------------
    # Re-run the BAM pass with fine-grained timing to profile
    # parse_bam, Fragment.from_reads, resolve_fragment vs rest
    print("--- Sub-profiling scan_and_buffer internals ---")

    bamfh = pysam.AlignmentFile(bam_path, "rb")
    bam_stats2 = PipelineStats().as_bam_stats_dict()

    t_parse = 0.0
    t_fragment_build = 0.0
    t_resolve = 0.0
    t_model_train = 0.0
    t_buffer_append = 0.0
    t_mm_pair = 0.0
    n_frags = 0
    n_resolved = 0
    n_intergenic = 0

    strand_models2 = StrandModels()
    frag_length_models2 = FragmentLengthModels(max_size=1000)
    buffer2 = FragmentBuffer(
        t_to_g_arr=index.t_to_g_arr,
        chunk_size=1_000_000,
        max_memory_bytes=2 * 1024**3,
    )

    _t0 = time.perf_counter()
    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        bamfh.fetch(until_eof=True),
        bam_stats2,
        skip_duplicates=True,
        include_multimap=False,
    ):
        _t1 = time.perf_counter()
        t_parse += _t1 - _t0

        # Multimapper pairing
        secondary_pairs = []
        if sec_r1 or sec_r2:
            _ta = time.perf_counter()
            secondary_pairs = pair_multimapper_reads(
                sec_r1, sec_r2, index,
                sj_strand_tag=sj_strand_tag,
            )
            t_mm_pair += time.perf_counter() - _ta

        all_hits = list(hits) + secondary_pairs
        num_hits = max(nh, len(all_hits))
        is_unique_mapper = num_hits == 1

        for r1_reads, r2_reads in all_hits:
            _ta = time.perf_counter()
            frag = Fragment.from_reads(r1_reads, r2_reads, sj_strand_tag=sj_strand_tag)
            _tb = time.perf_counter()
            t_fragment_build += _tb - _ta
            n_frags += 1

            if not frag.exons:
                continue

            _ta = time.perf_counter()
            result = resolve_fragment(frag, index)
            _tb = time.perf_counter()
            t_resolve += _tb - _ta

            if result is None:
                n_intergenic += 1
                continue

            n_resolved += 1
            result.num_hits = num_hits

            # Model training (simplified timing)
            _ta = time.perf_counter()
            # ... (just time the overhead, don't re-train)
            t_model_train += time.perf_counter() - _ta

            _ta = time.perf_counter()
            buffer2.append(result, frag_id=n_frags)
            t_buffer_append += time.perf_counter() - _ta

        _t0 = time.perf_counter()

    bamfh.close()
    buffer2.finalize()

    t_scan_total = t_parse + t_fragment_build + t_resolve + t_mm_pair + t_buffer_append
    print(f"  parse_bam_file:        {t_parse:.3f}s ({t_parse/t_scan_total*100:.1f}%)")
    print(f"  Fragment.from_reads:   {t_fragment_build:.3f}s ({t_fragment_build/t_scan_total*100:.1f}%)")
    print(f"  resolve_fragment:      {t_resolve:.3f}s ({t_resolve/t_scan_total*100:.1f}%)")
    print(f"  MM pairing:            {t_mm_pair:.3f}s ({t_mm_pair/t_scan_total*100:.1f}%)")
    print(f"  buffer.append:         {t_buffer_append:.3f}s ({t_buffer_append/t_scan_total*100:.1f}%)")
    print(f"  Total accounted:       {t_scan_total:.3f}s")
    print(f"  Fragments: {n_frags:,}, Resolved: {n_resolved:,}, Intergenic: {n_intergenic:,}")
    print()
    buffer2.cleanup()

    # ---------------------------------------------------------------
    # Stage 3: Finalize models
    # ---------------------------------------------------------------
    with Timer("Finalize models") as t3:
        strand_models.finalize()
        frag_length_models.finalize()
    print(f"[Stage 3] Finalize models: {t3.elapsed:.3f}s")
    print()

    # ---------------------------------------------------------------
    # Stage 4: count_from_buffer (EM + assignment)
    # ---------------------------------------------------------------
    # We need to break this down further. Replicate the internals:

    # 4a: Compute geometry
    with Timer("Compute geometry") as t4a:
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
    print(f"[Stage 4a] Compute geometry: {t4a.elapsed:.3f}s")

    # 4b: Create AbundanceEstimator
    with Timer("Create estimator") as t4b:
        gdna_splice_penalties = dict(_GDNA_SPLICE_PENALTIES)
        gdna_splice_penalties[_SPLICE_UNANNOT] = _DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
        overhang_log_penalty = overhang_alpha_to_log_penalty(DEFAULT_OVERHANG_ALPHA)
        mismatch_log_penalty = overhang_alpha_to_log_penalty(DEFAULT_MISMATCH_ALPHA)

        counter = AbundanceEstimator(
            index.num_transcripts, index.num_genes,
            seed=None, alpha=0.5,
            effective_lengths=effective_lengths,
            exonic_lengths=exonic_lengths,
            t_to_g=index.t_to_g_arr,
            gene_spans=gene_spans,
            mean_frag=mean_frag,
            intronic_spans=intronic_spans,
            transcript_spans=transcript_spans,
        )
    print(f"[Stage 4b] Create estimator: {t4b.elapsed:.3f}s")

    # 4c: _scan_and_build_em_data
    with Timer("scan_and_build_em_data") as t4c:
        em_data = _scan_and_build_em_data(
            buffer, index, strand_models, frag_length_models,
            counter, stats, log_every=10_000,
            gdna_splice_penalties=gdna_splice_penalties,
            overhang_log_penalty=overhang_log_penalty,
            mismatch_log_penalty=mismatch_log_penalty,
        )
    print(f"[Stage 4c] scan_and_build_em_data: {t4c.elapsed:.3f}s")
    print(f"  EM units:      {em_data.n_units:,}")
    print(f"  EM candidates: {em_data.n_candidates:,}")

    # 4d: Compute nRNA init
    with Timer("Compute nRNA init") as t4d:
        nrna_init = _compute_nrna_init(
            counter.transcript_intronic_sense,
            counter.transcript_intronic_antisense,
            transcript_spans,
            exonic_lengths,
            mean_frag,
            strand_models,
        )
        counter.nrna_init = nrna_init
    print(f"[Stage 4d] Compute nRNA init: {t4d.elapsed:.3f}s")

    # 4e: Build loci
    with Timer("Build loci") as t4e:
        loci = _build_loci(em_data, index) if em_data.n_units > 0 else []
    print(f"[Stage 4e] Build loci: {t4e.elapsed:.3f}s")
    if loci:
        max_t = max(len(l.transcript_indices) for l in loci)
        max_u = max(len(l.unit_indices) for l in loci)
        print(f"  Loci: {len(loci)}, largest: {max_t} transcripts, {max_u} units")

    # 4f: Compute EB gDNA priors
    with Timer("EB gDNA priors") as t4f:
        if loci:
            gdna_inits = _compute_eb_gdna_priors(
                loci, em_data, counter, index, strand_models,
            )
        else:
            gdna_inits = []
    print(f"[Stage 4f] EB gDNA priors: {t4f.elapsed:.3f}s")

    # 4g: Per-locus EM
    with Timer("Per-locus EM") as t4g:
        t_build_em = 0.0
        t_run_em = 0.0
        t_assign = 0.0
        n_loci_processed = 0
        large_loci = []
        for i, locus in enumerate(loci):
            if len(locus.unit_indices) == 0:
                continue
            _ta = time.perf_counter()
            locus_em = _build_locus_em_data(
                locus, em_data, counter, index, mean_frag,
                gdna_init=gdna_inits[i],
            )
            _tb = time.perf_counter()
            t_build_em += _tb - _ta

            _ta = time.perf_counter()
            theta, alpha = counter.run_locus_em(
                locus_em, em_iterations=1000,
            )
            _tb = time.perf_counter()
            t_run_em += _tb - _ta

            _ta = time.perf_counter()
            gdna_assigned = counter.assign_locus_ambiguous(
                locus_em, theta,
                confidence_threshold=0.95,
            )
            _tb = time.perf_counter()
            t_assign += _tb - _ta

            n_loci_processed += 1
            if len(locus.transcript_indices) > 20:
                large_loci.append((
                    i,
                    len(locus.transcript_indices),
                    len(locus.unit_indices),
                    _tb - _ta + (time.perf_counter() - _ta),
                ))

    print(f"[Stage 4g] Per-locus EM: {t4g.elapsed:.3f}s")
    print(f"  Build locus EM data: {t_build_em:.3f}s")
    print(f"  Run locus EM:        {t_run_em:.3f}s")
    print(f"  Assign posteriors:   {t_assign:.3f}s")
    print(f"  Loci processed:      {n_loci_processed:,}")
    print()

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    total = t0.elapsed + t1.elapsed + t2.elapsed + t3.elapsed + t4a.elapsed + t4b.elapsed + t4c.elapsed + t4d.elapsed + t4e.elapsed + t4f.elapsed + t4g.elapsed
    print("=" * 72)
    print("STAGE SUMMARY")
    print("=" * 72)
    stages = [
        ("Load index",              t0.elapsed),
        ("Detect SJ tag",           t1.elapsed),
        ("scan_and_buffer",         t2.elapsed),
        ("Finalize models",         t3.elapsed),
        ("Compute geometry",        t4a.elapsed),
        ("Create estimator",        t4b.elapsed),
        ("scan_and_build_em_data",  t4c.elapsed),
        ("Compute nRNA init",       t4d.elapsed),
        ("Build loci",              t4e.elapsed),
        ("EB gDNA priors",          t4f.elapsed),
        ("Per-locus EM",            t4g.elapsed),
    ]
    for name, dur in stages:
        pct = dur / total * 100 if total > 0 else 0
        bar = "#" * int(pct / 2)
        print(f"  {name:<28s} {dur:7.3f}s  {pct:5.1f}%  {bar}")
    print(f"  {'TOTAL':<28s} {total:7.3f}s")
    print()

    # ---------------------------------------------------------------
    # cProfile the full pipeline if requested
    # ---------------------------------------------------------------
    if enable_cprofile:
        print("=" * 72)
        print("cProfile: Full pipeline run")
        print("=" * 72)
        pr = cProfile.Profile()
        pr.enable()

        # Re-run the full pipeline
        import pysam as _pysam
        bamfh = _pysam.AlignmentFile(bam_path, "rb")
        try:
            stats3, sm3, flm3, buf3 = scan_and_buffer(
                bamfh.fetch(until_eof=True), index,
                skip_duplicates=True, include_multimap=False,
                sj_strand_tag=sj_strand_tag,
            )
        finally:
            bamfh.close()
        sm3.finalize()
        flm3.finalize()

        gdna_sp3 = dict(_GDNA_SPLICE_PENALTIES)
        gdna_sp3[_SPLICE_UNANNOT] = _DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
        counter3 = count_from_buffer(
            buf3, index, sm3, flm3, stats3,
            gdna_splice_penalties=gdna_sp3,
            overhang_log_penalty=overhang_log_penalty,
            mismatch_log_penalty=mismatch_log_penalty,
        )
        buf3.cleanup()

        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s)
        ps.sort_stats("cumulative")
        ps.print_stats(60)
        print(s.getvalue())

        s2 = io.StringIO()
        ps2 = pstats.Stats(pr, stream=s2)
        ps2.sort_stats("tottime")
        ps2.print_stats(60)
        print("\n--- Sorted by total time ---")
        print(s2.getvalue())

    buffer.cleanup()


def main():
    parser = argparse.ArgumentParser(description="Profile hulkrna pipeline")
    parser.add_argument("--bam", required=True, help="Path to name-sorted BAM")
    parser.add_argument("--index", required=True, help="Path to hulkrna index dir")
    parser.add_argument("--cprofile", action="store_true", help="Enable cProfile")
    args = parser.parse_args()
    profile_stages(args.bam, args.index, enable_cprofile=args.cprofile)


if __name__ == "__main__":
    main()
