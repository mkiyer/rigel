#!/usr/bin/env python3
"""Profiling tool for rigel performance analysis.

Runs the rigel pipeline with detailed timing, memory tracking,
and optional cProfile instrumentation.  Outputs machine-readable
JSON alongside a human-readable text report.

Modes
-----
**Simple** (default) — Times ``run_pipeline()`` as a whole with RSS
memory sampling in a background thread.

**Stage decomposition** (``--stages``) — Separately times:
    1. Index load
    2. ``scan_and_buffer`` (C++ BAM scan)
    3. Model finalization
    4. ``quant_from_buffer`` (scoring, routing, locus-level EM)

**cProfile** (``--cprofile``) — Enables Python function-level
profiling for the pipeline run.  Writes a ``.prof`` file that can
be visualized with ``snakeviz`` or ``py-spy``.

**Comparison** — Supply multiple ``rigel_configs`` in the YAML to
profile the same BAM with different parameter sets side by side.

Usage
-----
::

    python scripts/profiler.py --bam reads.bam --index rigel_index/
    python scripts/profiler.py --config scripts/profile_example.yaml
    python scripts/profiler.py --bam reads.bam --index rigel_index/ \\
        --cprofile --stages --outdir profile_results/

Output
------
::

    <outdir>/
        profile_summary.json   — machine-readable results
        profile_report.txt     — human-readable text report
        memory_timeline.csv    — RSS samples over time
        profile_<name>.prof    — cProfile data (if --cprofile)
"""
from __future__ import annotations

import argparse
import cProfile
import csv
import gc
import io
import json
import logging
import os
import platform
import pstats
import resource
import sys
import threading
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]

from rigel.config import (
    BamScanConfig,
    EMConfig,
    FragmentScoringConfig,
    PipelineConfig,
    TranscriptGeometry,
)
from rigel.estimator import AbundanceEstimator
from rigel.index import TranscriptIndex
from rigel.locus import build_loci, build_locus_em_data, compute_eb_gdna_priors, compute_nrna_init
from rigel.estimator import compute_global_gdna_density, compute_hybrid_nrna_frac_priors
from rigel.pipeline import quant_from_buffer, run_pipeline, scan_and_buffer, _compute_intergenic_density
from rigel.scan import FragmentRouter
from rigel.scoring import (
    GDNA_SPLICE_PENALTIES,
    SPLICE_UNANNOT,
    FragmentScorer,
    overhang_alpha_to_log_penalty,
)
from rigel.stats import PipelineStats
from rigel.strand_model import StrandModels

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
# Memory tracking
# ═══════════════════════════════════════════════════════════════════


def _get_rss_mb() -> float:
    """Current peak RSS of this process in MB."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    if sys.platform == "darwin":
        return usage.ru_maxrss / (1024 * 1024)  # bytes on macOS
    return usage.ru_maxrss / 1024  # KB on Linux


class MemoryTimeline:
    """Background thread that samples RSS at regular intervals.

    The timeline records (elapsed_sec, rss_mb) tuples which can be
    written to CSV for visualization.
    """

    def __init__(self, interval_sec: float = 0.1):
        self.interval_sec = interval_sec
        self.samples: list[tuple[float, float]] = []
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._t0 = 0.0

    def start(self) -> None:
        self._t0 = time.monotonic()
        self._stop.clear()
        self.samples = [(_snap_rss_current(), 0.0)]
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)
        self.samples.append((_snap_rss_current(), time.monotonic() - self._t0))

    def _run(self) -> None:
        while not self._stop.is_set():
            elapsed = time.monotonic() - self._t0
            rss = _snap_rss_current()
            self.samples.append((rss, elapsed))
            self._stop.wait(self.interval_sec)

    @property
    def peak_mb(self) -> float:
        return max((s[0] for s in self.samples), default=0.0)

    def write_csv(self, path: Path) -> None:
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["elapsed_sec", "rss_mb"])
            for rss, elapsed in self.samples:
                w.writerow([f"{elapsed:.3f}", f"{rss:.1f}"])
        logger.info("Wrote memory timeline: %s (%d samples)", path, len(self.samples))


# Pre-initialize macOS mach task_info for fast RSS snapshots
_MACH_TASK_INFO_READY = False
_mach_libc = None
_MachTaskBasicInfo = None
_MACH_TASK_BASIC_INFO = 20

if sys.platform == "darwin":
    try:
        import ctypes
        import ctypes.util

        class _MachTaskBasicInfoStruct(ctypes.Structure):
            _fields_ = [
                ("suspend_count", ctypes.c_int32),
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("user_time", ctypes.c_uint64),
                ("system_time", ctypes.c_uint64),
                ("policy", ctypes.c_int32),
            ]

        _mach_libc = ctypes.CDLL(ctypes.util.find_library("c"))
        _mach_libc.mach_task_self.restype = ctypes.c_uint32
        _MachTaskBasicInfo = _MachTaskBasicInfoStruct
        _MACH_TASK_INFO_READY = True
    except Exception:
        pass


def _snap_rss_current() -> float:
    """Snapshot current RSS in MB.

    Uses /proc/self/status on Linux (VmRSS) or task_info on macOS via
    ctypes mach API.  Falls back to peak RSS if current is unavailable.
    """
    # Linux: read VmRSS from proc
    if sys.platform == "linux":
        try:
            with open("/proc/self/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        return int(line.split()[1]) / 1024  # KB → MB
        except FileNotFoundError:
            pass

    # macOS: use pre-initialized mach task_info
    if _MACH_TASK_INFO_READY:
        try:
            info = _MachTaskBasicInfo()
            count = ctypes.c_uint32(ctypes.sizeof(info) // 4)
            task = _mach_libc.mach_task_self()
            ret = _mach_libc.task_info(
                task, _MACH_TASK_BASIC_INFO, ctypes.byref(info), ctypes.byref(count),
            )
            if ret == 0:
                return info.resident_size / (1024 * 1024)
        except Exception:
            pass

    # Fallback: peak RSS
    return _get_rss_mb()


# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════


@dataclass
class HulkrnaParams:
    """rigel parameter set (same keys as benchmark.py HulkrnaConfig)."""

    name: str = "default"
    params: dict = field(default_factory=dict)


@dataclass
class ProfileConfig:
    """Profiling configuration."""

    bam: str = ""
    index: str = ""
    outdir: str = "profile_output"

    rigel_configs: list[HulkrnaParams] = field(default_factory=list)

    stages: bool = False
    enable_cprofile: bool = False
    cprofile_top_n: int = 60
    memory_sample_interval_ms: int = 100
    verbose: bool = True
    tmpdir: str | None = None


def parse_yaml_config(path: str | Path) -> ProfileConfig:
    """Parse a YAML config file into ProfileConfig."""
    if yaml is None:
        raise ImportError("PyYAML required: pip install pyyaml")
    with open(path) as f:
        raw = yaml.safe_load(f)

    cfg = ProfileConfig()
    cfg.bam = raw.get("bam", "")
    cfg.index = raw.get("index", "")
    cfg.outdir = raw.get("outdir", "profile_output")
    cfg.stages = bool(raw.get("stages", False))
    cfg.enable_cprofile = bool(raw.get("enable_cprofile", False))
    cfg.cprofile_top_n = int(raw.get("cprofile_top_n", 60))
    cfg.memory_sample_interval_ms = int(raw.get("memory_sample_interval_ms", 100))
    cfg.verbose = bool(raw.get("verbose", True))

    for name, params in raw.get("rigel_configs", {"default": {}}).items():
        cfg.rigel_configs.append(HulkrnaParams(name=name, params=params or {}))

    if not cfg.rigel_configs:
        cfg.rigel_configs.append(HulkrnaParams())

    return cfg


# ═══════════════════════════════════════════════════════════════════
# Pipeline config builder (shared with benchmark.py)
# ═══════════════════════════════════════════════════════════════════


def _build_pipeline_config(
    rigel_params: HulkrnaParams | None = None,
    tmpdir: str | None = None,
) -> PipelineConfig:
    """Build PipelineConfig from profile config."""
    raw: dict = {}
    if rigel_params is not None:
        raw.update(rigel_params.params)

    em_kw: dict = {}
    _EM_ALIASES = {
        "em_convergence_delta": "convergence_delta",
        "em_iterations": "iterations",
        "em_prior_alpha": "prior_alpha",
        "em_prior_gamma": "prior_gamma",
        "em_mode": "mode",
        "prune_threshold": "prune_threshold",
        "confidence_threshold": "confidence_threshold",
        "n_threads": "n_threads",
    }
    for raw_key, cfg_key in _EM_ALIASES.items():
        if raw_key in raw:
            em_kw[cfg_key] = raw.pop(raw_key)

    scoring_kw: dict = {}
    ov_alpha = raw.pop("overhang_alpha", None)
    if ov_alpha is not None:
        scoring_kw["overhang_log_penalty"] = overhang_alpha_to_log_penalty(ov_alpha)
    mm_alpha = raw.pop("mismatch_alpha", None)
    if mm_alpha is not None:
        scoring_kw["mismatch_log_penalty"] = overhang_alpha_to_log_penalty(mm_alpha)
    gdna_pen = raw.pop("gdna_splice_penalty_unannot", None)
    if gdna_pen is not None:
        penalties = dict(GDNA_SPLICE_PENALTIES)
        penalties[SPLICE_UNANNOT] = gdna_pen
        scoring_kw["gdna_splice_penalties"] = penalties

    scan_kw: dict = {"sj_strand_tag": "auto", "include_multimap": True}
    if tmpdir is not None:
        from pathlib import Path as _Path
        scan_kw["spill_dir"] = _Path(tmpdir)
    if "n_threads" in em_kw:
        scan_kw["n_scan_threads"] = em_kw["n_threads"]

    return PipelineConfig(
        em=EMConfig(**em_kw),
        scan=BamScanConfig(**scan_kw),
        scoring=FragmentScoringConfig(**scoring_kw),
    )


# ═══════════════════════════════════════════════════════════════════
# Timer helper
# ═══════════════════════════════════════════════════════════════════


class Timer:
    """Context-manager wall-clock timer."""

    def __init__(self, label: str):
        self.label = label
        self.elapsed = 0.0

    def __enter__(self):
        self._start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self._start


# ═══════════════════════════════════════════════════════════════════
# Profile run result
# ═══════════════════════════════════════════════════════════════════


@dataclass
class StageTimings:
    """Per-stage wall-clock timings (seconds)."""

    index_load: float = 0.0
    scan_and_buffer: float = 0.0
    finalize_models: float = 0.0
    quant_from_buffer: float = 0.0

    # Sub-stages of quant_from_buffer (if --stages)
    compute_geometry: float = 0.0
    create_estimator: float = 0.0
    fragment_scorer: float = 0.0
    fragment_router_scan: float = 0.0
    compute_nrna_init: float = 0.0
    build_loci: float = 0.0
    compute_eb_gdna_priors: float = 0.0
    locus_em: float = 0.0
    locus_em_build: float = 0.0
    locus_em_run: float = 0.0
    locus_em_assign: float = 0.0


@dataclass
class ProfileResult:
    """Complete profile run result."""

    config_name: str
    bam_path: str
    index_path: str

    wall_time_sec: float
    peak_rss_mb: float
    rss_before_mb: float
    rss_after_mb: float

    n_fragments: int
    n_buffered: int
    n_unique: int
    n_multimapping: int
    n_intergenic: int
    n_gdna_total: int
    n_loci: int
    throughput_frags_per_sec: float

    stages: StageTimings
    pipeline_stats: dict

    # Index metadata
    n_transcripts: int
    n_genes: int

    # Locus summary
    max_locus_transcripts: int = 0
    max_locus_units: int = 0


# ═══════════════════════════════════════════════════════════════════
# Simple profile: run_pipeline() as a whole
# ═══════════════════════════════════════════════════════════════════


def profile_simple(
    bam_path: str,
    index: TranscriptIndex,
    config_name: str,
    rigel_params: HulkrnaParams | None = None,
    enable_cprofile: bool = False,
    tmpdir: str | None = None,
) -> tuple[ProfileResult, cProfile.Profile | None]:
    """Profile run_pipeline() as a single timed call."""

    pcfg = _build_pipeline_config(rigel_params, tmpdir=tmpdir)
    profiler = cProfile.Profile() if enable_cprofile else None

    rss_before = _get_rss_mb()

    t0 = time.perf_counter()
    if profiler:
        profiler.enable()
    pipe = run_pipeline(bam_path, index, config=pcfg)
    if profiler:
        profiler.disable()
    wall = time.perf_counter() - t0

    rss_after = _get_rss_mb()

    stats = pipe.stats
    result = ProfileResult(
        config_name=config_name,
        bam_path=bam_path,
        index_path="",
        wall_time_sec=wall,
        peak_rss_mb=rss_after,
        rss_before_mb=rss_before,
        rss_after_mb=rss_after,
        n_fragments=stats.total,
        n_buffered=stats.n_fragments,
        n_unique=stats.unique,
        n_multimapping=stats.multimapping,
        n_intergenic=stats.n_intergenic,
        n_gdna_total=stats.n_gdna_total,
        n_loci=0,
        throughput_frags_per_sec=stats.total / wall if wall > 0 else 0,
        stages=StageTimings(),
        pipeline_stats=stats.to_dict(),
        n_transcripts=index.num_transcripts,
        n_genes=index.num_genes,
    )
    return result, profiler


# ═══════════════════════════════════════════════════════════════════
# Stage-decomposed profile
# ═══════════════════════════════════════════════════════════════════


def profile_stages(
    bam_path: str,
    index: TranscriptIndex,
    config_name: str,
    rigel_params: HulkrnaParams | None = None,
    enable_cprofile: bool = False,
    tmpdir: str | None = None,
) -> tuple[ProfileResult, cProfile.Profile | None]:
    """Profile pipeline with per-stage timing decomposition."""

    pcfg = _build_pipeline_config(rigel_params, tmpdir=tmpdir)
    profiler = cProfile.Profile() if enable_cprofile else None
    timings = StageTimings()

    rss_before = _get_rss_mb()

    # ── Stage 1: scan_and_buffer ────────────────────────────
    if profiler:
        profiler.enable()
    with Timer("scan_and_buffer") as t_scan:
        stats, strand_models, frag_length_models, buffer = scan_and_buffer(
            bam_path, index, pcfg.scan,
        )
    timings.scan_and_buffer = t_scan.elapsed

    # ── Stage 2: Finalize models ────────────────────────────
    with Timer("finalize_models") as t_fin:
        strand_models.finalize()
        frag_length_models.finalize()
    timings.finalize_models = t_fin.elapsed

    # ── Stage 3: quant_from_buffer (decomposed) ─────────────

    em_config = pcfg.em
    scoring = pcfg.scoring

    # 3a: Compute geometry
    with Timer("compute_geometry") as t_geom:
        exonic_lengths = index.t_df["length"].values.astype(np.float64)
        if frag_length_models.global_model.n_observations > 0:
            mean_frag = frag_length_models.global_model.mean
        else:
            mean_frag = 200.0

        if frag_length_models.global_model.n_observations > 0:
            effective_lengths = (
                frag_length_models.global_model.compute_all_transcript_eff_lens(
                    exonic_lengths.astype(np.int64),
                )
            )
        else:
            effective_lengths = np.maximum(exonic_lengths - mean_frag + 1.0, 1.0)

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
    timings.compute_geometry = t_geom.elapsed

    # 3b: TSS groups + create estimator
    with Timer("create_estimator") as t_est:
        if index.t_to_tss_group is None:
            index.compute_and_set_tss_groups(tss_window=em_config.tss_window)
        estimator = AbundanceEstimator(
            index.num_transcripts,
            em_config=em_config,
            geometry=geometry,
        )
    timings.create_estimator = t_est.elapsed

    # 3c: FragmentScorer
    with Timer("fragment_scorer") as t_scorer:
        ctx = FragmentScorer.from_models(
            strand_models, frag_length_models, index, estimator,
            overhang_log_penalty=scoring.overhang_log_penalty,
            mismatch_log_penalty=scoring.mismatch_log_penalty,
            gdna_splice_penalties=scoring.gdna_splice_penalties,
        )
    timings.fragment_scorer = t_scorer.elapsed

    # 3d: FragmentRouter.scan (scoring + routing)
    with Timer("fragment_router_scan") as t_route:
        builder = FragmentRouter(ctx, estimator, stats, index, strand_models)
        em_data = builder.scan(buffer, log_every=1_000_000)
    timings.fragment_router_scan = t_route.elapsed

    # -- Free scanner accumulators + buffer after scan --
    del builder, ctx          # release ~1.3 GB of array.array accumulators
    buffer.release()
    gc.collect()

    # 3e: nRNA init
    with Timer("compute_nrna_init") as t_nrna:
        nrna_init = compute_nrna_init(
            estimator.transcript_intronic_sense,
            estimator.transcript_intronic_antisense,
            geometry.transcript_spans,
            geometry.exonic_lengths,
            geometry.mean_frag,
            strand_models,
        )
        estimator.nrna_init = nrna_init
    timings.compute_nrna_init = t_nrna.elapsed

    # 3f–3h: Locus-level EM
    n_loci = 0
    max_locus_t = 0
    max_locus_u = 0

    if em_data.n_units > 0:
        with Timer("build_loci") as t_loci:
            loci = build_loci(em_data, index)
        timings.build_loci = t_loci.elapsed
        n_loci = len(loci) if loci else 0

        if loci:
            max_locus_t = max(len(l.transcript_indices) for l in loci)
            max_locus_u = max(len(l.unit_indices) for l in loci)

            # Assign locus_id to every transcript
            for locus in loci:
                for t_idx in locus.transcript_indices:
                    estimator.locus_id_per_transcript[int(t_idx)] = locus.locus_id

            with Timer("compute_eb_gdna_priors") as t_gdna:
                intergenic_density = _compute_intergenic_density(stats, index)
                gdna_inits = compute_eb_gdna_priors(
                    loci, em_data, estimator, index, strand_models,
                    intergenic_density=intergenic_density,
                    kappa_chrom=em_config.gdna_kappa_chrom,
                    kappa_locus=em_config.gdna_kappa_locus,
                    mom_min_evidence_chrom=em_config.gdna_mom_min_evidence_chrom,
                    mom_min_evidence_locus=em_config.gdna_mom_min_evidence_locus,
                    kappa_min=em_config.nrna_frac_kappa_min,
                    kappa_max=em_config.nrna_frac_kappa_max,
                    kappa_fallback=em_config.nrna_frac_kappa_fallback,
                    kappa_min_obs=em_config.nrna_frac_kappa_min_obs,
                )

                # gDNA density + nrna_frac priors (needed before EM)
                gdna_density = compute_global_gdna_density(
                    estimator, strand_models.strand_specificity,
                )
                compute_hybrid_nrna_frac_priors(
                    estimator,
                    t_to_tss_group=index.t_to_tss_group,
                    t_to_strand=index.t_to_strand_arr,
                    locus_id_per_transcript=estimator.locus_id_per_transcript,
                    strand_specificity=strand_models.strand_specificity,
                    gdna_density=gdna_density,
                    kappa_global=em_config.nrna_frac_kappa_global,
                    kappa_locus=em_config.nrna_frac_kappa_locus,
                    kappa_tss=em_config.nrna_frac_kappa_tss,
                    mom_min_evidence_global=em_config.nrna_frac_mom_min_evidence_global,
                    mom_min_evidence_locus=em_config.nrna_frac_mom_min_evidence_locus,
                    mom_min_evidence_tss=em_config.nrna_frac_mom_min_evidence_tss,
                    kappa_min=em_config.nrna_frac_kappa_min,
                    kappa_max=em_config.nrna_frac_kappa_max,
                    kappa_fallback=em_config.nrna_frac_kappa_fallback,
                    kappa_min_obs=em_config.nrna_frac_kappa_min_obs,
                )
            timings.compute_eb_gdna_priors = t_gdna.elapsed

            with Timer("locus_em") as t_em:
                # Use batch C++ path (same as quant_from_buffer default)
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
            timings.locus_em = t_em.elapsed
            timings.locus_em_build = 0.0
            timings.locus_em_run = timings.locus_em
            timings.locus_em_assign = 0.0

    if profiler:
        profiler.disable()

    # -- Phase 5B: Free ScoredFragments CSR arrays --
    del em_data
    gc.collect()

    timings.quant_from_buffer = (
        t_geom.elapsed + t_est.elapsed + t_scorer.elapsed + t_route.elapsed
        + t_nrna.elapsed + timings.build_loci + timings.compute_eb_gdna_priors
        + timings.locus_em
    )

    rss_after = _get_rss_mb()
    total_wall = timings.scan_and_buffer + timings.finalize_models + timings.quant_from_buffer

    result = ProfileResult(
        config_name=config_name,
        bam_path=bam_path,
        index_path="",
        wall_time_sec=total_wall,
        peak_rss_mb=rss_after,
        rss_before_mb=rss_before,
        rss_after_mb=rss_after,
        n_fragments=stats.total,
        n_buffered=stats.n_fragments,
        n_unique=stats.unique,
        n_multimapping=stats.multimapping,
        n_intergenic=stats.n_intergenic,
        n_gdna_total=stats.n_gdna_total,
        n_loci=n_loci,
        throughput_frags_per_sec=stats.total / total_wall if total_wall > 0 else 0,
        stages=timings,
        pipeline_stats=stats.to_dict(),
        n_transcripts=index.num_transcripts,
        n_genes=index.num_genes,
        max_locus_transcripts=max_locus_t,
        max_locus_units=max_locus_u,
    )
    return result, profiler


# ═══════════════════════════════════════════════════════════════════
# Report formatting
# ═══════════════════════════════════════════════════════════════════


def _format_bar(pct: float, width: int = 40) -> str:
    filled = int(pct / 100 * width)
    return "█" * filled + "░" * (width - filled)


def format_report(results: list[ProfileResult], stage_mode: bool) -> str:
    """Format profile results as a human-readable text report."""
    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("HULKRNA PROFILE REPORT")
    lines.append("=" * 72)
    lines.append(f"Platform:   {platform.platform()}")
    lines.append(f"Python:     {sys.version.split()[0]}")
    lines.append(f"Date:       {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    for r in results:
        lines.append("-" * 72)
        lines.append(f"Config: {r.config_name}")
        lines.append("-" * 72)
        lines.append(f"  BAM:             {r.bam_path}")
        lines.append(f"  Index:           {r.index_path}")
        lines.append(f"  Transcripts:     {r.n_transcripts:,}")
        lines.append(f"  Genes:           {r.n_genes:,}")
        lines.append("")
        lines.append(f"  Wall time:       {r.wall_time_sec:.2f}s")
        lines.append(f"  Peak RSS:        {r.peak_rss_mb:.0f} MB")
        lines.append(f"  RSS before:      {r.rss_before_mb:.0f} MB")
        lines.append(f"  RSS after:       {r.rss_after_mb:.0f} MB")
        lines.append(f"  RSS delta:       {r.rss_after_mb - r.rss_before_mb:.0f} MB")
        lines.append("")
        lines.append(f"  Fragments:       {r.n_fragments:,}")
        lines.append(f"  Buffered:        {r.n_buffered:,}")
        lines.append(f"  Unique align:    {r.n_unique:,}")
        lines.append(f"  Multimapping:    {r.n_multimapping:,}")
        lines.append(f"  Intergenic:      {r.n_intergenic:,}")
        lines.append(f"  gDNA total:      {r.n_gdna_total:,}")
        lines.append(f"  Throughput:      {r.throughput_frags_per_sec:,.0f} frags/s")
        lines.append("")

        if stage_mode:
            s = r.stages
            stages = [
                ("scan_and_buffer", s.scan_and_buffer),
                ("finalize_models", s.finalize_models),
                ("compute_geometry", s.compute_geometry),
                ("create_estimator", s.create_estimator),
                ("fragment_scorer", s.fragment_scorer),
                ("fragment_router_scan", s.fragment_router_scan),
                ("compute_nrna_init", s.compute_nrna_init),
                ("build_loci", s.build_loci),
                ("eb_gdna_priors", s.compute_eb_gdna_priors),
                ("locus_em", s.locus_em),
            ]
            total = sum(t for _, t in stages)
            lines.append("  Stage Breakdown:")
            lines.append(f"  {'Stage':<24s} {'Time':>8s}  {'%':>6s}  Bar")
            lines.append(f"  {'-'*24} {'-'*8}  {'-'*6}  {'-'*40}")
            for name, dur in stages:
                pct = dur / total * 100 if total > 0 else 0
                bar = _format_bar(pct)
                lines.append(f"  {name:<24s} {dur:7.3f}s  {pct:5.1f}%  {bar}")
            lines.append(f"  {'TOTAL':<24s} {total:7.3f}s")
            lines.append("")

            if s.locus_em > 0:
                lines.append("  Locus EM sub-stages:")
                lines.append(f"    build_locus_em_data: {s.locus_em_build:.3f}s")
                lines.append(f"    run_locus_em:        {s.locus_em_run:.3f}s")
                lines.append(f"    assign_posteriors:   {s.locus_em_assign:.3f}s")
                lines.append(f"    Loci:                {r.n_loci:,}")
                lines.append(f"    Max transcripts:     {r.max_locus_transcripts:,}")
                lines.append(f"    Max units:           {r.max_locus_units:,}")
                lines.append("")

            # Top-3 hotspots
            ranked = sorted(stages, key=lambda x: -x[1])
            lines.append("  Top-3 hotspots:")
            for i, (name, dur) in enumerate(ranked[:3]):
                pct = dur / total * 100 if total > 0 else 0
                lines.append(f"    {i + 1}. {name}: {dur:.3f}s ({pct:.1f}%)")
            lines.append("")

    # Comparison table (if multiple configs)
    if len(results) > 1:
        lines.append("=" * 72)
        lines.append("COMPARISON")
        lines.append("=" * 72)
        lines.append(f"  {'Config':<20s} {'Wall (s)':>10s} {'RSS (MB)':>10s} "
                      f"{'Frags':>12s} {'Throughput':>15s}")
        lines.append(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*12} {'-'*15}")
        for r in results:
            lines.append(
                f"  {r.config_name:<20s} {r.wall_time_sec:10.2f} "
                f"{r.peak_rss_mb:10.0f} {r.n_fragments:12,} "
                f"{r.throughput_frags_per_sec:15,.0f}"
            )
        lines.append("")

    return "\n".join(lines)


def format_cprofile(profiler: cProfile.Profile, top_n: int = 60) -> str:
    """Format cProfile stats as text."""
    lines: list[str] = []

    s_cum = io.StringIO()
    ps = pstats.Stats(profiler, stream=s_cum)
    ps.sort_stats("cumulative")
    ps.print_stats(top_n)
    lines.append(f"{'='*72}")
    lines.append(f"Top {top_n} by cumulative time")
    lines.append(f"{'='*72}")
    lines.append(s_cum.getvalue())

    s_tot = io.StringIO()
    ps2 = pstats.Stats(profiler, stream=s_tot)
    ps2.sort_stats("tottime")
    ps2.print_stats(top_n)
    lines.append(f"{'='*72}")
    lines.append(f"Top {top_n} by self time")
    lines.append(f"{'='*72}")
    lines.append(s_tot.getvalue())

    total_calls = sum(v[0] for v in ps.stats.values())
    lines.append(f"Total function calls: {total_calls:,}")

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════


def run_profile(cfg: ProfileConfig) -> list[ProfileResult]:
    """Run profiling for all configured parameter sets."""
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam_path = cfg.bam
    index_dir = cfg.index

    if not Path(bam_path).exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    if not Path(index_dir).is_dir():
        raise FileNotFoundError(f"Index not found: {index_dir}")

    # ── Load index (timed separately) ───────────────────────
    print(f"Loading index from {index_dir} ...", flush=True)
    t0 = time.perf_counter()
    index = TranscriptIndex.load(str(index_dir))
    index_load_time = time.perf_counter() - t0
    print(f"  {index.num_transcripts:,} transcripts, "
          f"{index.num_genes:,} genes ({index_load_time:.1f}s)", flush=True)

    results: list[ProfileResult] = []
    all_profilers: dict[str, cProfile.Profile] = {}

    for hc in cfg.rigel_configs:
        print(f"\n{'='*72}", flush=True)
        print(f"Profiling: {hc.name}", flush=True)
        print(f"  BAM: {Path(bam_path).name}", flush=True)
        if hc.params:
            print(f"  Params: {hc.params}", flush=True)
        print(f"{'='*72}\n", flush=True)

        # Start memory timeline
        mem_interval = cfg.memory_sample_interval_ms / 1000.0
        mem_timeline = MemoryTimeline(interval_sec=mem_interval)
        mem_timeline.start()

        try:
            if cfg.stages:
                result, profiler = profile_stages(
                    bam_path, index, hc.name, hc, cfg.enable_cprofile,
                    tmpdir=cfg.tmpdir,
                )
            else:
                result, profiler = profile_simple(
                    bam_path, index, hc.name, hc, cfg.enable_cprofile,
                    tmpdir=cfg.tmpdir,
                )
        finally:
            mem_timeline.stop()

        result.index_path = str(index_dir)
        result.stages.index_load = index_load_time
        result.peak_rss_mb = max(result.peak_rss_mb, mem_timeline.peak_mb)

        results.append(result)
        if profiler:
            all_profilers[hc.name] = profiler

        # Write memory timeline per config
        mem_csv = outdir / f"memory_timeline_{hc.name}.csv"
        mem_timeline.write_csv(mem_csv)

        # Print inline summary
        print(f"\nWall time: {result.wall_time_sec:.2f}s", flush=True)
        print(f"Fragments: {result.n_fragments:,}", flush=True)
        print(f"Throughput: {result.throughput_frags_per_sec:,.0f} frags/s", flush=True)
        print(f"Peak RSS: {result.peak_rss_mb:.0f} MB", flush=True)

    # ── Write reports ────────────────────────────────────────
    report_txt = format_report(results, stage_mode=cfg.stages)

    # Append cProfile output
    for name, profiler in all_profilers.items():
        report_txt += f"\n\n{'='*72}\n"
        report_txt += f"cProfile: {name}\n"
        report_txt += f"{'='*72}\n"
        report_txt += format_cprofile(profiler, cfg.cprofile_top_n)
        # Save .prof file
        prof_path = outdir / f"profile_{name}.prof"
        profiler.dump_stats(str(prof_path))
        print(f"Saved cProfile data: {prof_path}", flush=True)

    # Write text report
    report_path = outdir / "profile_report.txt"
    with open(report_path, "w") as f:
        f.write(report_txt)
    print(f"\nReport: {report_path}", flush=True)

    # Print report to stdout
    print(f"\n{report_txt}", flush=True)

    # Write JSON summary
    json_data = {
        "platform": platform.platform(),
        "python": sys.version,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "index_load_sec": index_load_time,
        "profiles": [
            {
                "config_name": r.config_name,
                "bam_path": r.bam_path,
                "index_path": r.index_path,
                "wall_time_sec": r.wall_time_sec,
                "peak_rss_mb": r.peak_rss_mb,
                "rss_before_mb": r.rss_before_mb,
                "rss_after_mb": r.rss_after_mb,
                "n_fragments": r.n_fragments,
                "n_buffered": r.n_buffered,
                "n_unique": r.n_unique,
                "n_multimapping": r.n_multimapping,
                "n_intergenic": r.n_intergenic,
                "n_gdna_total": r.n_gdna_total,
                "n_loci": r.n_loci,
                "throughput_frags_per_sec": r.throughput_frags_per_sec,
                "n_transcripts": r.n_transcripts,
                "n_genes": r.n_genes,
                "stages": asdict(r.stages),
                "pipeline_stats": r.pipeline_stats,
            }
            for r in results
        ],
    }
    json_path = outdir / "profile_summary.json"
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"JSON:   {json_path}", flush=True)

    return results


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Profile rigel pipeline performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/profiler.py --bam reads.bam --index rigel_index/\n"
            "  python scripts/profiler.py --config scripts/profile_example.yaml --stages\n"
            "  python scripts/profiler.py --bam reads.bam --index rigel_index/ "
            "--cprofile --stages\n"
        ),
    )
    p.add_argument("--config", help="YAML configuration file")
    p.add_argument("--bam", help="Name-sorted BAM file (overrides YAML)")
    p.add_argument("--index", help="rigel index directory (overrides YAML)")
    p.add_argument("--outdir", help="Output directory (overrides YAML)")
    p.add_argument("--stages", action="store_true", default=None,
                    help="Enable per-stage timing decomposition")
    p.add_argument("--cprofile", action="store_true", default=None,
                    help="Enable cProfile function-level profiling")
    p.add_argument("--memory-interval", type=int, default=None,
                    help="RSS sampling interval in ms (default: 100)")
    p.add_argument("--verbose", action="store_true", default=None,
                    help="Verbose logging")
    p.add_argument("--threads", type=int, default=None,
                    help="Number of threads for parallel locus EM "
                         "(0=all cores, 1=sequential)")
    p.add_argument("--tmpdir", default=None,
                    help="Directory for temporary buffer spill files "
                         "(default: system temp directory)")
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    # Build config from YAML or defaults
    if args.config:
        cfg = parse_yaml_config(args.config)
    else:
        cfg = ProfileConfig()
        if not cfg.rigel_configs:
            cfg.rigel_configs.append(HulkrnaParams())

    # CLI overrides
    if args.bam:
        cfg.bam = args.bam
    if args.index:
        cfg.index = args.index
    if args.outdir:
        cfg.outdir = args.outdir
    if args.stages is not None:
        cfg.stages = args.stages
    if args.cprofile is not None:
        cfg.enable_cprofile = args.cprofile
    if args.memory_interval is not None:
        cfg.memory_sample_interval_ms = args.memory_interval
    if args.verbose is not None:
        cfg.verbose = args.verbose
    if args.threads is not None:
        for hc in cfg.rigel_configs:
            hc.params["n_threads"] = args.threads
    if args.tmpdir is not None:
        cfg.tmpdir = args.tmpdir

    # Validate
    if not cfg.bam:
        print("Error: --bam is required (or specify bam in YAML)", file=sys.stderr)
        return 1
    if not cfg.index:
        print("Error: --index is required (or specify index in YAML)", file=sys.stderr)
        return 1

    level = logging.DEBUG if cfg.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    print("rigel profiler", flush=True)
    print(f"  BAM:             {cfg.bam}", flush=True)
    print(f"  Index:           {cfg.index}", flush=True)
    print(f"  Output:          {cfg.outdir}", flush=True)
    print(f"  Stages:          {cfg.stages}", flush=True)
    print(f"  cProfile:        {cfg.enable_cprofile}", flush=True)
    print(f"  Memory interval: {cfg.memory_sample_interval_ms}ms", flush=True)
    print(f"  Configs:         {[h.name for h in cfg.rigel_configs]}", flush=True)

    results = run_profile(cfg)

    print(f"\nProfiling complete. {len(results)} config(s) profiled.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
