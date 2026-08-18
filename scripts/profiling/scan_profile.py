#!/usr/bin/env python3
"""WHAT DOES THE SCAN ALONE COST, AND HOW DOES IT ANSWER TO ITS THREAD AND CHUNK KNOBS?

This harness loads a Rigel index once, then runs only ``scan_and_buffer`` across a grid of
scan/decompression thread budgets, reporting wall time, throughput and RSS for each. It is also the
command to launch under macOS Instruments/xctrace or py-spy, because native samples then land almost
entirely inside ``rigel._bam_impl`` and htslib rather than in a harness.

⭐ It is the only instrument here that SWEEPS a knob, and its shape is the one a memory-bounded
chunked solve wants: name a knob, sweep it, and put wall clock and peak RSS in the same row.
⚠ Its knobs are the scanner's (``--threads``, ``--scan-bgzf-threads``, ``--scan-fragments-per-chunk``,
``--scan-buffer-size``); it reaches no calibration or EM knob, and ``profiler.py`` is the instrument
that sees those phases at all.

⛔⛔ **PROFILE A HIGH-DEPTH REAL RNA-seq LIBRARY, NOT cfRNA AND NOT A PANEL CONDITION** (owner,
2026-08-17). `docs/TESTING.md` §7 is the ruling.

⛔⛔ **TWO FILES ARE CALLED ``scan_profile.py`` AND THEY MEASURE DIFFERENT THINGS.**
``scripts/design/scan_profile.py`` regresses the ACCUMULATOR's cost per fragment over several BAMs
(wired vs unwired, slope and intercept). This one measures the SCAN itself. Neither is a version of
the other; a reader who greps the basename finds both, and the gate that covers them prints the tree
in the case id so a failure is never ambiguous.
"""

from __future__ import annotations

import argparse
import csv
import gc
import itertools
import json
import logging
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

if str(Path(__file__).resolve().parent) not in sys.path:  # `profiler.py` is a sibling, not a package
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from profiler import MemoryTimeline, _peak_rss_mb, _snap_rss_current  # noqa: E402
from rigel.config import BamScanConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.native import detect_sj_strand_tag  # noqa: E402
from rigel.pipeline import scan_and_buffer  # noqa: E402

logger = logging.getLogger("rigel_scan_profile")


@dataclass
class ScanRun:
    """One scan-only profiling result."""

    name: str
    bam_path: str
    index_path: str
    wall_time_sec: float
    #: ⭐ the largest CURRENT-RSS sample inside THIS run's span — falls when memory is released
    peak_rss_mb: float
    #: ⚠ ``ru_maxrss``: exact, but a LIFETIME high-water mark, so in a sweep it can only rise. Reported
    #: beside the sampled peak rather than instead of it, because the two answer different questions.
    process_peak_rss_mb: float
    rss_before_mb: float
    rss_after_scan_mb: float
    rss_after_release_mb: float
    threads: int
    scan_worker_threads: int
    scan_bgzf_threads: int
    requested_scan_bgzf_threads: int
    scan_fragments_per_chunk: int
    scan_read_name_batch_size: int
    scan_buffer_size_bytes: int
    n_bam_records: int
    n_read_names: int
    n_physical_fragments: int
    n_buffered_fragments: int
    n_unique: int
    n_multimapping: int
    n_intergenic: int
    n_chimeric: int
    buffer_memory_mb: float
    buffer_chunks: int
    buffer_spilled_chunks: int
    bam_records_per_sec: float
    read_names_per_sec: float
    buffered_fragments_per_sec: float
    memory_timeline_csv: str


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Profile Rigel scan_and_buffer without running quantification.",
    )
    parser.add_argument("--bam", required=True, type=Path, help="Name-sorted input BAM")
    parser.add_argument("--index", required=True, type=Path, help="Rigel index directory")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--name-prefix", default="scan", help="Prefix for run names")
    parser.add_argument(
        "--threads",
        nargs="+",
        type=int,
        default=[8],
        help="One or more total scan thread budgets",
    )
    parser.add_argument(
        "--scan-bgzf-threads",
        nargs="+",
        type=int,
        default=[4],
        help="One or more BGZF decompression thread counts reserved from --threads",
    )
    parser.add_argument("--repeat", type=int, default=1, help="Repeat every configuration")
    parser.add_argument("--scan-fragments-per-chunk", type=int, default=1_000_000)
    parser.add_argument("--scan-read-name-batch-size", type=int, default=512)
    parser.add_argument("--scan-buffer-size", type=float, default=2.0)
    parser.add_argument("--spill-dir", type=Path, default=None)
    parser.add_argument("--max-frag-length", type=int, default=1000)
    parser.add_argument("--sj-strand-tag", default="auto")
    parser.add_argument("--splicing-anchor-tolerance", type=int, default=3)
    parser.add_argument("--include-multimap", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--skip-duplicates", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--memory-interval-ms", type=int, default=250)
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def _run_one(
    *,
    name: str,
    bam_path: Path,
    index_path: Path,
    index: TranscriptIndex,
    outdir: Path,
    scan_cfg: BamScanConfig,
    memory_interval_sec: float,
) -> ScanRun:
    run_dir = outdir / name
    run_dir.mkdir(parents=True, exist_ok=True)

    gc.collect()
    rss_before = _snap_rss_current()
    timeline = MemoryTimeline(interval_sec=memory_interval_sec)

    scan_workers, bgzf_threads = scan_cfg.resolved_scan_threads()
    logger.info(
        "Starting %s: threads=%d scan_workers=%d scan_bgzf_threads=%d "
        "scan_fragments_per_chunk=%d scan_read_name_batch_size=%d "
        "scan_buffer_size=%.2f GiB",
        name,
        scan_cfg.resolved_total_threads(),
        scan_workers,
        bgzf_threads,
        scan_cfg.fragments_per_chunk,
        scan_cfg.read_name_batch_size,
        scan_cfg.buffer_size_bytes / 1024**3,
    )
    timeline.start()
    if not timeline.available:
        # ⛔ Say it, per run, rather than printing `nan` columns that read as a measurement.
        logger.warning(
            "no CURRENT-RSS reader on this platform (%s): every rss_* column in %s is `nan`. "
            "Wall time and throughput are unaffected.",
            platform.system(),
            name,
        )
    t0 = time.perf_counter()
    stats, strand_models, buffer, cal_payload = scan_and_buffer(
        str(bam_path),
        index,
        scan_cfg,
    )
    wall = time.perf_counter() - t0
    timeline.stop()

    memory_csv = run_dir / "memory_timeline.csv"
    timeline.write_csv(memory_csv)

    buffer_memory_mb = buffer.memory_bytes / 1024**2
    n_buffered = buffer.total_fragments
    n_chunks = buffer.n_chunks
    n_spilled = buffer.n_spilled
    rss_after_scan = _snap_rss_current()
    stats_dict = stats.to_dict()

    # Release large scan objects before the next sweep point.
    buffer.release()
    del buffer, cal_payload, strand_models
    gc.collect()
    rss_after_release = _snap_rss_current()

    run = ScanRun(
        name=name,
        bam_path=str(bam_path),
        index_path=str(index_path),
        wall_time_sec=wall,
        peak_rss_mb=timeline.peak_mb,
        process_peak_rss_mb=_peak_rss_mb(),
        rss_before_mb=rss_before,
        rss_after_scan_mb=rss_after_scan,
        rss_after_release_mb=rss_after_release,
        threads=scan_cfg.total_threads,
        scan_worker_threads=scan_workers,
        scan_bgzf_threads=bgzf_threads,
        requested_scan_bgzf_threads=scan_cfg.bgzf_threads,
        scan_fragments_per_chunk=scan_cfg.fragments_per_chunk,
        scan_read_name_batch_size=scan_cfg.read_name_batch_size,
        scan_buffer_size_bytes=scan_cfg.buffer_size_bytes,
        n_bam_records=int(stats_dict["total"]),
        n_read_names=int(stats_dict["n_read_names"]),
        n_physical_fragments=int(stats_dict["n_fragments"]),
        n_buffered_fragments=int(n_buffered),
        n_unique=int(stats_dict["unique"]),
        n_multimapping=int(stats_dict["multimapping"]),
        n_intergenic=int(stats_dict["n_intergenic"]),
        n_chimeric=int(stats_dict["n_chimeric"]),
        buffer_memory_mb=buffer_memory_mb,
        buffer_chunks=int(n_chunks),
        buffer_spilled_chunks=int(n_spilled),
        bam_records_per_sec=stats_dict["total"] / wall if wall > 0 else 0.0,
        read_names_per_sec=stats_dict["n_read_names"] / wall if wall > 0 else 0.0,
        buffered_fragments_per_sec=n_buffered / wall if wall > 0 else 0.0,
        memory_timeline_csv=str(memory_csv),
    )
    (run_dir / "scan_result.json").write_text(json.dumps(asdict(run), indent=2))
    logger.info(
        "Finished %s: %.2fs, %.0f read-names/s, peak RSS %.1f MB, spilled %d chunks",
        name,
        wall,
        run.read_names_per_sec,
        run.peak_rss_mb,
        run.buffer_spilled_chunks,
    )
    return run


def _write_summary(outdir: Path, summary: dict, runs: list[ScanRun]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "scan_profile_summary.json").write_text(json.dumps(summary, indent=2))
    if not runs:
        return
    with open(outdir / "scan_profile_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(runs[0]).keys()))
        writer.writeheader()
        for run in runs:
            writer.writerow(asdict(run))


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )
    args.outdir.mkdir(parents=True, exist_ok=True)

    index_t0 = time.perf_counter()
    index = TranscriptIndex.load(args.index)
    index_load_sec = time.perf_counter() - index_t0

    sj_tag = args.sj_strand_tag
    if sj_tag == "auto":
        sj_tag = detect_sj_strand_tag(str(args.bam))
        logger.info("Detected splice-junction strand tag: %s", sj_tag)

    scan_buffer_size_bytes = int(args.scan_buffer_size * 1024**3)
    runs: list[ScanRun] = []
    memory_interval_sec = max(args.memory_interval_ms, 1) / 1000.0

    for repeat_index in range(args.repeat):
        for threads, bgzf_threads in itertools.product(
            args.threads,
            args.scan_bgzf_threads,
        ):
            suffix = f"t{threads}_bgzf{bgzf_threads}"
            if args.repeat > 1:
                suffix = f"{suffix}_r{repeat_index + 1}"
            name = f"{args.name_prefix}_{suffix}"
            spill_base = args.spill_dir or (args.outdir / name / "spill")
            scan_cfg = BamScanConfig(
                skip_duplicates=args.skip_duplicates,
                include_multimap=args.include_multimap,
                max_frag_length=args.max_frag_length,
                sj_strand_tag=sj_tag,
                total_threads=threads,
                bgzf_threads=bgzf_threads,
                fragments_per_chunk=args.scan_fragments_per_chunk,
                read_name_batch_size=args.scan_read_name_batch_size,
                buffer_size_bytes=scan_buffer_size_bytes,
                spill_dir=spill_base,
                splicing_anchor_tolerance=args.splicing_anchor_tolerance,
            )
            runs.append(
                _run_one(
                    name=name,
                    bam_path=args.bam,
                    index_path=args.index,
                    index=index,
                    outdir=args.outdir,
                    scan_cfg=scan_cfg,
                    memory_interval_sec=memory_interval_sec,
                )
            )

    summary = {
        "platform": platform.platform(),
        "python": sys.version,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "bam": str(args.bam),
        "index": str(args.index),
        "index_load_sec": index_load_sec,
        "runs": [asdict(run) for run in runs],
    }
    _write_summary(args.outdir, summary, runs)
    logger.info("Wrote scan profile summary to %s", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())