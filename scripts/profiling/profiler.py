#!/usr/bin/env python
"""Rigel pipeline profiler — per-stage wall-time + memory, with optional cProfile.

Profiles ``run_pipeline(bam, index)`` and reports where the time goes across the three pipeline
stages:

  1. **scan**      — ``scan_and_buffer`` (C++ htslib BAM scan, fragment resolution, model
                      training, the fractional accumulator).
  2. **calibrate** — the acyclic gDNA/RNA deconvolution + library hyperparameter fit.
  3. **quant**     — ``quant_from_buffer`` (fragment scoring, routing, per-locus EM).

It works by **wrapping those three functions and running the real `run_pipeline`** — it does not
re-implement the pipeline's internals, so it always profiles the production flow and cannot drift
from the actual stage wiring. Use it to find bottlenecks before optimizing on real RNA-seq data.

Usage:
    python scripts/profiling/profiler.py --bam sample.bam --index index/ \\
        [--threads N] [--repeat 3] [--cprofile] [--out profile.json]

Inside the activated `rigel` conda env. The BAM must be name-sorted with NH tags (as `rigel quant`
requires). For a quick synthetic input, generate an oracle BAM + index with `rigel sim` /
`simulate_suite.py` (see docs/CARRY_FORWARD.md).
"""

from __future__ import annotations

import argparse
import cProfile
import gc
import json
import logging
import platform
import pstats
import resource
import statistics
import sys
import time
from contextlib import contextmanager
from dataclasses import replace
from io import StringIO
from pathlib import Path

import rigel.calibration as _cal_pkg
import rigel.pipeline as _pipeline
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline

# ru_maxrss is bytes on macOS, kilobytes on Linux → MB.
_RSS_TO_MB = (1.0 / 1024**2) if sys.platform == "darwin" else (1.0 / 1024)

# The three pipeline stages, in execution order: (module, attribute, label).
# scan_and_buffer / quant_from_buffer are module-level in rigel.pipeline; calibrate is imported
# locally inside run_pipeline via `from .calibration import calibrate`, so patch it on the package.
_STAGES = (
    (_pipeline, "scan_and_buffer", "scan"),
    (_cal_pkg, "calibrate", "calibrate"),
    (_pipeline, "quant_from_buffer", "quant"),
)


def _peak_rss_mb() -> float:
    """Process peak resident set size (high-water mark) in MB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_TO_MB


@contextmanager
def _timed_stages(record: dict[str, dict]):
    """Wrap each stage function with a timer + peak-RSS snapshot; restore on exit.

    Records ``{label: {seconds, peak_rss_mb}}`` for each stage as the real ``run_pipeline`` calls
    it. Wrappers forward ``*args, **kwargs`` so they are immune to stage-signature changes.
    """
    originals = [(mod, attr, getattr(mod, attr)) for mod, attr, _ in _STAGES]

    def make_timed(orig, label):
        def timed(*args, **kwargs):
            start = time.perf_counter()
            try:
                return orig(*args, **kwargs)
            finally:
                record[label] = {
                    "seconds": time.perf_counter() - start,
                    "peak_rss_mb": _peak_rss_mb(),
                }

        return timed

    try:
        for (mod, attr, _), (_, _, orig) in zip(_STAGES, originals):
            label = next(lbl for m, a, lbl in _STAGES if m is mod and a is attr)
            setattr(mod, attr, make_timed(orig, label))
        yield
    finally:
        for mod, attr, orig in originals:
            setattr(mod, attr, orig)


def profile_once(bam_path: str, index: TranscriptIndex, config: PipelineConfig,
                 *, cprofile: bool = False) -> dict:
    """Run the pipeline once; return a result dict with per-stage timings + totals."""
    gc.collect()
    stages: dict[str, dict] = {}
    profiler = cProfile.Profile() if cprofile else None
    rss_before = _peak_rss_mb()
    start = time.perf_counter()
    if profiler:
        profiler.enable()
    with _timed_stages(stages):
        result = run_pipeline(bam_path, index, config=config)
    if profiler:
        profiler.disable()
    total = time.perf_counter() - start

    staged = sum(s["seconds"] for s in stages.values())
    n_frag = int(getattr(result.stats, "n_fragments", 0) or 0)
    return {
        "total_seconds": total,
        "other_seconds": max(0.0, total - staged),  # setup/teardown outside the three stages
        "stages": stages,
        "n_fragments": n_frag,
        "fragments_per_sec": (n_frag / total) if total > 0 else 0.0,
        "peak_rss_mb": _peak_rss_mb(),
        "rss_before_mb": rss_before,
        "_cprofile": profiler,
    }


def _aggregate(runs: list[dict]) -> dict:
    """Mean (± stdev) across repeated runs, per stage + totals."""
    labels = [lbl for _, _, lbl in _STAGES]

    def ms(vals):
        return {"mean": statistics.fmean(vals),
                "stdev": statistics.stdev(vals) if len(vals) > 1 else 0.0}

    agg = {
        "n_runs": len(runs),
        "total_seconds": ms([r["total_seconds"] for r in runs]),
        "other_seconds": ms([r["other_seconds"] for r in runs]),
        "peak_rss_mb": ms([r["peak_rss_mb"] for r in runs]),
        "fragments_per_sec": ms([r["fragments_per_sec"] for r in runs]),
        "n_fragments": runs[-1]["n_fragments"],
        "stages": {},
    }
    for lbl in labels:
        secs = [r["stages"].get(lbl, {}).get("seconds", 0.0) for r in runs]
        rss = [r["stages"].get(lbl, {}).get("peak_rss_mb", 0.0) for r in runs]
        agg["stages"][lbl] = {"seconds": ms(secs), "peak_rss_mb_cumulative": ms(rss)}
    return agg


def _bar(frac: float, width: int = 36) -> str:
    fill = int(round(max(0.0, min(1.0, frac)) * width))
    return "█" * fill + "·" * (width - fill)


def format_report(agg: dict, label: str) -> str:
    total = agg["total_seconds"]["mean"]
    lines = [
        "",
        "=" * 74,
        f"  RIGEL PROFILE — {label}   ({agg['n_runs']} run(s))",
        "=" * 74,
        f"  fragments: {agg['n_fragments']:,}    "
        f"throughput: {agg['fragments_per_sec']['mean']:,.0f} frag/s    "
        f"peak RSS: {agg['peak_rss_mb']['mean']:,.0f} MB",
        "",
        f"  {'stage':<12} {'seconds':>10} {'± stdev':>9} {'% total':>8}  {'':<36}",
        "  " + "-" * 70,
    ]
    for _, _, lbl in _STAGES:
        st = agg["stages"][lbl]
        sec, sd = st["seconds"]["mean"], st["seconds"]["stdev"]
        pct = (sec / total) if total > 0 else 0.0
        lines.append(f"  {lbl:<12} {sec:>10.3f} {sd:>9.3f} {100*pct:>7.1f}%  {_bar(pct)}")
    other = agg["other_seconds"]["mean"]
    opct = (other / total) if total > 0 else 0.0
    lines.append(f"  {'other':<12} {other:>10.3f} {'':>9} {100*opct:>7.1f}%  {_bar(opct)}")
    lines.append("  " + "-" * 70)
    lines.append(f"  {'TOTAL':<12} {total:>10.3f} {agg['total_seconds']['stdev']:>9.3f} {100.0:>7.1f}%")
    lines.append("=" * 74)
    return "\n".join(lines)


def _build_config(args: argparse.Namespace) -> PipelineConfig:
    """Production-default PipelineConfig, with an optional thread override."""
    cfg = PipelineConfig()
    if args.threads is not None:
        cfg = replace(
            cfg,
            scan=replace(cfg.scan, total_threads=args.threads),
            em=replace(cfg.em, n_threads=args.threads),
        )
    return cfg


def _cprofile_text(prof: cProfile.Profile, top: int) -> str:
    buf = StringIO()
    pstats.Stats(prof, stream=buf).sort_stats("cumulative").print_stats(top)
    return buf.getvalue()


def main() -> int:
    ap = argparse.ArgumentParser(description="Profile the Rigel quantification pipeline.")
    ap.add_argument("--bam", required=True, help="Name-sorted BAM with NH tags.")
    ap.add_argument("--index", required=True, help="Rigel index directory.")
    ap.add_argument("--label", default=None, help="Run label for the report (default: BAM name).")
    ap.add_argument("--threads", type=int, default=None, help="Thread budget (scan + EM).")
    ap.add_argument("--repeat", type=int, default=1, help="Repeat N times and average.")
    ap.add_argument("--cprofile", action="store_true",
                    help="Also collect cProfile (Python hotspots) on the first run.")
    ap.add_argument("--cprofile-top", type=int, default=30, help="Top-N cProfile rows to print.")
    ap.add_argument("--out", default=None, help="Write the JSON report to this path.")
    args = ap.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    bam_path = str(Path(args.bam).resolve())
    label = args.label or Path(args.bam).name
    config = _build_config(args)

    index = TranscriptIndex.load(Path(args.index))
    runs: list[dict] = []
    cprofile_text = None
    for i in range(max(1, args.repeat)):
        res = profile_once(bam_path, index, config, cprofile=(args.cprofile and i == 0))
        if res.get("_cprofile") is not None:
            cprofile_text = _cprofile_text(res.pop("_cprofile"), args.cprofile_top)
        else:
            res.pop("_cprofile", None)
        runs.append(res)

    agg = _aggregate(runs)
    print(format_report(agg, label))
    if cprofile_text:
        print(f"\n  cProfile — top {args.cprofile_top} by cumulative time:\n")
        print(cprofile_text)

    if args.out:
        report = {
            "label": label,
            "bam": bam_path,
            "index": str(Path(args.index).resolve()),
            "threads": args.threads,
            "system": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "cpu_count": __import__("os").cpu_count(),
            },
            "aggregate": agg,
            "runs": runs,
        }
        Path(args.out).write_text(json.dumps(report, indent=2))
        print(f"\n  JSON report → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
