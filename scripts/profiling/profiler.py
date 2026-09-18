#!/usr/bin/env python
"""Where do the time and the memory go? The whole pipeline on one library, as a tree of named stages.

Runs the real `run_pipeline` on one BAM and attributes wall clock and memory to named stages: the index
load, the scan, the second pass, each index-derived geometry build, the fragment-length models, every
calibration stage down to each sweep and each directional pass, and every quant stage down to the locus
EM (or, with `--scan-only`, the scan alone, which is the loop for tuning the scanner's thread and chunk
knobs on a deep library without waiting for calibration). Stages are production functions wrapped at runtime (`PROBES`), patched in every loaded module that
holds a reference, so the profile always follows the production flow and never re-implements it; a probe
whose target has moved is reported as missing, never silently dropped. Each stage reports its calls,
inclusive and self seconds, the peak current RSS sampled while it ran (a reading that can fall, unlike
the process high-water mark), and the RSS it left held on exit. `--compare` joins two reports stage by
stage, which is how a speed-up is judged; whether it changed any number is `rename_identity.py --bam`'s
question, not this one's. Profile the real libraries (the cfRNA ones as smoke tests, a deep cell-line
library for the regime that matters), never a panel condition or a toy, which rank hotspots differently.

Usage::

    python scripts/profiling/profiler.py --bam lib.bam --index idx/ --threads 8 --out lib.json
    python scripts/profiling/profiler.py --bam lib.bam --index idx/ --set calibration.sweep_logodds_step=0.2
    python scripts/profiling/profiler.py --bam lib.bam --index idx/ --scan-only --set scan.total_threads=4
    python scripts/profiling/profiler.py --bam lib.bam --index idx/ --cprofile lib.prof
    python scripts/profiling/profiler.py --compare before.json after.json
    python scripts/profiling/profiler.py --self-test
    sudo py-spy record --native -o lib.svg -- python scripts/profiling/profiler.py --bam lib.bam --index idx/
"""

from __future__ import annotations

import argparse
import cProfile
import ctypes
import dataclasses
import importlib
import inspect
import io
import json
import os
import platform
import pstats
import resource
import subprocess
import sys
import threading
import time
import types
from pathlib import Path

# the one SECTION.FIELD=VALUE parser, shared with the design instruments' ``--set``
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "design"))
from _shared import set_field  # noqa: E402

#: every stage worth naming, as (label, module, attribute path). Nesting is not declared: it is whatever
#: the call stack says at run time, so a stage called from two places appears under both parents. Probe
#: stages, not per-slot kernels: each probed call reads the RSS twice, about 5 µs, which is noise on a
#: stage and a distortion on a function called millions of times.
PROBES: tuple[tuple[str, str, str], ...] = (
    ("index load", "rigel.index", "TranscriptIndex.load"),
    ("run_pipeline", "rigel.pipeline", "run_pipeline"),
    ("scan", "rigel.pipeline", "scan_and_buffer"),
    ("second pass (drain)", "rigel.pipeline", "_drain_side_buffer"),
    ("score held fragments", "rigel.second_pass", "score_held_fragments"),
    ("choose hypotheses", "rigel.second_pass", "choose_hypotheses"),
    ("drain held fragments", "rigel.second_pass", "drain"),
    ("region arrays", "rigel.calibration.region_arrays", "RegionArrays.from_index"),
    ("boundary flags", "rigel.calibration.splice_graph", "build_boundary_flags_array"),
    ("mature wall distances", "rigel.calibration.splice_graph", "build_mature_wall_distances"),
    ("contiguous reach", "rigel.calibration.splice_graph", "build_contiguous_boundary_reach_arrays"),
    ("sj geometry", "rigel.calibration.splice_graph", "build_sj_geometry_arrays"),
    ("region partition", "rigel.calibration.splice_graph", "build_region_partition_arrays"),
    ("sj opportunity", "rigel.calibration.sj_opportunity", "crossing_probability_from_index"),
    ("gdna opportunity", "rigel.calibration.gdna_opportunity", "gdna_opportunity_from_index"),
    ("fl models", "rigel.calibration.fl", "build_fl_models"),
    ("calibrate", "rigel.calibration.calibrate", "calibrate"),
    ("region chain", "rigel.calibration.region_chain", "build_region_chain"),
    ("region geometry", "rigel.calibration.region_geometry", "build_region_geometry"),
    ("region statics", "rigel.calibration.region_geometry", "build_region_statics"),
    ("sweep", "rigel.calibration.sweep", "solve_chain"),
    ("locus block", "rigel.calibration.sweep", "_solve_block"),
    ("own claims (region init)", "rigel.calibration.region_init", "build_region_init"),
    ("policy prepare", "rigel.calibration.messages.transfer", "TransferPolicy.prepare"),
    ("directional pass", "rigel.calibration.sweep", "_pass"),
    ("policy solve", "rigel.calibration.messages.transfer", "_PreparedTransfer.solve"),
    ("backbone checks", "rigel.calibration.sweep", "_check_message"),
    ("psi grid solve", "rigel.calibration.simplex_logodds", "_solve_regions_logodds_all"),
    ("landscape fit", "rigel.calibration.landscape", "fit_landscape"),
    ("region deconv", "rigel.calibration.sweep", "chain_region_deconv"),
    ("boundary deconv", "rigel.calibration.sweep", "chain_boundary_deconv"),
    ("quant", "rigel.pipeline", "quant_from_buffer"),
    ("setup geometry + estimator", "rigel.pipeline", "_setup_geometry_and_estimator"),
    ("capture effective lengths", "rigel.calibration.capture_eff_length", "transcript_capture_eff_lengths"),
    ("score fragments", "rigel.pipeline", "_score_fragments"),
    ("build loci", "rigel.locus", "build_multi_loci"),
    ("assemble priors", "rigel.calibration.priors", "assemble_priors"),
    ("region locus shares", "rigel.calibration.priors", "_region_locus_shares"),
    ("boundary locus shares", "rigel.calibration.priors", "_boundary_locus_shares"),
    ("sum by locus", "rigel.calibration.priors", "_sum_by_locus"),
    ("global reference density", "rigel.calibration.capture_eff_length", "_global_reference_density"),
    ("partition", "rigel.locus_partition", "partition_and_free"),
    ("locus EM", "rigel.pipeline", "_run_locus_em_partitioned"),
    ("gDNA track", "rigel.calibration.track", "build_gdna_track"),
)

_RSS_TO_MB = 1.0 / (1024 * 1024) if sys.platform == "darwin" else 1.0 / 1024


# ── memory ─────────────────────────────────────────────────────────────────────────────────────────


class _ProcTaskInfo(ctypes.Structure):
    """macOS ``struct proc_taskinfo``. The trailing fields make ``sizeof`` match what the kernel expects;
    ``proc_pidinfo`` returns 0 for a buffer smaller than the flavour's own size."""

    _fields_ = [
        ("pti_virtual_size", ctypes.c_uint64),
        ("pti_resident_size", ctypes.c_uint64),
        ("pti_total_user", ctypes.c_uint64),
        ("pti_total_system", ctypes.c_uint64),
        ("pti_threads_user", ctypes.c_uint64),
        ("pti_threads_system", ctypes.c_uint64),
    ] + [(f"pti_i{i}", ctypes.c_int32) for i in range(12)]


_PROC_PIDTASKINFO = 4  # <sys/proc_info.h>
_LIBPROC = None
if sys.platform == "darwin":
    try:
        _LIBPROC = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        _LIBPROC.proc_pidinfo.restype = ctypes.c_int
    except OSError:  # pragma: no cover
        _LIBPROC = None


def current_rss_mb() -> float:
    """Current resident set size in MB, or ``nan`` where the platform offers no reader.

    Not the ``getrusage`` high-water mark, which never falls and so cannot attribute memory to a stage.
    macOS reads ``proc_pidinfo`` through ctypes (no fork: a ``ps`` subprocess would fork the very process
    being measured); Linux reads ``/proc/self/statm``.
    """
    if _LIBPROC is not None:
        info = _ProcTaskInfo()
        n = _LIBPROC.proc_pidinfo(os.getpid(), _PROC_PIDTASKINFO, 0, ctypes.byref(info), ctypes.sizeof(info))
        return info.pti_resident_size / 1024**2 if n == ctypes.sizeof(info) else float("nan")
    try:
        with open("/proc/self/statm") as handle:
            return int(handle.read().split()[1]) * resource.getpagesize() / 1024**2
    except (OSError, IndexError, ValueError):
        return float("nan")


def process_peak_mb() -> float:
    """The process's lifetime high-water mark in MB — exact for a single run, useless per stage."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_TO_MB


class MemoryTimeline:
    """Samples current RSS on a daemon thread between ``start()`` and ``stop()``; ``stop()`` joins the
    thread and takes one final sample so a short span's peak is never read below a later reading."""

    def __init__(self, interval_sec: float = 0.1) -> None:
        self.interval_sec = max(float(interval_sec), 0.001)
        self.samples: list[tuple[float, float]] = []
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _loop(self) -> None:
        while not self._stop.is_set():
            self.samples.append((time.perf_counter(), current_rss_mb()))
            self._stop.wait(self.interval_sec)

    def start(self) -> None:
        self.samples = [(time.perf_counter(), current_rss_mb())]
        self._stop.clear()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5.0)
            self._thread = None
        self.samples.append((time.perf_counter(), current_rss_mb()))

    def peak_between(self, t0: float, t1: float) -> float:
        vals = [mb for t, mb in self.samples if t0 <= t <= t1 and mb == mb]
        return max(vals) if vals else float("nan")


# ── the stage tree ─────────────────────────────────────────────────────────────────────────────────


class StageRecorder:
    """Records every probed call as ``(path, t0, t1, rss_in, rss_out)``, where ``path`` is the tuple of
    probe labels on the call stack. Calls on a thread other than the one that installed the recorder are
    recorded under that thread's own stack."""

    def __init__(self) -> None:
        self.calls: list[tuple[tuple[str, ...], float, float, float, float]] = []
        self._local = threading.local()

    def _stack(self) -> list[str]:
        if not hasattr(self._local, "stack"):
            self._local.stack = []
        return self._local.stack

    def wrap(self, label: str, fn):
        recorder = self

        def probed(*args, **kwargs):
            stack = recorder._stack()
            stack.append(label)
            path, t0, rss_in = tuple(stack), time.perf_counter(), current_rss_mb()
            try:
                return fn(*args, **kwargs)
            finally:
                recorder.calls.append((path, t0, time.perf_counter(), rss_in, current_rss_mb()))
                stack.pop()

        probed.__wrapped__ = fn
        probed.__name__ = getattr(fn, "__name__", label)
        probed.__qualname__ = getattr(fn, "__qualname__", label)
        probed.__doc__ = getattr(fn, "__doc__", None)
        return probed

    def summarise(self, timeline: MemoryTimeline | None, run_seconds: float) -> list[dict]:
        """One row per stage path: calls, inclusive seconds, self seconds, peak MB, held MB."""
        by_path: dict[tuple[str, ...], dict] = {}
        for path, t0, t1, rss_in, rss_out in self.calls:
            row = by_path.setdefault(path, dict(calls=0, total=0.0, windows=[], held=0.0))
            row["calls"] += 1
            row["total"] += t1 - t0
            row["windows"].append((t0, t1))
            if rss_in == rss_in and rss_out == rss_out:
                row["held"] += rss_out - rss_in
        rows = []
        for path, row in by_path.items():
            children = sum(r["total"] for p, r in by_path.items() if len(p) == len(path) + 1 and p[:-1] == path)
            peak = float("nan")
            if timeline is not None:
                peaks = [timeline.peak_between(a, b) for a, b in row["windows"]]
                peaks = [p for p in peaks if p == p]
                peak = max(peaks) if peaks else float("nan")
            rows.append(dict(
                path=list(path), label=path[-1], depth=len(path) - 1, calls=row["calls"],
                total=row["total"], self=row["total"] - children, peak_mb=peak, held_mb=row["held"],
                share=row["total"] / run_seconds if run_seconds > 0 else float("nan"),
            ))
        return _tree_order(rows)


def _tree_order(rows: list[dict]) -> list[dict]:
    """Depth-first, children after their parent, siblings by first-seen order of the parent's calls."""
    first_seen = {tuple(r["path"]): i for i, r in enumerate(rows)}
    return sorted(rows, key=lambda r: tuple(first_seen.get(tuple(r["path"][: k + 1]), 0) for k in range(len(r["path"]))))


def _resolve(module_name: str, attr_path: str):
    """``(owner, attribute name, original object, static attribute)`` for ``module.attr_path``."""
    module = importlib.import_module(module_name)
    owner = module
    *parents, name = attr_path.split(".")
    for part in parents:
        owner = getattr(owner, part)
    return owner, name, getattr(owner, name), inspect.getattr_static(owner, name)


def install_probes(recorder: StageRecorder, probes=PROBES) -> tuple[list, list[str]]:
    """Patch every probe target in place; return ``(undo list, missing probe descriptions)``.

    A module-level function is replaced in every loaded ``rigel`` module that holds it, which covers
    ``from x import f`` bindings; a class/static method is replaced on its class.
    """
    undo: list = []
    missing: list[str] = []
    for label, module_name, attr_path in probes:
        try:
            owner, name, original, static = _resolve(module_name, attr_path)
        except (ImportError, AttributeError) as exc:
            missing.append(f"{label}: {module_name}.{attr_path} ({type(exc).__name__})")
            continue
        if isinstance(owner, type):
            if isinstance(static, classmethod):
                replacement = classmethod(recorder.wrap(label, static.__func__))
            elif isinstance(static, staticmethod):
                replacement = staticmethod(recorder.wrap(label, static.__func__))
            else:
                replacement = recorder.wrap(label, static)
            undo.append((owner, name, static))
            setattr(owner, name, replacement)
            continue
        wrapped = recorder.wrap(label, original)
        for mod in list(sys.modules.values()):
            if not isinstance(mod, types.ModuleType):
                continue
            mod_name = getattr(mod, "__name__", "") or ""
            if not (mod_name == module_name or mod_name.startswith("rigel") or mod_name == "__main__"):
                continue
            for attr, value in list(vars(mod).items()):
                if value is original:
                    undo.append((mod, attr, original))
                    setattr(mod, attr, wrapped)
    return undo, missing


def remove_probes(undo: list) -> None:
    for owner, name, original in reversed(undo):
        setattr(owner, name, original)


# ── one run ────────────────────────────────────────────────────────────────────────────────────────


def _git_sha() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                             cwd=Path(__file__).resolve().parents[2], check=True)
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                               text=True, cwd=Path(__file__).resolve().parents[2]).stdout.strip()
        return out.stdout.strip() + ("+dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def profile_run(bam: str, index_dir: str, *, threads: int | None, knobs: list[str], label: str,
                cprofile_path: str | None, scan_only: bool = False) -> dict:
    """Load the index and run the pipeline once under the probes; return the report dict."""
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    cfg = PipelineConfig()
    if threads is not None:
        cfg = dataclasses.replace(cfg, scan=dataclasses.replace(cfg.scan, total_threads=threads),
                                  em=dataclasses.replace(cfg.em, n_threads=threads),
                                  calibration=dataclasses.replace(cfg.calibration, n_threads=threads))
    for knob in knobs:
        cfg = set_field(cfg, knob)

    import rigel.pipeline as pipeline  # noqa: F401 — imported so its bindings exist to be patched

    recorder = StageRecorder()
    undo, missing = install_probes(recorder)
    timeline = MemoryTimeline()
    profiler = cProfile.Profile() if cprofile_path else None
    timeline.start()
    t0 = time.perf_counter()
    try:
        if profiler is not None:
            profiler.enable()
        index = TranscriptIndex.load(index_dir)
        if scan_only:
            scan_cfg = cfg.scan
            if scan_cfg.sj_strand_tag == "auto":
                scan_cfg = dataclasses.replace(scan_cfg, sj_strand_tag=pipeline._native_detect_sj_tag(bam))
            stats = pipeline.scan_and_buffer(bam, index, scan_cfg)[0]
        else:
            stats = pipeline.run_pipeline(bam, index, cfg).stats
    finally:
        if profiler is not None:
            profiler.disable()
        wall = time.perf_counter() - t0
        timeline.stop()
        remove_probes(undo)
    if profiler is not None:
        profiler.dump_stats(cprofile_path)

    n_frag = int(getattr(stats, "n_fragments", 0) or 0)
    return dict(
        label=label, bam=str(bam), index=str(index_dir), git=_git_sha(), host=platform.node(),
        machine=platform.machine(), python=platform.python_version(), threads=threads, knobs=knobs,
        scan_only=scan_only, fragments=n_frag, wall_seconds=wall, process_peak_mb=process_peak_mb(),
        sampled_peak_mb=max((mb for _, mb in timeline.samples if mb == mb), default=float("nan")),
        missing_probes=missing, stages=recorder.summarise(timeline, wall),
    )


# ── reporting ──────────────────────────────────────────────────────────────────────────────────────


def format_report(rep: dict) -> str:
    lines = []
    rate = rep["fragments"] / rep["wall_seconds"] if rep["wall_seconds"] > 0 else float("nan")
    lines.append(f"  RIGEL PROFILE — {rep['label']}   git {rep['git']}   threads {rep['threads']}")
    lines.append(f"  {rep['fragments']:,} fragments   {rep['wall_seconds']:.1f} s   {rate:,.0f} frag/s   "
                 f"peak {rep['process_peak_mb']:,.0f} MB" + (f"   knobs {rep['knobs']}" if rep["knobs"] else ""))
    if rep["missing_probes"]:
        lines.append("  ⚠ probes whose target moved (not measured): " + "; ".join(rep["missing_probes"]))
    head = f"  {'stage':<40s} {'calls':>7s} {'total s':>9s} {'self s':>9s} {'% run':>7s} {'peak MB':>9s} {'held MB':>9s}"
    lines += ["", head, "  " + "-" * (len(head) - 2)]
    for r in rep["stages"]:
        name = "  " * r["depth"] + r["label"]
        lines.append(f"  {name:<40s} {r['calls']:>7d} {r['total']:>9.2f} {r['self']:>9.2f} "
                     f"{100 * r['share']:>6.1f}% {r['peak_mb']:>9,.0f} {r['held_mb']:>+9,.0f}")
    return "\n".join(lines)


def format_compare(a: dict, b: dict) -> str:
    lines = [f"  COMPARE  A = {a['label']} ({a['git']})   B = {b['label']} ({b['git']})",
             f"  wall {a['wall_seconds']:.1f} s -> {b['wall_seconds']:.1f} s   "
             f"peak {a['process_peak_mb']:,.0f} -> {b['process_peak_mb']:,.0f} MB", ""]
    head = f"  {'stage':<40s} {'A total':>9s} {'B total':>9s} {'B/A':>7s} {'A peak':>9s} {'B peak':>9s}"
    lines += [head, "  " + "-" * (len(head) - 2)]
    rows_b = {tuple(r["path"]): r for r in b["stages"]}
    seen = set()
    for r in a["stages"]:
        path = tuple(r["path"])
        seen.add(path)
        s = rows_b.get(path)
        name = "  " * r["depth"] + r["label"]
        if s is None:
            lines.append(f"  {name:<40s} {r['total']:>9.2f} {'—':>9s}")
            continue
        ratio = s["total"] / r["total"] if r["total"] > 0 else float("nan")
        lines.append(f"  {name:<40s} {r['total']:>9.2f} {s['total']:>9.2f} {ratio:>7.2f} "
                     f"{r['peak_mb']:>9,.0f} {s['peak_mb']:>9,.0f}")
    for path, s in rows_b.items():
        if path not in seen:
            lines.append(f"  {'  ' * s['depth'] + s['label']:<40s} {'—':>9s} {s['total']:>9.2f}   (only in B)")
    return "\n".join(lines)


def cprofile_top(path: str, top: int = 40) -> str:
    buf = io.StringIO()
    pstats.Stats(path, stream=buf).sort_stats("tottime").print_stats(top)
    return buf.getvalue()


# ── self-test: the probe machinery falsified without a BAM ─────────────────────────────────────────


def self_test() -> int:
    import numpy as np

    mod = types.ModuleType("rigel_profiler_selftest")
    sys.modules[mod.__name__] = mod

    def inner(n):
        time.sleep(0.05)
        return np.ones(n)

    def outer():
        time.sleep(0.05)
        a = mod.inner(25_000_000)  # 200 MB of float64
        b = mod.inner(10)
        time.sleep(0.2)
        return float(a.sum() + b.sum())

    class Holder:
        @classmethod
        def build(cls, x):
            return x + 1

    mod.inner, mod.outer, mod.Holder = inner, outer, Holder
    probes = (("outer", mod.__name__, "outer"), ("inner", mod.__name__, "inner"),
              ("classmethod", mod.__name__, "Holder.build"), ("ghost", mod.__name__, "not_there"))
    checks = []
    recorder = StageRecorder()
    undo, missing = install_probes(recorder, probes)
    timeline = MemoryTimeline(interval_sec=0.01)
    timeline.start()
    t0 = time.perf_counter()
    try:
        mod.outer()
        holder_result = mod.Holder.build(1)
    finally:
        wall = time.perf_counter() - t0
        timeline.stop()
        remove_probes(undo)
    rows = {tuple(r["path"]): r for r in recorder.summarise(timeline, wall)}

    checks.append(("a missing target is reported, not dropped", missing and "ghost" in missing[0]))
    checks.append(("nesting follows the call stack", ("outer", "inner") in rows))
    checks.append(("calls are counted", rows.get(("outer", "inner"), {}).get("calls") == 2))
    o = rows.get(("outer",), {})
    checks.append(("self = inclusive minus children", abs(o.get("self", -1) - (o.get("total", 0) - rows.get(("outer", "inner"), {}).get("total", 0))) < 1e-9))
    checks.append(("inclusive covers the sleeps", o.get("total", 0) >= 0.35))
    checks.append(("a 200 MB allocation shows in the stage's peak",
                   o.get("peak_mb", float("nan")) - timeline.samples[0][1] > 150))
    checks.append(("a classmethod probe records and still returns", ("classmethod",) in rows and holder_result == 2))
    checks.append(("probes are removed afterwards", mod.inner is inner and mod.outer is outer
                   and not hasattr(inspect.getattr_static(mod.Holder, "build").__func__, "__wrapped__")))
    rep = dict(label="t", git="g", threads=1, knobs=[], fragments=0, wall_seconds=wall, process_peak_mb=0.0,
               missing_probes=missing, stages=recorder.summarise(timeline, wall))
    checks.append(("the report renders every stage", all(r["label"] in format_report(rep) for r in rep["stages"])))
    checks.append(("compare joins identical reports at ratio 1.00", "1.00" in format_compare(rep, rep)))

    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}s}")
    passed = sum(bool(ok) for _, ok in checks)
    print(f"\n{passed}/{len(checks)} self-test checks pass")
    return 0 if passed == len(checks) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", help="name-sorted BAM with NH tags")
    ap.add_argument("--index", help="rigel index directory")
    ap.add_argument("--label", default=None, help="report label (default: the library directory)")
    ap.add_argument("--threads", type=int, default=None, help="the scan's, the EM's and calibration's thread budget")
    ap.add_argument("--set", dest="knobs", action="append", default=[], metavar="SECTION.FIELD=VALUE",
                    help="override one PipelineConfig field; repeatable")
    ap.add_argument("--cprofile", default=None, metavar="OUT.prof", help="also write a cProfile dump")
    ap.add_argument("--scan-only", action="store_true", help="run the scan alone, not the pipeline")
    ap.add_argument("--out", default=None, help="write the JSON report here")
    ap.add_argument("--compare", nargs=2, metavar=("A.json", "B.json"), help="compare two reports")
    ap.add_argument("--self-test", action="store_true", help="falsify the probe machinery; no BAM")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.compare:
        a, b = (json.loads(Path(p).read_text()) for p in args.compare)
        print(format_compare(a, b))
        return 0
    if not (args.bam and args.index):
        ap.error("--bam and --index are required (or --compare / --self-test)")
    label = args.label or Path(args.bam).resolve().parent.parent.name
    rep = profile_run(args.bam, args.index, threads=args.threads, knobs=args.knobs, label=label,
                      cprofile_path=args.cprofile, scan_only=args.scan_only)
    print(format_report(rep))
    if args.cprofile:
        print("\n  cProfile — top 40 by self time:\n")
        print(cprofile_top(args.cprofile))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(rep, indent=2))
        print(f"\n  report -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
