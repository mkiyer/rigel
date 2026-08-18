#!/usr/bin/env python
"""WHERE DOES THE TIME AND THE MEMORY GO? — the whole pipeline, per phase, on one library.

Profiles ``run_pipeline(bam, index)`` and attributes both **wall clock** and **peak RSS** to each of
the three pipeline phases:

  1. **scan**      — ``scan_and_buffer`` (C++ htslib BAM scan, fragment resolution, model
                      training, the fractional accumulator).
  2. **calibrate** — the gDNA/RNA deconvolution + library hyperparameter fit. ⭐ This is the phase
                      the memory-bounded chunked grid solve lives in, so it is the one to watch.
  3. **quant**     — ``quant_from_buffer`` (fragment scoring, routing, per-locus EM).

It works by **wrapping those three functions and running the real `run_pipeline`** — it does not
re-implement the pipeline's internals, so it always profiles the production flow and cannot drift
from the actual phase wiring.

⛔⛔ **PROFILE A HIGH-DEPTH REAL RNA-seq LIBRARY, NOT cfRNA AND NOT A PANEL CONDITION** (owner,
2026-08-17, reversing the older *"profile on real cfRNA"* rule). The cfRNA on disk is sparse and
small; it under-represents exactly the regime the performance work targets, and a toy ranks hotspots
backwards outright (`TRAPS: toys-rank-hotspots-backwards`). `docs/TESTING.md` §7 is the ruling.

Usage::

    python scripts/profiling/profiler.py --bam sample.bam --index index/ \\
        [--threads N] [--repeat 3] [--cprofile] [--out profile.json]

    # ⭐ COMPUTE AGAINST MEMORY: put any PipelineConfig field on the x-axis.
    python scripts/profiling/profiler.py --bam s.bam --index index/ \\
        --sweep calibration.sweep_n_grid=30,60,120

    # native + OpenMP frames, which cProfile cannot see:
    py-spy record -o pipeline.svg -- python scripts/profiling/profiler.py --bam s.bam --index index/

    # the instrument's own falsification, no BAM and no index:
    python scripts/profiling/profiler.py --self-test

Inside the activated `rigel` conda env. The BAM must be name-sorted with NH tags (as `rigel quant`
requires).

⭐⭐ **WHAT `--sweep` IS FOR, AND IT ALREADY HAS A KNOB WORTH TURNING.** Measured on 10 M fragments at
4 threads, `calibration.sweep_n_grid` **30 → 120**: the calibrate phase goes **5.2 s → 16.2 s (3.1×)**
and process peak RSS **7,788 → 8,925 MB**. That is the compute-vs-memory curve the grid solve has to be
argued on, and nothing here could draw it before 2026-08-17.

⭐⭐ **THIS MODULE IS ALSO THE RSS TOOLKIT FOR `scan_profile.py`** — :func:`_snap_rss_current` and
:class:`MemoryTimeline`. ⚠ They are here rather than in `scan_profile.py` because the CURRENT-vs-PEAK
distinction below is one fact and belongs in one place.

⭐ Gated by `tests/test_scripts_index.py`, which reaches this tree since 2026-08-17. Before that no
gate did, and the directory rotted **twice** with the one defect class that gate already catches
(`TRAPS: a-green-suite-hid-five-dead-instruments`).
"""

from __future__ import annotations

import argparse
import cProfile
import csv
import ctypes
import gc
import json
import logging
import os
import platform
import pstats
import resource
import statistics
import sys
import threading
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
    """Process peak resident set size (high-water mark) in MB.

    ⛔ **CUMULATIVE, and that is the whole reason :class:`MemoryTimeline` exists.** ``ru_maxrss`` is a
    high-water mark for the LIFETIME of the process, so in a sweep it can only rise: the second
    configuration's number contains the first one's peak and cannot fall below it. It is exact and it is
    the right number for a single run; it cannot answer "did this run give the memory back".
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_TO_MB


class _ProcTaskInfo(ctypes.Structure):
    """macOS ``struct proc_taskinfo`` (``<sys/proc_info.h>``) — only the first two fields are read.

    ⚠ The trailing fields are declared so ``sizeof`` matches what the kernel expects: ``proc_pidinfo``
    returns 0 when handed a buffer smaller than the flavour's own size, so a truncated struct silently
    reads as "unavailable" rather than as a short read.
    """

    _fields_ = [
        ("pti_virtual_size", ctypes.c_uint64),
        ("pti_resident_size", ctypes.c_uint64),
        ("pti_total_user", ctypes.c_uint64),
        ("pti_total_system", ctypes.c_uint64),
        ("pti_threads_user", ctypes.c_uint64),
        ("pti_threads_system", ctypes.c_uint64),
    ] + [(f"pti_i{i}", ctypes.c_int32) for i in range(12)]


#: ``PROC_PIDTASKINFO`` from ``<sys/proc_info.h>`` — an OS ABI constant, not a tunable.
_PROC_PIDTASKINFO = 4

_LIBPROC = None
if sys.platform == "darwin":
    try:
        _LIBPROC = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        _LIBPROC.proc_pidinfo.restype = ctypes.c_int
    except OSError:  # pragma: no cover — a macOS without libproc is not a supported platform
        _LIBPROC = None


def _snap_rss_current() -> float:
    """**CURRENT** resident set size in MB, or ``nan`` where the platform offers no reader.

    ⛔⛔ **NOT ``_peak_rss_mb``, AND SUBSTITUTING IT WOULD BE A LIE.** ``resource`` exposes only the
    high-water mark, which by construction never falls — so a "before / after / after-release" triple
    built on it reports memory as never released no matter what the code does, and a per-PHASE peak
    built on it reports every phase's peak as at least the previous phase's. Attributing memory to a
    phase is the whole point of the memory work, so it needs a reading that CAN go down.

    ⭐ One reader per platform, no third-party dependency (`psutil` is not in this environment):

    * Linux — ``/proc/self/statm`` field 2 (resident pages) × the page size;
    * macOS — ``proc_pidinfo(PROC_PIDTASKINFO).pti_resident_size``, through ``ctypes``.

    ⛔⛔ **THE macOS READER WAS ``ps -o rss=`` IN A SUBPROCESS UNTIL 2026-08-17, AND THAT FORKS THE
    PROCESS BEING PROFILED.** ``subprocess`` cannot take the ``posix_spawn`` path while ``close_fds``
    is true (its default), so every sample was a real ``fork()`` — whose cost scales with the parent's
    page-table size, which is precisely the quantity under measurement. Measured on this machine at
    **1.49 GB** resident: ``ps`` **6,985.7 µs**/sample against ``libproc`` **2.4 µs**/sample, a factor
    of **2,925**. At the shipped 250 ms interval that is 2.8 % of wall clock spent forking at 1.5 GB
    and worse as the library grows — an instrument that perturbs its own subject. ⭐ The two agree:
    15.77 MB vs 15.86 MB read back to back, so this is a cost repair and not a definition change.

    ⚠ Returns ``nan`` rather than ``0.0`` when the platform has no reader. A zero would be
    indistinguishable from a real reading of an empty process and would silently poison every RSS
    column; callers report ``available`` from it (:class:`MemoryTimeline`).

    ⭐ A standing finding, and it belongs to the ALLOCATOR rather than to this function: a freed 200 MB
    numpy array moved this reading by **0.0 MB** on macOS 25.6, and ``buffer.release()`` did not move
    it either. ⛔ That is not evidence the two readers are interchangeable — ``ru_maxrss`` is *defined*
    as non-decreasing, so it could not have shown a release had one happened; this one is unconstrained
    and showed none. Reproduce with ``scan_profile.py --repeat 2``, comparing ``peak_rss_mb`` against
    ``process_peak_rss_mb``.
    """
    if _LIBPROC is not None:
        info = _ProcTaskInfo()
        written = _LIBPROC.proc_pidinfo(
            os.getpid(), _PROC_PIDTASKINFO, 0, ctypes.byref(info), ctypes.sizeof(info)
        )
        if written == ctypes.sizeof(info):
            return info.pti_resident_size / 1024**2
        return float("nan")
    try:
        with open("/proc/self/statm") as handle:
            resident_pages = int(handle.read().split()[1])
        return resident_pages * resource.getpagesize() / 1024**2
    except (OSError, IndexError, ValueError):
        return float("nan")


class MemoryTimeline:
    """Sample **current** RSS on a background thread for the span between ``start()`` and ``stop()``.

    ⭐ **Per-SPAN, which is what ``_peak_rss_mb`` cannot be.** One timeline wraps one run, so a sweep's
    second configuration is measured on its own, and a release WOULD be visible as a fall (⚠ none was
    observed on the workload measured — see :func:`_snap_rss_current`; the allocator kept the pages).

    ⛔ **An unavailable reader is REPORTED, never faked.** If :func:`_snap_rss_current` cannot read the
    platform, ``available`` is ``False``, ``peak_mb`` is ``nan`` and the CSV is written with its one
    header row — an instrument must not print an inert column as a measurement
    (`TRAPS: an-ablation-that-never-ran`, the same failure in the reporting direction).

    ⚠ The sampler is a ``daemon`` thread and ``stop()`` JOINS it, so a run that raises cannot leave a
    thread sampling into the next configuration's numbers.
    """

    def __init__(self, interval_sec: float = 0.25) -> None:
        self.interval_sec = max(float(interval_sec), 0.001)
        self.samples: list[tuple[float, float]] = []
        self.available = False
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._t0 = 0.0

    @property
    def peak_mb(self) -> float:
        """The largest CURRENT-RSS sample in the span, or ``nan`` if nothing was sampled."""
        values = [mb for _, mb in self.samples if mb == mb]  # `mb == mb` drops nan
        return max(values) if values else float("nan")

    def _loop(self) -> None:
        while not self._stop.is_set():
            self.samples.append((time.perf_counter() - self._t0, _snap_rss_current()))
            self._stop.wait(self.interval_sec)

    def start(self) -> None:
        # ⚠ ONE reading, used for both the availability probe and the span's first sample — on macOS
        # each call is a `ps` subprocess, and a second start() must not append to a finished span.
        first = _snap_rss_current()
        self.available = first == first  # False iff nan
        self.samples = [(0.0, first)]
        self._t0 = time.perf_counter()
        self._stop.clear()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        """End the span, joining the sampler, and take one FINAL sample at the boundary.

        ⛔ **The final sample is not cosmetic.** Without it the span's last ``interval_sec`` is unsampled,
        so on a short run ``peak_mb`` can read BELOW the RSS measured immediately afterwards — a "peak"
        smaller than a later reading, which a reader is right to distrust. Measured on a 0.46 s scan at a
        100 ms interval: sampled peak **246.66 MB** against **252.89 MB** read one line later.
        """
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5.0)
            self._thread = None
        self.samples.append((time.perf_counter() - self._t0, _snap_rss_current()))

    def write_csv(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["seconds", "rss_mb"])
            for seconds, mb in self.samples:
                writer.writerow([f"{seconds:.4f}", f"{mb:.3f}"])


@contextmanager
def _timed_stages(record: dict[str, dict], *, memory_interval_sec: float):
    """Wrap each stage function with a timer + a per-stage memory timeline; restore on exit.

    Records ``{label: {seconds, peak_rss_mb, rss_entry_mb, rss_exit_mb, process_peak_rss_mb}}`` for
    each stage as the real ``run_pipeline`` calls it. Wrappers forward ``*args, **kwargs`` so they are
    immune to stage-signature changes.

    ⛔⛔ **``peak_rss_mb`` IS SAMPLED PER SPAN, AND THE ``ru_maxrss`` IT REPLACED COULD NOT ATTRIBUTE
    MEMORY TO A PHASE AT ALL.** ``ru_maxrss`` is a process-LIFETIME high-water mark: it never falls, so
    the calibrate row inherited scan's peak and the quant row inherited both, and the three numbers were
    a running maximum wearing three labels. To price compute against memory you need the peak *of a
    phase*, which needs a reading that can go down — hence :class:`MemoryTimeline` here.

    ⚠ ``process_peak_rss_mb`` is kept beside it, honestly named, because it is EXACT where the sampled
    peak is a sample: a spike shorter than ``memory_interval_sec`` is invisible to the sampler and the
    high-water mark still catches it. Read them together — a large gap means the phase spiked between
    samples. ⚠ ``rss_entry_mb`` / ``rss_exit_mb`` bracket the span, so memory a phase HANDS ON (the
    fragment buffer, the payload) is separable from memory it merely borrowed.

    ⛔ **ON A PIPELINE THAT ONLY GROWS THE TWO PEAK COLUMNS COINCIDE, AND THAT IS NOT A REDUNDANCY.**
    Measured on the smoke run: 3,254 / 3,802 / 7,912 from both, identically, because each phase's own
    peak *was* the running maximum. They separate the moment a phase releases — which is what the
    memory work is trying to make happen, so the column has to be able to show it before the work
    starts, not after (`TRAPS: an-ablation-that-never-ran`, in its reporting form).
    """
    originals = [(mod, attr, getattr(mod, attr)) for mod, attr, _ in _STAGES]

    def make_timed(orig, label):
        def timed(*args, **kwargs):
            timeline = MemoryTimeline(interval_sec=memory_interval_sec)
            entry = _snap_rss_current()
            timeline.start()
            start = time.perf_counter()
            try:
                return orig(*args, **kwargs)
            finally:
                seconds = time.perf_counter() - start
                timeline.stop()
                record[label] = {
                    "seconds": seconds,
                    "peak_rss_mb": timeline.peak_mb,
                    "rss_entry_mb": entry,
                    "rss_exit_mb": _snap_rss_current(),
                    "process_peak_rss_mb": _peak_rss_mb(),
                    "memory_samples": len(timeline.samples),
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
                 *, cprofile: bool = False, memory_interval_sec: float = 0.25) -> dict:
    """Run the pipeline once; return a result dict with per-stage timings, per-stage memory + totals."""
    gc.collect()
    stages: dict[str, dict] = {}
    profiler = cProfile.Profile() if cprofile else None
    rss_before = _snap_rss_current()
    start = time.perf_counter()
    if profiler:
        profiler.enable()
    with _timed_stages(stages, memory_interval_sec=memory_interval_sec):
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
        "process_peak_rss_mb": _peak_rss_mb(),
        "rss_before_mb": rss_before,
        "rss_after_mb": _snap_rss_current(),
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
        "process_peak_rss_mb": ms([r["process_peak_rss_mb"] for r in runs]),
        "fragments_per_sec": ms([r["fragments_per_sec"] for r in runs]),
        "n_fragments": runs[-1]["n_fragments"],
        "stages": {},
    }
    for lbl in labels:
        rows = [r["stages"].get(lbl, {}) for r in runs]
        agg["stages"][lbl] = {
            # ⭐ `peak_rss_mb` is SAMPLED PER SPAN, so it is attributable to this phase; the old
            # `ru_maxrss` column could only be a running maximum. See `_timed_stages`.
            name: ms([row.get(name, float("nan")) for row in rows])
            for name in ("seconds", "peak_rss_mb", "rss_entry_mb", "rss_exit_mb")
        }
    return agg


def _bar(frac: float, width: int = 20) -> str:
    fill = int(round(max(0.0, min(1.0, frac)) * width))
    return "█" * fill + "·" * (width - fill)


def format_report(agg: dict, label: str) -> str:
    """⭐ TIME AND MEMORY SIDE BY SIDE, because the choice the tuning has to make is between them.

    ⛔⛔ **READ `peak` AS "THE PROCESS HIGH-WATER MARK WHILE THIS PHASE RAN", NOT AS "WHAT THIS PHASE
    NEEDS".** It contains everything its predecessors handed it — on the smoke run below `calibrate`
    reads 3,802 MB while *entering* at 3,285 MB, so its own contribution is 517 MB and not 3.8 GB.
    ⭐ `held` is the separating column: `rss_exit − rss_entry`, the memory a phase hands ON rather than
    merely borrows. **A large `peak − exit` gap with a small `held` is a chunking candidate** — a
    transient the phase could be made to take in pieces; a large `held` is a lifetime problem and
    chunking will not touch it.

    ⭐ Measured once end to end, 10 M fragments at 4 threads (⚠ a PANEL condition, which is the wrong
    substrate for a real measurement — `docs/TESTING.md` §7 — and is used here only to show the columns
    are not inert): scan 38.5 s / peak 3,254 / held **+3,072** (the fragment buffer, a lifetime cost);
    calibrate 6.8 s / peak 3,802 / held +517; quant 51.0 s / peak 7,912 / exit 6,018 — ⭐ **a 1.9 GB
    transient the process peak is set by and nothing keeps.** That gap is the number a memory bound has
    to be argued against, and no column in this report before 2026-08-17 could express it.
    """
    total = agg["total_seconds"]["mean"]
    boundaries = [
        "",
        "=" * 92,
        f"  RIGEL PROFILE — {label}   ({agg['n_runs']} run(s))",
        "=" * 92,
        f"  fragments: {agg['n_fragments']:,}    "
        f"throughput: {agg['fragments_per_sec']['mean']:,.0f} frag/s    "
        f"process peak RSS: {agg['process_peak_rss_mb']['mean']:,.0f} MB",
        "",
        f"  {'phase':<12} {'seconds':>10} {'± stdev':>9} {'% total':>8}  {'peak MB':>10} {'held MB':>10}"
        f"  {'':<20}",
        "  " + "-" * 88,
    ]
    for _, _, lbl in _STAGES:
        st = agg["stages"][lbl]
        sec, sd = st["seconds"]["mean"], st["seconds"]["stdev"]
        pct = (sec / total) if total > 0 else 0.0
        peak = st["peak_rss_mb"]["mean"]
        held = st["rss_exit_mb"]["mean"] - st["rss_entry_mb"]["mean"]
        boundaries.append(
            f"  {lbl:<12} {sec:>10.3f} {sd:>9.3f} {100*pct:>7.1f}%  {peak:>10,.0f} {held:>+10,.0f}"
            f"  {_bar(pct, 20)}"
        )
    other = agg["other_seconds"]["mean"]
    opct = (other / total) if total > 0 else 0.0
    boundaries.append(
        f"  {'other':<12} {other:>10.3f} {'':>9} {100*opct:>7.1f}%  {'':>10} {'':>10}  {_bar(opct, 20)}"
    )
    boundaries.append("  " + "-" * 88)
    boundaries.append(
        f"  {'TOTAL':<12} {total:>10.3f} {agg['total_seconds']['stdev']:>9.3f} {100.0:>7.1f}%"
    )
    boundaries.append("=" * 92)
    return "\n".join(boundaries)


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


def _set_knob(cfg: PipelineConfig, dotted: str, raw: str) -> PipelineConfig:
    """Return ``cfg`` with one ``<section>.<field>`` replaced, the value parsed from the field's TYPE.

    ⭐ **THE SWEEP KNOB IS THE CONFIG ITSELF, so this needs no per-knob code and gains new knobs for
    free.** Pricing compute against memory means putting a dial on the x-axis and wall clock + peak RSS
    on the y — and every dial that matters is already a `PipelineConfig` field
    (`calibration.sweep_n_grid`, `calibration.n_grid_tilt`, `scan.fragments_per_chunk`,
    `scan.buffer_size_bytes`, `em.n_threads`). A knob added later is sweepable the day it lands.

    ⛔ The field must EXIST and its current value must be an `int`, `float` or `bool` — an unknown name
    is an error rather than a silent no-op, which is the failure `TRAPS: an-ablation-that-never-ran`
    names: a sweep whose knob does nothing prints a flat line that reads exactly like "this dial does
    not matter".
    """
    section, _, field = dotted.partition(".")
    if not field or not hasattr(cfg, section):
        raise SystemExit(f"--sweep needs <section>.<field>; `{dotted}` names no section of PipelineConfig")
    sub = getattr(cfg, section)
    if not hasattr(sub, field):
        raise SystemExit(f"{type(sub).__name__} has no field `{field}`")
    current = getattr(sub, field)
    if isinstance(current, bool):
        value = raw.lower() in ("1", "true", "yes")
    elif isinstance(current, int):
        value = int(raw)
    elif isinstance(current, float):
        value = float(raw)
    else:
        raise SystemExit(f"{dotted} is {type(current).__name__}; only int/float/bool are sweepable")
    return replace(cfg, **{section: replace(sub, **{field: value})})


def _cprofile_text(prof: cProfile.Profile, top: int) -> str:
    buf = StringIO()
    pstats.Stats(prof, stream=buf).sort_stats("cumulative").print_stats(top)
    return buf.getvalue()


def _self_test() -> int:
    """⭐ FALSIFY THE MEMORY INSTRUMENT WITH NO BAM AND NO INDEX — a 2 s check that it can SEE memory.

    ⛔ The failure this exists for is `TRAPS: an-ablation-that-never-ran` in its reporting form: an RSS
    column that is inert reads exactly like a phase that used no memory. Three gates, each perturbed:

      ① the reader is AVAILABLE and returns a finite, positive MB;
      ② allocating and TOUCHING a known block MOVES the reading by at least most of that block —
        this is what a broken reader (`nan`, a constant, a fork that reads the child) fails;
      ③ a :class:`MemoryTimeline` spanning that allocation reports a `peak_mb` that CONTAINS it, and
        does so from its background thread — ⛔ `stop()` takes a final sample, so a peak alone would
        pass with a dead sampler; the gate therefore also requires samples the LOOP produced.

    ⛔⛔ **ONE HOLE, PINNED AS A HOLE RATHER THAN LEFT TO READ AS COVERAGE.** Substituting
    :func:`_peak_rss_mb` for the reader passes **3/3** — measured, in a fresh process. It has to: this
    workload only ever grows, and on a monotone span a high-water mark and a current reading are the
    same number. The substitution is wrong for the reason the reader's own docstring gives (a per-phase
    peak built on `ru_maxrss` is a running maximum), and no gate here can see it, because the allocator
    on this platform does not hand freed pages back — so the falling reading that WOULD discriminate
    them is not observable. ⭐ Perturbations that DO fire, each measured in a fresh process: a `nan`
    reader 0/3, a constant reader 1/3, a sampler whose loop never runs 2/3 (③ alone, at exactly the
    2 samples `start()` and `stop()` write).
    """
    block_mb = 256
    checks: list[tuple[str, bool, str]] = []

    before = _snap_rss_current()
    checks.append(("① reader available and finite", before == before and before > 0.0, f"{before:.1f} MB"))

    timeline = MemoryTimeline(interval_sec=0.01)
    timeline.start()
    import numpy as np

    block = np.ones(block_mb * 1024**2 // 8, dtype=np.float64)  # touched by `ones`, so it is resident
    block[0] = 1.0
    time.sleep(0.1)
    after = _snap_rss_current()
    timeline.stop()
    grew = after - before
    # ⚠ 0.75 of the block, not the whole of it: the allocator may already have had free pages mapped.
    checks.append(("② a touched 256 MB block moves the reading", grew > 0.75 * block_mb,
                   f"+{grew:.1f} MB of {block_mb} MB"))
    # ⛔ `start()` writes one sample and `stop()` writes one more, so `> 2` is what proves the LOOP ran.
    sampled_the_span = timeline.peak_mb >= before + 0.75 * block_mb and len(timeline.samples) > 2
    checks.append(("③ its background sampler saw the span", sampled_the_span,
                   f"peak {timeline.peak_mb:.1f} MB over {len(timeline.samples)} samples"))
    del block

    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}   {detail}")
    n_ok = sum(ok for _, ok, _ in checks)
    print(f"\n  {n_ok}/{len(checks)} — memory instrumentation on {platform.system()}")
    return 0 if n_ok == len(checks) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description="Profile the Rigel quantification pipeline.")
    ap.add_argument("--bam", help="Name-sorted BAM with NH tags.")
    ap.add_argument("--index", help="Rigel index directory.")
    ap.add_argument("--label", default=None, help="Run label for the report (default: BAM name).")
    ap.add_argument("--threads", type=int, default=None, help="Thread budget (scan + EM).")
    ap.add_argument("--repeat", type=int, default=1, help="Repeat N times and average.")
    ap.add_argument("--cprofile", action="store_true",
                    help="Also collect cProfile (Python hotspots) on the first run.")
    ap.add_argument("--cprofile-top", type=int, default=30, help="Top-N cProfile rows to print.")
    ap.add_argument("--memory-interval-ms", type=int, default=250,
                    help="Per-phase RSS sampling interval. A spike shorter than this is only visible "
                         "in the process high-water mark.")
    ap.add_argument("--sweep", default=None, metavar="SECTION.FIELD=V1,V2,…",
                    help="Re-run once per value of one PipelineConfig field and print time against "
                         "peak RSS. e.g. `--sweep calibration.sweep_n_grid=30,60,120`.")
    ap.add_argument("--out", default=None, help="Write the JSON report to this path.")
    ap.add_argument("--self-test", action="store_true",
                    help="Falsify the memory instrumentation with no BAM and no index; then exit.")
    args = ap.parse_args()

    if args.self_test:
        return _self_test()
    if not args.bam or not args.index:
        ap.error("--bam and --index are required (or pass --self-test)")

    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    bam_path = str(Path(args.bam).resolve())
    label = args.label or Path(args.bam).name
    base_config = _build_config(args)

    knob, values = None, [None]
    if args.sweep:
        knob, _, joined = args.sweep.partition("=")
        values = [v for v in joined.split(",") if v]
        if not values:
            ap.error(f"--sweep {args.sweep} lists no values")

    index = TranscriptIndex.load(Path(args.index))
    memory_interval_sec = max(args.memory_interval_ms, 1) / 1000.0
    points: list[dict] = []
    cprofile_text = None
    for value in values:
        config = base_config if value is None else _set_knob(base_config, knob, value)
        point_label = label if value is None else f"{label}  [{knob}={value}]"
        runs: list[dict] = []
        for i in range(max(1, args.repeat)):
            res = profile_once(bam_path, index, config, cprofile=(args.cprofile and i == 0),
                               memory_interval_sec=memory_interval_sec)
            if res.get("_cprofile") is not None:
                cprofile_text = _cprofile_text(res.pop("_cprofile"), args.cprofile_top)
            else:
                res.pop("_cprofile", None)
            runs.append(res)
        agg = _aggregate(runs)
        print(format_report(agg, point_label))
        points.append({"knob": knob, "value": value, "aggregate": agg, "runs": runs})

    if knob is not None:
        print(f"\n  {knob:<34} {'seconds':>10} {'peak RSS MB':>13}")
        print("  " + "-" * 60)
        secs, rss = [], []
        for point in points:
            a = point["aggregate"]
            # ⛔⛔ THE SAMPLED PEAK, NEVER `ru_maxrss`. A sweep is exactly the case `_peak_rss_mb`'s
            # docstring warns about: it is a process-LIFETIME high-water mark, so the second point's
            # number contains the first point's and can only rise. Measured while building this: the
            # `ru_maxrss` column read 7,788 → 8,925 where the second point's own phases peaked at
            # 8,925 — the rise was real there, but the column could not have shown a FALL had the
            # bigger grid been cheaper, which makes it useless as a comparison.
            secs.append(a["total_seconds"]["mean"])
            rss.append(max(s["peak_rss_mb"]["mean"] for s in a["stages"].values()))
            print(f"  {str(point['value']):<34} {secs[-1]:>10.1f} {rss[-1]:>13,.0f}")
        if len(points) > 1:
            print(
                "\n  ⚠ ONE PROCESS RUNS EVERY POINT, and the allocator does not hand freed pages back "
                "(see `_snap_rss_current`), so a later point ENTERS higher than it would alone: the "
                "memory column is an UPPER bound on each point after the first. For a number to quote, "
                "run one point per invocation."
            )
        # ⛔ `TRAPS: an-ablation-that-never-ran` — a knob the pipeline ignores draws a flat line that
        # reads exactly like "this dial is free". Say which it was rather than leaving it to the eye.
        if len(points) > 1 and max(secs) - min(secs) == 0.0 and max(rss) - min(rss) == 0.0:
            print(f"\n  ⛔ `{knob}` moved NEITHER time nor memory across {len(points)} values. Either "
                  f"the pipeline does not read it, or the values were too close together — this is "
                  f"NOT evidence that the dial is free.")

    if cprofile_text:
        print(f"\n  cProfile — top {args.cprofile_top} by cumulative time:\n")
        print(cprofile_text)

    if args.out:
        report = {
            "label": label,
            "bam": bam_path,
            "index": str(Path(args.index).resolve()),
            "threads": args.threads,
            "sweep": args.sweep,
            "system": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "cpu_count": os.cpu_count(),
            },
            "points": points,
        }
        Path(args.out).write_text(json.dumps(report, indent=2))
        print(f"\n  JSON report → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
