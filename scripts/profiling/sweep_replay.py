#!/usr/bin/env python
"""Replay one calibration sweep in isolation, on a real library, and prove a faster sweep changes nothing.

On a deep library the calibration sweep is most of the run, and a whole pipeline run is far too slow a
loop to optimise it in. `capture` runs the pipeline once and pickles every `solve_chain` call's inputs and
its returned belief; `replay` loads one captured call, runs the current `solve_chain` on it, times it, and
compares every array of the result with the captured one bit for bit, so a change to the sweep is proven a
numeric no-op on real data in the time of one sweep. `--cprofile` ranks the Python inside that one sweep
without the rest of the pipeline diluting it. The capture is a statement about the code that captured it:
re-capture after any change that is allowed to move numbers.

Usage::

    python scripts/profiling/sweep_replay.py capture --bam lib.bam --index idx/ --out DIR [--threads 8]
    python scripts/profiling/sweep_replay.py replay --dir DIR [--call 0] [--cprofile out.prof] [--block-slots N]
"""

from __future__ import annotations

import argparse
import cProfile
import dataclasses
import io
import os
import pickle
import pstats
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


class FrozenRows:
    """Unpickling shim for captures taken before 2026-09-11, whose policy carried a ``rows_at`` closure
    frozen as this table; the intron factory's rows now travel on the sweep's context and a capture
    pickles the policy as it is."""

    def __init__(self, table: dict):
        self.table = table

    def __call__(self, n_grid, window):
        return self.table[(int(n_grid), float(window))]


def capture(bam: str, index_dir: str, out: Path, threads: int | None) -> int:
    """Run the pipeline once, pickling each sweep's (args, kwargs) and returned belief into ``out``."""
    import importlib

    # the package re-exports the function under the module's own name, so ``import ... as`` would bind
    # the function; the module that calls ``solve_chain`` is the one in sys.modules
    cal = importlib.import_module("rigel.calibration.calibrate")
    sweep = importlib.import_module("rigel.calibration.sweep")
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import run_pipeline

    out.mkdir(parents=True, exist_ok=True)
    original = sweep.solve_chain
    calls = {"n": 0}

    def recording(*args, **kwargs):
        k = calls["n"]
        calls["n"] += 1
        result = original(*args, **kwargs)
        with open(out / f"sweep_{k}.in.pkl", "wb") as fh:
            pickle.dump((args, kwargs), fh, protocol=5)
        with open(out / f"sweep_{k}.out.pkl", "wb") as fh:
            pickle.dump(result, fh, protocol=5)
        return result

    cal.solve_chain = recording
    try:
        cfg = PipelineConfig()
        if threads is not None:
            cfg = dataclasses.replace(cfg, scan=dataclasses.replace(cfg.scan, total_threads=threads),
                                      em=dataclasses.replace(cfg.em, n_threads=threads))
        run_pipeline(bam, TranscriptIndex.load(index_dir), cfg)
    finally:
        cal.solve_chain = original
    print(f"captured {calls['n']} sweep call(s) -> {out}")
    return 0 if calls["n"] else 1


def _compare(a, b, path="belief") -> list[str]:
    """Every field that differs between two results, recursing into dataclasses, dicts and sequences."""
    if dataclasses.is_dataclass(a) and dataclasses.is_dataclass(b):
        diffs = []
        for f in dataclasses.fields(a):
            diffs += _compare(getattr(a, f.name), getattr(b, f.name), f"{path}.{f.name}")
        return diffs
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        same = a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b, equal_nan=a.dtype.kind == "f")
        return [] if same else [f"{path}: {a.dtype}{a.shape} vs {b.dtype}{b.shape}"]
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return [f"{path}: keys differ"]
        return [d for k in a for d in _compare(a[k], b[k], f"{path}[{k!r}]")]
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in _compare(x, y, f"{path}[{i}]")]
    if isinstance(a, float) and isinstance(b, float) and a != a and b != b:
        return []
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def replay(directory: Path, call: int, cprofile_path: str | None, block_slots: str | None = None) -> int:
    """Run the current ``solve_chain`` on one captured call and compare with the captured result.

    ``block_slots`` overrides the sweep's locus-block size (``"none"`` for the whole chain as one block):
    the answer must not depend on it, so a replay at any value against a capture taken at another is the
    chunk-exactness of the whole sweep, proven on real data."""
    import rigel.calibration.sweep as sweep

    with open(directory / f"sweep_{call}.in.pkl", "rb") as fh:
        args, kwargs = pickle.load(fh)
    if block_slots is not None:
        kwargs = dict(kwargs, block_slots=None if block_slots.lower() == "none" else int(block_slots))
    with open(directory / f"sweep_{call}.out.pkl", "rb") as fh:
        expected = pickle.load(fh)
    profiler = cProfile.Profile() if cprofile_path else None
    t0 = time.perf_counter()
    if profiler is not None:
        profiler.enable()
    result = sweep.solve_chain(*args, **kwargs)
    if profiler is not None:
        profiler.disable()
        profiler.dump_stats(cprofile_path)
    seconds = time.perf_counter() - t0
    diffs = _compare(expected, result)
    tag = "" if block_slots is None else f" [block_slots={block_slots}]"
    print(f"  sweep {call}:{tag} {seconds:.2f} s   " + ("BIT-IDENTICAL" if not diffs else f"{len(diffs)} FIELD(S) DIFFER"))
    for d in diffs[:12]:
        print(f"     {d}")
    if profiler is not None:
        buf = io.StringIO()
        pstats.Stats(cprofile_path, stream=buf).sort_stats("tottime").print_stats(35)
        print(buf.getvalue())
    return 0 if not diffs else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("capture", help="run the pipeline once and pickle every sweep's inputs and result")
    c.add_argument("--bam", required=True)
    c.add_argument("--index", required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--threads", type=int, default=None)
    r = sub.add_parser("replay", help="replay one captured sweep and compare its result bit for bit")
    r.add_argument("--dir", type=Path, required=True)
    r.add_argument("--call", type=int, default=0)
    r.add_argument("--cprofile", default=None, metavar="OUT.prof")
    r.add_argument("--block-slots", default=None, metavar="N|none",
                   help="override the locus-block size the replayed sweep solves in; the answer must not move")
    args = ap.parse_args()
    if args.cmd == "capture":
        return capture(args.bam, args.index, args.out, args.threads)
    return replay(args.dir, args.call, args.cprofile, args.block_slots)


if __name__ == "__main__":
    sys.exit(main())
