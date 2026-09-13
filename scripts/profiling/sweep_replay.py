#!/usr/bin/env python
"""Replay one calibration sweep in isolation, on a real library, and prove a faster sweep changes nothing.

On a deep library the calibration sweep is most of the run, and a whole pipeline run is far too slow a
loop to optimise it in. `capture` runs the pipeline once and pickles every `solve_chain` call's inputs and
its returned belief; `replay` loads one captured call, runs the current `solve_chain` on it, times it, and
compares every array of the result with the captured one bit for bit, so a change to the sweep is proven a
numeric no-op on real data in the time of one sweep. `--cprofile` ranks the Python inside that one sweep
without the rest of the pipeline diluting it. The capture is a statement about the code that captured it:
re-capture after any change that is allowed to move numbers.

`--tolerance` adds the second verdict a change that MAY move numbers needs (a compiled port, a solver
rewritten in one precision): per output array the slots moved, the largest absolute and relative move, and
the DERIVED budget beside them as a sanity bound. The budget is the read-out's conditioning: ψ's answer is a
posterior mean over the grid, so an error ``δ_k`` in cell ``k``'s log-density moves a fraction by at most
``½·max|δ_k|`` (σ lies in [0, 1], so its mean absolute deviation about its mean is ≤ ½) and a log-variance by
at most ``L̃²·max|δ_k|`` (``L̃ = log(1 + e^L)``, the range of ``log σ`` on the window); an implementation
that keeps the solver's expressions perturbs each of the ``T`` terms of a cell by one rounding unit of its
magnitude, and every term is bounded by ``A = max(c_κ·N, K²/2)`` — the strand Gaussian's scale
``c_κ·n`` with ``c_κ = 1/(2κ(1−κ))`` and ``N`` the largest slot count on the chain (strand
overdispersion above zero caps it at ``c_κ/od`` however deep the slot, so the bound is looser there), and
the fitted arms' kernel range on the ``K``-cell grid — so ``|δ_k| ≤ ε·T·A``. The rounding unit is the
solve's (float32 on the AMBIG cube while it has one). The bound is loose by construction — achievable
rounding lands orders inside it, since neighbouring cells' errors are incoherent and the read-out averages
them — and it is a bound on the MEAN; the median the fraction is read from shares it where the posterior
is unimodal, while a balanced bimodal posterior has none (its median jumps between the modes under any
perturbation). So the report always shows the actual move beside the budget, and a slot flagged beyond it
is a slot to look at, not a verdict. Gate: ``tests/test_sweep_replay_tolerance.py``.

Usage::

    python scripts/profiling/sweep_replay.py capture --bam lib.bam --index idx/ --out DIR [--threads 8]
    python scripts/profiling/sweep_replay.py replay --dir DIR [--call 0] [--cprofile out.prof] [--block-slots N] [--tolerance]
    python scripts/profiling/sweep_replay.py self-test
"""

from __future__ import annotations

import argparse
import cProfile
import dataclasses
import inspect
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


#: the terms summed into one ψ cell: the strand Gaussian, the gDNA arm, the RNA arm, the factory rows,
#: the message rows, and the cube row at an AMBIG slot — counted from `simplex_logodds`
PSI_TERMS = 6
EPS64 = float(np.finfo(np.float64).eps) / 2.0  # the rounding unit, half an ulp at 1
EPS32 = float(np.finfo(np.float32).eps) / 2.0


def budget(count_max: float, kappa: float, window: float, eps, n_grid: int = 1) -> tuple:
    """The derived per-slot budget ``(on a fraction, on a log-variance)`` — see the module docstring:
    ``½·ε·T·A`` and ``L̃²·ε·T·A`` with ``A = max(c_κ·N, K²/2)``. ``eps`` may be a scalar or a per-slot
    array."""
    c_kappa = 1.0 / (2.0 * kappa * (1.0 - kappa))
    amplitude = max(c_kappa * float(count_max), 0.5 * float(n_grid) ** 2)
    cell = np.asarray(eps, np.float64) * PSI_TERMS * amplitude
    L_tilde = float(np.log1p(np.exp(float(window))))
    return 0.5 * cell, L_tilde * L_tilde * cell


def moves(a, b) -> dict:
    """One output array against its capture: the slots moved, the largest absolute and relative move.
    Refuses a shape or dtype change — that is not a move, it is a different quantity."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        raise ValueError(f"not comparable: {a.dtype}{a.shape} vs {b.dtype}{b.shape}")
    if a.dtype.kind != "f":
        moved = a != b
        return dict(moved=int(moved.sum()), where=np.flatnonzero(moved), abs=0.0, rel=0.0)
    same = (a == b) | (np.isnan(a) & np.isnan(b))
    where = np.flatnonzero(~same)
    d = np.abs(a[where] - b[where])
    scale = np.maximum(np.abs(a[where]), np.abs(b[where]))
    rel = np.where(scale > 0.0, d / np.where(scale > 0.0, scale, 1.0), np.inf)
    return dict(
        moved=int(where.size),
        where=where,
        abs=float(d.max()) if where.size else 0.0,
        rel=float(rel.max()) if where.size else 0.0,
        per_slot=d,
    )


def tolerance_report(expected, result, args, kwargs) -> list[str]:
    """The report lines for one replayed sweep: every belief array, its moves, and the budget the slots
    that moved are held to. The budget's inputs are read off the captured call: the largest slot count,
    the strand model's κ, the window, and each slot's rounding unit (float32 where both strands are
    live — the AMBIG cube — else float64)."""
    _chain, statics, geometry, _belief, _ra = args
    counts = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1) + np.asarray(
        geometry.spliced_count, np.float64
    ).sum(axis=1)
    ambig = np.asarray(statics.free_pos, bool) & np.asarray(statics.free_neg, bool)
    eps = np.where(ambig, EPS32, EPS64)
    kappa = float(kwargs["rna_sense_frac"])
    window = float(kwargs.get("logodds_window", 10.0))
    n_grid = int(kwargs["n_grid"])
    b_frac, b_var = budget(float(counts.max()) if counts.size else 0.0, kappa, window, eps, n_grid)
    lines = [
        f"     tolerance: N = {counts.max():.0f}, κ = {kappa:.4f}, L = {window:g}, K = {n_grid}; the budget "
        f"is ½·ε·T·max(c_κ·N, K²/2) on a fraction ({b_frac.min():.2e} at float64, {b_frac.max():.2e} at "
        f"float32) and L̃² of it on a log-variance"
    ]
    for f in dataclasses.fields(expected):
        a, b = getattr(expected, f.name), getattr(result, f.name)
        if a is None or b is None:
            lines.append(f"     {f.name:<16} {'absent' if a is None and b is None else 'PRESENT IN ONE'}")
            continue
        m = moves(a, b)
        if m["moved"] == 0:
            lines.append(f"     {f.name:<16} 0 moved")
            continue
        bound = (b_var if f.name.startswith("var") else b_frac)[m["where"]]
        ratio = float((m["per_slot"] / bound).max()) if "per_slot" in m else float("nan")
        verdict = "within the budget" if ratio <= 1.0 else "BEYOND the budget — look at these slots"
        lines.append(
            f"     {f.name:<16} {m['moved']} moved   max |Δ| {m['abs']:.3e}   max rel Δ {m['rel']:.3e}   "
            f"max |Δ|/budget {ratio:.3e}   {verdict}"
        )
    return lines


def self_test() -> int:
    """The comparator's own falsification: identical reads zero moved; a one-ulp nudge in one slot reads
    one moved at one ulp, inside the budget; a nudge past the budget reads EXCEEDS; a NaN that stays a
    NaN is not a move; a shape change is refused."""
    checks = []
    rng = np.random.default_rng(3021)
    a = rng.uniform(0.0, 1.0, 50)
    a[7] = np.nan
    m = moves(a, a.copy())
    checks.append(("identical arrays read zero moved (a NaN kept is not a move)", m["moved"] == 0))
    b = a.copy()
    b[3] = np.nextafter(b[3], 1.0)
    m = moves(a, b)
    checks.append(("a one-ulp nudge reads one moved at one ulp", m["moved"] == 1 and m["where"][0] == 3
                   and m["abs"] == abs(b[3] - a[3]) and 0 < m["rel"] < 3 * EPS64))
    bf, bv = budget(1.0e6, 0.99, 10.0, EPS64)
    checks.append(("the float64 budget on a fraction at N = 1e6, κ = 0.99 is ½·ε·6·50.5·1e6",
                   np.isclose(bf, 0.5 * EPS64 * 6 * (1 / (2 * 0.99 * 0.01)) * 1e6, rtol=1e-12) and bv > bf))
    checks.append(("a one-ulp move sits inside the budget; a move of 1e-6 exceeds it",
                   m["abs"] <= bf and 1e-6 > bf))
    try:
        moves(a, a[:-1])
        checks.append(("a shape change is refused", False))
    except ValueError:
        checks.append(("a shape change is refused", True))
    bf32, _ = budget(1.0e6, 0.99, 10.0, EPS32)
    checks.append(("the float32 budget is 2^29 times the float64 one", abs(bf32 / bf - 2.0**29) < 1e-6))
    n_ok = 0
    for name, ok in checks:
        n_ok += bool(ok)
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\n{n_ok}/{len(checks)} comparator gates fire")
    return 0 if n_ok == len(checks) else 1


def replay(
    directory: Path,
    call: int,
    cprofile_path: str | None,
    block_slots: str | None = None,
    tolerance: bool = False,
) -> int:
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
    # a capture outlives the signature it was taken under: a keyword `solve_chain` no longer takes
    # is dropped and named, so an old capture still replays (the tolerance report says what moved)
    accepted = set(inspect.signature(sweep.solve_chain).parameters)
    dropped = sorted(k for k in kwargs if k not in accepted)
    if dropped:
        print(f"     replay: the capture carries {dropped}, which solve_chain no longer takes; dropped")
    result = sweep.solve_chain(*args, **{k: v for k, v in kwargs.items() if k in accepted})
    if profiler is not None:
        profiler.disable()
        profiler.dump_stats(cprofile_path)
    seconds = time.perf_counter() - t0
    diffs = _compare(expected, result)
    tag = "" if block_slots is None else f" [block_slots={block_slots}]"
    print(f"  sweep {call}:{tag} {seconds:.2f} s   " + ("BIT-IDENTICAL" if not diffs else f"{len(diffs)} FIELD(S) DIFFER"))
    for d in diffs[:12]:
        print(f"     {d}")
    if tolerance:
        for line in tolerance_report(expected, result, args, kwargs):
            print(line)
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
    r.add_argument("--tolerance", action="store_true",
                   help="beside the bit verdict, per output array: slots moved, max |Δ|, max rel Δ, the derived budget")
    sub.add_parser("self-test", help="the comparator's own falsification")
    args = ap.parse_args()
    if args.cmd == "capture":
        return capture(args.bam, args.index, args.out, args.threads)
    if args.cmd == "self-test":
        return self_test()
    return replay(args.dir, args.call, args.cprofile, args.block_slots, args.tolerance)


if __name__ == "__main__":
    sys.exit(main())
