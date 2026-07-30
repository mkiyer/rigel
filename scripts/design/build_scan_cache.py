"""Scan a suite ONCE and cache what calibration needs, so calibration can be iterated on without rescanning.

    TODO item 2 (the cached substrate)   ·   `docs/testing/testing_plan.md`   ·   Library: `rigel.scan_cache`

⭐ **THE POINT.** Calibration is the phase under development and it is the expensive one — index load
~8 s, BAM scan ~2 s, **calibration ~66 s** on a real cfRNA library (`CARRY_FORWARD.md` §1 fact 22) — while
a 5 M-fragment simulated condition costs far more than that to scan. Caching the scan took a 24-condition
sweep from ~13 min to ~9 s on the old path. One cache per condition, each independently keyed and valid,
so a partial suite is usable and one condition can be rebuilt without touching the others.

    python scripts/design/build_scan_cache.py --index IDX --suite SUITE [--conditions A B ...]
    python scripts/design/build_scan_cache.py --index IDX --bam X.bam --out CACHE_DIR

⚠ **A cache is refused at load unless it describes the index it is loaded against** — `graph_hash`, a
**reach digest** (which neither `partition_hash` nor `graph_hash` covers, and a rebuild once moved ~38 %
of human contiguous reaches with both byte-identical), and the scan config. See `rigel.scan_cache`.
"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from rigel.config import BamScanConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer  # noqa: E402
from rigel.scan_cache import read_scan_cache, write_scan_cache  # noqa: E402

#: Where a simulated condition keeps its oracle BAM, relative to the condition directory.
ORACLE_BAM = "sim_oracle.bam"


def cache_one(index: TranscriptIndex, bam: Path, out_dir: Path) -> tuple[float, float, int]:
    """Scan *bam* against *index* and write its cache. Returns (scan seconds, load seconds, bytes)."""
    scan = BamScanConfig(sj_strand_tag=_native_detect_sj_tag(str(bam)))

    start = time.perf_counter()
    _stats, strand_model, fl_models, _buffer, payload = scan_and_buffer(str(bam), index, scan)
    scan_seconds = time.perf_counter() - start

    write_scan_cache(
        out_dir,
        payload=payload,
        strand_model=strand_model,
        frag_length_models=fl_models,
        index=index,
        bam=str(bam),
        scan_config=scan,
    )
    # ⚠ Read it straight back: a cache that cannot be loaded against the index it was just written
    # from is worse than no cache, and the failure should surface here rather than three days later.
    start = time.perf_counter()
    read_scan_cache(out_dir, index)
    load_seconds = time.perf_counter() - start

    size = sum(p.stat().st_size for p in out_dir.iterdir() if p.is_file())
    return scan_seconds, load_seconds, size


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--suite", type=Path, default=None, help="Suite dir; caches every condition in it")
    ap.add_argument("--bam", type=Path, default=None, help="A single BAM instead of a suite")
    ap.add_argument("--out", type=Path, default=None, help="Cache dir (single-BAM mode)")
    ap.add_argument("--conditions", nargs="*", default=None, help="Subset of condition names")
    ap.add_argument("--force", action="store_true", help="Rebuild caches that already exist")
    args = ap.parse_args()

    if (args.suite is None) == (args.bam is None):
        raise SystemExit("give exactly one of --suite or --bam")
    if args.bam is not None and args.out is None:
        raise SystemExit("--bam needs --out")

    start = time.perf_counter()
    index = TranscriptIndex.load(str(args.index))
    print(f"index loaded in {time.perf_counter() - start:.1f} s  ({args.index})")

    if args.bam is not None:
        jobs = [(args.bam.stem, args.bam, args.out)]
    else:
        names = args.conditions or sorted(
            p.name for p in args.suite.iterdir() if (p / ORACLE_BAM).exists()
        )
        jobs = [(name, args.suite / name / ORACLE_BAM, args.suite / "scan_cache" / name)
                for name in names]
        if not jobs:
            raise SystemExit(f"no condition under {args.suite} has a {ORACLE_BAM}")

    print(f"\n{'condition':<48} {'scan':>9} {'load':>9} {'size':>10}")
    total_scan = 0.0
    total_load = 0.0
    for name, bam, out_dir in jobs:
        if out_dir.exists() and not args.force:
            print(f"{name:<48} {'cached':>9} {'':>9} {'skip':>10}")
            continue
        if not bam.exists():
            print(f"{name:<48} {'MISSING':>9} {'':>9} {'':>10}")
            continue
        scan_seconds, load_seconds, size = cache_one(index, bam, out_dir)
        total_scan += scan_seconds
        total_load += load_seconds
        print(f"{name:<48} {scan_seconds:>8.1f}s {load_seconds:>8.2f}s {size / 1e6:>8.1f} MB")

    if total_scan > 0:
        print(
            f"\n⭐ scanned once in {total_scan:.1f} s; every later calibration run reloads in "
            f"{total_load:.2f} s — {total_scan / max(total_load, 1e-9):.0f}x"
        )


if __name__ == "__main__":
    main()
