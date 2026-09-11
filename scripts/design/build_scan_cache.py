"""Scan once, calibrate many times: cache what calibration reads from a BAM scan, one cache per condition.

Calibration is the phase under development and re-running it must cost seconds, not a re-scan, so this
script runs the single-pass scan on each condition's oracle BAM (or on one BAM) and writes the payload
and strand model to a `rigel.scan_cache` directory. Each cache is independently keyed against the index
it was built from — `graph_hash`, a reach digest (which neither `partition_hash` nor `graph_hash`
covers) and the scan config — and is refused at load if the index does not match, so a partial suite is
usable and one condition can be rebuilt without touching the others. Every cache is read straight back
after it is written so a cache that cannot be loaded fails here rather than in a later session. The
cache key hashes the accumulator's deposit rule and not the fragment construction, so a change to which
fragments are offered needs `--force`. Nothing is calibrated or scored.

Usage::

    python scripts/design/build_scan_cache.py --index IDX --suite SUITE                    # every condition with an oracle BAM
    python scripts/design/build_scan_cache.py --index IDX --suite SUITE --conditions A B   # a subset
    python scripts/design/build_scan_cache.py --index IDX --suite SUITE --force            # rebuild existing caches
    python scripts/design/build_scan_cache.py --index IDX --bam X.bam --out CACHE_DIR      # a single BAM
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
    _stats, strand_model, _buffer, payload = scan_and_buffer(str(bam), index, scan)
    scan_seconds = time.perf_counter() - start

    write_scan_cache(
        out_dir,
        payload=payload,
        strand_model=strand_model,
        index=index,
        bam=str(bam),
        scan_config=scan,
    )
    # Read it straight back: a cache that cannot be loaded against the index it was just written
    # from is worse than no cache, and the failure should surface here.
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
