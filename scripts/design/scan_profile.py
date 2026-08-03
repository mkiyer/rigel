"""What the accumulator costs per fragment, measured on a real BAM.

    

⚠ WHY THIS EXISTS. There was no scan-profiling harness in this repository, and the only recorded figure
— *"the BAM scan is ~2 % of runtime"* — was measured on a 39 MB BAM with 155,352 read groups. That
percentage does not generalise: at 10^8 fragments a 124 ns deposit is 12.4 s of single-thread CPU on its
own. **The accumulator's budget must be stated per fragment, not per run.**

HOW THE DEPOSIT IS ISOLATED. The scanner only accumulates when ``set_regions`` has been called
(``pipeline._wire_calibration_regions``). So the same BAM is scanned twice — once with the accumulator
wired and once without — and the difference, divided by the fragment count, is the per-fragment deposit
cost. Everything else (htslib parsing, fragment assembly, overlap resolution, model training, columnar
buffering) is common to both runs and cancels.

⚠ TWO COSTS WITH DIFFERENT SCALING, AND THEY MUST NOT BE AVERAGED TOGETHER. Wiring the accumulator adds
work that is **O(fragments)** — the deposit itself — and work that is **O(partition)**: allocating and
zeroing a per-worker accumulator, merging the workers, and copying the per-reference structs into the
flat payload. On a shallow BAM against the 1.04 M-node human partition the second term dominates, so a
single-BAM ``delta / n_fragments`` is not a per-fragment cost.

So pass **two or more BAMs of different depth against the same index** and the harness regresses::

    delta(n)  =  fixed_partition_cost  +  per_fragment_cost * n

the slope being the number that generalises and the intercept the price of the partition itself.

    OMP_NUM_THREADS=1 python scripts/design/scan_profile.py INDEX BAM [BAM ...] [--repeats N]

⚠ Run it on a quiet machine and with ``OMP_NUM_THREADS=1`` unless you are deliberately measuring the
parallel path: the deposit is per-worker with a serial merge, so a threaded number answers a different
question.
"""

from __future__ import annotations

import argparse
import statistics
import time
from pathlib import Path

from rigel import pipeline
from rigel.calibration.splice_graph import build_node_partition_arrays
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex


def _timed_scan(bam: str, index: TranscriptIndex, scan: BamScanConfig, accumulate: bool):
    """One scan, returning ``(seconds, stats, payload)``. With ``accumulate=False`` the accumulator is
    never wired, so the scan does everything except deposit."""
    original = pipeline._wire_calibration_regions
    if not accumulate:
        pipeline._wire_calibration_regions = lambda *a, **k: None
    try:
        t0 = time.perf_counter()
        stats, _strand, _fl, _buffer, payload = pipeline.scan_and_buffer(bam, index, scan)
        return time.perf_counter() - t0, stats, payload
    finally:
        pipeline._wire_calibration_regions = original


def _profile_one(bam: str, index, scan, repeats: int):
    """``(n_fragments, scan_seconds, accumulator_seconds)`` for one BAM, minima over ``repeats``."""
    arms = {}
    n_fragments = 0
    for accumulate in (False, True):
        times = []
        for _ in range(repeats):
            seconds, stats, payload = _timed_scan(bam, index, scan, accumulate)
            times.append(seconds)
            if accumulate:
                n_fragments = int(stats.n_fragments)
                if payload is None:
                    raise SystemExit(f"{bam}: no payload; set_regions did not run")
        arms[accumulate] = min(times)
    return n_fragments, arms[True], arms[True] - arms[False]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index")
    ap.add_argument("bams", nargs="+")
    ap.add_argument("--repeats", type=int, default=3, help="scans per arm; the MINIMUM is used")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    scan = BamScanConfig(total_threads=args.threads)
    cuts, _offsets, _types = build_node_partition_arrays(index)

    print(f"index      {args.index}")
    print(f"partition  {len(index.nodes_df):,} nodes, {cuts.size:,} cut positions")
    print(f"threads    {args.threads}   repeats {args.repeats} (minimum used)\n")
    print(f"{'BAM':<26} {'fragments':>12} {'scan':>9} {'accumulator':>12} {'delta/frag':>12}")
    print("-" * 76)

    points = []
    for bam in args.bams:
        n, scan_s, acc_s = _profile_one(bam, index, scan, args.repeats)
        if n == 0:
            print(f"{Path(bam).name[:26]:<26} {'no fragments':>12}")
            continue
        points.append((n, acc_s))
        print(
            f"{Path(bam).name[:26]:<26} {n:>12,} {scan_s:>8.3f}s {acc_s:>11.3f}s "
            f"{1e9 * acc_s / n:>10.1f}ns"
        )

    if len(points) < 2:
        print(
            "\n⚠ ONE BAM ONLY, so the two costs cannot be separated. The ns/fragment above still "
            "contains the O(partition) allocation, merge and payload copy, and is therefore an UPPER "
            "BOUND on the true per-fragment deposit. Pass a second BAM of different depth."
        )
        return

    # least squares on delta(n) = fixed + c*n
    xs = [float(n) for n, _ in points]
    ys = [s for _, s in points]
    mean_x, mean_y = statistics.fmean(xs), statistics.fmean(ys)
    var_x = sum((x - mean_x) ** 2 for x in xs)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / var_x
    intercept = mean_y - slope * mean_x

    print(f"\nregression over {len(points)} BAMs:  delta(n) = fixed + c*n")
    print(
        f"  ⭐ c      = {1e9 * slope:8.1f} ns/fragment      <- the deposit, the number that generalises"
    )
    print(f"  ⭐ fixed  = {intercept:8.3f} s                 <- the price of the partition itself")
    if slope <= 0:
        print(
            "\n⚠ NEGATIVE SLOPE: the per-fragment cost is below this harness's noise floor at these "
            "depths. Use deeper BAMs or raise --repeats; do not quote a per-fragment budget from this."
        )
        return
    print("\nprojected single-thread DEPOSIT cost, partition cost excluded:")
    for depth in (1e6, 1e7, 1e8, 1e9):
        print(f"    {depth:.0e} fragments   {slope * depth:8.1f} s")


if __name__ == "__main__":
    main()
