#!/usr/bin/env python
"""Re-measure FRAGMENT_LENGTH_AUDIT.md §2's gap table against the CURRENT anchor.

⭐ **The falsification this area never had.** On a zero-gDNA condition every fragment is RNA, so the
unconditional anchor and the RNA pool describe **one population** and any gap between them is bias.
§2 measured that gap at **+11.6 % mean / +71.1 % sd** with the scanner's histogram as the anchor; C1
built the accumulator's own histogram and measured **+7.7 % / +32.0 %**, isolating the residual as the
junction-opportunity tilt (C3's target). This script says which of those the *shipped* code is on.

⚠ It reads the anchor **through ``build_fl_models``**, not off the payload directly — the question is
what the tool is wired to use, not what is available to it.

Usage::

    python scripts/design/fl_anchor_gap.py [--index DIR] [--pilot DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def moments(counts: np.ndarray) -> tuple[float, float, int]:
    """Mean, sd and support ceiling of a raw 1-bp histogram."""
    c = np.asarray(counts, dtype=np.float64)
    total = c.sum()
    if total <= 0:
        return float("nan"), float("nan"), 0
    bins = np.arange(c.size, dtype=np.float64)
    mean = float((bins * c).sum() / total)
    var = float((bins * bins * c).sum() / total) - mean * mean
    nz = np.nonzero(c)[0]
    return mean, float(np.sqrt(max(var, 0.0))), int(nz[-1]) if nz.size else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    args = ap.parse_args()

    if not args.pilot.is_dir():
        print(f"no pilot scan-cache dir at {args.pilot}", file=sys.stderr)
        return 2

    from rigel.calibration.fl import build_fl_models, rna_fl_mass
    from rigel.index import TranscriptIndex
    from rigel.scan_cache import read_scan_cache

    # ⚠ A cache is REFUSED unless it describes the index it is loaded against (graph_hash, reach_digest
    # and payload_schema_digest). That refusal is the point of the keys, so it is reported, not caught.
    index = TranscriptIndex.load(str(args.index))

    conditions = sorted(p for p in args.pilot.iterdir() if p.is_dir())
    print(f"{'condition':<44} {'anchor mean':>11} {'RNA mean':>9} {'d%':>7} "
          f"{'anchor sd':>10} {'RNA sd':>8} {'d%':>7} {'ceilings':>12}")
    print("-" * 116)

    for cond in conditions:
        try:
            cache = read_scan_cache(cond, index)
        except Exception as exc:  # noqa: BLE001 — a refused cache is a result, not a crash
            print(f"{cond.name:<44} REFUSED: {exc}")
            continue

        fl = build_fl_models(cache.payload)
        a_mean, a_sd, a_top = moments(fl.global_counts)
        r_mean, r_sd, r_top = moments(rna_fl_mass(cache.payload))
        d_mean = 100.0 * (r_mean - a_mean) / a_mean
        d_sd = 100.0 * (r_sd - a_sd) / a_sd
        print(f"{cond.name:<44} {a_mean:11.1f} {r_mean:9.1f} {d_mean:+6.1f}% "
              f"{a_sd:10.1f} {r_sd:8.1f} {d_sd:+6.1f}% {a_top:5d} vs {r_top:<5d}")

    print()
    print("⭐ The ZERO-gDNA rows (`gdna_none_*`) are the falsification: every fragment there is RNA,")
    print("   so anchor and RNA pool describe ONE population and the gap is bias, not composition.")
    print("   §2 (scanner anchor): +11.6 % / +71.1 %.   C1 (accumulator anchor): +7.7 % / +32.0 %.")
    print("   The residual is the junction-opportunity tilt — C3's target, not C2's.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
