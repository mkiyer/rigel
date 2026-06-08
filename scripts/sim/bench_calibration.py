#!/usr/bin/env python3
"""Deprecated shim — the canonical gDNA-vs-RNA benchmark now lives in ``rigel.sim.analysis``.

The calibration-era pool-soft + per-fragment-hard confusion this script used to compute has
been superseded by the **net fragment-flow** deconvolution analysis (the primary, transcript-
and locus-level metric that cancels unrecoverable symmetric misassignment). That analysis,
plus the pool fraction, the gross per-fragment confusion, transcript-abundance accuracy and the
calibration scalars, is produced by ``rigel.sim.analysis.main`` (a.k.a.
``scripts/sim/evaluate_suite.py``).

This wrapper preserves the old ``--sim-base/--run/--force`` CLI and simply forwards to the
canonical workflow so existing call sites keep working.

    python scripts/sim/bench_calibration.py --sim-base <suite> --run --force
    # ≡ python scripts/sim/evaluate_suite.py --sim-base <suite>

Outputs land in the suite dir: ``analysis_report.txt``, ``net_flow_per_transcript.tsv``,
``net_flow_per_locus.tsv``, ``condition_report.md``, ``condition_metrics.tsv``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim-base", type=Path, required=True, help="Suite dir (has manifest.json + conditions)")
    ap.add_argument("--run", action="store_true", help="Build index + run rigel quant before evaluating")
    ap.add_argument("--force", action="store_true", help="Re-quant even if rigel_out exists (with --run)")
    ap.add_argument("--conditions", nargs="*", default=None, help="Optional subset of condition names")
    args = ap.parse_args()

    from rigel.sim.analysis import discover_conditions, main as analysis_main

    print(
        "[bench_calibration] deprecated → forwarding to rigel.sim.analysis "
        "(net fragment-flow is now the primary metric).",
        file=sys.stderr,
    )

    # --force: drop cached quant so run_quant re-quants (analysis.run_quant skips when present).
    if args.run and args.force:
        for cond in discover_conditions(args.sim_base, args.conditions):
            qf = args.sim_base / cond / "rigel_out" / "quant.feather"
            if qf.exists():
                qf.unlink()

    forwarded = ["analysis", "--sim-base", str(args.sim_base)]
    if not args.run:
        forwarded.append("--skip-quant")
    if args.conditions:
        forwarded += ["--conditions", *args.conditions]

    saved_argv = sys.argv
    try:
        sys.argv = forwarded
        analysis_main()
    finally:
        sys.argv = saved_argv


if __name__ == "__main__":
    main()
