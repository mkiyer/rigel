"""CLI dispatcher for the benchmarking package.

Usage::

    python -m scripts.benchmarking run     -c benchmark.yaml [--force] [--conditions ...]
    python -m scripts.benchmarking analyze -c benchmark.yaml [-o results/report]
    python -m scripts.benchmarking status  -c benchmark.yaml
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

from .config import load_config


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-c", "--config", required=True,
        help="Path to benchmark YAML config file",
    )
    parser.add_argument(
        "--conditions", nargs="*", default=None,
        help="Subset of conditions to process (default: all)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose logging",
    )


# ─── align ───────────────────────────────────────────────────────


def cmd_align(args: argparse.Namespace) -> int:
    from .aligner import run_all_alignments

    cfg = load_config(args.config)
    t0 = time.monotonic()
    summary = run_all_alignments(
        cfg,
        aligners=args.aligners,
        force=args.force,
        dry_run=args.dry_run,
        condition_filter=args.conditions,
    )
    elapsed = time.monotonic() - t0
    print(f"\nTotal elapsed: {elapsed:.1f}s")
    return 1 if summary.n_failed > 0 else 0


# ─── run ─────────────────────────────────────────────────────────


def cmd_run(args: argparse.Namespace) -> int:
    from .runner import run_all

    cfg = load_config(args.config)
    t0 = time.monotonic()
    summary = run_all(
        cfg,
        force=args.force,
        dry_run=args.dry_run,
        config_names=args.config_names,
        condition_filter=args.conditions,
    )
    elapsed = time.monotonic() - t0
    print(f"\nTotal elapsed: {elapsed:.1f}s")
    return 1 if summary.n_failed > 0 else 0


# ─── analyze ─────────────────────────────────────────────────────


def cmd_run_tools(args: argparse.Namespace) -> int:
    from .tools import run_all_tools

    cfg = load_config(args.config)
    t0 = time.monotonic()
    summary = run_all_tools(
        cfg,
        tools=args.tools,
        force=args.force,
        dry_run=args.dry_run,
        condition_filter=args.conditions,
    )
    elapsed = time.monotonic() - t0
    print(f"\nTotal elapsed: {elapsed:.1f}s")
    return 1 if summary.n_failed > 0 else 0


# ─── analyze ─────────────────────────────────────────────────────


def cmd_analyze(args: argparse.Namespace) -> int:
    from .analysis import run_analysis

    cfg = load_config(args.config)
    outdir = Path(args.output) if args.output else None
    t0 = time.monotonic()
    run_analysis(cfg, outdir=outdir, condition_filter=args.conditions)
    elapsed = time.monotonic() - t0
    print(f"\nCompleted in {elapsed:.1f}s")
    return 0


# ─── status ──────────────────────────────────────────────────────


def cmd_status(args: argparse.Namespace) -> int:
    cfg = load_config(args.config)

    print(f"Benchmark directory: {cfg.benchmark_dir}")
    print(f"Rigel index: {cfg.rigel_index}")
    print(f"BAM relative path: {cfg.bam_relpath}")
    print(f"Threads: {cfg.threads}  Seed: {cfg.seed}")
    print()

    # Named configs
    if cfg.rigel_configs:
        print(f"Rigel configs ({len(cfg.rigel_configs)}):")
        for name, nc in sorted(cfg.rigel_configs.items()):
            params_str = ", ".join(f"{k}={v}" for k, v in nc.params.items())
            print(f"  {name}: {params_str or '(defaults)'}")
    else:
        print("No rigel configs defined.")
    print()

    # Conditions
    try:
        all_conditions = cfg.discover_conditions()
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    selected = args.conditions or cfg.get_conditions()
    print(f"Conditions: {len(selected)} selected / {len(all_conditions)} available")
    print()

    # Status matrix
    from .analysis import discover_tools

    header = ["condition", "bam"]
    config_names = sorted(cfg.rigel_configs) if cfg.rigel_configs else []
    header.extend(f"rigel/{n}" for n in config_names)
    header.append("other_tools")

    rows = []
    for cond in selected:
        cond_dir = cfg.condition_dir(cond)
        bam_ok = cfg.bam_path(cond).exists()

        # Check each named config output
        config_status = []
        for cn in config_names:
            outdir = cfg.output_dir(cond, cn)
            if (outdir / "quant.feather").exists():
                config_status.append("done")
            else:
                config_status.append("—")

        # Discover other tools
        tools = discover_tools(cond_dir)
        other = [t.name for t in tools if not t.name.startswith("rigel/")]
        other_str = ", ".join(other) if other else "—"

        row = [cond, "ok" if bam_ok else "MISSING"]
        row.extend(config_status)
        row.append(other_str)
        rows.append(row)

    # Print table
    col_widths = [
        max(len(str(row[i])) for row in [header] + rows)
        for i in range(len(header))
    ]
    fmt = "  ".join(f"{{:<{w}}}" for w in col_widths)
    print(fmt.format(*header))
    print(fmt.format(*["─" * w for w in col_widths]))
    for row in rows:
        print(fmt.format(*row))

    # Summary counts
    n_done = sum(1 for r in rows for s in r[2:-1] if s == "done")
    n_pending = sum(1 for r in rows for s in r[2:-1] if s == "—")
    n_missing_bam = sum(1 for r in rows if r[1] == "MISSING")
    print()
    print(f"Done: {n_done}  Pending: {n_pending}  Missing BAM: {n_missing_bam}")

    return 0


# ─── main ────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="python -m scripts.benchmarking",
        description="Rigel full-scale simulation benchmarking",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # align
    align_p = subparsers.add_parser("align", help="Align FASTQ reads (minimap2)")
    _add_common_args(align_p)
    align_p.add_argument("--force", action="store_true", help="Re-align even if BAM exists")
    align_p.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    align_p.add_argument(
        "--aligners", nargs="*", default=None,
        help="Aligners to run (default: minimap2). Options: minimap2, oracle",
    )
    align_p.set_defaults(func=cmd_align)

    # run
    run_p = subparsers.add_parser("run", help="Run rigel quant with named configs")
    _add_common_args(run_p)
    run_p.add_argument("--force", action="store_true", help="Re-run even if output exists")
    run_p.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    run_p.add_argument(
        "--configs", dest="config_names", nargs="*", default=None,
        help="Subset of named configs to run (default: all)",
    )
    run_p.set_defaults(func=cmd_run)

    # run-tools
    tools_p = subparsers.add_parser("run-tools", help="Run external tools (salmon, kallisto)")
    _add_common_args(tools_p)
    tools_p.add_argument("--force", action="store_true", help="Re-run even if output exists")
    tools_p.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    tools_p.add_argument(
        "--tools", nargs="*", default=None,
        help="Tools to run (default: auto-detect from config). Options: salmon, kallisto",
    )
    tools_p.set_defaults(func=cmd_run_tools)

    # analyze
    analyze_p = subparsers.add_parser("analyze", help="Analyze benchmark results")
    _add_common_args(analyze_p)
    analyze_p.add_argument("-o", "--output", default=None, help="Output directory for reports")
    analyze_p.set_defaults(func=cmd_analyze)

    # status
    status_p = subparsers.add_parser("status", help="Show benchmark status")
    _add_common_args(status_p)
    status_p.set_defaults(func=cmd_status)

    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
