"""Runner: execute ``rigel quant`` across conditions with named configs."""
from __future__ import annotations

import logging
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from .config import BenchmarkConfig, RigelNamedConfig

logger = logging.getLogger(__name__)


@dataclass
class RunResult:
    """Result of a single rigel quant invocation."""

    condition: str
    config_name: str
    success: bool
    elapsed: float = 0.0
    skipped: bool = False
    returncode: int = 0
    error: str = ""


@dataclass
class RunSummary:
    """Aggregate results from a batch of runs."""

    results: list[RunResult] = field(default_factory=list)

    @property
    def n_total(self) -> int:
        return len(self.results)

    @property
    def n_success(self) -> int:
        return sum(1 for r in self.results if r.success and not r.skipped)

    @property
    def n_skipped(self) -> int:
        return sum(1 for r in self.results if r.skipped)

    @property
    def n_failed(self) -> int:
        return sum(1 for r in self.results if not r.success and not r.skipped)

    def print_summary(self) -> None:
        """Print a human-readable summary."""
        print(f"\n{'='*60}")
        print(f"Run Summary: {self.n_success} succeeded, "
              f"{self.n_skipped} skipped, {self.n_failed} failed "
              f"(of {self.n_total} total)")
        print(f"{'='*60}")

        if self.n_failed > 0:
            print("\nFailed runs:")
            for r in self.results:
                if not r.success and not r.skipped:
                    print(f"  {r.condition} / {r.config_name}: "
                          f"exit={r.returncode} — {r.error}")

        if self.n_success > 0:
            elapsed_runs = [r for r in self.results if r.success and not r.skipped]
            total_time = sum(r.elapsed for r in elapsed_runs)
            print(f"\nTotal runtime: {total_time:.1f}s "
                  f"(avg {total_time / len(elapsed_runs):.1f}s per run)")


def _check_output_exists(outdir: Path) -> bool:
    """Check if a completed rigel output already exists."""
    return (outdir / "quant.feather").exists()


def run_rigel_config(
    cfg: BenchmarkConfig,
    condition: str,
    named_config: RigelNamedConfig,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> RunResult:
    """Run ``rigel quant`` for a single condition + config combination.

    Parameters
    ----------
    cfg : BenchmarkConfig
        Top-level benchmark configuration.
    condition : str
        Condition name (e.g. ``gdna_none_ss_1.00_nrna_none``).
    named_config : RigelNamedConfig
        The rigel configuration to run.
    force : bool
        Re-run even if output already exists.
    dry_run : bool
        Print the command but don't execute it.

    Returns
    -------
    RunResult
    """
    outdir = cfg.output_dir(condition, named_config.name)
    bam_path = cfg.bam_path(condition, named_config)

    # Check if already done
    if not force and _check_output_exists(outdir):
        logger.info("Skipping %s / %s — output exists", condition, named_config.name)
        return RunResult(
            condition=condition,
            config_name=named_config.name,
            success=True,
            skipped=True,
        )

    # Validate BAM exists
    if not bam_path.exists():
        msg = f"BAM not found: {bam_path}"
        logger.error(msg)
        return RunResult(
            condition=condition,
            config_name=named_config.name,
            success=False,
            error=msg,
        )

    # Build command
    cmd = [
        "rigel", "quant",
        "--bam", str(bam_path),
        "--index", str(cfg.rigel_index),
        "-o", str(outdir),
    ]

    # Add global defaults
    if "threads" not in named_config.params:
        cmd.extend(["--threads", str(cfg.threads)])
    if "seed" not in named_config.params:
        cmd.extend(["--seed", str(cfg.seed)])

    # Add config-specific params
    cmd.extend(named_config.to_cli_args())

    cmd_str = shlex.join(cmd)

    if dry_run:
        print(f"[dry-run] {cmd_str}")
        return RunResult(
            condition=condition,
            config_name=named_config.name,
            success=True,
            skipped=True,
        )

    # Execute
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"  Running: {condition} / {named_config.name}", flush=True)
    logger.info("Command: %s", cmd_str)

    t0 = time.monotonic()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )
        elapsed = time.monotonic() - t0

        if result.returncode != 0:
            stderr_tail = result.stderr[-500:] if result.stderr else ""
            logger.error(
                "%s / %s failed (exit %d): %s",
                condition, named_config.name, result.returncode, stderr_tail,
            )
            return RunResult(
                condition=condition,
                config_name=named_config.name,
                success=False,
                elapsed=elapsed,
                returncode=result.returncode,
                error=stderr_tail,
            )

        logger.info(
            "%s / %s completed in %.1fs", condition, named_config.name, elapsed
        )
        return RunResult(
            condition=condition,
            config_name=named_config.name,
            success=True,
            elapsed=elapsed,
        )

    except FileNotFoundError:
        return RunResult(
            condition=condition,
            config_name=named_config.name,
            success=False,
            error="'rigel' command not found. Is the conda env activated?",
        )


def run_all(
    cfg: BenchmarkConfig,
    *,
    force: bool = False,
    dry_run: bool = False,
    config_names: list[str] | None = None,
    condition_filter: list[str] | None = None,
) -> RunSummary:
    """Run all condition × config combinations.

    Parameters
    ----------
    cfg : BenchmarkConfig
        Benchmark configuration.
    force : bool
        Re-run even if output exists.
    dry_run : bool
        Print commands without executing.
    config_names : list[str] or None
        Subset of named configs to run (default: all).
    condition_filter : list[str] or None
        Subset of conditions (default: all or from config).
    """
    summary = RunSummary()

    conditions = condition_filter or cfg.get_conditions()
    configs_to_run = cfg.rigel_configs
    if config_names:
        configs_to_run = {
            k: v for k, v in configs_to_run.items() if k in config_names
        }
        unknown = set(config_names) - set(configs_to_run)
        if unknown:
            logger.warning("Unknown config names: %s", sorted(unknown))

    if not configs_to_run:
        print("No rigel configs defined — nothing to run.", file=sys.stderr)
        return summary

    print(f"Running {len(configs_to_run)} config(s) × "
          f"{len(conditions)} condition(s) = "
          f"{len(configs_to_run) * len(conditions)} total runs\n")

    for config_name, named_config in sorted(configs_to_run.items()):
        print(f"\n[config: {config_name}]")
        for condition in conditions:
            result = run_rigel_config(
                cfg, condition, named_config,
                force=force, dry_run=dry_run,
            )
            summary.results.append(result)

    summary.print_summary()
    return summary
