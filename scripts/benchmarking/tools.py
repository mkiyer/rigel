"""Runners for external quantification tools (salmon, kallisto)."""
from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

from .config import BenchmarkConfig

logger = logging.getLogger(__name__)


@dataclass
class ToolRunResult:
    """Result of a single external tool invocation."""

    condition: str
    tool: str
    success: bool
    elapsed: float = 0.0
    skipped: bool = False
    error: str = ""


@dataclass
class ToolRunSummary:
    """Aggregate results from external tool runs."""

    results: list[ToolRunResult] = field(default_factory=list)

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
        print(f"\n{'='*60}")
        print(f"Tool Run Summary: {self.n_success} succeeded, "
              f"{self.n_skipped} skipped, {self.n_failed} failed "
              f"(of {self.n_total} total)")
        print(f"{'='*60}")
        if self.n_failed > 0:
            print("\nFailed:")
            for r in self.results:
                if not r.success and not r.skipped:
                    print(f"  {r.condition} / {r.tool}: {r.error}")
        if self.n_success > 0:
            runs = [r for r in self.results if r.success and not r.skipped]
            total_time = sum(r.elapsed for r in runs)
            print(f"\nTotal tool runtime: {total_time:.1f}s")


def _find_tool(name: str) -> str:
    """Resolve tool name to an absolute path."""
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(f"Tool '{name}' not found in PATH")


def _salmon_output_exists(out_dir: Path) -> bool:
    """Check if salmon output is complete."""
    return (out_dir / "quant.sf.gz").exists() or (out_dir / "quant.sf").exists()


def _compress_salmon_output(out_dir: Path) -> None:
    """Gzip salmon output files for consistency with benchmark expectations."""
    import gzip as _gzip

    for name in ("quant.sf", "quant.genes.sf"):
        path = out_dir / name
        gz_path = out_dir / f"{name}.gz"
        if path.exists() and not gz_path.exists():
            with open(path, "rb") as f_in, _gzip.open(gz_path, "wb") as f_out:
                f_out.writelines(f_in)
            path.unlink()
            logger.debug("Compressed %s → %s", path.name, gz_path.name)


def run_salmon(
    cfg: BenchmarkConfig,
    condition: str,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> ToolRunResult:
    """Run salmon quant for a condition.

    Requires ``salmon_index`` to be set in the benchmark config.
    """
    if not cfg.salmon_index:
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=False,
            error="No salmon_index configured",
        )

    cond_dir = cfg.condition_dir(condition)
    out_dir = cond_dir / "salmon"

    if not force and _salmon_output_exists(out_dir):
        logger.info("Skipping %s / salmon — output exists", condition)
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=True, skipped=True,
        )

    r1 = cond_dir / "sim_R1.fq.gz"
    r2 = cond_dir / "sim_R2.fq.gz"
    if not r1.exists() or not r2.exists():
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=False,
            error=f"FASTQ not found: {r1}",
        )

    quant_dir = out_dir / "_quant_tmp"
    cmd = [
        _find_tool("salmon"), "quant",
        "-i", str(cfg.salmon_index),
        "-l", "A",
        "-1", str(r1), "-2", str(r2),
        "-o", str(quant_dir),
        "-p", str(cfg.threads),
        "--validateMappings",
    ]

    cmd_str = shlex.join(cmd)

    if dry_run:
        print(f"[dry-run] {cmd_str}")
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=True, skipped=True,
        )

    print(f"  Running: {condition} / salmon", flush=True)
    logger.info("Command: %s", cmd_str)

    quant_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.monotonic() - t0

        if result.returncode != 0:
            stderr_tail = result.stderr[-500:] if result.stderr else ""
            return ToolRunResult(
                condition=condition, tool="salmon",
                success=False, elapsed=elapsed,
                error=f"exit {result.returncode}: {stderr_tail}",
            )

        # Move quant.sf to final location and compress
        out_dir.mkdir(parents=True, exist_ok=True)
        for name in ("quant.sf", "quant.genes.sf"):
            src = quant_dir / name
            dst = out_dir / name
            if src.exists():
                shutil.move(str(src), str(dst))

        _compress_salmon_output(out_dir)

        # Clean up temp dir
        shutil.rmtree(quant_dir, ignore_errors=True)

        logger.info("%s / salmon completed in %.1fs", condition, elapsed)
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=True, elapsed=elapsed,
        )

    except FileNotFoundError:
        return ToolRunResult(
            condition=condition, tool="salmon",
            success=False,
            error="'salmon' command not found",
        )


def run_kallisto(
    cfg: BenchmarkConfig,
    condition: str,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> ToolRunResult:
    """Run kallisto quant for a condition.

    Requires ``kallisto_index`` to be set in the benchmark config.
    """
    if not cfg.kallisto_index:
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=False,
            error="No kallisto_index configured",
        )

    cond_dir = cfg.condition_dir(condition)
    out_dir = cond_dir / "kallisto"

    if not force and (out_dir / "abundance.tsv").exists():
        logger.info("Skipping %s / kallisto — output exists", condition)
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=True, skipped=True,
        )

    r1 = cond_dir / "sim_R1.fq.gz"
    r2 = cond_dir / "sim_R2.fq.gz"
    if not r1.exists() or not r2.exists():
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=False,
            error=f"FASTQ not found: {r1}",
        )

    # Parse strand specificity from condition name to set --rf-stranded
    from .analysis import parse_condition
    meta = parse_condition(condition)
    ss = meta.get("strand_specificity", 1.0)

    cmd = [
        _find_tool("kallisto"), "quant",
        "-i", str(cfg.kallisto_index),
        "-o", str(out_dir),
        "-t", str(cfg.threads),
    ]
    if ss >= 0.9:
        cmd.append("--rf-stranded")
    cmd.extend([str(r1), str(r2)])

    cmd_str = shlex.join(cmd)

    if dry_run:
        print(f"[dry-run] {cmd_str}")
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=True, skipped=True,
        )

    print(f"  Running: {condition} / kallisto", flush=True)
    logger.info("Command: %s", cmd_str)

    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.monotonic() - t0

        if result.returncode != 0:
            stderr_tail = result.stderr[-500:] if result.stderr else ""
            return ToolRunResult(
                condition=condition, tool="kallisto",
                success=False, elapsed=elapsed,
                error=f"exit {result.returncode}: {stderr_tail}",
            )

        logger.info("%s / kallisto completed in %.1fs", condition, elapsed)
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=True, elapsed=elapsed,
        )

    except FileNotFoundError:
        return ToolRunResult(
            condition=condition, tool="kallisto",
            success=False,
            error="'kallisto' command not found",
        )


def run_all_tools(
    cfg: BenchmarkConfig,
    *,
    tools: list[str] | None = None,
    force: bool = False,
    dry_run: bool = False,
    condition_filter: list[str] | None = None,
) -> ToolRunSummary:
    """Run external tools (salmon, kallisto) for all conditions.

    Parameters
    ----------
    tools : list[str] or None
        Tools to run. Default: auto-detect from config (salmon if
        salmon_index set, kallisto if kallisto_index set).
    """
    summary = ToolRunSummary()
    conditions = condition_filter or cfg.get_conditions()

    if tools is None:
        tools = []
        if cfg.salmon_index:
            tools.append("salmon")
        if cfg.kallisto_index:
            tools.append("kallisto")

    if not tools:
        print("No external tools configured (set salmon_index / kallisto_index).")
        return summary

    print(f"Running {len(tools)} tool(s) × {len(conditions)} condition(s)\n")

    for tool in tools:
        print(f"\n[tool: {tool}]")
        for condition in conditions:
            if tool == "salmon":
                result = run_salmon(cfg, condition, force=force, dry_run=dry_run)
            elif tool == "kallisto":
                result = run_kallisto(cfg, condition, force=force, dry_run=dry_run)
            else:
                result = ToolRunResult(
                    condition=condition, tool=tool,
                    success=False, error=f"Unknown tool: {tool}",
                )
            summary.results.append(result)

    summary.print_summary()
    return summary
