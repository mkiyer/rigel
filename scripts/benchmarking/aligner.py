"""Alignment: run minimap2 (and oracle symlink) for benchmark conditions."""
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
class AlignResult:
    """Result of a single alignment job."""

    condition: str
    aligner: str
    success: bool
    elapsed: float = 0.0
    skipped: bool = False
    error: str = ""
    bam_path: str = ""


@dataclass
class AlignSummary:
    """Aggregate alignment results."""

    results: list[AlignResult] = field(default_factory=list)

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
        print(f"Alignment Summary: {self.n_success} succeeded, "
              f"{self.n_skipped} skipped, {self.n_failed} failed "
              f"(of {self.n_total} total)")
        print(f"{'='*60}")
        if self.n_failed > 0:
            print("\nFailed:")
            for r in self.results:
                if not r.success and not r.skipped:
                    print(f"  {r.condition} / {r.aligner}: {r.error}")
        if self.n_success > 0:
            runs = [r for r in self.results if r.success and not r.skipped]
            total_time = sum(r.elapsed for r in runs)
            print(f"\nTotal alignment time: {total_time:.1f}s "
                  f"(avg {total_time / len(runs):.1f}s per run)")


def _find_tool(name: str) -> str:
    """Resolve tool name to an absolute path."""
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(f"Tool '{name}' not found in PATH")


def _build_minimap2_bed(
    gtf_path: str,
    output_path: Path,
) -> Path:
    """Build a BED12 splice junction file from GTF using paftools.js.

    This file tells minimap2 about annotated splice sites, improving
    alignment of reads near exon boundaries and enabling alignment of
    reads spanning small exons.
    """
    if output_path.exists():
        logger.info("BED12 already exists: %s", output_path)
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Building minimap2 BED12 from GTF: %s", gtf_path)

    try:
        k8 = _find_tool("k8")
        paftools = _find_tool("paftools.js")
        with open(output_path, "w") as f:
            subprocess.run(
                [k8, paftools, "gff2bed", str(gtf_path)],
                stdout=f, stderr=subprocess.PIPE, check=True,
            )
        logger.info("Wrote BED12: %s", output_path)
        return output_path
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        logger.warning("paftools.js gff2bed failed: %s", exc)
        raise RuntimeError(
            "Cannot build minimap2 splice BED: paftools.js unavailable or failed. "
            "Install with: conda install k8"
        ) from exc


def align_minimap2(
    genome_fasta: str,
    fastq_r1: Path,
    fastq_r2: Path,
    out_bam: Path,
    *,
    bed_path: Path | None = None,
    threads: int = 8,
) -> None:
    """Align paired-end short reads with minimap2 → name-sorted BAM.

    Uses ``-ax splice:sr`` for short RNA-seq reads with these key options:

    - ``--secondary=yes -N 20``: report up to 20 secondary alignments per
      read. Critical for genes with many processed pseudogenes (e.g.
      GAPDH, EEF1A1 with 60+ pseudogene copies in the human genome).
    - ``-j FILE``: annotated splice junctions (BED12) to improve alignment
      at annotated exon boundaries and enable tiny-exon recovery.

    Output is name-sorted (``samtools sort -n``) as required by rigel.
    """
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    mm2_cmd = [
        _find_tool("minimap2"),
        "-ax", "splice:sr",
        "--secondary=yes", "-N", "20",
        "-t", str(threads),
    ]
    if bed_path is not None:
        mm2_cmd.extend(["-j", str(bed_path)])
    mm2_cmd.extend([str(genome_fasta), str(fastq_r1), str(fastq_r2)])

    sort_cmd = [
        _find_tool("samtools"), "sort", "-n",
        "-@", str(max(1, threads // 2)),
        "-o", str(out_bam),
    ]

    logger.info("minimap2 | samtools sort -n → %s", out_bam)
    logger.info("minimap2 cmd: %s", shlex.join(mm2_cmd))

    p1 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stderr=subprocess.PIPE)
    p1.stdout.close()
    _, stderr2 = p2.communicate()
    mm2_stderr = p1.stderr.read().decode(errors="replace") if p1.stderr else ""
    p1.wait()

    if mm2_stderr:
        for line in mm2_stderr.strip().split("\n"):
            logger.debug("  minimap2: %s", line)

    if p1.returncode != 0:
        raise RuntimeError(f"minimap2 failed (exit {p1.returncode})")
    if p2.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed (exit {p2.returncode}): "
            f"{stderr2.decode(errors='replace')[-500:]}"
        )


def _oracle_bam_path(cond_dir: Path) -> Path:
    """Return the oracle BAM path for a condition directory."""
    return cond_dir / "sim_oracle.bam"


def run_alignment(
    cfg: BenchmarkConfig,
    condition: str,
    aligner: str,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> AlignResult:
    """Run alignment for a single condition.

    Parameters
    ----------
    aligner : str
        ``"minimap2"`` or ``"oracle"``.
    """
    cond_dir = cfg.condition_dir(condition)
    manifest = cfg.load_manifest()

    genome_fasta = manifest.get("genome", "")
    gtf = manifest.get("gtf", "")

    if aligner == "oracle":
        # Oracle: just check the oracle BAM exists
        oracle_bam = _oracle_bam_path(cond_dir)
        if oracle_bam.exists():
            return AlignResult(
                condition=condition, aligner="oracle",
                success=True, skipped=True,
                bam_path=str(oracle_bam),
            )
        return AlignResult(
            condition=condition, aligner="oracle",
            success=False,
            error=f"Oracle BAM not found: {oracle_bam}",
        )

    elif aligner == "minimap2":
        # Output BAM path
        out_dir = cond_dir / "minimap2"
        out_bam = out_dir / "aligned.bam"

        if not force and out_bam.exists():
            logger.info("Skipping %s / minimap2 — BAM exists", condition)
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=True, skipped=True,
                bam_path=str(out_bam),
            )

        # Find FASTQ files
        r1 = cond_dir / "sim_R1.fq.gz"
        r2 = cond_dir / "sim_R2.fq.gz"
        if not r1.exists() or not r2.exists():
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=False,
                error=f"FASTQ not found: {r1}",
            )

        if not genome_fasta:
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=False,
                error="No genome FASTA specified in manifest",
            )

        # Build BED12 for splice junctions (cached)
        bed_path = None
        if gtf:
            try:
                bed_path = _build_minimap2_bed(
                    gtf, cfg.benchmark_dir / "minimap2_junctions.bed"
                )
            except RuntimeError as exc:
                logger.warning("BED12 build failed, aligning without junctions: %s", exc)

        cmd_str = (
            f"minimap2 -ax splice:sr --secondary=yes -N 20 "
            f"-t {cfg.threads} "
            f"{'-j ' + str(bed_path) + ' ' if bed_path else ''}"
            f"{genome_fasta} {r1} {r2} | samtools sort -n -o {out_bam}"
        )

        if dry_run:
            print(f"[dry-run] {cmd_str}")
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=True, skipped=True,
            )

        print(f"  Aligning: {condition} / minimap2", flush=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        t0 = time.monotonic()
        try:
            align_minimap2(
                genome_fasta, r1, r2, out_bam,
                bed_path=bed_path, threads=cfg.threads,
            )
            elapsed = time.monotonic() - t0
            logger.info("%s / minimap2 completed in %.1fs", condition, elapsed)
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=True, elapsed=elapsed,
                bam_path=str(out_bam),
            )
        except Exception as exc:
            elapsed = time.monotonic() - t0
            logger.error("%s / minimap2 failed: %s", condition, exc)
            return AlignResult(
                condition=condition, aligner="minimap2",
                success=False, elapsed=elapsed,
                error=str(exc),
            )

    else:
        return AlignResult(
            condition=condition, aligner=aligner,
            success=False,
            error=f"Unknown aligner: {aligner}",
        )


def run_all_alignments(
    cfg: BenchmarkConfig,
    *,
    aligners: list[str] | None = None,
    force: bool = False,
    dry_run: bool = False,
    condition_filter: list[str] | None = None,
) -> AlignSummary:
    """Run alignment for all conditions × aligners.

    Parameters
    ----------
    aligners : list[str] or None
        Aligners to run (default: ``["minimap2"]``).
    """
    summary = AlignSummary()

    if aligners is None:
        aligners = ["minimap2"]

    conditions = condition_filter or cfg.get_conditions()

    print(f"Aligning {len(aligners)} aligner(s) × "
          f"{len(conditions)} condition(s)\n")

    for aligner in aligners:
        print(f"\n[aligner: {aligner}]")
        for condition in conditions:
            result = run_alignment(
                cfg, condition, aligner,
                force=force, dry_run=dry_run,
            )
            summary.results.append(result)

    summary.print_summary()
    return summary
