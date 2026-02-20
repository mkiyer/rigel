#!/usr/bin/env python3
"""Regional benchmark: hulkrna (two modes) vs salmon vs kallisto vs htseq.

Supports a combinatorial grid of gDNA contamination levels and strand
specificity values.  Per-region setup (extract, abundance assignment, index
build) is done once; the condition sweep re-simulates reads and re-runs all
tools for each (gDNA level, strand specificity) pair.

Tools
-----
- hulkrna       : include_multimap=False (default mode)
- hulkrna_mm    : include_multimap=True  (multimapping mode)
- salmon        : salmon quant --validateMappings
- kallisto      : kallisto quant (--rf-stranded when strand_spec >= 0.9)
- htseq         : htseq-count (gene-level only, optional, via --include-htseq)

gDNA levels (computed from assigned transcript abundance distribution)
---------------------------------------------------------------------
- none     : 0
- low      : 10th percentile
- moderate : 20th percentile
- high     : 50th percentile (median)

Output per seed directory
─────────────────────────
    summary.json        – full structured results
    summary.csv         – tidy long-form CSV
    summary.md          – human-readable Markdown tables
    diagnostics.md      – aggregate diagnostic notes
    <region>/
        region/         – extracted FASTA & GTF
        transcripts.fa
        hulkrna_index/
        <condition>/          e.g. gdna_none_ss_1.00
            reads/
            align/
            salmon/
            kallisto/
            htseq/
            per_transcript_counts.csv
            per_gene_counts.csv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import os
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from hulkrna.gtf import GTF
from hulkrna.index import HulkIndex, write_bed12
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, ReadSimulator, SimConfig, reverse_complement
from hulkrna.transcript import Transcript
from hulkrna.types import Interval, Strand

logger = logging.getLogger(__name__)

# ── Tool identifiers ────────────────────────────────────────────────

TRANSCRIPT_TOOLS = ("hulkrna", "hulkrna_mm", "salmon", "kallisto")
GENE_TOOLS = ("hulkrna", "hulkrna_mm", "salmon", "kallisto", "htseq")

# ── Default gDNA levels and strand specificities ────────────────────

DEFAULT_GDNA_LEVELS = ("none", "low", "moderate", "high")
DEFAULT_STRAND_SPECIFICITIES = (0.95, 0.99, 1.0)

GDNA_PERCENTILES = {
    "none": 0.0,
    "low": 10.0,
    "moderate": 20.0,
    "high": 50.0,
}


# ── Dataclasses ─────────────────────────────────────────────────────


@dataclass
class RegionSpec:
    label: str
    chrom: str
    start0: int  # 0-based inclusive
    end0: int  # 0-based exclusive


@dataclass
class ToolMetrics:
    tool: str
    elapsed_sec: float
    total_truth: float
    total_observed: float
    total_abs_error: float
    mean_abs_error: float
    rmse: float
    pearson: float
    spearman: float


@dataclass
class ConditionResult:
    """Results for one region under one (gDNA, strand specificity) condition."""

    region: str
    n_transcripts: int
    n_genes: int
    n_fragments: int
    n_gdna_actual: int
    gdna_label: str
    gdna_abundance: float
    strand_specificity: float
    transcript_metrics: dict  # tool name → ToolMetrics (as dict)
    gene_metrics: dict  # tool name → ToolMetrics (as dict)


@dataclass
class StringGenome:
    """Minimal genome interface expected by ReadSimulator."""

    seq: str
    name: str

    def __getitem__(self, key: slice | int) -> str:
        """Support slicing via genome[start:end]."""
        return self.seq[key]

    def __len__(self) -> int:
        return len(self.seq)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        assert chrom == self.name
        return self.seq[start:end]

    def get_reference_length(self, chrom: str) -> int:
        assert chrom == self.name
        return len(self.seq)

    def references(self) -> list[str]:
        return [self.name]

    @property
    def lengths(self) -> dict[str, int]:
        return {self.name: len(self.seq)}


# ── Region parsing ──────────────────────────────────────────────────


def parse_region(region_str: str) -> RegionSpec:
    m = re.match(r"(.+):(\d+)-(\d+)$", region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}")
    chrom = m.group(1)
    start1 = int(m.group(2))
    end1 = int(m.group(3))
    label = f"{chrom}_{start1}_{end1}"
    return RegionSpec(label=label, chrom=chrom, start0=start1 - 1, end0=end1)


def parse_regions(args: argparse.Namespace) -> list[RegionSpec]:
    regions: list[RegionSpec] = []
    for r in args.region:
        regions.append(parse_region(r))
    if args.region_file is not None:
        with open(args.region_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                start1 = int(parts[1])
                end1 = int(parts[2])
                label = parts[3] if len(parts) > 3 else f"{chrom}_{start1}_{end1}"
                regions.append(
                    RegionSpec(label=label, chrom=chrom, start0=start1 - 1, end0=end1)
                )
    if not regions:
        raise ValueError("No regions specified (use --region or --region-file)")
    return regions


# ── GTF / FASTA helpers ─────────────────────────────────────────────


def open_textmaybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def extract_region(
    genome_fa: Path,
    gtf_path: Path,
    region: RegionSpec,
    out_prefix: Path,
) -> tuple[Path, Path, list[Transcript], str]:
    """Extract region FASTA, GTF, transcripts, and raw sequence."""
    out_fa = Path(f"{out_prefix}.fa")
    out_gtf = Path(f"{out_prefix}.gtf")
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Extracting region %s (%s:%d-%d)",
        region.label,
        region.chrom,
        region.start0 + 1,
        region.end0,
    )

    result = subprocess.run(
        [
            "samtools",
            "faidx",
            str(genome_fa),
            f"{region.chrom}:{region.start0 + 1}-{region.end0}",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    # Rewrite with region label as contig name
    fasta_lines = result.stdout.strip().split("\n")
    region_seq = "".join(fasta_lines[1:]).upper()
    with open(out_fa, "w") as fh:
        fh.write(f">{region.label}\n")
        for i in range(0, len(region_seq), 80):
            fh.write(region_seq[i : i + 80] + "\n")
    subprocess.run(
        ["samtools", "faidx", str(out_fa)],
        check=True,
        capture_output=True,
    )

    # Parse GTF transcripts in region
    tx_exons: dict[str, list] = defaultdict(list)
    tx_meta: dict[str, dict] = {}

    with open_textmaybe_gzip(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[0] != region.chrom:
                continue
            feature = parts[2]
            if feature != "exon":
                continue
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            strand_str = parts[6]
            attrs_str = parts[8]

            attrs = {}
            for a in attrs_str.split(";"):
                a = a.strip()
                if not a:
                    continue
                m = re.match(r'(\w+)\s+"([^"]*)"', a)
                if m:
                    attrs[m.group(1)] = m.group(2)

            tid = attrs.get("transcript_id", "")
            gid = attrs.get("gene_id", "")
            gname = attrs.get("gene_name", gid)
            gtype = attrs.get("gene_type", attrs.get("gene_biotype", "unknown"))
            if not tid or not gid:
                continue

            # basic tag filter
            has_basic = "basic" in attrs_str or "tag" not in attrs_str
            if not has_basic:
                tag_val = attrs.get("tag", "")
                has_basic = "basic" in tag_val or "tag" not in attrs_str

            tx_exons[tid].append((start0, end0))
            tx_meta[tid] = {
                "g_id": gid,
                "g_name": gname,
                "g_type": gtype,
                "strand": Strand.from_str(strand_str),
            }

    # Filter to transcripts fully within region
    transcripts: list[Transcript] = []
    t_index = 0
    for tid in sorted(tx_exons):
        exons = tx_exons[tid]
        all_inside = all(s >= region.start0 and e <= region.end0 for s, e in exons)
        if not all_inside:
            continue
        meta = tx_meta[tid]
        rel_exons = sorted(
            [Interval(s - region.start0, e - region.start0) for s, e in exons],
            key=lambda x: x.start,
        )
        t = Transcript(
            t_id=tid,
            t_index=t_index,
            g_id=meta["g_id"],
            g_name=meta["g_name"],
            g_type=meta["g_type"],
            ref=region.label,
            strand=meta["strand"],
            exons=rel_exons,
        )
        transcripts.append(t)
        t_index += 1

    logger.info(
        "Region %s: %d transcripts from %d genes",
        region.label,
        len(transcripts),
        len({t.g_id for t in transcripts}),
    )

    # Write initial GTF (will be rewritten with abundances later)
    with open(out_gtf, "w") as out:
        for t in transcripts:
            for exon in t.exons:
                gtf_obj = GTF(
                    seqname=region.label,
                    source="hulkrna_region_bench",
                    feature="exon",
                    start=exon.start,
                    end=exon.end,
                    score=1.0,
                    strand=Strand.to_str(t.strand),
                    phase=".",
                    attrs={
                        "gene_id": t.g_id,
                        "transcript_id": t.t_id,
                        "gene_name": t.g_name,
                        "gene_type": t.g_type,
                    },
                    tags={"basic"},
                )
                out.write(str(gtf_obj) + "\n")

    return out_fa, out_gtf, transcripts, region_seq


# ── Abundance assignment ────────────────────────────────────────────


def assign_abundances(
    transcripts: list[Transcript],
    *,
    mode: str = "random",
    rng: np.random.Generator | None = None,
    min_abund: float = 1.0,
    max_abund: float = 1000.0,
    abundance_file: Path | None = None,
) -> None:
    """Assign expression abundances to transcripts.

    Modes:
      random  – log-uniform in [min_abund, max_abund], rescaled so median ≈ 100
      uniform – all transcripts get the same abundance (100)
      file    – load from CSV/TSV with columns transcript_id, abundance
    """
    if mode == "uniform":
        for t in transcripts:
            t.abundance = 100.0
        return

    if mode == "random":
        if rng is None:
            rng = np.random.default_rng(42)
        log_ab = rng.uniform(
            low=np.log(max(min_abund, 1e-9)),
            high=np.log(max(max_abund, 1.0)),
            size=len(transcripts),
        )
        raw = np.exp(log_ab)
        # Rescale so the median ≈ 100
        med = float(np.median(raw))
        if med > 0:
            raw *= 100.0 / med
        for t, ab in zip(transcripts, raw):
            t.abundance = max(1e-9, float(ab))
        return

    if mode == "file":
        if abundance_file is None:
            raise ValueError("abundance_file required for mode='file'")
        sep = "\t" if str(abundance_file).endswith(".tsv") else ","
        df = pd.read_csv(abundance_file, sep=sep)
        if "transcript_id" not in df.columns or "abundance" not in df.columns:
            raise ValueError(
                "Abundance file must contain columns: transcript_id, abundance"
            )
        abund_map = {
            str(tid): float(ab)
            for tid, ab in zip(df["transcript_id"], df["abundance"], strict=False)
        }
        missing = 0
        for t in transcripts:
            if t.t_id in abund_map:
                t.abundance = max(1e-9, float(abund_map[t.t_id]))
            else:
                missing += 1
                t.abundance = 1.0
        logger.info(
            "Loaded abundances from %s: matched=%d missing=%d",
            abundance_file,
            len(transcripts) - missing,
            missing,
        )
        return

    raise ValueError(f"Unsupported abundance mode: {mode}")


# ── Sequence helpers ────────────────────────────────────────────────


def write_transcript_fasta(
    transcripts: list[Transcript], region_seq: str, out_fa: Path
) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fa, "w") as out:
        for t in transcripts:
            parts = [region_seq[e.start : e.end] for e in t.exons]
            seq = "".join(parts)
            if t.strand == Strand.NEG:
                seq = reverse_complement(seq)
            out.write(f">{t.t_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


# ── Alignment ───────────────────────────────────────────────────────


def align_minimap2(
    region_fa: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_bam: Path,
    bed_path: Path | None = None,
) -> None:
    """Align PE reads with minimap2 and produce a name-sorted BAM."""
    minimap2_cmd = [
        "minimap2",
        "-ax",
        "splice:sr",
        "--secondary=yes",
    ]
    if bed_path is not None:
        minimap2_cmd.extend(["-j", str(bed_path)])
    minimap2_cmd.extend([
        str(region_fa),
        str(fastq_r1),
        str(fastq_r2),
    ])
    sort_cmd = ["samtools", "sort", "-n", "-o", str(out_bam), "-"]

    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(
        sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    p1.stdout.close()

    _, p2_stderr = p2.communicate()
    p1_stderr = p1.stderr.read()
    p1.wait()

    if p2.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {p2_stderr.decode()}")
    if p1.returncode != 0:
        logger.warning(
            "minimap2 returned rc=%d: %s", p1.returncode, p1_stderr.decode()
        )


def coord_sort_bam(input_bam: Path, output_bam: Path) -> None:
    """Re-sort a BAM file to coordinate order and index it."""
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["samtools", "sort", "-o", str(output_bam), str(input_bam)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["samtools", "index", str(output_bam)],
        check=True,
        capture_output=True,
    )


# ── Truth parsing ───────────────────────────────────────────────────


def parse_truth_from_fastq(r1_fastq: Path) -> tuple[dict[str, int], int]:
    counts: Counter[str] = Counter()
    n_gdna = 0
    with open(r1_fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            qname = line[1:].strip()
            if qname.endswith("/1"):
                qname = qname[:-2]
            tid = qname.split(":", 1)[0]
            if tid == "gdna":
                n_gdna += 1
            else:
                counts[tid] += 1
    return dict(counts), n_gdna


# ── Tool runners ────────────────────────────────────────────────────


def run_hulkrna_tool(
    bam_path: Path,
    index: HulkIndex,
    include_multimap: bool,
    args: argparse.Namespace,
) -> tuple[dict[str, float], float]:
    """Run hulkrna pipeline and return (transcript_counts, elapsed_sec)."""
    t0 = time.monotonic()
    pipe = run_pipeline(
        bam_path,
        index,
        seed=args.pipeline_seed,
        sj_strand_tag="auto",
        include_multimap=include_multimap,
        gdna_threshold=args.gdna_threshold,
        gdna_splice_penalty_unannot=args.gdna_splice_penalty_unannot,
    )
    elapsed = time.monotonic() - t0

    counts_df = pipe.counter.get_counts_df(index)
    counts = {
        row.transcript_id: float(row.count)
        for row in counts_df.itertuples(index=False)
    }
    return counts, elapsed


def run_salmon(
    transcript_fa: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_dir: Path,
    threads: int,
) -> tuple[dict[str, float], float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    idx = out_dir / "index"
    quant_dir = out_dir / "quant"

    t0 = time.monotonic()
    subprocess.run(
        ["salmon", "index", "-t", str(transcript_fa), "-i", str(idx), "-k", "23"],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [
            "salmon",
            "quant",
            "-i",
            str(idx),
            "-l",
            "A",
            "-1",
            str(fastq_r1),
            "-2",
            str(fastq_r2),
            "-o",
            str(quant_dir),
            "-p",
            str(threads),
            "--validateMappings",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    elapsed = time.monotonic() - t0

    quant_sf = quant_dir / "quant.sf"
    df = pd.read_csv(quant_sf, sep="\t")
    return dict(zip(df["Name"], df["NumReads"])), elapsed


def run_kallisto(
    transcript_fa: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_dir: Path,
    threads: int,
    strand_specificity: float,
) -> tuple[dict[str, float], float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    idx = out_dir / "index.kallisto"
    quant_dir = out_dir / "quant"

    t0 = time.monotonic()
    subprocess.run(
        ["kallisto", "index", "-i", str(idx), str(transcript_fa)],
        check=True,
        capture_output=True,
        text=True,
    )

    cmd = [
        "kallisto",
        "quant",
        "-i",
        str(idx),
        "-o",
        str(quant_dir),
        "-t",
        str(threads),
    ]
    # Simulator uses dUTP-like orientation (R1 reverse of transcript) when stranded.
    if strand_specificity >= 0.9:
        cmd.append("--rf-stranded")
    cmd.extend([str(fastq_r1), str(fastq_r2)])

    subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.monotonic() - t0

    abund_tsv = quant_dir / "abundance.tsv"
    df = pd.read_csv(abund_tsv, sep="\t")
    return dict(zip(df["target_id"], df["est_counts"])), elapsed


def run_htseq(
    bam_path: Path,
    gtf_path: Path,
    out_dir: Path,
    strand_specificity: float,
    conda_env: str = "htseq",
) -> tuple[dict[str, int], float]:
    """Run htseq-count on a coordinate-sorted BAM.

    Returns (gene_id → count, elapsed_sec).  Special counters
    (__no_feature, __ambiguous, etc.) are excluded.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "htseq_counts.txt"

    stranded = "reverse" if strand_specificity >= 0.9 else "no"

    cmd = [
        "conda",
        "run",
        "-n",
        conda_env,
        "htseq-count",
        "-f",
        "bam",
        "-r",
        "pos",
        f"--stranded={stranded}",
        "-t",
        "exon",
        "-i",
        "gene_id",
        str(bam_path),
        str(gtf_path),
    ]

    t0 = time.monotonic()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    elapsed = time.monotonic() - t0

    # Save raw output
    with open(out_file, "w") as f:
        f.write(result.stdout)

    counts: dict[str, int] = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        gene_id = parts[0]
        if gene_id.startswith("__"):
            continue
        counts[gene_id] = int(parts[1])

    return counts, elapsed


# ── Scoring ─────────────────────────────────────────────────────────


def score_tool(
    tool: str,
    truth: dict[str, int | float],
    observed: dict[str, float],
    elapsed_sec: float,
) -> ToolMetrics:
    """Score a tool's predictions against ground truth."""
    tids = sorted(truth)
    truth_arr = np.array([truth.get(t, 0) for t in tids], dtype=float)
    obs_arr = np.array([float(observed.get(t, 0.0)) for t in tids], dtype=float)

    diff = obs_arr - truth_arr
    abs_diff = np.abs(diff)
    total_abs_error = float(abs_diff.sum())
    mean_abs_error = float(abs_diff.mean()) if len(abs_diff) else 0.0
    rmse = float(np.sqrt(np.mean(diff * diff))) if len(diff) else 0.0

    if len(tids) > 1:
        pearson = (
            float(np.corrcoef(truth_arr, obs_arr)[0, 1])
            if np.std(truth_arr) > 0 and np.std(obs_arr) > 0
            else 0.0
        )
        s_truth = pd.Series(truth_arr)
        s_obs = pd.Series(obs_arr)
        spearman = float(s_truth.corr(s_obs, method="spearman"))
        if math.isnan(spearman):
            spearman = 0.0
    else:
        pearson = 0.0
        spearman = 0.0

    return ToolMetrics(
        tool=tool,
        elapsed_sec=round(elapsed_sec, 3),
        total_truth=float(truth_arr.sum()),
        total_observed=float(obs_arr.sum()),
        total_abs_error=round(total_abs_error, 3),
        mean_abs_error=round(mean_abs_error, 3),
        rmse=round(rmse, 3),
        pearson=round(pearson, 5),
        spearman=round(spearman, 5),
    )


# ── Gene-level aggregation ──────────────────────────────────────────


def aggregate_to_genes(
    tx_counts: dict[str, float | int],
    gene_map: dict[str, str],
) -> dict[str, float]:
    """Sum transcript-level counts to gene level.

    Parameters
    ----------
    tx_counts : dict
        transcript_id → count
    gene_map : dict
        transcript_id → gene_id
    """
    gene_counts: dict[str, float] = defaultdict(float)
    for tid, count in tx_counts.items():
        gid = gene_map.get(tid)
        if gid is not None:
            gene_counts[gid] += float(count)
    return dict(gene_counts)


# ── gDNA level computation ──────────────────────────────────────────


def compute_gdna_levels(
    transcripts: list[Transcript],
    requested_levels: tuple[str, ...] = DEFAULT_GDNA_LEVELS,
) -> dict[str, float]:
    """Compute gDNA abundance values from transcript abundance percentiles.

    Returns a dict mapping level label to actual abundance value.
    """
    abundances = np.array([t.abundance for t in transcripts])
    levels: dict[str, float] = {}
    for label in requested_levels:
        pctl = GDNA_PERCENTILES.get(label)
        if pctl is None:
            raise ValueError(
                f"Unknown gDNA level '{label}'. "
                f"Valid levels: {', '.join(GDNA_PERCENTILES)}"
            )
        if pctl == 0.0:
            levels[label] = 0.0
        else:
            levels[label] = float(np.percentile(abundances, pctl))
    return levels


# ── Condition directory naming ───────────────────────────────────────


def condition_dir_name(gdna_label: str, strand_specificity: float) -> str:
    return f"gdna_{gdna_label}_ss_{strand_specificity:.2f}"


# ── Main benchmark orchestration ─────────────────────────────────────


def run_region_benchmark(
    region: RegionSpec,
    args: argparse.Namespace,
    out_root: Path,
    gdna_levels: dict[str, float],
    strand_specificities: tuple[float, ...],
    include_htseq: bool,
    htseq_conda_env: str,
) -> list[ConditionResult]:
    """Run benchmark for one region across all conditions.

    Returns a list of ConditionResult (one per condition combination).
    """
    reg_dir = out_root / region.label
    reg_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[REGION] %s:%d-%d", region.chrom, region.start0 + 1, region.end0)

    # ── Per-region setup (done once) ────────────────────────────────

    region_fa, region_gtf, transcripts, region_seq = extract_region(
        args.genome,
        args.gtf,
        region,
        reg_dir / "region",
    )
    if not transcripts:
        logger.warning(
            "No fully-contained transcripts in region %s; skipping",
            region.label,
        )
        return []

    rng = np.random.default_rng(
        args.abundance_seed + abs(hash(region.label)) % 10_000
    )
    assign_abundances(
        transcripts,
        mode=args.abundance_mode,
        rng=rng,
        min_abund=args.abundance_min,
        max_abund=args.abundance_max,
        abundance_file=args.abundance_file,
    )

    # Rewrite region GTF with abundance in score field
    with open(region_gtf, "w") as out:
        for t in transcripts:
            for exon in t.exons:
                gtf_obj = GTF(
                    seqname=region.label,
                    source="hulkrna_region_bench",
                    feature="exon",
                    start=exon.start,
                    end=exon.end,
                    score=t.abundance,
                    strand=Strand.to_str(t.strand),
                    phase=".",
                    attrs={
                        "gene_id": t.g_id,
                        "transcript_id": t.t_id,
                        "gene_name": t.g_name,
                        "gene_type": t.g_type,
                    },
                    tags={"basic"},
                )
                out.write(str(gtf_obj) + "\n")

    transcript_fa = reg_dir / "transcripts.fa"
    write_transcript_fasta(transcripts, region_seq, transcript_fa)

    # BED12 for minimap2 annotation-guided alignment (reused across conditions)
    bed_path = reg_dir / "annotation.bed"
    write_bed12(transcripts, bed_path)

    # hulkrna index (reused across conditions)
    index_dir = reg_dir / "hulkrna_index"
    HulkIndex.build(region_fa, region_gtf, index_dir, write_tsv=False)
    index = HulkIndex.load(index_dir)

    # Gene mapping for gene-level aggregation
    gene_map = {t.t_id: t.g_id for t in transcripts}
    abundance_map = {t.t_id: t.abundance for t in transcripts}
    gene_abundance = defaultdict(float)
    for t in transcripts:
        gene_abundance[t.g_id] += t.abundance

    n_transcripts = len(transcripts)
    n_genes = len({t.g_id for t in transcripts})

    # ── Condition sweep ─────────────────────────────────────────────

    results: list[ConditionResult] = []

    for gdna_label, gdna_abundance in gdna_levels.items():
        for strand_spec in strand_specificities:
            cond_name = condition_dir_name(gdna_label, strand_spec)
            cond_dir = reg_dir / cond_name
            cond_dir.mkdir(parents=True, exist_ok=True)

            logger.info(
                "  [COND] %s | gDNA=%s (%.2f) | strand_spec=%.2f",
                region.label,
                gdna_label,
                gdna_abundance,
                strand_spec,
            )

            # Simulate reads
            sim_cfg = SimConfig(
                frag_mean=args.frag_mean,
                frag_std=args.frag_std,
                frag_min=args.frag_min,
                frag_max=args.frag_max,
                read_length=args.read_length,
                error_rate=args.error_rate,
                strand_specificity=strand_spec,
                seed=args.sim_seed,
            )
            gdna_cfg = None
            if gdna_abundance > 0:
                gdna_cfg = GDNAConfig(
                    abundance=gdna_abundance,
                    frag_mean=args.gdna_frag_mean,
                    frag_std=args.gdna_frag_std,
                    frag_min=args.gdna_frag_min,
                    frag_max=args.gdna_frag_max,
                )

            sim_genome = StringGenome(region_seq, region.label)
            simulator = ReadSimulator(
                sim_genome, transcripts, config=sim_cfg, gdna_config=gdna_cfg
            )
            fastq_dir = cond_dir / "reads"
            fastq_r1, fastq_r2 = simulator.write_fastq(
                fastq_dir, args.n_fragments, prefix="sim"
            )
            truth_counts, n_gdna_actual = parse_truth_from_fastq(fastq_r1)

            # Align (name-sorted for hulkrna)
            bam_ns = cond_dir / "align" / "reads_namesort.bam"
            bam_ns.parent.mkdir(parents=True, exist_ok=True)
            align_minimap2(region_fa, fastq_r1, fastq_r2, bam_ns, bed_path=bed_path)

            # ── Run transcript-level tools ──────────────────────────

            hulkrna_counts, hulkrna_elapsed = run_hulkrna_tool(
                bam_ns, index, include_multimap=False, args=args
            )
            hulkrna_mm_counts, hulkrna_mm_elapsed = run_hulkrna_tool(
                bam_ns, index, include_multimap=True, args=args
            )
            salmon_counts, salmon_elapsed = run_salmon(
                transcript_fa,
                fastq_r1,
                fastq_r2,
                cond_dir / "salmon",
                threads=args.threads,
            )
            kallisto_counts, kallisto_elapsed = run_kallisto(
                transcript_fa,
                fastq_r1,
                fastq_r2,
                cond_dir / "kallisto",
                threads=args.threads,
                strand_specificity=strand_spec,
            )

            tx_tool_counts = {
                "hulkrna": (hulkrna_counts, hulkrna_elapsed),
                "hulkrna_mm": (hulkrna_mm_counts, hulkrna_mm_elapsed),
                "salmon": (salmon_counts, salmon_elapsed),
                "kallisto": (kallisto_counts, kallisto_elapsed),
            }

            # ── Run htseq (gene-level only) ─────────────────────────

            htseq_gene_counts: dict[str, int] = {}
            htseq_elapsed = 0.0
            if include_htseq:
                bam_cs = cond_dir / "align" / "reads_coordsort.bam"
                coord_sort_bam(bam_ns, bam_cs)
                htseq_gene_counts, htseq_elapsed = run_htseq(
                    bam_cs,
                    region_gtf,
                    cond_dir / "htseq",
                    strand_specificity=strand_spec,
                    conda_env=htseq_conda_env,
                )

            # ── Score transcript-level ──────────────────────────────

            transcript_metrics = {}
            for tool, (counts, elapsed) in tx_tool_counts.items():
                transcript_metrics[tool] = asdict(
                    score_tool(tool, truth_counts, counts, elapsed)
                )

            # ── Score gene-level ────────────────────────────────────

            truth_gene = aggregate_to_genes(truth_counts, gene_map)
            gene_metrics = {}
            for tool, (counts, elapsed) in tx_tool_counts.items():
                gene_counts = aggregate_to_genes(counts, gene_map)
                gene_metrics[tool] = asdict(
                    score_tool(tool, truth_gene, gene_counts, elapsed)
                )
            if include_htseq:
                gene_metrics["htseq"] = asdict(
                    score_tool("htseq", truth_gene, htseq_gene_counts, htseq_elapsed)
                )

            # ── Write per-transcript CSV ────────────────────────────

            per_tx_csv = cond_dir / "per_transcript_counts.csv"
            tids = sorted(truth_counts)
            with open(per_tx_csv, "w", newline="") as f:
                writer = csv.writer(f)
                header = [
                    "transcript_id",
                    "gene_id",
                    "truth",
                    "hulkrna",
                    "hulkrna_mm",
                    "salmon",
                    "kallisto",
                    "abundance",
                ]
                writer.writerow(header)
                for tid in tids:
                    writer.writerow([
                        tid,
                        gene_map.get(tid, ""),
                        int(truth_counts.get(tid, 0)),
                        float(hulkrna_counts.get(tid, 0.0)),
                        float(hulkrna_mm_counts.get(tid, 0.0)),
                        float(salmon_counts.get(tid, 0.0)),
                        float(kallisto_counts.get(tid, 0.0)),
                        float(abundance_map.get(tid, 0.0)),
                    ])

            # ── Write per-gene CSV ──────────────────────────────────

            per_gene_csv = cond_dir / "per_gene_counts.csv"
            gene_ids = sorted(truth_gene)
            with open(per_gene_csv, "w", newline="") as f:
                writer = csv.writer(f)
                header = [
                    "gene_id",
                    "truth",
                    "hulkrna",
                    "hulkrna_mm",
                    "salmon",
                    "kallisto",
                ]
                if include_htseq:
                    header.append("htseq")
                header.append("abundance")
                writer.writerow(header)
                for gid in gene_ids:
                    row = [
                        gid,
                        int(truth_gene.get(gid, 0)),
                        float(aggregate_to_genes(hulkrna_counts, gene_map).get(gid, 0.0)),
                        float(aggregate_to_genes(hulkrna_mm_counts, gene_map).get(gid, 0.0)),
                        float(aggregate_to_genes(salmon_counts, gene_map).get(gid, 0.0)),
                        float(aggregate_to_genes(kallisto_counts, gene_map).get(gid, 0.0)),
                    ]
                    if include_htseq:
                        row.append(int(htseq_gene_counts.get(gid, 0)))
                    row.append(float(gene_abundance.get(gid, 0.0)))
                    writer.writerow(row)

            results.append(
                ConditionResult(
                    region=region.label,
                    n_transcripts=n_transcripts,
                    n_genes=n_genes,
                    n_fragments=args.n_fragments,
                    n_gdna_actual=n_gdna_actual,
                    gdna_label=gdna_label,
                    gdna_abundance=round(gdna_abundance, 4),
                    strand_specificity=strand_spec,
                    transcript_metrics=transcript_metrics,
                    gene_metrics=gene_metrics,
                )
            )

            logger.info(
                "    hulkrna MAE=%.2f | hulkrna_mm MAE=%.2f | "
                "salmon MAE=%.2f | kallisto MAE=%.2f%s",
                transcript_metrics["hulkrna"]["mean_abs_error"],
                transcript_metrics["hulkrna_mm"]["mean_abs_error"],
                transcript_metrics["salmon"]["mean_abs_error"],
                transcript_metrics["kallisto"]["mean_abs_error"],
                f" | htseq gene-MAE={gene_metrics['htseq']['mean_abs_error']:.2f}"
                if include_htseq
                else "",
            )

    return results


# ── Output functions ─────────────────────────────────────────────────


def write_summary(results: list[ConditionResult], out_root: Path, include_htseq: bool) -> None:
    """Write summary.json, summary.csv, and summary.md."""
    out_json = out_root / "summary.json"
    out_csv = out_root / "summary.csv"
    out_md = out_root / "summary.md"

    # ── JSON ────────────────────────────────────────────────────────

    with open(out_json, "w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)

    # ── CSV (tidy long-form) ────────────────────────────────────────

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "region",
            "n_transcripts",
            "n_genes",
            "n_fragments",
            "n_gdna_actual",
            "gdna_label",
            "gdna_abundance",
            "strand_specificity",
            "level",
            "tool",
            "elapsed_sec",
            "total_abs_error",
            "mean_abs_error",
            "rmse",
            "pearson",
            "spearman",
        ])
        for r in results:
            base = [
                r.region,
                r.n_transcripts,
                r.n_genes,
                r.n_fragments,
                r.n_gdna_actual,
                r.gdna_label,
                r.gdna_abundance,
                r.strand_specificity,
            ]
            for tool, tm in r.transcript_metrics.items():
                writer.writerow(
                    base
                    + [
                        "transcript",
                        tm["tool"],
                        tm["elapsed_sec"],
                        tm["total_abs_error"],
                        tm["mean_abs_error"],
                        tm["rmse"],
                        tm["pearson"],
                        tm["spearman"],
                    ]
                )
            for tool, tm in r.gene_metrics.items():
                writer.writerow(
                    base
                    + [
                        "gene",
                        tm["tool"],
                        tm["elapsed_sec"],
                        tm["total_abs_error"],
                        tm["mean_abs_error"],
                        tm["rmse"],
                        tm["pearson"],
                        tm["spearman"],
                    ]
                )

    # ── Markdown ────────────────────────────────────────────────────

    tx_tools = list(TRANSCRIPT_TOOLS)
    gene_tools = list(GENE_TOOLS) if include_htseq else list(TRANSCRIPT_TOOLS)

    lines = [
        "# Regional Benchmark Summary",
        "",
        "## Transcript-Level Metrics",
        "",
        "| Region | gDNA | SS | Tool | MAE | RMSE | Pearson | Spearman | Time (s) |",
        "| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in results:
        for tool in tx_tools:
            tm = r.transcript_metrics.get(tool)
            if tm is None:
                continue
            lines.append(
                f"| {r.region} | {r.gdna_label} | {r.strand_specificity:.2f} | "
                f"{tm['tool']} | {tm['mean_abs_error']:.2f} | {tm['rmse']:.2f} | "
                f"{tm['pearson']:.4f} | {tm['spearman']:.4f} | {tm['elapsed_sec']:.2f} |"
            )

    lines.extend([
        "",
        "## Gene-Level Metrics",
        "",
        "| Region | gDNA | SS | Tool | MAE | RMSE | Pearson | Spearman | Time (s) |",
        "| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for r in results:
        for tool in gene_tools:
            tm = r.gene_metrics.get(tool)
            if tm is None:
                continue
            lines.append(
                f"| {r.region} | {r.gdna_label} | {r.strand_specificity:.2f} | "
                f"{tm['tool']} | {tm['mean_abs_error']:.2f} | {tm['rmse']:.2f} | "
                f"{tm['pearson']:.4f} | {tm['spearman']:.4f} | {tm['elapsed_sec']:.2f} |"
            )

    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diagnostics(
    results: list[ConditionResult], out_root: Path, include_htseq: bool
) -> None:
    """Write aggregate diagnostic report."""
    tx_tools = list(TRANSCRIPT_TOOLS)

    # ── Aggregate transcript-level metrics by tool ──────────────────

    rows = []
    for r in results:
        for tool in tx_tools:
            tm = r.transcript_metrics.get(tool)
            if tm is None:
                continue
            rows.append({
                "region": r.region,
                "gdna_label": r.gdna_label,
                "strand_specificity": r.strand_specificity,
                "tool": tool,
                "mae": tm["mean_abs_error"],
                "rmse": tm["rmse"],
                "pearson": tm["pearson"],
                "spearman": tm["spearman"],
            })
    df = pd.DataFrame(rows)

    lines = [
        "# Diagnostic Notes",
        "",
        "## Overall transcript-level metrics (mean across all regions & conditions)",
        "",
        "| Tool | MAE | RMSE | Pearson | Spearman |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    if len(df):
        agg = (
            df.groupby("tool")[["mae", "rmse", "pearson", "spearman"]]
            .mean()
            .sort_values("mae")
        )
        for tool, row in agg.iterrows():
            lines.append(
                f"| {tool} | {row['mae']:.3f} | {row['rmse']:.3f} | "
                f"{row['pearson']:.4f} | {row['spearman']:.4f} |"
            )

    # ── By gDNA level ───────────────────────────────────────────────

    lines.extend([
        "",
        "## Transcript-level MAE by gDNA level",
        "",
    ])
    if len(df):
        header = "| gDNA Level | " + " | ".join(tx_tools) + " |"
        sep = "| --- | " + " | ".join(["---:"] * len(tx_tools)) + " |"
        lines.extend([header, sep])
        for gdna_label in DEFAULT_GDNA_LEVELS:
            sub = df[df["gdna_label"] == gdna_label]
            if len(sub) == 0:
                continue
            vals = []
            for tool in tx_tools:
                tsub = sub[sub["tool"] == tool]
                vals.append(f"{tsub['mae'].mean():.3f}" if len(tsub) else "—")
            lines.append(f"| {gdna_label} | " + " | ".join(vals) + " |")

    # ── By strand specificity ───────────────────────────────────────

    lines.extend([
        "",
        "## Transcript-level MAE by strand specificity",
        "",
    ])
    if len(df):
        header = "| Strand Spec | " + " | ".join(tx_tools) + " |"
        sep = "| --- | " + " | ".join(["---:"] * len(tx_tools)) + " |"
        lines.extend([header, sep])
        for ss in sorted(df["strand_specificity"].unique()):
            sub = df[df["strand_specificity"] == ss]
            vals = []
            for tool in tx_tools:
                tsub = sub[sub["tool"] == tool]
                vals.append(f"{tsub['mae'].mean():.3f}" if len(tsub) else "—")
            lines.append(f"| {ss:.2f} | " + " | ".join(vals) + " |")

    # ── Per-transcript error by abundance quartile ──────────────────

    per_tx_parts: list[pd.DataFrame] = []
    for r in results:
        cond_name = condition_dir_name(r.gdna_label, r.strand_specificity)
        f = out_root / r.region / cond_name / "per_transcript_counts.csv"
        if f.exists():
            part = pd.read_csv(f)
            part["region"] = r.region
            part["gdna_label"] = r.gdna_label
            part["strand_specificity"] = r.strand_specificity
            per_tx_parts.append(part)

    if per_tx_parts:
        alltx = pd.concat(per_tx_parts, ignore_index=True)
        for tool in tx_tools:
            if tool in alltx.columns:
                alltx[f"ae_{tool}"] = (alltx[tool] - alltx["truth"]).abs()
        alltx["ab_q"] = pd.qcut(
            alltx["abundance"].rank(method="first"),
            4,
            labels=["Q1_low", "Q2", "Q3", "Q4_high"],
        )

        lines.extend([
            "",
            "## Mean absolute error by abundance quartile (all conditions pooled)",
            "",
            "| Quartile | " + " | ".join(tx_tools) + " |",
            "| --- | " + " | ".join(["---:"] * len(tx_tools)) + " |",
        ])
        for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
            sub = alltx[alltx["ab_q"] == q]
            if len(sub) == 0:
                continue
            vals = []
            for tool in tx_tools:
                col = f"ae_{tool}"
                vals.append(f"{sub[col].mean():.3f}" if col in sub.columns else "—")
            lines.append(f"| {q} | " + " | ".join(vals) + " |")

        lines.extend([
            "",
            "## Dropout rate (truth > 0 but predicted ≤ 0)",
            "",
            "| Tool | Dropout rate |",
            "| --- | ---: |",
        ])
        for tool in tx_tools:
            if tool in alltx.columns:
                rate = float(
                    ((alltx["truth"] > 0) & (alltx[tool] <= 0)).mean()
                )
                lines.append(f"| {tool} | {rate:.4f} |")

    out = out_root / "diagnostics.md"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")


# ── CLI ──────────────────────────────────────────────────────────────


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Regional benchmark: hulkrna vs salmon vs kallisto vs htseq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--genome", type=Path, required=True, help="Genome FASTA path")
    p.add_argument(
        "--gtf", type=Path, required=True, help="Gene annotation GTF/GTF.GZ path"
    )
    p.add_argument(
        "--region",
        action="append",
        default=[],
        help="Region: chr:start-end (1-based)",
    )
    p.add_argument(
        "--region-file",
        type=Path,
        default=None,
        help="TSV: chrom<TAB>start1<TAB>end1<TAB>optional_label",
    )
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")

    p.add_argument("--n-fragments", type=int, default=20000)
    p.add_argument("--threads", type=int, default=1)

    p.add_argument("--sim-seed", type=int, default=42)
    p.add_argument("--pipeline-seed", type=int, default=42)

    # Fragment parameters
    p.add_argument("--frag-mean", type=float, default=250.0)
    p.add_argument("--frag-std", type=float, default=50.0)
    p.add_argument("--frag-min", type=int, default=50)
    p.add_argument("--frag-max", type=int, default=1000)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--error-rate", type=float, default=0.0)

    # gDNA contamination parameters
    p.add_argument(
        "--gdna-levels",
        type=str,
        default="none,low,moderate,high",
        help="Comma-separated gDNA levels: none,low,moderate,high (default: all four)",
    )
    p.add_argument("--gdna-frag-mean", type=float, default=350.0)
    p.add_argument("--gdna-frag-std", type=float, default=100.0)
    p.add_argument("--gdna-frag-min", type=int, default=100)
    p.add_argument("--gdna-frag-max", type=int, default=1000)
    p.add_argument("--gdna-threshold", type=float, default=0.5)
    p.add_argument("--gdna-splice-penalty-unannot", type=float, default=0.01)

    # Strand specificity sweep
    p.add_argument(
        "--strand-specificities",
        type=str,
        default="0.95,0.99,1.0",
        help="Comma-separated strand specificity values (default: 0.95,0.99,1.0)",
    )

    # Abundance parameters
    p.add_argument(
        "--abundance-mode",
        choices=["random", "uniform", "file"],
        default="random",
    )
    p.add_argument("--abundance-seed", type=int, default=123)
    p.add_argument("--abundance-min", type=float, default=1.0)
    p.add_argument("--abundance-max", type=float, default=1000.0)
    p.add_argument(
        "--abundance-file",
        type=Path,
        default=None,
        help="CSV/TSV with columns: transcript_id, abundance",
    )

    # HTSeq
    p.add_argument(
        "--include-htseq",
        action="store_true",
        help="Include htseq-count (gene-level) in the benchmark",
    )
    p.add_argument(
        "--htseq-conda-env",
        type=str,
        default="htseq",
        help="Conda environment name where htseq-count is installed (default: htseq)",
    )

    # Runtime
    p.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue with next region on failure",
    )
    p.add_argument("--verbose", action="store_true")

    return p


def ensure_tools(include_htseq: bool = False, htseq_conda_env: str = "htseq") -> None:
    required = ["minimap2", "samtools", "salmon", "kallisto"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(f"Missing required tools in PATH: {', '.join(missing)}")

    if include_htseq:
        # Check htseq-count is available in the specified conda env
        try:
            subprocess.run(
                ["conda", "run", "-n", htseq_conda_env, "htseq-count", "--help"],
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise RuntimeError(
                f"htseq-count not available in conda env '{htseq_conda_env}'. "
                f"Install with: conda create -n {htseq_conda_env} -c bioconda htseq"
            ) from e


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    include_htseq = args.include_htseq
    htseq_conda_env = args.htseq_conda_env

    ensure_tools(include_htseq=include_htseq, htseq_conda_env=htseq_conda_env)

    if args.abundance_mode == "file" and args.abundance_file is None:
        raise ValueError("--abundance-file is required when --abundance-mode=file")

    regions = parse_regions(args)

    # Parse gDNA levels and strand specificities from CLI
    gdna_level_names = tuple(args.gdna_levels.split(","))
    strand_specificities = tuple(float(s) for s in args.strand_specificities.split(","))

    logger.info(
        "Benchmark grid: %d gDNA levels × %d strand specs × %d regions = %d conditions",
        len(gdna_level_names),
        len(strand_specificities),
        len(regions),
        len(gdna_level_names) * len(strand_specificities) * len(regions),
    )

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_results: list[ConditionResult] = []

    for region in regions:
        try:
            # Extract & assign abundances to compute gDNA percentiles.
            # We do a lightweight pre-pass for gDNA level computation,
            # then the full benchmark reuses those abundances.
            # Actually, run_region_benchmark does everything internally,
            # but we need gDNA levels before calling it.

            # Pre-compute gDNA levels for this region
            reg_dir = args.outdir / region.label
            reg_dir.mkdir(parents=True, exist_ok=True)

            _, _, transcripts_pre, _ = extract_region(
                args.genome,
                args.gtf,
                region,
                reg_dir / "region",
            )
            if not transcripts_pre:
                logger.warning(
                    "No fully-contained transcripts in region %s; skipping",
                    region.label,
                )
                continue

            rng_pre = np.random.default_rng(
                args.abundance_seed + abs(hash(region.label)) % 10_000
            )
            assign_abundances(
                transcripts_pre,
                mode=args.abundance_mode,
                rng=rng_pre,
                min_abund=args.abundance_min,
                max_abund=args.abundance_max,
                abundance_file=args.abundance_file,
            )
            gdna_levels = compute_gdna_levels(transcripts_pre, gdna_level_names)
            logger.info(
                "  gDNA levels for %s: %s",
                region.label,
                {k: f"{v:.2f}" for k, v in gdna_levels.items()},
            )

            region_results = run_region_benchmark(
                region,
                args,
                args.outdir,
                gdna_levels=gdna_levels,
                strand_specificities=strand_specificities,
                include_htseq=include_htseq,
                htseq_conda_env=htseq_conda_env,
            )
            all_results.extend(region_results)

            logger.info("[DONE] %s: %d condition results", region.label, len(region_results))

        except Exception as exc:
            logger.exception("Region failed: %s", region.label)
            if not args.keep_going:
                raise
            logger.error("Continuing after failure: %s", exc)

    if not all_results:
        logger.warning("No successful region runs.")
        return 1

    write_summary(all_results, args.outdir, include_htseq=include_htseq)
    write_diagnostics(all_results, args.outdir, include_htseq=include_htseq)
    logger.info("Wrote summaries to %s", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
