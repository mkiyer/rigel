#!/usr/bin/env python3
"""
Regional head-to-head benchmarking for hulkrna, salmon, and kallisto.

This script enables focused benchmarking on real genome regions:

1) Extract region genome FASTA + transcript annotations from a full genome+GTF
2) Build tool indexes (hulkrna, salmon, kallisto) for that region
3) Simulate reads with configurable strand specificity + gDNA contamination
4) Quantify with hulkrna/salmon/kallisto
5) Compare against simulation ground truth

Example
-------
PYTHONPATH=src conda run -n hulkrna python scripts/benchmark_region_competition.py \
  --genome /path/to/genome.fa \
  --gtf /path/to/genes.gtf.gz \
  --region chr1:155180000-155260000 \
  --region chr17:43043000-43126000 \
  --outdir regional_bench \
  --n-fragments 20000 \
  --strand-specificity 0.85 \
  --gdna-abundance 20 \
  --abundance-mode random
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import math
import shutil
import subprocess
import time
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from hulkrna.gtf import GTF
from hulkrna.index import HulkIndex, write_bed12
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, ReadSimulator, SimConfig, reverse_complement
from hulkrna.transcript import Transcript
from hulkrna.types import Interval, Strand


logger = logging.getLogger(__name__)


@dataclass(slots=True)
class RegionSpec:
    label: str
    chrom: str
    start0: int
    end0: int

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass(slots=True)
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


@dataclass(slots=True)
class RegionRunResult:
    region: str
    n_transcripts: int
    n_genes: int
    n_fragments: int
    n_gdna_expected: int
    hulkrna: ToolMetrics
    salmon: ToolMetrics
    kallisto: ToolMetrics


class StringGenome:
    """Minimal genome object compatible with ReadSimulator."""

    def __init__(self, seq: str, name: str):
        self._seq = seq
        self.name = name

    def __len__(self) -> int:
        return len(self._seq)

    def __getitem__(self, key):
        return self._seq[key]


def parse_region(text: str) -> RegionSpec:
    """Parse region string 'chr:start-end' (1-based inclusive coordinates)."""
    if ":" not in text or "-" not in text:
        raise ValueError(f"Invalid region format: {text}")
    chrom, rest = text.split(":", 1)
    start_s, end_s = rest.replace(",", "").split("-", 1)
    start1 = int(start_s)
    end1 = int(end_s)
    if start1 < 1 or end1 < start1:
        raise ValueError(f"Invalid region coordinates: {text}")
    label = f"{chrom}_{start1}_{end1}"
    return RegionSpec(label=label, chrom=chrom, start0=start1 - 1, end0=end1)


def parse_regions(args: argparse.Namespace) -> list[RegionSpec]:
    regions: list[RegionSpec] = []
    for reg in args.region:
        regions.append(parse_region(reg))
    if args.region_file is not None:
        with open(args.region_file) as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                start1 = int(fields[1])
                end1 = int(fields[2])
                label = fields[3] if len(fields) > 3 else f"{chrom}_{start1}_{end1}"
                regions.append(
                    RegionSpec(label=label, chrom=chrom, start0=start1 - 1, end0=end1)
                )
    if not regions:
        raise ValueError("No regions provided. Use --region and/or --region-file.")
    return regions


def open_textmaybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def extract_region(
    genome_fa: Path,
    gtf_path: Path,
    region: RegionSpec,
    out_dir: Path,
) -> tuple[Path, Path, list[Transcript], str]:
    """Extract region FASTA and transcript annotations fully contained in region."""
    out_dir.mkdir(parents=True, exist_ok=True)

    with pysam.FastaFile(str(genome_fa)) as fa:
        region_seq = fa.fetch(region.chrom, region.start0, region.end0)

    region_fa = out_dir / "region.fa"
    with open(region_fa, "w") as out:
        out.write(f">{region.label}\n")
        for i in range(0, len(region_seq), 80):
            out.write(region_seq[i : i + 80] + "\n")
    pysam.faidx(str(region_fa))

    # Parse exon rows and collect transcripts fully contained in region.
    # We only keep transcripts where ALL exons are in the window.
    by_tid: dict[str, dict] = {}
    with open_textmaybe_gzip(gtf_path) as fh:
        for rec in GTF.parse(fh):
            if rec.feature != "exon":
                continue
            if rec.seqname != region.chrom:
                continue

            tid = rec.attrs.get("transcript_id")
            gid = rec.attrs.get("gene_id")
            if tid is None or gid is None:
                continue

            exon_start = rec.start
            exon_end = rec.end

            data = by_tid.setdefault(
                tid,
                {
                    "g_id": gid,
                    "g_name": rec.attrs.get("gene_name", gid),
                    "g_type": rec.attrs.get("gene_type", "protein_coding"),
                    "strand": rec.strand,
                    "exons": [],
                    "all_inside": True,
                },
            )

            data["exons"].append((exon_start, exon_end))
            if exon_start < region.start0 or exon_end > region.end0:
                data["all_inside"] = False

    transcripts: list[Transcript] = []
    for tid, data in by_tid.items():
        if not data["all_inside"]:
            continue

        exons = sorted(data["exons"])
        if not exons:
            continue

        strand = Strand.from_str(data["strand"])
        rel_exons = [Interval(s - region.start0, e - region.start0) for s, e in exons]

        t = Transcript(
            ref=region.label,
            strand=strand,
            exons=rel_exons,
            t_id=tid,
            g_id=data["g_id"],
            g_name=data["g_name"],
            g_type=data["g_type"],
            is_basic=True,
        )
        t.compute_length()
        transcripts.append(t)

    transcripts.sort(key=lambda t: (t.g_id, t.t_id))

    # Assign indices
    g_id_to_idx: dict[str, int] = {}
    for t_idx, t in enumerate(transcripts):
        if t.g_id not in g_id_to_idx:
            g_id_to_idx[t.g_id] = len(g_id_to_idx)
        t.t_index = t_idx
        t.g_index = g_id_to_idx[t.g_id]

    region_gtf = out_dir / "region.gtf"
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

    return region_fa, region_gtf, transcripts, region_seq


def assign_abundances(
    transcripts: list[Transcript],
    mode: str,
    rng: np.random.Generator,
    min_abund: float,
    max_abund: float,
    abundance_file: Path | None = None,
) -> None:
    if mode == "uniform":
        for t in transcripts:
            t.abundance = 100.0
        return

    if mode == "random":
        lo = math.log10(min_abund)
        hi = math.log10(max_abund)
        vals = rng.uniform(lo, hi, size=len(transcripts))
        abund = np.power(10.0, vals)
        # Rescale around ~100 for readability
        abund = abund / np.median(abund) * 100.0
        for t, a in zip(transcripts, abund):
            t.abundance = float(a)
        return

    if mode == "file":
        if abundance_file is None:
            raise ValueError("--abundance-file is required when --abundance-mode=file")
        if not abundance_file.exists():
            raise FileNotFoundError(f"Abundance file not found: {abundance_file}")

        df = pd.read_csv(abundance_file, sep=None, engine="python")
        needed = {"transcript_id", "abundance"}
        if not needed.issubset(df.columns):
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


def write_transcript_fasta(transcripts: list[Transcript], region_seq: str, out_fa: Path) -> None:
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


def align_minimap2(
    region_fa: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_bam: Path,
    bed_path: Path | None = None,
) -> None:
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
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()

    _, p2_stderr = p2.communicate()
    p1_stderr = p1.stderr.read()
    p1.wait()

    if p2.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {p2_stderr.decode()}")
    if p1.returncode != 0:
        logger.warning("minimap2 returned rc=%d: %s", p1.returncode, p1_stderr.decode())


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


def run_salmon(transcript_fa: Path, fastq_r1: Path, fastq_r2: Path, out_dir: Path, threads: int) -> tuple[dict[str, float], float]:
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


def score_tool(tool: str, truth: dict[str, int], observed: dict[str, float], elapsed_sec: float) -> ToolMetrics:
    tids = sorted(truth)
    truth_arr = np.array([truth.get(t, 0) for t in tids], dtype=float)
    obs_arr = np.array([float(observed.get(t, 0.0)) for t in tids], dtype=float)

    diff = obs_arr - truth_arr
    abs_diff = np.abs(diff)
    total_abs_error = float(abs_diff.sum())
    mean_abs_error = float(abs_diff.mean()) if len(abs_diff) else 0.0
    rmse = float(np.sqrt(np.mean(diff * diff))) if len(diff) else 0.0

    if len(tids) > 1:
        # guard constant vectors
        pearson = float(np.corrcoef(truth_arr, obs_arr)[0, 1]) if np.std(truth_arr) > 0 and np.std(obs_arr) > 0 else 0.0
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


def run_region_benchmark(region: RegionSpec, args: argparse.Namespace, out_root: Path) -> RegionRunResult | None:
    reg_dir = out_root / region.label
    reg_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[REGION] %s:%d-%d", region.chrom, region.start0 + 1, region.end0)

    region_fa, region_gtf, transcripts, region_seq = extract_region(
        args.genome,
        args.gtf,
        region,
        reg_dir / "region",
    )
    if not transcripts:
        logger.warning("No fully-contained transcripts in region %s; skipping", region.label)
        return None

    rng = np.random.default_rng(args.abundance_seed + abs(hash(region.label)) % 10_000)
    assign_abundances(
        transcripts,
        mode=args.abundance_mode,
        rng=rng,
        min_abund=args.abundance_min,
        max_abund=args.abundance_max,
        abundance_file=args.abundance_file,
    )

    # Rewrite region GTF with abundance in score
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

    # Simulate reads
    sim_cfg = SimConfig(
        frag_mean=args.frag_mean,
        frag_std=args.frag_std,
        frag_min=args.frag_min,
        frag_max=args.frag_max,
        read_length=args.read_length,
        error_rate=args.error_rate,
        strand_specificity=args.strand_specificity,
        seed=args.sim_seed,
    )
    gdna_cfg = None
    if args.gdna_abundance > 0:
        gdna_cfg = GDNAConfig(
            abundance=args.gdna_abundance,
            frag_mean=args.gdna_frag_mean,
            frag_std=args.gdna_frag_std,
            frag_min=args.gdna_frag_min,
            frag_max=args.gdna_frag_max,
        )

    sim_genome = StringGenome(region_seq, region.label)
    simulator = ReadSimulator(sim_genome, transcripts, config=sim_cfg, gdna_config=gdna_cfg)
    fastq_dir = reg_dir / "reads"
    fastq_r1, fastq_r2 = simulator.write_fastq(fastq_dir, args.n_fragments, prefix="sim")
    truth_counts, n_gdna_expected = parse_truth_from_fastq(fastq_r1)

    # hulkrna path: align -> index -> count
    bam_path = reg_dir / "align" / "reads.bam"
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    # Generate BED12 for minimap2 -j annotation-guided alignment
    bed_path = reg_dir / "align" / "annotation.bed"
    write_bed12(transcripts, bed_path)
    align_minimap2(region_fa, fastq_r1, fastq_r2, bam_path, bed_path=bed_path)

    index_dir = reg_dir / "hulkrna_index"
    HulkIndex.build(region_fa, region_gtf, index_dir, write_tsv=False)
    index = HulkIndex.load(index_dir)

    t0 = time.monotonic()
    pipe = run_pipeline(
        bam_path,
        index,
        seed=args.pipeline_seed,
        sj_strand_tag="auto",
        include_multimap=True,
        gdna_threshold=args.gdna_threshold,
        gdna_splice_penalty_unannot=args.gdna_splice_penalty_unannot,
    )
    hulkrna_elapsed = time.monotonic() - t0

    hulkrna_counts_df = pipe.counter.get_counts_df(index)
    hulkrna_counts = {
        row.transcript_id: float(row.count)
        for row in hulkrna_counts_df.itertuples(index=False)
    }

    salmon_counts, salmon_elapsed = run_salmon(
        transcript_fa,
        fastq_r1,
        fastq_r2,
        reg_dir / "salmon",
        threads=args.threads,
    )
    kallisto_counts, kallisto_elapsed = run_kallisto(
        transcript_fa,
        fastq_r1,
        fastq_r2,
        reg_dir / "kallisto",
        threads=args.threads,
        strand_specificity=args.strand_specificity,
    )

    hulkrna_metrics = score_tool("hulkrna", truth_counts, hulkrna_counts, hulkrna_elapsed)
    salmon_metrics = score_tool("salmon", truth_counts, salmon_counts, salmon_elapsed)
    kallisto_metrics = score_tool("kallisto", truth_counts, kallisto_counts, kallisto_elapsed)

    per_tx_csv = reg_dir / "per_transcript_counts.csv"
    tids = sorted(truth_counts)
    with open(per_tx_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["transcript_id", "truth", "hulkrna", "salmon", "kallisto", "abundance"])
        abundance_map = {t.t_id: t.abundance for t in transcripts}
        for tid in tids:
            writer.writerow([
                tid,
                int(truth_counts.get(tid, 0)),
                float(hulkrna_counts.get(tid, 0.0)),
                float(salmon_counts.get(tid, 0.0)),
                float(kallisto_counts.get(tid, 0.0)),
                float(abundance_map.get(tid, 0.0)),
            ])

    return RegionRunResult(
        region=region.label,
        n_transcripts=len(transcripts),
        n_genes=len({t.g_id for t in transcripts}),
        n_fragments=args.n_fragments,
        n_gdna_expected=n_gdna_expected,
        hulkrna=hulkrna_metrics,
        salmon=salmon_metrics,
        kallisto=kallisto_metrics,
    )


def write_summary(results: list[RegionRunResult], out_root: Path) -> None:
    out_json = out_root / "summary.json"
    out_csv = out_root / "summary.csv"
    out_md = out_root / "summary.md"

    with open(out_json, "w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "region",
                "n_transcripts",
                "n_genes",
                "n_fragments",
                "n_gdna_expected",
                "tool",
                "elapsed_sec",
                "total_abs_error",
                "mean_abs_error",
                "rmse",
                "pearson",
                "spearman",
            ]
        )
        for r in results:
            for tm in (r.hulkrna, r.salmon, r.kallisto):
                writer.writerow(
                    [
                        r.region,
                        r.n_transcripts,
                        r.n_genes,
                        r.n_fragments,
                        r.n_gdna_expected,
                        tm.tool,
                        tm.elapsed_sec,
                        tm.total_abs_error,
                        tm.mean_abs_error,
                        tm.rmse,
                        tm.pearson,
                        tm.spearman,
                    ]
                )

    lines = [
        "# Regional Benchmark Summary",
        "",
        "| Region | Tool | Total Abs Error | Mean Abs Error | RMSE | Pearson | Spearman | Time (s) |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in results:
        for tm in (r.hulkrna, r.salmon, r.kallisto):
            lines.append(
                f"| {r.region} | {tm.tool} | {tm.total_abs_error:.2f} | {tm.mean_abs_error:.2f} | "
                f"{tm.rmse:.2f} | {tm.pearson:.4f} | {tm.spearman:.4f} | {tm.elapsed_sec:.2f} |"
            )
    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diagnostics(results: list[RegionRunResult], out_root: Path) -> None:
    """Write simple aggregate diagnostic report highlighting performance gaps."""
    rows = []
    for r in results:
        for tool_name in ("hulkrna", "salmon", "kallisto"):
            tm = getattr(r, tool_name)
            rows.append(
                {
                    "region": r.region,
                    "tool": tool_name,
                    "n_transcripts": r.n_transcripts,
                    "n_fragments": r.n_fragments,
                    "n_gdna_expected": r.n_gdna_expected,
                    "mae": tm.mean_abs_error,
                    "rmse": tm.rmse,
                    "pearson": tm.pearson,
                    "spearman": tm.spearman,
                }
            )
    df = pd.DataFrame(rows)
    agg = (
        df.groupby("tool")[["mae", "rmse", "pearson", "spearman"]]
        .mean()
        .sort_values("mae")
    )

    per_tx_parts: list[pd.DataFrame] = []
    for r in results:
        f = out_root / r.region / "per_transcript_counts.csv"
        if f.exists():
            part = pd.read_csv(f)
            part["region"] = r.region
            per_tx_parts.append(part)

    lines = [
        "# Diagnostic Notes",
        "",
        "## Aggregate tool metrics (mean across regions)",
        "",
        "| Tool | MAE | RMSE | Pearson | Spearman |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for tool, row in agg.iterrows():
        lines.append(
            f"| {tool} | {row['mae']:.3f} | {row['rmse']:.3f} | {row['pearson']:.4f} | {row['spearman']:.4f} |"
        )

    if per_tx_parts:
        alltx = pd.concat(per_tx_parts, ignore_index=True)
        for tool in ("hulkrna", "salmon", "kallisto"):
            alltx[f"ae_{tool}"] = (alltx[tool] - alltx["truth"]).abs()
        alltx["ab_q"] = pd.qcut(
            alltx["abundance"].rank(method="first"),
            4,
            labels=["Q1_low", "Q2", "Q3", "Q4_high"],
        )

        lines.extend(
            [
                "",
                "## Mean absolute error by abundance quartile",
                "",
                "| Quartile | hulkrna | salmon | kallisto |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
            sub = alltx[alltx["ab_q"] == q]
            if len(sub) == 0:
                continue
            lines.append(
                f"| {q} | {sub['ae_hulkrna'].mean():.3f} | {sub['ae_salmon'].mean():.3f} | {sub['ae_kallisto'].mean():.3f} |"
            )

        lines.extend(
            [
                "",
                "## Dropout rate (truth > 0 but predicted <= 0)",
                "",
                "| Tool | Dropout rate |",
                "| --- | ---: |",
            ]
        )
        for tool in ("hulkrna", "salmon", "kallisto"):
            rate = float(((alltx["truth"] > 0) & (alltx[tool] <= 0)).mean())
            lines.append(f"| {tool} | {rate:.4f} |")

    out = out_root / "diagnostics.md"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Regional hulkrna/salmon/kallisto benchmarking")
    p.add_argument("--genome", type=Path, required=True, help="Genome FASTA path")
    p.add_argument("--gtf", type=Path, required=True, help="Gene annotation GTF/GTF.GZ path")
    p.add_argument("--region", action="append", default=[], help="Region: chr:start-end (1-based)")
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

    p.add_argument("--frag-mean", type=float, default=250.0)
    p.add_argument("--frag-std", type=float, default=50.0)
    p.add_argument("--frag-min", type=int, default=50)
    p.add_argument("--frag-max", type=int, default=1000)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--error-rate", type=float, default=0.0)

    p.add_argument("--strand-specificity", type=float, default=1.0)
    p.add_argument("--gdna-abundance", type=float, default=0.0)
    p.add_argument("--gdna-frag-mean", type=float, default=350.0)
    p.add_argument("--gdna-frag-std", type=float, default=100.0)
    p.add_argument("--gdna-frag-min", type=int, default=100)
    p.add_argument("--gdna-frag-max", type=int, default=1000)

    p.add_argument("--gdna-threshold", type=float, default=0.5)
    p.add_argument("--gdna-splice-penalty-unannot", type=float, default=0.01)

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
        help="CSV/TSV with columns: transcript_id, abundance (used when --abundance-mode file)",
    )

    p.add_argument("--keep-going", action="store_true", help="Continue with next region on failure")
    p.add_argument("--verbose", action="store_true")
    return p


def ensure_tools() -> None:
    required = ["minimap2", "samtools", "salmon", "kallisto"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(f"Missing required tools in PATH: {', '.join(missing)}")


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    ensure_tools()
    if args.abundance_mode == "file" and args.abundance_file is None:
        raise ValueError("--abundance-file is required when --abundance-mode=file")
    regions = parse_regions(args)

    args.outdir.mkdir(parents=True, exist_ok=True)

    results: list[RegionRunResult] = []
    for region in regions:
        try:
            res = run_region_benchmark(region, args, args.outdir)
            if res is not None:
                results.append(res)
                logger.info(
                    "[DONE] %s | hulkrna MAE=%.2f | salmon MAE=%.2f | kallisto MAE=%.2f",
                    res.region,
                    res.hulkrna.mean_abs_error,
                    res.salmon.mean_abs_error,
                    res.kallisto.mean_abs_error,
                )
        except Exception as exc:
            logger.exception("Region failed: %s", region.label)
            if not args.keep_going:
                raise
            logger.error("Continuing after failure: %s", exc)

    if not results:
        logger.warning("No successful region runs.")
        return 1

    write_summary(results, args.outdir)
    write_diagnostics(results, args.outdir)
    logger.info("Wrote summaries to %s", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
