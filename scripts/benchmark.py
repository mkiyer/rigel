#!/usr/bin/env python3
"""Benchmark rigel (and optionally salmon, kallisto, htseq).

.. note:: Supports gzipped FASTQ files (.fq.gz) transparently.

Works with two data sources:

1. **Simulated data** — Point ``sim_dir`` to the output of ``sim.py``
   (reads ``manifest.json`` to auto-discover conditions, FASTQ, oracle
   BAM, and truth).

2. **Explicit datasets** — Provide a list of datasets in the YAML
   config, each with FASTQ or BAM paths and an optional truth file.

For simulated data the benchmark scores tool estimates against
per-transcript ground truth parsed from FASTQ read names.  For real
data with no truth file, tool outputs are reported without scoring.

Usage
-----
::

    python scripts/benchmark.py --config scripts/benchmark_example.yaml

Output
------
::

    <outdir>/
        rigel_index/
        annotation.bed12
        <condition>/
            align_<aligner>/
            per_transcript_counts_<aligner>.csv
        summary.json
        summary.csv
        summary.md
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]

import pysam

from rigel.config import (
    BamScanConfig,
    EMConfig,
    FragmentScoringConfig,
    PipelineConfig,
)
from rigel.index import TranscriptIndex, write_bed12
from rigel.pipeline import run_pipeline
from rigel.scoring import (
    GDNA_SPLICE_PENALTIES,
    SPLICE_UNANNOT,
    overhang_alpha_to_log_penalty,
)
from rigel.transcript import Transcript
from rigel.types import Strand

logger = logging.getLogger(__name__)

# ── Utilities ───────────────────────────────────────────────────────

_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _find_tool(name: str) -> str:
    """Resolve *name* to an absolute path, checking the active env first."""
    env_bin = Path(sys.executable).resolve().parent / name
    if env_bin.is_file():
        return str(env_bin)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"Tool '{name}' not found in PATH or {env_bin.parent}"
    )


def reverse_complement(seq: str) -> str:
    return seq.translate(_DNA_COMPLEMENT)[::-1]


def _get_peak_rss_mb() -> float:
    """Peak resident set size of this process in MB."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    if sys.platform == "darwin":
        return usage.ru_maxrss / (1024 * 1024)  # bytes on macOS
    return usage.ru_maxrss / 1024  # KB on Linux


# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════


@dataclass
class HulkrnaConfig:
    """Named rigel parameterization for benchmarking."""

    name: str
    params: dict = field(default_factory=dict)


@dataclass
class AlignerConfig:
    """Named aligner configuration."""

    name: str
    type: str  # "oracle", "minimap2", "hisat2"
    params: dict = field(default_factory=dict)


@dataclass
class Dataset:
    """A single dataset to benchmark — one simulation condition or
    one real-data sample."""

    name: str
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None
    oracle_bam: Path | None = None
    truth_path: Path | None = None
    # Metadata (from sim manifest or user-supplied)
    gdna_label: str = ""
    gdna_rate: float = 0.0
    strand_specificity: float = 1.0
    n_rna: int = 0
    n_gdna: int = 0

    # Pipeline seed (from simulation, for EM reproducibility)
    pipeline_seed: int = 42


@dataclass
class BenchmarkConfig:
    """Top-level benchmark configuration."""

    genome: str = ""
    gtf: str = ""
    outdir: str = "benchmark_output"
    threads: int = 8
    transcript_filter: str = "all"

    sim_dir: str | None = None
    datasets: list[Dataset] = field(default_factory=list)

    aligners: list[AlignerConfig] = field(default_factory=list)
    rigel_configs: list[HulkrnaConfig] = field(default_factory=list)

    salmon_enabled: bool = False
    salmon_index: str | None = None
    kallisto_enabled: bool = False
    kallisto_index: str | None = None
    htseq_enabled: bool = False

    # Minimap2 pre-built index (.mmi) and splice junction BED
    minimap2_index: str | None = None
    minimap2_bed: str | None = None

    keep_going: bool = True
    verbose: bool = True


def parse_yaml_config(path: str | Path) -> BenchmarkConfig:
    """Parse a YAML config file into a BenchmarkConfig."""
    if yaml is None:
        raise ImportError("PyYAML is required: pip install pyyaml")

    with open(path) as f:
        raw = yaml.safe_load(f)

    cfg = BenchmarkConfig()
    cfg.genome = raw.get("genome", "")
    cfg.gtf = raw.get("gtf", "")
    cfg.outdir = raw.get("outdir", "benchmark_output")
    cfg.threads = int(raw.get("threads", 8))
    cfg.transcript_filter = raw.get("transcript_filter", "all")

    cfg.sim_dir = raw.get("sim_dir", None)

    # Explicit datasets
    for ds_raw in raw.get("datasets", []):
        ds = Dataset(name=ds_raw["name"])
        if "fastq_r1" in ds_raw:
            ds.fastq_r1 = Path(ds_raw["fastq_r1"])
        if "fastq_r2" in ds_raw:
            ds.fastq_r2 = Path(ds_raw["fastq_r2"])
        if "bam" in ds_raw:
            ds.oracle_bam = Path(ds_raw["bam"])
        if "oracle_bam" in ds_raw:
            ds.oracle_bam = Path(ds_raw["oracle_bam"])
        if "truth" in ds_raw:
            ds.truth_path = Path(ds_raw["truth"])
        ds.gdna_label = ds_raw.get("gdna_label", "")
        ds.gdna_rate = float(ds_raw.get("gdna_rate", 0.0))
        ds.strand_specificity = float(ds_raw.get("strand_specificity", 1.0))
        ds.n_rna = int(ds_raw.get("n_rna", 0))
        ds.n_gdna = int(ds_raw.get("n_gdna", 0))
        cfg.datasets.append(ds)

    # Aligners
    aligners_raw = raw.get("aligners", {"oracle": {}})
    for name, params in aligners_raw.items():
        params = params or {}
        atype = params.pop("type", name)
        if atype not in ("oracle", "minimap2", "hisat2"):
            atype = name
        cfg.aligners.append(AlignerConfig(name=name, type=atype, params=params))

    # HulkRNA configs
    rigel_raw = raw.get("rigel_configs", {"default": {}})
    for name, params in rigel_raw.items():
        cfg.rigel_configs.append(HulkrnaConfig(name=name, params=params or {}))

    # Optional tools
    salmon_raw = raw.get("salmon", {})
    cfg.salmon_enabled = bool(salmon_raw.get("enabled", False)) if salmon_raw else False
    if salmon_raw:
        cfg.salmon_index = salmon_raw.get("index", None)
    kallisto_raw = raw.get("kallisto", {})
    cfg.kallisto_enabled = bool(kallisto_raw.get("enabled", False)) if kallisto_raw else False
    if kallisto_raw:
        cfg.kallisto_index = kallisto_raw.get("index", None)
    htseq_raw = raw.get("htseq", {})
    cfg.htseq_enabled = bool(htseq_raw.get("enabled", False)) if htseq_raw else False

    # Minimap2 pre-built index and BED
    mm2_raw = raw.get("minimap2", {})
    if mm2_raw:
        cfg.minimap2_index = mm2_raw.get("index", None)
        cfg.minimap2_bed = mm2_raw.get("bed", None)

    # Runtime
    rt_raw = raw.get("runtime", {})
    cfg.keep_going = bool(rt_raw.get("keep_going", True))
    cfg.verbose = bool(raw.get("verbose", True))

    return cfg


# ═══════════════════════════════════════════════════════════════════
# Sim manifest loading
# ═══════════════════════════════════════════════════════════════════


def load_sim_manifest(sim_dir: str | Path) -> tuple[dict, list[Dataset]]:
    """Load manifest.json from a sim.py output directory.

    Returns
    -------
    manifest : dict
        Full manifest contents (genome, gtf, etc.).
    datasets : list[Dataset]
        One Dataset per simulation condition.
    """
    sim_dir = Path(sim_dir)
    manifest_path = sim_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    with open(manifest_path) as f:
        manifest = json.load(f)

    truth_path = sim_dir / manifest.get("truth_abundances", "truth_abundances.tsv")
    sim_params = manifest.get("simulation", {})
    pipeline_seed = sim_params.get("pipeline_seed", sim_params.get("sim_seed", 42))

    datasets: list[Dataset] = []
    for cond in manifest.get("conditions", []):
        ds = Dataset(name=cond["name"])
        ds.gdna_label = cond.get("gdna_label", "")
        ds.gdna_rate = float(cond.get("gdna_rate", 0.0))
        ds.strand_specificity = float(cond.get("strand_specificity", 1.0))
        ds.n_rna = int(cond.get("n_rna", 0))
        ds.n_gdna = int(cond.get("n_gdna", 0))
        ds.pipeline_seed = pipeline_seed

        if cond.get("fastq_r1"):
            ds.fastq_r1 = sim_dir / cond["fastq_r1"]
        if cond.get("fastq_r2"):
            ds.fastq_r2 = sim_dir / cond["fastq_r2"]
        if cond.get("oracle_bam"):
            ds.oracle_bam = sim_dir / cond["oracle_bam"]
        if truth_path.exists():
            ds.truth_path = truth_path

        datasets.append(ds)

    return manifest, datasets


# ═══════════════════════════════════════════════════════════════════
# Transcript loading
# ═══════════════════════════════════════════════════════════════════


def load_transcripts(
    gtf_path: str | Path,
    *,
    transcript_filter: str = "all",
) -> list[Transcript]:
    """Load transcripts from a GTF with optional filtering."""
    logger.info("Loading transcripts from %s (filter=%s)", gtf_path, transcript_filter)
    transcripts = Transcript.read_gtf(str(gtf_path), parse_mode="warn-skip")
    logger.info("Read %d transcripts from GTF", len(transcripts))

    if transcript_filter == "basic":
        transcripts = [t for t in transcripts if t.is_basic]
    elif transcript_filter == "mane":
        transcripts = [t for t in transcripts if t.is_mane]
    elif transcript_filter == "ccds":
        transcripts = [t for t in transcripts if t.is_ccds]

    for i, t in enumerate(transcripts):
        t.t_index = i

    logger.info("Final: %d transcripts", len(transcripts))
    return transcripts


# ═══════════════════════════════════════════════════════════════════
# Index building
# ═══════════════════════════════════════════════════════════════════


def build_rigel_index(
    genome_path: str | Path,
    gtf_path: str | Path,
    index_dir: Path,
) -> TranscriptIndex:
    """Build rigel index from full genome + GTF."""
    required = ["transcripts.feather", "intervals.feather", "ref_lengths.feather"]
    if all((index_dir / f).exists() for f in required):
        logger.info("rigel index already exists at %s, loading...", index_dir)
    else:
        logger.info("Building rigel index → %s", index_dir)
        TranscriptIndex.build(
            fasta_file=genome_path,
            gtf_file=gtf_path,
            output_dir=index_dir,
            write_tsv=True,
            gtf_parse_mode="warn-skip",
        )
    idx = TranscriptIndex.load(str(index_dir))
    logger.info("Loaded rigel index: %d transcripts, %d genes",
                idx.num_transcripts, idx.num_genes)
    return idx


def write_transcript_fasta(
    fasta_path: str | Path,
    transcripts: list[Transcript],
    out_path: Path,
) -> Path:
    """Write spliced transcript FASTA for salmon/kallisto."""
    fasta_fh = pysam.FastaFile(str(fasta_path))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with open(out_path, "w") as out:
        for t in transcripts:
            exon_seqs = []
            for e in t.exons:
                seq = fasta_fh.fetch(t.ref, e.start, e.end)
                exon_seqs.append(seq.upper())
            seq = "".join(exon_seqs)
            if t.strand == Strand.NEG:
                seq = reverse_complement(seq)
            if not seq:
                continue
            out.write(f">{t.t_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
            n_written += 1
    fasta_fh.close()
    logger.info("Wrote %d transcript sequences to %s", n_written, out_path)
    return out_path


def build_salmon_index(
    transcript_fa: Path,
    out_dir: Path,
    threads: int = 8,
    prebuilt: str | None = None,
) -> Path:
    """Return a salmon index directory, using *prebuilt* if provided."""
    if prebuilt:
        p = Path(prebuilt)
        if p.is_dir():
            logger.info("Using pre-built salmon index: %s", p)
            return p
        raise FileNotFoundError(f"Pre-built salmon index not found: {p}")
    idx_dir = out_dir / "salmon_index"
    if idx_dir.exists():
        logger.info("Salmon index already exists at %s", idx_dir)
        return idx_dir
    logger.info("Building salmon index → %s", idx_dir)
    subprocess.run(
        [_find_tool("salmon"), "index", "-t", str(transcript_fa), "-i", str(idx_dir),
         "-k", "23", "-p", str(threads)],
        check=True, capture_output=True, text=True,
    )
    return idx_dir


def build_kallisto_index(
    transcript_fa: Path,
    out_dir: Path,
    prebuilt: str | None = None,
) -> Path:
    """Return a kallisto index file, using *prebuilt* if provided."""
    if prebuilt:
        p = Path(prebuilt)
        if p.is_file() or p.is_dir():
            logger.info("Using pre-built kallisto index: %s", p)
            return p
        raise FileNotFoundError(f"Pre-built kallisto index not found: {p}")
    idx_path = out_dir / "kallisto_index.idx"
    if idx_path.exists():
        logger.info("Kallisto index already exists at %s", idx_path)
        return idx_path
    logger.info("Building kallisto index → %s", idx_path)
    subprocess.run(
        [_find_tool("kallisto"), "index", "-i", str(idx_path), str(transcript_fa)],
        check=True, capture_output=True, text=True,
    )
    return idx_path


def build_minimap2_index(
    genome_fa: Path,
    out_dir: Path,
    threads: int = 8,
    prebuilt: str | None = None,
) -> Path:
    """Return a minimap2 .mmi index, using *prebuilt* if provided."""
    if prebuilt:
        p = Path(prebuilt)
        if p.is_file():
            logger.info("Using pre-built minimap2 index: %s", p)
            return p
        raise FileNotFoundError(f"Pre-built minimap2 index not found: {p}")
    idx_path = out_dir / "genome.mmi"
    if idx_path.exists():
        logger.info("minimap2 index already exists at %s", idx_path)
        return idx_path
    logger.info("Building minimap2 index → %s (this may take a while)", idx_path)
    subprocess.run(
        [_find_tool("minimap2"), "-x", "splice:sr", "-d", str(idx_path),
         "-t", str(threads), str(genome_fa)],
        check=True, capture_output=True, text=True,
    )
    return idx_path


def build_minimap2_bed(
    gtf_path: Path,
    out_dir: Path,
    prebuilt: str | None = None,
) -> Path:
    """Return a splice-junction BED for minimap2 ``--junc-bed``.

    Uses minimap2's bundled ``paftools.js gff2bed`` which produces the
    12-column BED that minimap2's ``--junc-bed`` expects.  Falls back
    to rigel's ``write_bed12`` if paftools.js is not available.
    """
    if prebuilt:
        p = Path(prebuilt)
        if p.is_file():
            logger.info("Using pre-built minimap2 BED: %s", p)
            return p
        raise FileNotFoundError(f"Pre-built minimap2 BED not found: {p}")
    bed_path = out_dir / "minimap2_annotation.bed"
    if bed_path.exists():
        logger.info("minimap2 BED already exists at %s", bed_path)
        return bed_path
    logger.info("Building minimap2 splice BED via paftools.js → %s", bed_path)
    with open(bed_path, "w") as fh:
        subprocess.run(
            [_find_tool("k8"), _find_tool("paftools.js"), "gff2bed", str(gtf_path)],
            check=True, stdout=fh, stderr=subprocess.PIPE,
        )
    return bed_path


# ═══════════════════════════════════════════════════════════════════
# Alignment
# ═══════════════════════════════════════════════════════════════════


def align_minimap2(
    genome_or_index: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_bam: Path,
    *,
    bed_path: Path | None = None,
    threads: int = 8,
) -> None:
    """Align paired-end reads with minimap2 → name-sorted BAM.

    *genome_or_index* may be a raw FASTA or a pre-built ``.mmi`` index.
    """
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    mm2_cmd = [
        _find_tool("minimap2"), "-a", "-x", "splice:sr",
        "--secondary=yes", "-N", "10", "-t", str(threads),
    ]
    if bed_path is not None:
        mm2_cmd.extend(["--junc-bed", str(bed_path)])
    mm2_cmd.extend([str(genome_or_index), str(fastq_r1), str(fastq_r2)])

    sort_cmd = [
        _find_tool("samtools"), "sort", "-n",
        "-@", str(max(1, threads // 2)),
        "-o", str(out_bam),
    ]

    logger.info("Running minimap2 | samtools sort -n → %s", out_bam)
    p1 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    _, stderr2 = p2.communicate()
    p1.wait()
    if p1.returncode != 0 or p2.returncode != 0:
        raise RuntimeError(f"minimap2 alignment failed: mm2={p1.returncode}, sort={p2.returncode}")


def coord_sort_bam(
    input_bam: Path, output_bam: Path, threads: int = 4,
) -> Path:
    """Re-sort a BAM to coordinate order and index."""
    subprocess.run(
        [_find_tool("samtools"), "sort", "-@", str(threads), "-o", str(output_bam), str(input_bam)],
        check=True, capture_output=True,
    )
    subprocess.run([_find_tool("samtools"), "index", str(output_bam)], check=True, capture_output=True)
    return output_bam


# ═══════════════════════════════════════════════════════════════════
# Tool runners
# ═══════════════════════════════════════════════════════════════════


def _build_pipeline_config(
    pipeline_seed: int = 42,
    rigel_config: HulkrnaConfig | None = None,
    annotated_bam_path: Path | None = None,
) -> PipelineConfig:
    """Build PipelineConfig from benchmark config + overrides."""
    raw: dict = {}
    if rigel_config is not None:
        raw.update(rigel_config.params)

    em_kw: dict = {"seed": pipeline_seed}
    _EM_ALIASES = {
        "em_convergence_delta": "convergence_delta",
        "em_iterations": "iterations",
        "em_prior_alpha": "prior_alpha",
        "em_prior_gamma": "prior_gamma",
        "em_mode": "mode",
        "prune_threshold": "prune_threshold",
        "confidence_threshold": "confidence_threshold",
    }
    for raw_key, cfg_key in _EM_ALIASES.items():
        if raw_key in raw:
            em_kw[cfg_key] = raw.pop(raw_key)

    scoring_kw: dict = {}
    ov_alpha = raw.pop("overhang_alpha", None)
    if ov_alpha is not None:
        scoring_kw["overhang_log_penalty"] = overhang_alpha_to_log_penalty(ov_alpha)
    mm_alpha = raw.pop("mismatch_alpha", None)
    if mm_alpha is not None:
        scoring_kw["mismatch_log_penalty"] = overhang_alpha_to_log_penalty(mm_alpha)
    gdna_pen = raw.pop("gdna_splice_penalty_unannot", None)
    if gdna_pen is not None:
        penalties = dict(GDNA_SPLICE_PENALTIES)
        penalties[SPLICE_UNANNOT] = gdna_pen
        scoring_kw["gdna_splice_penalties"] = penalties

    scan_kw: dict = {"sj_strand_tag": "auto", "include_multimap": True}

    return PipelineConfig(
        em=EMConfig(**em_kw),
        scan=BamScanConfig(**scan_kw),
        scoring=FragmentScoringConfig(**scoring_kw),
        annotated_bam_path=str(annotated_bam_path) if annotated_bam_path else None,
    )


def run_rigel_tool(
    bam_path: Path,
    index: TranscriptIndex,
    pipeline_seed: int = 42,
    rigel_config: HulkrnaConfig | None = None,
    annotated_bam_path: Path | None = None,
) -> tuple[dict[str, float], float, dict[str, float], int]:
    """Run rigel pipeline → (transcript_counts, elapsed, pool_counts, n_fragments)."""
    t0 = time.monotonic()
    cfg = _build_pipeline_config(pipeline_seed, rigel_config, annotated_bam_path)
    pipe = run_pipeline(bam_path, index, config=cfg)
    elapsed = time.monotonic() - t0

    n_fragments = pipe.stats.total

    counts_df = pipe.estimator.get_counts_df(index)
    counts = {
        row.transcript_id: float(row.mrna)
        for row in counts_df.itertuples(index=False)
    }
    mature_pred = float(sum(counts.values()))
    nascent_pred = float(pipe.estimator.nrna_em_counts.sum())
    genomic_pred = float(pipe.stats.n_gdna_total)
    pool_counts = {
        "mature_rna": mature_pred,
        "nascent_rna": nascent_pred,
        "genomic_dna": genomic_pred,
    }
    return counts, elapsed, pool_counts, n_fragments


def run_salmon(
    salmon_index: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_dir: Path,
    threads: int,
) -> tuple[dict[str, float], float]:
    quant_dir = out_dir / "quant"
    quant_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    subprocess.run(
        [_find_tool("salmon"), "quant", "-i", str(salmon_index), "-l", "A",
         "-1", str(fastq_r1), "-2", str(fastq_r2),
         "-o", str(quant_dir), "-p", str(threads), "--validateMappings"],
        check=True, capture_output=True, text=True,
    )
    elapsed = time.monotonic() - t0
    df = pd.read_csv(quant_dir / "quant.sf", sep="\t")
    return dict(zip(df["Name"], df["NumReads"])), elapsed


def run_kallisto(
    kallisto_index: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_dir: Path,
    threads: int,
    strand_specificity: float,
) -> tuple[dict[str, float], float]:
    quant_dir = out_dir / "quant"
    quant_dir.mkdir(parents=True, exist_ok=True)
    cmd = [_find_tool("kallisto"), "quant", "-i", str(kallisto_index),
           "-o", str(quant_dir), "-t", str(threads)]
    if strand_specificity >= 0.9:
        cmd.append("--rf-stranded")
    cmd.extend([str(fastq_r1), str(fastq_r2)])
    t0 = time.monotonic()
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.monotonic() - t0
    df = pd.read_csv(quant_dir / "abundance.tsv", sep="\t")
    # Kallisto preserves full FASTA headers (pipe-delimited); extract transcript ID
    ids = df["target_id"].str.split("|").str[0]
    return dict(zip(ids, df["est_counts"])), elapsed


def run_htseq(
    bam_path: Path,
    gtf_path: str | Path,
    out_dir: Path,
    strand_specificity: float,
    conda_env: str = "htseq",
) -> tuple[dict[str, int], float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    stranded = "reverse" if strand_specificity >= 0.9 else "no"
    cmd = ["conda", "run", "-n", conda_env,
           "htseq-count", "-f", "bam", "-r", "pos",
           f"--stranded={stranded}", "-t", "exon", "-i", "gene_id",
           str(bam_path), str(gtf_path)]
    t0 = time.monotonic()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    elapsed = time.monotonic() - t0
    counts: dict[str, int] = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2 or parts[0].startswith("__"):
            continue
        counts[parts[0]] = int(parts[1])
    return counts, elapsed


# ═══════════════════════════════════════════════════════════════════
# Truth parsing
# ═══════════════════════════════════════════════════════════════════


def parse_truth_from_fastq(
    r1_path: Path,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Parse ground-truth from FASTQ read names.

    Returns
    -------
    mrna_counts : dict[transcript_id, count]
    nrna_counts : dict[transcript_id, count]
    n_gdna : int
    """
    mrna_counts: Counter[str] = Counter()
    nrna_counts: Counter[str] = Counter()
    n_gdna = 0

    opener = gzip.open if str(r1_path).endswith('.gz') else open
    with opener(r1_path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            qname = line[1:].strip()
            if qname.endswith("/1"):
                qname = qname[:-2]
            t_id = qname.split(":")[0]
            if t_id.startswith("gdna"):
                n_gdna += 1
            elif t_id.startswith("nrna_"):
                nrna_counts[t_id[5:]] += 1
            else:
                mrna_counts[t_id] += 1

    return dict(mrna_counts), dict(nrna_counts), n_gdna


# ═══════════════════════════════════════════════════════════════════
# Scoring
# ═══════════════════════════════════════════════════════════════════


@dataclass
class ToolMetrics:
    """Aggregate accuracy metrics for one tool in one condition."""

    tool: str
    elapsed_sec: float
    peak_rss_mb: float
    throughput_frags_per_sec: float
    n_transcripts: int
    total_truth: float
    total_observed: float
    total_abs_error: float
    mean_abs_error: float
    rmse: float
    pearson: float
    spearman: float


def _spearman_r(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation using numpy (no scipy)."""
    n = len(x)
    if n < 2:
        return 0.0
    rx = np.empty(n, dtype=np.float64)
    ry = np.empty(n, dtype=np.float64)
    for src, dst in ((x, rx), (y, ry)):
        order = np.argsort(src)
        # Average ranks for ties
        i = 0
        while i < n:
            j = i + 1
            while j < n and src[order[j]] == src[order[i]]:
                j += 1
            avg_rank = (i + j - 1) / 2.0  # 0-based average
            for k in range(i, j):
                dst[order[k]] = avg_rank
            i = j
    # Pearson on ranks
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def score_tool(
    tool: str,
    truth: dict[str, float],
    observed: dict[str, float],
    elapsed: float,
    peak_rss_mb: float = 0.0,
    throughput: float = 0.0,
) -> ToolMetrics:
    """Compare observed counts against truth → ToolMetrics."""
    all_tids = sorted(set(truth) | set(observed))
    n = len(all_tids)
    if n == 0:
        return ToolMetrics(tool=tool, elapsed_sec=elapsed,
                           peak_rss_mb=peak_rss_mb,
                           throughput_frags_per_sec=throughput,
                           n_transcripts=0,
                           total_truth=0, total_observed=0, total_abs_error=0,
                           mean_abs_error=0, rmse=0, pearson=0, spearman=0)

    truth_arr = np.array([truth.get(t, 0.0) for t in all_tids])
    obs_arr = np.array([observed.get(t, 0.0) for t in all_tids])

    abs_err = np.abs(truth_arr - obs_arr)
    mae = float(abs_err.mean())
    rmse = float(np.sqrt(np.mean((truth_arr - obs_arr) ** 2)))

    if truth_arr.std() > 0 and obs_arr.std() > 0:
        pearson = float(np.corrcoef(truth_arr, obs_arr)[0, 1])
        spearman = _spearman_r(truth_arr, obs_arr)
    else:
        pearson = 0.0
        spearman = 0.0

    return ToolMetrics(
        tool=tool, elapsed_sec=elapsed,
        peak_rss_mb=peak_rss_mb,
        throughput_frags_per_sec=throughput,
        n_transcripts=n,
        total_truth=float(truth_arr.sum()),
        total_observed=float(obs_arr.sum()),
        total_abs_error=float(abs_err.sum()),
        mean_abs_error=mae, rmse=rmse,
        pearson=pearson, spearman=spearman,
    )


def aggregate_to_genes(
    tx_counts: dict[str, float],
    transcripts: list[Transcript],
) -> dict[str, float]:
    """Sum transcript counts to gene-level counts."""
    gene_counts: dict[str, float] = defaultdict(float)
    tid_to_gid = {t.t_id: t.g_id for t in transcripts}
    for tid, cnt in tx_counts.items():
        gid = tid_to_gid.get(tid, tid)
        gene_counts[gid] += cnt
    return dict(gene_counts)


# ═══════════════════════════════════════════════════════════════════
# Result dataclass
# ═══════════════════════════════════════════════════════════════════


@dataclass
class ConditionResult:
    """Results for one (dataset × aligner) combination."""

    dataset_name: str
    gdna_label: str
    gdna_rate: float
    n_gdna: int
    strand_specificity: float
    n_rna_fragments: int
    n_total_fragments: int
    n_gdna_truth: int
    n_mrna_truth: int
    n_nrna_truth: int
    aligner: str
    transcript_metrics: dict  # tool → ToolMetrics (as dict)
    gene_metrics: dict
    pool_counts: dict  # tool → {mature_rna, nascent_rna, genomic_dna}
    elapsed: dict  # tool → elapsed_sec
    peak_rss: dict  # tool → peak_rss_mb
    throughput: dict  # tool → frags/sec
    alignment_elapsed: float  # seconds for this aligner


# ═══════════════════════════════════════════════════════════════════
# Output writers
# ═══════════════════════════════════════════════════════════════════


def write_summary_json(results: list[ConditionResult], path: Path) -> None:
    data = [asdict(r) for r in results]
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    logger.info("Wrote %s", path)


def write_summary_csv(results: list[ConditionResult], path: Path) -> None:
    """Tidy long-form CSV: one row per (condition × tool × level)."""
    rows: list[dict] = []
    for r in results:
        base = {
            "dataset": r.dataset_name,
            "gdna_label": r.gdna_label,
            "gdna_rate": r.gdna_rate,
            "strand_specificity": r.strand_specificity,
            "aligner": r.aligner,
            "n_rna_fragments": r.n_rna_fragments,
            "n_gdna": r.n_gdna,
            "alignment_elapsed_sec": r.alignment_elapsed,
        }
        for tool, m in r.transcript_metrics.items():
            row = {**base, "level": "transcript", "tool": tool}
            if isinstance(m, dict):
                row.update({k: v for k, v in m.items() if k != "tool"})
            rows.append(row)
        for tool, m in r.gene_metrics.items():
            row = {**base, "level": "gene", "tool": tool}
            if isinstance(m, dict):
                row.update({k: v for k, v in m.items() if k != "tool"})
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    logger.info("Wrote %s (%d rows)", path, len(df))


def write_summary_md(results: list[ConditionResult], path: Path) -> None:
    with open(path, "w") as f:
        f.write("# Benchmark Results\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Transcript-level table
        f.write("## Transcript-Level Metrics\n\n")
        f.write("| Dataset | Aligner | Tool | MAE | RMSE | Pearson | Spearman "
                "| Time (s) | RSS (MB) | Throughput (frag/s) |\n")
        f.write("|---------|---------|------|-----|------|---------|----------"
                "|----------|----------|---------------------|\n")
        for r in results:
            for tool, m in r.transcript_metrics.items():
                if isinstance(m, dict):
                    rss = m.get('peak_rss_mb', 0)
                    rss_str = f"{rss:.0f}" if rss > 0 else "\u2014"
                    tp = m.get('throughput_frags_per_sec', 0)
                    tp_str = f"{tp:,.0f}" if tp > 0 else "\u2014"
                    f.write(
                        f"| {r.dataset_name} | {r.aligner} | {tool} | "
                        f"{m.get('mean_abs_error', 0):.2f} | "
                        f"{m.get('rmse', 0):.2f} | "
                        f"{m.get('pearson', 0):.4f} | "
                        f"{m.get('spearman', 0):.4f} | "
                        f"{m.get('elapsed_sec', 0):.1f} | "
                        f"{rss_str} | {tp_str} |\n"
                    )

        # Pool-level summary
        f.write("\n## Pool-Level Counts\n\n")
        f.write("| Dataset | Aligner | Tool | mRNA (truth) | mRNA (pred) | "
                "nRNA (truth) | nRNA (pred) | gDNA (truth) | gDNA (pred) |\n")
        f.write("|---------|---------|------|-------------|-------------|"
                "-------------|-------------|-------------|-------------|\n")
        for r in results:
            for tool, pc in r.pool_counts.items():
                if isinstance(pc, dict):
                    f.write(
                        f"| {r.dataset_name} | {r.aligner} | {tool} | "
                        f"{r.n_mrna_truth} | {pc.get('mature_rna', 0):.0f} | "
                        f"{r.n_nrna_truth} | {pc.get('nascent_rna', 0):.0f} | "
                        f"{r.n_gdna_truth} | {pc.get('genomic_dna', 0):.0f} |\n"
                    )

        # Performance summary (always shown, even without truth)
        f.write("\n## Performance\n\n")
        f.write("| Dataset | Aligner | Align (s) | Tool | Quant (s) "
                "| RSS (MB) | Throughput (frag/s) |\n")
        f.write("|---------|---------|-----------|------|----------"
                "|----------|---------------------|\n")
        for r in results:
            for tool, elapsed_t in r.elapsed.items():
                rss = r.peak_rss.get(tool, 0)
                rss_str = f"{rss:.0f}" if rss > 0 else "\u2014"
                tp = r.throughput.get(tool, 0)
                tp_str = f"{tp:,.0f}" if tp > 0 else "\u2014"
                f.write(
                    f"| {r.dataset_name} | {r.aligner} | "
                    f"{r.alignment_elapsed:.1f} | {tool} | "
                    f"{elapsed_t:.1f} | {rss_str} | {tp_str} |\n"
                )

    logger.info("Wrote %s", path)


def write_per_transcript_csv(
    transcripts: list[Transcript],
    mrna_truth: dict[str, int],
    nrna_truth: dict[str, int],
    tool_tx_counts: dict[str, dict[str, float]],
    path: Path,
) -> None:
    """Per-transcript truth vs predicted CSV."""
    tools = sorted(tool_tx_counts.keys())
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        header = [
            "transcript_id", "gene_id", "gene_name",
            "mrna_abundance", "nrna_abundance",
            "mrna_truth", "nrna_truth",
        ] + tools
        writer.writerow(header)

        for t in transcripts:
            row = [
                t.t_id, t.g_id, t.g_name,
                f"{t.abundance or 0:.4f}", f"{t.nrna_abundance:.4f}",
                mrna_truth.get(t.t_id, 0),
                nrna_truth.get(t.t_id, 0),
            ]
            for tool in tools:
                row.append(f"{tool_tx_counts[tool].get(t.t_id, 0.0):.2f}")
            writer.writerow(row)
    logger.info("Wrote %s", path)


# ═══════════════════════════════════════════════════════════════════
# Tool naming helpers
# ═══════════════════════════════════════════════════════════════════


def _rigel_tool_name(
    aligner: AlignerConfig,
    rigel_config: HulkrnaConfig,
    multi_rigel: bool = False,
) -> str:
    if multi_rigel:
        return f"rigel_{rigel_config.name}_{aligner.name}"
    return f"rigel_{aligner.name}"


# ═══════════════════════════════════════════════════════════════════
# Main benchmark orchestration
# ═══════════════════════════════════════════════════════════════════


def run_benchmark(cfg: BenchmarkConfig) -> list[ConditionResult]:
    """Run the full benchmark.

    1. Resolve datasets (from sim_dir manifest and/or explicit datasets)
    2. Build indexes
    3. Per dataset × aligner: align → quantify → score
    """
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = Path(cfg.genome)
    gtf_path = Path(cfg.gtf)

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome not found: {genome_path}")
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    # ── 1. Resolve datasets ──────────────────────────────────────
    datasets: list[Dataset] = list(cfg.datasets)

    if cfg.sim_dir:
        manifest, sim_datasets = load_sim_manifest(cfg.sim_dir)
        datasets.extend(sim_datasets)
        # Inherit genome/gtf/filter from manifest if not overridden
        if not cfg.genome and manifest.get("genome"):
            genome_path = Path(manifest["genome"])
            cfg.genome = str(genome_path)
        if not cfg.gtf and manifest.get("gtf"):
            gtf_path = Path(manifest["gtf"])
            cfg.gtf = str(gtf_path)
        if cfg.transcript_filter == "all" and manifest.get("transcript_filter"):
            cfg.transcript_filter = manifest["transcript_filter"]

    if not datasets:
        raise ValueError("No datasets to benchmark — provide sim_dir or datasets in config")

    # ── 2. Load transcripts for scoring / FASTA writing ──────────
    transcripts = load_transcripts(
        gtf_path, transcript_filter=cfg.transcript_filter,
    )

    # ── 3. Build indexes ─────────────────────────────────────────
    rigel_index = build_rigel_index(
        genome_path, gtf_path, outdir / "rigel_index",
    )

    # Minimap2 index + splice BED (build-if-absent)
    has_minimap2 = any(a.type == "minimap2" for a in cfg.aligners)
    mm2_index_path: Path | None = None
    mm2_bed_path: Path | None = None
    if has_minimap2:
        mm2_index_path = build_minimap2_index(
            genome_path, outdir, cfg.threads, prebuilt=cfg.minimap2_index,
        )
        mm2_bed_path = build_minimap2_bed(
            gtf_path, outdir, prebuilt=cfg.minimap2_bed,
        )

    # BED12 for rigel (always needed)
    bed_path = outdir / "annotation.bed12"
    if not bed_path.exists():
        write_bed12(transcripts, bed_path)

    # Optional tool indexes
    transcript_fa = outdir / "transcripts.fa"
    salmon_idx = None
    kallisto_idx = None

    needs_transcript_fa = (
        (cfg.salmon_enabled and not cfg.salmon_index)
        or (cfg.kallisto_enabled and not cfg.kallisto_index)
    )
    if needs_transcript_fa and not transcript_fa.exists():
        write_transcript_fasta(cfg.genome, transcripts, transcript_fa)

    if cfg.salmon_enabled:
        salmon_idx = build_salmon_index(
            transcript_fa, outdir, cfg.threads, prebuilt=cfg.salmon_index,
        )
    if cfg.kallisto_enabled:
        kallisto_idx = build_kallisto_index(
            transcript_fa, outdir, prebuilt=cfg.kallisto_index,
        )

    # ── 4. Benchmark loop ────────────────────────────────────────
    multi_rigel = len(cfg.rigel_configs) > 1
    results: list[ConditionResult] = []

    for ds_idx, ds in enumerate(datasets):
        n_total = ds.n_rna + ds.n_gdna
        print(
            f"\n[{ds_idx + 1}/{len(datasets)}] {ds.name}: "
            f"RNA={ds.n_rna:,} gDNA={ds.n_gdna:,} SS={ds.strand_specificity:.2f}",
            flush=True,
        )

        ds_dir = outdir / ds.name
        ds_dir.mkdir(parents=True, exist_ok=True)

        # Parse truth if FASTQ available
        mrna_truth: dict[str, int] = {}
        nrna_truth: dict[str, int] = {}
        n_gdna_truth = 0
        has_truth = False

        if ds.fastq_r1 and ds.fastq_r1.exists():
            mrna_truth, nrna_truth, n_gdna_truth = parse_truth_from_fastq(ds.fastq_r1)
            has_truth = True

        n_mrna_truth = sum(mrna_truth.values())
        n_nrna_truth = sum(nrna_truth.values())

        # Per-aligner
        for ac in cfg.aligners:
            align_dir = ds_dir / f"align_{ac.name}"
            align_dir.mkdir(parents=True, exist_ok=True)
            bam_ns = align_dir / "reads_namesort.bam"

            # Produce BAM
            if bam_ns.exists():
                print(f"  BAM ({ac.name}) exists, skipping", flush=True)
            elif ac.type == "oracle":
                if ds.oracle_bam and ds.oracle_bam.exists():
                    # Symlink oracle BAM from sim output
                    bam_ns.symlink_to(ds.oracle_bam.resolve())
                    print(f"  Oracle BAM → {ds.oracle_bam}", flush=True)
                else:
                    print(f"  Oracle BAM not available for {ds.name}, skipping aligner", flush=True)
                    continue
            elif ac.type == "minimap2":
                if ds.fastq_r1 and ds.fastq_r2:
                    print(f"  minimap2...", end="", flush=True)
                    t0 = time.monotonic()
                    align_minimap2(
                        mm2_index_path or genome_path,
                        ds.fastq_r1, ds.fastq_r2, bam_ns,
                        bed_path=mm2_bed_path, threads=cfg.threads,
                    )
                    alignment_elapsed = time.monotonic() - t0
                    print(f" done ({alignment_elapsed:.1f}s)", flush=True)
                else:
                    print(f"  No FASTQ for minimap2, skipping", flush=True)
                    continue
            else:
                print(f"  Unsupported aligner: {ac.type}", flush=True)
                continue

            # Run tools
            tool_tx_counts: dict[str, dict[str, float]] = {}
            tool_pool: dict[str, dict[str, float]] = {}
            tool_elapsed: dict[str, float] = {}
            tool_peak_rss: dict[str, float] = {}
            tool_throughput: dict[str, float] = {}
            alignment_elapsed = 0.0

            for hc in cfg.rigel_configs:
                hn = _rigel_tool_name(ac, hc, multi_rigel)
                print(f"    rigel ({hn})...", end="", flush=True)
                try:
                    counts, elapsed, pools, n_frags = run_rigel_tool(
                        bam_ns, rigel_index, ds.pipeline_seed, hc,
                    )
                    tool_tx_counts[hn] = counts
                    tool_pool[hn] = pools
                    tool_elapsed[hn] = elapsed
                    tool_peak_rss[hn] = _get_peak_rss_mb()
                    tool_throughput[hn] = n_frags / elapsed if elapsed > 0 else 0.0
                    print(f" done ({elapsed:.1f}s, {tool_peak_rss[hn]:.0f}MB RSS, "
                          f"{tool_throughput[hn]:,.0f} frags/s)", flush=True)
                except Exception as e:
                    print(f" FAILED: {e}", flush=True)
                    if not cfg.keep_going:
                        raise

            # Salmon (once per dataset, on first aligner)
            if cfg.salmon_enabled and salmon_idx and ac == cfg.aligners[0] and ds.fastq_r1:
                salmon_dir = ds_dir / "salmon"
                print("    salmon...", end="", flush=True)
                try:
                    counts, elapsed = run_salmon(
                        salmon_idx, ds.fastq_r1, ds.fastq_r2, salmon_dir, cfg.threads,
                    )
                    tool_tx_counts["salmon"] = counts
                    tool_elapsed["salmon"] = elapsed
                    tool_peak_rss["salmon"] = 0.0
                    tool_throughput["salmon"] = n_total / elapsed if elapsed > 0 else 0.0
                    print(f" done ({elapsed:.1f}s)", flush=True)
                except Exception as e:
                    print(f" FAILED: {e}", flush=True)
                    if not cfg.keep_going:
                        raise

            # Kallisto
            if cfg.kallisto_enabled and kallisto_idx and ac == cfg.aligners[0] and ds.fastq_r1:
                kallisto_dir = ds_dir / "kallisto"
                print("    kallisto...", end="", flush=True)
                try:
                    counts, elapsed = run_kallisto(
                        kallisto_idx, ds.fastq_r1, ds.fastq_r2,
                        kallisto_dir, cfg.threads, ds.strand_specificity,
                    )
                    tool_tx_counts["kallisto"] = counts
                    tool_elapsed["kallisto"] = elapsed
                    tool_peak_rss["kallisto"] = 0.0
                    tool_throughput["kallisto"] = n_total / elapsed if elapsed > 0 else 0.0
                    print(f" done ({elapsed:.1f}s)", flush=True)
                except Exception as e:
                    print(f" FAILED: {e}", flush=True)
                    if not cfg.keep_going:
                        raise

            # HTSeq
            if cfg.htseq_enabled and ac == cfg.aligners[0]:
                htseq_dir = ds_dir / "htseq"
                bam_cs = align_dir / "reads_coordsort.bam"
                print("    htseq...", end="", flush=True)
                try:
                    if not bam_cs.exists():
                        coord_sort_bam(bam_ns, bam_cs, cfg.threads)
                    gene_counts, elapsed = run_htseq(
                        bam_cs, cfg.gtf, htseq_dir, ds.strand_specificity,
                    )
                    tool_elapsed["htseq"] = elapsed
                    tool_peak_rss["htseq"] = 0.0
                    tool_throughput["htseq"] = n_total / elapsed if elapsed > 0 else 0.0
                    print(f" done ({elapsed:.1f}s)", flush=True)
                except Exception as e:
                    print(f" FAILED: {e}", flush=True)
                    if not cfg.keep_going:
                        raise

            # Score + performance metrics
            tx_metrics: dict = {}
            gene_metrics: dict = {}

            mrna_truth_f = {k: float(v) for k, v in mrna_truth.items()} if has_truth else {}
            for tool, counts in tool_tx_counts.items():
                elapsed_t = tool_elapsed.get(tool, 0.0)
                rss_t = tool_peak_rss.get(tool, 0.0)
                tp_t = tool_throughput.get(tool, 0.0)
                if has_truth:
                    tm = score_tool(tool, mrna_truth_f, counts, elapsed_t, rss_t, tp_t)
                    tx_metrics[tool] = asdict(tm)
                    gene_truth = aggregate_to_genes(mrna_truth_f, transcripts)
                    gene_obs = aggregate_to_genes(counts, transcripts)
                    gm = score_tool(tool, gene_truth, gene_obs, elapsed_t, rss_t, tp_t)
                    gene_metrics[tool] = asdict(gm)
                else:
                    tx_metrics[tool] = asdict(ToolMetrics(
                        tool=tool, elapsed_sec=elapsed_t,
                        peak_rss_mb=rss_t, throughput_frags_per_sec=tp_t,
                        n_transcripts=0, total_truth=0, total_observed=0,
                        total_abs_error=0, mean_abs_error=0, rmse=0,
                        pearson=0, spearman=0,
                    ))

            # Per-transcript CSV
            if has_truth and tool_tx_counts:
                write_per_transcript_csv(
                    transcripts, mrna_truth, nrna_truth, tool_tx_counts,
                    ds_dir / f"per_transcript_counts_{ac.name}.csv",
                )

            results.append(ConditionResult(
                dataset_name=ds.name,
                gdna_label=ds.gdna_label,
                gdna_rate=ds.gdna_rate,
                n_gdna=ds.n_gdna,
                strand_specificity=ds.strand_specificity,
                n_rna_fragments=ds.n_rna,
                n_total_fragments=n_total,
                n_gdna_truth=n_gdna_truth,
                n_mrna_truth=n_mrna_truth,
                n_nrna_truth=n_nrna_truth,
                aligner=ac.name,
                transcript_metrics=tx_metrics,
                gene_metrics=gene_metrics,
                pool_counts={k: v for k, v in tool_pool.items()},
                elapsed=tool_elapsed,
                peak_rss=tool_peak_rss,
                throughput=tool_throughput,
                alignment_elapsed=alignment_elapsed,
            ))

    return results


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Benchmark rigel (and optional tools) on simulated or real data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", required=True, help="YAML configuration file")
    p.add_argument("--sim-dir", help="sim.py output directory (overrides YAML)")
    p.add_argument("--genome", help="Genome FASTA (overrides YAML)")
    p.add_argument("--gtf", help="Gene annotation GTF (overrides YAML)")
    p.add_argument("--outdir", help="Output directory (overrides YAML)")
    p.add_argument("--threads", type=int, help="Number of threads (overrides YAML)")
    p.add_argument("--verbose", action="store_true", default=None, help="Verbose logging")
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    cfg = parse_yaml_config(args.config)

    # CLI overrides
    if args.sim_dir:
        cfg.sim_dir = args.sim_dir
    if args.genome:
        cfg.genome = args.genome
    if args.gtf:
        cfg.gtf = args.gtf
    if args.outdir:
        cfg.outdir = args.outdir
    if args.threads is not None:
        cfg.threads = args.threads
    if args.verbose is not None:
        cfg.verbose = args.verbose

    level = logging.DEBUG if cfg.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    if not cfg.genome and not cfg.sim_dir:
        print("Error: genome FASTA not specified (provide --genome or --sim-dir)", file=sys.stderr)
        return 1
    if not cfg.gtf and not cfg.sim_dir:
        print("Error: GTF not specified (provide --gtf or --sim-dir)", file=sys.stderr)
        return 1

    print("Benchmark", flush=True)
    print(f"  Genome:          {cfg.genome}", flush=True)
    print(f"  GTF:             {cfg.gtf}", flush=True)
    print(f"  Output:          {cfg.outdir}", flush=True)
    print(f"  Threads:         {cfg.threads}", flush=True)
    if cfg.sim_dir:
        print(f"  Sim dir:         {cfg.sim_dir}", flush=True)
    print(f"  Aligners:        {[a.name for a in cfg.aligners]}", flush=True)
    print(f"  HulkRNA configs: {[h.name for h in cfg.rigel_configs]}", flush=True)
    print(f"  Salmon:          {cfg.salmon_enabled}"
          f"{' (index: ' + cfg.salmon_index + ')' if cfg.salmon_index else ''}",
          flush=True)
    print(f"  Kallisto:        {cfg.kallisto_enabled}"
          f"{' (index: ' + cfg.kallisto_index + ')' if cfg.kallisto_index else ''}",
          flush=True)
    print(f"  HTSeq:           {cfg.htseq_enabled}", flush=True)
    if cfg.minimap2_index:
        print(f"  minimap2 index:  {cfg.minimap2_index}", flush=True)
    if cfg.minimap2_bed:
        print(f"  minimap2 BED:    {cfg.minimap2_bed}", flush=True)

    t0 = time.monotonic()
    results = run_benchmark(cfg)
    total_elapsed = time.monotonic() - t0

    outdir = Path(cfg.outdir)
    write_summary_json(results, outdir / "summary.json")
    write_summary_csv(results, outdir / "summary.csv")
    write_summary_md(results, outdir / "summary.md")

    print(f"\nBenchmark complete in {total_elapsed:.1f}s", flush=True)
    print(f"Results: {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
