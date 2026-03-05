#!/usr/bin/env python3
"""Regional benchmark: hulkrna multimap vs salmon vs kallisto vs htseq.

Supports a combinatorial grid of gDNA contamination rates and strand
specificity values. Per-region setup (extract, abundance assignment, index
build) is done once; the condition sweep re-simulates reads and re-runs all
tools for each (gDNA rate, strand specificity) pair.

Multiple aligners can be evaluated simultaneously (e.g. ``--aligner minimap2
oracle``).  Each aligner produces a separate hulkrna result tagged as
``hulkrna_<aligner>`` (e.g. ``hulkrna_minimap2``, ``hulkrna_oracle``).
FASTQ-based tools (salmon, kallisto) are run once per condition and shared
across aligners.

YAML configuration
------------------
All settings can be specified in a YAML config file (``--config``).

**Named regions** – region entries can be plain coordinate strings or
single-key dicts that supply a human-readable label::

    region:
      - HBB: chr11:5225000-5310000
      - chr7:55019017-55211628        # auto-labelled from coordinates

**hulkrna parameterizations** – define multiple configurations to benchmark
hulkrna against itself (e.g. different EM thresholds)::

    hulkrna_configs:
      default: {}
      strict_em:
        em_convergence_delta: 1.0e-5

When there are multiple configs, tool names become
``hulkrna_<config>_<aligner>`` (e.g. ``hulkrna_strict_em_oracle``).

**Aligner parameterizations** – define named aligner configurations with
explicit parameters (replaces the ``aligner:`` list)::

    aligners:
      oracle: {}
      hisat2_plain:
        type: hisat2
        index_type: plain
      hisat2_graph:
        type: hisat2
        index_type: graph
        splice_sites_at_build: true
        exons_at_build: true

Tools
-----
- hulkrna_<aligner>        : one result per aligner (single hulkrna config)
- hulkrna_<cfg>_<aligner>  : one result per config × aligner (multiple configs)
- salmon                   : salmon quant --validateMappings
- kallisto                 : kallisto quant (--rf-stranded when strand_spec >= 0.9)
- htseq_<aligner>          : htseq-count (gene-level only), one result per aligner

Genome aligners for hulkrna/htseq inputs
----------------------------------------
- minimap2      : splice:sr with secondary alignments + annotation BED12
- hisat2        : hisat2/hisat2-build with secondary alignments
- oracle        : perfect BAM bypassing alignment entirely

gDNA model
----------
gDNA abundance is computed from total transcript abundance:
    gDNA_abundance = sum(transcript_abundances) * gDNA_rate
where gDNA_rate is in [0, 1].

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
        <condition>/          e.g. gdna_r0.25_ss_1.00
            reads/
            align/
            salmon/
            kallisto/
            htseq/
            per_transcript_counts.csv
            per_gene_counts.csv
            per_pool_counts.csv
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
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import yaml
except ImportError:  # pragma: no cover - handled at runtime when config is used
    yaml = None

from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig, FragmentScoringConfig
from hulkrna.gtf import GTFRecord
from hulkrna.index import TranscriptIndex, write_bed12
from hulkrna.pipeline import run_pipeline
from hulkrna.scoring import (
    GDNA_SPLICE_PENALTIES,
    SPLICE_UNANNOT,
    overhang_alpha_to_log_penalty,
)
from hulkrna.sim import OracleBamSimulator, GDNAConfig, ReadSimulator, SimConfig, reverse_complement
from hulkrna.transcript import Transcript
from hulkrna.types import Interval, Strand

logger = logging.getLogger(__name__)

# ── Tool identifiers ────────────────────────────────────────────────

# Legacy constants kept for aggregate_benchmarks backward-compatibility
TRANSCRIPT_TOOLS = ("hulkrna_mm", "salmon", "kallisto")
GENE_TOOLS = ("hulkrna_mm", "salmon", "kallisto", "htseq")

ALIGNER_CHOICES = ("minimap2", "hisat2", "oracle")


# ── Configuration dataclasses ───────────────────────────────────────


@dataclass
class HulkrnaConfig:
    """Named hulkrna parameterization for benchmarking.

    Parameters in *params* are forwarded as keyword arguments to
    ``run_pipeline()``.  Anything not listed is left at its
    pipeline default.
    """

    name: str
    params: dict = field(default_factory=dict)


# Default minimap2 / hisat2 parameter dicts used when the YAML does
# not provide an explicit ``aligners:`` section.

MINIMAP2_DEFAULTS: dict = {
    "preset": "splice:sr",
    "secondary": True,
    "bed_guided": True,
}

HISAT2_DEFAULTS: dict = {
    "k": 50,
    "secondary": True,
    "no_unal": True,
    "index_type": "plain",
    "splice_sites_at_align": True,
    "splice_sites_at_build": False,
    "exons_at_build": False,
}


@dataclass
class AlignerConfig:
    """Named aligner configuration.

    *type* must be one of ``ALIGNER_CHOICES``.  *params* holds
    aligner-specific knobs (preset, k, index_type, …).
    """

    name: str
    type: str  # "oracle", "minimap2", "hisat2"
    params: dict = field(default_factory=dict)


# ── Tool-naming helpers ─────────────────────────────────────────────


def _hulkrna_tool_name(
    aligner: str | AlignerConfig,
    hulkrna_config: HulkrnaConfig | None = None,
    *,
    multi_hulk: bool = False,
) -> str:
    """Return hulkrna tool identifier.

    When there is only one hulkrna config (``multi_hulk=False``), the
    name is ``hulkrna_<aligner>``.  With multiple configs it becomes
    ``hulkrna_<config>_<aligner>`` so results are distinguishable.
    """
    aname = aligner.name if isinstance(aligner, AlignerConfig) else aligner
    if multi_hulk and hulkrna_config is not None:
        return f"hulkrna_{hulkrna_config.name}_{aname}"
    return f"hulkrna_{aname}"


def _htseq_tool_name(aligner: str | AlignerConfig) -> str:
    """Return htseq tool identifier for a given aligner."""
    aname = aligner.name if isinstance(aligner, AlignerConfig) else aligner
    return f"htseq_{aname}"


def _get_transcript_tools(
    aligner_configs: list[AlignerConfig],
    hulkrna_configs: list[HulkrnaConfig],
) -> tuple[str, ...]:
    """Build transcript-level tool list from active configs."""
    multi = len(hulkrna_configs) > 1
    tools: list[str] = []
    for ac in aligner_configs:
        for hc in hulkrna_configs:
            tools.append(_hulkrna_tool_name(ac, hc, multi_hulk=multi))
    tools.extend(["salmon", "kallisto"])
    return tuple(tools)


def _get_gene_tools(
    aligner_configs: list[AlignerConfig],
    hulkrna_configs: list[HulkrnaConfig],
    include_htseq: bool = False,
) -> tuple[str, ...]:
    """Build gene-level tool list from active configs."""
    multi = len(hulkrna_configs) > 1
    tools: list[str] = []
    for ac in aligner_configs:
        for hc in hulkrna_configs:
            tools.append(_hulkrna_tool_name(ac, hc, multi_hulk=multi))
    tools.extend(["salmon", "kallisto"])
    if include_htseq:
        for ac in aligner_configs:
            tools.append(_htseq_tool_name(ac))
    return tuple(tools)

# ── Default gDNA levels and strand specificities ────────────────────

DEFAULT_STRAND_SPECIFICITIES = (0.95, 0.99, 1.0)
DEFAULT_GDNA_RATES = (0.0, 0.1, 0.25, 0.5)
DEFAULT_NRNA_RATES = (0.0,)


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
    gdna_rate: float
    gdna_abundance: float
    nrna_label: str
    nrna_rate: float
    nrna_abundance: float
    strand_specificity: float
    aligner: str
    transcript_metrics: dict  # tool name → ToolMetrics (as dict)
    gene_metrics: dict  # tool name → ToolMetrics (as dict)
    pool_metrics: list[dict]  # rows with pool-level truth/observed/error per tool


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


def parse_region(region_str: str, name: str | None = None) -> RegionSpec:
    """Parse a ``chr:start-end`` string into a :class:`RegionSpec`.

    *name*, when given, is used as the human-friendly label; otherwise
    a label is auto-generated from the coordinates.
    """
    m = re.match(r"(.+):(\d+)-(\d+)$", region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}")
    chrom = m.group(1)
    start1 = int(m.group(2))
    end1 = int(m.group(3))
    label = name if name else f"{chrom}_{start1}_{end1}"
    return RegionSpec(label=label, chrom=chrom, start0=start1 - 1, end0=end1)


def _parse_region_entry(entry) -> RegionSpec:
    """Parse a single region entry (string or ``{name: coord}`` dict)."""
    if isinstance(entry, str):
        return parse_region(entry)
    if isinstance(entry, dict):
        if len(entry) != 1:
            raise ValueError(
                f"Named region dict must have exactly one key: value pair, got {entry}"
            )
        name, coord = next(iter(entry.items()))
        return parse_region(str(coord), name=str(name))
    raise TypeError(f"Unsupported region entry type: {type(entry).__name__}")


def parse_regions(args: argparse.Namespace) -> list[RegionSpec]:
    regions: list[RegionSpec] = []
    for r in args.region:
        regions.append(_parse_region_entry(r))
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
                gtf_obj = GTFRecord(
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
    *,
    aligner_config: AlignerConfig | None = None,
) -> None:
    """Align PE reads with minimap2 and produce a name-sorted BAM.

    When *aligner_config* is provided its ``params`` dict may contain:

    - **preset** (str) – minimap2 preset, e.g. ``splice:sr`` (default).
    - **secondary** (bool) – emit secondary alignments (default True).
    - **bed_guided** (bool) – use BED12 annotation for splice guidance
      (default True).  When False, *bed_path* is ignored.
    """
    params = dict(MINIMAP2_DEFAULTS)
    if aligner_config is not None:
        params.update(aligner_config.params)

    preset = params.get("preset", "splice:sr")
    secondary = params.get("secondary", True)
    bed_guided = params.get("bed_guided", True)

    minimap2_cmd = [
        "minimap2",
        "-ax",
        str(preset),
    ]
    if secondary:
        minimap2_cmd.append("--secondary=yes")
    if bed_guided and bed_path is not None:
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


def write_hisat2_splice_sites(
    transcripts: list[Transcript],
    out_path: Path,
) -> Path:
    """Write a HISAT2 splice-site file from transcript exon coordinates.

    Format (tab-separated, 0-based)::

        chrom  left_flanking_base  right_flanking_base  strand

    where *left_flanking_base* is the last base of the upstream exon and
    *right_flanking_base* is the first base of the downstream exon.
    Coordinates use the 0-based convention expected by ``hisat2-build --ss``
    and ``hisat2 --known-splicesite-infile``.
    """
    seen: set[tuple[str, int, int, str]] = set()
    with open(out_path, "w") as fh:
        for t in transcripts:
            ref = t.ref or ""
            strand_str = "+" if t.strand == Strand.POS else "-"
            for i in range(1, len(t.exons)):
                left = t.exons[i - 1].end - 1   # last base of left exon (0-based)
                right = t.exons[i].start         # first base of right exon (0-based)
                key = (ref, left, right, strand_str)
                if key not in seen:
                    seen.add(key)
                    fh.write(f"{ref}\t{left}\t{right}\t{strand_str}\n")
    logger.info("Wrote %d unique splice sites to %s", len(seen), out_path)
    return out_path


def write_hisat2_exons(
    transcripts: list[Transcript],
    out_path: Path,
) -> Path:
    """Write a HISAT2 exon file from transcript exon coordinates.

    Format (tab-separated, 0-based inclusive)::

        chrom  left_position  right_position

    Coordinates are 0-based, inclusive on both ends, matching the format
    expected by ``hisat2-build --exon``.
    """
    seen: set[tuple[str, int, int]] = set()
    with open(out_path, "w") as fh:
        for t in transcripts:
            ref = t.ref or ""
            for e in t.exons:
                key = (ref, e.start, e.end - 1)
                if key not in seen:
                    seen.add(key)
                    fh.write(f"{ref}\t{e.start}\t{e.end - 1}\n")
    logger.info("Wrote %d unique exons to %s", len(seen), out_path)
    return out_path


def align_hisat2(
    region_fa: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    out_bam: Path,
    *,
    threads: int,
    hisat2_exe: str,
    hisat2_build_exe: str,
    hisat2_index_dir: Path,
    splice_sites_file: Path | None = None,
    exons_file: Path | None = None,
    aligner_config: AlignerConfig | None = None,
) -> None:
    """Align PE reads with HISAT2 and produce a name-sorted BAM.

    When *aligner_config* is provided its ``params`` dict may contain:

    - **k** (int) – max alignments per read (default 50).
    - **secondary** (bool) – emit secondary alignments (default True).
    - **no_unal** (bool) – suppress unaligned reads (default True).
    - **index_type** (``"plain"`` | ``"graph"``) – build a plain BWT
      index (default, avoids graph-bias) or an HGFM graph index with
      ``--ss`` / ``--exon`` annotations baked in.
    - **splice_sites_at_align** (bool) – pass known splice-sites at
      alignment time via ``--known-splicesite-infile`` (default True).
    - **splice_sites_at_build** (bool) – pass ``--ss`` to
      ``hisat2-build`` when building a graph index (default False;
      only used when index_type=graph).
    - **exons_at_build** (bool) – pass ``--exon`` to ``hisat2-build``
      when building a graph index (default False; only used when
      index_type=graph).
    """
    params = dict(HISAT2_DEFAULTS)
    if aligner_config is not None:
        params.update(aligner_config.params)

    k = int(params.get("k", 50))
    secondary = params.get("secondary", True)
    no_unal = params.get("no_unal", True)
    index_type = str(params.get("index_type", "plain"))
    ss_at_align = params.get("splice_sites_at_align", True)
    ss_at_build = params.get("splice_sites_at_build", False)
    exons_at_build = params.get("exons_at_build", False)

    # Use a sub-directory per index_type so plain and graph indexes can
    # coexist when running multiple hisat2 configs on the same region.
    idx_subdir = hisat2_index_dir / index_type
    idx_subdir.mkdir(parents=True, exist_ok=True)
    index_prefix = idx_subdir / "idx"

    if not (idx_subdir / "idx.1.ht2").exists():
        build_cmd = [
            hisat2_build_exe,
            "-p", str(max(1, threads)),
        ]
        if index_type == "graph":
            if ss_at_build and splice_sites_file is not None:
                build_cmd.extend(["--ss", str(splice_sites_file)])
            if exons_at_build and exons_file is not None:
                build_cmd.extend(["--exon", str(exons_file)])
            logger.info(
                "Building HISAT2 index (graph): %s",
                " ".join(build_cmd + [str(region_fa), str(index_prefix)]),
            )
        else:
            # Plain BWT index — avoids graph-bias that suppresses long-
            # intron junctions in paralogous-gene regions (66× improvement
            # for ENST642908 in the HBB locus).
            logger.info(
                "Building HISAT2 index (plain): %s",
                " ".join(build_cmd + [str(region_fa), str(index_prefix)]),
            )
        build_cmd.extend([str(region_fa), str(index_prefix)])
        subprocess.run(build_cmd, check=True, capture_output=True, text=True)

    hisat2_cmd = [
        hisat2_exe,
        "-x",
        str(index_prefix),
        "-1",
        str(fastq_r1),
        "-2",
        str(fastq_r2),
        "-p",
        str(max(1, threads)),
        "-k",
        str(k),
    ]
    if secondary:
        hisat2_cmd.append("--secondary")
    if no_unal:
        hisat2_cmd.append("--no-unal")

    # Provide known splice sites at alignment time for maximum sensitivity
    if ss_at_align and splice_sites_file is not None:
        hisat2_cmd.extend(["--known-splicesite-infile", str(splice_sites_file)])
    sort_cmd = ["samtools", "sort", "-n", "-o", str(out_bam), "-"]

    p1 = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
        logger.warning("hisat2 returned rc=%d: %s", p1.returncode, p1_stderr.decode())


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


def parse_truth_from_bam(bam_path: Path) -> tuple[dict[str, int], int]:
    """Parse ground-truth counts from read names in a BAM file.

    Read-name format is identical to the FASTQ simulator:
    ``{tid}:{start}-{end}:{strand}:{idx}`` for RNA, ``gdna:...`` for gDNA,
    ``nrna_{tid}:...`` for nascent RNA.

    Only counts R1 reads (flag 0x40) to avoid double-counting.
    """
    import pysam

    counts: Counter[str] = Counter()
    n_gdna = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not (read.flag & 0x40):  # only count R1 (first in pair)
                continue
            qname = read.query_name
            tid = qname.split(":", 1)[0]
            if tid == "gdna":
                n_gdna += 1
            else:
                counts[tid] += 1
    return dict(counts), n_gdna


def compute_truth_pool_counts(
    truth_counts: dict[str, int],
    n_gdna_actual: int,
) -> dict[str, float]:
    mature_truth = float(
        sum(count for tid, count in truth_counts.items() if not str(tid).startswith("nrna_"))
    )
    nascent_truth = float(
        sum(count for tid, count in truth_counts.items() if str(tid).startswith("nrna_"))
    )
    genomic_truth = float(n_gdna_actual)
    return {
        "mature_rna": mature_truth,
        "nascent_rna": nascent_truth,
        "genomic_dna": genomic_truth,
        "rna": mature_truth + nascent_truth,
        "dna": genomic_truth,
    }


def compute_pool_metric_rows(
    tool: str,
    pool_truth: dict[str, float],
    pool_observed: dict[str, float],
) -> list[dict]:
    rows: list[dict] = []
    for pool in ("mature_rna", "nascent_rna", "genomic_dna", "rna", "dna"):
        truth = float(pool_truth.get(pool, 0.0))
        observed = float(pool_observed.get(pool, 0.0))
        signed_error = observed - truth
        rows.append({
            "tool": tool,
            "pool": pool,
            "truth": round(truth, 6),
            "observed": round(observed, 6),
            "signed_error": round(signed_error, 6),
            "abs_error": round(abs(signed_error), 6),
        })
    return rows


# ── Annotated BAM truth comparison ──────────────────────────────────

# Verdict codes for per-read accuracy.
VERDICT_CORRECT_TRANSCRIPT = "correct_transcript"
VERDICT_CORRECT_GENE = "correct_gene"
VERDICT_CORRECT_POOL = "correct_pool"
VERDICT_WRONG_POOL = "wrong_pool"
VERDICT_INTERGENIC = "intergenic"

# Mapping from ZP tag values to canonical truth pool names.
_POOL_NORM = {
    "mRNA": "mRNA",
    "nRNA": "nRNA",
    "gDNA": "gDNA",
    "intergenic": "intergenic",
    "chimeric": "chimeric",
}


def _parse_truth_from_qname(qname: str) -> tuple[str, str]:
    """Extract ground-truth transcript ID and pool from a read name.

    Read-name format (from the simulator):
      ``{tid}:{start}-{end}:{strand}:{idx}``

    Returns ``(truth_tid, truth_pool)`` where:
    - For gDNA reads: ``("gdna", "gDNA")``
    - For nRNA reads: ``("nrna_ENST...", "nRNA")``
    - For mRNA reads: ``("ENST...", "mRNA")``
    """
    tid = qname.split(":", 1)[0]
    if tid == "gdna":
        return "gdna", "gDNA"
    if tid.startswith("nrna_"):
        return tid, "nRNA"
    return tid, "mRNA"


@dataclass
class ReadAccuracyRow:
    """Per-read accuracy record for a single fragment."""

    qname: str
    truth_tid: str
    truth_pool: str
    assigned_tid: str
    assigned_gid: str
    assigned_pool: str
    posterior: float
    n_candidates: int
    frag_class: str
    splice_type: str
    verdict: str


@dataclass
class ReadAccuracySummary:
    """Aggregate read-level accuracy statistics."""

    total: int = 0
    correct_transcript: int = 0
    correct_gene: int = 0
    correct_pool: int = 0
    wrong_pool: int = 0
    intergenic: int = 0

    # Per truth-pool breakdown.
    by_truth_pool: dict = field(default_factory=lambda: defaultdict(Counter))

    def add(self, truth_pool: str, verdict: str) -> None:
        self.total += 1
        if verdict == VERDICT_CORRECT_TRANSCRIPT:
            self.correct_transcript += 1
        elif verdict == VERDICT_CORRECT_GENE:
            self.correct_gene += 1
        elif verdict == VERDICT_CORRECT_POOL:
            self.correct_pool += 1
        elif verdict == VERDICT_WRONG_POOL:
            self.wrong_pool += 1
        elif verdict == VERDICT_INTERGENIC:
            self.intergenic += 1
        self.by_truth_pool[truth_pool][verdict] += 1

    def pct(self, n: int) -> float:
        return 100.0 * n / self.total if self.total > 0 else 0.0

    def summary_rows(self) -> list[dict]:
        """Return tidy rows for CSV output."""
        rows: list[dict] = []
        for pool in ("mRNA", "nRNA", "gDNA"):
            cts = self.by_truth_pool.get(pool, Counter())
            pool_total = sum(cts.values())
            for v in (
                VERDICT_CORRECT_TRANSCRIPT,
                VERDICT_CORRECT_GENE,
                VERDICT_CORRECT_POOL,
                VERDICT_WRONG_POOL,
                VERDICT_INTERGENIC,
            ):
                count = cts.get(v, 0)
                rows.append({
                    "truth_pool": pool,
                    "verdict": v,
                    "count": count,
                    "pct": 100.0 * count / pool_total if pool_total > 0 else 0.0,
                    "pool_total": pool_total,
                })
        # Overall row.
        for v in (
            VERDICT_CORRECT_TRANSCRIPT,
            VERDICT_CORRECT_GENE,
            VERDICT_CORRECT_POOL,
            VERDICT_WRONG_POOL,
            VERDICT_INTERGENIC,
        ):
            n = getattr(self, v, 0)
            rows.append({
                "truth_pool": "ALL",
                "verdict": v,
                "count": n,
                "pct": self.pct(n),
                "pool_total": self.total,
            })
        return rows


def compare_annotated_bam_to_truth(
    annotated_bam_path: Path,
    gene_map: dict[str, str],
) -> tuple[list[ReadAccuracyRow], ReadAccuracySummary]:
    """Compare annotated BAM assignments against simulated read-name truth.

    Parameters
    ----------
    annotated_bam_path
        Path to the annotated BAM produced by ``run_pipeline`` with
        ``annotated_bam_path`` set.
    gene_map
        Mapping ``{transcript_id: gene_id}`` from the reference index.

    Returns
    -------
    tuple[list[ReadAccuracyRow], ReadAccuracySummary]
        Per-read rows and aggregate summary.
    """
    import pysam

    rows: list[ReadAccuracyRow] = []
    summary = ReadAccuracySummary()

    seen_qnames: set[str] = set()

    with pysam.AlignmentFile(str(annotated_bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            # Only count R1 (first in pair) to avoid double-counting.
            if not (rec.flag & 0x40):
                continue
            # Only count the primary hit (ZH=1).
            try:
                zh = rec.get_tag("ZH")
            except KeyError:
                continue
            if zh != 1:
                continue
            # De-duplicate by qname (in case of supplementary records).
            qname = rec.query_name
            if qname in seen_qnames:
                continue
            seen_qnames.add(qname)

            truth_tid, truth_pool = _parse_truth_from_qname(qname)

            # Read assignment tags.
            assigned_tid = str(rec.get_tag("ZT"))
            assigned_gid = str(rec.get_tag("ZG"))
            assigned_pool = str(rec.get_tag("ZP"))
            posterior = float(rec.get_tag("ZW"))
            n_candidates = int(rec.get_tag("ZN"))
            frag_class = str(rec.get_tag("ZC"))
            splice_type = str(rec.get_tag("ZS"))

            # Classify the assignment.
            verdict = _classify_verdict(
                truth_tid, truth_pool, assigned_tid, assigned_gid,
                assigned_pool, gene_map,
            )

            summary.add(truth_pool, verdict)

            rows.append(ReadAccuracyRow(
                qname=qname,
                truth_tid=truth_tid,
                truth_pool=truth_pool,
                assigned_tid=assigned_tid,
                assigned_gid=assigned_gid,
                assigned_pool=assigned_pool,
                posterior=posterior,
                n_candidates=n_candidates,
                frag_class=frag_class,
                splice_type=splice_type,
                verdict=verdict,
            ))

    return rows, summary


def _classify_verdict(
    truth_tid: str,
    truth_pool: str,
    assigned_tid: str,
    assigned_gid: str,
    assigned_pool: str,
    gene_map: dict[str, str],
) -> str:
    """Classify a single read as correct/error.

    Verdict hierarchy (best → worst):
    1. ``correct_transcript`` — pool matches AND transcript matches
    2. ``correct_gene`` — pool mRNA/nRNA AND gene matches but
       transcript differs (isoform confusion)
    3. ``correct_pool`` — pool matches but gene/transcript wrong
    4. ``intergenic`` — assigned to intergenic when truth was RNA/gDNA
    5. ``wrong_pool`` — assigned to a different active pool
    """
    norm_pool = _POOL_NORM.get(assigned_pool, assigned_pool)

    # Intergenic assignment is its own category.
    if norm_pool == "intergenic":
        return VERDICT_INTERGENIC

    # --- gDNA reads ---
    if truth_pool == "gDNA":
        if norm_pool == "gDNA":
            # gDNA has no specific transcript; pool match is the best
            # we can do.
            return VERDICT_CORRECT_TRANSCRIPT
        return VERDICT_WRONG_POOL

    # --- nRNA reads ---
    if truth_pool == "nRNA":
        if norm_pool == "nRNA":
            # nRNA truth_tid is "nrna_ENST…"; assigned_tid is
            # the bare transcript ID ("ENST…").
            bare_tid = truth_tid[5:] if truth_tid.startswith("nrna_") else truth_tid
            if assigned_tid == bare_tid:
                return VERDICT_CORRECT_TRANSCRIPT
            truth_gid = gene_map.get(bare_tid, "")
            assigned_gene = gene_map.get(assigned_tid, assigned_gid)
            if truth_gid and truth_gid == assigned_gene:
                return VERDICT_CORRECT_GENE
            return VERDICT_CORRECT_POOL
        return VERDICT_WRONG_POOL

    # --- mRNA reads ---
    if truth_pool == "mRNA":
        if norm_pool == "mRNA":
            if assigned_tid == truth_tid:
                return VERDICT_CORRECT_TRANSCRIPT
            truth_gid = gene_map.get(truth_tid, "")
            assigned_gene = gene_map.get(assigned_tid, assigned_gid)
            if truth_gid and truth_gid == assigned_gene:
                return VERDICT_CORRECT_GENE
            return VERDICT_CORRECT_POOL
        return VERDICT_WRONG_POOL

    # Fallback for unexpected pools.
    return VERDICT_WRONG_POOL


def write_read_accuracy_csv(
    rows: list[ReadAccuracyRow],
    out_path: Path,
) -> None:
    """Write per-read accuracy detail CSV."""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "qname", "truth_tid", "truth_pool",
            "assigned_tid", "assigned_gid", "assigned_pool",
            "posterior", "n_candidates", "frag_class", "splice_type",
            "verdict",
        ])
        for r in rows:
            writer.writerow([
                r.qname, r.truth_tid, r.truth_pool,
                r.assigned_tid, r.assigned_gid, r.assigned_pool,
                f"{r.posterior:.6f}", r.n_candidates, r.frag_class,
                r.splice_type, r.verdict,
            ])


def write_read_accuracy_summary_csv(
    summary: ReadAccuracySummary,
    out_path: Path,
) -> None:
    """Write aggregate read accuracy summary CSV."""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["truth_pool", "verdict", "count", "pct", "pool_total"])
        for row in summary.summary_rows():
            writer.writerow([
                row["truth_pool"],
                row["verdict"],
                row["count"],
                f"{row['pct']:.2f}",
                row["pool_total"],
            ])


# ── Tool runners ────────────────────────────────────────────────────


def _build_pipeline_config(
    args: argparse.Namespace,
    hulkrna_config: HulkrnaConfig | None = None,
    annotated_bam_path: Path | None = None,
) -> PipelineConfig:
    """Build a :class:`PipelineConfig` from benchmark args + overrides."""
    # ── Collect raw param dict (base from args, then overlay) ────
    raw: dict = {}
    if hulkrna_config is not None:
        raw.update(hulkrna_config.params)

    # ── EM config ────────────────────────────────────────────────
    em_kw: dict = {"seed": args.pipeline_seed}
    _EM_ALIASES = {
        "em_convergence_delta": "convergence_delta",
        "em_iterations": "iterations",
        "em_prior_alpha": "prior_alpha",
        "em_pseudocount": "prior_alpha",
        "em_prior_gamma": "prior_gamma",
        "em_mode": "mode",
        "prune_threshold": "prune_threshold",
        "confidence_threshold": "confidence_threshold",
    }
    for raw_key, cfg_key in _EM_ALIASES.items():
        if raw_key in raw:
            em_kw[cfg_key] = raw.pop(raw_key)

    # ── Scoring config ───────────────────────────────────────────
    scoring_kw: dict = {}
    # overhang_alpha — from override or from top-level args
    ov_alpha = raw.pop("overhang_alpha", None) or getattr(args, "overhang_alpha", None)
    if ov_alpha is not None:
        scoring_kw["overhang_log_penalty"] = overhang_alpha_to_log_penalty(ov_alpha)

    mm_alpha = raw.pop("mismatch_alpha", None)
    if mm_alpha is not None:
        scoring_kw["mismatch_log_penalty"] = overhang_alpha_to_log_penalty(mm_alpha)

    gdna_pen = raw.pop("gdna_splice_penalty_unannot", None)
    if gdna_pen is None:
        gdna_pen = getattr(args, "gdna_splice_penalty_unannot", None)
    if gdna_pen is not None:
        penalties = dict(GDNA_SPLICE_PENALTIES)
        penalties[SPLICE_UNANNOT] = gdna_pen
        scoring_kw["gdna_splice_penalties"] = penalties

    # ── Scan config ──────────────────────────────────────────────
    scan_kw: dict = {
        "sj_strand_tag": "auto",
        "include_multimap": True,
    }

    return PipelineConfig(
        em=EMConfig(**em_kw),
        scan=BamScanConfig(**scan_kw),
        scoring=FragmentScoringConfig(**scoring_kw),
        annotated_bam_path=annotated_bam_path,
    )


def run_hulkrna_tool(
    bam_path: Path,
    index: TranscriptIndex,
    args: argparse.Namespace,
    hulkrna_config: HulkrnaConfig | None = None,
    annotated_bam_path: Path | None = None,
) -> tuple[dict[str, float], float, dict[str, float]]:
    """Run hulkrna pipeline and return (transcript_counts, elapsed_sec, pool_counts).

    When *hulkrna_config* is provided, its ``params`` are merged on
    top of the base arguments drawn from *args*.  Supported params:
    ``em_convergence_delta``, ``em_iterations``, ``em_prior_alpha``,
    ``em_prior_gamma``,
    ``overhang_alpha``, ``mismatch_alpha``,
    ``gdna_splice_penalty_unannot``, etc.

    When *annotated_bam_path* is provided, an annotated BAM with
    per-fragment assignment tags (ZT, ZG, ZP, ZW, …) is written to
    that path.
    """
    t0 = time.monotonic()
    cfg = _build_pipeline_config(args, hulkrna_config, annotated_bam_path)
    pipe = run_pipeline(bam_path, index, config=cfg)
    elapsed = time.monotonic() - t0

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
        "rna": mature_pred + nascent_pred,
        "dna": genomic_pred,
    }
    return counts, elapsed, pool_counts


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
    requested_rates: tuple[float, ...],
    requested_labels: tuple[str, ...] | None = None,
) -> dict[str, tuple[float, float]]:
    """Compute gDNA abundance from total transcript abundance and gDNA rates.

    gDNA abundance is defined as:
      sum(transcript abundances) * gdna_rate

    Returns dict[label] = (rate, abundance).
    """
    total_tx_abundance = float(sum(float(t.abundance) for t in transcripts))
    levels: dict[str, tuple[float, float]] = {}
    for idx, rate in enumerate(requested_rates):
        if rate < 0.0 or rate > 1.0:
            raise ValueError(f"Invalid gDNA rate {rate}; expected in [0.0, 1.0]")
        if requested_labels is not None and idx < len(requested_labels):
            label = requested_labels[idx].strip()
        else:
            label = f"r{rate:g}"
        abundance = total_tx_abundance * float(rate)
        levels[label] = (float(rate), float(abundance))
    return levels


def compute_nrna_levels(
    transcripts: list[Transcript],
    requested_rates: tuple[float, ...],
    requested_labels: tuple[str, ...] | None = None,
) -> dict[str, tuple[float, float]]:
    """Compute nRNA abundance from total transcript abundance and nRNA rates.

    nRNA abundance is defined as:
      nRNA_abundance_total = sum(transcript abundances) * nrna_rate

    Individual transcript nRNA abundance is set per condition as:
      t.nrna_abundance = t.abundance * nrna_rate

    Returns dict[label] = (rate, total_nrna_abundance).
    """
    total_tx_abundance = float(sum(float(t.abundance) for t in transcripts))
    levels: dict[str, tuple[float, float]] = {}
    for idx, rate in enumerate(requested_rates):
        if rate < 0.0 or rate > 1.0:
            raise ValueError(f"Invalid nRNA rate {rate}; expected in [0.0, 1.0]")
        if requested_labels is not None and idx < len(requested_labels):
            label = requested_labels[idx].strip()
        else:
            label = f"r{rate:g}"
        total_abundance = total_tx_abundance * float(rate)
        levels[label] = (float(rate), float(total_abundance))
    return levels


# ── Condition directory naming ───────────────────────────────────────


def condition_dir_name(gdna_label: str, nrna_label: str, strand_specificity: float) -> str:
    return f"gdna_{gdna_label}_nrna_{nrna_label}_ss_{strand_specificity:.2f}"


# ── Main benchmark orchestration ─────────────────────────────────────


def run_region_benchmark(
    region: RegionSpec,
    args: argparse.Namespace,
    out_root: Path,
    gdna_levels: dict[str, tuple[float, float]],
    nrna_levels: dict[str, tuple[float, float]],
    strand_specificities: tuple[float, ...],
    include_htseq: bool,
    htseq_conda_env: str,
    aligner_configs: list[AlignerConfig] | None = None,
    hulkrna_configs: list[HulkrnaConfig] | None = None,
) -> list[ConditionResult]:
    """Run benchmark for one region across all conditions.

    Parameters
    ----------
    aligner_configs : list of AlignerConfig, optional
        Named aligner configurations.  When *None*, falls back to
        building default configs from ``args.aligner``.
    hulkrna_configs : list of HulkrnaConfig, optional
        Named hulkrna parameterizations.  When *None* a single
        ``HulkrnaConfig(name="default")`` is used.

    Returns a list of ConditionResult (one per condition combination).
    """
    # ── Resolve configs ─────────────────────────────────────────────

    if aligner_configs is None:
        aligner_configs = [
            AlignerConfig(name=a, type=a) for a in args.aligner
        ]
    if hulkrna_configs is None:
        hulkrna_configs = [HulkrnaConfig(name="default")]

    multi_hulk = len(hulkrna_configs) > 1
    has_hisat2 = any(ac.type == "hisat2" for ac in aligner_configs)

    reg_dir = out_root / region.label
    reg_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[REGION] %s:%d-%d", region.chrom, region.start0 + 1, region.end0)
    print(f"\n>>> REGION: {region.label} ({region.chrom}:{region.start0+1}-{region.end0})", flush=True)

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
                gtf_obj = GTFRecord(
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
    TranscriptIndex.build(region_fa, region_gtf, index_dir, write_tsv=False)
    index = TranscriptIndex.load(index_dir)

    hisat2_index_dir = reg_dir / "hisat2_index"

    # HISAT2 splice-site and exon annotation files for RNA-Seq-aware indexing
    hisat2_ss_file = reg_dir / "hisat2_splice_sites.txt"
    hisat2_exon_file = reg_dir / "hisat2_exons.txt"
    if has_hisat2:
        write_hisat2_splice_sites(transcripts, hisat2_ss_file)
        write_hisat2_exons(transcripts, hisat2_exon_file)

    # Gene mapping for gene-level aggregation
    gene_map = {t.t_id: t.g_id for t in transcripts}
    abundance_map = {t.t_id: t.abundance for t in transcripts}
    gene_abundance = defaultdict(float)
    for t in transcripts:
        gene_abundance[t.g_id] += t.abundance

    n_transcripts = len(transcripts)
    n_genes = len({t.g_id for t in transcripts})

    tx_tools = _get_transcript_tools(aligner_configs, hulkrna_configs)
    gene_tools = _get_gene_tools(aligner_configs, hulkrna_configs, include_htseq=include_htseq)

    # ── Condition sweep ─────────────────────────────────────────────

    results: list[ConditionResult] = []

    for gdna_label, (gdna_rate, gdna_abundance) in gdna_levels.items():
        for nrna_label, (nrna_rate, nrna_total_abundance) in nrna_levels.items():
            for strand_spec in strand_specificities:
                cond_name = condition_dir_name(gdna_label, nrna_label, strand_spec)
                cond_dir = reg_dir / cond_name
                cond_dir.mkdir(parents=True, exist_ok=True)

                logger.info(
                    "  [COND] %s | gDNA=%s (rate=%.3f, abundance=%.2f) | "
                    "nRNA=%s (rate=%.3f, abundance=%.2f) | strand_spec=%.2f",
                    region.label,
                    gdna_label,
                    gdna_rate,
                    gdna_abundance,
                    nrna_label,
                    nrna_rate,
                    nrna_total_abundance,
                    strand_spec,
                )
                print(
                    f"  CONDITION: gDNA={gdna_label} nRNA={nrna_label} SS={strand_spec:.2f}",
                    flush=True,
                )

                # Set per-transcript nRNA abundances for this condition.
                # Mature mRNA abundance is not reduced.
                for t in transcripts:
                    t.nrna_abundance = float(t.abundance) * float(nrna_rate)

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

                # Always write FASTQ (needed for salmon/kallisto and
                # for real aligners; oracle also uses it for truth)
                simulator = ReadSimulator(
                    sim_genome, transcripts,
                    config=sim_cfg, gdna_config=gdna_cfg,
                )
                fastq_dir = cond_dir / "reads"
                fastq_r1, fastq_r2 = simulator.write_fastq(
                    fastq_dir, args.n_fragments, prefix="sim"
                )
                truth_counts, n_gdna_actual = parse_truth_from_fastq(fastq_r1)

                # ── Per-aligner: BAM + hulkrna + htseq ──────────────

                tx_tool_counts: dict[str, tuple[dict[str, float], float]] = {}
                tool_pool_counts: dict[str, dict[str, float]] = {}
                htseq_gene_by_aligner: dict[str, tuple[dict[str, int], float]] = {}

                for ac in aligner_configs:
                    align_dir = cond_dir / f"align_{ac.name}"
                    align_dir.mkdir(parents=True, exist_ok=True)
                    bam_ns = align_dir / "reads_namesort.bam"

                    if ac.type == "oracle":
                        print(f"    Aligner: {ac.name} — generating oracle BAM...", flush=True)
                        oracle_sim = OracleBamSimulator(
                            sim_genome, transcripts,
                            config=sim_cfg, gdna_config=gdna_cfg,
                            ref_name=region.label,
                        )
                        oracle_sim.write_bam(bam_ns, args.n_fragments)
                    elif ac.type == "minimap2":
                        print(f"    Aligner: {ac.name} — running minimap2...", flush=True)
                        align_minimap2(
                            region_fa, fastq_r1, fastq_r2, bam_ns,
                            bed_path=bed_path,
                            aligner_config=ac,
                        )
                    elif ac.type == "hisat2":
                        align_hisat2(
                            region_fa, fastq_r1, fastq_r2, bam_ns,
                            threads=args.threads,
                            hisat2_exe=args.hisat2_exe,
                            hisat2_build_exe=args.hisat2_build_exe,
                            hisat2_index_dir=hisat2_index_dir,
                            splice_sites_file=hisat2_ss_file,
                            exons_file=hisat2_exon_file,
                            aligner_config=ac,
                        )
                    else:
                        raise RuntimeError(f"Unsupported aligner type: {ac.type}")

                    # Run hulkrna on this aligner's BAM — once per config
                    for hc in hulkrna_configs:
                        hn = _hulkrna_tool_name(ac, hc, multi_hulk=multi_hulk)
                        print(f"    Running hulkrna ({hn})...", end="", flush=True)
                        ann_bam = align_dir / f"annotated_{hc.name}.bam"
                        hk_counts, hk_elapsed, hk_pools = run_hulkrna_tool(
                            bam_ns, index, args=args, hulkrna_config=hc,
                            annotated_bam_path=ann_bam,
                        )
                        tx_tool_counts[hn] = (hk_counts, hk_elapsed)
                        tool_pool_counts[hn] = hk_pools
                        print(f" done ({hk_elapsed:.1f}s)", flush=True)

                        # ── Coord-sort and index annotated BAM for genome browsers ──
                        if ann_bam.exists():
                            ann_bam_sorted = align_dir / f"annotated_{hc.name}.sorted.bam"
                            coord_sort_bam(ann_bam, ann_bam_sorted)

                        # ── Per-read accuracy from annotated BAM ────
                        if ann_bam.exists():
                            acc_rows, acc_summary = compare_annotated_bam_to_truth(
                                ann_bam, gene_map,
                            )
                            acc_dir = align_dir / f"accuracy_{hc.name}"
                            acc_dir.mkdir(parents=True, exist_ok=True)
                            write_read_accuracy_csv(
                                acc_rows, acc_dir / "per_read_accuracy.csv",
                            )
                            write_read_accuracy_summary_csv(
                                acc_summary, acc_dir / "read_accuracy_summary.csv",
                            )
                            n_total = acc_summary.total
                            n_correct = acc_summary.correct_transcript
                            pct = acc_summary.pct(n_correct)
                            logger.info(
                                "    [ACCURACY] %s: %d/%d (%.1f%%) correct transcript, "
                                "%d correct gene, %d correct pool, "
                                "%d wrong pool, %d intergenic",
                                hn, n_correct, n_total, pct,
                                acc_summary.correct_gene,
                                acc_summary.correct_pool,
                                acc_summary.wrong_pool,
                                acc_summary.intergenic,
                            )

                    # htseq (gene-level, needs coord-sorted BAM)
                    if include_htseq:
                        ht_name = _htseq_tool_name(ac)
                        bam_cs = align_dir / "reads_coordsort.bam"
                        coord_sort_bam(bam_ns, bam_cs)
                        htseq_counts, htseq_elapsed = run_htseq(
                            bam_cs,
                            region_gtf,
                            cond_dir / ht_name,
                            strand_specificity=strand_spec,
                            conda_env=htseq_conda_env,
                        )
                        htseq_gene_by_aligner[ht_name] = (
                            htseq_counts, htseq_elapsed,
                        )

            # ── Run transcript-level tools (FASTQ-based, once) ──────

                print(f"    Running salmon...", end="", flush=True)
                salmon_counts, salmon_elapsed = run_salmon(
                    transcript_fa,
                    fastq_r1,
                    fastq_r2,
                    cond_dir / "salmon",
                    threads=args.threads,
                )
                print(f" done ({salmon_elapsed:.1f}s)", flush=True)
                print(f"    Running kallisto...", end="", flush=True)
                kallisto_counts, kallisto_elapsed = run_kallisto(
                    transcript_fa,
                    fastq_r1,
                    fastq_r2,
                    cond_dir / "kallisto",
                    threads=args.threads,
                    strand_specificity=strand_spec,
                )
                print(f" done ({kallisto_elapsed:.1f}s)", flush=True)

                tx_tool_counts["salmon"] = (salmon_counts, salmon_elapsed)
                tx_tool_counts["kallisto"] = (kallisto_counts, kallisto_elapsed)

                tool_pool_counts["salmon"] = {
                    "mature_rna": float(sum(salmon_counts.values())),
                    "nascent_rna": 0.0,
                    "genomic_dna": 0.0,
                    "rna": float(sum(salmon_counts.values())),
                    "dna": 0.0,
                }
                tool_pool_counts["kallisto"] = {
                    "mature_rna": float(sum(kallisto_counts.values())),
                    "nascent_rna": 0.0,
                    "genomic_dna": 0.0,
                    "rna": float(sum(kallisto_counts.values())),
                    "dna": 0.0,
                }

                pool_truth = compute_truth_pool_counts(truth_counts, n_gdna_actual)
                pool_metrics: list[dict] = []
                for tool_name, observed_pools in tool_pool_counts.items():
                    pool_metrics.extend(
                        compute_pool_metric_rows(tool_name, pool_truth, observed_pools)
                    )

            # ── Score transcript-level ──────────────────────────────

                # Use mature-only truth for per-transcript metrics.
                # nRNA truth entries (nrna_ENST...) are evaluated
                # exclusively via pool-level metrics.
                # Include ALL annotated transcripts so that false
                # positives on truth=0 transcripts are captured.
                mature_truth_counts = {
                    tid: 0 for tid in gene_map
                }
                mature_truth_counts.update({
                    tid: count
                    for tid, count in truth_counts.items()
                    if not str(tid).startswith("nrna_")
                })

                transcript_metrics = {}
                for tool, (counts, elapsed) in tx_tool_counts.items():
                    transcript_metrics[tool] = asdict(
                        score_tool(tool, mature_truth_counts, counts, elapsed)
                    )

            # ── Score gene-level ────────────────────────────────────

                truth_gene = aggregate_to_genes(mature_truth_counts, gene_map)
                # Ensure all annotated genes appear in truth_gene
                for gid in gene_abundance:
                    truth_gene.setdefault(gid, 0)
                gene_metrics = {}
                for tool, (counts, elapsed) in tx_tool_counts.items():
                    gene_counts = aggregate_to_genes(counts, gene_map)
                    gene_metrics[tool] = asdict(
                        score_tool(tool, truth_gene, gene_counts, elapsed)
                    )
                for htseq_tool, (htseq_counts, htseq_elapsed) in htseq_gene_by_aligner.items():
                    gene_metrics[htseq_tool] = asdict(
                        score_tool(htseq_tool, truth_gene, htseq_counts, htseq_elapsed)
                    )

            # ── Write per-transcript CSV ────────────────────────────

                per_tx_csv = cond_dir / "per_transcript_counts.csv"
                tids = sorted(mature_truth_counts)
                with open(per_tx_csv, "w", newline="") as f:
                    writer = csv.writer(f)
                    header = ["transcript_id", "gene_id", "truth"]
                    header.extend(tx_tools)
                    header.extend(["abundance", "nrna_abundance"])
                    writer.writerow(header)
                    for tid in tids:
                        row: list = [
                            tid,
                            gene_map.get(tid, ""),
                            int(truth_counts.get(tid, 0)),
                        ]
                        for tool in tx_tools:
                            row.append(float(tx_tool_counts[tool][0].get(tid, 0.0)))
                        row.extend([
                            float(abundance_map.get(tid, 0.0)),
                            float(abundance_map.get(tid, 0.0) * nrna_rate),
                        ])
                        writer.writerow(row)

            # ── Write per-gene CSV ──────────────────────────────────

                per_gene_csv = cond_dir / "per_gene_counts.csv"
                gene_ids = sorted(truth_gene)
                # Aggregate each tool's counts to gene-level
                gene_level_counts: dict[str, dict[str, float]] = {}
                for tool in tx_tools:
                    gene_level_counts[tool] = aggregate_to_genes(
                        tx_tool_counts[tool][0], gene_map,
                    )
                for htseq_tool, (htseq_counts, _) in htseq_gene_by_aligner.items():
                    gene_level_counts[htseq_tool] = {
                        gid: float(c) for gid, c in htseq_counts.items()
                    }
                with open(per_gene_csv, "w", newline="") as f:
                    writer = csv.writer(f)
                    header = ["gene_id", "truth"]
                    header.extend(gene_tools)
                    header.extend(["abundance", "nrna_abundance"])
                    writer.writerow(header)
                    for gid in gene_ids:
                        row = [
                            gid,
                            int(truth_gene.get(gid, 0)),
                        ]
                        for tool in gene_tools:
                            row.append(float(gene_level_counts.get(tool, {}).get(gid, 0.0)))
                        row.extend([
                            float(gene_abundance.get(gid, 0.0)),
                            float(gene_abundance.get(gid, 0.0) * nrna_rate),
                        ])
                        writer.writerow(row)

                per_pool_csv = cond_dir / "per_pool_counts.csv"
                with open(per_pool_csv, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([
                        "tool",
                        "pool",
                        "truth",
                        "observed",
                        "signed_error",
                        "abs_error",
                    ])
                    for pm_row in pool_metrics:
                        writer.writerow([
                            pm_row["tool"],
                            pm_row["pool"],
                            pm_row["truth"],
                            pm_row["observed"],
                            pm_row["signed_error"],
                            pm_row["abs_error"],
                        ])

                aligner_str = ",".join(ac.name for ac in aligner_configs)
                results.append(
                    ConditionResult(
                        region=region.label,
                        n_transcripts=n_transcripts,
                        n_genes=n_genes,
                        n_fragments=args.n_fragments,
                        n_gdna_actual=n_gdna_actual,
                        gdna_label=gdna_label,
                        gdna_rate=round(gdna_rate, 6),
                        gdna_abundance=round(gdna_abundance, 4),
                        nrna_label=nrna_label,
                        nrna_rate=round(nrna_rate, 6),
                        nrna_abundance=round(nrna_total_abundance, 4),
                        strand_specificity=strand_spec,
                        aligner=aligner_str,
                        transcript_metrics=transcript_metrics,
                        gene_metrics=gene_metrics,
                        pool_metrics=pool_metrics,
                    )
                )

                hulkrna_parts = " | ".join(
                    f"{_hulkrna_tool_name(ac, hc, multi_hulk=multi_hulk)} MAE="
                    f"{transcript_metrics[_hulkrna_tool_name(ac, hc, multi_hulk=multi_hulk)]['mean_abs_error']:.2f}"
                    for ac in aligner_configs
                    for hc in hulkrna_configs
                )
                htseq_parts = ""
                if include_htseq:
                    htseq_parts = " | " + " | ".join(
                        f"{ht} gene-MAE={gene_metrics[ht]['mean_abs_error']:.2f}"
                        for ht in htseq_gene_by_aligner
                    )
                logger.info(
                    "    aligners=%s | %s | "
                    "salmon MAE=%.2f | kallisto MAE=%.2f%s",
                    aligner_str,
                    hulkrna_parts,
                    transcript_metrics["salmon"]["mean_abs_error"],
                    transcript_metrics["kallisto"]["mean_abs_error"],
                    htseq_parts,
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
            "aligner",
            "n_transcripts",
            "n_genes",
            "n_fragments",
            "n_gdna_actual",
            "gdna_label",
            "gdna_rate",
            "gdna_abundance",
            "nrna_label",
            "nrna_rate",
            "nrna_abundance",
            "strand_specificity",
            "level",
            "tool",
            "elapsed_sec",
            "total_abs_error",
            "mean_abs_error",
            "rmse",
            "pearson",
            "spearman",
            "pool",
            "truth",
            "observed",
            "signed_error",
        ])
        for r in results:
            base = [
                r.region,
                r.aligner,
                r.n_transcripts,
                r.n_genes,
                r.n_fragments,
                r.n_gdna_actual,
                r.gdna_label,
                r.gdna_rate,
                r.gdna_abundance,
                r.nrna_label,
                r.nrna_rate,
                r.nrna_abundance,
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
                        "",
                        "",
                        "",
                        "",
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
                        "",
                        "",
                        "",
                        "",
                    ]
                )
            for pm in r.pool_metrics:
                writer.writerow(
                    base
                    + [
                        "pool",
                        pm["tool"],
                        "",
                        pm["abs_error"],
                        pm["abs_error"],
                        "",
                        "",
                        "",
                        pm["pool"],
                        pm["truth"],
                        pm["observed"],
                        pm["signed_error"],
                    ]
                )

    # ── Markdown ────────────────────────────────────────────────────

    # Derive tool lists dynamically from the first result
    if results:
        tx_tools = list(results[0].transcript_metrics.keys())
        gene_tools = list(results[0].gene_metrics.keys())
    else:
        tx_tools = list(TRANSCRIPT_TOOLS)
        gene_tools = list(GENE_TOOLS) if include_htseq else list(TRANSCRIPT_TOOLS)

    lines = [
        "# Regional Benchmark Summary",
        "",
        "## Transcript-Level Metrics",
        "",
        "| Region | Aligner | gDNA | SS | Tool | MAE | RMSE | Pearson | Spearman | Time (s) |",
        "| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in results:
        for tool in tx_tools:
            tm = r.transcript_metrics.get(tool)
            if tm is None:
                continue
            lines.append(
                f"| {r.region} | {r.aligner} | {r.gdna_label} | {r.strand_specificity:.2f} | "
                f"{tm['tool']} | {tm['mean_abs_error']:.2f} | {tm['rmse']:.2f} | "
                f"{tm['pearson']:.4f} | {tm['spearman']:.4f} | {tm['elapsed_sec']:.2f} |"
            )

    lines.extend([
        "",
        "## Gene-Level Metrics",
        "",
        "| Region | Aligner | gDNA | SS | Tool | MAE | RMSE | Pearson | Spearman | Time (s) |",
        "| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for r in results:
        for tool in gene_tools:
            tm = r.gene_metrics.get(tool)
            if tm is None:
                continue
            lines.append(
                f"| {r.region} | {r.aligner} | {r.gdna_label} | {r.strand_specificity:.2f} | "
                f"{tm['tool']} | {tm['mean_abs_error']:.2f} | {tm['rmse']:.2f} | "
                f"{tm['pearson']:.4f} | {tm['spearman']:.4f} | {tm['elapsed_sec']:.2f} |"
            )

    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_diagnostics(
    results: list[ConditionResult], out_root: Path, include_htseq: bool
) -> None:
    """Write aggregate diagnostic report."""
    # Derive tool lists dynamically from results
    if results:
        tx_tools = list(results[0].transcript_metrics.keys())
    else:
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
        for gdna_label in sorted(df["gdna_label"].unique()):
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
        cond_name = condition_dir_name(r.gdna_label, r.nrna_label, r.strand_specificity)
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

    pool_rows = []
    for r in results:
        for pm in r.pool_metrics:
            pool_rows.append({
                "region": r.region,
                "gdna_label": r.gdna_label,
                "nrna_label": r.nrna_label,
                "strand_specificity": r.strand_specificity,
                "tool": pm["tool"],
                "pool": pm["pool"],
                "abs_error": pm["abs_error"],
                "signed_error": pm["signed_error"],
            })
    pool_df = pd.DataFrame(pool_rows)
    if len(pool_df):
        pool_tool_header = "| Pool | " + " | ".join(tx_tools) + " |"
        pool_tool_sep = "| --- | " + " | ".join(["---:"] * len(tx_tools)) + " |"
        lines.extend([
            "",
            "## Pool-level absolute error (mature_rna, nascent_rna, genomic_dna, rna, dna)",
            "",
            pool_tool_header,
            pool_tool_sep,
        ])
        for pool in ["mature_rna", "nascent_rna", "genomic_dna", "rna", "dna"]:
            sub = pool_df[pool_df["pool"] == pool]
            vals = []
            for tool in tx_tools:
                tsub = sub[sub["tool"] == tool]
                vals.append(f"{tsub['abs_error'].mean():.3f}" if len(tsub) else "—")
            lines.append(f"| {pool} | " + " | ".join(vals) + " |")

    out = out_root / "diagnostics.md"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")


# ── CLI ──────────────────────────────────────────────────────────────


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Regional benchmark: hulkrna vs salmon vs kallisto vs htseq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--genome", type=Path, required=False, default=None, help="Genome FASTA path")
    p.add_argument(
        "--gtf", type=Path, required=False, default=None, help="Gene annotation GTF/GTF.GZ path"
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
    p.add_argument("--outdir", type=Path, required=False, default=None, help="Output directory")
    p.add_argument(
        "--config",
        type=Path,
        default=None,
        help="YAML configuration file (CLI flags override config values)",
    )

    p.add_argument("--n-fragments", type=int, default=20000)
    p.add_argument("--threads", type=int, default=1)

    p.add_argument(
        "--aligner",
        nargs="+",
        choices=list(ALIGNER_CHOICES),
        default=["minimap2"],
        help="Genome aligner(s) for hulkrna/htseq inputs; specify one or more "
             "(e.g. --aligner minimap2 oracle). oracle: bypass alignment with perfect BAM",
    )
    p.add_argument(
        "--hisat2-exe",
        type=str,
        default="hisat2",
        help="HISAT2 executable name/path (used with --aligner hisat2)",
    )
    p.add_argument(
        "--hisat2-build-exe",
        type=str,
        default="hisat2-build",
        help="hisat2-build executable name/path (used with --aligner hisat2)",
    )

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
        "--gdna-rates",
        type=str,
        default=",".join(str(x) for x in DEFAULT_GDNA_RATES),
        help="Comma-separated total gDNA rates in [0,1]; abundance is sum(tx_abundance)*rate",
    )
    p.add_argument(
        "--gdna-rate-labels",
        type=str,
        default="",
        help="Optional comma-separated labels matching --gdna-rates (same length)",
    )
    p.add_argument(
        "--nrna-rates",
        type=str,
        default=",".join(str(x) for x in DEFAULT_NRNA_RATES),
        help="Comma-separated total nRNA rates in [0,1]; per transcript nRNA = abundance*rate",
    )
    p.add_argument(
        "--nrna-rate-labels",
        type=str,
        default="",
        help="Optional comma-separated labels matching --nrna-rates (same length)",
    )
    p.add_argument(
        "--gdna-levels",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )
    p.add_argument("--gdna-frag-mean", type=float, default=350.0)
    p.add_argument("--gdna-frag-std", type=float, default=100.0)
    p.add_argument("--gdna-frag-min", type=int, default=100)
    p.add_argument("--gdna-frag-max", type=int, default=1000)
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
        "--no-htseq",
        action="store_true",
        help="Disable htseq-count (enabled by default)",
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

    # Overhang alpha
    p.add_argument(
        "--overhang-alpha",
        dest="overhang_alpha",
        type=float,
        default=None,
        help="Overhang alpha for hulkrna (default: use pipeline default 0.01). "
             "0 = hard binary, 1 = off.",
    )

    return p


def ensure_tools(
    *,
    aligner_configs: list[AlignerConfig],
    hisat2_exe: str,
    hisat2_build_exe: str,
    include_htseq: bool = True,
    htseq_conda_env: str = "htseq",
) -> None:
    required = ["samtools", "salmon", "kallisto"]
    for ac in aligner_configs:
        if ac.type == "minimap2":
            required.append("minimap2")
        elif ac.type == "hisat2":
            required.extend([hisat2_exe, hisat2_build_exe])
    # Deduplicate while preserving order
    required = list(dict.fromkeys(required))
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


def _flatten_config_dict(d: dict, prefix: str = "") -> dict:
    """Flatten nested YAML dicts for scalar option mapping.

    Structured sections (``hulkrna_configs``, ``aligners``) are
    preserved as-is and NOT flattened — they are handled separately
    by :func:`apply_yaml_config`.
    """
    # Sections that should be kept as structured dicts
    _STRUCTURED_SECTIONS = {"hulkrna_configs", "aligners"}

    out: dict = {}
    for key, value in d.items():
        joined = f"{prefix}_{key}" if prefix else str(key)
        if not prefix and key in _STRUCTURED_SECTIONS:
            out[key] = value
            continue
        if isinstance(value, dict):
            out.update(_flatten_config_dict(value, joined))
        else:
            out[joined] = value
    return out


def _cli_dest_overrides(argv: list[str]) -> set[str]:
    dests: set[str] = set()
    for tok in argv:
        if not tok.startswith("--"):
            continue
        flag = tok[2:].split("=", 1)[0]
        dests.add(flag.replace("-", "_"))
    return dests


def apply_yaml_config(args: argparse.Namespace, argv: list[str]) -> argparse.Namespace:
    """Apply YAML configuration to ``args``, respecting CLI overrides.

    Handles three YAML‑only sections that have no CLI equivalent:

    - **region** entries can be plain strings (``chr:start-end``) or
      single-key dicts (``{Name: "chr:start-end"}``) for named regions.
    - **hulkrna_configs** – mapping of ``name → {param: value, …}``.
    - **aligners** – mapping of ``name → {type: …, param: value, …}``.
    """
    if args.config is None:
        return args
    if yaml is None:
        raise RuntimeError(
            "YAML config requested but PyYAML is not installed. Install with `pip install pyyaml`."
        )
    with open(args.config) as fh:
        cfg_raw = yaml.safe_load(fh) or {}
    if not isinstance(cfg_raw, dict):
        raise ValueError("YAML config root must be a mapping/object")

    cfg = _flatten_config_dict(cfg_raw)
    cli_overrides = _cli_dest_overrides(argv)

    # ── Structured sections (not in argparse) ───────────────────────

    # hulkrna_configs
    if "hulkrna_configs" in cfg and "hulkrna_configs" not in cli_overrides:
        raw_hc = cfg.pop("hulkrna_configs")
        if isinstance(raw_hc, dict):
            args.hulkrna_configs = [
                HulkrnaConfig(name=str(name), params=dict(params) if isinstance(params, dict) else {})
                for name, params in raw_hc.items()
            ]
        else:
            raise ValueError("hulkrna_configs must be a mapping of name → {params}")

    # aligners
    if "aligners" in cfg and "aligner" not in cli_overrides:
        raw_al = cfg.pop("aligners")
        if isinstance(raw_al, dict):
            aligner_cfgs: list[AlignerConfig] = []
            for aname, aval in raw_al.items():
                aname = str(aname)
                if aval is None:
                    aval = {}
                if not isinstance(aval, dict):
                    raise ValueError(
                        f"aligners.{aname} must be a mapping (or empty), got {type(aval).__name__}"
                    )
                atype = str(aval.pop("type", aname))
                if atype not in ALIGNER_CHOICES:
                    raise ValueError(
                        f"aligners.{aname}: unsupported type '{atype}'; "
                        f"choose from {ALIGNER_CHOICES}"
                    )
                aligner_cfgs.append(AlignerConfig(name=aname, type=atype, params=dict(aval)))
            args.aligner_configs = aligner_cfgs
        else:
            raise ValueError("aligners must be a mapping of name → {type, params…}")

    # ── aliases for nested config sections ──────────────────────────

    alias_map = {
        "simulation_n_fragments": "n_fragments",
        "simulation_sim_seed": "sim_seed",
        "simulation_pipeline_seed": "pipeline_seed",
        "simulation_frag_mean": "frag_mean",
        "simulation_frag_std": "frag_std",
        "simulation_frag_min": "frag_min",
        "simulation_frag_max": "frag_max",
        "simulation_read_length": "read_length",
        "simulation_error_rate": "error_rate",
        "gdna_rates": "gdna_rates",
        "gdna_rate_labels": "gdna_rate_labels",
        "nrna_rates": "nrna_rates",
        "nrna_rate_labels": "nrna_rate_labels",
        "gdna_frag_mean": "gdna_frag_mean",
        "gdna_frag_std": "gdna_frag_std",
        "gdna_frag_min": "gdna_frag_min",
        "gdna_frag_max": "gdna_frag_max",
        "abundance_mode": "abundance_mode",
        "abundance_seed": "abundance_seed",
        "abundance_min": "abundance_min",
        "abundance_max": "abundance_max",
        "htseq_enabled": "no_htseq",
        "htseq_conda_env": "htseq_conda_env",
        "runtime_keep_going": "keep_going",
    }

    for raw_key, value in cfg.items():
        # Skip already-handled structured sections
        if raw_key in {"hulkrna_configs", "aligners"}:
            continue
        dest = raw_key.replace("-", "_")
        dest = alias_map.get(dest, dest)
        if not hasattr(args, dest):
            continue
        if dest in cli_overrides:
            continue
        if dest in {"genome", "gtf", "outdir", "region_file", "abundance_file", "config"} and isinstance(value, str):
            setattr(args, dest, Path(value))
            continue
        # Region: support list of strings, dicts, or a bare string
        if dest == "region" and isinstance(value, str):
            setattr(args, dest, [value])
            continue
        if dest == "region" and isinstance(value, list):
            # Preserve dicts for named regions; stringify bare values
            entries: list = []
            for item in value:
                if isinstance(item, dict):
                    entries.append(item)
                else:
                    entries.append(str(item))
            setattr(args, dest, entries)
            continue
        if dest == "aligner" and isinstance(value, str):
            setattr(args, dest, [value])
            continue
        if dest == "aligner" and isinstance(value, list):
            setattr(args, dest, [str(x) for x in value])
            continue
        if dest in {"gdna_rates", "gdna_rate_labels", "nrna_rates", "nrna_rate_labels", "strand_specificities"} and isinstance(value, list):
            setattr(args, dest, ",".join(str(x) for x in value))
            continue
        if dest == "no_htseq" and raw_key.replace("-", "_") == "htseq_enabled":
            setattr(args, dest, not bool(value))
            continue
        setattr(args, dest, value)
    return args


def main() -> int:
    parser = build_arg_parser()
    parser.add_argument(
        "--include-htseq",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    args = apply_yaml_config(args, argv)

    # Force unbuffered logging to stderr so progress is visible
    # even when running under 'conda run' (which may buffer stdout).
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        handlers=[handler],
    )

    if args.genome is None:
        raise ValueError("--genome is required (or set genome in --config)")
    if args.gtf is None:
        raise ValueError("--gtf is required (or set gtf in --config)")
    if args.outdir is None:
        raise ValueError("--outdir is required (or set outdir in --config)")

    include_htseq = not args.no_htseq
    if getattr(args, "include_htseq", False):
        include_htseq = True
    htseq_conda_env = args.htseq_conda_env

    # ── Build AlignerConfig list ────────────────────────────────────
    # Priority:  YAML ``aligners:`` section  >  CLI/YAML ``aligner:``

    aligner_configs: list[AlignerConfig] = getattr(args, "aligner_configs", None) or []
    if not aligner_configs:
        # Normalise aligner to a list (YAML may set a bare string)
        if isinstance(args.aligner, str):
            args.aligner = [args.aligner]
        aligner_configs = [
            AlignerConfig(name=a, type=a) for a in args.aligner
        ]

    # ── Build HulkrnaConfig list ────────────────────────────────────

    hulkrna_configs: list[HulkrnaConfig] = getattr(args, "hulkrna_configs", None) or []
    if not hulkrna_configs:
        hulkrna_configs = [HulkrnaConfig(name="default")]

    # Keep args.aligner in sync for any code that still reads it
    args.aligner = [ac.name for ac in aligner_configs]

    ensure_tools(
        aligner_configs=aligner_configs,
        hisat2_exe=args.hisat2_exe,
        hisat2_build_exe=args.hisat2_build_exe,
        include_htseq=include_htseq,
        htseq_conda_env=htseq_conda_env,
    )

    if args.abundance_mode == "file" and args.abundance_file is None:
        raise ValueError("--abundance-file is required when --abundance-mode=file")

    if args.gdna_levels:
        legacy_rate_map = {
            "none": 0.0,
            "low": 0.1,
            "moderate": 0.25,
            "high": 0.5,
        }
        tokens = [t.strip() for t in str(args.gdna_levels).split(",") if t.strip()]
        if tokens:
            mapped: list[float] = []
            for token in tokens:
                if token not in legacy_rate_map:
                    raise ValueError(
                        f"Unsupported --gdna-levels token '{token}'. Use --gdna-rates instead."
                    )
                mapped.append(legacy_rate_map[token])
            args.gdna_rates = ",".join(str(x) for x in mapped)
            if not args.gdna_rate_labels:
                args.gdna_rate_labels = ",".join(tokens)
            logger.warning("--gdna-levels is deprecated; use --gdna-rates")

    regions = parse_regions(args)

    # Parse gDNA rates and strand specificities from CLI/config
    gdna_rates = tuple(
        float(x.strip())
        for x in str(args.gdna_rates).split(",")
        if x.strip()
    )
    if not gdna_rates:
        raise ValueError("At least one gDNA rate is required")
    gdna_rate_labels = tuple(
        x.strip()
        for x in str(args.gdna_rate_labels).split(",")
        if x.strip()
    )
    if gdna_rate_labels and len(gdna_rate_labels) != len(gdna_rates):
        raise ValueError("--gdna-rate-labels must have the same number of entries as --gdna-rates")

    nrna_rates = tuple(
        float(x.strip())
        for x in str(args.nrna_rates).split(",")
        if x.strip()
    )
    if not nrna_rates:
        raise ValueError("At least one nRNA rate is required")
    nrna_rate_labels = tuple(
        x.strip()
        for x in str(args.nrna_rate_labels).split(",")
        if x.strip()
    )
    if nrna_rate_labels and len(nrna_rate_labels) != len(nrna_rates):
        raise ValueError("--nrna-rate-labels must have the same number of entries as --nrna-rates")

    strand_specificities = tuple(
        float(s.strip())
        for s in str(args.strand_specificities).split(",")
        if s.strip()
    )

    n_aligner = len(aligner_configs)
    n_hulk = len(hulkrna_configs)
    n_conditions = len(gdna_rates) * len(nrna_rates) * len(strand_specificities) * len(regions)
    logger.info(
        "Benchmark grid: %d gDNA rates × %d nRNA rates × %d strand specs "
        "× %d regions × %d aligners × %d hulkrna configs = %d conditions",
        len(gdna_rates),
        len(nrna_rates),
        len(strand_specificities),
        len(regions),
        n_aligner,
        n_hulk,
        n_conditions,
    )
    print(
        f"\n{'='*60}\n"
        f"BENCHMARK: {n_conditions} conditions, "
        f"{len(regions)} regions, {n_aligner} aligners\n"
        f"{'='*60}",
        flush=True,
    )
    logger.info(
        "Aligners: %s",
        ", ".join(f"{ac.name} ({ac.type})" for ac in aligner_configs),
    )
    if n_hulk > 1:
        logger.info(
            "HulkRNA configs: %s",
            ", ".join(
                f"{hc.name}" + (f" ({hc.params})" if hc.params else "")
                for hc in hulkrna_configs
            ),
        )

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_results: list[ConditionResult] = []

    for region in regions:
        try:
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
            gdna_levels = compute_gdna_levels(
                transcripts_pre,
                gdna_rates,
                gdna_rate_labels if gdna_rate_labels else None,
            )
            nrna_levels = compute_nrna_levels(
                transcripts_pre,
                nrna_rates,
                nrna_rate_labels if nrna_rate_labels else None,
            )
            logger.info(
                "  gDNA levels for %s: %s",
                region.label,
                {k: f"rate={v[0]:.3f}, abundance={v[1]:.2f}" for k, v in gdna_levels.items()},
            )
            logger.info(
                "  nRNA levels for %s: %s",
                region.label,
                {k: f"rate={v[0]:.3f}, abundance={v[1]:.2f}" for k, v in nrna_levels.items()},
            )

            region_results = run_region_benchmark(
                region,
                args,
                args.outdir,
                gdna_levels=gdna_levels,
                nrna_levels=nrna_levels,
                strand_specificities=strand_specificities,
                include_htseq=include_htseq,
                htseq_conda_env=htseq_conda_env,
                aligner_configs=aligner_configs,
                hulkrna_configs=hulkrna_configs,
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
    print(f"\n{'='*60}\nBENCHMARK COMPLETE — results in {args.outdir}\n{'='*60}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
