#!/usr/bin/env python3
"""Whole-genome RNA-seq read simulator.

Generates paired-end RNA-seq reads across the full transcriptome with
configurable mRNA/nRNA abundances, gDNA contamination, and strand
specificity.  Outputs FASTQ files and optional oracle BAM files.

Abundance model
---------------
Two modes:

**random** — For each transcript, sample whether it is expressed
(Bernoulli with probability ``frac_expressed``).  For expressed
transcripts, draw total RNA abundance from a **log-uniform**
distribution:

    log_total ~ Uniform(log(min), log(max))
    total = exp(log_total)

Then draw a per-transcript nascent RNA fraction:

    nrna_frac ~ Uniform(nrna_frac_min, nrna_frac_max)
    mRNA = total × (1 − nrna_frac)
    nRNA = total × nrna_frac

Single-exon transcripts always get nRNA = 0.

**file** — Load from a TSV with columns ``transcript_id``,
``mrna_abundance``, ``nrna_abundance``.

Fragment allocation
-------------------
The number of RNA fragments (``n_rna_fragments``) is fixed across
all conditions.  gDNA fragments are *added on top*:

    n_gdna = round(gdna_rate × n_rna_fragments)
    n_total = n_rna_fragments + n_gdna

Condition grid
--------------
Sweeps: ``nrna_fracs × gdna_rates × strand_specificities``.

When the abundance file provides explicit nRNA data, the nRNA sweep is
skipped (single condition using the file's nRNA values).

Usage
-----
::

    python scripts/sim.py --config scripts/sim_example.yaml

Output
------
::

    <outdir>/
        manifest.json
        truth_abundances.tsv
        <condition>/               e.g. gdna_none_ss_1.00
            sim_R1.fq.gz, sim_R2.fq.gz
            sim_oracle.bam         (if oracle_bam enabled)
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
import copy
from collections import defaultdict
from typing import Generator

import numpy as np
import pysam

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]

from rigel.transcript import Transcript
from rigel.types import Strand

logger = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────

_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_DNA_BASES = np.array(list("ACGT"), dtype="<U1")

# SAM flag bits
_FLAG_PAIRED = 0x1
_FLAG_PROPER_PAIR = 0x2
_FLAG_REVERSE = 0x10
_FLAG_MATE_REVERSE = 0x20
_FLAG_READ1 = 0x40
_FLAG_READ2 = 0x80
_BASE_R1_FLAG = _FLAG_PAIRED | _FLAG_PROPER_PAIR | _FLAG_READ1
_BASE_R2_FLAG = _FLAG_PAIRED | _FLAG_PROPER_PAIR | _FLAG_READ2


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_DNA_COMPLEMENT)[::-1]


# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════


@dataclass
class SimulationParams:
    """Read simulation parameters."""

    n_rna_fragments: int = 1_000_000
    sim_seed: int = 42
    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0


@dataclass
class AbundanceConfig:
    """Abundance assignment configuration."""

    mode: str = "random"  # "random" or "file"
    seed: int = 42
    min: float = 0.1
    max: float = 10000.0
    frac_expressed: float = 0.6
    file: str | None = None


@dataclass
class GDNASimConfig:
    """gDNA contamination configuration."""

    rates: list[float] = field(default_factory=lambda: [0.0])
    rate_labels: list[str] | None = None
    frag_mean: float = 350.0
    frag_std: float = 100.0
    frag_min: int = 100
    frag_max: int = 1000


@dataclass
class NRNAConfig:
    """Nascent RNA spike-in sweep configuration.

    Each entry in ``fracs`` is a (min, max) range for the per-transcript
    nRNA fraction.  A separate simulation condition is produced for each
    entry (crossed with gDNA rates and strand specificities).

    When abundances come from a file that already contains explicit nRNA
    data, the ``fracs`` sweep is ignored — the file's nRNA values are
    used as-is in a single nRNA condition.
    """

    fracs: list[tuple[float, float]] = field(
        default_factory=lambda: [(0.0, 0.5)]
    )
    frac_labels: list[str] | None = None


@dataclass
class SimConfig:
    """Top-level simulation configuration."""

    genome: str = ""
    gtf: str = ""
    outdir: str = "sim_output"
    transcript_filter: str = "all"  # "all", "basic", "mane", "ccds"

    simulation: SimulationParams = field(default_factory=SimulationParams)
    abundance: AbundanceConfig = field(default_factory=AbundanceConfig)
    gdna: GDNASimConfig = field(default_factory=GDNASimConfig)
    nrna: NRNAConfig = field(default_factory=NRNAConfig)
    strand_specificities: list[float] = field(
        default_factory=lambda: [1.0]
    )

    oracle_bam: bool = True
    verbose: bool = True


def parse_yaml_config(path: str | Path) -> SimConfig:
    """Parse a YAML config file into a SimConfig."""
    if yaml is None:
        raise ImportError("PyYAML is required: pip install pyyaml")

    with open(path) as f:
        raw = yaml.safe_load(f)

    cfg = SimConfig()
    cfg.genome = raw.get("genome", "")
    cfg.gtf = raw.get("gtf", "")
    cfg.outdir = raw.get("outdir", "sim_output")
    cfg.transcript_filter = raw.get("transcript_filter", "all")

    # Simulation params
    sim_raw = raw.get("simulation", {})
    sim = cfg.simulation
    sim.n_rna_fragments = int(sim_raw.get("n_rna_fragments", 1_000_000))
    sim.sim_seed = int(sim_raw.get("sim_seed", 42))
    sim.frag_mean = float(sim_raw.get("frag_mean", 250.0))
    sim.frag_std = float(sim_raw.get("frag_std", 50.0))
    sim.frag_min = int(sim_raw.get("frag_min", 50))
    sim.frag_max = int(sim_raw.get("frag_max", 1000))
    sim.read_length = int(sim_raw.get("read_length", 150))
    sim.error_rate = float(sim_raw.get("error_rate", 0.0))

    # Abundance
    ab_raw = raw.get("abundance", {})
    ab = cfg.abundance
    ab.mode = ab_raw.get("mode", "random")
    ab.seed = int(ab_raw.get("seed", 42))
    ab.min = float(ab_raw.get("min", 0.1))
    ab.max = float(ab_raw.get("max", 10000.0))
    ab.frac_expressed = float(ab_raw.get("frac_expressed", 0.6))
    ab.file = ab_raw.get("file", None)

    # nRNA spike-in sweep — top-level "nrna:" section is canonical
    nrna_raw = raw.get("nrna", {})
    nrna = cfg.nrna
    raw_fracs = nrna_raw.get("fracs", None)
    if raw_fracs is not None:
        nrna.fracs = [(float(p[0]), float(p[1])) for p in raw_fracs]
    nrna.frac_labels = nrna_raw.get("frac_labels", None)

    # gDNA
    gd_raw = raw.get("gdna", {})
    gd = cfg.gdna
    gd.rates = [float(r) for r in gd_raw.get("rates", [0.0])]
    gd.rate_labels = gd_raw.get("rate_labels", None)
    gd.frag_mean = float(gd_raw.get("frag_mean", 350.0))
    gd.frag_std = float(gd_raw.get("frag_std", 100.0))
    gd.frag_min = int(gd_raw.get("frag_min", 100))
    gd.frag_max = int(gd_raw.get("frag_max", 1000))

    # Strand specificities
    cfg.strand_specificities = [
        float(s) for s in raw.get("strand_specificities", [1.0])
    ]

    # Misc
    cfg.oracle_bam = bool(raw.get("oracle_bam", True))
    cfg.verbose = bool(raw.get("verbose", True))

    return cfg


# ═══════════════════════════════════════════════════════════════════
# Transcript loading and filtering
# ═══════════════════════════════════════════════════════════════════


def load_transcripts(
    gtf_path: str | Path,
    *,
    transcript_filter: str = "all",
) -> list[Transcript]:
    """Load transcripts from a GTF file with optional filtering.

    Parameters
    ----------
    gtf_path : path
        GTF annotation file (may be gzipped).
    transcript_filter : str
        One of ``"all"``, ``"basic"``, ``"mane"``, ``"ccds"``.

    Returns
    -------
    list[Transcript]
        With ``t_index`` assigned sequentially.
    """
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

    n_genes = len({t.g_id for t in transcripts})
    logger.info("Final: %d transcripts from %d genes", len(transcripts), n_genes)
    return transcripts


# ═══════════════════════════════════════════════════════════════════
# Abundance assignment
# ═══════════════════════════════════════════════════════════════════


def assign_random_abundances(
    transcripts: list[Transcript],
    config: AbundanceConfig,
) -> None:
    """Assign random total-RNA abundances using log-uniform sampling.

    All abundance is assigned to ``t.abundance`` (mRNA).  ``nrna_abundance``
    is left at zero.  The nRNA spike-in is applied per-condition in
    ``run_simulation`` via ``_spike_in_nrna``.

    1. Bernoulli(frac_expressed) → expressed flag.
    2. For expressed: total_RNA ~ LogUniform(min, max).
    3. Assign total RNA to ``t.abundance``, nRNA = 0.
    """
    if config.min <= 0:
        raise ValueError(f"abundance min must be > 0 for log-uniform, got {config.min}")

    rng = np.random.default_rng(config.seed)
    n = len(transcripts)

    # Step 1: expressed?
    expressed = rng.random(n) < config.frac_expressed

    # Step 2: log-uniform total RNA for expressed transcripts
    total_rna = np.zeros(n)
    n_expr = int(expressed.sum())
    log_min = np.log(config.min)
    log_max = np.log(config.max)
    total_rna[expressed] = np.exp(rng.uniform(log_min, log_max, size=n_expr))

    # Step 3: assign total RNA as mRNA, nRNA = 0
    for i, t in enumerate(transcripts):
        t.abundance = float(total_rna[i])
        t.nrna_abundance = 0.0

    total_mrna = sum(t.abundance for t in transcripts)
    logger.info(
        "Log-uniform abundances: %d/%d expressed, total mRNA=%.1f",
        n_expr, n, total_mrna,
    )


# ═══════════════════════════════════════════════════════════════════
# Abundance file loaders (salmon, kallisto, sim TSV)
# ═══════════════════════════════════════════════════════════════════


def _detect_abundance_format(
    path: str | Path,
) -> tuple[str, "pd.DataFrame"]:
    """Auto-detect file format and return (format_name, dataframe).

    Accepts plain-text or gzip-compressed (``.gz``) files.

    Supported formats
    -----------------
    salmon   — ``quant.sf`` with columns ``Name``, ``TPM``.
    kallisto — ``abundance.tsv`` with columns ``target_id``, ``tpm``.
    sim      — sim.py TSV with ``transcript_id``, ``mrna_abundance``
               (and optionally ``nrna_abundance``).
    """
    import pandas as pd

    # Detect gzip by magic bytes (robust even without .gz extension)
    compression: str = "infer"
    try:
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                compression = "gzip"
    except OSError:
        pass

    df = pd.read_csv(path, sep="\t", compression=compression)
    cols = set(df.columns)

    # salmon quant.sf
    if {"Name", "TPM"}.issubset(cols):
        return "salmon", df
    # kallisto abundance.tsv
    if {"target_id", "tpm"}.issubset(cols):
        return "kallisto", df
    # sim TSV (transcript_id + mrna_abundance required)
    if {"transcript_id", "mrna_abundance"}.issubset(cols):
        return "sim", df
    # Try case-insensitive fallback for kallisto (some versions capitalize)
    low_cols = {c.lower() for c in df.columns}
    if {"target_id", "tpm"}.issubset(low_cols):
        df.columns = [c.lower() for c in df.columns]
        return "kallisto", df

    raise ValueError(
        f"Unrecognized abundance file format.  Columns: {sorted(cols)}.\n"
        "Expected one of:\n"
        "  salmon quant.sf  — Name, TPM\n"
        "  kallisto abundance.tsv — target_id, tpm\n"
        "  sim TSV — transcript_id, mrna_abundance [, nrna_abundance]"
    )


def _load_abundance_map(
    path: str | Path,
) -> tuple[dict[str, tuple[float, float | None]], str]:
    """Load transcript abundances → {transcript_id: (total_rna, nrna|None)}.

    Returns
    -------
    abund_map : dict[str, (float, float | None)]
        Keys are transcript IDs.  The second element is ``None`` when
        nRNA information is not available (salmon, kallisto, or sim TSV
        without ``nrna_abundance``).
    fmt : str
        Detected format ("salmon", "kallisto", or "sim").
    """
    fmt, df = _detect_abundance_format(path)

    abund_map: dict[str, tuple[float, float | None]] = {}

    if fmt == "salmon":
        for _, row in df.iterrows():
            abund_map[str(row["Name"])] = (float(row["TPM"]), None)
    elif fmt == "kallisto":
        for _, row in df.iterrows():
            abund_map[str(row["target_id"])] = (float(row["tpm"]), None)
    elif fmt == "sim":
        has_nrna = "nrna_abundance" in df.columns
        for _, row in df.iterrows():
            mrna = float(row["mrna_abundance"])
            nrna: float | None = float(row["nrna_abundance"]) if has_nrna else None
            abund_map[str(row["transcript_id"])] = (mrna, nrna)

    return abund_map, fmt


def _spike_in_nrna(
    transcripts: list[Transcript],
    nrna_frac_range: tuple[float, float],
    seed: int = 42,
) -> None:
    """Add random nRNA as a fraction of each transcript's total abundance.

    For each transcript with abundance > 0:
        nrna_frac ~ Uniform(nrna_frac_min, nrna_frac_max)
        nRNA = total × nrna_frac
        mRNA = total × (1 − nrna_frac)

    Single-exon transcripts always get nRNA = 0.
    """
    nrna_frac_min, nrna_frac_max = nrna_frac_range
    rng = np.random.default_rng(seed)
    n_spiked = 0
    n_single = 0

    for t in transcripts:
        total = (t.abundance or 0.0) + t.nrna_abundance
        if total <= 0 or len(t.exons) <= 1:
            if len(t.exons) <= 1 and total > 0:
                t.abundance = total
                t.nrna_abundance = 0.0
                n_single += 1
            continue
        frac = rng.uniform(nrna_frac_min, nrna_frac_max)
        t.nrna_abundance = total * frac
        t.abundance = total * (1.0 - frac)
        n_spiked += 1

    total_nrna = sum(t.nrna_abundance for t in transcripts)
    logger.info(
        "Spiked nRNA: %d transcripts, total nRNA=%.1f "
        "(frac range [%.2f, %.2f]), %d single-exon zeroed",
        n_spiked, total_nrna,
        nrna_frac_min, nrna_frac_max, n_single,
    )


def assign_file_abundances(
    transcripts: list[Transcript],
    tsv_path: str | Path,
) -> bool:
    """Load abundances from a file (salmon, kallisto, or sim TSV).

    Auto-detects file format:

    - **salmon** ``quant.sf`` — uses ``TPM`` as total RNA abundance.
    - **kallisto** ``abundance.tsv`` — uses ``tpm`` as total RNA abundance.
    - **sim TSV** — uses ``mrna_abundance`` (and optionally
      ``nrna_abundance``).

    Returns ``True`` if the file provided explicit nRNA data (sim TSV
    with ``nrna_abundance`` column), ``False`` otherwise.  When the file
    does not supply nRNA, the caller applies the ``nrna.fracs`` sweep.
    """
    abund_map, fmt = _load_abundance_map(tsv_path)
    logger.info("Detected abundance format: %s (%d entries)", fmt, len(abund_map))

    matched = 0
    has_nrna_data = False

    for t in transcripts:
        if t.t_id in abund_map:
            total_or_mrna, nrna = abund_map[t.t_id]
            if nrna is not None:
                # File provided both mRNA and nRNA
                t.abundance = total_or_mrna
                t.nrna_abundance = nrna
                has_nrna_data = True
            else:
                # File provided total RNA only (salmon/kallisto TPM)
                t.abundance = total_or_mrna
                t.nrna_abundance = 0.0
            matched += 1
        else:
            t.abundance = 0.0
            t.nrna_abundance = 0.0

    logger.info(
        "File abundances (%s): matched=%d/%d, has_nrna=%s",
        fmt, matched, len(transcripts), has_nrna_data,
    )

    return has_nrna_data


def write_truth_abundances(
    transcripts: list[Transcript],
    path: Path,
) -> None:
    """Write ground-truth abundances to a TSV file."""
    with open(path, "w") as f:
        f.write("transcript_id\tgene_id\tgene_name\tref\tstrand\t"
                "mrna_abundance\tnrna_abundance\ttotal_rna\tn_exons\t"
                "spliced_length\tgenomic_span\n")
        for t in transcripts:
            total = (t.abundance or 0.0) + t.nrna_abundance
            genomic_span = t.end - t.start if t.end and t.start else 0
            strand_str = t.strand.to_str()
            f.write(
                f"{t.t_id}\t{t.g_id}\t{t.g_name}\t{t.ref}\t"
                f"{strand_str}\t"
                f"{t.abundance:.4f}\t{t.nrna_abundance:.4f}\t"
                f"{total:.4f}\t{len(t.exons)}\t"
                f"{t.length or t.compute_length()}\t{genomic_span}\n"
            )
    logger.info("Wrote truth abundances to %s", path)


# ═══════════════════════════════════════════════════════════════════
# Whole-genome read simulator
# ═══════════════════════════════════════════════════════════════════


class WholeGenomeSimulator:
    """Whole-genome paired-end RNA-seq read simulator.

    Generates gzipped FASTQ and collated oracle BAM in a single pass.

    Architecture
    ------------
    1. Pre-compute abundance weight vectors (numpy arrays).
    2. Sample from fragment-length distribution, count fragments per
       unique length.
    3. For each fragment length, sample mRNA/nRNA transcripts
       proportional to abundance × effective_length.  Accumulate
       per-transcript fragment counts.
    4. Iterate transcripts: pull sequence ONCE, generate all
       fragments, write FASTQ.gz + BAM simultaneously.
    """

    _OVERSAMPLE_RATIO = 1.5
    _OVERSAMPLE_EXTRA = 10

    def __init__(
        self,
        fasta_path: str | Path,
        transcripts: list[Transcript],
        sim_params: SimulationParams,
        gdna_config: GDNASimConfig,
        *,
        strand_specificity: float = 1.0,
        seed: int | None = None,
    ):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.transcripts = transcripts
        self.sim_params = sim_params
        self.gdna_config = gdna_config
        self.strand_specificity = strand_specificity
        self._rng = np.random.default_rng(
            seed if seed is not None else sim_params.sim_seed,
        )

        N = len(transcripts)

        # Pre-extract spliced mRNA sequences
        logger.info("Pre-extracting %d mRNA sequences...", N)
        self._t_seqs: list[str] = []
        self._t_lengths: list[int] = []
        for t in transcripts:
            seq = self._extract_mrna_seq(t)
            self._t_seqs.append(seq)
            self._t_lengths.append(len(seq))

        # Vectorised abundance / length arrays for weight computation
        self._mrna_abund = np.array(
            [t.abundance or 0.0 for t in transcripts],
        )
        self._nrna_abund = np.array(
            [t.nrna_abundance for t in transcripts],
        )
        self._t_lengths_arr = np.array(self._t_lengths, dtype=np.int64)
        self._premrna_lengths: list[int] = [
            t.end - t.start for t in transcripts
        ]
        self._premrna_lengths_arr = np.array(
            self._premrna_lengths, dtype=np.int64,
        )

        # Clamp frag_max to longest transcript
        max_mrna = int(self._t_lengths_arr.max()) if N > 0 else 0
        self._frag_max = (
            min(sim_params.frag_max, max_mrna) if max_mrna > 0
            else sim_params.frag_max
        )

        # gDNA setup: annotated chromosomes weighted by length
        annotated_chroms = {t.ref for t in transcripts}
        self._gdna_chroms: list[str] = []
        self._gdna_chrom_lengths: list[int] = []
        for ref in self.fasta.references:
            if ref in annotated_chroms:
                self._gdna_chroms.append(ref)
                self._gdna_chrom_lengths.append(
                    self.fasta.get_reference_length(ref),
                )

        # BAM header and reference mapping (built once, reused)
        self._ref_names = list(self.fasta.references)
        self._ref_lengths = [
            self.fasta.get_reference_length(r) for r in self._ref_names
        ]
        self._ref_name_to_id: dict[str, int] = {
            name: i for i, name in enumerate(self._ref_names)
        }
        self._bam_header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
            "SQ": [
                {"SN": name, "LN": length}
                for name, length in zip(self._ref_names, self._ref_lengths)
            ],
            "PG": [{
                "ID": "rigel_sim",
                "PN": "rigel_sim",
                "VN": "1.0",
                "CL": "simulated",
            }],
        })

        logger.info(
            "Simulator ready: %d transcripts, %d gDNA chromosomes, "
            "max mRNA length=%d",
            N, len(self._gdna_chroms), max_mrna,
        )

    # -- Sequence extraction ------------------------------------------------

    def _extract_mrna_seq(self, t: Transcript) -> str:
        """Extract spliced mRNA sequence (5'-to-3' oriented)."""
        exon_seqs = []
        for e in t.exons:
            seq = self.fasta.fetch(t.ref, e.start, e.end)
            exon_seqs.append(seq.upper())
        seq = "".join(exon_seqs)
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    def _fetch_premrna_seq(self, t: Transcript) -> str:
        """Fetch unspliced pre-mRNA sequence on demand."""
        seq = self.fasta.fetch(t.ref, t.start, t.end).upper()
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    # -- Fragment length sampling -------------------------------------------

    def _sample_frag_lengths(
        self, n: int, mean: float, std: float,
        frag_min: int, frag_max: int,
    ) -> np.ndarray:
        """Sample *n* fragment lengths from a truncated normal."""
        rng = self._rng
        result = np.empty(n, dtype=int)
        filled = 0
        while filled < n:
            needed = n - filled
            size = int(needed * self._OVERSAMPLE_RATIO) + self._OVERSAMPLE_EXTRA
            raw = rng.normal(mean, std, size).astype(int)
            valid = raw[(raw >= frag_min) & (raw <= frag_max)]
            nkeep = min(len(valid), needed)
            result[filled : filled + nkeep] = valid[:nkeep]
            filled += nkeep
        return result

    def _sample_rna_frag_lengths(self, n: int) -> np.ndarray:
        p = self.sim_params
        return self._sample_frag_lengths(
            n, p.frag_mean, p.frag_std, p.frag_min, self._frag_max,
        )

    def _sample_gdna_frag_lengths(self, n: int) -> np.ndarray:
        g = self.gdna_config
        return self._sample_frag_lengths(
            n, g.frag_mean, g.frag_std, g.frag_min, g.frag_max,
        )

    # -- Error introduction -------------------------------------------------

    def _introduce_errors(self, seq: str) -> str:
        if self.sim_params.error_rate <= 0:
            return seq
        arr = np.array(list(seq), dtype="<U1")
        mask = self._rng.random(len(seq)) < self.sim_params.error_rate
        if not np.any(mask):
            return seq
        arr[mask] = _DNA_BASES[self._rng.integers(4, size=int(mask.sum()))]
        return "".join(arr)

    # -- Steps 1-3: accumulate per-transcript fragment counts ---------------

    def _accumulate_rna_counts(
        self, n_rna: int,
    ) -> tuple[dict[int, dict[int, int]], dict[int, dict[int, int]]]:
        """Sample fragment lengths and transcript assignments.

        Returns
        -------
        mrna_counts : dict[t_idx, dict[frag_len, count]]
        nrna_counts : dict[t_idx, dict[frag_len, count]]
        """
        if n_rna <= 0:
            return {}, {}

        N = len(self.transcripts)
        rng = self._rng

        # Step 2: sample all RNA fragment lengths up front
        frag_lengths = self._sample_rna_frag_lengths(n_rna)
        unique_lengths, length_counts = np.unique(
            frag_lengths, return_counts=True,
        )

        mrna_counts: dict[int, dict[int, int]] = defaultdict(
            lambda: defaultdict(int),
        )
        nrna_counts: dict[int, dict[int, int]] = defaultdict(
            lambda: defaultdict(int),
        )

        # Step 3: for each fragment length, sample transcripts
        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)

            # Vectorised weight computation
            mrna_eff = np.maximum(0, self._t_lengths_arr - fl + 1)
            nrna_eff = np.maximum(0, self._premrna_lengths_arr - fl + 1)
            weights = np.concatenate([
                self._mrna_abund * mrna_eff,
                self._nrna_abund * nrna_eff,
            ])
            total_w = weights.sum()
            if total_w <= 0:
                continue

            probs = weights / total_w
            indices = rng.choice(2 * N, size=fc, p=probs)
            unique_idx, idx_counts = np.unique(indices, return_counts=True)

            for idx, cnt in zip(unique_idx, idx_counts):
                idx, cnt = int(idx), int(cnt)
                if idx < N:
                    mrna_counts[idx][fl] += cnt
                else:
                    nrna_counts[idx - N][fl] += cnt

        return dict(mrna_counts), dict(nrna_counts)

    # -- Step 4: generate reads and write FASTQ + BAM ----------------------

    def _write_mrna_bam_pair(
        self,
        bam_fh: pysam.AlignmentFile,
        qname: str,
        r1_seq: str,
        r2_seq: str,
        t: Transcript,
        frag_start: int,
        frag_end: int,
        read_len: int,
        flipped: bool,
        ref_id: int,
    ) -> None:
        """Write one mRNA R1+R2 BAM record pair (collated)."""
        r2_t_end = min(frag_start + read_len, frag_end)
        r2_blocks = _transcript_to_genomic_blocks(frag_start, r2_t_end, t)
        r1_t_start = max(frag_end - read_len, frag_start)
        r1_blocks = _transcript_to_genomic_blocks(r1_t_start, frag_end, t)
        if not r2_blocks or not r1_blocks:
            return

        r2_cigar = _blocks_to_cigar(r2_blocks)
        r1_cigar = _blocks_to_cigar(r1_blocks)
        r2_start = r2_blocks[0][0]
        r1_start = r1_blocks[0][0]

        if t.strand == Strand.POS:
            r2_is_rev, r1_is_rev = False, True
        else:
            r2_is_rev, r1_is_rev = True, False
        if flipped:
            r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev

        leftmost = min(r1_blocks[0][0], r2_blocks[0][0])
        rightmost = max(r1_blocks[-1][1], r2_blocks[-1][1])
        tlen = rightmost - leftmost

        tags: list = [("NH", 1)]
        if len(r2_blocks) > 1 or len(r1_blocks) > 1:
            xs_strand = "-" if t.strand == Strand.NEG else "+"
            tags.append(("XS", xs_strand))

        r1_flag = _BASE_R1_FLAG
        if r1_is_rev:
            r1_flag |= _FLAG_REVERSE
        if r2_is_rev:
            r1_flag |= _FLAG_MATE_REVERSE
        r2_flag = _BASE_R2_FLAG
        if r2_is_rev:
            r2_flag |= _FLAG_REVERSE
        if r1_is_rev:
            r2_flag |= _FLAG_MATE_REVERSE

        if r1_start <= r2_start:
            r1_tlen, r2_tlen = tlen, -tlen
        else:
            r1_tlen, r2_tlen = -tlen, tlen

        bam_fh.write(_make_bam_record(
            self._bam_header, qname, r1_seq, r1_flag, ref_id,
            r1_start, r1_cigar, ref_id, r2_start, r1_tlen, tags=tags,
        ))
        bam_fh.write(_make_bam_record(
            self._bam_header, qname, r2_seq, r2_flag, ref_id,
            r2_start, r2_cigar, ref_id, r1_start, r2_tlen, tags=tags,
        ))

    def _write_nrna_bam_pair(
        self,
        bam_fh: pysam.AlignmentFile,
        qname: str,
        r1_seq: str,
        r2_seq: str,
        t: Transcript,
        frag_start: int,
        frag_end: int,
        read_len: int,
        flipped: bool,
        ref_id: int,
    ) -> None:
        """Write one nRNA R1+R2 BAM record pair (collated)."""
        g_start, g_end = _premrna_to_genomic_interval(
            frag_start, frag_end, t,
        )
        r2_g_start = g_start
        r2_g_end = min(g_start + read_len, g_end)
        r1_g_end = g_end
        r1_g_start = max(g_end - read_len, g_start)

        if t.strand == Strand.POS:
            r2_is_rev, r1_is_rev = False, True
        else:
            r2_is_rev, r1_is_rev = True, False
        if flipped:
            r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev

        tlen = g_end - g_start
        tags: list = [("NH", 1)]

        r1_flag = _BASE_R1_FLAG
        if r1_is_rev:
            r1_flag |= _FLAG_REVERSE
        if r2_is_rev:
            r1_flag |= _FLAG_MATE_REVERSE
        r2_flag = _BASE_R2_FLAG
        if r2_is_rev:
            r2_flag |= _FLAG_REVERSE
        if r1_is_rev:
            r2_flag |= _FLAG_MATE_REVERSE

        r1_read_len = r1_g_end - r1_g_start
        r2_read_len = r2_g_end - r2_g_start

        if r1_g_start <= r2_g_start:
            r1_tlen, r2_tlen = tlen, -tlen
        else:
            r1_tlen, r2_tlen = -tlen, tlen

        bam_fh.write(_make_bam_record(
            self._bam_header, qname, r1_seq, r1_flag, ref_id,
            r1_g_start, [(pysam.CMATCH, r1_read_len)],
            ref_id, r2_g_start, r1_tlen, tags=tags,
        ))
        bam_fh.write(_make_bam_record(
            self._bam_header, qname, r2_seq, r2_flag, ref_id,
            r2_g_start, [(pysam.CMATCH, r2_read_len)],
            ref_id, r1_g_start, r2_tlen, tags=tags,
        ))

    def _write_mrna_reads(
        self,
        t_idx: int,
        len_counts: dict[int, int],
        r1_fh,
        r2_fh,
        bam_fh: pysam.AlignmentFile | None,
    ) -> int:
        """Generate and write all mRNA reads for one transcript."""
        t = self.transcripts[t_idx]
        t_seq = self._t_seqs[t_idx]
        t_len = self._t_lengths[t_idx]
        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh else None
        rng = self._rng
        ss = self.strand_specificity
        strand_char = "r" if t.strand == Strand.NEG else "f"
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            read_len = min(self.sim_params.read_length, frag_len)
            eff_len = t_len - frag_len + 1
            if eff_len <= 0:
                continue

            frag_starts = rng.integers(0, eff_len, size=count)
            flip_mask = rng.random(count) >= ss if ss < 1.0 else None
            quals = "I" * read_len

            for i in range(count):
                frag_start = int(frag_starts[i])
                frag_end = frag_start + frag_len
                frag_seq = t_seq[frag_start:frag_end]

                r1_seq = reverse_complement(frag_seq[-read_len:])
                r2_seq = frag_seq[:read_len]

                flipped = flip_mask is not None and flip_mask[i]
                if flipped:
                    r1_seq, r2_seq = r2_seq, r1_seq

                r1_seq = self._introduce_errors(r1_seq)
                r2_seq = self._introduce_errors(r2_seq)

                qname = f"{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{n_written}"

                r1_fh.write(f"@{qname}/1\n{r1_seq}\n+\n{quals}\n")
                r2_fh.write(f"@{qname}/2\n{r2_seq}\n+\n{quals}\n")

                if bam_fh is not None and ref_id is not None:
                    self._write_mrna_bam_pair(
                        bam_fh, qname, r1_seq, r2_seq, t,
                        frag_start, frag_end, read_len, flipped, ref_id,
                    )

                n_written += 1

        return n_written

    def _write_nrna_reads(
        self,
        t_idx: int,
        len_counts: dict[int, int],
        r1_fh,
        r2_fh,
        bam_fh: pysam.AlignmentFile | None,
    ) -> int:
        """Generate and write all nRNA reads for one transcript.

        Fetches the pre-mRNA sequence ONCE for this transcript.
        """
        t = self.transcripts[t_idx]
        premrna_len = self._premrna_lengths[t_idx]
        premrna_seq = self._fetch_premrna_seq(t)
        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh else None
        rng = self._rng
        ss = self.strand_specificity
        strand_char = "r" if t.strand == Strand.NEG else "f"
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            read_len = min(self.sim_params.read_length, frag_len)
            eff_len = premrna_len - frag_len + 1
            if eff_len <= 0:
                continue

            frag_starts = rng.integers(0, eff_len, size=count)
            flip_mask = rng.random(count) >= ss if ss < 1.0 else None
            quals = "I" * read_len

            for i in range(count):
                frag_start = int(frag_starts[i])
                frag_end = frag_start + frag_len
                frag_seq = premrna_seq[frag_start:frag_end]

                r1_seq = reverse_complement(frag_seq[-read_len:])
                r2_seq = frag_seq[:read_len]

                flipped = flip_mask is not None and flip_mask[i]
                if flipped:
                    r1_seq, r2_seq = r2_seq, r1_seq

                r1_seq = self._introduce_errors(r1_seq)
                r2_seq = self._introduce_errors(r2_seq)

                qname = f"nrna_{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{n_written}"

                r1_fh.write(f"@{qname}/1\n{r1_seq}\n+\n{quals}\n")
                r2_fh.write(f"@{qname}/2\n{r2_seq}\n+\n{quals}\n")

                if bam_fh is not None and ref_id is not None:
                    self._write_nrna_bam_pair(
                        bam_fh, qname, r1_seq, r2_seq, t,
                        frag_start, frag_end, read_len, flipped, ref_id,
                    )

                n_written += 1

        return n_written

    def _write_gdna_reads(
        self,
        n_gdna: int,
        r1_fh,
        r2_fh,
        bam_fh: pysam.AlignmentFile | None,
    ) -> int:
        """Generate and write all gDNA reads."""
        if n_gdna == 0 or not self._gdna_chroms:
            return 0

        frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
        unique_lengths, length_counts = np.unique(
            frag_lengths, return_counts=True,
        )
        rng = self._rng
        read_len_cfg = self.sim_params.read_length
        n_written = 0

        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            chrom_eff = np.array(
                [max(0, cl - fl + 1) for cl in self._gdna_chrom_lengths],
                dtype=float,
            )
            total_eff = chrom_eff.sum()
            if total_eff <= 0:
                continue
            chrom_probs = chrom_eff / total_eff
            chrom_indices = rng.choice(
                len(self._gdna_chroms), size=fc, p=chrom_probs,
            )
            unique_chroms, chrom_counts = np.unique(
                chrom_indices, return_counts=True,
            )

            for ci, cc in zip(unique_chroms, chrom_counts):
                ci, cc = int(ci), int(cc)
                chrom = self._gdna_chroms[ci]
                chrom_len = self._gdna_chrom_lengths[ci]
                ref_id = self._ref_name_to_id.get(chrom)
                eff_len = chrom_len - fl + 1
                starts = rng.integers(0, eff_len, size=cc)
                strands = rng.integers(0, 2, size=cc)

                for j in range(cc):
                    start = int(starts[j])
                    end = start + fl
                    read_len = min(read_len_cfg, fl)
                    seq = self.fasta.fetch(chrom, start, end).upper()

                    is_neg = bool(strands[j])
                    if is_neg:
                        seq = reverse_complement(seq)
                    strand_char = "r" if is_neg else "f"

                    r1_seq = reverse_complement(seq[-read_len:])
                    r2_seq = seq[:read_len]

                    r1_seq = self._introduce_errors(r1_seq)
                    r2_seq = self._introduce_errors(r2_seq)

                    qname = f"gdna:{chrom}:{start}-{end}:{strand_char}:{n_written}"
                    quals = "I" * read_len

                    r1_fh.write(f"@{qname}/1\n{r1_seq}\n+\n{quals}\n")
                    r2_fh.write(f"@{qname}/2\n{r2_seq}\n+\n{quals}\n")

                    if bam_fh is not None and ref_id is not None:
                        r1_flag = _BASE_R1_FLAG
                        r2_flag = _BASE_R2_FLAG
                        if is_neg:
                            r1_flag |= _FLAG_REVERSE
                            r2_flag |= _FLAG_MATE_REVERSE
                        else:
                            r1_flag |= _FLAG_MATE_REVERSE
                            r2_flag |= _FLAG_REVERSE

                        tlen = end - start
                        r2_start_pos = start
                        r1_start_pos = end - read_len

                        if r1_start_pos <= r2_start_pos:
                            r1_tlen, r2_tlen = tlen, -tlen
                        else:
                            r1_tlen, r2_tlen = -tlen, tlen

                        bam_fh.write(_make_bam_record(
                            self._bam_header, qname, r1_seq,
                            r1_flag, ref_id, r1_start_pos,
                            [(pysam.CMATCH, read_len)],
                            ref_id, r2_start_pos, r1_tlen,
                        ))
                        bam_fh.write(_make_bam_record(
                            self._bam_header, qname, r2_seq,
                            r2_flag, ref_id, r2_start_pos,
                            [(pysam.CMATCH, read_len)],
                            ref_id, r1_start_pos, r2_tlen,
                        ))

                    n_written += 1

        return n_written

    # -- Main entry point ---------------------------------------------------

    def simulate_and_write(
        self,
        output_dir: Path,
        n_rna: int,
        n_gdna: int = 0,
        *,
        oracle_bam: bool = True,
        prefix: str = "sim",
    ) -> tuple[Path, Path, Path | None]:
        """Single-pass simulation: accumulate counts, generate, write.

        Steps
        -----
        1-3. Accumulate per-transcript (mRNA/nRNA) fragment counts.
        4.   Iterate transcripts: pull sequence once, generate all
             fragments, write gzipped FASTQ + collated BAM.

        Returns
        -------
        (r1_path, r2_path, bam_path | None)
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        r1_path = output_dir / f"{prefix}_R1.fq.gz"
        r2_path = output_dir / f"{prefix}_R2.fq.gz"
        bam_path = output_dir / f"{prefix}_oracle.bam" if oracle_bam else None

        # Steps 1-3: accumulate per-transcript fragment counts
        logger.info("Accumulating RNA fragment counts...")
        mrna_counts, nrna_counts = self._accumulate_rna_counts(n_rna)

        n_mrna_tx = len(mrna_counts)
        n_nrna_tx = len(nrna_counts)
        total_mrna = sum(sum(d.values()) for d in mrna_counts.values())
        total_nrna = sum(sum(d.values()) for d in nrna_counts.values())
        logger.info(
            "Fragment counts: %d mRNA frags across %d transcripts, "
            "%d nRNA frags across %d transcripts",
            total_mrna, n_mrna_tx, total_nrna, n_nrna_tx,
        )

        # Step 4: generate reads and write outputs
        n_written = 0
        bam_fh: pysam.AlignmentFile | None = None

        with (
            gzip.open(r1_path, "wt") as r1_fh,
            gzip.open(r2_path, "wt") as r2_fh,
        ):
            if oracle_bam:
                bam_fh = pysam.AlignmentFile(
                    str(bam_path), "wb", header=self._bam_header,
                )
            try:
                # mRNA reads (sequence already pre-extracted)
                for t_idx in sorted(mrna_counts):
                    n_written += self._write_mrna_reads(
                        t_idx, mrna_counts[t_idx],
                        r1_fh, r2_fh, bam_fh,
                    )

                # nRNA reads (pre-mRNA fetched once per transcript)
                for t_idx in sorted(nrna_counts):
                    n_written += self._write_nrna_reads(
                        t_idx, nrna_counts[t_idx],
                        r1_fh, r2_fh, bam_fh,
                    )

                # gDNA reads
                n_written += self._write_gdna_reads(
                    n_gdna, r1_fh, r2_fh, bam_fh,
                )
            finally:
                if bam_fh is not None:
                    bam_fh.close()

        logger.info(
            "Wrote %d read pairs -> %s (RNA=%d, gDNA=%d, oracle=%s)",
            n_written, output_dir, n_rna, n_gdna, oracle_bam,
        )
        return r1_path, r2_path, bam_path

    def close(self) -> None:
        """Close the underlying FASTA file handle."""
        self.fasta.close()


# ═══════════════════════════════════════════════════════════════════
# BAM record helpers
# ═══════════════════════════════════════════════════════════════════


def _transcript_to_genomic_blocks(
    frag_start: int,
    frag_end: int,
    transcript: Transcript,
) -> list[tuple[int, int]]:
    """Map transcript-space interval to genomic exon blocks.

    For negative-strand transcripts, mirrors coordinates before
    mapping (the transcript sequence is reverse-complemented so
    position 0 is the mRNA 5′ end / rightmost genomic coordinate).
    """
    exons = transcript.exons
    if transcript.strand == Strand.NEG:
        t_len = sum(e.end - e.start for e in exons)
        frag_start, frag_end = t_len - frag_end, t_len - frag_start

    blocks: list[tuple[int, int]] = []
    consumed = 0
    for exon in exons:
        exon_len = exon.end - exon.start
        exon_tx_start = consumed
        exon_tx_end = consumed + exon_len
        overlap_start = max(frag_start, exon_tx_start)
        overlap_end = min(frag_end, exon_tx_end)
        if overlap_start < overlap_end:
            offset_start = overlap_start - exon_tx_start
            offset_end = overlap_end - exon_tx_start
            blocks.append((exon.start + offset_start, exon.start + offset_end))
        consumed += exon_len
        if consumed >= frag_end:
            break
    return blocks


def _premrna_to_genomic_interval(
    frag_start: int,
    frag_end: int,
    transcript: Transcript,
) -> tuple[int, int]:
    """Map pre-mRNA-space interval to genomic coordinates."""
    premrna_len = transcript.end - transcript.start
    if transcript.strand == Strand.NEG:
        frag_start, frag_end = premrna_len - frag_end, premrna_len - frag_start
    return (transcript.start + frag_start, transcript.start + frag_end)


def _blocks_to_cigar(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Convert genomic blocks to pysam CIGAR tuples."""
    cigar: list[tuple[int, int]] = []
    for i, (bstart, bend) in enumerate(blocks):
        if i > 0:
            prev_end = blocks[i - 1][1]
            intron_len = bstart - prev_end
            if intron_len > 0:
                cigar.append((pysam.CREF_SKIP, intron_len))
        match_len = bend - bstart
        if match_len > 0:
            cigar.append((pysam.CMATCH, match_len))
    return cigar


def _make_bam_record(
    header: pysam.AlignmentHeader,
    query_name: str,
    query_sequence: str,
    flag: int,
    reference_id: int,
    reference_start: int,
    cigar: list[tuple[int, int]],
    mate_reference_id: int,
    mate_reference_start: int,
    template_length: int,
    mapping_quality: int = 255,
    tags: list | None = None,
) -> pysam.AlignedSegment:
    """Build a pysam AlignedSegment."""
    a = pysam.AlignedSegment(header)
    a.query_name = query_name
    a.query_sequence = query_sequence
    a.flag = flag
    a.reference_id = reference_id
    a.reference_start = reference_start
    a.cigar = cigar
    a.mapping_quality = mapping_quality
    a.query_qualities = pysam.qualitystring_to_array("I" * len(query_sequence))
    a.next_reference_id = mate_reference_id
    a.next_reference_start = mate_reference_start
    a.template_length = template_length
    if tags:
        a.set_tags(tags)
    return a


# ═══════════════════════════════════════════════════════════════════
# Condition naming helpers
# ═══════════════════════════════════════════════════════════════════


def condition_dir_name(
    gdna_label: str,
    strand_specificity: float,
    nrna_label: str,
) -> str:
    return f"gdna_{gdna_label}_ss_{strand_specificity:.2f}_nrna_{nrna_label}"


def gdna_label_for_rate(rate: float, labels: list[str] | None, idx: int) -> str:
    if labels and idx < len(labels):
        return labels[idx].strip()
    return f"r{rate:g}"


def nrna_label_for_frac(
    frac_range: tuple[float, float],
    labels: list[str] | None,
    idx: int,
) -> str:
    if labels and idx < len(labels):
        return labels[idx].strip()
    lo, hi = frac_range
    if lo == hi == 0.0:
        return "none"
    return f"{lo:.2f}_{hi:.2f}"


def _build_nrna_pairs(
    cfg: SimConfig,
    has_file_nrna: bool,
) -> list[tuple[str, tuple[float, float] | None]]:
    """Build nRNA sweep pairs.

    When the abundance file supplied explicit nRNA data, returns a
    single entry ``("file", None)`` — no spike-in.  Otherwise returns
    one entry per ``cfg.nrna.fracs``.
    """
    if has_file_nrna:
        return [("file", None)]
    pairs: list[tuple[str, tuple[float, float] | None]] = []
    for i, frac_range in enumerate(cfg.nrna.fracs):
        label = nrna_label_for_frac(
            frac_range, cfg.nrna.frac_labels, i,
        )
        pairs.append((label, frac_range))
    return pairs


# ═══════════════════════════════════════════════════════════════════
# Manifest
# ═══════════════════════════════════════════════════════════════════


def write_manifest(
    outdir: Path,
    cfg: SimConfig,
    conditions: list[dict],
) -> None:
    """Write manifest.json summarizing all simulation outputs."""
    manifest = {
        "version": 1,
        "genome": str(Path(cfg.genome).resolve()),
        "gtf": str(Path(cfg.gtf).resolve()),
        "transcript_filter": cfg.transcript_filter,
        "truth_abundances": "truth_abundances.tsv",
        "simulation": asdict(cfg.simulation),
        "abundance": asdict(cfg.abundance),
        "conditions": conditions,
    }
    path = outdir / "manifest.json"
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("Wrote manifest to %s", path)


# ═══════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════


def run_simulation(cfg: SimConfig) -> list[dict]:
    """Run full simulation.

    1. Load genome + GTF → transcripts
    2. Assign base abundances (total RNA) once
    3. Sweep nrna_fracs × gdna_rates × strand_specificities
    4. Write manifest

    When the abundance file provides explicit nRNA data, the nRNA
    sweep is skipped — file values are used as-is.  Otherwise each
    entry in ``nrna.fracs`` produces a separate nRNA condition via
    ``_spike_in_nrna``.

    Returns list of condition dicts (for manifest).
    """
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = Path(cfg.genome)
    gtf_path = Path(cfg.gtf)

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome not found: {genome_path}")
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    # 1. Load transcripts
    transcripts = load_transcripts(gtf_path, transcript_filter=cfg.transcript_filter)
    if not transcripts:
        raise RuntimeError("No transcripts loaded from GTF")

    # 2. Assign base abundances (total RNA, nRNA = 0)
    ab = cfg.abundance
    has_file_nrna = False
    if ab.mode == "random":
        assign_random_abundances(transcripts, ab)
    elif ab.mode == "file":
        if not ab.file:
            raise ValueError("abundance.file must be specified for mode='file'")
        has_file_nrna = assign_file_abundances(transcripts, ab.file)
    else:
        raise ValueError(f"Unknown abundance mode: {ab.mode}")

    if has_file_nrna:
        logger.info(
            "Abundance file provided nRNA data — skipping nRNA spike-in sweep"
        )

    # Save base abundances
    base_abundances = [
        (t.abundance or 0.0, t.nrna_abundance) for t in transcripts
    ]

    # 3. Build condition grid: nrna × gdna × strand_specificities
    sim = cfg.simulation

    gdna_pairs: list[tuple[str, float]] = []
    for i, rate in enumerate(cfg.gdna.rates):
        label = gdna_label_for_rate(rate, cfg.gdna.rate_labels, i)
        gdna_pairs.append((label, rate))

    nrna_pairs = _build_nrna_pairs(cfg, has_file_nrna)

    total_conditions = (
        len(nrna_pairs) * len(gdna_pairs) * len(cfg.strand_specificities)
    )
    cond_num = 0
    conditions: list[dict] = []

    for nrna_label, nrna_frac_range in nrna_pairs:
        # Deep-copy transcripts for this nRNA configuration
        cond_transcripts = copy.deepcopy(transcripts)

        # Restore base abundances
        for t, (base_mrna, base_nrna) in zip(
            cond_transcripts, base_abundances,
        ):
            t.abundance = base_mrna
            t.nrna_abundance = base_nrna

        # Apply nRNA spike-in (skipped when file provided nRNA)
        nrna_lo = nrna_hi = 0.0
        if nrna_frac_range is not None:
            nrna_lo, nrna_hi = nrna_frac_range
            if nrna_lo > 0.0 or nrna_hi > 0.0:
                _spike_in_nrna(
                    cond_transcripts, nrna_frac_range, seed=ab.seed,
                )

        # Write truth abundances for this nRNA configuration
        truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        write_truth_abundances(cond_transcripts, outdir / truth_name)

        for gdna_label, gdna_rate in gdna_pairs:
            for strand_spec in cfg.strand_specificities:
                cond_num += 1
                n_rna = sim.n_rna_fragments
                n_gdna = round(gdna_rate * n_rna)

                cond_name = condition_dir_name(
                    gdna_label, strand_spec, nrna_label,
                )
                cond_dir = outdir / cond_name

                print(
                    f"\n[{cond_num}/{total_conditions}] {cond_name}: "
                    f"RNA={n_rna:,} gDNA={n_gdna:,} SS={strand_spec:.2f} "
                    f"nRNA=[{nrna_lo:.2f},{nrna_hi:.2f}]",
                    flush=True,
                )

                cond_entry: dict = {
                    "name": cond_name,
                    "gdna_label": gdna_label,
                    "gdna_rate": gdna_rate,
                    "strand_specificity": strand_spec,
                    "nrna_label": nrna_label,
                    "nrna_frac_min": nrna_lo,
                    "nrna_frac_max": nrna_hi,
                    "n_rna": n_rna,
                    "n_gdna": n_gdna,
                    "truth_abundances": truth_name,
                    "fastq_r1": f"{cond_name}/sim_R1.fq.gz",
                    "fastq_r2": f"{cond_name}/sim_R2.fq.gz",
                }

                # Check if already done
                if (cond_dir / "sim_R1.fq.gz").exists():
                    print("  Output exists, skipping", flush=True)
                    cond_entry["oracle_bam"] = (
                        f"{cond_name}/sim_oracle.bam" if cfg.oracle_bam
                        else None
                    )
                else:
                    print("  Simulating...", end="", flush=True)
                    t0 = time.monotonic()
                    simulator = WholeGenomeSimulator(
                        genome_path, cond_transcripts, sim, cfg.gdna,
                        strand_specificity=strand_spec,
                    )
                    _, _, bam_path = simulator.simulate_and_write(
                        cond_dir, n_rna, n_gdna,
                        oracle_bam=cfg.oracle_bam, prefix="sim",
                    )
                    simulator.close()
                    cond_entry["oracle_bam"] = (
                        f"{cond_name}/sim_oracle.bam" if bam_path
                        else None
                    )
                    print(f" done ({time.monotonic() - t0:.1f}s)", flush=True)

                conditions.append(cond_entry)
    # 4. Write manifest
    write_manifest(outdir, cfg, conditions)

    return conditions


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Whole-genome RNA-seq read simulator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", required=True, help="YAML configuration file")
    p.add_argument("--genome", help="Genome FASTA (overrides YAML)")
    p.add_argument("--gtf", help="Gene annotation GTF (overrides YAML)")
    p.add_argument("--outdir", help="Output directory (overrides YAML)")
    p.add_argument("--n-rna", type=int, help="Number of RNA fragments (overrides YAML)")
    p.add_argument("--no-oracle", action="store_true", help="Skip oracle BAM generation")
    p.add_argument("--verbose", action="store_true", default=None, help="Verbose logging")
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    cfg = parse_yaml_config(args.config)

    # CLI overrides
    if args.genome:
        cfg.genome = args.genome
    if args.gtf:
        cfg.gtf = args.gtf
    if args.outdir:
        cfg.outdir = args.outdir
    if args.n_rna is not None:
        cfg.simulation.n_rna_fragments = args.n_rna
    if args.no_oracle:
        cfg.oracle_bam = False
    if args.verbose is not None:
        cfg.verbose = args.verbose

    level = logging.DEBUG if cfg.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    if not cfg.genome:
        print("Error: genome FASTA not specified", file=sys.stderr)
        return 1
    if not cfg.gtf:
        print("Error: GTF not specified", file=sys.stderr)
        return 1

    print("Whole-genome RNA-seq simulator", flush=True)
    print(f"  Genome:           {cfg.genome}", flush=True)
    print(f"  GTF:              {cfg.gtf}", flush=True)
    print(f"  Output:           {cfg.outdir}", flush=True)
    print(f"  RNA fragments:    {cfg.simulation.n_rna_fragments:,}", flush=True)
    print(f"  gDNA rates:       {cfg.gdna.rates}", flush=True)
    print(f"  Strand specs:     {cfg.strand_specificities}", flush=True)
    print(f"  nRNA fracs:       {cfg.nrna.fracs}", flush=True)
    print(f"  Transcript filter:{cfg.transcript_filter}", flush=True)
    print(f"  Oracle BAM:       {cfg.oracle_bam}", flush=True)

    t0 = time.monotonic()
    conditions = run_simulation(cfg)
    elapsed = time.monotonic() - t0

    print(f"\nSimulation complete in {elapsed:.1f}s", flush=True)
    print(f"  {len(conditions)} conditions generated → {cfg.outdir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
