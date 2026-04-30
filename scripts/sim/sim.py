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
import multiprocessing
import shutil
import sys
import time
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pysam

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]

try:
    import pgzip
except ImportError:
    pgzip = None  # type: ignore[assignment]

from rigel.transcript import Transcript
from rigel.types import Strand

logger = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────

_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Byte-level complement table for vectorized reverse-complement
_BYTE_COMPLEMENT = np.zeros(256, dtype=np.uint8)
for _c, _rc in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _BYTE_COMPLEMENT[_c] = _rc

_FASTQ_BUFFER_SIZE = 100_000

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


def _seq_to_bytes(seq: str) -> np.ndarray:
    """Convert a DNA string to a numpy uint8 array (writable copy)."""
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()


def _batch_extract_reads(
    seq_bytes: np.ndarray,
    frag_starts: np.ndarray,
    frag_len: int,
    read_len: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized extraction of R1 and R2 read sequences.

    Given a source sequence as a byte array, extract all fragments at once
    using stride tricks, then slice R2 (first read_len bases) and R1
    (reverse complement of last read_len bases).

    Returns (r2_seqs, r1_seqs) as 2D uint8 arrays of shape (count, read_len).
    """
    offsets = np.arange(read_len, dtype=np.int64)

    # R2: first read_len bases of each fragment
    r2_indices = frag_starts[:, np.newaxis] + offsets[np.newaxis, :]
    r2_seqs = seq_bytes[r2_indices]

    # R1: reverse complement of last read_len bases
    r1_start = frag_starts + (frag_len - read_len)
    r1_indices = r1_start[:, np.newaxis] + offsets[np.newaxis, :]
    r1_seqs = _BYTE_COMPLEMENT[seq_bytes[r1_indices][:, ::-1]]

    return r2_seqs, r1_seqs


def _open_gzip_write(path: Path, threads: int = 4):
    """Open a gzip file for text writing, using pgzip if available."""
    if pgzip is not None:
        return pgzip.open(path, "wt", thread=threads, compresslevel=3)
    return gzip.open(path, "wt", compresslevel=3)


class _FastqBuffer:
    """Buffered FASTQ writer that accumulates records before flushing.

    Reduces the number of gzip write calls from ~8 per fragment to
    ~8 per flush (every ``buf_size`` fragments).
    """

    __slots__ = ("_fh", "_parts", "_size", "_count")

    def __init__(self, fh, buf_size: int = _FASTQ_BUFFER_SIZE):
        self._fh = fh
        self._parts: list[str] = []
        self._size = buf_size
        self._count = 0

    def write_record(self, qname: str, seq: str, quals: str) -> None:
        self._parts.append(f"@{qname}\n{seq}\n+\n{quals}\n")
        self._count += 1
        if self._count >= self._size:
            self.flush()

    def write_records_batch(
        self,
        qnames: list[str],
        seqs: np.ndarray,
        quals: str,
        suffix: str,
    ) -> None:
        """Write a batch of FASTQ records from numpy byte arrays.

        Parameters
        ----------
        qnames : list of query names (without /1 or /2 suffix)
        seqs : 2D uint8 array (count x read_len)
        quals : quality string (same for all reads)
        suffix : "/1" or "/2"
        """
        parts = self._parts
        for i in range(len(qnames)):
            seq_str = seqs[i].tobytes().decode("ascii")
            parts.append(f"@{qnames[i]}{suffix}\n{seq_str}\n+\n{quals}\n")
        self._count += len(qnames)
        if self._count >= self._size:
            self.flush()

    def flush(self) -> None:
        if self._parts:
            self._fh.write("".join(self._parts))
            self._parts.clear()
            self._count = 0

    def close(self) -> None:
        self.flush()


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
    n_workers: int = 1


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
    sim.n_workers = int(sim_raw.get("n_workers", 1))

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
# Genome cache — pre-load chromosome sequences as numpy byte arrays
# ═══════════════════════════════════════════════════════════════════


class GenomeCache:
    """Caches full chromosome sequences as numpy byte arrays.

    Eliminates per-fragment ``pysam.FastaFile.fetch()`` overhead which
    dominates gDNA simulation time.  Chromosomes are loaded lazily on
    first access and stored as uppercase uint8 arrays.
    """

    __slots__ = ("_fasta", "_cache")

    def __init__(self, fasta: pysam.FastaFile):
        self._fasta = fasta
        self._cache: dict[str, np.ndarray] = {}

    def get(self, ref: str) -> np.ndarray:
        arr = self._cache.get(ref)
        if arr is None:
            seq = self._fasta.fetch(ref).upper()
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            self._cache[ref] = arr
        return arr

    def preload(self, refs: list[str]) -> None:
        for ref in refs:
            self.get(ref)

    def clear(self) -> None:
        self._cache.clear()


# ═══════════════════════════════════════════════════════════════════
# Whole-genome read simulator
# ═══════════════════════════════════════════════════════════════════


class WholeGenomeSimulator:
    """Whole-genome paired-end RNA-seq read simulator.

    Generates gzipped FASTQ and collated oracle BAM in a single pass.

    Performance strategy
    --------------------
    1. Pre-cache all chromosome sequences as numpy byte arrays at init.
    2. Pre-extract all mRNA and pre-mRNA sequences as byte arrays.
    3. Sample fragment lengths in bulk, then sample transcript
       assignments per unique fragment length.
    4. Vectorized read extraction via numpy fancy indexing.
    5. Multi-threaded gzip compression via pgzip (if available).
    6. Unified RNA read writer handles both mRNA and nRNA.

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
        self._genome = GenomeCache(self.fasta)

        # Pre-extract spliced mRNA sequences as byte arrays
        logger.info("Pre-extracting %d mRNA sequences...", N)
        self._t_seq_bytes: list[np.ndarray] = []
        self._t_lengths = np.empty(N, dtype=np.int64)
        for i, t in enumerate(transcripts):
            seq = self._extract_mrna_seq(t)
            self._t_seq_bytes.append(_seq_to_bytes(seq))
            self._t_lengths[i] = len(seq)

        # Pre-mRNA lengths and sequences
        self._premrna_lengths = np.array(
            [t.end - t.start for t in transcripts], dtype=np.int64,
        )

        # Pre-extract pre-mRNA sequences for nRNA-eligible transcripts
        logger.info("Pre-extracting pre-mRNA sequences for nRNA...")
        self._premrna_seq_bytes: list[np.ndarray | None] = [None] * N
        for i, t in enumerate(transcripts):
            if t.nrna_abundance > 0 and len(t.exons) > 1:
                self._premrna_seq_bytes[i] = self._fetch_premrna_bytes(t)

        # Strand chars and abundance arrays
        self._t_strand_chars: list[str] = [
            "r" if t.strand == Strand.NEG else "f" for t in transcripts
        ]
        self._mrna_abund = np.array([t.abundance or 0.0 for t in transcripts])
        self._nrna_abund = np.array([t.nrna_abundance for t in transcripts])

        # Quality string cache
        self._quals_cache: dict[int, str] = {}

        # Fast-path flag for error introduction
        self._has_errors = sim_params.error_rate > 0

        # Clamp frag_max to longest transcript
        max_mrna = int(self._t_lengths.max()) if N > 0 else 0
        self._frag_max = min(sim_params.frag_max, max_mrna) if max_mrna > 0 else sim_params.frag_max

        # gDNA setup: annotated references weighted by length
        annotated_refs = {t.ref for t in transcripts}
        self._gdna_refs: list[str] = []
        self._gdna_ref_lengths: list[int] = []
        for ref in self.fasta.references:
            if ref in annotated_refs:
                self._gdna_refs.append(ref)
                self._gdna_ref_lengths.append(self.fasta.get_reference_length(ref))

        # Pre-load gDNA chromosomes into cache (major perf win)
        if self._gdna_refs:
            logger.info("Pre-loading %d chromosomes for gDNA...", len(self._gdna_refs))
            self._genome.preload(self._gdna_refs)

        # BAM header and reference mapping (built once, reused)
        self._ref_names = list(self.fasta.references)
        self._ref_lengths = [self.fasta.get_reference_length(r) for r in self._ref_names]
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

        # Pre-built shared tag lists for BAM records
        self._nh1_tags = [("NH", 1)]

        logger.info(
            "Simulator ready: %d transcripts, %d gDNA refs, max mRNA len=%d",
            N, len(self._gdna_refs), max_mrna,
        )

    # -- Sequence extraction ------------------------------------------------

    def _extract_mrna_seq(self, t: Transcript) -> str:
        """Extract spliced mRNA sequence (5'-to-3' oriented)."""
        exon_seqs = []
        for e in t.exons:
            exon_seqs.append(self.fasta.fetch(t.ref, e.start, e.end).upper())
        seq = "".join(exon_seqs)
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    def _fetch_premrna_bytes(self, t: Transcript) -> np.ndarray:
        """Fetch unspliced pre-mRNA sequence as uint8 array."""
        seq = self.fasta.fetch(t.ref, t.start, t.end).upper()
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return _seq_to_bytes(seq)

    def _get_premrna_bytes(self, t_idx: int) -> np.ndarray:
        """Get pre-mRNA bytes, fetching on demand if not pre-cached."""
        arr = self._premrna_seq_bytes[t_idx]
        if arr is None:
            arr = self._fetch_premrna_bytes(self.transcripts[t_idx])
            self._premrna_seq_bytes[t_idx] = arr
        return arr

    def _get_quals(self, read_len: int) -> str:
        """Return cached quality string of the given length."""
        q = self._quals_cache.get(read_len)
        if q is None:
            q = "I" * read_len
            self._quals_cache[read_len] = q
        return q

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

    def _introduce_errors_batch(self, seqs: np.ndarray) -> np.ndarray:
        """Apply random substitution errors to a batch of sequences in-place.

        Parameters
        ----------
        seqs : 2D uint8 array (count x read_len)
        """
        if not self._has_errors:
            return seqs
        count, read_len = seqs.shape
        mask = self._rng.random((count, read_len)) < self.sim_params.error_rate
        n_errors = int(mask.sum())
        if n_errors > 0:
            bases = np.array([65, 67, 71, 84], dtype=np.uint8)  # A, C, G, T
            seqs[mask] = bases[self._rng.integers(4, size=n_errors)]
        return seqs

    # -- Fragment count accumulation ----------------------------------------

    def _accumulate_pool(
        self,
        n_frags: int,
        abundances: np.ndarray,
        lengths: np.ndarray,
    ) -> dict[int, dict[int, int]]:
        """Sample fragment counts from a single abundance/length pool.

        Returns dict[t_idx, dict[frag_len, count]].
        """
        if n_frags <= 0:
            return {}

        rng = self._rng
        frag_lengths = self._sample_rna_frag_lengths(n_frags)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            eff = np.maximum(0, lengths - fl + 1)
            weights = abundances * eff
            total_w = weights.sum()
            if total_w <= 0:
                continue
            probs = weights / total_w
            indices = rng.choice(len(abundances), size=fc, p=probs)
            unique_idx, idx_counts = np.unique(indices, return_counts=True)
            for idx, cnt in zip(unique_idx, idx_counts):
                counts[int(idx)][fl] += int(cnt)

        return dict(counts)

    def _accumulate_rna_counts(
        self,
        n_rna: int,
        *,
        n_mrna: int | None = None,
        n_nrna: int | None = None,
    ) -> tuple[dict[int, dict[int, int]], dict[int, dict[int, int]]]:
        """Sample fragment lengths and transcript assignments.

        If *n_mrna* and *n_nrna* are given, samples each pool
        independently (fragment-count control).  Otherwise falls back
        to the combined-pool approach where *n_rna* fragments are
        drawn jointly from mRNA + nRNA weighted by abundance × length.

        Returns (mrna_counts, nrna_counts) where each is
        dict[t_idx, dict[frag_len, count]].
        """
        # Separate-pool mode: precise fragment-level control
        if n_mrna is not None and n_nrna is not None:
            mrna_counts = self._accumulate_pool(
                n_mrna, self._mrna_abund, self._t_lengths,
            )
            nrna_counts = self._accumulate_pool(
                n_nrna, self._nrna_abund, self._premrna_lengths,
            )
            return mrna_counts, nrna_counts

        # Combined-pool mode (original behaviour)
        if n_rna <= 0:
            return {}, {}

        N = len(self.transcripts)
        rng = self._rng

        frag_lengths = self._sample_rna_frag_lengths(n_rna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        mrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        nrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)

            mrna_eff = np.maximum(0, self._t_lengths - fl + 1)
            nrna_eff = np.maximum(0, self._premrna_lengths - fl + 1)
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

    # -- Unified RNA read writer --------------------------------------------

    def _write_rna_reads(
        self,
        t_idx: int,
        len_counts: dict[int, int],
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
        *,
        is_nrna: bool,
    ) -> int:
        """Generate and write all mRNA or nRNA reads for one transcript.

        When ``is_nrna=False``, uses the pre-extracted spliced mRNA sequence.
        When ``is_nrna=True``, uses the pre-mRNA (unspliced) sequence.
        """
        t = self.transcripts[t_idx]
        rng = self._rng
        ss = self.strand_specificity
        strand_char = self._t_strand_chars[t_idx]
        t_id = t.t_id

        if is_nrna:
            seq_bytes = self._get_premrna_bytes(t_idx)
            seq_len = int(self._premrna_lengths[t_idx])
            qname_prefix = f"nrna_{t_id}"
        else:
            seq_bytes = self._t_seq_bytes[t_idx]
            seq_len = int(self._t_lengths[t_idx])
            qname_prefix = t_id

        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh else None
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            read_len = min(self.sim_params.read_length, frag_len)
            eff_len = seq_len - frag_len + 1
            if eff_len <= 0:
                continue

            frag_starts = rng.integers(0, eff_len, size=count)
            flip_mask = rng.random(count) >= ss if ss < 1.0 else None
            quals = self._get_quals(read_len)

            # Vectorized read extraction
            r2_seqs, r1_seqs = _batch_extract_reads(seq_bytes, frag_starts, frag_len, read_len)

            # Strand flipping
            if flip_mask is not None:
                flip_idx = np.where(flip_mask)[0]
                if len(flip_idx) > 0:
                    r1_seqs[flip_idx], r2_seqs[flip_idx] = (
                        r2_seqs[flip_idx].copy(),
                        r1_seqs[flip_idx].copy(),
                    )

            # Batch error introduction
            if self._has_errors:
                self._introduce_errors_batch(r1_seqs)
                self._introduce_errors_batch(r2_seqs)

            # Build query names
            base_n = n_written
            qnames = [
                f"{qname_prefix}:{int(frag_starts[i])}-{int(frag_starts[i]) + frag_len}:{strand_char}:{base_n + i}"
                for i in range(count)
            ]

            # Batch FASTQ writes
            r1_buf.write_records_batch(qnames, r1_seqs, quals, "/1")
            r2_buf.write_records_batch(qnames, r2_seqs, quals, "/2")

            # BAM records (pysam API requires per-record construction)
            if bam_fh is not None and ref_id is not None:
                if is_nrna:
                    self._write_nrna_bam_batch(
                        bam_fh, qnames, r1_seqs, r2_seqs, t,
                        frag_starts, frag_len, read_len, flip_mask, ref_id,
                    )
                else:
                    self._write_mrna_bam_batch(
                        bam_fh, qnames, r1_seqs, r2_seqs, t,
                        frag_starts, frag_len, read_len, flip_mask, ref_id,
                    )

            n_written += count

        return n_written

    # -- BAM batch writers --------------------------------------------------

    def _write_mrna_bam_batch(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        r1_seqs: np.ndarray,
        r2_seqs: np.ndarray,
        t: Transcript,
        frag_starts: np.ndarray,
        frag_len: int,
        read_len: int,
        flip_mask: np.ndarray | None,
        ref_id: int,
    ) -> None:
        for i in range(len(qnames)):
            frag_start = int(frag_starts[i])
            frag_end = frag_start + frag_len
            flipped = flip_mask is not None and flip_mask[i]

            r2_t_end = min(frag_start + read_len, frag_end)
            r2_blocks = _transcript_to_genomic_blocks(frag_start, r2_t_end, t)
            r1_t_start = max(frag_end - read_len, frag_start)
            r1_blocks = _transcript_to_genomic_blocks(r1_t_start, frag_end, t)
            if not r2_blocks or not r1_blocks:
                continue

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

            tags = self._nh1_tags
            if len(r2_blocks) > 1 or len(r1_blocks) > 1:
                xs_strand = "-" if t.strand == Strand.NEG else "+"
                tags = [("NH", 1), ("XS", xs_strand)]

            r1_flag = _BASE_R1_FLAG
            r2_flag = _BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= _FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= _FLAG_MATE_REVERSE
                r2_flag |= _FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= _FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_start <= r2_start else -tlen
            r2_tlen = -r1_tlen

            r1_seq_str = r1_seqs[i].tobytes().decode("ascii")
            r2_seq_str = r2_seqs[i].tobytes().decode("ascii")

            bam_fh.write(_make_bam_record(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_start, r1_cigar, ref_id, r2_start, r1_tlen, tags=tags,
            ))
            bam_fh.write(_make_bam_record(
                self._bam_header, qnames[i], r2_seq_str, r2_flag, ref_id,
                r2_start, r2_cigar, ref_id, r1_start, r2_tlen, tags=tags,
            ))

    def _write_nrna_bam_batch(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        r1_seqs: np.ndarray,
        r2_seqs: np.ndarray,
        t: Transcript,
        frag_starts: np.ndarray,
        frag_len: int,
        read_len: int,
        flip_mask: np.ndarray | None,
        ref_id: int,
    ) -> None:
        for i in range(len(qnames)):
            frag_start = int(frag_starts[i])
            frag_end = frag_start + frag_len
            flipped = flip_mask is not None and flip_mask[i]

            g_start, g_end = _premrna_to_genomic_interval(frag_start, frag_end, t)
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

            r1_flag = _BASE_R1_FLAG
            r2_flag = _BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= _FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= _FLAG_MATE_REVERSE
                r2_flag |= _FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= _FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_g_start <= r2_g_start else -tlen
            r2_tlen = -r1_tlen

            r1_read_len = r1_g_end - r1_g_start
            r2_read_len = r2_g_end - r2_g_start

            r1_seq_str = r1_seqs[i].tobytes().decode("ascii")
            r2_seq_str = r2_seqs[i].tobytes().decode("ascii")

            bam_fh.write(_make_bam_record(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_g_start, [(pysam.CMATCH, r1_read_len)],
                ref_id, r2_g_start, r1_tlen, tags=self._nh1_tags,
            ))
            bam_fh.write(_make_bam_record(
                self._bam_header, qnames[i], r2_seq_str, r2_flag, ref_id,
                r2_g_start, [(pysam.CMATCH, r2_read_len)],
                ref_id, r1_g_start, r2_tlen, tags=self._nh1_tags,
            ))

    # -- gDNA reads (vectorized per chromosome) -----------------------------

    def _accumulate_gdna_counts(self, n_gdna: int) -> dict[tuple[int, int], int]:
        """Sample gDNA fragment lengths + chromosomes.

        Returns dict[(ref_idx, frag_len)] = count.
        """
        if n_gdna == 0 or not self._gdna_refs:
            return {}

        rng = self._rng
        frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)
        gdna_ref_lengths_arr = np.array(self._gdna_ref_lengths, dtype=np.float64)

        counts: dict[tuple[int, int], int] = {}
        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            chrom_eff = np.maximum(0, gdna_ref_lengths_arr - fl + 1)
            total_eff = chrom_eff.sum()
            if total_eff <= 0:
                continue
            chrom_probs = chrom_eff / total_eff
            chrom_indices = rng.choice(len(self._gdna_refs), size=fc, p=chrom_probs)
            unique_chroms, chrom_counts = np.unique(chrom_indices, return_counts=True)
            for ci, cc in zip(unique_chroms, chrom_counts):
                counts[(int(ci), fl)] = int(cc)
        return counts

    def _write_gdna_chunk(
        self,
        ref_idx: int,
        fl: int,
        count: int,
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
        n_offset: int = 0,
    ) -> int:
        """Generate ``count`` gDNA fragments at chromosome ref_idx with frag length fl.

        Returns number of fragments written.
        """
        if count <= 0:
            return 0
        rng = self._rng
        ref = self._gdna_refs[ref_idx]
        chrom_len = self._gdna_ref_lengths[ref_idx]
        ref_id = self._ref_name_to_id.get(ref)
        chrom_bytes = self._genome.get(ref)
        eff_len = chrom_len - fl + 1
        if eff_len <= 0:
            return 0

        read_len = min(self.sim_params.read_length, fl)
        quals = self._get_quals(read_len)

        starts = rng.integers(0, eff_len, size=count)
        chrom_strands = rng.integers(0, 2, size=count)

        r2_seqs, r1_seqs = _batch_extract_reads(chrom_bytes, starts, fl, read_len)

        # Negative strand: swap R1/R2 (fragment = revcomp of genomic)
        neg_mask = chrom_strands.astype(bool)
        if neg_mask.any():
            neg_idx = np.where(neg_mask)[0]
            r1_seqs[neg_idx], r2_seqs[neg_idx] = (
                r2_seqs[neg_idx].copy(),
                r1_seqs[neg_idx].copy(),
            )

        if self._has_errors:
            self._introduce_errors_batch(r1_seqs)
            self._introduce_errors_batch(r2_seqs)

        qnames = [
            f"gdna:{ref}:{int(starts[j])}-{int(starts[j]) + fl}:"
            f"{'r' if chrom_strands[j] else 'f'}:{n_offset + j}"
            for j in range(count)
        ]

        r1_buf.write_records_batch(qnames, r1_seqs, quals, "/1")
        r2_buf.write_records_batch(qnames, r2_seqs, quals, "/2")

        if bam_fh is not None and ref_id is not None:
            for j in range(count):
                start = int(starts[j])
                end = start + fl
                is_neg = bool(chrom_strands[j])

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
                r1_tlen = tlen if r1_start_pos <= r2_start_pos else -tlen
                r2_tlen = -r1_tlen

                r1_seq_str = r1_seqs[j].tobytes().decode("ascii")
                r2_seq_str = r2_seqs[j].tobytes().decode("ascii")

                bam_fh.write(_make_bam_record(
                    self._bam_header, qnames[j], r1_seq_str,
                    r1_flag, ref_id, r1_start_pos,
                    [(pysam.CMATCH, read_len)],
                    ref_id, r2_start_pos, r1_tlen,
                    tags=self._nh1_tags,
                ))
                bam_fh.write(_make_bam_record(
                    self._bam_header, qnames[j], r2_seq_str,
                    r2_flag, ref_id, r2_start_pos,
                    [(pysam.CMATCH, read_len)],
                    ref_id, r1_start_pos, r2_tlen,
                    tags=self._nh1_tags,
                ))
        return count

    def _write_gdna_from_counts(
        self,
        gdna_counts: dict[tuple[int, int], int],
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
    ) -> int:
        """Write gDNA reads from a pre-sampled (ref_idx, fl) -> count map."""
        n_written = 0
        for (ref_idx, fl), count in gdna_counts.items():
            n_written += self._write_gdna_chunk(
                ref_idx, fl, count, r1_buf, r2_buf, bam_fh, n_offset=n_written,
            )
        return n_written

    # -- Main entry point ---------------------------------------------------

    def simulate_and_write(
        self,
        output_dir: Path,
        n_rna: int,
        n_gdna: int = 0,
        *,
        n_mrna: int | None = None,
        n_nrna: int | None = None,
        oracle_bam: bool = True,
        prefix: str = "sim",
        n_workers: int = 1,
    ) -> tuple[Path, Path, Path | None]:
        """Single-pass simulation: accumulate counts, generate, write.

        When *n_mrna* and *n_nrna* are provided, the mRNA and nRNA
        fragment pools are sampled independently (precise fragment-count
        control).  *n_rna* is ignored in this mode but still used for
        logging.  Otherwise *n_rna* fragments are drawn from the
        combined mRNA + nRNA pool (original behaviour).

        When *n_workers* > 1 the per-transcript and per-(chrom, frag-len)
        gDNA work is sharded across worker processes (fork-based) and the
        resulting FASTQ.gz / BAM shard files are concatenated in the
        parent.

        Returns (r1_path, r2_path, bam_path | None).
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        r1_path = output_dir / f"{prefix}_R1.fq.gz"
        r2_path = output_dir / f"{prefix}_R2.fq.gz"
        bam_path = output_dir / f"{prefix}_oracle.bam" if oracle_bam else None

        # Phase 1: accumulate per-transcript / per-chunk counts (main process).
        # Uses self._rng so determinism follows the configured seed.
        logger.info("Accumulating RNA fragment counts...")
        mrna_counts, nrna_counts = self._accumulate_rna_counts(
            n_rna, n_mrna=n_mrna, n_nrna=n_nrna,
        )
        gdna_counts = self._accumulate_gdna_counts(n_gdna)

        total_mrna = sum(sum(d.values()) for d in mrna_counts.values())
        total_nrna = sum(sum(d.values()) for d in nrna_counts.values())
        total_gdna = sum(gdna_counts.values())
        logger.info(
            "Fragment counts: %d mRNA (%d txs), %d nRNA (%d txs), %d gDNA chunks",
            total_mrna, len(mrna_counts), total_nrna, len(nrna_counts), len(gdna_counts),
        )

        if n_workers <= 1:
            self._write_all_serial(
                r1_path, r2_path, bam_path, mrna_counts, nrna_counts, gdna_counts,
                oracle_bam,
            )
        else:
            self._write_all_parallel(
                output_dir, prefix, r1_path, r2_path, bam_path,
                mrna_counts, nrna_counts, gdna_counts,
                oracle_bam, n_workers,
            )

        n_written = total_mrna + total_nrna + total_gdna
        logger.info(
            "Wrote %d read pairs -> %s (RNA=%d, gDNA=%d, oracle=%s, workers=%d)",
            n_written, output_dir, n_rna, n_gdna, oracle_bam, n_workers,
        )
        return r1_path, r2_path, bam_path

    # -- Serial / parallel write paths --------------------------------------

    def _write_all_serial(
        self,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_counts: dict[int, dict[int, int]],
        nrna_counts: dict[int, dict[int, int]],
        gdna_counts: dict[tuple[int, int], int],
        oracle_bam: bool,
    ) -> None:
        bam_fh: pysam.AlignmentFile | None = None
        with (
            _open_gzip_write(r1_path) as r1_fh,
            _open_gzip_write(r2_path) as r2_fh,
        ):
            r1_buf = _FastqBuffer(r1_fh)
            r2_buf = _FastqBuffer(r2_fh)
            if oracle_bam:
                bam_fh = pysam.AlignmentFile(
                    str(bam_path), "wb", header=self._bam_header,
                )
            try:
                for t_idx in sorted(mrna_counts):
                    self._write_rna_reads(
                        t_idx, mrna_counts[t_idx], r1_buf, r2_buf, bam_fh,
                        is_nrna=False,
                    )
                for t_idx in sorted(nrna_counts):
                    self._write_rna_reads(
                        t_idx, nrna_counts[t_idx], r1_buf, r2_buf, bam_fh,
                        is_nrna=True,
                    )
                self._write_gdna_from_counts(gdna_counts, r1_buf, r2_buf, bam_fh)
                r1_buf.close()
                r2_buf.close()
            finally:
                if bam_fh is not None:
                    bam_fh.close()

    def _write_all_parallel(
        self,
        output_dir: Path,
        prefix: str,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_counts: dict[int, dict[int, int]],
        nrna_counts: dict[int, dict[int, int]],
        gdna_counts: dict[tuple[int, int], int],
        oracle_bam: bool,
        n_workers: int,
    ) -> None:
        # Build per-shard task lists via greedy LPT bin-packing on fragment counts.
        mrna_items = [(t, dict(d)) for t, d in mrna_counts.items()]
        nrna_items = [(t, dict(d)) for t, d in nrna_counts.items()]
        gdna_items = list(gdna_counts.items())  # [((ref_idx, fl), count), ...]

        mrna_shards = _shard_by_count(
            mrna_items, n_workers, weight=lambda x: sum(x[1].values()),
        )
        nrna_shards = _shard_by_count(
            nrna_items, n_workers, weight=lambda x: sum(x[1].values()),
        )
        gdna_shards = _shard_by_count(
            gdna_items, n_workers, weight=lambda x: x[1],
        )

        shard_dir = output_dir / f".{prefix}_shards"
        shard_dir.mkdir(parents=True, exist_ok=True)

        base_seed = int(self._rng.integers(0, 2**31 - 1))

        tasks = []
        for k in range(n_workers):
            tasks.append((
                k,
                mrna_shards[k],
                nrna_shards[k],
                gdna_shards[k],
                str(shard_dir),
                prefix,
                oracle_bam,
                base_seed + k,
            ))

        # Bind self into a module global so fork'd children inherit it without
        # pickling.  Passing self via Pool.map would fail (pysam handles).
        global _WORKER_SIM
        _WORKER_SIM = self

        ctx = multiprocessing.get_context("fork")
        results: list[tuple[str, str, str | None]]
        try:
            with ctx.Pool(n_workers) as pool:
                results = pool.map(_shard_task, tasks)
        finally:
            _WORKER_SIM = None

        shard_r1 = [Path(r[0]) for r in results]
        shard_r2 = [Path(r[1]) for r in results]
        shard_bams = [Path(r[2]) for r in results if r[2] is not None]

        logger.info("Concatenating %d FASTQ shards...", len(shard_r1))
        _concat_files_binary(shard_r1, r1_path)
        _concat_files_binary(shard_r2, r2_path)

        if oracle_bam and shard_bams:
            logger.info("Concatenating %d BAM shards...", len(shard_bams))
            pysam.cat("-o", str(bam_path), *[str(p) for p in shard_bams])

        # Cleanup shard files
        for p in shard_r1 + shard_r2 + shard_bams:
            try:
                p.unlink()
            except OSError:
                pass
        try:
            shard_dir.rmdir()
        except OSError:
            pass

    def _run_shard(
        self,
        shard_id: int,
        mrna_items: list[tuple[int, dict[int, int]]],
        nrna_items: list[tuple[int, dict[int, int]]],
        gdna_items: list[tuple[tuple[int, int], int]],
        shard_dir: str,
        prefix: str,
        oracle_bam: bool,
        seed: int,
    ) -> tuple[str, str, str | None]:
        """Worker entry: write one shard's FASTQ + BAM files."""
        # Independent RNG per shard
        self._rng = np.random.default_rng(seed)

        shard_path = Path(shard_dir)
        r1_path = shard_path / f"{prefix}.shard{shard_id:03d}.R1.fq.gz"
        r2_path = shard_path / f"{prefix}.shard{shard_id:03d}.R2.fq.gz"
        bam_path = shard_path / f"{prefix}.shard{shard_id:03d}.bam" if oracle_bam else None

        bam_fh: pysam.AlignmentFile | None = None
        with (
            _open_gzip_write(r1_path, threads=1) as r1_fh,
            _open_gzip_write(r2_path, threads=1) as r2_fh,
        ):
            r1_buf = _FastqBuffer(r1_fh)
            r2_buf = _FastqBuffer(r2_fh)
            if oracle_bam:
                bam_fh = pysam.AlignmentFile(
                    str(bam_path), "wb", header=self._bam_header,
                )
            try:
                for t_idx, len_counts in mrna_items:
                    self._write_rna_reads(
                        t_idx, len_counts, r1_buf, r2_buf, bam_fh, is_nrna=False,
                    )
                for t_idx, len_counts in nrna_items:
                    self._write_rna_reads(
                        t_idx, len_counts, r1_buf, r2_buf, bam_fh, is_nrna=True,
                    )
                n_offset = 0
                for (ref_idx, fl), count in gdna_items:
                    n_offset += self._write_gdna_chunk(
                        ref_idx, fl, count, r1_buf, r2_buf, bam_fh, n_offset=n_offset,
                    )
                r1_buf.close()
                r2_buf.close()
            finally:
                if bam_fh is not None:
                    bam_fh.close()

        return (str(r1_path), str(r2_path), str(bam_path) if bam_path else None)

    def close(self) -> None:
        """Release resources."""
        self._genome.clear()
        self.fasta.close()


# ═══════════════════════════════════════════════════════════════════
# Parallel-shard helpers (module-level for fork-based multiprocessing)
# ═══════════════════════════════════════════════════════════════════


# Set in the parent process before Pool creation; inherited by fork()
# children so the pre-extracted sequence cache is shared via copy-on-write
# without pickling pysam handles.
_WORKER_SIM: "WholeGenomeSimulator | None" = None


def _shard_task(args):
    """Pool worker entry. Dispatches to the inherited simulator's _run_shard."""
    if _WORKER_SIM is None:
        raise RuntimeError("Worker simulator not initialized (fork inheritance failed)")
    return _WORKER_SIM._run_shard(*args)


def _shard_by_count(items, n_shards, *, weight):
    """Greedy LPT bin-packing: distribute items across shards by descending weight."""
    if n_shards <= 1:
        return [list(items)]
    sorted_items = sorted(items, key=lambda x: -weight(x))
    shards: list[list] = [[] for _ in range(n_shards)]
    loads = [0] * n_shards
    for item in sorted_items:
        idx = min(range(n_shards), key=lambda i: loads[i])
        shards[idx].append(item)
        loads[idx] += weight(item)
    return shards


def _concat_files_binary(srcs: list[Path], dst: Path) -> None:
    """Concatenate files byte-wise. Valid for gzip streams (per RFC 1952)."""
    with open(dst, "wb") as out:
        for s in srcs:
            with open(s, "rb") as f:
                shutil.copyfileobj(f, out, length=4 * 1024 * 1024)


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
        "gdna": asdict(cfg.gdna),
        "nrna": asdict(cfg.nrna),
        "strand_specificities": cfg.strand_specificities,
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

    1. Load genome + GTF -> transcripts
    2. Assign base abundances (total RNA) once
    3. Sweep nrna_fracs x gdna_rates x strand_specificities
    4. Write manifest

    Returns list of condition dicts (for manifest).
    """
    import copy
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
                        n_workers=sim.n_workers,
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
    p.add_argument("-j", "--n-workers", type=int, default=None,
                   help="Worker processes for parallel read generation (overrides YAML)")
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
    if args.n_workers is not None:
        cfg.simulation.n_workers = args.n_workers
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
    print(f"  Workers:          {cfg.simulation.n_workers}", flush=True)
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
