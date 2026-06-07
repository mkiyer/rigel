#!/usr/bin/env python3
"""Whole-genome RNA-seq read simulator.

Generates paired-end RNA-seq reads across the full transcriptome with
configurable mRNA/nRNA abundances, gDNA contamination, and strand
specificity.  Outputs FASTQ files and optional oracle BAM files.

Abundance model
---------------
Two abundance sources:

**random** — For each transcript, sample whether it is expressed
(Bernoulli with probability ``frac_expressed``).  For expressed
transcripts, draw total RNA abundance from a **log-uniform**
distribution:

    log_total ~ Uniform(log(min), log(max))
    total = exp(log_total)

Nascent RNA is controlled separately by the ``nrna`` section. The canonical
benchmark mode is additive:

    mRNA = base abundance
    nRNA = base abundance × nrna_ratio

Single-exon transcripts always get nRNA = 0.

**file** — Load from a TSV with columns ``transcript_id``,
``mrna_abundance``, ``nrna_abundance``.

Fragment allocation
-------------------
In additive-ratio mode, ``n_rna_fragments`` is the mature RNA depth. nRNA and
gDNA fragments are added on top:

    n_mrna = n_rna_fragments
    n_nrna = round(n_mrna × nrna_ratio)
    n_gdna = round(n_mrna × gdna_rate)
    n_total = n_mrna + n_nrna + n_gdna

Condition grid
--------------
Sweeps: ``nrna ratios × gdna_rates × strand_specificities``.

When the abundance file provides explicit nRNA data, the nRNA sweep is
skipped (single condition using the file's nRNA values).

Usage
-----
::

    python scripts/sim/simulate_reads.py --config scripts/sim/configs/sim_example.yaml

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
import logging
import multiprocessing
import shutil
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
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
from rigel.sim.bam import (
    BASE_R1_FLAG,
    BASE_R2_FLAG,
    FLAG_MATE_REVERSE,
    FLAG_REVERSE,
    blocks_to_cigar,
    make_aligned_segment,
    premrna_to_genomic_interval,
    transcript_to_genomic_blocks,
)
from rigel.sim.capture import CaptureConfig, CaptureScenario, CaptureSampler
from rigel.sim.genome import reverse_complement
from rigel.sim.sampling import truncated_normal_frag_lengths
from rigel.sim.manifest import (
    gdna_label_for_rate,
    write_manifest as write_manifest_file,
)

logger = logging.getLogger(__name__)

# Byte-level complement table for vectorized reverse-complement
_BYTE_COMPLEMENT = np.zeros(256, dtype=np.uint8)
for _c, _rc in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _BYTE_COMPLEMENT[_c] = _rc

_FASTQ_BUFFER_SIZE = 100_000


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


def _bam_seq_from_fastq_bytes(seq: np.ndarray, is_reverse: bool) -> str:
    """Return the SAM/BAM stored SEQ for a FASTQ-oriented read sequence."""
    if is_reverse:
        seq = _BYTE_COMPLEMENT[seq[::-1]]
    return seq.tobytes().decode("ascii")


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
    #: gDNA strand intra-class correlation in [0, 1) — region-to-region overdispersion of the
    #: sense/antisense split around ½. 0 ⇒ exact 50/50 (Binomial). See GDNAConfig in reads.py.
    #: This is the *per-condition* value the simulator reads; the suite sweep axis is below.
    strand_overdispersion: float = 0.0
    #: Suite sweep axis over gDNA strand overdispersion (one condition per value). ``None`` ⇒ a
    #: single condition at ``strand_overdispersion`` (backward-compatible: no name change).
    strand_overdispersions: list[float] | None = None
    strand_overdispersion_labels: list[str] | None = None


@dataclass
class NRNAConfig:
    """Nascent RNA spike-in sweep configuration.

    The only generated sweep mode is ``additive_ratio``. Each entry in
    ``ratios`` adds nascent RNA independently of mature RNA:

        nrna_abundance = mrna_abundance * ratio

    When abundances come from a file that already contains explicit nRNA
    data, the configured sweep is ignored and the file's nRNA values are
    used as-is in a single nRNA condition.
    """

    mode: str = "additive_ratio"
    ratios: list[float] = field(default_factory=lambda: [0.0])
    ratio_ranges: list[tuple[float, float]] | None = None
    ratio_labels: list[str] | None = None
    eligible_fraction: float = 1.0
    seed: int = 42


@dataclass
class WholeGenomeSimConfig:
    """Top-level simulation configuration."""

    genome: str = ""
    gtf: str = ""
    outdir: str = "sim_output"
    transcript_filter: str = "all"  # "all", "basic", "mane", "ccds"

    simulation: SimulationParams = field(default_factory=SimulationParams)
    abundance: AbundanceConfig = field(default_factory=AbundanceConfig)
    gdna: GDNASimConfig = field(default_factory=GDNASimConfig)
    nrna: NRNAConfig = field(default_factory=NRNAConfig)
    capture: CaptureConfig = field(default_factory=CaptureConfig)
    capture_configs: list[CaptureScenario] = field(default_factory=list)
    strand_specificities: list[float] = field(
        default_factory=lambda: [1.0]
    )

    oracle_bam: bool = True
    verbose: bool = True


def _capture_config_from_mapping(
    raw: dict,
    defaults: dict | None = None,
    *,
    label: str = "capture",
    require_probes_when_enabled: bool = False,
) -> CaptureConfig:
    """Build a CaptureConfig from a YAML mapping plus optional defaults."""
    merged = dict(defaults or {})
    merged.update(raw)
    enabled = bool(merged.get("enabled", True))
    if not enabled:
        return CaptureConfig()

    probes = merged.get("probes", None)
    if not probes:
        if require_probes_when_enabled and raw.get("enabled") is True:
            raise ValueError(f"capture config '{label}' is enabled but has no probes")
        return CaptureConfig()

    return CaptureConfig(
        probes=str(probes),
        probe_format=str(merged.get("probe_format", merged.get("format", "auto"))),
        off_target_weight=float(merged.get("off_target_weight", 1.0)),
        binding_per_base=float(merged.get("binding_per_base", 10.0)),
        gdna_split_penalty=float(merged.get("gdna_split_penalty", 0.2)),
        min_overlap=int(merged.get("min_overlap", 1)),
    )


def _capture_label_text(value: object) -> str:
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value)


def _capture_scenarios_from_mapping(raw: dict) -> list[CaptureScenario]:
    """Parse capture.configs/capture.scenarios from a YAML mapping."""
    raw_configs = raw.get("configs", raw.get("scenarios"))
    if raw_configs is None:
        return []
    if not isinstance(raw_configs, list) or not raw_configs:
        raise ValueError("capture.configs must be a non-empty list")

    defaults = {
        key: value
        for key, value in raw.items()
        if key not in {"configs", "scenarios"}
    }
    scenarios: list[CaptureScenario] = []
    seen_labels: set[str] = set()
    for index, item in enumerate(raw_configs, start=1):
        if not isinstance(item, dict):
            raise ValueError("each capture.configs entry must be a mapping")
        label = _capture_label_text(item.get("label", item.get("name", f"c{index}")))
        if not label:
            raise ValueError("capture config labels must be non-empty")
        if label in seen_labels:
            raise ValueError(f"duplicate capture config label: {label}")
        seen_labels.add(label)
        scenarios.append(CaptureScenario(
            label=label,
            config=_capture_config_from_mapping(
                item,
                defaults,
                label=label,
                require_probes_when_enabled=True,
            ),
        ))
    return scenarios


def _capture_grid(cfg: WholeGenomeSimConfig) -> tuple[list[CaptureScenario], bool]:
    """Return capture scenarios and whether condition names should include labels."""
    if cfg.capture_configs:
        return cfg.capture_configs, True
    label = "on" if cfg.capture.probes else "off"
    return [CaptureScenario(label=label, config=cfg.capture)], False


def parse_yaml_config(path: str | Path) -> WholeGenomeSimConfig:
    """Parse a YAML config file into a WholeGenomeSimConfig."""
    if yaml is None:
        raise ImportError("PyYAML is required: pip install pyyaml")

    with open(path) as f:
        raw = yaml.safe_load(f)

    cfg = WholeGenomeSimConfig()
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
    nrna.mode = str(nrna_raw.get("mode", "additive_ratio"))
    if nrna.mode not in {"additive_ratio", "random_fraction"}:
        raise ValueError("nrna.mode must be 'additive_ratio' or 'random_fraction'")
    if "fracs" in nrna_raw or "frac_labels" in nrna_raw:
        raise ValueError("nrna.fracs is no longer supported; use nrna.ratios")
    raw_ratios = nrna_raw.get("ratios", None)
    if raw_ratios is not None:
        nrna.ratios = [float(r) for r in raw_ratios]
    raw_ratio_ranges = nrna_raw.get("ratio_ranges", None)
    if raw_ratio_ranges is not None:
        nrna.ratio_ranges = [
            (float(pair[0]), float(pair[1])) for pair in raw_ratio_ranges
        ]
    nrna.ratio_labels = nrna_raw.get("ratio_labels", None)
    nrna.eligible_fraction = float(nrna_raw.get("eligible_fraction", 1.0))
    nrna.seed = int(nrna_raw.get("seed", 42))
    if not 0.0 <= nrna.eligible_fraction <= 1.0:
        raise ValueError("nrna.eligible_fraction must be between 0 and 1")
    if nrna.mode == "additive_ratio":
        expected_len = len(nrna.ratios)
    else:
        if nrna.ratio_ranges is None:
            raise ValueError("nrna.ratio_ranges is required for mode='random_fraction'")
        for lo, hi in nrna.ratio_ranges:
            if lo < 0 or hi < 0 or hi < lo:
                raise ValueError("nrna.ratio_ranges entries must satisfy 0 <= min <= max")
        expected_len = len(nrna.ratio_ranges)
    if nrna.ratio_labels is not None and len(nrna.ratio_labels) != expected_len:
        raise ValueError("nrna.ratio_labels must match the number of nRNA scenarios")

    # gDNA
    gd_raw = raw.get("gdna", {})
    gd = cfg.gdna
    gd.rates = [float(r) for r in gd_raw.get("rates", [0.0])]
    gd.rate_labels = gd_raw.get("rate_labels", None)
    gd.frag_mean = float(gd_raw.get("frag_mean", 350.0))
    gd.frag_std = float(gd_raw.get("frag_std", 100.0))
    gd.frag_min = int(gd_raw.get("frag_min", 100))
    gd.frag_max = int(gd_raw.get("frag_max", 1000))
    gd.strand_overdispersion = float(gd_raw.get("strand_overdispersion", 0.0))
    _ods = gd_raw.get("strand_overdispersions", None)
    gd.strand_overdispersions = [float(o) for o in _ods] if _ods is not None else None
    gd.strand_overdispersion_labels = gd_raw.get("strand_overdispersion_labels", None)

    # Hybrid capture
    cap_raw = raw.get("capture", {}) or {}
    if not isinstance(cap_raw, dict):
        raise ValueError("capture must be a YAML mapping")
    cfg.capture_configs = _capture_scenarios_from_mapping(cap_raw)
    if cfg.capture_configs:
        cfg.capture = cfg.capture_configs[0].config
    else:
        cfg.capture = _capture_config_from_mapping(cap_raw)

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
    is left at zero. The additive nRNA ratio is applied per-condition in
    ``run_simulation`` via ``apply_nrna_ratio``.

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
) -> tuple[str, object]:
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


def apply_nrna_ratio(
    transcripts: list[Transcript],
    ratio: float,
) -> None:
    """Set nRNA abundance as an additive ratio of mature RNA abundance."""
    n_spiked = 0
    n_single = 0

    for t in transcripts:
        mrna = t.abundance or 0.0
        if mrna <= 0:
            t.nrna_abundance = 0.0
            continue
        if len(t.exons) <= 1:
            t.nrna_abundance = 0.0
            n_single += 1
            continue
        t.nrna_abundance = mrna * ratio
        if ratio > 0:
            n_spiked += 1

    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    total_nrna = sum(t.nrna_abundance for t in transcripts)
    logger.info(
        "Set additive nRNA ratio: %.3g (%d transcripts, mRNA=%.1f, nRNA=%.1f, "
        "%d single-exon zeroed)",
        ratio, n_spiked, total_mrna, total_nrna, n_single,
    )


def apply_random_nrna_fraction(
    transcripts: list[Transcript],
    ratio_range: tuple[float, float],
    *,
    eligible_fraction: float,
    seed: int,
) -> float:
    """Assign nRNA to a random subset of expressed multi-exon transcripts.

    Returns the realized total nRNA:mRNA molecular abundance ratio.
    """
    lo, hi = ratio_range
    if lo < 0 or hi < lo:
        raise ValueError("ratio_range must satisfy 0 <= min <= max")
    if not 0.0 <= eligible_fraction <= 1.0:
        raise ValueError("eligible_fraction must be between 0 and 1")

    rng = np.random.default_rng(seed)
    n_eligible = 0
    n_spiked = 0
    n_single = 0

    for t in transcripts:
        mrna = t.abundance or 0.0
        t.nrna_abundance = 0.0
        if mrna <= 0:
            continue
        if len(t.exons) <= 1:
            n_single += 1
            continue
        n_eligible += 1
        if rng.random() >= eligible_fraction:
            continue
        ratio = float(rng.uniform(lo, hi)) if hi > lo else lo
        t.nrna_abundance = mrna * ratio
        if ratio > 0:
            n_spiked += 1

    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    total_nrna = sum(t.nrna_abundance for t in transcripts)
    realized_ratio = total_nrna / total_mrna if total_mrna > 0 else 0.0
    logger.info(
        "Set random nRNA fractions: range=[%.3g, %.3g], eligible_fraction=%.3g, "
        "spiked=%d/%d expressed multi-exon, mRNA=%.1f, nRNA=%.1f, "
        "realized_ratio=%.4g, %d single-exon zeroed",
        lo, hi, eligible_fraction, n_spiked, n_eligible,
        total_mrna, total_nrna, realized_ratio, n_single,
    )
    return realized_ratio


def total_nrna_to_mrna_ratio(transcripts: list[Transcript]) -> float:
    """Return total nRNA abundance divided by total mature RNA abundance."""
    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    if total_mrna <= 0:
        return 0.0
    return sum(t.nrna_abundance for t in transcripts) / total_mrna


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
    does not supply nRNA, the caller applies the configured nRNA sweep.
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

    def __init__(
        self,
        fasta_path: str | Path,
        transcripts: list[Transcript],
        sim_params: SimulationParams,
        gdna_config: GDNASimConfig,
        *,
        strand_specificity: float = 1.0,
        seed: int | None = None,
        capture_config: CaptureConfig | None = None,
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

        # gDNA strand overdispersion: per-ref exon-derived regions each with a shared +strand
        # rate ~ Beta(a, a), so a region's gDNA strand split is Beta-Binomial(½, overdispersion)
        # — matching the calibrator's seed-region partition. 0 ⇒ uniform 50/50 (no regions built).
        self._gdna_strand_regions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        gdna_od = float(getattr(self.gdna_config, "strand_overdispersion", 0.0) or 0.0)
        if gdna_od > 0.0 and self._gdna_refs:
            self._init_gdna_strand_regions(gdna_od)

        # BAM header and reference mapping (built once, reused)
        self._ref_names = list(self.fasta.references)
        self._ref_lengths = [self.fasta.get_reference_length(r) for r in self._ref_names]
        self._ref_name_to_id: dict[str, int] = {
            name: i for i, name in enumerate(self._ref_names)
        }
        self.capture = CaptureSampler.from_config(
            capture_config,
            transcripts,
            dict(zip(self._ref_names, self._ref_lengths)),
        )
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

    def _sample_rna_frag_lengths(self, n: int) -> np.ndarray:
        p = self.sim_params
        return truncated_normal_frag_lengths(
            self._rng, n, p.frag_mean, p.frag_std, p.frag_min, self._frag_max
        )

    def _init_gdna_strand_regions(self, overdispersion: float) -> None:
        """Build per-gDNA-ref exon-derived regions, each with a shared +strand rate ~ Beta(a, a).

        ``Beta(a, a)`` has intra-class correlation ``1/(2a + 1)``, so ``a = ½(1 − od)/od`` gives a
        region-to-region strand overdispersion of ``od``. All gDNA fragments in a region share its
        rate, so the per-region sense/antisense count is Beta-Binomial(½, od) — the distribution
        the calibrator's gDNA strand model fits.
        """
        alpha = 0.5 * (1.0 - overdispersion) / overdispersion
        exons_by_ref: dict[str, set[int]] = defaultdict(set)
        for t in self.transcripts:
            if t.ref is None:
                continue
            for e in t.exons:
                exons_by_ref[t.ref].add(int(e.start))
                exons_by_ref[t.ref].add(int(e.end))
        for ref, length in zip(self._gdna_refs, self._gdna_ref_lengths):
            bounds = {b for b in exons_by_ref.get(ref, set()) if 0 < b < length}
            boundaries = np.array([0, *sorted(bounds), int(length)], dtype=np.int64)
            n_regions = max(len(boundaries) - 1, 1)
            p_plus = self._rng.beta(alpha, alpha, size=n_regions)
            self._gdna_strand_regions[ref] = (boundaries, p_plus)

    def _sample_gdna_frag_lengths(self, n: int) -> np.ndarray:
        g = self.gdna_config
        return truncated_normal_frag_lengths(
            self._rng, n, g.frag_mean, g.frag_std, g.frag_min, g.frag_max
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
        *,
        space: str,
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
            eff = self.capture.partition_array(space, range(len(abundances)), lengths, fl)
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
                n_mrna, self._mrna_abund, self._t_lengths, space="mrna",
            )
            nrna_counts = self._accumulate_pool(
                n_nrna, self._nrna_abund, self._premrna_lengths, space="nrna",
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

            mrna_eff = self.capture.partition_array("mrna", range(N), self._t_lengths, fl)
            nrna_eff = self.capture.partition_array("nrna", range(N), self._premrna_lengths, fl)
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

            capture_space = "nrna" if is_nrna else "mrna"
            frag_starts = self.capture.sample_starts(
                capture_space, t_idx, seq_len, frag_len, count, rng,
            )
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
            r2_blocks = transcript_to_genomic_blocks(frag_start, r2_t_end, t)
            r1_t_start = max(frag_end - read_len, frag_start)
            r1_blocks = transcript_to_genomic_blocks(r1_t_start, frag_end, t)
            if not r2_blocks or not r1_blocks:
                continue

            if t.strand == Strand.POS:
                r2_is_rev, r1_is_rev = False, True
            else:
                r2_is_rev, r1_is_rev = True, False
            if flipped:
                r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev
                r1_blocks, r2_blocks = r2_blocks, r1_blocks

            r2_cigar = blocks_to_cigar(r2_blocks)
            r1_cigar = blocks_to_cigar(r1_blocks)
            r2_start = r2_blocks[0][0]
            r1_start = r1_blocks[0][0]

            leftmost = min(r1_blocks[0][0], r2_blocks[0][0])
            rightmost = max(r1_blocks[-1][1], r2_blocks[-1][1])
            tlen = rightmost - leftmost

            tags = self._nh1_tags
            if len(r2_blocks) > 1 or len(r1_blocks) > 1:
                xs_strand = "-" if t.strand == Strand.NEG else "+"
                tags = [("NH", 1), ("XS", xs_strand)]

            r1_flag = BASE_R1_FLAG
            r2_flag = BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= FLAG_MATE_REVERSE
                r2_flag |= FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_start <= r2_start else -tlen
            r2_tlen = -r1_tlen

            r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[i], r1_is_rev)
            r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[i], r2_is_rev)

            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_start, r1_cigar, ref_id, r2_start, r1_tlen, tags=tags,
            ))
            bam_fh.write(make_aligned_segment(
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

            g_start, g_end = premrna_to_genomic_interval(frag_start, frag_end, t)
            if t.strand == Strand.POS:
                r2_g_start = g_start
                r2_g_end = min(g_start + read_len, g_end)
                r1_g_start = max(g_end - read_len, g_start)
                r1_g_end = g_end
                r2_is_rev, r1_is_rev = False, True
            else:
                r1_g_start = g_start
                r1_g_end = min(g_start + read_len, g_end)
                r2_g_start = max(g_end - read_len, g_start)
                r2_g_end = g_end
                r2_is_rev, r1_is_rev = True, False
            if flipped:
                r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev
                r1_g_start, r2_g_start = r2_g_start, r1_g_start
                r1_g_end, r2_g_end = r2_g_end, r1_g_end

            tlen = g_end - g_start

            r1_flag = BASE_R1_FLAG
            r2_flag = BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= FLAG_MATE_REVERSE
                r2_flag |= FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_g_start <= r2_g_start else -tlen
            r2_tlen = -r1_tlen

            r1_read_len = r1_g_end - r1_g_start
            r2_read_len = r2_g_end - r2_g_start

            r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[i], r1_is_rev)
            r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[i], r2_is_rev)

            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_g_start, [(pysam.CMATCH, r1_read_len)],
                ref_id, r2_g_start, r1_tlen, tags=self._nh1_tags,
            ))
            bam_fh.write(make_aligned_segment(
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
        counts: dict[tuple[int, int], int] = {}
        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            chrom_eff = np.array(
                [
                    self.capture.partition("gdna", ref, ref_len, fl)
                    for ref, ref_len in zip(self._gdna_refs, self._gdna_ref_lengths)
                ],
                dtype=np.float64,
            )
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

        starts = self.capture.sample_starts("gdna", ref, chrom_len, fl, count, rng)
        # Strand: uniform 50/50, or — with overdispersion — by the fragment's exon-derived region's
        # shared Beta(a, a) +strand rate (0 = +/forward, 1 = −/reverse, matching the qname encoding).
        strand_regions = self._gdna_strand_regions.get(ref)
        if strand_regions is None:
            chrom_strands = rng.integers(0, 2, size=count)
        else:
            boundaries, p_plus = strand_regions
            reg_idx = np.clip(
                np.searchsorted(boundaries, starts, side="right") - 1, 0, len(p_plus) - 1
            )
            chrom_strands = (rng.random(count) >= p_plus[reg_idx]).astype(int)

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

                if is_neg:
                    r1_is_rev = False
                    r2_is_rev = True
                    r1_start_pos = start
                    r2_start_pos = end - read_len
                else:
                    r1_is_rev = True
                    r2_is_rev = False
                    r1_start_pos = end - read_len
                    r2_start_pos = start

                r1_flag = BASE_R1_FLAG
                r2_flag = BASE_R2_FLAG
                if r1_is_rev:
                    r1_flag |= FLAG_REVERSE
                if r2_is_rev:
                    r1_flag |= FLAG_MATE_REVERSE
                    r2_flag |= FLAG_REVERSE
                if r1_is_rev:
                    r2_flag |= FLAG_MATE_REVERSE

                tlen = end - start
                r1_tlen = tlen if r1_start_pos <= r2_start_pos else -tlen
                r2_tlen = -r1_tlen

                r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[j], r1_is_rev)
                r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[j], r2_is_rev)

                bam_fh.write(make_aligned_segment(
                    self._bam_header, qnames[j], r1_seq_str,
                    r1_flag, ref_id, r1_start_pos,
                    [(pysam.CMATCH, read_len)],
                    ref_id, r2_start_pos, r1_tlen,
                    tags=self._nh1_tags,
                ))
                bam_fh.write(make_aligned_segment(
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


def nrna_label_for_ratio(
    ratio: float | tuple[float, float],
    labels: list[str] | None,
    idx: int,
) -> str:
    if labels and idx < len(labels):
        return labels[idx].strip()
    if isinstance(ratio, tuple):
        lo, hi = ratio
        if lo == hi == 0.0:
            return "none"
        return f"range_{lo:g}_{hi:g}"
    if ratio == 0:
        return "none"
    return f"ratio_{ratio:g}"


def _build_nrna_pairs(
    cfg: WholeGenomeSimConfig,
    has_file_nrna: bool,
) -> list[tuple[str, str, float | tuple[float, float] | None, int]]:
    """Build nRNA sweep pairs.

    When the abundance file supplied explicit nRNA data, returns a
    single entry ``("file", None)`` — no spike-in.  Otherwise returns
    one entry per configured additive ratio.
    """
    if has_file_nrna:
        return [("file", "file", None, 0)]

    mode = cfg.nrna.mode
    pairs: list[tuple[str, str, float | tuple[float, float] | None, int]] = []
    if mode == "additive_ratio":
        for i, ratio in enumerate(cfg.nrna.ratios):
            label = nrna_label_for_ratio(ratio, cfg.nrna.ratio_labels, i)
            pairs.append((label, mode, ratio, i))
        return pairs
    if mode == "random_fraction":
        if cfg.nrna.ratio_ranges is None:
            raise ValueError("nrna.ratio_ranges is required for mode='random_fraction'")
        for i, ratio_range in enumerate(cfg.nrna.ratio_ranges):
            label = nrna_label_for_ratio(ratio_range, cfg.nrna.ratio_labels, i)
            pairs.append((label, mode, ratio_range, i))
        return pairs

    raise ValueError(f"Unknown nRNA simulation mode: {mode}")


# ═══════════════════════════════════════════════════════════════════
# Manifest
# ═══════════════════════════════════════════════════════════════════


def write_manifest(
    outdir: Path,
    cfg: WholeGenomeSimConfig,
    conditions: list[dict],
) -> None:
    """Write manifest.json summarizing all simulation outputs."""
    path = write_manifest_file(outdir, cfg, conditions)
    logger.info("Wrote manifest to %s", path)


# ═══════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════


def run_simulation(cfg: WholeGenomeSimConfig) -> list[dict]:
    """Run full simulation.

    1. Load genome + GTF -> transcripts
    2. Assign base abundances (total RNA) once
    3. Sweep nRNA settings x gdna_rates x strand_specificities
    4. Write manifest

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

    # 3. Build condition grid: nrna × gdna × strand_specificities × capture
    sim = cfg.simulation

    gdna_pairs: list[tuple[str, float]] = []
    for i, rate in enumerate(cfg.gdna.rates):
        label = gdna_label_for_rate(rate, cfg.gdna.rate_labels, i)
        gdna_pairs.append((label, rate))

    # gDNA strand-overdispersion sweep axis (one condition per value). None ⇒ a single value
    # (= strand_overdispersion), which keeps condition names unchanged when it is 0.
    _od_values = (
        cfg.gdna.strand_overdispersions
        if cfg.gdna.strand_overdispersions is not None
        else [cfg.gdna.strand_overdispersion]
    )
    _od_labels = cfg.gdna.strand_overdispersion_labels
    gdna_od_pairs: list[tuple[str, float]] = []
    for i, od in enumerate(_od_values):
        od = float(od)
        if not (0.0 <= od < 1.0):
            raise ValueError(f"gdna strand_overdispersion must be in [0, 1); got {od}")
        label = (
            _od_labels[i]
            if _od_labels is not None and i < len(_od_labels)
            else ("od00" if od == 0.0 else f"od{od:g}".replace(".", "p"))
        )
        gdna_od_pairs.append((label, od))

    nrna_pairs = _build_nrna_pairs(cfg, has_file_nrna)
    capture_scenarios, include_capture_in_names = _capture_grid(cfg)

    from .orchestrator import run_condition_grid  # local import avoids an import cycle

    conditions = run_condition_grid(
        outdir=outdir,
        genome_path=genome_path,
        transcripts=transcripts,
        base_abundances=base_abundances,
        sim=sim,
        gdna=cfg.gdna,
        nrna=cfg.nrna,
        gdna_pairs=gdna_pairs,
        gdna_od_pairs=gdna_od_pairs,
        strand_specificities=cfg.strand_specificities,
        nrna_pairs=nrna_pairs,
        capture_scenarios=capture_scenarios,
        include_capture_in_names=include_capture_in_names,
        base_seed=sim.sim_seed,
        oracle_bam=cfg.oracle_bam,
        skip_existing=True,
    )

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
    if cfg.nrna.mode == "additive_ratio":
        print(f"  nRNA ratios:      {cfg.nrna.ratios}", flush=True)
    elif cfg.nrna.mode == "random_fraction":
        print(f"  nRNA ranges:      {cfg.nrna.ratio_ranges}", flush=True)
        print(f"  nRNA eligible:    {cfg.nrna.eligible_fraction}", flush=True)
    else:
        print("  nRNA:             explicit file values", flush=True)
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
