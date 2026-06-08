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

    python scripts/sim/simulate_reads.py --config scripts/sim/configs/example_existing_reference.yaml

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
import logging
import sys
import time
from pathlib import Path

import numpy as np

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]

try:
    import pgzip
except ImportError:
    pgzip = None  # type: ignore[assignment]

from rigel.transcript import Transcript
# Config dataclasses live in wgs_config (data layer); re-exported so existing
# `from rigel.sim.whole_genome import SimulationParams` call sites keep working.
from .wgs_config import (  # noqa: F401
    AbundanceConfig,
    GDNASimConfig,
    NRNAConfig,
    SimulationParams,
    WholeGenomeSimConfig,
)
from .wgs_engine import GenomeCache, WholeGenomeSimulator  # noqa: F401  re-export
from rigel.sim.capture import CaptureConfig, CaptureScenario
from rigel.sim.manifest import (
    gdna_label_for_rate,
    write_manifest as write_manifest_file,
)

logger = logging.getLogger(__name__)

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
