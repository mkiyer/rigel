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

**file** — Load from a TSV with columns ``transcript_id``,
``mrna_abundance``, ``nrna_abundance``.

Nascent RNA
-----------
Nascent RNA is controlled separately by the ``nrna`` section, and it lives on the
rigel index's NASCENT-RNA ENTITIES (owner, 2026-08-19): every multi-exon transcript
links (``nrna_t_index``) to one single-exon entity over its TSS/TES-clustered span —
a synthetic transcript, or an annotated single-exon transcript that already covers
the span. Annotated multi-exon transcripts carry ``nrna_abundance = 0`` always.

Whatever sets those numbers, sampling is ONE multinomial over every RNA row — mature
rows and entity rows alike — with probability ∝ abundance × effective length, so the
nascent FRAGMENT share follows from the molecules and their lengths and is never
imposed. Three modes set them (``NRNAConfig``):

**sparse** (``abundance_ranges`` + ``on_fraction``) — ⭐ the mode the benchmark panels
use, and the one to read first. Nascent RNA is ABSENT from most gene spans and present
in a minority: each ENTITY is switched on with probability ``on_fraction``, and where
it is on its abundance is drawn LOG-UNIFORMLY over ``(lo, hi)``. That level is
**ABSOLUTE** and independent of the mature level, so ``nascent > mature`` occurs.
The fragment share is EMERGENT — priceable in advance with :func:`expected_rna_weights`
and recorded per condition by the orchestrator. See :func:`apply_sparse_nrna`.

**additive_ratio** (``ratios``) and **fragment_share** (``shares``) — the two RATIO
modes, where the entity's molecules are pooled from its contributors:

    entity.nrna_abundance = Σ over the entity's contributors of (contributor.abundance × nrna_ratio)

``additive_ratio`` states ``nrna_ratio`` directly; ``fragment_share`` states the
nascent share of RNA FRAGMENTS in the uncaptured library and SOLVES for the ratio that
produces it (:func:`apply_nrna_fragment_share`). ⛔ Under both, nascent mass tracks
mature abundance and can never exceed it, which is why the panels no longer use them.

Fragment allocation
-------------------
``n_rna_fragments`` is the whole RNA budget and gDNA is added on top; or, when
``n_total_fragments`` is set, the total is fixed and ``gdna_rate`` decides only the
split (:func:`rigel.sim.orchestrator.resolve_depths`, which the panels use):

    n_gdna = round(n_rna × gdna_rate)              # n_total_fragments unset
    n_rna  = round(total / (1 + gdna_rate));  n_gdna = total − n_rna

⛔ Nascent comes OUT of ``n_rna``, never on top of it, in every mode: the entities are
rows of the one RNA multinomial, so the mature/nascent split inside ``n_rna`` is
realised and read off the origin counts afterwards rather than allocated.

Condition grid
--------------
Sweeps: ``nrna × gdna_rates × gdna strand overdispersions × strand_specificities ×
capture scenarios``. The ``nrna`` axis is one condition per configured entry of
whichever key the mode reads — ``ratios``, ``shares`` or ``abundance_ranges``.

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
from rigel.types import Interval, Strand

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

    defaults = {key: value for key, value in raw.items() if key not in {"configs", "scenarios"}}
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
        scenarios.append(
            CaptureScenario(
                label=label,
                config=_capture_config_from_mapping(
                    item,
                    defaults,
                    label=label,
                    require_probes_when_enabled=True,
                ),
            )
        )
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
    cfg.shadow_gtf = raw.get("shadow_gtf", None)
    cfg.index = raw.get("index", None)
    cfg.outdir = raw.get("outdir", "sim_output")
    cfg.transcript_filter = raw.get("transcript_filter", "all")

    # Simulation params
    sim_raw = raw.get("simulation", {})
    sim = cfg.simulation
    sim.n_rna_fragments = int(sim_raw.get("n_rna_fragments", 1_000_000))
    # ⭐ Optional FIXED-TOTAL budget: when present, `gdna.rates` decides only the RNA/gDNA
    # SPLIT instead of adding gDNA on top (`orchestrator.resolve_depths`). Absent ⇒ legacy.
    _total = sim_raw.get("n_total_fragments")
    sim.n_total_fragments = None if _total is None else int(_total)
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
    if nrna.mode not in {"additive_ratio", "sparse", "fragment_share"}:
        raise ValueError("nrna.mode must be 'additive_ratio', 'fragment_share' or 'sparse'")
    raw_shares = nrna_raw.get("shares", None)
    if raw_shares is not None:
        nrna.shares = [float(x) for x in raw_shares]
    if nrna.mode == "fragment_share" and nrna.shares is None:
        raise ValueError("nrna.shares is required for mode='fragment_share'")
    if "fracs" in nrna_raw or "frac_labels" in nrna_raw:
        raise ValueError("nrna.fracs is no longer supported; use nrna.ratios")
    raw_ratios = nrna_raw.get("ratios", None)
    if raw_ratios is not None:
        nrna.ratios = [float(r) for r in raw_ratios]
    raw_abundance_ranges = nrna_raw.get("abundance_ranges", None)
    if raw_abundance_ranges is not None:
        for pair in raw_abundance_ranges:
            # ⛔ a range is a PAIR. A three-element entry used to be silently TRUNCATED to its first
            # two, so `[1, 10, 100]` ran as `(1, 10)` and the panel was quietly not the one configured.
            if not hasattr(pair, "__len__") or len(pair) != 2:
                raise ValueError(
                    f"nrna.abundance_ranges entries must be [lo, hi] pairs; got {pair!r}"
                )
        nrna.abundance_ranges = [(float(pair[0]), float(pair[1])) for pair in raw_abundance_ranges]
    nrna.ratio_labels = nrna_raw.get("ratio_labels", None)
    nrna.on_fraction = float(nrna_raw.get("on_fraction", 1.0))
    nrna.seed = int(nrna_raw.get("seed", 42))
    if nrna.mode == "sparse" and nrna.abundance_ranges is None:
        raise ValueError("nrna.abundance_ranges is required for mode='sparse'")
    # ⛔⛔ A FIELD THAT THE SELECTED MODE CANNOT READ IS A CONFIG THE AUTHOR DID NOT WRITE. Without
    # this, `abundance_ranges` + `on_fraction` with `mode:` omitted parsed CLEAN and ran
    # `additive_ratio` at ratio 0.0 — a NASCENT-FREE panel whose conditions were still named for the
    # nascent label. Silence is the whole defect: the numbers look like a panel, and are another one.
    _MODE_FIELDS = {
        "additive_ratio": ("ratios",),
        "fragment_share": ("shares",),
        "sparse": ("abundance_ranges", "on_fraction"),
    }
    _ignored = sorted(
        f
        for mode, fields in _MODE_FIELDS.items()
        if mode != nrna.mode
        for f in fields
        if f in nrna_raw and f not in _MODE_FIELDS[nrna.mode]
    )
    if _ignored:
        raise ValueError(
            f"nrna.mode={nrna.mode!r} cannot read {_ignored} — those fields belong to another mode "
            f"and would be silently ignored. Remove them, or set the mode they belong to."
        )
    if not 0.0 <= nrna.on_fraction <= 1.0:
        raise ValueError("nrna.on_fraction must be between 0 and 1")
    if nrna.mode == "additive_ratio":
        expected_len = len(nrna.ratios)
    elif nrna.mode == "fragment_share":
        for share in nrna.shares or []:
            if not 0.0 <= share < 1.0:
                raise ValueError("nrna.shares entries must satisfy 0 <= share < 1")
        expected_len = len(nrna.shares or [])
    else:
        for lo, hi in nrna.abundance_ranges or []:
            # ⛔ a LOG-uniform draw has no zero end: express "no nascent" with on_fraction = 0
            if lo <= 0 or hi < lo:
                raise ValueError("nrna.abundance_ranges entries must satisfy 0 < min <= max")
        expected_len = len(nrna.abundance_ranges or [])
    if nrna.ratio_labels is not None and len(nrna.ratio_labels) != expected_len:
        raise ValueError("nrna.ratio_labels must match the number of nRNA scenarios")

    # gDNA
    gd_raw = raw.get("gdna", {})
    gd = cfg.gdna
    gd.rates = [float(r) for r in gd_raw.get("rates", [0.0])]
    gd.rate_labels = gd_raw.get("rate_labels", None)
    _genomic_refs = gd_raw.get("genomic_refs", None)
    gd.genomic_refs = [str(ref) for ref in _genomic_refs] if _genomic_refs is not None else None
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
    cfg.strand_specificities = [float(s) for s in raw.get("strand_specificities", [1.0])]

    # Misc
    cfg.oracle_bam = bool(raw.get("oracle_bam", True))
    # ⭐ Drop sim_R1/R2.fq.gz after each condition's truth is written. No calibration
    # instrument reads a FASTQ and they are ~half a panel's on-disk size; ``scripts/benchmarking/``
    # DOES read them, so a panel built this way cannot be compared against another tool
    # without re-simulating. `sim.orchestrator.run_condition_grid`.
    cfg.emit_fastq = bool(raw.get("emit_fastq", True))
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


def merge_shadow_transcripts(
    transcripts: list[Transcript], shadow_gtf: str | Path
) -> list[Transcript]:
    """Append SHADOW transcripts — unannotated transcription the simulator draws from and the index
    never sees (owner design, 2026-08-29).

    ``transcripts`` is the index's list (annotated + nascent entities). The shadow GTF's transcripts are
    loaded with the ordinary GTF loader and appended with ``t_index`` continuing after the index's rows,
    ``is_nrna = is_synthetic = False`` and ``nrna_t_index = -1`` (a shadow carries no nascent: it is
    itself the unannotated RNA). ⛔ A shadow whose ``t_id`` the index already knows is REFUSED — it
    would be annotated, and the whole point is that the tool cannot know about it. Their fragments are
    named like any RNA fragment (``{t_id}:…``), so the oracle split files them as ``mrna`` and the
    certified per-slot truth shows RNA exactly where the annotation says there is none.
    Gate: ``tests/test_sim_shadow_transcripts.py``."""
    shadows = load_transcripts(shadow_gtf, transcript_filter="all")
    known = {t.t_id for t in transcripts}
    clash = sorted(t.t_id for t in shadows if t.t_id in known)
    if clash:
        raise ValueError(
            f"{len(clash)} shadow transcript(s) are KNOWN to the index and are therefore not shadows: "
            f"{clash[:5]}{'…' if len(clash) > 5 else ''} — remove them from the shadow GTF or the annotation"
        )
    base = max((t.t_index for t in transcripts), default=-1) + 1
    for k, t in enumerate(shadows):
        t.t_index = base + k
        t.is_nrna = False
        t.is_synthetic = False
        t.nrna_t_index = -1
        t.nrna_abundance = 0.0
    logger.info(
        "Shadow transcripts: %d appended from %s (unknown to the index by construction)",
        len(shadows),
        shadow_gtf,
    )
    return transcripts + shadows


def load_transcripts_from_index(index_dir: str | Path) -> list[Transcript]:
    """⭐ The simulation's transcriptome IS the rigel index's (owner, 2026-08-19).

    Rebuilds one :class:`Transcript` per index row — annotated transcripts AND the synthetic nascent
    entities, with ``is_nrna`` / ``is_synthetic`` / ``nrna_t_index`` / ``nrna_n_contributors`` exactly as
    the index carries them and ``t_index`` equal to the index row — so the simulator, the oracle and
    `rigel quant` all speak about the same transcript set. ⛔ The GTF loader is not used here because
    it cannot know what the index did (duplicate collapse, nascent consolidation); feeding the simulator
    the raw GTF is how its nascent model diverged from the tool's.
    """
    from ..index import TranscriptIndex

    index = TranscriptIndex.load(index_dir, retain_test_structures=True)
    t_df = index.t_df
    transcripts: list[Transcript] = []
    for row in t_df.itertuples(index=False):
        exons = index.get_exon_intervals(int(row.t_index))
        if exons is None or len(exons) == 0:
            raise ValueError(f"index row {row.t_index} ({row.t_id}) has no exon intervals")
        t = Transcript(
            ref=str(row.ref),
            strand=Strand(int(row.strand)),
            exons=[Interval(int(s), int(e)) for s, e in np.asarray(exons).tolist()],
            t_id=str(row.t_id),
            g_id=str(row.g_id),
            g_name=str(row.g_name) if row.g_name is not None else None,
            g_type=str(row.g_type) if row.g_type is not None else None,
            t_index=int(row.t_index),
            g_index=int(row.g_index),
            is_basic=bool(row.is_basic),
            is_mane=bool(row.is_mane),
            is_ccds=bool(row.is_ccds),
            is_nrna=bool(row.is_nrna),
            is_synthetic=bool(row.is_synthetic),
            nrna_t_index=int(row.nrna_t_index),
            nrna_n_contributors=int(row.nrna_n_contributors),
        )
        t.length = t.compute_length()
        transcripts.append(t)
    if [t.t_index for t in transcripts] != list(range(len(transcripts))):
        raise ValueError("index rows are not contiguous t_index 0..n-1")
    n_syn = sum(1 for t in transcripts if t.is_synthetic)
    logger.info(
        "Loaded %d transcripts from index %s (%d synthetic nascent entities)",
        len(transcripts),
        index_dir,
        n_syn,
    )
    return transcripts


def ensure_index(cfg: "WholeGenomeSimConfig") -> Path:
    """The index the simulation reads its transcriptome from: ``cfg.index`` if set, else built from
    ``genome`` + ``gtf`` into ``<outdir>/rigel_index`` (duplicates collapsed, as `panel.py build`
    does) — and reused on the next run."""
    from ..index import TranscriptIndex

    if cfg.index:
        index_dir = Path(cfg.index)
        if not index_dir.is_dir():
            raise FileNotFoundError(f"index not found: {index_dir}")
        return index_dir
    index_dir = Path(cfg.outdir) / "rigel_index"
    if not index_dir.is_dir():
        logger.info(
            "No index configured — building one from %s + %s into %s",
            cfg.genome,
            cfg.gtf,
            index_dir,
        )
        TranscriptIndex.build(
            cfg.genome,
            cfg.gtf,
            index_dir,
            write_tsv=False,
            collapse_duplicate_transcripts=True,
        )
    return index_dir


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
    # ⛔ A synthetic nascent entity is not a mature molecule: it never draws a mature abundance. The
    # draw runs over the ANNOTATED rows only, in index order, so the synthetic rows do not shift it.
    annotated = np.array([not t.is_synthetic for t in transcripts], dtype=bool)

    # Step 1: expressed?
    expressed = np.zeros(n, dtype=bool)
    expressed[annotated] = rng.random(int(annotated.sum())) < config.frac_expressed

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
        n_expr,
        n,
        total_mrna,
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


def assign_nrna_to_entities(transcripts: list[Transcript], per_contributor: np.ndarray) -> int:
    """⭐ Pool each contributor's nascent molecules onto its nRNA ENTITY (owner, 2026-08-19).

    ``per_contributor[i]`` is transcript ``i``'s nascent molecular abundance (zero for anything that
    is not an expressed multi-exon transcript). Every multi-exon transcript links to one entity through
    ``nrna_t_index`` — the index's synthetic single-exon transcript over the clustered span, or the
    annotated single-exon transcript that already covers it — and the entity's ``nrna_abundance`` is
    the SUM over its contributors. Annotated multi-exon transcripts end with ``nrna_abundance = 0``:
    their nascent molecules are the entity's, and are sampled on the entity's template.

    ⛔ A contributor with no entity (``nrna_t_index < 0``) is an index defect, not a case to skip.
    Returns the number of entities that received nascent.
    """
    by_index = {t.t_index: t for t in transcripts}
    for t in transcripts:
        t.nrna_abundance = 0.0
    n_entities = 0
    for t, amount in zip(transcripts, per_contributor):
        if amount <= 0.0:
            continue
        if len(t.exons) <= 1:
            raise ValueError(f"{t.t_id} is single-exon and cannot contribute nascent RNA")
        if t.nrna_t_index < 0 or t.nrna_t_index not in by_index:
            raise ValueError(
                f"{t.t_id} is multi-exon but links to no nascent entity (nrna_t_index="
                f"{t.nrna_t_index}); the transcript list did not come from a rigel index"
            )
        entity = by_index[t.nrna_t_index]
        if entity.nrna_abundance == 0.0:
            n_entities += 1
        entity.nrna_abundance += float(amount)
    return n_entities


def apply_nrna_ratio(
    transcripts: list[Transcript],
    ratio: float,
) -> None:
    """Nascent RNA at a MOLECULAR ratio of mature RNA: ``ratio`` nascent molecules per mature molecule
    of every expressed multi-exon transcript, pooled onto its nRNA entity
    (:func:`assign_nrna_to_entities`)."""
    per = np.array(
        [
            (t.abundance or 0.0) * ratio if (t.abundance or 0.0) > 0 and len(t.exons) > 1 else 0.0
            for t in transcripts
        ]
    )
    n_entities = assign_nrna_to_entities(transcripts, per)
    n_contrib = int(np.sum(per > 0))
    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    total_nrna = sum(t.nrna_abundance for t in transcripts)
    logger.info(
        "Set nRNA molecular ratio %.3g: %d contributors onto %d entities (mRNA=%.1f, nRNA=%.1f)",
        ratio,
        n_contrib,
        n_entities,
        total_mrna,
        total_nrna,
    )


def fl_pmf(sim: "SimulationParams") -> tuple[np.ndarray, np.ndarray]:
    """The fragment-length pmf the engine draws from — a Normal truncated to ``[frag_min, frag_max]``,
    evaluated exactly rather than sampled, so a quantity derived from it is deterministic.
    Matches `sampling.truncated_normal_frag_lengths`'s rejection loop in distribution."""
    widths = np.arange(int(sim.frag_min), int(sim.frag_max) + 1, dtype=np.int64)
    logp = -0.5 * ((widths - float(sim.frag_mean)) / float(sim.frag_std)) ** 2
    p = np.exp(logp - logp.max())
    return widths, p / p.sum()


def expected_rna_weights(
    transcripts: list[Transcript], sim: "SimulationParams"
) -> tuple[float, float]:
    """``(mature, nascent)`` expected sampling weight in the UNCAPTURED library:
    ``Σ_t abundance_t · E_w[max(0, L_t − w + 1)]``, the exact expectation over the fragment-length
    pmf rather than the effective length at the mean width (they differ wherever a transcript is
    short relative to the fragment length, which is most of them).
    """
    widths, p = fl_pmf(sim)
    L = np.array([float(t.length or t.compute_length()) for t in transcripts])
    eff = (np.maximum(0.0, L[:, None] - widths[None, :] + 1.0) * p[None, :]).sum(axis=1)
    am = np.array([float(t.abundance or 0.0) for t in transcripts])
    an = np.array([float(t.nrna_abundance) for t in transcripts])
    return float((am * eff).sum()), float((an * eff).sum())


def apply_nrna_fragment_share(
    transcripts: list[Transcript], share: float, sim: "SimulationParams"
) -> float:
    """⭐ Set nascent molecules so nascent takes ``share`` of the RNA **FRAGMENTS** in the uncaptured
    library, and return the molecular ratio that achieves it.

    ⛔ **Why a panel states the FRAGMENT share and not the molecular ratio.** A nascent entity spans a
    whole gene and a mature transcript is spliced — on the ladder's index, mean 40,667 bp against
    1,708 bp — so a molecular ratio of 0.25 puts **86 %** of RNA fragments in nascent RNA. The ratio
    that gives the 20 % the panel has always meant is **0.0100**, and it is a property of the
    annotation, not a number anyone should hand-write into a config
    (`TRAPS: no-magic-numbers`).

    Each expressed multi-exon transcript contributes ``ratio × abundance`` nascent molecules to its
    entity, so ``W_nascent`` is linear in the ratio and the solve is exact::

        share = c·W_n1 / (W_m + c·W_n1)   ⇒   c = (share / (1 − share)) · W_m / W_n1

    with ``W_n1`` the nascent weight at ratio 1. ⚠ Solved on UNCAPTURED lengths: it fixes the
    LIBRARY's molecular composition, and the realised share then moves under capture, which is
    physically right — capture acts on molecules that already exist.
    """
    if not 0.0 <= share < 1.0:
        raise ValueError(f"nrna share must be in [0, 1); got {share}")
    if share == 0.0:
        assign_nrna_to_entities(transcripts, np.zeros(len(transcripts)))
        return 0.0
    apply_nrna_ratio(transcripts, 1.0)
    w_mature, w_nascent_unit = expected_rna_weights(transcripts, sim)
    if w_nascent_unit <= 0.0:
        raise ValueError(
            "no nascent opportunity: no expressed multi-exon transcript has an nRNA entity, so a "
            f"nascent fragment share of {share} is unreachable"
        )
    ratio = (share / (1.0 - share)) * w_mature / w_nascent_unit
    apply_nrna_ratio(transcripts, ratio)
    logger.info(
        "nRNA fragment share %.4g ⇒ molecular ratio %.6g (W_mature=%.4g, W_nascent@1=%.4g)",
        share,
        ratio,
        w_mature,
        w_nascent_unit,
    )
    return ratio


def apply_sparse_nrna(
    transcripts: list[Transcript],
    abundance_range: tuple[float, float],
    *,
    on_fraction: float,
    seed: int,
) -> float:
    """⭐⭐⭐ **NASCENT RNA IS SPARSE: ABSENT FROM MOST GENE SPANS, PRESENT AND MEASURABLE IN A
    MINORITY** (owner, 2026-08-22). Draw, per nascent ENTITY, whether it is transcribed at all
    (Bernoulli ``on_fraction``) and, if it is, its **ABSOLUTE** molecular abundance LOG-UNIFORMLY over
    ``abundance_range``. Returns the realised nascent:mature molecular ratio.

    ⛔⛔ **THE LEVEL IS INDEPENDENT OF THE MATURE LEVEL, AND THAT IS THE POINT** (owner's ruling; it is
    what this mode changes). The retired ratio modes set ``nascent = mature x ratio``, so nascent mass
    tracked mature abundance and nascent could never exceed it. The steady state says the opposite:
    mature is synthesis/degradation and nascent is synthesis, so the nascent:mature ratio is a
    STABILITY parameter, not an expression one — an abundant stable transcript shows almost no nascent
    signal and an unstable rare one can show more nascent than mature. Drawing the two independently
    makes ``nascent > mature`` a real case the tool must survive.

    ⭐ **LOG-UNIFORM, because the levels span decades where nascent is present** — a linear draw on a
    range like (0.05, 2) puts 97 % of its mass in the top decade and cannot express "very low in some,
    high in others".

    ⛔⛔ **THE UNIT OF SPARSITY IS THE ENTITY, NOT THE TRANSCRIPT, and the difference is measurable.**
    An entity is one TSS/TES-clustered gene span and several isoforms share it, so drawing per
    contributor would give a gene with 5 isoforms ``1 - (1 - 0.1)^5 = 41 %`` chance of carrying nascent
    at ``on_fraction = 0.1`` — the INTRON slots, which is what calibration reads, would be four times
    less sparse than configured. Per entity, the configured fraction is the fraction of gene spans and
    therefore of intron slots. Pre-mRNA is a property of a locus being transcribed, which is the same
    unit.

    ⚠ The nascent FRAGMENT share is EMERGENT here rather than solved (owner: acceptable, since a sparse
    subset and a bounded range control it). It is not a free parameter — the caller can price it
    exactly with :func:`expected_rna_weights`, and the orchestrator records it per condition.
    """
    lo, hi = float(abundance_range[0]), float(abundance_range[1])
    if lo <= 0.0 or hi < lo:
        raise ValueError(
            f"abundance_range must satisfy 0 < lo <= hi for a LOG-uniform draw; got ({lo}, {hi})"
        )
    if not 0.0 <= on_fraction <= 1.0:
        raise ValueError(f"on_fraction must be between 0 and 1; got {on_fraction}")

    rng = np.random.default_rng(seed)
    by_index = {t.t_index: t for t in transcripts}
    for t in transcripts:
        t.nrna_abundance = 0.0

    # an entity is ELIGIBLE iff at least one EXPRESSED multi-exon transcript names it: a silent gene
    # is not being transcribed, so it has no pre-mRNA (`TRAPS: starved-is-not-depleted` — this is the
    # "biology puts nothing there" side, and it must read as an exact zero rather than a small number)
    eligible: list[Transcript] = []
    seen: set[int] = set()
    for t in transcripts:
        if (t.abundance or 0.0) <= 0.0 or len(t.exons) <= 1:
            continue
        entity = by_index.get(t.nrna_t_index)
        if entity is None:
            raise ValueError(
                f"transcript {t.t_id} has no nascent entity (nrna_t_index={t.nrna_t_index}); a "
                f"transcript list without entities is not a rigel index's"
            )
        if entity.t_index not in seen:
            seen.add(entity.t_index)
            eligible.append(entity)

    log_lo, log_hi = np.log(lo), np.log(hi)
    n_on = 0
    for entity in eligible:
        if rng.random() >= on_fraction:
            continue
        entity.nrna_abundance = float(np.exp(rng.uniform(log_lo, log_hi)) if hi > lo else lo)
        n_on += 1

    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    total_nrna = sum(t.nrna_abundance for t in transcripts)
    realized_ratio = total_nrna / total_mrna if total_mrna > 0 else 0.0
    logger.info(
        "Sparse nRNA: %d/%d gene spans ON (on_fraction=%.3g), level ~ logU(%.3g, %.3g), "
        "mRNA=%.1f, nRNA=%.1f, realised molar ratio=%.4g",
        n_on,
        len(eligible),
        on_fraction,
        lo,
        hi,
        total_mrna,
        total_nrna,
        realized_ratio,
    )
    return realized_ratio


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
    # ⭐ A file names ANNOTATED transcripts, and its nRNA column is that transcript's nascent
    # molecules — which live on its ENTITY (owner, 2026-08-19). Collect them per contributor and pool
    # in one pass, exactly as the ratio modes do; a file row naming an entity directly is honoured too.
    per_contributor = np.zeros(len(transcripts))

    for i, t in enumerate(transcripts):
        if t.t_id in abund_map:
            total_or_mrna, nrna = abund_map[t.t_id]
            t.abundance = total_or_mrna
            if nrna is not None:
                has_nrna_data = True
                if len(t.exons) > 1:
                    per_contributor[i] = nrna
                else:
                    t.nrna_abundance = nrna  # already an entity / single-exon row
            matched += 1
        else:
            t.abundance = 0.0
    if per_contributor.any():
        direct = {
            i: transcripts[i].nrna_abundance
            for i in range(len(transcripts))
            if transcripts[i].nrna_abundance > 0
        }
        assign_nrna_to_entities(transcripts, per_contributor)
        for i, amount in direct.items():  # keep any nascent named on an entity row directly
            transcripts[i].nrna_abundance += amount

    logger.info(
        "File abundances (%s): matched=%d/%d, has_nrna=%s",
        fmt,
        matched,
        len(transcripts),
        has_nrna_data,
    )

    return has_nrna_data


def write_truth_abundances(
    transcripts: list[Transcript],
    path: Path,
) -> None:
    """Write ground-truth abundances to a TSV file."""
    with open(path, "w") as f:
        f.write(
            "transcript_id\tgene_id\tgene_name\tref\tstrand\t"
            "mrna_abundance\tnrna_abundance\ttotal_rna\tn_exons\t"
            "spliced_length\tgenomic_span\n"
        )
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
# nRNA sweep-pair construction
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
    """Build nRNA sweep pairs — one ``(label, mode, value, index)`` per nRNA condition.

    When the abundance file supplied explicit nRNA data, returns the single entry
    ``("file", "file", None, 0)``: the TSV is the one source of nascent weight and the sweep
    is skipped entirely. Otherwise one entry per configured value of whichever key the mode
    reads — ``ratios`` for ``additive_ratio``, ``shares`` for ``fragment_share``, and an
    ``(lo, hi)`` pair from ``abundance_ranges`` for ``sparse``. ``index`` is the entry's
    position, which ``sparse`` also folds into its per-condition seed.
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
    if mode == "fragment_share":
        if cfg.nrna.shares is None:
            raise ValueError("nrna.shares is required for mode='fragment_share'")
        for i, share in enumerate(cfg.nrna.shares):
            label = nrna_label_for_ratio(share, cfg.nrna.ratio_labels, i)
            pairs.append((label, mode, share, i))
        return pairs
    if mode == "sparse":
        if cfg.nrna.abundance_ranges is None:
            raise ValueError("nrna.abundance_ranges is required for mode='sparse'")
        for i, abundance_range in enumerate(cfg.nrna.abundance_ranges):
            label = nrna_label_for_ratio(abundance_range, cfg.nrna.ratio_labels, i)
            pairs.append((label, mode, abundance_range, i))
        return pairs

    raise ValueError(f"Unknown nRNA simulation mode: {mode}")


# ═══════════════════════════════════════════════════════════════════
# Manifest
# ═══════════════════════════════════════════════════════════════════


def write_manifest(
    outdir: Path,
    cfg: WholeGenomeSimConfig,
    conditions: list[dict],
    *,
    n_shadow_transcripts: int = 0,
) -> None:
    """Write manifest.json summarizing all simulation outputs."""
    path = write_manifest_file(outdir, cfg, conditions, n_shadow_transcripts=n_shadow_transcripts)
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

    # 1. Load transcripts — from the rigel INDEX, so the simulated transcriptome (annotated transcripts
    #    plus the synthetic nascent entities) is exactly the one `rigel quant` reads.
    if cfg.transcript_filter != "all":
        raise ValueError(
            "transcript_filter is not applicable when the transcriptome comes from a rigel index: the "
            "index IS the transcript set. Build the index from a filtered GTF instead."
        )
    index_dir = ensure_index(cfg)
    transcripts = load_transcripts_from_index(index_dir)
    if not transcripts:
        raise RuntimeError(f"No transcripts loaded from index {index_dir}")
    n_shadow = 0
    if cfg.shadow_gtf:
        if not Path(cfg.shadow_gtf).exists():
            raise FileNotFoundError(f"shadow GTF not found: {cfg.shadow_gtf}")
        n_before = len(transcripts)
        transcripts = merge_shadow_transcripts(transcripts, cfg.shadow_gtf)
        n_shadow = len(transcripts) - n_before

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
        logger.info("Abundance file provided nRNA data — skipping nRNA spike-in sweep")

    # Save base abundances
    base_abundances = [(t.abundance or 0.0, t.nrna_abundance) for t in transcripts]

    # 3. Build condition grid: nrna × gdna × strand_specificities × capture
    sim = cfg.simulation

    gdna_pairs: list[tuple[str, float]] = []
    for i, rate in enumerate(cfg.gdna.rates):
        label = gdna_label_for_rate(rate, cfg.gdna.rate_labels, i)
        gdna_pairs.append((label, rate))

    # ⭐ The genomic/RNA-only split is a SCENARIO input. A config that asks for gDNA without stating
    # which references are genomic is rejected rather than guessed at, because the guess that used to
    # live in the engine — "has an annotation" — put gDNA on RNA-only spike-ins.
    genomic_refs = cfg.gdna.genomic_refs
    if genomic_refs is None:
        if any(rate > 0 for rate in cfg.gdna.rates):
            raise ValueError(
                "gdna.genomic_refs must list the references genomic DNA is drawn from whenever any "
                "gdna.rates entry is non-zero. Every reference is genomic or RNA-only and the "
                "classification is an input, not an inference from the annotation."
            )
        genomic_refs = []

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
        genomic_refs=genomic_refs,
        gdna_pairs=gdna_pairs,
        gdna_od_pairs=gdna_od_pairs,
        strand_specificities=cfg.strand_specificities,
        nrna_pairs=nrna_pairs,
        capture_scenarios=capture_scenarios,
        include_capture_in_names=include_capture_in_names,
        base_seed=sim.sim_seed,
        oracle_bam=cfg.oracle_bam,
        skip_existing=True,
        emit_fastq=getattr(cfg, "emit_fastq", True),
    )

    # 4. Write manifest
    write_manifest(outdir, cfg, conditions, n_shadow_transcripts=n_shadow)
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
    p.add_argument(
        "-j",
        "--n-workers",
        type=int,
        default=None,
        help="Worker processes for parallel read generation (overrides YAML)",
    )
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
    # ⚠ Print the CONFIGURED nRNA mode, whichever it is. The `fragment_share` branch used to fall
    # through to "explicit file values" — right for a config whose abundance TSV supplies nascent
    # weights and wrong for one that does not, and the two are not distinguishable here: whether the
    # sweep is skipped is decided later, by whether the loaded file carried an `nrna_abundance` column.
    if cfg.nrna.mode == "additive_ratio":
        print(f"  nRNA ratios:      {cfg.nrna.ratios}", flush=True)
    elif cfg.nrna.mode == "fragment_share":
        print(
            f"  nRNA frag shares: {cfg.nrna.shares} (molecular ratio SOLVED per share)", flush=True
        )
    elif cfg.nrna.mode == "sparse":
        print(
            f"  nRNA ranges:      {cfg.nrna.abundance_ranges} (LOG-uniform, absolute)", flush=True
        )
        print(f"  nRNA on_fraction: {cfg.nrna.on_fraction} of gene spans", flush=True)
    if cfg.abundance.mode == "file":
        print(
            "  nRNA:             the sweep above is SKIPPED if the abundance file supplies "
            "nrna_abundance",
            flush=True,
        )
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
