"""Whole-genome read simulation from one YAML config: parsing, the transcript set, abundances, and the run.

:func:`run_simulation` loads the transcript set from the rigel index (:func:`load_transcripts_from_index`: annotated
transcripts plus the synthetic nascent entities, ``t_index`` equal to the index row), appends any shadow transcripts,
assigns the base abundances (a log-uniform draw, or a salmon / kallisto / sim file) and hands the condition grid to
:func:`rigel.sim.orchestrator.run_condition_grid`. It writes ``manifest.json`` last.

Nascent abundance lives on the entities: every multi-exon transcript links through ``span_t_index`` to the entity
over its TSS/TES cluster and keeps ``nrna_abundance = 0`` itself. Sampling is one multinomial over all RNA rows,
transcript and entity alike, in proportion to abundance x effective length, so the nascent fragment share follows
from molecules and lengths and is never imposed.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import numpy as np
import yaml
from scipy.special import ndtr

from rigel.transcript import Transcript
from rigel.types import Interval, Strand

from .wgs_config import AbundanceConfig, SimulationParams, WholeGenomeSimConfig
from rigel.sim.capture import CaptureConfig, CaptureScenario
from rigel.sim.manifest import (
    gdna_label_for_rate,
    write_manifest,
)

logger = logging.getLogger(__name__)


#: Every key a capture mapping may carry: the scenario's own (``label``, ``enabled``) and the
#: ``CaptureConfig`` fields.
_CAPTURE_KEYS = frozenset(
    {
        "label",
        "enabled",
        "probes",
        "probe_format",
        "off_target_weight",
        "binding_per_base",
        "min_overlap",
    }
)


#: Every top-level key :func:`parse_yaml_config` reads.
_TOP_KEYS = frozenset(
    {
        "genome",
        "gtf",
        "shadow_gtf",
        "index",
        "outdir",
        "simulation",
        "abundance",
        "nrna",
        "gdna",
        "capture",
        "strand_specificities",
        "oracle_bam",
        "emit_fastq",
        "verbose",
    }
)
#: Every key :func:`parse_yaml_config` reads in each plain section (``capture`` checks its own).
_SECTION_KEYS = {
    "simulation": frozenset(
        {
            "n_rna_fragments",
            "n_total_fragments",
            "sim_seed",
            "frag_mean",
            "frag_std",
            "frag_min",
            "frag_max",
            "read_length",
            "error_rate",
            "n_workers",
        }
    ),
    "abundance": frozenset({"mode", "seed", "min", "max", "frac_expressed", "file"}),
    "nrna": frozenset(
        {"mode", "ratios", "abundance_ranges", "ratio_labels", "on_fraction", "seed"}
    ),
    "gdna": frozenset(
        {
            "rates",
            "rate_labels",
            "genomic_refs",
            "frag_mean",
            "frag_std",
            "frag_min",
            "frag_max",
            "strand_overdispersion",
            "strand_overdispersions",
            "strand_overdispersion_labels",
        }
    ),
}


def _refuse_unknown_keys(raw: dict, known: frozenset[str], where: str) -> None:
    """Raise on a key the parser does not read, since a misspelt key would otherwise run on its default."""
    unknown = sorted(set(raw) - known)
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {unknown}; known: {sorted(known)}")


#: A variable reference in a path, ``$NAME`` or ``${NAME}``, as :func:`os.path.expandvars` reads it (an unclosed
#: ``${`` included), so the check sees every reference the expansion would leave in place.
_PATH_VARIABLE = re.compile(r"\$(\w+|\{[^}]*\}?)", re.ASCII)


def expand_config_path(value, key: str):
    """A config path with ``~`` and every ``$NAME`` / ``${NAME}`` expanded, so a config can name its data root by a
    variable. An unset or empty variable is refused by name: left in place it would make a literal ``$NAME`` path,
    and expanded to nothing it would move the path to ``/``. A non-string (an absent key's ``None``) is returned
    as is."""
    if not isinstance(value, str):
        return value
    names = {ref.strip("{}") for ref in _PATH_VARIABLE.findall(value)}
    unset = sorted(name for name in names if not os.environ.get(name))
    if unset:
        raise ValueError(
            f"{key}: {value!r} names an environment variable that is unset or empty "
            f"({', '.join('$' + name for name in unset)}); export it, or write the path out"
        )
    return os.path.expanduser(os.path.expandvars(value))


def _refuse_a_run_that_keeps_no_reads(cfg: WholeGenomeSimConfig) -> None:
    """Raise when ``oracle_bam`` and ``emit_fastq`` are both false: no condition would keep a read, and without
    either file a condition has no resume key."""
    if not cfg.oracle_bam and not cfg.emit_fastq:
        raise ValueError(
            "oracle_bam and emit_fastq are both false (--no-oracle sets oracle_bam false): no condition "
            "would keep a read, since no oracle BAM is written and the FASTQs are deleted once the truth "
            "is, and every run would simulate every condition again. Keep oracle_bam true (the reads "
            "every instrument scores) or emit_fastq true."
        )


def _capture_config_from_mapping(
    raw: dict,
    defaults: dict | None = None,
    *,
    label: str = "capture",
    require_probes_when_enabled: bool = False,
) -> CaptureConfig:
    """Build a CaptureConfig from a YAML mapping plus optional defaults. An unknown key is refused, since ignored it
    would read as applied."""
    merged = dict(defaults or {})
    merged.update(raw)
    unknown = sorted(set(merged) - _CAPTURE_KEYS)
    if unknown:
        raise ValueError(
            f"capture config '{label}': unknown key(s) {unknown}; known: {sorted(_CAPTURE_KEYS)}"
        )
    enabled = bool(merged.get("enabled", True))
    if not enabled:
        return CaptureConfig()

    probes = merged.get("probes")
    if not probes:
        if require_probes_when_enabled and raw.get("enabled") is True:
            raise ValueError(f"capture config '{label}' is enabled but has no probes")
        return CaptureConfig()

    return CaptureConfig(
        probes=expand_config_path(str(probes), f"capture config '{label}': probes"),
        probe_format=str(merged.get("probe_format", CaptureConfig.probe_format)),
        off_target_weight=float(merged.get("off_target_weight", CaptureConfig.off_target_weight)),
        binding_per_base=float(merged.get("binding_per_base", CaptureConfig.binding_per_base)),
        min_overlap=int(merged.get("min_overlap", CaptureConfig.min_overlap)),
    )


def _capture_label_text(value: object) -> str:
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value)


def _capture_scenarios_from_mapping(raw: dict) -> list[CaptureScenario]:
    """Parse capture.configs from a YAML mapping."""
    raw_configs = raw.get("configs")
    if raw_configs is None:
        return []
    if not isinstance(raw_configs, list) or not raw_configs:
        raise ValueError("capture.configs must be a non-empty list")

    defaults = {key: value for key, value in raw.items() if key != "configs"}
    scenarios: list[CaptureScenario] = []
    seen_labels: set[str] = set()
    for index, item in enumerate(raw_configs, start=1):
        if not isinstance(item, dict):
            raise ValueError("each capture.configs entry must be a mapping")
        label = _capture_label_text(item.get("label", f"c{index}"))
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
    with open(path) as f:
        raw = yaml.safe_load(f)
    _refuse_unknown_keys(raw, _TOP_KEYS, "sim config")
    for section, known in _SECTION_KEYS.items():
        _refuse_unknown_keys(raw.get(section, {}), known, section)

    cfg = WholeGenomeSimConfig()
    cfg.genome = expand_config_path(raw.get("genome", cfg.genome), "genome")
    cfg.gtf = expand_config_path(raw.get("gtf", cfg.gtf), "gtf")
    cfg.shadow_gtf = expand_config_path(raw.get("shadow_gtf", cfg.shadow_gtf), "shadow_gtf")
    cfg.index = expand_config_path(raw.get("index", cfg.index), "index")
    cfg.outdir = expand_config_path(raw.get("outdir", cfg.outdir), "outdir")

    # Simulation params
    sim_raw = raw.get("simulation", {})
    sim = cfg.simulation
    sim.n_rna_fragments = int(sim_raw.get("n_rna_fragments", sim.n_rna_fragments))
    _total = sim_raw.get("n_total_fragments")
    sim.n_total_fragments = None if _total is None else int(_total)
    sim.sim_seed = int(sim_raw.get("sim_seed", sim.sim_seed))
    sim.frag_mean = float(sim_raw.get("frag_mean", sim.frag_mean))
    sim.frag_std = float(sim_raw.get("frag_std", sim.frag_std))
    sim.frag_min = int(sim_raw.get("frag_min", sim.frag_min))
    sim.frag_max = int(sim_raw.get("frag_max", sim.frag_max))
    sim.read_length = int(sim_raw.get("read_length", sim.read_length))
    sim.error_rate = float(sim_raw.get("error_rate", sim.error_rate))
    sim.n_workers = int(sim_raw.get("n_workers", sim.n_workers))

    # Abundance
    ab_raw = raw.get("abundance", {})
    ab = cfg.abundance
    ab.mode = ab_raw.get("mode", ab.mode)
    ab.seed = int(ab_raw.get("seed", ab.seed))
    ab.min = float(ab_raw.get("min", ab.min))
    ab.max = float(ab_raw.get("max", ab.max))
    ab.frac_expressed = float(ab_raw.get("frac_expressed", ab.frac_expressed))
    ab.file = expand_config_path(ab_raw.get("file", ab.file), "abundance.file")

    # Nascent RNA sweep
    nrna_raw = raw.get("nrna", {})
    nrna = cfg.nrna
    nrna.mode = str(nrna_raw.get("mode", nrna.mode))
    if nrna.mode not in {"additive_ratio", "sparse"}:
        raise ValueError("nrna.mode must be 'additive_ratio' or 'sparse'")
    raw_ratios = nrna_raw.get("ratios", None)
    if raw_ratios is not None:
        nrna.ratios = [float(r) for r in raw_ratios]
    raw_abundance_ranges = nrna_raw.get("abundance_ranges", None)
    if raw_abundance_ranges is not None:
        for pair in raw_abundance_ranges:
            # a range must be a pair: truncating a longer entry would run another scenario
            if not hasattr(pair, "__len__") or len(pair) != 2:
                raise ValueError(
                    f"nrna.abundance_ranges entries must be [lo, hi] pairs; got {pair!r}"
                )
        nrna.abundance_ranges = [(float(pair[0]), float(pair[1])) for pair in raw_abundance_ranges]
    nrna.ratio_labels = nrna_raw.get("ratio_labels", nrna.ratio_labels)
    nrna.on_fraction = float(nrna_raw.get("on_fraction", nrna.on_fraction))
    nrna.seed = int(nrna_raw.get("seed", nrna.seed))
    if nrna.mode == "sparse" and nrna.abundance_ranges is None:
        raise ValueError("nrna.abundance_ranges is required for mode='sparse'")
    # A field the selected mode cannot read is refused: `abundance_ranges` + `on_fraction` with `mode:` omitted
    # would otherwise run `additive_ratio` at ratio 0.0, a nascent-free run named like the configured one.
    _MODE_FIELDS = {
        "additive_ratio": ("ratios",),
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
    else:
        for lo, hi in nrna.abundance_ranges or []:
            # a log-uniform draw has no zero end: express "no nascent" with on_fraction = 0
            if lo <= 0 or hi < lo:
                raise ValueError("nrna.abundance_ranges entries must satisfy 0 < min <= max")
        expected_len = len(nrna.abundance_ranges or [])
    if nrna.ratio_labels is not None and len(nrna.ratio_labels) != expected_len:
        raise ValueError("nrna.ratio_labels must match the number of nRNA scenarios")

    # gDNA
    gd_raw = raw.get("gdna", {})
    gd = cfg.gdna
    gd.rates = [float(r) for r in gd_raw.get("rates", gd.rates)]
    gd.rate_labels = gd_raw.get("rate_labels", gd.rate_labels)
    _genomic_refs = gd_raw.get("genomic_refs", None)
    gd.genomic_refs = [str(ref) for ref in _genomic_refs] if _genomic_refs is not None else None
    gd.frag_mean = float(gd_raw.get("frag_mean", gd.frag_mean))
    gd.frag_std = float(gd_raw.get("frag_std", gd.frag_std))
    gd.frag_min = int(gd_raw.get("frag_min", gd.frag_min))
    gd.frag_max = int(gd_raw.get("frag_max", gd.frag_max))
    gd.strand_overdispersion = float(gd_raw.get("strand_overdispersion", gd.strand_overdispersion))
    _ods = gd_raw.get("strand_overdispersions", None)
    gd.strand_overdispersions = [float(o) for o in _ods] if _ods is not None else None
    gd.strand_overdispersion_labels = gd_raw.get(
        "strand_overdispersion_labels", gd.strand_overdispersion_labels
    )

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
    # False drops each condition's FASTQs once its truth is written, saving most of a panel's disk; a panel
    # without FASTQs cannot be re-aligned or given to another quantifier without re-simulating.
    cfg.emit_fastq = bool(raw.get("emit_fastq", True))
    _refuse_a_run_that_keeps_no_reads(cfg)
    cfg.verbose = bool(raw.get("verbose", True))

    return cfg


# ═══════════════════════════════════════════════════════════════════
# Transcript loading
# ═══════════════════════════════════════════════════════════════════


def load_transcripts(gtf_path: str | Path) -> list[Transcript]:
    """Load transcripts from a GTF file.

    Parameters
    ----------
    gtf_path : path
        GTF annotation file (may be gzipped).

    Returns
    -------
    list[Transcript]
        With ``t_index`` assigned sequentially.
    """
    logger.info("Loading transcripts from %s", gtf_path)
    transcripts = Transcript.read_gtf(str(gtf_path), parse_mode="warn-skip")
    logger.info("Read %d transcripts from GTF", len(transcripts))

    for i, t in enumerate(transcripts):
        t.t_index = i

    n_genes = len({t.g_id for t in transcripts})
    logger.info("Final: %d transcripts from %d genes", len(transcripts), n_genes)
    return transcripts


def merge_shadow_transcripts(
    transcripts: list[Transcript], shadow_gtf: str | Path
) -> list[Transcript]:
    """Append shadow transcripts, unannotated transcription the simulator draws from and the index never sees.

    The shadow GTF's rows follow the index's (``transcripts``) in ``t_index``, with ``is_synthetic = False``,
    ``span_t_index = -1`` and no nascent. A shadow whose ``t_id`` the index knows is refused, since the tool could
    know it. Their fragments are named like any RNA fragment, so the truth shows RNA where the annotation has none.
    """
    shadows = load_transcripts(shadow_gtf)
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
        t.is_synthetic = False
        t.span_t_index = -1
        t.nrna_abundance = 0.0
    logger.info(
        "Shadow transcripts: %d appended from %s (unknown to the index by construction)",
        len(shadows),
        shadow_gtf,
    )
    return transcripts + shadows


def load_transcripts_from_index(index_dir: str | Path) -> list[Transcript]:
    """One :class:`Transcript` per row of the rigel index, annotated transcripts and nascent entities, with the
    index's structure flags and ``t_index`` equal to the row, so the simulator and `rigel quant` see the same set.

    The index, not the GTF, is read because only it knows the duplicate collapse and the nascent entities. Each
    transcript carries the ids of the duplicates folded into it (``aliases``), so a probe panel designed on the raw
    annotation still binds. Raises if the rows are not contiguous ``t_index`` 0..n-1 or a row has no exons.
    """
    import pandas as pd

    from ..index import TRANSCRIPT_ALIASES_FEATHER, TranscriptIndex

    index = TranscriptIndex.load(index_dir, retain_lookup_maps=True)
    t_df = index.t_df
    aliases: dict[str, list[str]] = {}
    alias_path = Path(index_dir) / TRANSCRIPT_ALIASES_FEATHER
    if alias_path.exists():
        for dropped, kept in pd.read_feather(alias_path)[["dropped_t_id", "kept_t_id"]].itertuples(
            index=False
        ):
            aliases.setdefault(str(kept), []).append(str(dropped))
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
            is_single_exon=bool(row.is_single_exon),
            is_synthetic=bool(row.is_synthetic),
            span_t_index=int(row.span_t_index),
            span_n_contributors=int(row.span_n_contributors),
            aliases=tuple(sorted(aliases.get(str(row.t_id), ()))),
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
    """The index the transcriptome comes from: ``cfg.index`` if set, else one built from ``genome`` + ``gtf`` into
    ``<outdir>/rigel_index`` (duplicates collapsed) and reused on the next run."""
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


def _index_provenance(index_dir: Path, gtf: str | Path) -> dict:
    """The index record ``manifest.json`` carries: its resolved path and its own manifest verbatim.

    Refuses a ``gtf`` whose sha256 is not the one the index recorded, so a config cannot name one annotation and
    simulate another; a byte-identical copy elsewhere passes. Gate: ``tests/test_whole_genome_sim_config.py``."""
    from ..index import load_manifest, source_record

    index_manifest = load_manifest(index_dir)
    built_from = index_manifest["sources"]["gtf"]
    given = source_record(gtf)
    if given["sha256"] != built_from["sha256"]:
        raise ValueError(
            f"gtf {given['path']} is not the GTF the index {Path(index_dir).resolve()} was built from "
            f"({built_from['path']}): sha256 {given['sha256']} against {built_from['sha256']}. The "
            "simulated transcriptome is the index's, so gtf: must name that GTF; or build the index "
            "from this one."
        )
    return {"path": str(Path(index_dir).resolve()), "manifest": index_manifest}


# ═══════════════════════════════════════════════════════════════════
# Abundance assignment
# ═══════════════════════════════════════════════════════════════════


def assign_random_abundances(
    transcripts: list[Transcript],
    config: AbundanceConfig,
) -> None:
    """Draw each annotated transcript's ``abundance``: expressed with probability ``frac_expressed``, then
    LogUniform(min, max). ``nrna_abundance`` is set to 0; the nascent mode is applied per condition.
    """
    if config.min <= 0:
        raise ValueError(f"abundance min must be > 0 for log-uniform, got {config.min}")

    rng = np.random.default_rng(config.seed)
    n = len(transcripts)
    # The draw runs over the annotated rows only, in index order: a nascent entity has no mature abundance, and
    # the synthetic rows must not shift the stream.
    annotated = np.array([not t.is_synthetic for t in transcripts], dtype=bool)

    expressed = np.zeros(n, dtype=bool)
    expressed[annotated] = rng.random(int(annotated.sum())) < config.frac_expressed

    total_rna = np.zeros(n)
    n_expr = int(expressed.sum())
    log_min = np.log(config.min)
    log_max = np.log(config.max)
    total_rna[expressed] = np.exp(rng.uniform(log_min, log_max, size=n_expr))

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
    """Detect an abundance file's format, plain or gzipped, and return ``(format_name, dataframe)``:

    salmon   — ``quant.sf``: ``Name``, ``TPM``.
    kallisto — ``abundance.tsv``: ``target_id``, ``tpm``.
    sim      — the molar truth table (:func:`rigel.sim.truth.write_molar_truth`): ``transcript_id``,
               ``mrna_abundance`` and optionally ``nrna_abundance``.

    A condition's truth table is refused: it holds realised fragments, not molecules.
    """
    import pandas as pd

    # gzip is detected by its magic bytes, whatever the extension
    compression: str = "infer"
    try:
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                compression = "gzip"
    except OSError:
        pass

    df = pd.read_csv(path, sep="\t", compression=compression)
    cols = set(df.columns)

    # A condition's table is known by its realised column, whatever else it carries.
    if "observed_mrna_fragments" in cols:
        raise ValueError(
            f"{path} is a condition's truth table, the fragments one condition realised, not an "
            "abundance input; the molar table is the run's truth_abundances_nrna_<label>.tsv"
        )
    if {"Name", "TPM"}.issubset(cols):
        return "salmon", df
    if {"target_id", "tpm"}.issubset(cols):
        return "kallisto", df
    if {"transcript_id", "mrna_abundance"}.issubset(cols):
        return "sim", df
    # kallisto, case-insensitively
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
    """Load ``{transcript_id: (abundance, nrna or None)}`` and the detected format; the nRNA value is ``None`` when
    the file has none."""
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
    """Pool each contributor's nascent molecules onto its nascent entity; return how many entities received any.

    ``per_contributor[i]`` is transcript ``i``'s nascent molecular abundance. Each multi-exon transcript links
    through ``span_t_index`` to one single-exon entity (synthetic, or an annotated transcript covering the span),
    whose ``nrna_abundance`` becomes the sum over its contributors; every other row ends at 0. Raises when a
    contributor is single-exon or has no entity: the list did not come from a rigel index.
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
        if t.span_t_index < 0 or t.span_t_index not in by_index:
            raise ValueError(
                f"{t.t_id} is multi-exon but links to no nascent entity (span_t_index="
                f"{t.span_t_index}); the transcript list did not come from a rigel index"
            )
        entity = by_index[t.span_t_index]
        if entity.nrna_abundance == 0.0:
            n_entities += 1
        entity.nrna_abundance += float(amount)
    return n_entities


def apply_nrna_ratio(
    transcripts: list[Transcript],
    ratio: float,
) -> None:
    """Nascent RNA at a molecular ratio: ``ratio`` nascent molecules per molecule of every expressed multi-exon
    transcript, pooled onto its entity (:func:`assign_nrna_to_entities`)."""
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
    """The exact fragment-length pmf the engine draws from: `sampling.truncated_normal_frag_lengths` rounds a Normal
    draw and keeps it in ``[frag_min, frag_max]``, so length ``k`` has the Normal's mass on ``[k − ½, k + ½)``,
    renormalised over the window. Each mass is taken on the side of the mean where it is a difference of two small
    tails, never of two numbers near 1."""
    widths = np.arange(int(sim.frag_min), int(sim.frag_max) + 1, dtype=np.int64)
    lo = (widths - 0.5 - float(sim.frag_mean)) / float(sim.frag_std)
    hi = lo + 1.0 / float(sim.frag_std)
    p = np.where(lo > 0.0, ndtr(-lo) - ndtr(-hi), ndtr(hi) - ndtr(lo))
    return widths, p / p.sum()


def apply_sparse_nrna(
    transcripts: list[Transcript],
    abundance_range: tuple[float, float],
    *,
    on_fraction: float,
    seed: int,
) -> float:
    """Sparse nascent RNA: each eligible entity is transcribed with probability ``on_fraction`` and, if so, draws an
    absolute molecular abundance log-uniformly over ``abundance_range``. Returns the realised nascent:mature
    molecular ratio. Requires ``0 < lo <= hi`` and ``0 <= on_fraction <= 1``.

    The level is independent of the mature level, because the nascent:mature ratio reflects stability, not
    expression: a stable abundant transcript shows little nascent signal and an unstable rare one can show more
    nascent than mature. The draw is log-uniform because present levels span decades. The unit is the entity, not
    the transcript, so ``on_fraction`` is the fraction of gene spans: drawn per isoform, a gene with ``n`` isoforms
    would carry nascent with probability ``1 - (1 - p)^n``.
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

    # an entity is eligible iff an expressed multi-exon transcript names it: a silent gene's span must be exactly 0
    eligible: list[Transcript] = []
    seen: set[int] = set()
    for t in transcripts:
        if (t.abundance or 0.0) <= 0.0 or len(t.exons) <= 1:
            continue
        entity = by_index.get(t.span_t_index)
        if entity is None:
            raise ValueError(
                f"transcript {t.t_id} has no nascent entity (span_t_index={t.span_t_index}); a "
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
    """Load abundances from a salmon, kallisto or sim file (:func:`_detect_abundance_format`); a transcript the file
    does not name gets 0. Returns whether the file supplied nRNA values; when it did not, the caller applies the
    configured nascent sweep.
    """
    abund_map, fmt = _load_abundance_map(tsv_path)
    logger.info("Detected abundance format: %s (%d entries)", fmt, len(abund_map))

    matched = 0
    has_nrna_data = False
    # A file's nRNA value for a multi-exon transcript is pooled onto its entity, as in the ratio mode; a value on an
    # entity or single-exon row is kept directly.
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
                    t.nrna_abundance = nrna
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
        for i, amount in direct.items():
            transcripts[i].nrna_abundance += amount

    logger.info(
        "File abundances (%s): matched=%d/%d, has_nrna=%s",
        fmt,
        matched,
        len(transcripts),
        has_nrna_data,
    )

    return has_nrna_data


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
        return f"range_{lo:g}_{hi:g}"
    if ratio == 0:
        return "none"
    return f"ratio_{ratio:g}"


def _build_nrna_pairs(
    cfg: WholeGenomeSimConfig,
    has_file_nrna: bool,
) -> list[tuple[str, str, float | tuple[float, float] | None, int]]:
    """One ``(label, mode, value, index)`` per nascent condition: per ``ratios`` entry for ``additive_ratio``, per
    ``abundance_ranges`` pair for ``sparse`` (which folds ``index`` into its seed), or the single
    ``("file", "file", None, 0)`` when the abundance file supplied nRNA values.
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
    if mode == "sparse":
        if cfg.nrna.abundance_ranges is None:
            raise ValueError("nrna.abundance_ranges is required for mode='sparse'")
        for i, abundance_range in enumerate(cfg.nrna.abundance_ranges):
            label = nrna_label_for_ratio(abundance_range, cfg.nrna.ratio_labels, i)
            pairs.append((label, mode, abundance_range, i))
        return pairs

    raise ValueError(f"Unknown nRNA simulation mode: {mode}")


# ═══════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════


def run_simulation(cfg: WholeGenomeSimConfig) -> list[dict]:
    """Run the whole simulation: load the index's transcripts and any shadows, assign the base abundances once, run
    the condition grid and write the manifest. Returns the per-condition manifest entries.
    """
    # checked again: a caller can set `oracle_bam` after the config is read (--no-oracle)
    _refuse_a_run_that_keeps_no_reads(cfg)
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = Path(cfg.genome)
    gtf_path = Path(cfg.gtf)

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome not found: {genome_path}")
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    index_dir = ensure_index(cfg)
    transcripts = load_transcripts_from_index(index_dir)
    if not transcripts:
        raise RuntimeError(f"No transcripts loaded from index {index_dir}")
    index_record = _index_provenance(index_dir, gtf_path)
    n_shadow = 0
    if cfg.shadow_gtf:
        if not Path(cfg.shadow_gtf).exists():
            raise FileNotFoundError(f"shadow GTF not found: {cfg.shadow_gtf}")
        n_before = len(transcripts)
        transcripts = merge_shadow_transcripts(transcripts, cfg.shadow_gtf)
        n_shadow = len(transcripts) - n_before

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
        logger.info("Abundance file provided nascent abundances — skipping the nascent sweep")

    base_abundances = [(t.abundance or 0.0, t.nrna_abundance) for t in transcripts]

    sim = cfg.simulation

    gdna_pairs: list[tuple[str, float]] = []
    for i, rate in enumerate(cfg.gdna.rates):
        label = gdna_label_for_rate(rate, cfg.gdna.rate_labels, i)
        gdna_pairs.append((label, rate))

    # Which references are genomic is an input, never inferred: an annotated reference can be an RNA-only spike-in.
    genomic_refs = cfg.gdna.genomic_refs
    if genomic_refs is None:
        if any(rate > 0 for rate in cfg.gdna.rates):
            raise ValueError(
                "gdna.genomic_refs must list the references genomic DNA is drawn from whenever any "
                "gdna.rates entry is non-zero. Every reference is genomic or RNA-only and the "
                "classification is an input, not an inference from the annotation."
            )
        genomic_refs = []

    # The strand-overdispersion axis: one condition per value, or one at `strand_overdispersion`.
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
        emit_fastq=cfg.emit_fastq,
    )

    manifest_path = write_manifest(
        outdir, cfg, conditions, index=index_record, n_shadow_transcripts=n_shadow
    )
    logger.info("Wrote manifest to %s", manifest_path)
    return conditions
