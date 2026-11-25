"""The whole-genome simulation config: the YAML `simulate_reads.py` and `panel.py` read, its dataclasses and its parser.
A helper module, not a script.

:func:`parse_yaml_config` refuses any key it does not read, since a misspelt key would otherwise run on its default,
and expands ``~`` and ``$NAME`` in every path (:func:`expand_config_path`). The engine's own parameters
(``SimulationParams``, ``GDNASimConfig``) and the capture weight model (``CaptureConfig``) are the shipped
simulator's; this module adds the sweep axes and the run's inputs around them.
"""

from __future__ import annotations

import os
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from rigel.sim.capture import CaptureConfig
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams


# ── the config ───────────────────────────────────────────────────────────────────────────────────


@dataclass
class AbundanceConfig:
    """How each transcript's base abundance is set: a log-uniform draw, or a salmon / kallisto / sim file."""

    mode: str = "random"  # "random" or "file"
    seed: int = 42
    min: float = 0.1
    max: float = 10000.0
    frac_expressed: float = 0.6
    file: str | None = None


@dataclass
class NascentConfig:
    """Nascent RNA sweep: the mode that sets the entity levels, and its sweep values. A field of the other mode is
    refused, not ignored.

    - ``sparse`` (``abundance_ranges`` + ``on_fraction``): each nascent entity is on with probability
      ``on_fraction`` and, if on, draws an absolute abundance log-uniformly from ``(lo, hi)``, independent of the
      mature level, so ``nascent > mature`` occurs and an intron can be exactly nascent-free.
    - ``additive_ratio`` (``ratios``): ``nascent_abundance = annotated_abundance * ratio`` on every expressed multi-exon
      span, so nascent tracks mature and never exceeds it.

    An abundance file with its own nascent RNA values replaces the sweep with one condition.
    """

    mode: str = "additive_ratio"
    ratios: list[float] = field(default_factory=lambda: [0.0])
    #: per-condition (lo, hi) for the log-uniform absolute nascent abundance, for ``mode="sparse"``
    abundance_ranges: list[tuple[float, float]] | None = None
    #: labels the condition whatever the quantity swept (ratios or abundance ranges)
    ratio_labels: list[str] | None = None
    #: fraction of eligible gene spans that carry nascent RNA at all, for ``mode="sparse"``
    on_fraction: float = 1.0
    seed: int = 42


@dataclass
class CaptureScenario:
    """One labeled hybrid-capture setting on the grid's capture axis."""

    label: str
    config: CaptureConfig


@dataclass
class WholeGenomeSimConfig:
    """Top-level simulation configuration."""

    genome: str = ""
    #: The GTF the index is built from; a GTF whose sha256 differs from the index's record is refused.
    gtf: str = ""
    #: The rigel index whose transcript list (annotated transcripts plus its synthetic nascent entities) is
    #: simulated, so the simulator and `rigel quant` see the same rows. ``None`` builds one from ``genome`` +
    #: ``gtf`` into ``<outdir>/rigel_index`` on first use.
    index: str | None = None
    #: Shadow transcripts: a GTF of unannotated transcription the simulator draws from and the index never sees.
    #: Their ids must be unknown to the index. They take abundances like any transcript, carry no nascent, and are
    #: named like any RNA fragment, so the truth sees RNA the tool cannot know about.
    shadow_gtf: str | None = None
    outdir: str = "sim_output"

    simulation: SimulationParams = field(default_factory=SimulationParams)
    abundance: AbundanceConfig = field(default_factory=AbundanceConfig)
    gdna: GDNASimConfig = field(default_factory=GDNASimConfig)
    nascent: NascentConfig = field(default_factory=NascentConfig)
    capture: CaptureConfig = field(default_factory=CaptureConfig)
    capture_configs: list[CaptureScenario] = field(default_factory=list)
    strand_specificities: list[float] = field(default_factory=lambda: [1.0])

    #: The two read outputs; at least one must be kept, since a run keeping neither keeps no reads.
    oracle_bam: bool = True
    emit_fastq: bool = True
    verbose: bool = True


# ── paths ────────────────────────────────────────────────────────────────────────────────────────

#: A variable reference in a path, ``$NAME`` or ``${NAME}``, as :func:`os.path.expandvars` reads it (an unclosed
#: ``${`` included), so the check sees every reference the expansion would leave in place.
_PATH_VARIABLE = re.compile(r"\$(\w+|\{[^}]*\}?)", re.ASCII)

#: The top-level keys that name a path.
_PATH_KEYS = ("genome", "gtf", "shadow_gtf", "index", "outdir")


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


# ── the plain sections ───────────────────────────────────────────────────────────────────────────


def _optional(cast: Callable) -> Callable:
    """*cast*, passing ``None`` through."""
    return lambda value: None if value is None else cast(value)


def _list_of(cast: Callable) -> Callable:
    return lambda values: [cast(value) for value in values]


def _as_is(value):
    return value


#: Each plain section's fields, in the order they are read, with the cast its YAML value takes. A field the YAML
#: does not name keeps its default.
_SIMULATION_FIELDS: dict[str, Callable] = {
    "n_rna_fragments": int,
    "n_total_fragments": _optional(int),
    "sim_seed": int,
    "frag_mean": float,
    "frag_std": float,
    "frag_min": int,
    "frag_max": int,
    "read_length": int,
    "error_rate": float,
    "n_workers": int,
}
_ABUNDANCE_FIELDS: dict[str, Callable] = {
    "mode": _as_is,
    "seed": int,
    "min": float,
    "max": float,
    "frac_expressed": float,
    "file": lambda value: expand_config_path(value, "abundance.file"),
}
_GDNA_FIELDS: dict[str, Callable] = {
    "rates": _list_of(float),
    "rate_labels": _as_is,
    "genomic_refs": _optional(_list_of(str)),
    "frag_mean": float,
    "frag_std": float,
    "frag_min": int,
    "frag_max": int,
    "strand_overdispersion": float,
    "strand_overdispersions": _optional(_list_of(float)),
    "strand_overdispersion_labels": _as_is,
}


def _read_fields(raw: dict, section, casts: dict[str, Callable]) -> None:
    for name, cast in casts.items():
        if name in raw:
            setattr(section, name, cast(raw[name]))


# ── the nascent section ──────────────────────────────────────────────────────────────────────────

#: The fields each nascent mode reads. A field of the other mode is refused: ``abundance_ranges`` + ``on_fraction``
#: with ``mode:`` omitted would otherwise run ``additive_ratio`` at ratio 0.0, a nascent-free run named like the
#: configured one.
_NASCENT_MODE_FIELDS = {
    "additive_ratio": ("ratios",),
    "sparse": ("abundance_ranges", "on_fraction"),
}


def _read_nascent(raw: dict, nascent: NascentConfig) -> None:
    nascent.mode = str(raw.get("mode", nascent.mode))
    if nascent.mode not in _NASCENT_MODE_FIELDS:
        raise ValueError("nascent.mode must be 'additive_ratio' or 'sparse'")
    if raw.get("ratios") is not None:
        nascent.ratios = [float(r) for r in raw["ratios"]]
    ranges = raw.get("abundance_ranges")
    if ranges is not None:
        for pair in ranges:
            # a range must be a pair: truncating a longer entry would run another scenario
            if not hasattr(pair, "__len__") or len(pair) != 2:
                raise ValueError(f"nascent.abundance_ranges entries must be [lo, hi] pairs; got {pair!r}")
        nascent.abundance_ranges = [(float(pair[0]), float(pair[1])) for pair in ranges]
    nascent.ratio_labels = raw.get("ratio_labels", nascent.ratio_labels)
    nascent.on_fraction = float(raw.get("on_fraction", nascent.on_fraction))
    nascent.seed = int(raw.get("seed", nascent.seed))
    _check_nascent(raw, nascent)


def _check_nascent(raw: dict, nascent: NascentConfig) -> None:
    """Refuse a nascent section its mode cannot run as written."""
    if nascent.mode == "sparse" and nascent.abundance_ranges is None:
        raise ValueError("nascent.abundance_ranges is required for mode='sparse'")
    own = _NASCENT_MODE_FIELDS[nascent.mode]
    ignored = sorted(
        f
        for mode, fields in _NASCENT_MODE_FIELDS.items()
        if mode != nascent.mode
        for f in fields
        if f in raw and f not in own
    )
    if ignored:
        raise ValueError(
            f"nascent.mode={nascent.mode!r} cannot read {ignored} — those fields belong to another mode "
            f"and would be silently ignored. Remove them, or set the mode they belong to."
        )
    if not 0.0 <= nascent.on_fraction <= 1.0:
        raise ValueError("nascent.on_fraction must be between 0 and 1")
    if nascent.mode == "additive_ratio":
        n_conditions = len(nascent.ratios)
    else:
        for lo, hi in nascent.abundance_ranges:
            # a log-uniform draw has no zero end: express "no nascent" with on_fraction = 0
            if lo <= 0 or hi < lo:
                raise ValueError("nascent.abundance_ranges entries must satisfy 0 < min <= max")
        n_conditions = len(nascent.abundance_ranges)
    if nascent.ratio_labels is not None and len(nascent.ratio_labels) != n_conditions:
        raise ValueError("nascent.ratio_labels must match the number of nascent RNA scenarios")


# ── the capture section ──────────────────────────────────────────────────────────────────────────

#: Every key a capture mapping may carry: the scenario's own (``label``, ``enabled``) and the ``CaptureConfig``
#: fields.
_CAPTURE_KEYS = frozenset(
    {
        "label",
        "enabled",
        "probes",
        "probe_format",
        "off_target_weight",
        "binding_per_base",
        "min_overlap",
        "gc_optimum",
        "gc_half_distance",
    }
)


def _capture_config(
    raw: dict,
    defaults: dict | None = None,
    *,
    label: str = "capture",
    require_probes_when_enabled: bool = False,
) -> CaptureConfig:
    """A CaptureConfig from one YAML mapping over *defaults*. An unknown key is refused, since ignored it would read
    as applied. With *require_probes_when_enabled*, a mapping that says ``enabled: true`` must name its probes."""
    merged = {**(defaults or {}), **raw}
    unknown = sorted(set(merged) - _CAPTURE_KEYS)
    if unknown:
        raise ValueError(
            f"capture config '{label}': unknown key(s) {unknown}; known: {sorted(_CAPTURE_KEYS)}"
        )
    if not bool(merged.get("enabled", True)):
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
        gc_optimum=float(merged.get("gc_optimum", CaptureConfig.gc_optimum)),
        gc_half_distance=float(merged.get("gc_half_distance", CaptureConfig.gc_half_distance)),
    )


def _capture_label(value: object) -> str:
    """A capture label as text; YAML reads a bare ``on`` / ``off`` as a bool."""
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value)


def _capture_scenarios(raw: dict) -> list[CaptureScenario]:
    """``capture.configs``, each entry over the section's other keys as its defaults."""
    entries = raw.get("configs")
    if entries is None:
        return []
    if not isinstance(entries, list) or not entries:
        raise ValueError("capture.configs must be a non-empty list")
    defaults = {key: value for key, value in raw.items() if key != "configs"}
    scenarios: list[CaptureScenario] = []
    seen: set[str] = set()
    for number, entry in enumerate(entries, start=1):
        if not isinstance(entry, dict):
            raise ValueError("each capture.configs entry must be a mapping")
        label = _capture_label(entry.get("label", f"c{number}"))
        if not label:
            raise ValueError("capture config labels must be non-empty")
        if label in seen:
            raise ValueError(f"duplicate capture config label: {label}")
        seen.add(label)
        config = _capture_config(entry, defaults, label=label, require_probes_when_enabled=True)
        scenarios.append(CaptureScenario(label=label, config=config))
    return scenarios


def _read_capture(raw, cfg: WholeGenomeSimConfig) -> None:
    raw = raw or {}
    if not isinstance(raw, dict):
        raise ValueError("capture must be a YAML mapping")
    cfg.capture_configs = _capture_scenarios(raw)
    cfg.capture = cfg.capture_configs[0].config if cfg.capture_configs else _capture_config(raw)


# ── the file ─────────────────────────────────────────────────────────────────────────────────────

#: Every top-level key :func:`parse_yaml_config` reads.
_TOP_KEYS = frozenset(
    {
        *_PATH_KEYS,
        "simulation",
        "abundance",
        "nascent",
        "gdna",
        "capture",
        "strand_specificities",
        "oracle_bam",
        "emit_fastq",
        "verbose",
    }
)
#: Every key each plain section reads (``capture`` checks its own).
_SECTION_KEYS = {
    "simulation": frozenset(_SIMULATION_FIELDS),
    "abundance": frozenset(_ABUNDANCE_FIELDS),
    "nascent": frozenset({"mode", "ratios", "abundance_ranges", "ratio_labels", "on_fraction", "seed"}),
    "gdna": frozenset(_GDNA_FIELDS),
}


def _refuse_unknown_keys(raw: dict, known: frozenset[str], where: str) -> None:
    """Raise on a key the parser does not read, since a misspelt key would otherwise run on its default."""
    unknown = sorted(set(raw) - known)
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {unknown}; known: {sorted(known)}")


def refuse_a_run_that_keeps_no_reads(cfg: WholeGenomeSimConfig) -> None:
    """Raise when ``oracle_bam`` and ``emit_fastq`` are both false: no condition would keep a read, and without
    either file a condition has no resume key."""
    if not cfg.oracle_bam and not cfg.emit_fastq:
        raise ValueError(
            "oracle_bam and emit_fastq are both false (--no-oracle sets oracle_bam false): no condition "
            "would keep a read, since no oracle BAM is written and the FASTQs are deleted once the truth "
            "is, and every run would simulate every condition again. Keep oracle_bam true (the reads "
            "every instrument scores) or emit_fastq true."
        )


def parse_yaml_config(path: str | Path) -> WholeGenomeSimConfig:
    """Parse a YAML config file into a WholeGenomeSimConfig."""
    with open(path) as f:
        raw = yaml.safe_load(f)
    _refuse_unknown_keys(raw, _TOP_KEYS, "sim config")
    for section, known in _SECTION_KEYS.items():
        _refuse_unknown_keys(raw.get(section, {}), known, section)

    cfg = WholeGenomeSimConfig()
    for key in _PATH_KEYS:
        setattr(cfg, key, expand_config_path(raw.get(key, getattr(cfg, key)), key))
    _read_fields(raw.get("simulation", {}), cfg.simulation, _SIMULATION_FIELDS)
    _read_fields(raw.get("abundance", {}), cfg.abundance, _ABUNDANCE_FIELDS)
    _read_nascent(raw.get("nascent", {}), cfg.nascent)
    _read_fields(raw.get("gdna", {}), cfg.gdna, _GDNA_FIELDS)
    _read_capture(raw.get("capture", {}), cfg)
    cfg.strand_specificities = [float(s) for s in raw.get("strand_specificities", [1.0])]
    cfg.oracle_bam = bool(raw.get("oracle_bam", True))
    # False drops each condition's FASTQs once its truth is written, saving most of a panel's disk; a panel
    # without FASTQs cannot be re-aligned or given to another quantifier without re-simulating.
    cfg.emit_fastq = bool(raw.get("emit_fastq", True))
    refuse_a_run_that_keeps_no_reads(cfg)
    cfg.verbose = bool(raw.get("verbose", True))
    return cfg
