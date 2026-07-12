"""Whole-genome simulator configuration dataclasses.

The data layer for the simulator: per-run + sweep parameters consumed by the engine
(:mod:`wgs_engine`) and the suite frontend (:mod:`whole_genome`). Kept dependency-light
(only the capture config) so both can import it without cycles.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .capture import CaptureConfig, CaptureScenario


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

    Two sweep modes are supported:

    - ``additive_ratio`` (via ``ratios``): each entry adds nascent RNA at a
      fixed ratio of mature RNA, ``nrna_abundance = mrna_abundance * ratio``.
    - ``random_fraction`` (via ``ratio_ranges`` + ``eligible_fraction``): each
      entry draws per-transcript ratios uniformly from a ``(lo, hi)`` range,
      applied to a random ``eligible_fraction`` of transcripts (the rest get
      zero nascent).

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
    strand_specificities: list[float] = field(default_factory=lambda: [1.0])

    oracle_bam: bool = True
    verbose: bool = True
