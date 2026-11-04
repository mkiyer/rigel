"""Whole-genome simulator configuration: the per-run and sweep parameters :mod:`wgs_engine` and :mod:`whole_genome`
read. It imports only the capture config, so both can import it without a cycle."""

from __future__ import annotations

from dataclasses import dataclass, field

from .capture import CaptureConfig, CaptureScenario


@dataclass
class SimulationParams:
    """Read simulation parameters."""

    n_rna_fragments: int = 1_000_000
    #: Fixed total fragment budget: `gdna.rates` then decides only the RNA/gDNA split
    #: (`rigel.sim.orchestrator.resolve_depths`), which reaches high gDNA fractions at a bounded depth.
    #: `None` fixes the RNA depth and adds gDNA to it.
    n_total_fragments: int | None = None
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
    #: The references genomic DNA is drawn from: an explicit input, never inferred from the annotation.
    #: Every other reference is RNA-only. A config that asks for gDNA with this ``None`` is refused.
    genomic_refs: list[str] | None = None
    frag_mean: float = 350.0
    frag_std: float = 100.0
    frag_min: int = 100
    frag_max: int = 1000
    #: The per-condition gDNA strand intra-class correlation in [0, 1): region-to-region overdispersion of the
    #: sense/antisense split around ½, 0 being an exact Binomial 50/50.
    strand_overdispersion: float = 0.0
    #: The sweep axis over strand overdispersion, one condition per value; ``None`` runs one condition at
    #: ``strand_overdispersion``.
    strand_overdispersions: list[float] | None = None
    strand_overdispersion_labels: list[str] | None = None


@dataclass
class NRNAConfig:
    """Nascent RNA sweep: the mode that sets the entity levels, and its sweep values. A field of the other mode is
    refused, not ignored.

    - ``sparse`` (``abundance_ranges`` + ``on_fraction``): each nascent entity is on with probability
      ``on_fraction`` and, if on, draws an absolute abundance log-uniformly from ``(lo, hi)``, independent of the
      mature level, so ``nascent > mature`` occurs and an intron can be exactly nascent-free.
    - ``additive_ratio`` (``ratios``): ``nrna_abundance = mrna_abundance * ratio`` on every expressed multi-exon
      span, so nascent tracks mature and never exceeds it.

    An abundance file with its own nRNA values replaces the sweep with one condition.
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
    nrna: NRNAConfig = field(default_factory=NRNAConfig)
    capture: CaptureConfig = field(default_factory=CaptureConfig)
    capture_configs: list[CaptureScenario] = field(default_factory=list)
    strand_specificities: list[float] = field(default_factory=lambda: [1.0])

    #: The two read outputs; at least one must be kept, since a run keeping neither keeps no reads.
    oracle_bam: bool = True
    emit_fastq: bool = True
    verbose: bool = True
