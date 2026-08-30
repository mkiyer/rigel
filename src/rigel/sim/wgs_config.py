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
    #: ⭐ Fixed TOTAL budget. When set, `gdna.rates` decides only the SPLIT between RNA and gDNA
    #: (`rigel.sim.orchestrator.resolve_depths`) rather than adding gDNA on top of a fixed RNA depth
    #: — the only way to reach the high-gDNA end of the spectrum at a simulatable depth. `None`
    #: keeps the legacy behaviour, so every pre-existing config is unaffected.
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
    #: ⭐ The references genomic DNA is drawn from — an EXPLICIT input, never inferred from the
    #: annotation. Every reference in the FASTA is genomic or RNA-only; naming the genomic ones
    #: classifies both, since the complement is exactly the RNA-only set. ``None`` is not a default
    #: for "guess": a config that asks for gDNA without stating this is rejected.
    #: ⚠ At least TWO real genomic references are wanted, so the gDNA deposit path is exercised over
    #: a non-trivial reference-id space; see `tests/test_sim_genomic_refs.py`.
    genomic_refs: list[str] | None = None
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
    """Nascent RNA sweep configuration.

    THREE sweep modes are supported, and ⭐ ``sparse`` is the one the benchmark panels use:

    - ⭐⭐ ``sparse`` (via ``abundance_ranges`` + ``on_fraction``): nascent RNA is
      ABSENT from most gene spans and present in a minority (owner, 2026-08-22).
      Each nascent ENTITY is switched on with probability ``on_fraction`` and, if
      on, given an ABSOLUTE molecular abundance drawn LOG-UNIFORMLY from a
      ``(lo, hi)`` range — INDEPENDENT of the mature level, so ``nascent >
      mature`` is a real case the tool must survive. The fragment share is
      EMERGENT rather than requested: it is priceable in advance with
      ``whole_genome.expected_rna_weights`` and the orchestrator records the
      realised molar ratio per condition. ⛔ ``on_fraction`` and the range TOGETHER
      set the nascent burden, so neither is readable on its own.
    - ``additive_ratio`` (via ``ratios``): each entry adds nascent RNA at a
      fixed ratio of mature RNA, ``nrna_abundance = mrna_abundance * ratio``.
    - ``fragment_share`` (via ``shares``): each entry states the nascent share
      of RNA **FRAGMENTS** in the uncaptured library, and the simulator SOLVES for
      the common molecular scale that produces it. ⭐ Why it exists at all is a
      factor of four: a nascent ENTITY spans a whole gene (mean 40,667 bp on the
      ladder's index) while a mature transcript is spliced (mean 1,708 bp), so a
      molecular ratio of 0.25 puts **86 %** of RNA fragments in nascent RNA, not
      25 % (measured 2026-08-19). The molecular ratio that gives a 20 % fragment
      share is **0.0100** — which is also the biologically sensible pre-mRNA:mRNA
      molar ratio. ⚠ The scale is solved on UNCAPTURED effective lengths, i.e. it
      fixes the LIBRARY's molecular composition; the realised fragment share then
      differs under capture, which is the physics (capture acts after the
      molecules exist).

    ⛔⛔ **BOTH RATIO MODES PUT NASCENT RNA ON EVERY EXPRESSED MULTI-EXON SPAN AT ONE
    LEVEL, AND THAT IS WHY THE PANELS LEFT THEM** (the owner's NASCENT SCOPE ruling,
    2026-08-22). Nascent mass tracks mature abundance there and can never exceed it,
    and no intron is ever exactly nascent-free — so a tool developed against them is
    designed AROUND nascent RNA rather than robust TO it. ⚠ They are not dead: the
    fl-gap side panels still use ``fragment_share``, so the two nascent models coexist
    across the panels on disk, and a config states which one it uses.

    When abundances come from a file that already contains explicit nRNA
    data, the configured sweep is ignored and the file's nRNA values are
    used as-is in a single nRNA condition.
    """

    mode: str = "additive_ratio"
    ratios: list[float] = field(default_factory=lambda: [0.0])
    #: per-condition (lo, hi) for the LOG-uniform ABSOLUTE nascent abundance, for ``mode="sparse"``
    abundance_ranges: list[tuple[float, float]] | None = None
    #: nascent share of RNA FRAGMENTS in the uncaptured library, for ``mode="fragment_share"``
    shares: list[float] | None = None
    #: labels the CONDITION whatever the quantity swept (ratios, shares or abundance ranges)
    ratio_labels: list[str] | None = None
    #: fraction of eligible gene spans that carry nascent RNA at all, for ``mode="sparse"``
    on_fraction: float = 1.0
    seed: int = 42


@dataclass
class WholeGenomeSimConfig:
    """Top-level simulation configuration."""

    genome: str = ""
    gtf: str = ""
    #: ⭐ The rigel index the simulation's TRANSCRIPTOME comes from (owner, 2026-08-19). The simulator
    #: takes the index's transcript list — annotated transcripts PLUS the synthetic nascent-RNA entities
    #: `rigel index` creates (`index.create_nrna_transcripts`, TSS/TES clustered within
    #: `NRNA_MERGE_TOLERANCE`) — so what it simulates is exactly what `rigel quant` sees. ``None`` ⇒ an
    #: index is built from ``genome`` + ``gtf`` into ``<outdir>/rigel_index`` on first use.
    index: str | None = None
    #: ⭐ SHADOW transcripts (owner design, 2026-08-29): a supplemental GTF the SIMULATOR draws fragments
    #: from but the INDEX never sees — unannotated transcription, simulated. Their ids must be unknown to
    #: the index (a shadow that is annotated is not a shadow — refused); they take abundances from the
    #: same abundance file, carry no nascent, and are tagged like any RNA fragment (``{t_id}:…``, origin
    #: ``mrna``), so the oracle split and the certified truth see them as RNA the tool cannot know about.
    shadow_gtf: str | None = None
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
