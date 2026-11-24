"""The read engine's parameters (:mod:`wgs_engine`): the read model and the gDNA draw. A whole-genome run's YAML
config reads them as its ``simulation`` and ``gdna`` sections, sweep axes included."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..config import CONSTANTS


@dataclass
class SimulationParams:
    """Read simulation parameters."""

    n_rna_fragments: int = 1_000_000
    #: Fixed total fragment budget: `gdna.rates` then decides only the RNA/gDNA split, which reaches high gDNA
    #: fractions at a bounded depth. `None` fixes the RNA depth and adds gDNA to it.
    n_total_fragments: int | None = None
    sim_seed: int = CONSTANTS.simulator.seed
    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0
    n_workers: int = 1


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
