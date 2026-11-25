"""The read engine's parameters (:mod:`wgs_engine`): the read model and the gDNA draw. A whole-genome run's YAML
config reads them as its ``simulation`` and ``gdna`` sections, sweep axes included."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..config import CONSTANTS
from .reads import GDNAConfig, ReadSimConfig


@dataclass
class SimulationParams:
    """Read simulation parameters."""

    n_rna_fragments: int = 1_000_000
    #: Fixed total fragment budget: `gdna.rates` then decides only the RNA/gDNA split, which reaches high gDNA
    #: fractions at a bounded depth. `None` fixes the RNA depth and adds gDNA to it.
    n_total_fragments: int | None = None
    sim_seed: int = CONSTANTS.simulator.seed
    # the read model's defaults are the scenario API's (`ReadSimConfig`), one set for both paths
    frag_mean: float = ReadSimConfig.frag_mean
    frag_std: float = ReadSimConfig.frag_std
    frag_min: int = ReadSimConfig.frag_min
    frag_max: int = ReadSimConfig.frag_max
    read_length: int = ReadSimConfig.read_length
    error_rate: float = ReadSimConfig.error_rate
    n_workers: int = 1


@dataclass
class GDNASimConfig:
    """gDNA contamination configuration."""

    rates: list[float] = field(default_factory=lambda: [0.0])
    rate_labels: list[str] | None = None
    #: The references genomic DNA is drawn from: an explicit input, never inferred from the annotation.
    #: Every other reference is RNA-only. A config that asks for gDNA with this ``None`` is refused.
    genomic_refs: list[str] | None = None
    # the gDNA length law's defaults are the scenario API's (`GDNAConfig`)
    frag_mean: float = GDNAConfig.frag_mean
    frag_std: float = GDNAConfig.frag_std
    frag_min: int = GDNAConfig.frag_min
    frag_max: int = GDNAConfig.frag_max
    #: The per-condition gDNA strand intra-class correlation in [0, 1): region-to-region overdispersion of the
    #: sense/antisense split around ½, 0 being an exact Binomial 50/50.
    strand_overdispersion: float = 0.0
    #: The sweep axis over strand overdispersion, one condition per value; ``None`` runs one condition at
    #: ``strand_overdispersion``.
    strand_overdispersions: list[float] | None = None
    strand_overdispersion_labels: list[str] | None = None

    def __post_init__(self) -> None:
        # the engine draws each region's strand rate from Beta(a, a), a from the correlation, so 1 or more has no law
        for od in [self.strand_overdispersion, *(self.strand_overdispersions or [])]:
            if not (0.0 <= float(od) < 1.0):
                raise ValueError(f"gdna strand_overdispersion must be in [0, 1); got {od!r}")
