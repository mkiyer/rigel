"""Single-condition read-simulation configuration: the dataclasses :class:`scenario.Scenario` reads and translates
into the engine's ``SimulationParams`` / ``GDNASimConfig``."""

from dataclasses import dataclass

__all__ = ["ReadSimConfig", "GDNAConfig"]


@dataclass
class ReadSimConfig:
    """Configuration for a single-condition read simulation.

    Attributes
    ----------
    frag_mean, frag_std, frag_min, frag_max : fragment-length normal distribution (bp).
    read_length : read length (bp); R1 and R2 share it.
    error_rate : per-base substitution error rate (0.0 = no errors).
    strand_specificity : probability an RNA fragment preserves correct read orientation
        (1.0 = perfectly stranded, 0.5 = no strand information); implemented as an R1↔R2 swap
        with probability ``1 − strand_specificity`` per fragment.
    r1_sense : the orientation the protocol targets — ``False`` (default) emits R1-antisense (dUTP),
        ``True`` emits R1-sense. This is a direction and ``strand_specificity`` is a fidelity about it;
        the two are independent.
    seed : random seed for reproducibility.
    """

    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0
    strand_specificity: float = 1.0
    r1_sense: bool = False
    seed: int = 42


@dataclass
class GDNAConfig:
    """Configuration for genomic DNA contamination.

    gDNA fragments are drawn from both genomic strands with their own fragment-length distribution.
    ``abundance`` is on the transcript-abundance scale: the gDNA fraction is ``abundance × genome_eff_len``
    over the total abundance × effective-length weight
    (:meth:`rigel.sim.wgs_engine.WholeGenomeSimulator.rna_gdna_split`).

    Attributes
    ----------
    abundance : relative gDNA abundance (transcript-abundance scale).
    frag_mean, frag_std, frag_min, frag_max : gDNA fragment-length distribution (bp).
    strand_overdispersion : float
        The intra-class correlation of the per-region sense/antisense split, in ``[0, 1)``, in the
        calibrator's units. ``0`` is an exact Binomial 50/50; larger is more region-to-region skew.
    """

    abundance: float = 10.0
    frag_mean: float = 350.0
    frag_std: float = 100.0
    frag_min: int = 100
    frag_max: int = 1000
    strand_overdispersion: float = 0.0

    def __post_init__(self) -> None:
        if not (0.0 <= float(self.strand_overdispersion) < 1.0):
            raise ValueError(
                "GDNAConfig.strand_overdispersion must be in [0, 1); "
                f"got {self.strand_overdispersion!r}."
            )
