"""rigel.sim.reads — read-simulation configuration dataclasses.

The read-generation *engine* is :class:`whole_genome.WholeGenomeSimulator` (one fast, vectorized,
parallel engine; see ``docs/CARRY_FORWARD.md``). This module holds the small
``ReadSimConfig`` / ``GDNAConfig`` dataclasses that :class:`scenario.Scenario` (and its tests) use
to describe a single-condition run; ``Scenario`` translates them into the engine's
``SimulationParams`` / ``GDNASimConfig``.

Read-name format encodes ground-truth origin (parsed by :mod:`read_name`):

    RNA:   {t_id}:{frag_start}-{frag_end}:{strand_char}:{index}/1
    gDNA:  gdna:{genomic_start}-{genomic_end}:{strand_char}:{index}/1
"""

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
    seed : random seed for reproducibility.
    """

    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0
    strand_specificity: float = 1.0
    seed: int = 42


@dataclass
class GDNAConfig:
    """Configuration for genomic DNA contamination.

    gDNA fragments are sampled from the genome (both strands) with an independent fragment-size
    distribution. ``abundance`` uses the same relative scale as transcript abundances — the gDNA
    fragment fraction is ``abundance × genome_eff_len`` over the total abundance × effective-length
    weight (see :meth:`whole_genome.WholeGenomeSimulator.pool_split`).

    Attributes
    ----------
    abundance : relative gDNA abundance (transcript-abundance scale).
    frag_mean, frag_std, frag_min, frag_max : gDNA fragment-length distribution (bp).
    gdna_strand_overdispersion : float or None
        gDNA strand overdispersion — the intra-class correlation of the per-region sense/antisense
        split in ``[0, 1)``, in the same units the calibrator fits. ``0`` ⇒ exact Binomial 50/50;
        larger ⇒ more region-to-region strand skew. Overrides ``strand_kappa`` when set.
    strand_kappa : float or None
        Legacy/low-level Beta concentration form (per-region rate ``Beta(kappa/2, kappa/2)`` with
        intra-class correlation ``1/(kappa + 1)``). Prefer ``gdna_strand_overdispersion``.
    """

    abundance: float = 10.0
    frag_mean: float = 350.0
    frag_std: float = 100.0
    frag_min: int = 100
    frag_max: int = 1000
    strand_kappa: float | None = None
    gdna_strand_overdispersion: float | None = None

    def __post_init__(self) -> None:
        # The clear-units knob (intra-class correlation) overrides the Beta concentration.
        if self.gdna_strand_overdispersion is not None:
            od = float(self.gdna_strand_overdispersion)
            if not (0.0 <= od < 1.0):
                raise ValueError(
                    "GDNAConfig.gdna_strand_overdispersion must be in [0, 1); "
                    f"got {self.gdna_strand_overdispersion!r}."
                )
            # od = 1/(kappa + 1) ⇒ kappa = (1 − od)/od; od == 0 ⇒ Binomial (no region partition).
            self.strand_kappa = None if od == 0.0 else (1.0 - od) / od
