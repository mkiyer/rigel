"""Hybrid-capture configuration dataclasses (data only).

``CaptureConfig`` parameterizes the capture-aware weighting model; ``CaptureScenario`` labels one
setting in a condition grid. The runtime sampler lives in :mod:`capture.sampler`.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CaptureConfig:
    """Hybrid-capture simulation configuration.

    Parameters
    ----------
    probes : str or None
        Probe file path.  Transcript-coordinate TSV and BED12 are supported.
        When unset, capture weighting is disabled and simulation remains
        uniform within each template.
    probe_format : {"auto", "transcript", "bed12"}
        Probe file format.  ``"auto"`` detects a BED12-looking row or a
        transcript-coordinate TSV/header.
    off_target_weight : float
        Baseline weight for every legal fragment start.  A positive value
        keeps off-target fragments possible.
    binding_per_base : float
        Additional weight per base of the best overlapping probe.  For example, with
        ``off_target_weight=1`` and ``binding_per_base=10``, a full 120 bp
        probe overlap has weight 1201 relative to weight 1 off target.
        Overlapping probes do not stack; the best single scaled overlap is used.
    gdna_split_penalty : float
        Multiplier applied to projected genomic/pre-mRNA blocks when a probe
        is split across exon-exon sj.  Mature RNA sees the contiguous
        transcript probe; unspliced molecules and gDNA only see separated
        genomic blocks and therefore get less binding weight.
    min_overlap : int
        Minimum overlap, in bases, required before a probe contributes weight.
    """

    probes: str | None = None
    probe_format: str = "auto"
    off_target_weight: float = 1.0
    binding_per_base: float = 10.0
    gdna_split_penalty: float = 0.2
    min_overlap: int = 1


@dataclass
class CaptureScenario:
    """One labeled hybrid-capture setting in a simulation condition grid."""

    label: str = "default"
    config: CaptureConfig = field(default_factory=CaptureConfig)
