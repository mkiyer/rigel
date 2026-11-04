"""Hybrid-capture configuration: ``CaptureConfig`` parameterises the weight model of :mod:`.sampler`."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CaptureConfig:
    """Hybrid-capture simulation configuration.

    Parameters
    ----------
    probes : str or None
        Probe file path, a transcript-coordinate TSV or BED12. Unset disables capture: starts are
        uniform within each template.
    probe_format : {"auto", "transcript", "bed12"}
        ``"auto"`` detects a BED12 row or a transcript-coordinate TSV.
    off_target_weight : float
        Baseline weight of every legal fragment start; positive keeps off-target fragments possible.
    binding_per_base : float
        Added weight per base of the fragment's best overlap with one contiguous probe part: with
        ``off_target_weight=1`` and ``binding_per_base=10``, a full 120 bp overlap weighs 1201.
        Overlapping probes do not stack, and neither do the parts of a probe split across a splice
        junction: a transcript holding the junction sees the probe whole, while gDNA, a nascent span
        and an isoform without the junction bind the better part.
    min_overlap : int
        Minimum overlap, in bases, before a probe contributes weight.
    """

    probes: str | None = None
    probe_format: str = "auto"
    off_target_weight: float = 1.0
    binding_per_base: float = 10.0
    min_overlap: int = 1
