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
        Added weight per base of overlap with the best-binding probe, one whose GC content is ``gc_optimum``:
        with ``off_target_weight=1`` and ``binding_per_base=10``, a full 120 bp overlap with it weighs 1201. A
        probe binds at its own strength (``gc_optimum``, ``gc_half_distance``) times this. A fragment binds the
        one contiguous probe part with the largest strength times overlap: overlapping probes do not stack, and
        neither do the parts of a probe split across a splice junction, so a transcript holding the junction
        sees the probe whole, while gDNA, a nascent span and an isoform without the junction bind the better part.
    min_overlap : int
        Minimum overlap, in bases, before a probe contributes weight.
    gc_optimum : float
        The GC fraction of the best-binding probe, where its strength is 1. 0.55 is the middle of the 50-60 % bait
        GC class that received the highest coverage in solution hybrid selection with long RNA baits, against
        the least at 30-40 % and 70-80 % (Gnirke et al. 2009, Nat Biotechnol 27:182, Results, "Effects of base
        composition", which attributes it mostly to hybridisation). A real captured exome gDNA library (VCaP), its
        copy number divided out, yields most at the same GC: the form fitted to it peaks at 0.54 to 0.57.
    gc_half_distance : float
        How far a probe's GC fraction lies from ``gc_optimum`` where it binds at half the strength: the strength
        is ``2 ** -((gc - gc_optimum) / gc_half_distance) ** 2``, one at the optimum, falling at both GC extremes.
        0.21 is measured: that VCaP library's yield falls to half its peak near GC 0.31 and 0.76, and a simulated
        library read the same way (yield per target against the target's GC, copy number divided out) matches
        its yield against GC most closely at 0.21 of the widths 0.17 to 0.23. The real fall is flatter at the top
        and steeper at the extremes than this symmetric form, so the match holds to about 0.1 in log2. ``inf`` is
        the form's own limit, where every probe binds alike whatever its sequence.
    """

    probes: str | None = None
    probe_format: str = "auto"
    off_target_weight: float = 1.0
    binding_per_base: float = 10.0
    min_overlap: int = 1
    gc_optimum: float = 0.55
    gc_half_distance: float = 0.21
