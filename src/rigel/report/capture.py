"""Capture-enrichment diagnostic for the report.

Thin wrapper over the core :func:`rigel.calibration.track.capture_summary` so
``summary.json`` and the report agree — here we just ask for the plottable KDE
curves (density by region count vs by gDNA mass) as well.
"""

from __future__ import annotations

from ..calibration.track import capture_summary


def capture_kde_from_track(track) -> dict | None:
    """Capture scalars + count-vs-mass KDE curves; ``None`` if uninformative."""
    return capture_summary(track, with_curve=True)
