"""rigel.report — build a self-contained HTML QC report from ``rigel quant`` outputs.

The report is a **separate, optional** step: ``rigel quant`` writes the report
*substrate* (``summary.json`` + the companion feather tables); ``rigel report``
turns that substrate into a shareable, offline HTML page at any later time. This
keeps the quant hot path dependency-lean — the report builder lives behind the
``rigel[report]`` extra.

Public entry point: :func:`build_report`.
"""

from .build import build_report

__all__ = ["build_report"]
