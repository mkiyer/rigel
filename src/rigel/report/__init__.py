"""Build a self-contained HTML QC report from ``rigel quant`` outputs (:func:`build_report`).

``rigel quant`` writes the report substrate (``summary.json`` and the companion feather tables); ``rigel report``
turns it into one offline HTML page at any later time, so quant carries none of the report's dependencies. Only the
charts need the ``rigel[report]`` extra; without it the report builds with the charts left out.
"""

from .build import build_report

__all__ = ["build_report"]
