"""The report pipeline in one function, :func:`build_report`: load the substrate (:mod:`.substrate`), build the
view model (:mod:`.model`) and the chart specs (:mod:`.specs`), and write one HTML file (:mod:`.html`).

The output directory must be one this version of ``rigel quant`` wrote (:func:`~.substrate.load_substrate` refuses
any other). Without ``vl-convert-python`` every chart is omitted, with a warning naming the extra to install.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .html import render_html, vega_bundle
from .model import build_view_model
from .specs import build_charts
from .substrate import load_substrate

logger = logging.getLogger(__name__)


def build_report(
    output_dir: str | Path,
    out_path: str | Path | None = None,
    title: str | None = None,
) -> Path:
    """Build a self-contained HTML QC report from a ``rigel quant`` output directory.

    Parameters
    ----------
    output_dir
        A ``rigel quant`` ``--output-dir`` (must contain ``summary.json``).
    out_path
        Destination HTML path. Defaults to ``<output_dir>/report.html``.
    title
        Document title. Defaults to ``"Rigel QC · <sample>"``.

    Returns
    -------
    Path
        The written HTML file.
    """
    sub = load_substrate(output_dir)

    model = build_view_model(sub)
    charts = build_charts(sub)

    if title is None:
        title = f"Rigel QC · {model['meta']['sample']}"
    out_path = Path(out_path) if out_path is not None else Path(output_dir) / "report.html"

    if vega_bundle() is None:
        logger.warning(
            "[report] vl-convert-python is not installed — every chart will be "
            "omitted. Install with: pip install 'rigel-rnaseq[report]'"
        )

    html = render_html(model, charts, title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    logger.info("[report] wrote %s (%.0f KB)", out_path, out_path.stat().st_size / 1024)
    return out_path
