"""``build_report`` — turn a ``rigel quant`` output directory into an HTML report."""

from __future__ import annotations

import logging
from pathlib import Path

from .html import render_html
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
    for w in sub.warnings:
        logger.warning("[report] %s", w)

    model = build_view_model(sub)
    charts = build_charts(sub)

    if title is None:
        title = f"Rigel QC · {model['meta']['sample']}"
    out_path = Path(out_path) if out_path is not None else Path(output_dir) / "report.html"

    try:
        import vl_convert  # noqa: F401
    except ImportError:
        logger.warning(
            "[report] vl-convert-python is not installed — the fragment-length "
            "charts will be omitted. Install with: pip install 'rigel-rnaseq[report]'"
        )

    html = render_html(model, charts, title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    logger.info("[report] wrote %s (%.0f KB)", out_path, out_path.stat().st_size / 1024)
    return out_path
