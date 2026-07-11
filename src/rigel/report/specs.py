"""Vega-Lite spec builders for the fragment-length distribution charts.

Charts are emitted as plain Vega-Lite JSON specs (no colours baked in — the
front-end driver injects the theme palette + axis colours at embed time and
re-embeds on theme change). Densities are pre-computed here so the same numbers
back the chart, its tooltip, and the summary table.

Only the fragment-length section uses Vega-Lite in v1; the composition / pool /
strand components are native HTML/SVG. The genome-track and capture-KDE charts
(Phase 2) will reuse this same spec + driver pipeline.
"""

from __future__ import annotations

import pandas as pd

_SCHEMA = "https://vega.github.io/schema/vega-lite/v5.json"

# Overlay compares the deconvolved pools; per-category small multiples show all.
_OVERLAY_CATEGORIES = ["rna", "gdna"]


def _binned_density(fl_df: pd.DataFrame, categories, binsize: int = 5) -> list[dict]:
    """Area-normalized density per category so pools of different N overlay fairly."""
    rows: list[dict] = []
    for cat in categories:
        sub = fl_df[fl_df["category"] == cat]
        if sub.empty:
            continue
        binned = sub.assign(bin=(sub["length"] // binsize) * binsize)
        agg = binned.groupby("bin")["count"].sum()
        total = float(agg.sum())
        if total <= 0:
            continue
        for b, c in agg.items():
            rows.append(
                {
                    "category": cat,
                    "length": int(b),
                    "density": float(c) / total / binsize,
                }
            )
    return rows


def fl_overlay_spec(fl_df: pd.DataFrame) -> dict | None:
    """Density overlay of the gDNA vs RNA fragment-length pools."""
    cats = [c for c in _OVERLAY_CATEGORIES if (fl_df["category"] == c).any()]
    if not cats:
        return None
    rows = _binned_density(fl_df, cats)
    if not rows:
        return None
    return {
        "$schema": _SCHEMA,
        "width": "container",
        "height": 250,
        "autosize": {"type": "fit", "contains": "padding"},
        "data": {"values": rows},
        "mark": {
            "type": "area",
            "line": {"strokeWidth": 2},
            "opacity": 0.22,
            "interpolate": "monotone",
        },
        "encoding": {
            "x": {
                "field": "length",
                "type": "quantitative",
                "title": "Fragment length (bp)",
                "scale": {"zero": False},
            },
            "y": {
                "field": "density",
                "type": "quantitative",
                "title": "Density",
                "stack": None,
                "axis": {"format": ".2~e"},
            },
            "color": {
                "field": "category",
                "type": "nominal",
                "title": None,
                "sort": _OVERLAY_CATEGORIES,
                "legend": {"orient": "top-right"},
            },
            "tooltip": [
                {"field": "category", "title": "pool"},
                {"field": "length", "title": "length (bp)"},
                {"field": "density", "title": "density", "format": ".3~e"},
            ],
        },
        "usermeta": {"embedOptions": {"actions": True}},
    }


def fl_small_multiples_spec(fl_df: pd.DataFrame, exclude=("global",)) -> dict | None:
    """Per-category histogram small multiples (independent y per facet)."""
    cats = [c for c in fl_df["category"].unique() if c not in exclude]
    if not cats:
        return None
    sub = fl_df[fl_df["category"].isin(cats)].copy()
    rows = [
        {"category": str(c), "length": int(x), "count": float(n)}
        for c, x, n in zip(sub["category"], sub["length"], sub["count"])
    ]
    return {
        "$schema": _SCHEMA,
        "data": {"values": rows},
        "columns": 3,
        "facet": {
            "field": "category",
            "type": "nominal",
            "title": None,
            "header": {"labelFontWeight": "bold"},
            "sort": cats,
        },
        "spec": {
            "width": 180,
            "height": 90,
            "mark": {"type": "bar", "binSpacing": 0},
            "encoding": {
                "x": {
                    "field": "length",
                    "type": "quantitative",
                    "bin": {"maxbins": 32},
                    "title": "bp",
                },
                "y": {"field": "count", "aggregate": "sum", "type": "quantitative", "title": None},
                "tooltip": [
                    {"field": "category", "title": "category"},
                    {"field": "length", "title": "length (bp)"},
                    {"field": "count", "aggregate": "sum", "title": "count"},
                ],
            },
        },
        "resolve": {"scale": {"y": "independent", "x": "independent"}},
    }


def build_fl_specs(fl_df: pd.DataFrame | None) -> dict:
    """Return the fragment-length Vega-Lite specs (``{}`` when no histogram data)."""
    if fl_df is None or fl_df.empty:
        return {}
    specs: dict = {}
    overlay = fl_overlay_spec(fl_df)
    if overlay is not None:
        specs["overlay"] = overlay
    small = fl_small_multiples_spec(fl_df)
    if small is not None:
        specs["small_multiples"] = small
    return specs
