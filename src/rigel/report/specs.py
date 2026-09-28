"""Vega-Lite spec builders for the report's charts: the fragment-length distributions and the
per-reference gDNA-density track.

Charts are emitted as plain Vega-Lite JSON specs (no colours baked in — the
front-end driver injects the theme palette + axis colours at embed time and
re-embeds on theme change). Densities are pre-computed here so the same numbers
back the chart, its tooltip, and the summary table.

The composition, pool and strand components are native HTML/SVG; every chart is built here and
embedded by the same driver.
"""

from __future__ import annotations

import numpy as np
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


def _bin_track(track: pd.DataFrame, bins_per_ref: int = 160) -> list[dict]:
    """Length-weighted per-bin gDNA density, binned within each reference.

    Binning happens here (Python) so a genome with millions of region rows still
    emits only ~bins_per_ref points per chromosome to the chart. One ``groupby``
    over the ref column keeps this O(regions) rather than O(refs × regions).
    """
    rows: list[dict] = []
    length = (track["end"] - track["start"]).clip(lower=1).to_numpy(dtype="float64")
    dens = track["gdna_density"].to_numpy(dtype="float64")
    starts = track["start"].to_numpy(dtype="float64")
    for ref, idx in track.groupby("ref", observed=True, sort=False).indices.items():
        s, ln, d = starts[idx], length[idx], dens[idx]
        span = float(s.max() + ln[s.argmax()])
        if span <= 0:
            continue
        width = span / bins_per_ref
        bin_idx = np.minimum((s / width).astype(int), bins_per_ref - 1)
        num = np.bincount(bin_idx, weights=d * ln, minlength=bins_per_ref)
        den = np.bincount(bin_idx, weights=ln, minlength=bins_per_ref)
        for b in np.nonzero(den > 0)[0]:
            rows.append(
                {
                    "ref": str(ref),
                    "pos_mb": round((b + 0.5) * width / 1e6, 4),
                    "density": float(num[b] / den[b]),
                }
            )
    return rows


def genome_track_spec(track: pd.DataFrame | None, top_n: int = 24) -> dict | None:
    """Per-reference gDNA-density track for the top references by gDNA mass.

    References are ranked by total gDNA mass (data-driven, organism-agnostic — no
    name-based filtering) and only the top ``top_n`` are faceted, so hundreds of
    small alt/haplotype contigs don't flood the panel (the full list is in the
    reference table). The y-axis is **log** and **independent per facet**, so a
    high-density outlier like chrM cannot flatten the other references' scales.
    """
    if track is None or track.empty:
        return None
    mass = track.groupby("ref", observed=True)["gdna_mass"].sum().sort_values(ascending=False)
    top = [str(r) for r in mass.index[:top_n]]
    sub = track[track["ref"].astype(str).isin(top)]
    # log axis needs strictly positive density
    rows = [r for r in _bin_track(sub) if r["density"] > 0]
    if not rows:
        return None
    n_shown = len({r["ref"] for r in rows})
    return {
        "$schema": _SCHEMA,
        "data": {"values": rows},
        "columns": 1 if n_shown <= 3 else 2,
        "facet": {
            "field": "ref",
            "type": "nominal",
            "title": None,
            "sort": top,  # biggest-mass references first
            "header": {"labelAnchor": "start", "labelFontWeight": "bold"},
        },
        "spec": {
            "width": "container",
            "height": 58,
            "mark": {
                "type": "line",
                "strokeWidth": 1.5,
                "interpolate": "monotone",
                "point": {"size": 6, "filled": True},
            },
            "encoding": {
                "x": {"field": "pos_mb", "type": "quantitative", "title": "position (Mb)"},
                "y": {
                    "field": "density",
                    "type": "quantitative",
                    "title": "gDNA density (log)",
                    "scale": {"type": "log"},
                },
                "tooltip": [
                    {"field": "ref", "title": "ref"},
                    {"field": "pos_mb", "title": "position (Mb)"},
                    {"field": "density", "title": "gDNA density", "format": ".3~g"},
                ],
            },
        },
        # Independent scales per facet — no reference (e.g. chrM) blows out another.
        "resolve": {"scale": {"x": "independent", "y": "independent"}},
    }


def build_charts(sub) -> dict:
    """All Vega-Lite charts for the report, keyed by container id (``vega-<key>``)."""
    charts = build_fl_specs(getattr(sub, "fragment_lengths", None))
    genome = genome_track_spec(getattr(sub, "calibration_track", None))
    if genome is not None:
        charts["genome"] = genome
    return charts
