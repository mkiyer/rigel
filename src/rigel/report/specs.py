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
    emits only ~bins_per_ref points per chromosome to the chart.
    """
    rows: list[dict] = []
    length = (track["end"] - track["start"]).clip(lower=1).to_numpy(dtype="float64")
    dens = track["gdna_density"].to_numpy(dtype="float64")
    starts = track["start"].to_numpy(dtype="float64")
    refs = track["ref"].astype(str).to_numpy()
    for ref in pd.unique(refs):
        m = refs == ref
        span = float(starts[m].max() + length[m][starts[m].argmax()])
        if span <= 0:
            continue
        width = span / bins_per_ref
        bin_idx = np.minimum((starts[m] / width).astype(int), bins_per_ref - 1)
        wl, wd = length[m], dens[m] * length[m]
        num = np.bincount(bin_idx, weights=wd, minlength=bins_per_ref)
        den = np.bincount(bin_idx, weights=wl, minlength=bins_per_ref)
        for b in range(bins_per_ref):
            if den[b] > 0:
                rows.append(
                    {
                        "ref": str(ref),
                        "pos_mb": round((b + 0.5) * width / 1e6, 4),
                        "density": float(num[b] / den[b]),
                    }
                )
    return rows


def genome_track_spec(track: pd.DataFrame | None) -> dict | None:
    """Whole-genome gDNA density overview — a binned area, faceted per reference."""
    if track is None or track.empty:
        return None
    rows = _bin_track(track)
    if not rows:
        return None
    n_refs = len({r["ref"] for r in rows})
    return {
        "$schema": _SCHEMA,
        "data": {"values": rows},
        "columns": 1 if n_refs <= 3 else 2,
        "facet": {
            "field": "ref",
            "type": "nominal",
            "title": None,
            "header": {"labelAnchor": "start", "labelFontWeight": "bold"},
        },
        "spec": {
            "width": "container",
            "height": 60,
            "mark": {
                "type": "area",
                "line": {"strokeWidth": 1},
                "opacity": 0.5,
                "interpolate": "monotone",
            },
            "encoding": {
                "x": {"field": "pos_mb", "type": "quantitative", "title": "position (Mb)"},
                "y": {"field": "density", "type": "quantitative", "title": "gDNA density"},
                "tooltip": [
                    {"field": "ref", "title": "ref"},
                    {"field": "pos_mb", "title": "position (Mb)"},
                    {"field": "density", "title": "gDNA density", "format": ".3~g"},
                ],
            },
        },
        "resolve": {"scale": {"x": "independent"}},
    }


_KIND_NAMES = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def capture_kde_spec(
    kde: pd.DataFrame | None, nodes: pd.DataFrame | None, capture: dict | None
) -> dict | None:
    """gDNA-density KDE: the ``P(log ρ_g)`` curve + labeled depleted/enriched modes + a node rug.

    Bimodality is the capture signal (a low off-target mode + a high on-target
    mode). No categorical verdict — the separation / enrichment numbers are shown
    and left to the analyst.
    """
    if kde is None or kde.empty:
        return None
    curve = [
        {"log_rho": float(x), "density": float(d)} for x, d in zip(kde["log_rho"], kde["density"])
    ]

    layers = [
        {
            "data": {"values": curve},
            "mark": {
                "type": "area",
                "line": {"strokeWidth": 2},
                "opacity": 0.18,
                "interpolate": "monotone",
            },
            "encoding": {
                "x": {
                    "field": "log_rho",
                    "type": "quantitative",
                    "title": "log gDNA density  ρg  (nats)",
                },
                "y": {"field": "density", "type": "quantitative", "title": "P(log ρg)"},
                "tooltip": [
                    {"field": "log_rho", "title": "log ρg", "format": ".2f"},
                    {"field": "density", "title": "density", "format": ".3~e"},
                ],
            },
        }
    ]

    # Node rug, coloured by node kind, along the bottom.
    if nodes is not None and not nodes.empty:
        rug = [
            {"log_rho": float(x), "kind": _KIND_NAMES.get(int(k), "?")}
            for x, k in zip(nodes["log_rho"], nodes["kind"])
        ]
        layers.append(
            {
                "data": {"values": rug},
                "mark": {
                    "type": "tick",
                    "opacity": 0.5,
                    "thickness": 1.5,
                    "size": 10,
                    "yOffset": 120,
                },
                "encoding": {
                    "x": {"field": "log_rho", "type": "quantitative"},
                    "color": {
                        "field": "kind",
                        "type": "nominal",
                        "title": "node",
                        "sort": ["intergenic", "intron", "exon", "boundary"],
                        "legend": {"orient": "top-right"},
                    },
                    "tooltip": [
                        {"field": "kind", "title": "node kind"},
                        {"field": "log_rho", "title": "log ρg", "format": ".2f"},
                    ],
                },
            }
        )

    # Depleted / enriched mode rules + labels.
    modes = []
    if capture:
        if capture.get("depleted_mode_log_rho") is not None:
            modes.append({"log_rho": capture["depleted_mode_log_rho"], "label": "depleted"})
        enr = capture.get("enriched_mode_log_rho")
        if enr is not None and capture.get("is_bimodal"):
            modes.append({"log_rho": enr, "label": "enriched"})
    if modes:
        layers.append(
            {
                "data": {"values": modes},
                "mark": {"type": "rule", "strokeDash": [3, 3], "opacity": 0.7},
                "encoding": {"x": {"field": "log_rho", "type": "quantitative"}},
            }
        )
        layers.append(
            {
                "data": {"values": modes},
                "mark": {"type": "text", "dy": -6, "fontWeight": "bold", "baseline": "bottom"},
                "encoding": {
                    "x": {"field": "log_rho", "type": "quantitative"},
                    "y": {"value": 0},
                    "text": {"field": "label"},
                },
            }
        )

    return {
        "$schema": _SCHEMA,
        "width": "container",
        "height": 250,
        "layer": layers,
        "autosize": {"type": "fit", "contains": "padding"},
    }


def _capture_meta(sub) -> dict | None:
    cal = getattr(sub, "summary", {}).get("calibration") or {}
    return cal.get("capture")


def build_charts(sub) -> dict:
    """All Vega-Lite charts for the report, keyed by container id (``vega-<key>``)."""
    charts = build_fl_specs(getattr(sub, "fragment_lengths", None))
    genome = genome_track_spec(getattr(sub, "calibration_track", None))
    if genome is not None:
        charts["genome"] = genome
    capture = capture_kde_spec(
        getattr(sub, "gdna_density_kde", None),
        getattr(sub, "gdna_density_nodes", None),
        _capture_meta(sub),
    )
    if capture is not None:
        charts["capture_kde"] = capture
    return charts
