"""Vega-Lite specs for the report's charts: the fragment-length distributions and the per-reference gDNA-density
track.

Specs carry no colours: the front end injects the theme palette at embed time and re-embeds on a theme change.
Densities are computed here, so the chart, its tooltip and the summary table show the same numbers. Bin widths,
reference caps and rounding are ``CONSTANTS.qc``'s.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import CONSTANTS

_SCHEMA = "https://vega.github.io/schema/vega-lite/v5.json"

# Overlay compares the deconvolved pools; per-category small multiples show all.
_OVERLAY_CATEGORIES = ["rna", "gdna"]

#: The label shown for each fragment-length category ``rigel quant`` writes: the anchor, the two modelled pools, and
#: the five structural pools they are built from. The summary table, the overlay legend and the facets all read it.
FL_CATEGORY_LABELS = {
    "global": "All fragments",
    "gdna": "gDNA",
    "rna": "RNA",
    "pool_gdna_intergenic": "gDNA pool · intergenic",
    "pool_gdna_intronic": "gDNA pool · intronic",
    "pool_gdna_intron_exon": "gDNA pool · intron|exon",
    "pool_gdna_intergenic_exon": "gDNA pool · intergenic|exon",
    "pool_rna_spliced": "RNA pool · spliced",
}


def fl_label(category: str) -> str:
    """A fragment-length category's label; a category the labels do not name shows as itself."""
    return FL_CATEGORY_LABELS.get(category, category)


def _binned_density(fl_df: pd.DataFrame, categories) -> list[dict]:
    """Area-normalized density per category, in bins ``CONSTANTS.qc.fl_overlay_bin_bp`` wide, so pools of
    different N overlay fairly; each row names its category by its label."""
    binsize = CONSTANTS.qc.fl_overlay_bin_bp
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
                    "category": fl_label(cat),
                    "length": int(b),
                    "density": float(c) / total / binsize,
                }
            )
    return rows


def fl_overlay_spec(fl_df: pd.DataFrame) -> dict | None:
    """Density overlay of the gDNA vs RNA fragment-length pools, named by their labels."""
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
                "sort": [fl_label(c) for c in _OVERLAY_CATEGORIES],
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


def fl_small_multiples_spec(fl_df: pd.DataFrame) -> dict | None:
    """Per-category histogram small multiples (independent y per facet), each headed by its label; the global
    histogram is left out."""
    cats = [c for c in fl_df["category"].unique() if c != "global"]
    if not cats:
        return None
    sub = fl_df[fl_df["category"].isin(cats)].copy()
    rows = [
        {"category": fl_label(str(c)), "length": int(x), "count": float(n)}
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
            "sort": [fl_label(str(c)) for c in cats],
        },
        "spec": {
            "width": 180,
            "height": 90,
            "mark": {"type": "bar", "binSpacing": 0},
            "encoding": {
                "x": {
                    "field": "length",
                    "type": "quantitative",
                    "bin": {"maxbins": CONSTANTS.qc.fl_small_multiples_bins},
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


def build_fl_specs(fl_df: pd.DataFrame) -> dict:
    """Return the fragment-length Vega-Lite specs (``{}`` when no histogram data)."""
    if fl_df.empty:
        return {}
    specs: dict = {}
    overlay = fl_overlay_spec(fl_df)
    if overlay is not None:
        specs["overlay"] = overlay
    small = fl_small_multiples_spec(fl_df)
    if small is not None:
        specs["small_multiples"] = small
    return specs


def _bin_track(track: pd.DataFrame) -> list[dict]:
    """Base-weighted gDNA density per bin, binned within each reference, so a genome with millions of regions still
    sends at most ``CONSTANTS.qc.genome_track_bins`` points per reference to the chart. A region counts in every bin
    it covers, by the bases it holds there, so one longer than a bin fills each of them. The regions are the
    partition, disjoint and in genomic order, so a bin's bases and density mass are differences of two cumulative
    sums at its edges: O(regions + bins)."""
    bins_per_ref = CONSTANTS.qc.genome_track_bins
    rows: list[dict] = []
    starts = track["start"].to_numpy(dtype="float64")
    ends = track["end"].to_numpy(dtype="float64")
    dens = track["gdna_density"].to_numpy(dtype="float64")
    for ref, idx in track.groupby("ref", observed=True, sort=False).indices.items():
        order = idx[np.argsort(starts[idx], kind="stable")]
        s, e, d = starts[order], ends[order], dens[order]
        span = float(e.max())
        if span <= 0:
            continue
        width = span / bins_per_ref
        edges = np.arange(bins_per_ref + 1, dtype=np.float64) * width
        bases = np.concatenate([[0.0], np.cumsum(e - s)])
        mass = np.concatenate([[0.0], np.cumsum(d * (e - s))])
        # below an edge: every region that ends at or before it, and the part of the next one that lies below it
        whole = np.searchsorted(e, edges, side="right")
        nxt = np.minimum(whole, s.size - 1)
        part = np.where(whole < s.size, np.clip(edges - s[nxt], 0.0, e[nxt] - s[nxt]), 0.0)
        den = np.diff(bases[whole] + part)
        num = np.diff(mass[whole] + d[nxt] * part)
        for b in np.nonzero(den > 0)[0]:
            rows.append(
                {
                    "ref": str(ref),
                    "pos_mb": round((b + 0.5) * width / 1e6, CONSTANTS.qc.position_mb_decimals),
                    "density": float(num[b] / den[b]),
                }
            )
    return rows


def genome_track_spec(track: pd.DataFrame) -> dict | None:
    """Per-reference gDNA-density track for the top ``CONSTANTS.qc.genome_track_refs`` references by gDNA count.

    Ranking by count, not by name, keeps small alt and haplotype contigs out of the panel on any organism; the
    reference table lists them all. The y-axis is log and independent per facet, so a high-density reference such
    as chrM cannot flatten the others.
    """
    if track.empty:
        return None
    gdna = track.groupby("ref", observed=True)["count_gdna"].sum().sort_values(ascending=False)
    top = [str(r) for r in gdna.index[: CONSTANTS.qc.genome_track_refs]]
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
            "sort": top,  # the references holding the most gDNA first
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
        "resolve": {"scale": {"x": "independent", "y": "independent"}},
    }


def build_charts(sub) -> dict:
    """All Vega-Lite charts for the report, keyed by container id (``vega-<key>``)."""
    charts = build_fl_specs(sub.fragment_lengths)
    genome = genome_track_spec(sub.calibration_track)
    if genome is not None:
        charts["genome"] = genome
    return charts
