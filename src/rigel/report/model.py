"""Transform a :class:`ReportSubstrate` into the render-ready view model: a JSON-serializable dict the front end
only lays out, so no logic runs in the browser.

Most values are raw numbers the front end formats, so tooltips, tables and charts share them. The headline captions
and a few text KPIs are finished strings, and some table cells are rounded here. Every tile cut, table cap and
rounding is ``CONSTANTS.qc``'s. The summary is read only through :class:`~rigel.report.substrate.SummaryReader`.
"""

from __future__ import annotations

from ..config import CONSTANTS
from .specs import fl_label
from .substrate import ReportSubstrate, SummaryReader

# Categorical colour slots of the palette; the front end resolves them to colours per theme.
_C = {i: f"f{i}" for i in range(1, 9)}

#: The ``gene_quant`` columns the gene table shows, in order. A column the table lacks reads as empty cells.
GENE_COLUMNS = (
    "gene_name",
    "gene_id",
    "gene_type",
    "ref",
    "strand",
    "tpm",
    "count",
    "n_transcripts",
)


def _pct(numer: float, denom: float) -> float:
    return (numer / denom) if denom else 0.0


def _capture(summary: SummaryReader) -> dict | None:
    """Calibration's capture answer: the gDNA landscape's bulk density (fragments per bp, the level every capture
    efficiency is read against) and its saturated level as a multiple of the bulk (1: nothing above it).
    ``bulk_density`` is ``None`` when no landscape was fitted; the answer is ``None`` when the summary has none."""
    if not summary.holds("calibration.gdna_bulk_density"):
        return None
    cal = summary["calibration"]
    return {
        "bulk_density": cal["gdna_bulk_density"],
        "saturation": cal["gdna_capture_saturation"],
    }


def _verdicts(summary: SummaryReader) -> list[dict]:
    """Headline QC tiles: mapping, strandedness, gDNA, usable fragments, + calibration's capture
    answer when the summary carries one."""
    tiles = [
        _mapping_tile(summary),
        _strandedness_tile(summary),
        _gdna_tile(summary),
        _usable_tile(summary),
    ]
    capture = _capture(summary)
    if capture is not None:
        tiles.append(_capture_tile(capture))
    return tiles


def _tile(k: str, icon: str, v, fmt: str, s: str, n: str) -> dict:
    """One headline tile: its label, icon, value and the value's format, status and note."""
    return {"k": k, "icon": icon, "v": v, "fmt": fmt, "s": s, "n": n}


def _grade(value: float, good: float, warning: float, *, higher_is_better: bool) -> str:
    """``good``, ``warning`` or ``critical``: *value* against the two cuts. A higher-is-better value passes a cut at
    or above it, a lower-is-better one strictly below it."""

    def passes(cut: float) -> bool:
        return value >= cut if higher_is_better else value < cut

    if passes(good):
        return "good"
    return "warning" if passes(warning) else "critical"


def _mapping_tile(summary: SummaryReader) -> dict:
    qc = CONSTANTS.qc
    aln = summary["alignment_stats"]
    total, mapped = aln["total_reads"], aln["mapped_reads"]
    rate = _pct(mapped, total)
    status = _grade(rate, qc.mapping_rate_good, qc.mapping_rate_warning, higher_is_better=True)
    return _tile(
        "Mapping rate", "check", rate, "pct", status, f"{mapped:,} of {total:,} reads mapped"
    )


def _strandedness_tile(summary: SummaryReader) -> dict:
    """The strand model's own verdict, the one calibration's strand channel reads."""
    strand = summary["strand_model"]
    if strand["strand_live"]:
        status, label = "good", "Reverse (R1−)" if not strand["read1_sense"] else "Forward (R1+)"
    else:
        status, label = "warning", "Unstranded"
    note = f"{strand['strand_specificity'] * 100:.1f}% strand-specific · {strand['protocol']}"
    return _tile("Strandedness", "arrow", label, "text", status, note)


def _gdna_tile(summary: SummaryReader) -> dict:
    qc = CONSTANTS.qc
    fraction = summary["quantification"]["fraction_gdna"]
    status = _grade(
        fraction, qc.gdna_fraction_good, qc.gdna_fraction_warning, higher_is_better=False
    )
    return _tile(
        "gDNA contamination",
        "warn",
        fraction,
        "pct",
        status,
        "EM-corrected; incl. unspliced intergenic",
    )


def _usable_tile(summary: SummaryReader) -> dict:
    """The genic fragments, the ones the quantifier can assign to RNA."""
    frag = summary["fragment_stats"]
    total, genic = frag["total"], frag["genic"]
    rate = _pct(genic, total)
    status = "good" if rate >= CONSTANTS.qc.genic_fraction_good else "warning"
    return _tile(
        "Usable fragments",
        "check",
        genic,
        "count",
        status,
        f"{rate * 100:.1f}% genic of {_si(total)}",
    )


def _capture_tile(capture: dict) -> dict:
    """Capture is descriptive: an enriched level is neither good nor bad, so the tile is neutral ("info")."""
    if capture["bulk_density"] is None:
        value, note = "None", "no gDNA landscape: lengths not corrected for capture"
    elif capture["saturation"] > 1.0:
        value = f"{capture['saturation']:,.0f}×"
        note = f"saturated gDNA level over the bulk ({capture['bulk_density']:.3g} fragments/bp)"
    else:
        value, note = "1×", "one gDNA level: nothing enriched above the bulk"
    return _tile("Capture", "target", value, "text", "info", note)


def _si(n: float) -> str:
    n = float(n)
    if n >= 1e6:
        return f"{n / 1e6:.1f}M"
    if n >= 1e3:
        return f"{n / 1e3:.1f}k"
    return f"{n:.0f}"


def _alignment(summary: SummaryReader) -> dict:
    a = summary["alignment_stats"]
    # Records (individual BAM alignment records).
    total = a["total_reads"]
    mapped = a["mapped_reads"]
    # Read pairs, one outcome each, summing to the total: the read fate's one denominator.
    pairs = summary["read_pairs"]
    n_pairs = pairs["total"]
    unique = summary["read_pairs.placed"]["unique"]
    multi = summary["read_pairs.placed"]["multimapping"]
    withheld = sum(summary["read_pairs.withheld"].values())
    excluded = pairs["excluded"]
    by_option = pairs["excluded_by_option"]

    fate = [
        {"label": "Placed, unique", "value": unique, "cls": _C[1]},
        {"label": "Placed, multimapping", "value": multi, "cls": _C[3]},
        {"label": "Withheld", "value": withheld, "cls": _C[6]},
        {"label": "Unmapped or filtered", "value": excluded, "cls": "fmuted"},
        {"label": "Excluded: multimapping", "value": by_option, "cls": _C[5]},
    ]
    fate = [seg for seg in fate if seg["value"] > 0]

    kpis = [
        {"l": "Read pairs", "v": n_pairs, "fmt": "count"},
        {"l": "Placed", "v": _pct(unique + multi, n_pairs), "fmt": "pct"},
        {"l": "Multimapping", "v": _pct(multi, n_pairs), "fmt": "pct"},
        {"l": "Alignment records", "v": total, "fmt": "count"},
        {"l": "Duplicates", "v": _pct(a["duplicate_reads"], total), "fmt": "pct"},
    ]
    # The table shows both units, read pairs and alignment records, each labelled.
    rows = [
        ["Read pairs", n_pairs, 1.0],
        ["  placed, unique", unique, _pct(unique, n_pairs)],
        ["  placed, multimapping", multi, _pct(multi, n_pairs)],
        ["  withheld", withheld, _pct(withheld, n_pairs)],
        ["  unmapped or filtered", excluded, _pct(excluded, n_pairs)],
        ["  excluded: multimapping", by_option, _pct(by_option, n_pairs)],
        ["Alignment records", total, 1.0],
        ["  mapped", mapped, _pct(mapped, total)],
        ["  secondary", a["secondary_reads"], _pct(a["secondary_reads"], total)],
        ["  supplementary", a["supplementary_reads"], _pct(a["supplementary_reads"], total)],
        ["  flagged proper pair", a["proper_pairs"], _pct(a["proper_pairs"], total)],
        ["  duplicates", a["duplicate_reads"], _pct(a["duplicate_reads"], total)],
        ["  QC-fail", a["qc_fail_reads"], _pct(a["qc_fail_reads"], total)],
    ]
    return {"kpis": kpis, "fate": fate, "table": rows}


def _fragments(summary: SummaryReader) -> dict:
    f = summary["fragment_stats"]
    total = f["total"]
    genic = f["genic"]
    inter = f["intergenic"]
    sp = summary["fragment_stats.splice"]
    spliced = sum(sp[k] for k in ("spliced_annotated", "spliced_implicit"))

    ctx = [
        {"label": "Genic", "value": genic, "cls": _C[1]},
        {"label": "Intergenic", "value": inter, "cls": _C[3]},
    ]
    splice_bar = [
        {"label": "Unspliced", "value": sp["unspliced"], "cls": "fmuted"},
        {"label": "Annotated", "value": sp["spliced_annotated"], "cls": _C[1]},
        {"label": "Implicit", "value": sp["spliced_implicit"], "cls": _C[2]},
    ]
    kpis = [
        {"l": "Fragments", "v": total, "fmt": "count"},
        {"l": "Genic", "v": _pct(genic, total), "fmt": "pct"},
        {"l": "Intergenic", "v": _pct(inter, total), "fmt": "pct"},
        {"l": "Spliced", "v": spliced, "fmt": "count"},
    ]
    return {"kpis": kpis, "ctx": ctx, "splice": splice_bar}


def _strand(summary: SummaryReader) -> dict:
    s = summary["strand_model"]
    d = summary["strand_model.diagnostics"]
    ci = s["ci_95"]
    eps = (ci[1] - ci[0]) / 2 if len(ci) == 2 else 0.0
    exonic = d["exonic_all_specificity"]
    return {
        "spec": s["strand_specificity"],
        "read1_sense": bool(s["read1_sense"]),
        "protocol": s["protocol"],
        "exonic_all_spec": exonic,
        "contamination_gap": d["contamination_gap"],
        "kpis": [
            {"l": "P(R1 sense)", "v": s["p_r1_sense"], "fmt": "float3"},
            {"l": "Specificity", "v": s["strand_specificity"], "fmt": "float3"},
            {"l": "95% CI", "v": f"±{eps:.4f}", "fmt": "text"},
            {"l": "Exonic (all)", "v": exonic if exonic is not None else "—", "fmt": "float3"},
        ],
    }


def _fragment_length(summary: SummaryReader) -> dict:
    """The summary table, one row per category in the order ``rigel quant`` wrote them: the category, its
    statistics, and the label it shows."""
    table = []
    for cat, d in summary["fragment_length.*"].items():
        table.append(
            [
                cat,
                d["n_observations"],
                d["mean"],
                d["std"],
                d["median"],
                d["mode"],
                round(d["overflow_fraction"] * 100, CONSTANTS.qc.overflow_percent_decimals),
                fl_label(cat),
            ]
        )
    return {"table": table}


def _quant(summary: SummaryReader) -> dict:
    q = summary["quantification"]
    annotated = q["count_annotated"]
    nascent = q["count_nascent"]
    gdna_em = q["count_gdna_em"]
    # every intergenic fragment is gDNA: unspliced, since a spliced one the scan keeps has a transcript. The pool is
    # summary.json's own, which holds the EM's share of the multimappers at their intergenic alignments beside the
    # fragments that reach no transcript, so the four pools add to its three
    inter = q["count_intergenic"]
    pools = [
        {"label": "Annotated transcripts", "value": annotated, "cls": _C[1]},
        {"label": "Nascent RNA", "value": nascent, "cls": _C[2]},
        {"label": "gDNA (EM)", "value": gdna_em, "cls": _C[8]},
        {"label": "Intergenic gDNA", "value": inter, "cls": _C[3]},
    ]
    table = [
        [
            "Annotated transcripts",
            annotated,
            q["fraction_annotated"],
            "every reference transcript, single-exon included",
        ],
        [
            "Nascent RNA",
            nascent,
            q["fraction_nascent"],
            "Rigel's synthetic spans, which model the unspliced RNA of multi-exon transcripts",
        ],
        [
            "gDNA (EM)",
            gdna_em,
            q["fraction_gdna_em"],
            "genic fragments the EM assigns to gDNA",
        ],
        [
            "Intergenic gDNA",
            inter,
            q["fraction_intergenic"],
            "no annotation overlap; gDNA by construction",
        ],
    ]
    rna_share = q["fraction_annotated"] + q["fraction_nascent"]
    kpis = [
        {"l": "Transcripts", "v": q["n_transcripts"], "fmt": "count"},
        {"l": "Genes", "v": q["n_genes"], "fmt": "count"},
        {"l": "Loci", "v": q["n_loci"], "fmt": "count"},
        {"l": "Annotated transcripts", "v": annotated, "fmt": "count"},
        {"l": "Nascent RNA", "v": nascent, "fmt": "count"},
    ]
    return {"pools": pools, "table": table, "rna_share": rna_share, "kpis": kpis}


def _genes(sub: ReportSubstrate) -> dict:
    """The gene table's rows: a gene's RNA is its annotated ``count``, since a nascent span belongs to no gene."""
    qc = CONSTANTS.qc
    gq = sub.gene_quant
    if gq is None or gq.empty:
        return {"rows": [], "total": 0, "shown": 0, "truncated": False}
    df = gq[[c for c in GENE_COLUMNS if c in gq.columns]]
    # Embed only expressed genes, capped to bound the HTML size; the rest stay in gene_quant.feather.
    if "count" in df.columns:
        df = df[df["count"] > 0]
    if "tpm" in df.columns:
        df = df.sort_values("tpm", ascending=False)
    total = len(df)
    truncated = total > qc.gene_table_rows
    df = df.head(qc.gene_table_rows)

    # Rows are built column-wise, not with iterrows.
    def col(name, cast, default):
        if name not in df.columns:
            return [default] * len(df)
        return [cast(v) for v in df[name].to_numpy()]

    # each cell as the front end shows it, and what a missing column shows: text, or a number to its decimals
    cell = {
        "tpm": (lambda v: round(float(v), qc.tpm_decimals), 0.0),
        "count": (lambda v: round(float(v), qc.count_decimals), 0.0),
        "n_transcripts": (int, 0),
    }
    cols = [col(name, *cell.get(name, (str, ""))) for name in GENE_COLUMNS]
    rows = [list(r) for r in zip(*cols)]
    return {"rows": rows, "total": total, "shown": len(rows), "truncated": truncated}


def _calibration(summary: SummaryReader, track) -> dict:
    """Two panels' worth of data: 'enrichment' (calibration's capture answer + its KPIs) and
    'density' (the genome track + per-reference table + its KPIs)."""
    cal = summary["calibration"]
    has_track = track is not None and len(track) > 0

    # Enrichment panel KPIs (capture).
    capture = _capture(summary)
    enrichment_kpis = []
    if capture is not None and capture["bulk_density"] is not None:
        enrichment_kpis = [
            {"l": "ρg bulk", "v": capture["bulk_density"], "fmt": "g4"},
            {"l": "Saturation × bulk", "v": capture["saturation"], "fmt": "g4"},
        ]
    elif capture is not None:
        enrichment_kpis = [{"l": "ρg bulk", "v": "None", "fmt": "text"}]

    # gDNA-density panel KPIs: calibration's library scalars, then the track's.
    density_kpis = [
        {"l": "Regions", "v": cal["n_regions"], "fmt": "count"},
        {"l": "ρg global", "v": cal["gdna_density_global"], "fmt": "g4"},
        {"l": "RNA sense", "v": cal["fraction_rna_sense"], "fmt": "float3"},
    ]
    if has_track:
        gf = track["fraction_gdna"].to_numpy(dtype="float64")
        cut = CONSTANTS.qc.region_gdna_fraction
        density_kpis.append({"l": "Mean gDNA fraction", "v": float(gf.mean()), "fmt": "pct"})
        density_kpis.append(
            {"l": f"Regions >{cut:.0%} gDNA", "v": int((gf > cut).sum()), "fmt": "count"}
        )

    ref_table = reference_table(track) if has_track else []
    return {
        "enrichment_kpis": enrichment_kpis,
        "density_kpis": density_kpis,
        "capture": capture,
        # an older-schema summary can hold a calibration block without a capture answer
        "calibrated": bool(summary[""]["calibration"]),
        "ref_table": ref_table,
        # the genome track's reference cap, which its caption states
        "genome_track_refs": CONSTANTS.qc.genome_track_refs,
    }


def reference_table(track) -> list:
    """Per-reference gDNA summary of every reference, by gDNA count descending (the track facets only the top few).

    Rows are ``[ref, n_regions, count_gdna, mean_density, fraction_gdna]``. ``mean_density`` is the reference's gDNA
    count per bp (Σ gDNA count / Σ region length): a per-region median is 0 on most references.
    """
    qc = CONSTANTS.qc
    length = (track["end"] - track["start"]).clip(lower=1)
    agg = (
        track.assign(_len=length)
        .groupby("ref", observed=True)
        .agg(
            n_regions=("start", "size"),
            count_gdna=("count_gdna", "sum"),
            count_rna=("count_rna", "sum"),
            len_bp=("_len", "sum"),
        )
    )
    agg["mean_density"] = (agg["count_gdna"] / agg["len_bp"]).fillna(0.0)
    agg["fraction_gdna"] = (agg["count_gdna"] / (agg["count_gdna"] + agg["count_rna"])).fillna(0.0)
    agg = agg.sort_values("count_gdna", ascending=False)
    return [
        [
            str(ref),
            int(v["n_regions"]),
            round(float(v["count_gdna"]), qc.count_decimals),
            float(v["mean_density"]),
            round(float(v["fraction_gdna"]), qc.fraction_decimals),
        ]
        for ref, v in agg.iterrows()
    ]


def _config(summary: SummaryReader) -> dict:
    cfg = summary[""]["configuration"] or {}
    groups: dict[str, dict] = {}
    # top-level scalars first
    top = {k: v for k, v in cfg.items() if not isinstance(v, dict)}
    if top:
        groups["General"] = {k: _fmt_val(v) for k, v in top.items()}
    for key, val in cfg.items():
        if isinstance(val, dict):
            groups[key.capitalize()] = {k: _fmt_val(v) for k, v in val.items()}
    return groups


def _fmt_val(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "—"
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v) if v else "—"
    return str(v)


def build_view_model(sub: ReportSubstrate) -> dict:
    """Assemble the complete JSON-serializable view model for the report."""
    s = SummaryReader(sub.summary)
    top = s[""]
    return {
        "meta": {
            "sample": sub.sample_name,
            "index": s["input"]["index_dir"],
            "created": top["timestamp"],
            "version": top["rigel_version"],
            "schema_version": top["schema_version"],
            "warnings": sub.warnings,
        },
        "verdicts": _verdicts(s),
        "alignment": _alignment(s),
        "fragments": _fragments(s),
        "strand": _strand(s),
        "fl": _fragment_length(s),
        "quant": _quant(s),
        "calibration": _calibration(s, sub.calibration_track),
        "genes": _genes(sub),
        "config": _config(s),
    }
