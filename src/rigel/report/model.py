"""Transform a :class:`ReportSubstrate` into the render-ready view model.

The view model is a plain, JSON-serializable dict consumed by the report's
front-end. Keeping this transform in Python (not JavaScript) means the HTML
layer only lays out already-computed values — no business logic in the browser.

Most values are emitted as raw numbers and formatted by the front-end
(thousands separators, percentages), so the same values drive tooltips, tables,
and chart data. The exceptions are finished here: the headline tiles' captions
and a few text KPIs are formatted strings, and some table cells are rounded (the
fragment-length overflow, as a percentage; the gene table's TPM and counts; the
per-reference gDNA count and fraction). Every cut a tile draws, the gene table's
cap, the genome track's reference cap its caption states and every rounding are
``CONSTANTS.qc``'s. The summary is read only through
:class:`~rigel.report.substrate.SummaryReader`, so only the keys it lists.
"""

from __future__ import annotations

from ..config import CONSTANTS
from .specs import fl_label
from .substrate import ReportSubstrate, SummaryReader

# Categorical colour slots (indices into the validated data-viz palette). The
# front-end resolves these to hex per theme; here we only assign roles.
_C = {i: f"f{i}" for i in range(1, 9)}

#: The ``gene_quant`` columns the gene table shows, in its row order. A column the table lacks reads as empty
#: cells, so ``tests/test_report.py`` holds each one to what the writer writes.
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
    """Calibration's capture answer: the bulk density of its gDNA landscape — the level, in fragments per
    bp, every capture efficiency is read against — and the landscape's saturated level as a multiple of
    that bulk (1: no level above it). ``bulk_density`` is ``None`` when no landscape was fitted; the
    answer is ``None`` when the summary carries none."""
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
    out: list[dict] = []
    qc = CONSTANTS.qc
    aln = summary["alignment_stats"]
    frag = summary["fragment_stats"]
    strand = summary["strand_model"]
    quant = summary["quantification"]

    # Mapping rate
    total = aln["total_reads"]
    mapped = aln["mapped_reads"]
    mrate = _pct(mapped, total)
    out.append(
        {
            "k": "Mapping rate",
            "icon": "check",
            "v": mrate,
            "fmt": "pct",
            "s": "good"
            if mrate >= qc.mapping_rate_good
            else "warning"
            if mrate >= qc.mapping_rate_warning
            else "critical",
            "n": f"{mapped:,} of {total:,} reads mapped",
        }
    )

    # Strandedness: the strand model's own verdict, the one calibration's strand channel reads
    spec = strand["strand_specificity"]
    protocol = strand["protocol"]
    if strand["strand_live"]:
        s, label = "good", "Reverse (R1−)" if not strand["read1_sense"] else "Forward (R1+)"
    else:
        s, label = "warning", "Unstranded"
    out.append(
        {
            "k": "Strandedness",
            "icon": "arrow",
            "v": label,
            "fmt": "text",
            "s": s,
            "n": f"{spec * 100:.1f}% strand-specific · {protocol}",
        }
    )

    # gDNA contamination
    gfrac = quant["fraction_gdna"]
    out.append(
        {
            "k": "gDNA contamination",
            "icon": "warn",
            "v": gfrac,
            "fmt": "pct",
            "s": "good"
            if gfrac < qc.gdna_fraction_good
            else "warning"
            if gfrac < qc.gdna_fraction_warning
            else "critical",
            "n": "EM-corrected; incl. unspliced intergenic",
        }
    )

    # Usable (genic) fragments
    ftot = frag["total"]
    genic = frag["genic"]
    grate = _pct(genic, ftot)
    out.append(
        {
            "k": "Usable fragments",
            "icon": "check",
            "v": genic,
            "fmt": "count",
            "s": "good" if grate >= qc.genic_fraction_good else "warning",
            "n": f"{grate * 100:.1f}% genic of {_si(ftot)}",
        }
    )

    # Capture — descriptive only: variable capture performance is expected, so an enriched level is
    # neither good nor bad. Styled neutral ("info").
    capture = _capture(summary)
    if capture is not None:
        if capture["bulk_density"] is None:
            value, note = "None", "no gDNA landscape: lengths not corrected for capture"
        elif capture["saturation"] > 1.0:
            value = f"{capture['saturation']:,.0f}×"
            note = (
                f"saturated gDNA level over the bulk ({capture['bulk_density']:.3g} fragments/bp)"
            )
        else:
            value, note = "1×", "one gDNA level: nothing enriched above the bulk"
        out.append(
            {"k": "Capture", "icon": "target", "v": value, "fmt": "text", "s": "info", "n": note}
        )
    return out


def _si(n: float) -> str:
    n = float(n)
    if n >= 1e6:
        return f"{n / 1e6:.1f}M"
    if n >= 1e3:
        return f"{n / 1e3:.1f}k"
    return f"{n:.0f}"


def _alignment(summary: SummaryReader) -> dict:
    a = summary["alignment_stats"]
    # Record-level (individual BAM alignment records).
    total = a["total_reads"]
    mapped = a["mapped_reads"]
    # Read-NAME-group level (one per read/pair) — the "read fate" denominator.
    unique = a["unique_reads"]
    multi = a["multimapping_reads"]
    groups = a["read_groups"]
    unmapped_grp = max(0, groups - unique - multi)

    # Fate over read groups: every segment shares the read-group denominator.
    fate = [
        {"label": "Uniquely mapped", "value": unique, "cls": _C[1]},
        {"label": "Multi-mapping", "value": multi, "cls": _C[3]},
    ]
    if unmapped_grp > 0:
        fate.append({"label": "Unmapped", "value": unmapped_grp, "cls": "fmuted"})

    kpis = [
        {"l": "Read groups", "v": groups, "fmt": "count"},
        {"l": "Uniquely mapped", "v": _pct(unique, groups), "fmt": "pct"},
        {"l": "Multi-mapping", "v": _pct(multi, groups), "fmt": "pct"},
        {"l": "Alignment records", "v": total, "fmt": "count"},
        {"l": "Duplicates", "v": _pct(a["duplicate_reads"], total), "fmt": "pct"},
    ]
    # Table keeps both unit families, each labeled so they are not conflated.
    rows = [
        ["Read groups (reads/pairs)", groups, 1.0],
        ["  uniquely mapped", unique, _pct(unique, groups)],
        ["  multi-mapping", multi, _pct(multi, groups)],
        ["  unmapped", unmapped_grp, _pct(unmapped_grp, groups)],
        ["Alignment records", total, 1.0],
        ["  mapped", mapped, _pct(mapped, total)],
        ["  secondary", a["secondary_reads"], _pct(a["secondary_reads"], total)],
        ["  supplementary", a["supplementary_reads"], _pct(a["supplementary_reads"], total)],
        ["  proper pairs", a["proper_pairs"], _pct(a["proper_pairs"], total)],
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
    spliced = sum(sp[k] for k in ("spliced_annotated", "spliced_implicit", "splice_artifact"))

    ctx = [
        {"label": "Genic", "value": genic, "cls": _C[1]},
        {"label": "Intergenic", "value": inter, "cls": _C[3]},
    ]
    splice_bar = [
        {"label": "Unspliced", "value": sp["unspliced"], "cls": "fmuted"},
        {"label": "Annotated", "value": sp["spliced_annotated"], "cls": _C[1]},
        {"label": "Implicit", "value": sp["spliced_implicit"], "cls": _C[2]},
        {"label": "Artifact", "value": sp["splice_artifact"], "cls": _C[6]},
    ]
    kpis = [
        {"l": "Fragments", "v": total, "fmt": "count"},
        {"l": "Genic", "v": _pct(genic, total), "fmt": "pct"},
        {"l": "Intergenic", "v": _pct(inter, total), "fmt": "pct"},
        {"l": "Spliced", "v": spliced, "fmt": "count"},
    ]
    # Blacklist provenance so a 0-artifact bar can be read correctly:
    #   loaded is False → no blacklist in the index; detection is OFF
    #   loaded is True  → detection ON; `size` junctions active
    return {
        "kpis": kpis,
        "ctx": ctx,
        "splice": splice_bar,
        "blacklisted": sp["sj_blacklisted"],
        "blacklist": {"size": sp["sj_blacklist_size"], "loaded": sp["sj_blacklist_loaded"]},
    }


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
    # every intergenic fragment is gDNA: unspliced, since a spliced one the scan keeps has a transcript
    inter = summary["fragment_stats"]["intergenic"]
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
    # Embed only expressed genes for lookup; the long tail of zero-count genes
    # stays in gene_quant.feather. Cap the embed to keep the HTML a sane size.
    if "count" in df.columns:
        df = df[df["count"] > 0]
    if "tpm" in df.columns:
        df = df.sort_values("tpm", ascending=False)
    total = len(df)
    truncated = total > qc.gene_table_rows
    df = df.head(qc.gene_table_rows)

    # Build the row arrays column-wise (vectorized) — far faster than iterrows.
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

    ref_table = _reference_table(track) if has_track else []
    return {
        "enrichment_kpis": enrichment_kpis,
        "density_kpis": density_kpis,
        "capture": capture,
        # a summary from before calibration's capture answer (0.7.1, schema 2) has a block but no answer
        "calibrated": bool(summary[""]["calibration"]),
        "ref_table": ref_table,
        # the genome track's reference cap, which its caption states
        "genome_track_refs": CONSTANTS.qc.genome_track_refs,
    }


def _reference_table(track) -> list:
    """Per-reference gDNA summary (all references), sorted by gDNA count desc.

    Rows: ``[ref, n_regions, count_gdna, mean_density, fraction_gdna]``. This is the
    complete, compact view of every reference — the track chart only facets the
    top few by gDNA count, so users never sift through hundreds of panels.

    ``mean_density`` is the reference-wide gDNA count per bp (Σ gDNA count / Σ region
    length) — NOT a median. On most references the median per-region density is 0
    (most regions carry no gDNA), which is a useless column; the count-per-bp
    average discriminates references and is what "density" means at this scale.
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
