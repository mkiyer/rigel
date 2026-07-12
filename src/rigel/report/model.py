"""Transform a :class:`ReportSubstrate` into the render-ready view model.

The view model is a plain, JSON-serializable dict consumed by the report's
front-end. Keeping this transform in Python (not JavaScript) means the HTML
layer only lays out already-computed values — no business logic in the browser.

All formatting (thousands separators, percentages) is deferred to the front-end;
this module emits raw numbers so the same values drive tooltips, tables, and
chart data without rounding drift.
"""

from __future__ import annotations

from .substrate import ReportSubstrate

# Categorical colour slots (indices into the validated data-viz palette). The
# front-end resolves these to hex per theme; here we only assign roles.
_C = {i: f"f{i}" for i in range(1, 9)}


def _pct(numer: float, denom: float) -> float:
    return (numer / denom) if denom else 0.0


def _verdicts(summary: dict, capture: dict | None = None) -> list[dict]:
    """Headline QC tiles: mapping, strandedness, gDNA, usable fragments, + capture
    enrichment (descriptive, mass-weighted) when a gDNA track is present."""
    out: list[dict] = []
    aln = summary.get("alignment_stats", {})
    frag = summary.get("fragment_stats", {})
    strand = summary.get("strand_model", {})
    quant = summary.get("quantification", {})

    # Mapping rate
    total = aln.get("total_reads", 0)
    mapped = aln.get("mapped_reads", 0)
    mrate = _pct(mapped, total)
    out.append(
        {
            "k": "Mapping rate",
            "icon": "check",
            "v": f"{mrate * 100:.1f}%",
            "s": "good" if mrate >= 0.75 else "warning" if mrate >= 0.5 else "critical",
            "n": f"{mapped:,} of {total:,} reads mapped",
        }
    )

    # Strandedness
    spec = strand.get("strand_specificity", 0.5)
    protocol = strand.get("protocol", "?")
    if spec >= 0.75:
        s, label = "good", "Reverse (R1−)" if not strand.get("read1_sense") else "Forward (R1+)"
    elif spec >= 0.6:
        s, label = "warning", "Weakly stranded"
    else:
        s, label = "warning", "Unstranded"
    out.append(
        {
            "k": "Strandedness",
            "icon": "arrow",
            "v": label,
            "s": s,
            "n": f"{spec * 100:.1f}% strand-specific · {protocol}",
        }
    )

    # gDNA contamination
    gfrac = quant.get("gdna_fraction", 0.0)
    out.append(
        {
            "k": "gDNA contamination",
            "icon": "warn",
            "v": f"{gfrac * 100:.1f}%",
            "s": "good" if gfrac < 0.1 else "warning" if gfrac < 0.3 else "critical",
            "n": "EM-corrected; incl. intergenic",
        }
    )

    # Usable (genic) fragments
    ftot = frag.get("total", 0)
    genic = frag.get("genic", 0)
    grate = _pct(genic, ftot)
    out.append(
        {
            "k": "Usable fragments",
            "icon": "check",
            "v": _si(genic),
            "s": "good" if grate >= 0.6 else "warning",
            "n": f"{grate * 100:.1f}% genic of {_si(ftot)}",
        }
    )

    # Capture enrichment — descriptive only. Variable capture performance is
    # expected and we do not interpret an enriched mode as good/bad; the tile
    # leads with the robust enriched-vs-median fold, the Calibration panel has the
    # full metric set. Styled neutral ("info").
    if capture:
        if capture.get("enriched"):
            out.append(
                {
                    "k": "Capture enrichment",
                    "icon": "target",
                    "v": _fold(capture.get("fold_vs_median", 1.0)),
                    "s": "info",
                    "n": f"on-target mode vs median; {capture.get('mass_frac_ontarget', 0) * 100:.0f}% of gDNA mass on-target",
                }
            )
        else:
            out.append(
                {
                    "k": "Capture enrichment",
                    "icon": "target",
                    "v": "None",
                    "s": "info",
                    "n": "no on-target gDNA density mode above the median",
                }
            )
    return out


def _fold(x: float) -> str:
    """Compact fold-change label (e.g. 45× / 4.4e4×)."""
    x = float(x)
    if x >= 1000:
        return f"{x:.0e}×".replace("e+0", "e").replace("e+", "e")
    return f"{x:.0f}×"


def _si(n: float) -> str:
    n = float(n)
    if n >= 1e6:
        return f"{n / 1e6:.1f}M"
    if n >= 1e3:
        return f"{n / 1e3:.1f}k"
    return f"{n:.0f}"


def _alignment(summary: dict) -> dict:
    a = summary.get("alignment_stats", {})
    # Record-level (individual BAM alignment records).
    total = a.get("total_reads", 0)
    mapped = a.get("mapped_reads", 0)
    # Read-NAME-group level (one per read/pair) — the "read fate" denominator.
    unique = a.get("unique_reads", 0)
    multi = a.get("multimapping_reads", 0)
    # read_groups may be absent on pre-fix summaries; fall back to unique+multi
    # (mapped groups only) rather than mixing in record-level unmapped counts.
    groups = a.get("read_groups", unique + multi)
    unmapped_grp = max(0, groups - unique - multi)

    # Fate over read groups — every segment shares the read-group denominator,
    # so it composes correctly (the earlier bug mixed records with groups).
    fate = [
        {"label": "Uniquely mapped", "value": unique, "cls": _C[1]},
        {"label": "Multi-mapping", "value": multi, "cls": _C[3]},
    ]
    if unmapped_grp > 0:
        fate.append({"label": "Unmapped", "value": unmapped_grp, "cls": "fmuted"})

    kpis = [
        {"l": "Read groups", "v": _si(groups)},
        {"l": "Uniquely mapped", "v": f"{_pct(unique, groups) * 100:.1f}", "u": "%"},
        {"l": "Multi-mapping", "v": f"{_pct(multi, groups) * 100:.1f}", "u": "%"},
        {"l": "Alignment records", "v": _si(total)},
        {"l": "Duplicates", "v": f"{_pct(a.get('duplicate_reads', 0), total) * 100:.1f}", "u": "%"},
    ]
    # Table keeps both unit families, each labeled so they are not conflated.
    rows = [
        ["Read groups (reads/pairs)", groups, 1.0],
        ["  uniquely mapped", unique, _pct(unique, groups)],
        ["  multi-mapping", multi, _pct(multi, groups)],
        ["  unmapped", unmapped_grp, _pct(unmapped_grp, groups)],
        ["Alignment records", total, 1.0],
        ["  mapped", mapped, _pct(mapped, total)],
        ["  secondary", a.get("secondary_reads", 0), _pct(a.get("secondary_reads", 0), total)],
        [
            "  supplementary",
            a.get("supplementary_reads", 0),
            _pct(a.get("supplementary_reads", 0), total),
        ],
        ["  proper pairs", a.get("proper_pairs", 0), _pct(a.get("proper_pairs", 0), total)],
        ["  duplicates", a.get("duplicate_reads", 0), _pct(a.get("duplicate_reads", 0), total)],
        ["  QC-fail", a.get("qc_fail_reads", 0), _pct(a.get("qc_fail_reads", 0), total)],
    ]
    return {"kpis": kpis, "fate": fate, "table": rows}


def _fragments(summary: dict) -> dict:
    f = summary.get("fragment_stats", {})
    total = f.get("total", 0)
    genic = f.get("genic", 0)
    inter = f.get("intergenic", 0)
    chim = f.get("chimeric", 0)
    sp = f.get("splice", {})
    spliced = sum(
        sp.get(k, 0)
        for k in ("spliced_annotated", "spliced_unannotated", "spliced_implicit", "splice_artifact")
    )

    ctx = [
        {"label": "Genic", "value": genic, "cls": _C[1]},
        {"label": "Intergenic", "value": inter, "cls": _C[3]},
        {"label": "Chimeric", "value": chim, "cls": _C[6]},
    ]
    chim_bar = [
        {"label": "Cis, same strand", "value": f.get("chimeric_cis_same", 0), "cls": _C[2]},
        {"label": "Trans", "value": f.get("chimeric_trans", 0), "cls": _C[5]},
        {"label": "Cis, opposite strand", "value": f.get("chimeric_cis_diff", 0), "cls": _C[6]},
    ]
    splice_bar = [
        {"label": "Unspliced", "value": sp.get("unspliced", 0), "cls": "fmuted"},
        {"label": "Annotated", "value": sp.get("spliced_annotated", 0), "cls": _C[1]},
        {"label": "Implicit", "value": sp.get("spliced_implicit", 0), "cls": _C[2]},
        {"label": "Unannotated", "value": sp.get("spliced_unannotated", 0), "cls": _C[3]},
        {"label": "Artifact", "value": sp.get("splice_artifact", 0), "cls": _C[6]},
    ]
    kpis = [
        {"l": "Fragments", "v": _si(total)},
        {"l": "Genic", "v": f"{_pct(genic, total) * 100:.1f}", "u": "%"},
        {"l": "Intergenic", "v": f"{_pct(inter, total) * 100:.1f}", "u": "%"},
        {"l": "Chimeric", "v": f"{_pct(chim, total) * 100:.1f}", "u": "%"},
        {"l": "Spliced", "v": _si(spliced)},
    ]
    return {
        "kpis": kpis,
        "ctx": ctx,
        "chim": chim_bar,
        "splice": splice_bar,
        "blacklisted": sp.get("sj_blacklisted", 0),
    }


def _strand(summary: dict) -> dict:
    s = summary.get("strand_model", {})
    d = s.get("diagnostics", {})
    ci = s.get("ci_95", [0.0, 1.0])
    eps = (ci[1] - ci[0]) / 2 if len(ci) == 2 else 0.0
    return {
        "spec": s.get("strand_specificity", 0.5),
        "p_r1_sense": s.get("p_r1_sense", 0.5),
        "read1_sense": bool(s.get("read1_sense", True)),
        "protocol": s.get("protocol", "?"),
        "n": s.get("n_training_fragments", 0),
        "ci": ci,
        "exonic_all_spec": d.get("exonic_all_specificity"),
        "contamination_gap": d.get("contamination_gap"),
        "kpis": [
            {"l": "P(R1 sense)", "v": f"{s.get('p_r1_sense', 0.5):.3f}"},
            {"l": "Specificity", "v": f"{s.get('strand_specificity', 0.5):.3f}"},
            {"l": "95% CI", "v": f"±{eps:.4f}"},
            {
                "l": "Exonic (all)",
                "v": f"{d.get('exonic_all_specificity', float('nan')):.3f}"
                if d.get("exonic_all_specificity") is not None
                else "—",
            },
        ],
    }


def _fragment_length(sub: ReportSubstrate) -> dict:
    fl = sub.summary.get("fragment_length", {})
    # summary table (one row per category)
    order = [
        "global",
        "gdna",
        "rna",
        "unspliced",
        "spliced_annot",
        "spliced_unannot",
        "spliced_implicit",
        "splice_artifact",
    ]
    present = [c for c in order if c in fl] + [c for c in fl if c not in order]
    table = []
    for cat in present:
        d = fl[cat]
        table.append(
            [
                cat,
                d.get("n_observations", 0),
                d.get("mean", 0),
                d.get("std", 0),
                d.get("median", 0),
                d.get("mode", 0),
                round(d.get("overflow_fraction", 0.0) * 100, 2),
            ]
        )
    return {"table": table, "categories": present, "has_hist": sub.fragment_lengths is not None}


def _quant(summary: dict) -> dict:
    q = summary.get("quantification", {})
    mrna = q.get("mrna_total", 0.0)
    nrna = q.get("nrna_total", 0.0)
    gdna_em = q.get("gdna_total", 0.0)
    inter = q.get("intergenic_total", 0.0)
    pools = [
        {"label": "Mature mRNA", "value": mrna, "cls": _C[1]},
        {"label": "Nascent RNA", "value": nrna, "cls": _C[2]},
        {"label": "gDNA (EM)", "value": gdna_em, "cls": _C[8]},
        {"label": "Intergenic gDNA", "value": inter, "cls": _C[3]},
    ]
    table = [
        ["Mature mRNA", mrna, q.get("mrna_fraction", 0.0), "spliced + EM-assigned transcript mass"],
        [
            "Nascent RNA",
            nrna,
            q.get("nrna_fraction", 0.0),
            "synthetic + annotated single-exon spans",
        ],
        [
            "gDNA — EM genic",
            gdna_em,
            q.get("gdna_em_fraction", 0.0),
            "re-attributed from ambiguous transcript overlap",
        ],
        [
            "gDNA — intergenic",
            inter,
            q.get("intergenic_fraction", 0.0),
            "no annotation overlap; gDNA by construction",
        ],
    ]
    rna_share = q.get("mrna_fraction", 0.0) + q.get("nrna_fraction", 0.0)
    kpis = [
        {"l": "Transcripts", "v": _si(q.get("n_transcripts", 0))},
        {"l": "Genes", "v": _si(q.get("n_genes", 0))},
        {"l": "Loci", "v": _si(q.get("n_loci", 0))},
        {"l": "mRNA", "v": _si(mrna)},
        {"l": "nRNA", "v": _si(nrna)},
    ]
    return {"pools": pools, "table": table, "rna_share": rna_share, "kpis": kpis}


def _genes(sub: ReportSubstrate, max_rows: int = 20000) -> dict:
    gq = sub.gene_quant
    if gq is None or gq.empty:
        return {"rows": [], "total": 0, "shown": 0, "truncated": False}
    want = [
        "gene_name",
        "gene_id",
        "gene_type",
        "ref",
        "strand",
        "tpm",
        "count",
        "mature_count",
        "nascent_count",
        "n_transcripts",
    ]
    df = gq[[c for c in want if c in gq.columns]].copy()
    # Embed only expressed genes for lookup; the long tail of zero-count genes
    # stays in gene_quant.feather. Cap the embed to keep the HTML a sane size.
    if "count" in df.columns:
        df = df[df["count"] > 0]
    if "tpm" in df.columns:
        df = df.sort_values("tpm", ascending=False)
    total = len(df)
    truncated = total > max_rows
    df = df.head(max_rows)
    rows = []
    for _, r in df.iterrows():
        rows.append(
            [
                str(r.get("gene_name", "")),
                str(r.get("gene_id", "")),
                str(r.get("gene_type", "")),
                str(r.get("ref", "")),
                str(r.get("strand", "")),
                round(float(r.get("tpm", 0)), 2),
                round(float(r.get("mature_count", 0)), 1),
                round(float(r.get("nascent_count", 0)), 1),
                int(r.get("n_transcripts", 0)),
            ]
        )
    return {"rows": rows, "total": total, "shown": len(rows), "truncated": truncated}


def _calibration(sub: ReportSubstrate, capture: dict | None = None) -> dict:
    cal = sub.summary.get("calibration") or {}
    track = sub.calibration_track
    kpis = [
        {"l": "Regions", "v": _si(cal.get("n_regions", 0))},
        {"l": "ρg global", "v": f"{cal.get('gdna_density_global', 0):.4g}"},
        {"l": "RNA sense", "v": f"{cal.get('rna_sense_frac', 0):.3f}"},
    ]
    has_track = track is not None and len(track) > 0
    if has_track:
        gf = track["gdna_frac"].to_numpy(dtype="float64")
        kpis.append({"l": "Mean gDNA frac", "v": f"{float(gf.mean()) * 100:.1f}", "u": "%"})
        kpis.append({"l": "Regions >50% gDNA", "v": _si(int((gf > 0.5).sum()))})
    if capture and capture.get("enriched"):
        kpis.append({"l": "Enr. vs median", "v": _fold(capture.get("fold_vs_median", 1.0))})
        kpis.append({"l": "Peak-to-peak", "v": _fold(capture.get("fold_peak_to_peak", 1.0))})
        kpis.append(
            {
                "l": "On-target mass",
                "v": f"{capture.get('mass_frac_ontarget', 0) * 100:.1f}",
                "u": "%",
            }
        )
    return {"kpis": kpis, "has_track": has_track, "capture": capture}


def _config(summary: dict) -> dict:
    cfg = summary.get("configuration", {})
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


def build_view_model(sub: ReportSubstrate, capture: dict | None = None) -> dict:
    """Assemble the complete JSON-serializable view model for the report."""
    s = sub.summary
    return {
        "meta": {
            "sample": sub.sample_name,
            "bam": s.get("input", {}).get("bam_file"),
            "index": s.get("input", {}).get("index_dir"),
            "created": s.get("timestamp"),
            "version": s.get("rigel_version"),
            "schema_version": s.get("schema_version"),
            "warnings": sub.warnings,
        },
        "verdicts": _verdicts(s, capture),
        "alignment": _alignment(s),
        "fragments": _fragments(s),
        "strand": _strand(s),
        "fl": _fragment_length(sub),
        "quant": _quant(s),
        "calibration": _calibration(sub, capture),
        "genes": _genes(sub),
        "config": _config(s),
    }
