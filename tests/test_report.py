"""The `rigel report` HTML builder, from a synthesized substrate rather than a pipeline run.

A minimal but realistic substrate — a v3 ``summary.json`` and its companion feathers — is written
to a temp directory, and the loader, the view model, the chart specifications and the full HTML
build run against it. The report must be self-contained, inlining its runtime, and must honour a
custom output path. The capture panel shows calibration's own answer, and the front end reads only the
keys and format tags the view model writes, and escapes every string a user or an annotation chose before it
becomes markup (read off ``report.js``'s source: the suite runs no JavaScript). Every chart the page embeds
must compile and read every mark property it sets. Vega-specific assertions are conditional on
``vl-convert-python``.
"""

import importlib.util
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import numpy as np

from rigel.report.build import build_report
from rigel.report.html import _asset
from rigel.report.model import _reference_table, build_view_model
from rigel.report.specs import build_charts, build_fl_specs, genome_track_spec
from rigel.report.substrate import SubstrateError, load_substrate

_HAS_VEGA = importlib.util.find_spec("vl_convert") is not None


def _write_substrate(d: Path) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    summary = {
        "schema_version": 3,
        "rigel_version": "0.7.0",
        "timestamp": "2026-07-11 09:42",
        "input": {"bam_file": "/data/SampleX.bam", "index_dir": "/refs/gencode"},
        "configuration": {"em": {"mode": "vbem", "n_threads": 8}, "seed": 42},
        "alignment_stats": {
            "total_reads": 1000,
            "mapped_reads": 950,
            "unique_reads": 800,
            "multimapping_reads": 150,
            "proper_pairs": 900,
            "duplicate_reads": 60,
            "qc_fail_reads": 5,
            "read_groups": 980,
        },
        "fragment_stats": {
            "total": 700,
            "genic": 600,
            "intergenic": 70,
            "chimeric": 30,
            "chimeric_trans": 10,
            "chimeric_cis_same": 15,
            "chimeric_cis_diff": 5,
            "with_annotated_sj": 400,
            "with_unannotated_sj": 20,
            "splice": {
                "unspliced": 250,
                "spliced_annotated": 400,
                "spliced_unannotated": 20,
                "spliced_implicit": 25,
                "splice_artifact": 5,
                "sj_blacklisted": 3,
            },
        },
        "strand_model": {
            "protocol": "R1-antisense",
            "strand_specificity": 0.98,
            "p_r1_sense": 0.02,
            "read1_sense": False,
            "n_training_fragments": 400,
            "posterior_variance": 0.0001,
            "ci_95": [0.975, 0.985],
            "diagnostics": {
                "exonic_all_specificity": 0.93,
                "exonic_all_p_r1_sense": 0.07,
                "exonic_all_n_fragments": 600,
                "contamination_gap": 0.05,
            },
        },
        "quantification": {
            "n_transcripts": 120,
            "n_genes": 40,
            "n_loci": 35,
            "mrna_total": 500.0,
            "nrna_total": 60.0,
            "gdna_total": 20.0,
            "intergenic_total": 70,
            "mrna_fraction": 0.77,
            "nrna_fraction": 0.09,
            "gdna_fraction": 0.14,
            "gdna_em_fraction": 0.03,
            "intergenic_fraction": 0.11,
        },
        "fragment_length": {
            "global": {
                "n_observations": 700,
                "mean": 200.0,
                "std": 40.0,
                "median": 195.0,
                "mode": 190,
                "max_size": 1000,
                "overflow_count": 0.0,
                "overflow_fraction": 0.0,
            },
            "rna": {
                "n_observations": 400,
                "mean": 205.0,
                "std": 35.0,
                "median": 200.0,
                "mode": 195,
                "max_size": 1000,
                "overflow_count": 0.0,
                "overflow_fraction": 0.0,
            },
            "gdna": {
                "n_observations": 120,
                "mean": 175.0,
                "std": 60.0,
                "median": 160.0,
                "mode": 150,
                "max_size": 1000,
                "overflow_count": 1.0,
                "overflow_fraction": 0.008,
            },
        },
    }
    (d / "summary.json").write_text(json.dumps(summary))

    # fragment_lengths.feather — a couple of gaussian-ish bumps
    rows = []
    for cat, mu in (("global", 200), ("rna", 205), ("gdna", 175)):
        for length in range(mu - 40, mu + 41, 5):
            rows.append((cat, length, max(0.0, 50.0 - abs(length - mu))))
    pd.DataFrame(rows, columns=["category", "length", "count"]).astype(
        {"length": "int32", "count": "float64"}
    ).to_feather(d / "fragment_lengths.feather")

    pd.DataFrame(
        {
            "gene_id": [f"ENSG{i:05d}" for i in range(5)],
            "gene_name": ["GAPDH", "ACTB", "TP53", "MYC", "MALAT1"],
            "count": [4200.0, 3800.0, 120.0, 300.0, 900.0],
            "count_spliced": [4100.0, 3700.0, 110.0, 280.0, 40.0],
            "tpm": [42000.0, 38000.0, 1200.0, 2900.0, 9000.0],
            "n_transcripts": [3, 4, 6, 2, 1],
        }
    ).to_feather(d / "gene_quant.feather")
    return d


def _enriched_track() -> pd.DataFrame:
    """A capture-like gDNA track: many low-density off-target regions carrying little gDNA mass, and a
    few high-density on-target regions carrying most of it."""
    rng = np.random.default_rng(0)
    low_d = np.exp(rng.normal(-9.0, 0.4, 4000))
    high_d = np.exp(rng.normal(0.0, 0.4, 120))
    dens = np.concatenate([low_d, high_d])
    gmass = np.concatenate([np.full(4000, 0.01), np.full(120, 50.0)])
    return pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1"] * len(dens)),
            "start": np.arange(len(dens)) * 100,
            "end": np.arange(len(dens)) * 100 + 50,
            "gdna_mass": gmass,
            "rna_mass": np.ones(len(dens)),
            "gdna_density": dens,
            "gdna_frac": np.clip(dens, 0, 1),
        }
    )


def _payload(html: str) -> dict:
    """The view model and chart specs a built report embeds."""
    m = re.search(r'<script id="rigel-data" type="application/json">(.*?)</script>', html, re.S)
    return json.loads(m.group(1).replace("<\\/", "</"))


def _mark_paths(node, path=()):
    """The path of every mark definition given as an object in a Vega-Lite spec."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "mark" and isinstance(value, dict):
                yield (*path, key)
            elif key != "data":
                yield from _mark_paths(value, (*path, key))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _mark_paths(value, (*path, i))


def _at(node, path):
    for key in path:
        node = node[key]
    return node


def test_load_substrate_missing(tmp_path):
    with pytest.raises(SubstrateError):
        load_substrate(tmp_path / "does_not_exist")
    (tmp_path / "empty").mkdir()
    with pytest.raises(SubstrateError):
        load_substrate(tmp_path / "empty")  # no summary.json


def test_view_model_shape(tmp_path):
    d = _write_substrate(tmp_path / "run")
    sub = load_substrate(d)
    assert sub.schema_version == 3
    assert sub.sample_name == "SampleX"
    vm = build_view_model(sub)

    assert len(vm["verdicts"]) == 4
    # formatting boundary: KPIs / verdicts carry raw values + a fmt tag (JS formats)
    for kpi in vm["alignment"]["kpis"]:
        assert "fmt" in kpi and "u" not in kpi
    numeric = next(k for k in vm["alignment"]["kpis"] if k["fmt"] in ("pct", "count"))
    assert isinstance(numeric["v"], (int, float))
    assert all("fmt" in v for v in vm["verdicts"])
    # alignment fate never double-counts past the total
    assert sum(seg["value"] for seg in vm["alignment"]["fate"]) <= 1000
    # splice bar carries implicit + artifact
    labels = {s["label"] for s in vm["fragments"]["splice"]}
    assert {"Implicit", "Artifact"} <= labels
    # strand contamination gap surfaced
    assert vm["strand"]["contamination_gap"] == 0.05
    # the fragment-length table keeps the categories in the order summary.json lists them
    assert [row[0] for row in vm["fl"]["table"]] == list(sub.summary["fragment_length"])
    # gene table populated and sorted by tpm desc
    assert vm["genes"]["rows"][0][0] == "GAPDH"
    assert vm["genes"]["total"] == 5


def test_fl_specs_present(tmp_path):
    d = _write_substrate(tmp_path / "run")
    sub = load_substrate(d)
    specs = build_fl_specs(sub.fragment_lengths)
    assert "overlay" in specs and "small_multiples" in specs


def test_genome_track_spec_bins_per_ref():
    track = pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1"] * 4 + ["chr2"] * 2),
            "start": [0, 1000, 2000, 3000, 0, 5000],
            "end": [1000, 2000, 3000, 4000, 5000, 10000],
            "gdna_mass": [1.0, 8.0, 2.0, 1.0, 30.0, 90.0],  # chr1 sum 12, chr2 sum 120
            "rna_mass": [9.0, 2.0, 8.0, 9.0, 7.0, 1.0],
            "gdna_density": [0.01, 0.09, 0.02, 0.01, 0.03, 0.10],
            "gdna_frac": [0.1, 0.8, 0.2, 0.1, 0.3, 0.9],
        }
    )
    spec = genome_track_spec(track)
    assert spec is not None
    refs = {row["ref"] for row in spec["data"]["values"]}
    assert refs == {"chr1", "chr2"}
    # log + independent y per facet (so a high-density ref can't flatten others)
    assert spec["spec"]["encoding"]["y"]["scale"]["type"] == "log"
    assert spec["resolve"]["scale"]["y"] == "independent"
    # top-N ranks by gDNA mass; chr2 (120) outranks chr1 (12)
    top1 = genome_track_spec(track, top_n=1)
    assert {row["ref"] for row in top1["data"]["values"]} == {"chr2"}
    # reference table: every reference, sorted by gDNA mass desc
    rt = _reference_table(track)
    assert [r[0] for r in rt] == ["chr2", "chr1"]
    assert rt[0][1] == 2  # chr2 n_regions
    # build_charts merges genome in when a track is present
    stub = SimpleNamespace(
        fragment_lengths=None,
        calibration_track=track,
        summary={},
    )
    assert set(build_charts(stub)) == {"genome"}
    empty = SimpleNamespace(
        fragment_lengths=None,
        calibration_track=None,
        summary={},
    )
    assert build_charts(empty) == {}


def _with_calibration(d: Path, **keys) -> Path:
    """Give the substrate's ``summary.json`` a calibration block carrying ``keys``."""
    path = d / "summary.json"
    summary = json.loads(path.read_text())
    summary["calibration"] = {"gdna_density_global": 0.03, "n_regions": 1760, **keys}
    path.write_text(json.dumps(summary))
    return d


def _js_function(name: str) -> str:
    """The body of one top-level function of the report's front end."""
    js = _asset("report.js")
    m = re.search(rf"\n  function {name}\(.*?\) \{{\n(.*?)\n  \}}\n", js, re.S)
    assert m, f"report.js has no function {name}"
    return m.group(1)


def test_the_escape_helper_escapes_every_html_special_character():
    m = re.search(r"\n  const esc = (.*?);\n", _asset("report.js"))
    assert m, "report.js has no esc helper"
    replaced = re.search(r"\.replace\(/\[([^\]]*)\]/g", m.group(1))
    assert replaced and set(replaced.group(1)) == set("&<>\"'"), "esc does not replace all five"
    entities = {ch: e for _, ch, e in re.findall(r"""(["'])(.)\1: "(&#?\w+;)\"""", m.group(1))}
    assert entities == {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}


@pytest.mark.parametrize(
    "sink, chosen",
    [
        ("hero", r"\bm\."),  # the sample, the index, the timestamp, the version, the warnings
        ("geneTable", r"\bv\b|\br\["),  # gene names, IDs and types, references
        ("refTable", r"\bv\b"),  # the reference names of the per-reference gDNA table
    ],
)
def test_a_string_a_user_or_an_annotation_chose_is_escaped_before_it_becomes_markup(sink, chosen):
    """``innerHTML`` parses what it is given, so a sample name, a path, a gene name or a reference holding ``<`` or
    ``&`` would become markup instead of text."""
    found = [e for e in re.findall(r"\$\{([^}]*)\}", _js_function(sink)) if re.search(chosen, e)]
    assert found, f"{sink} interpolates no chosen string, so this gate reads nothing"
    raw = [e for e in found if not (e.startswith("esc(") or ".map(esc)" in e)]
    assert not raw, f"{sink} puts {raw} into markup unescaped"


def test_the_configuration_panel_escapes_every_group_key_and_value():
    """The configuration holds the BAM, index and output paths the user typed."""
    shown = re.findall(r'H\("\w+", "[\w-]+", (.+?)\)\);', _js_function("config"))
    assert len(shown) == 3, shown
    assert all(a.startswith("esc(") for a in shown), shown


def test_the_capture_answer_is_calibrations_bulk_and_saturation(tmp_path):
    """The capture tile, KPIs and note read calibration's bulk density and saturated level and nothing else: an
    enriched landscape shows its saturation over the bulk, a one-level landscape says nothing is enriched, no
    landscape says so, and a summary without calibration shows no tile."""

    def view(name, **keys):
        return build_view_model(
            load_substrate(_with_calibration(_write_substrate(tmp_path / name), **keys))
        )

    enriched = view("on", gdna_bulk_density=1.2e-4, gdna_capture_saturation=852.0)
    assert enriched["calibration"]["capture"] == {"bulk_density": 1.2e-4, "saturation": 852.0}
    tile = next(v for v in enriched["verdicts"] if v["k"] == "Capture")
    assert (tile["v"], tile["fmt"], tile["s"]) == ("852×", "text", "info")
    assert "0.00012" in tile["n"]
    assert [k["v"] for k in enriched["calibration"]["enrichment_kpis"]] == [1.2e-4, 852.0]

    one = view("off", gdna_bulk_density=6e-3, gdna_capture_saturation=1.0)
    tile = next(v for v in one["verdicts"] if v["k"] == "Capture")
    assert (tile["v"], tile["fmt"]) == ("1×", "text")
    assert "one gDNA level" in tile["n"]

    none = view("none", gdna_bulk_density=None, gdna_capture_saturation=1.0)
    assert none["calibration"]["capture"] == {"bulk_density": None, "saturation": 1.0}
    tile = next(v for v in none["verdicts"] if v["k"] == "Capture")
    assert (tile["v"], tile["fmt"]) == ("None", "text")
    assert [k["v"] for k in none["calibration"]["enrichment_kpis"]] == ["None"]

    absent = build_view_model(load_substrate(_write_substrate(tmp_path / "absent")))
    assert absent["calibration"]["capture"] is None
    assert not any(v["k"] == "Capture" for v in absent["verdicts"])
    assert absent["calibration"]["enrichment_kpis"] == []


def test_the_capture_note_never_says_a_mode_was_located(tmp_path):
    """The reader locates no mode, so the note may not say it did — in any of its three cases — and it names
    the bulk the efficiencies are read against."""
    note = _js_function("captureNote")
    assert "mode was located" not in note and "located" not in note
    assert "bulk" in note


def test_the_capture_note_tells_an_older_summary_from_one_with_no_calibration(tmp_path):
    """A 0.7.1 or schema-2 summary has a calibration block but not calibration's capture answer: the note says
    the summary predates that answer, not that the run carries no calibration."""
    older = build_view_model(
        load_substrate(_with_calibration(_write_substrate(tmp_path / "older")))
    )
    assert older["calibration"]["capture"] is None
    assert older["calibration"]["calibrated"] is True
    bare = build_view_model(load_substrate(_write_substrate(tmp_path / "bare")))
    assert bare["calibration"]["calibrated"] is False
    note = _js_function("captureNote")
    assert "M.calibration.calibrated" in note and "predates" in note


def test_the_capture_note_reads_only_keys_the_view_model_writes(tmp_path):
    """A key the front end reads but the model never writes renders as ``NaN`` or ``undefined``, and
    nothing fails, so every key the note reads off the capture answer must be one the model writes."""
    d = _with_calibration(
        _write_substrate(tmp_path / "run"), gdna_bulk_density=1.2e-4, gdna_capture_saturation=852.0
    )
    written = set(build_view_model(load_substrate(d))["calibration"]["capture"])
    read = set(re.findall(r"\bc\.(\w+)", _js_function("captureNote")))
    assert read, "the capture note reads no key of the capture answer"
    assert read <= written, (
        f"report.js reads {sorted(read - written)}, which the model never writes"
    )


@pytest.mark.parametrize("bulk_density", [1.2e-4, None])
def test_the_report_shows_calibrations_rna_sense_fraction_whatever_the_capture_answer(
    tmp_path, bulk_density
):
    """The sense fraction is a library scalar, not a capture number, so it is shown whatever calibration says
    about capture."""
    d = _with_calibration(
        _write_substrate(tmp_path / "run"),
        rna_sense_frac=0.973,
        gdna_bulk_density=bulk_density,
        gdna_capture_saturation=1.0 if bulk_density is None else 852.0,
    )
    kpis = build_view_model(load_substrate(d))["calibration"]["density_kpis"]
    assert {"l": "RNA sense", "v": 0.973, "fmt": "float3"} in kpis


def test_every_format_tag_the_view_model_emits_is_one_the_front_end_formats(tmp_path):
    """``fmtValue`` falls through to the raw value on a tag it does not know, so a tag the model emits
    must be one of its cases (``text`` values are strings, which it passes through)."""
    d = _with_calibration(
        _write_substrate(tmp_path / "run"), gdna_bulk_density=1.2e-4, gdna_capture_saturation=852.0
    )
    _enriched_track().to_feather(d / "calibration_track.feather")
    vm = build_view_model(load_substrate(d))

    def tags(node):
        if isinstance(node, dict):
            if "fmt" in node:
                yield node["fmt"]
            for value in node.values():
                yield from tags(value)
        elif isinstance(node, list):
            for value in node:
                yield from tags(value)

    formatted = set(re.findall(r'case "(\w+)":', _js_function("fmtValue"))) | {"text"}
    emitted = set(tags(vm))
    assert emitted and emitted <= formatted, f"unformatted tags: {sorted(emitted - formatted)}"


def test_the_pools_book_spliced_intergenic_fragments_as_rna(tmp_path):
    """A spliced intergenic fragment is RNA no annotated transcript holds: the intergenic gDNA pool is the
    unspliced remainder, and the RNA share counts it."""
    d = _write_substrate(tmp_path / "run")
    path = d / "summary.json"
    summary = json.loads(path.read_text())
    summary["quantification"].update(
        intergenic_spliced=10, intergenic_rna_fraction=0.015, gdna_fraction=0.125
    )
    path.write_text(json.dumps(summary))
    quant = build_view_model(load_substrate(d))["quant"]
    pools = {p["label"]: p["value"] for p in quant["pools"]}
    assert pools["Intergenic gDNA"] == 60
    assert pools["Intergenic RNA"] == 10
    q = summary["quantification"]
    assert quant["rna_share"] == pytest.approx(
        q["mrna_fraction"] + q["nrna_fraction"] + q["intergenic_rna_fraction"]
    )


def test_build_report_self_contained(tmp_path):
    d = _write_substrate(tmp_path / "run")
    out = build_report(d)
    assert out.exists() and out.name == "report.html"
    html = out.read_text()

    # self-contained: no external scripts/styles/fetchable URLs
    assert not re.search(r"<script[^>]*\bsrc=", html)
    assert not re.search(r"<link[^>]*\bhref=", html)
    assert not re.search(r"(?:src|href)\s*=\s*[\"']https?://", html)

    # valid document + embedded payload that round-trips through JSON
    assert html.lstrip().startswith("<!doctype html>")
    payload = _payload(html)
    assert payload["model"]["meta"]["sample"] == "SampleX"
    assert set(payload["charts"]) == {"overlay", "small_multiples"}  # no track in this substrate


@pytest.mark.skipif(not _HAS_VEGA, reason="vl-convert-python not installed")
def test_build_report_inlines_vega_runtime(tmp_path):
    out = build_report(_write_substrate(tmp_path / "run"))
    html = out.read_text()
    assert "window.vegaEmbed" in html  # runtime bundled inline, offline-ready


@pytest.mark.skipif(not _HAS_VEGA, reason="vl-convert-python not installed")
def test_every_chart_the_report_builds_compiles_and_reads_every_mark_property(tmp_path):
    """Vega-Lite refuses a mark type it does not know but silently drops a mark property it does not
    know, so compiling is half the gate: removing any property of any mark definition must also change
    the compiled Vega. The charts are the ones the built page embeds, and every chart container on the
    page must carry one, so a chart added to the report is gated too."""
    import copy

    import vl_convert as vlc

    d = _write_substrate(tmp_path / "run")
    _enriched_track().to_feather(d / "calibration_track.feather")
    html = build_report(d).read_text()
    charts = _payload(html)["charts"]
    assert set(charts) == set(re.findall(r'id="vega-(\w+)"', html))

    for key, spec in charts.items():
        compiled = vlc.vegalite_to_vega(spec)
        for path in _mark_paths(spec):
            for prop in _at(spec, path):
                if prop == "type":
                    continue
                pruned = copy.deepcopy(spec)
                del _at(pruned, path)[prop]
                read = vlc.vegalite_to_vega(pruned) != compiled
                assert read, (
                    f"chart {key!r}: Vega-Lite ignores the mark property {prop!r} at {path}"
                )


def test_build_report_custom_output_path(tmp_path):
    d = _write_substrate(tmp_path / "run")
    dest = tmp_path / "reports" / "sample_x.html"
    out = build_report(d, out_path=dest, title="My QC")
    assert out == dest and dest.exists()
    assert "<title>My QC</title>" in dest.read_text()
    # the title is text, never markup
    build_report(d, out_path=dest, title="A & B </title><script>x</script>")
    assert (
        "<title>A &amp; B &lt;/title&gt;&lt;script&gt;x&lt;/script&gt;</title>" in dest.read_text()
    )
