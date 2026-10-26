"""The `rigel report` HTML builder, from a synthesized substrate rather than a pipeline run.

A minimal but realistic substrate — a schema-4 ``summary.json`` and its companion feathers — is written
to a temp directory, and the loader, the view model, the chart specifications and the full HTML build run
against it. The report must be self-contained, inlining its runtime, and must honour a custom output path.
The pools carry the output convention's five names, one spelling each, and the summary's own counts and
fractions; the strandedness tile is the strand model's own verdict; a gene's RNA is its count; every
fragment-length category shows one readable label; and every cut, cap and rounding the report applies is
``CONSTANTS.qc``'s. The capture panel shows calibration's own answer, and the front end reads only the keys,
cells and format tags the view model writes, and escapes every string a user or an annotation chose before it
becomes markup (read off ``report.js``'s source: the suite runs no JavaScript). Every chart the page embeds
must compile and read every mark property it sets. Vega-specific assertions are conditional on
``vl-convert-python``.

The report reads ``summary.json`` only through the keys ``SUMMARY_KEYS`` lists, and tolerantly, so a key
misspelled on the report's side would read its default without a sound — and the synthesized substrate, written
by the same hand, would agree with it. One gate therefore reads a REAL output directory, written by the writer
``rigel quant`` runs from a small simulated library: every listed key, every gene-table column and every
fragment-length category must be one that writer writes.
"""

import argparse
import importlib.util
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import numpy as np

from rigel import cli
from rigel.config import CONSTANTS, BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.report import model as report_model
from rigel.report import specs as report_specs
from rigel.report.build import build_report
from rigel.report.html import _SECTIONS, _asset
from rigel.report.model import GENE_COLUMNS, _reference_table, build_view_model
from rigel.report.specs import (
    FL_CATEGORY_LABELS,
    build_charts,
    build_fl_specs,
    fl_label,
    genome_track_spec,
)
from rigel.report.substrate import SUMMARY_KEYS, SubstrateError, SummaryReader, load_substrate
from rigel.sim import ReadSimConfig, Scenario

_HAS_VEGA = importlib.util.find_spec("vl_convert") is not None

#: The ``summary.json`` schema the synthesized substrate mirrors, and the one the real writer must write.
_SCHEMA = 6

#: The output convention's four pools, as the report names them.
_POOLS = ("Annotated transcripts", "Nascent RNA", "gDNA (EM)", "Intergenic gDNA")


def _write_substrate(d: Path) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    # The pools, as ``rigel quant`` books them: an intergenic fragment is gDNA.
    annotated, nascent, gdna_em, intergenic = 500.0, 60.0, 20.0, 70
    every = annotated + nascent + gdna_em + intergenic

    def share(x):
        return round(x / every, 6)

    summary = {
        "schema_version": _SCHEMA,
        "rigel_version": "0.7.0",
        "timestamp": "2026-07-11 09:42",
        "input": {"bam_file": "/data/SampleX.bam", "index_dir": "/refs/gencode"},
        "configuration": {"em": {"mode": "vbem", "n_threads": 8}, "seed": 42},
        "alignment_stats": {
            "total_reads": 1000,
            "mapped_reads": 950,
            "proper_pairs": 900,
            "duplicate_reads": 60,
            "qc_fail_reads": 5,
        },
        "read_pairs": {
            "total": 750,
            "excluded": 30,
            "withheld": {"unpaired": 5, "unannotated_splice": 10, "off_index": 0, "chimeric": 5},
            "excluded_by_option": 0,
            "placed": {"unique": 620, "multimapping": 80},
        },
        "fragment_stats": {
            "total": 700,
            "genic": 600,
            "intergenic": intergenic,
            "with_annotated_sj": 400,
            "splice": {
                "unspliced": 270,
                "spliced_annotated": 400,
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
            "strand_live": True,
            "diagnostics": {
                "exonic_all_specificity": 0.93,
                "exonic_all_p_r1_sense": 0.07,
                "exonic_all_n_fragments": 600,
                "contamination_gap": 0.05,
            },
        },
        "quantification": {
            "n_transcripts": 120,
            "n_nascent": 30,
            "n_genes": 40,
            "n_loci": 35,
            "count_annotated": annotated,
            "count_nascent": nascent,
            "count_gdna": gdna_em + intergenic,
            "count_gdna_em": gdna_em,
            "fraction_annotated": share(annotated),
            "fraction_nascent": share(nascent),
            "fraction_gdna": share(gdna_em + intergenic),
            "fraction_gdna_em": share(gdna_em),
            "fraction_intergenic": share(intergenic),
        },
        # every category ``rigel quant`` writes, in its order: the anchor, the two modelled pools, and the five
        # pools they are built from (``gdna`` is the two contained pools, ``rna`` the spliced one)
        "fragment_length": {
            cat: {
                "n_observations": n,
                "mean": mean,
                "std": std,
                "median": median,
                "mode": mode,
                "max_size": 1000,
                "overflow_count": overflow,
                "overflow_fraction": round(overflow / n, 6),
            }
            for cat, (n, mean, std, median, mode, overflow) in {
                "global": (700, 200.0, 40.0, 195.0, 190, 0.0),
                "gdna": (125, 175.0, 60.0, 160.0, 150, 1.0),
                "rna": (400, 205.0, 35.0, 200.0, 195, 0.0),
                "pool_gdna_intergenic": (60, 172.0, 58.0, 158.0, 150, 1.0),
                "pool_gdna_intronic": (65, 178.0, 61.0, 162.0, 150, 0.0),
                "pool_gdna_intron_exon": (20, 181.0, 55.0, 170.0, 160, 0.0),
                "pool_gdna_intergenic_exon": (10, 179.0, 57.0, 168.0, 155, 0.0),
                "pool_rna_spliced": (400, 205.0, 35.0, 200.0, 195, 0.0),
            }.items()
        },
    }
    (d / "summary.json").write_text(json.dumps(summary))

    # fragment_lengths.feather — a gaussian-ish bump per category
    rows = [
        (cat, length, max(0.0, 50.0 - abs(length - int(stats["mean"]))))
        for cat, stats in summary["fragment_length"].items()
        for length in range(int(stats["mean"]) - 40, int(stats["mean"]) + 41, 5)
    ]
    pd.DataFrame(rows, columns=["category", "length", "count"]).astype(
        {"category": "category", "length": "int32", "count": "float64"}
    ).to_feather(d / "fragment_lengths.feather")

    pd.DataFrame(
        {
            "gene_id": [f"ENSG{i:05d}" for i in range(5)],
            "gene_name": ["GAPDH", "ACTB", "TP53", "MYC", "MALAT1"],
            "count": [4200.37, 3800.0, 120.0, 300.0, 900.0],
            "count_spliced": [4100.0, 3700.0, 110.0, 280.0, 40.0],
            "tpm": [42000.456, 38000.0, 1200.0, 2900.0, 9000.0],
            "n_transcripts": [3, 4, 6, 2, 1],
        }
    ).to_feather(d / "gene_quant.feather")
    return d


def _enriched_track() -> pd.DataFrame:
    """A capture-like gDNA track: many low-density off-target regions carrying few gDNA fragments, and a few
    high-density on-target regions carrying most of them."""
    rng = np.random.default_rng(0)
    low_d = np.exp(rng.normal(-9.0, 0.4, 4000))
    high_d = np.exp(rng.normal(0.0, 0.4, 120))
    dens = np.concatenate([low_d, high_d])
    gdna = np.concatenate([np.full(4000, 0.01), np.full(120, 50.0)])
    return pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1"] * len(dens)),
            "start": np.arange(len(dens)) * 100,
            "end": np.arange(len(dens)) * 100 + 50,
            "count_gdna": gdna,
            "count_rna": np.ones(len(dens)),
            "gdna_density": dens,
            "fraction_gdna": np.clip(dens, 0, 1),
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
    assert sub.schema_version == _SCHEMA
    assert sub.sample_name == "SampleX"
    vm = build_view_model(sub)

    assert len(vm["verdicts"]) == 4
    # formatting boundary: KPIs / verdicts carry raw values + a fmt tag (JS formats)
    for kpi in vm["alignment"]["kpis"]:
        assert "fmt" in kpi and "u" not in kpi
    numeric = next(k for k in vm["alignment"]["kpis"] if k["fmt"] in ("pct", "count"))
    assert isinstance(numeric["v"], (int, float))
    assert all("fmt" in v for v in vm["verdicts"])
    # the read fate is every read pair, each counted once
    assert sum(seg["value"] for seg in vm["alignment"]["fate"]) == 750
    # the splice bar: the four splice types, one spliced (every junction annotated)
    labels = [s["label"] for s in vm["fragments"]["splice"]]
    assert labels == ["Unspliced", "Annotated", "Implicit", "Artifact"]
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


def test_genome_track_spec_bins_per_ref(monkeypatch):
    track = pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1"] * 4 + ["chr2"] * 2),
            "start": [0, 1000, 2000, 3000, 0, 5000],
            "end": [1000, 2000, 3000, 4000, 5000, 10000],
            "count_gdna": [1.0, 8.0, 2.0, 1.0, 30.0, 90.0],  # chr1 sum 12, chr2 sum 120
            "count_rna": [9.0, 2.0, 8.0, 9.0, 7.0, 1.0],
            "gdna_density": [0.01, 0.09, 0.02, 0.01, 0.03, 0.10],
            "fraction_gdna": [0.1, 0.8, 0.2, 0.1, 0.3, 0.9],
        }
    )
    spec = genome_track_spec(track)
    assert spec is not None
    refs = {row["ref"] for row in spec["data"]["values"]}
    assert refs == {"chr1", "chr2"}
    # log + independent y per facet (so a high-density ref can't flatten others)
    assert spec["spec"]["encoding"]["y"]["scale"]["type"] == "log"
    assert spec["resolve"]["scale"]["y"] == "independent"
    # top-N ranks by gDNA count; chr2 (120) outranks chr1 (12)
    _charts_with(monkeypatch, genome_track_refs=1)
    top1 = genome_track_spec(track)
    assert {row["ref"] for row in top1["data"]["values"]} == {"chr2"}
    monkeypatch.undo()
    # reference table: every reference, sorted by gDNA count desc
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


def _with(d: Path, section: str, **keys) -> Path:
    """Give one section of the substrate's ``summary.json`` these keys."""
    path = d / "summary.json"
    summary = json.loads(path.read_text())
    summary[section].update(keys)
    path.write_text(json.dumps(summary))
    return d


def _view(d: Path, monkeypatch, **qc) -> dict:
    """The view model with these ``CONSTANTS.qc`` fields changed, the way an experiment varies a constant."""
    constants = CONSTANTS
    for name, value in qc.items():
        constants = constants.replaced(f"qc.{name}", value)
    monkeypatch.setattr(report_model, "CONSTANTS", constants)
    return build_view_model(load_substrate(d))


def _charts_with(monkeypatch, **qc) -> None:
    """Have the chart builders read these ``CONSTANTS.qc`` fields changed, the way an experiment varies a
    constant."""
    constants = CONSTANTS
    for name, value in qc.items():
        constants = constants.replaced(f"qc.{name}", value)
    monkeypatch.setattr(report_specs, "CONSTANTS", constants)


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


@pytest.mark.parametrize(
    "table, rows",
    [
        ("flTable", lambda vm: vm["fl"]["table"]),
        ("geneTable", lambda vm: vm["genes"]["rows"]),
        ("refTable", lambda vm: vm["calibration"]["ref_table"]),
    ],
)
def test_each_table_reads_every_cell_the_view_model_writes_and_no_other(tmp_path, table, rows):
    """A row position the front end reads past the end of a row renders ``undefined``, and one it never reads is a
    column nobody sees; nothing fails in either case, so the positions each table reads must be exactly the ones
    its rows carry."""
    d = _with_calibration(_write_substrate(tmp_path / "run"))
    _enriched_track().to_feather(d / "calibration_track.feather")
    written = rows(build_view_model(load_substrate(d)))
    assert written, f"{table} has no rows in this substrate, so this gate reads nothing"
    (width,) = {len(r) for r in written}
    js = _js_function(table)
    read = {int(i) for i in re.findall(r"\bidx: (\d+)", js) + re.findall(r"\br\[(\d+)\]", js)}
    assert read == set(range(width)), f"{table} reads {sorted(read)} of a {width}-cell row"


@pytest.mark.parametrize("bulk_density", [1.2e-4, None])
def test_the_report_shows_calibrations_rna_sense_fraction_whatever_the_capture_answer(
    tmp_path, bulk_density
):
    """The sense fraction is a library scalar, not a capture number, so it is shown whatever calibration says
    about capture."""
    d = _with_calibration(
        _write_substrate(tmp_path / "run"),
        fraction_rna_sense=0.973,
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


def test_each_pool_is_named_once_and_shows_the_summarys_own_count_and_fraction(tmp_path):
    """The bar, its legend, the donut, the table and the KPIs name each pool the same way, in the output
    convention's words, and show the count and fraction ``summary.json`` holds for it. The words the convention
    retired name nothing the report shows."""
    d = _write_substrate(tmp_path / "run")
    summary = json.loads((d / "summary.json").read_text())
    q, intergenic = summary["quantification"], summary["fragment_stats"]["intergenic"]
    quant = build_view_model(load_substrate(d))["quant"]
    expected = {
        "Annotated transcripts": (q["count_annotated"], q["fraction_annotated"]),
        "Nascent RNA": (q["count_nascent"], q["fraction_nascent"]),
        "gDNA (EM)": (q["count_gdna_em"], q["fraction_gdna_em"]),
        "Intergenic gDNA": (intergenic, q["fraction_intergenic"]),
    }
    assert set(expected) == set(_POOLS)
    assert {p["label"]: p["value"] for p in quant["pools"]} == {
        k: n for k, (n, _) in expected.items()
    }
    assert {row[0]: (row[1], row[2]) for row in quant["table"]} == expected
    kpis = {k["l"]: k["v"] for k in quant["kpis"]}
    assert kpis["Annotated transcripts"] == q["count_annotated"]
    assert kpis["Nascent RNA"] == q["count_nascent"]
    assert kpis["Transcripts"] == q["n_transcripts"]  # annotated transcripts only, in schema 4
    shown = json.dumps(quant) + _SECTIONS + _asset("report.js")
    assert not re.findall(r"\b(?:[Mm]ature|mRNA|nRNA)\b", shown)


def test_the_rna_share_is_the_annotated_and_nascent_pools(tmp_path):
    """Every intergenic fragment is gDNA, so the RNA share is the annotated and nascent fractions."""
    d = _write_substrate(tmp_path / "run")
    quant = build_view_model(load_substrate(d))["quant"]
    q = json.loads((d / "summary.json").read_text())["quantification"]
    assert quant["rna_share"] == pytest.approx(q["fraction_annotated"] + q["fraction_nascent"])


@pytest.mark.parametrize(
    "live, read1_sense, specificity, shown, state",
    [
        (True, False, 0.98, "Reverse (R1−)", "good"),
        (True, True, 0.98, "Forward (R1+)", "good"),
        (True, False, 0.55, "Reverse (R1−)", "good"),
        (False, False, 0.98, "Unstranded", "warning"),
    ],
)
def test_the_strandedness_tile_is_the_strand_models_own_verdict(
    tmp_path, live, read1_sense, specificity, shown, state
):
    """The tile shows ``strand_live``, the verdict calibration's strand channel reads, and draws no cut of its own
    on the strand specificity. The two are set apart here so the test can tell which one the tile read: a live
    channel at a specificity of 0.55 reads stranded, a dead one at 0.98 unstranded. A stranded tile names the
    protocol's orientation."""
    d = _with(
        _write_substrate(tmp_path / "run"),
        "strand_model",
        strand_live=live,
        read1_sense=read1_sense,
        strand_specificity=specificity,
    )
    tile = next(
        v for v in build_view_model(load_substrate(d))["verdicts"] if v["k"] == "Strandedness"
    )
    assert (tile["v"], tile["s"]) == (shown, state)


def test_a_genes_rna_is_its_count(tmp_path):
    """The gene table shows a gene's annotated ``count``, and no mature or nascent column: under the output
    convention a single-exon transcript is annotated, and the nascent spans belong to no gene."""
    d = _write_substrate(tmp_path / "run")
    rows = build_view_model(load_substrate(d))["genes"]["rows"]
    counts = pd.read_feather(d / "gene_quant.feather").set_index("gene_name")["count"]
    assert {r[0]: r[6] for r in rows} == {
        g: round(c, CONSTANTS.qc.count_decimals) for g, c in counts.items()
    }
    assert not re.search(r"mature|nascent", _js_function("geneTable"))


def test_every_fragment_length_category_shows_one_readable_label(tmp_path):
    """The summary table, the overlay's legend and the per-category facets name each category by its one label,
    every category ``rigel quant`` writes has a label of its own, and a category the labels do not know shows as
    itself."""
    d = _write_substrate(tmp_path / "run")
    sub = load_substrate(d)
    cats = list(sub.summary["fragment_length"])
    labels = [fl_label(c) for c in cats]
    assert len(set(labels)) == len(labels) and not set(labels) & set(cats)
    assert fl_label("not_a_category") == "not_a_category"

    assert [row[-1] for row in build_view_model(sub)["fl"]["table"]] == labels
    specs = build_fl_specs(sub.fragment_lengths)
    overlay = {r["category"] for r in specs["overlay"]["data"]["values"]}
    assert overlay == {fl_label("gdna"), fl_label("rna")}
    assert set(specs["overlay"]["encoding"]["color"]["sort"]) == overlay
    facets = specs["small_multiples"]["facet"]["sort"]
    assert facets == [fl_label(c) for c in cats if c != "global"]
    assert {r["category"] for r in specs["small_multiples"]["data"]["values"]} == set(facets)


def test_the_tiles_cut_where_the_qc_constants_say(tmp_path, monkeypatch):
    """Each headline tile's cuts are ``CONSTANTS.qc``'s display thresholds, so moving one moves the tile. The
    substrate maps 95 % of its records and holds 13.8 % gDNA and 85.7 % genic fragments."""
    d = _write_substrate(tmp_path / "run")

    def states(**qc):
        vm = _view(d, monkeypatch, **qc)
        return tuple(
            next(v["s"] for v in vm["verdicts"] if v["k"] == k)
            for k in ("Mapping rate", "gDNA contamination", "Usable fragments")
        )

    assert states() == ("good", "warning", "good")
    assert states(mapping_rate_good=0.96, gdna_fraction_good=0.2, genic_fraction_good=0.9) == (
        "warning",
        "good",
        "warning",
    )
    assert states(
        mapping_rate_good=0.99, mapping_rate_warning=0.96, gdna_fraction_warning=0.12
    ) == ("critical", "critical", "good")


def test_the_gdna_region_count_reads_the_qc_cut_and_names_it(tmp_path, monkeypatch):
    """The density panel counts the track's regions above ``CONSTANTS.qc``'s gDNA-fraction cut, and its label
    states that cut."""
    d = _with_calibration(_write_substrate(tmp_path / "run"))
    track = _enriched_track()
    track.to_feather(d / "calibration_track.feather")
    for cut, label in ((None, "Regions >50% gDNA"), (0.99, "Regions >99% gDNA")):
        moved = {} if cut is None else {"region_gdna_fraction": cut}
        kpis = {
            k["l"]: k["v"] for k in _view(d, monkeypatch, **moved)["calibration"]["density_kpis"]
        }
        threshold = CONSTANTS.qc.region_gdna_fraction if cut is None else cut
        assert kpis[label] == int((track["fraction_gdna"] > threshold).sum())


def test_the_gene_table_embeds_at_most_the_qc_cap(tmp_path, monkeypatch):
    """Expressed genes, the highest TPM first, up to ``CONSTANTS.qc``'s cap."""
    d = _write_substrate(tmp_path / "run")
    genes = _view(d, monkeypatch)["genes"]
    assert (genes["shown"], genes["total"], genes["truncated"]) == (5, 5, False)
    genes = _view(d, monkeypatch, gene_table_rows=2)["genes"]
    assert (genes["shown"], genes["total"], genes["truncated"]) == (2, 5, True)
    assert [r[0] for r in genes["rows"]] == ["GAPDH", "ACTB"]


def test_every_cell_the_report_rounds_takes_its_decimals_from_qc(tmp_path, monkeypatch):
    """A gene's TPM and count, a reference's gDNA count and fraction, and a category's overflow percentage are
    rounded to ``CONSTANTS.qc``'s decimals."""
    d = _with_calibration(_write_substrate(tmp_path / "run"))
    pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1", "chr1"]),
            "start": [0, 1000],
            "end": [1000, 2000],
            "count_gdna": [10.26, 2.0],
            "count_rna": [5.0, 3.0],
            "gdna_density": [0.01, 0.002],
            "fraction_gdna": [0.67, 0.4],
        }
    ).to_feather(d / "calibration_track.feather")

    def cells(**qc):
        vm = _view(d, monkeypatch, **qc)
        gapdh, chr1 = vm["genes"]["rows"][0], vm["calibration"]["ref_table"][0]
        overflow = {row[0]: row[6] for row in vm["fl"]["table"]}["gdna"]
        return gapdh[5], gapdh[6], chr1[2], chr1[4], overflow

    assert cells() == (42000.46, 4200.4, 12.3, 0.6051, 0.8)
    moved = cells(
        tpm_decimals=0, count_decimals=0, fraction_decimals=1, overflow_percent_decimals=0
    )
    assert moved == (42000.0, 4200.0, 12.0, 0.6, 1.0)


def test_the_charts_bin_cap_and_round_where_the_qc_constants_say(monkeypatch):
    """The fragment-length overlay's bin width, the genome track's bins per reference, its reference cap and the
    rounding of its positions are ``CONSTANTS.qc``'s, so moving one moves the chart. chr1 spans 10 Mb in a
    thousand regions and holds the most gDNA, chr2 the next most."""
    fl = pd.DataFrame(
        {"category": ["rna"] * 20, "length": np.arange(100, 120), "count": np.ones(20)}
    )
    n = {"chr1": 1000, "chr2": 10, "chr3": 10}
    starts = np.concatenate([np.arange(k) * 10_000 for k in n.values()])
    track = pd.DataFrame(
        {
            "ref": pd.Categorical([ref for ref, k in n.items() for _ in range(k)]),
            "start": starts,
            "end": starts + 10_000,
            "count_gdna": np.repeat([1.0, 5.0, 1.0], list(n.values())),
            "count_rna": np.ones(len(starts)),
            "gdna_density": np.full(len(starts), 1e-4),
            "fraction_gdna": np.full(len(starts), 0.5),
        }
    )

    def charts():
        overlay = sorted({r["length"] for r in build_fl_specs(fl)["overlay"]["data"]["values"]})
        rows = genome_track_spec(track)["data"]["values"]
        chr1 = [r["pos_mb"] for r in rows if r["ref"] == "chr1"]
        return overlay, {r["ref"] for r in rows}, chr1

    overlay, refs, chr1 = charts()
    assert overlay == [100, 105, 110, 115]
    assert refs == {"chr1", "chr2", "chr3"}
    assert chr1 == [round((b + 0.5) * 10 / 160, 4) for b in range(160)]
    _charts_with(
        monkeypatch,
        fl_overlay_bin_bp=10,
        genome_track_refs=2,
        genome_track_bins=7,
        position_mb_decimals=1,
    )
    overlay, refs, chr1 = charts()
    assert overlay == [100, 110]
    assert refs == {"chr1", "chr2"}
    assert chr1 == [0.7, 2.1, 3.6, 5.0, 6.4, 7.9, 9.3]


def test_the_genome_caption_states_the_cap_the_view_model_passes(tmp_path, monkeypatch):
    """The caption over the genome track states its reference cap. The front end reads it off the view model,
    which writes ``CONSTANTS.qc``'s, so the caption and the track cannot disagree and the number is written down
    once: the caption's text holds no number of its own."""
    d = _with_calibration(_write_substrate(tmp_path / "run"))
    calibration = build_view_model(load_substrate(d))["calibration"]
    assert calibration["genome_track_refs"] == CONSTANTS.qc.genome_track_refs
    assert _view(d, monkeypatch, genome_track_refs=3)["calibration"]["genome_track_refs"] == 3
    js = _js_function("refTable")
    assert "const cal = M.calibration" in js
    read = set(re.findall(r"\bcal\.(\w+)", js))
    assert "genome_track_refs" in read and read <= set(calibration), sorted(read - set(calibration))
    (caption,) = re.findall(r"cap\.textContent = `(.*?)`;", js)
    assert "cal.genome_track_refs" in caption
    assert not re.search(r"\d", caption), f"the caption restates a number: {caption}"


def test_the_summary_reader_hands_out_only_the_keys_it_lists():
    """A listed key reads the summary's own value or its listed default, a missing or ``null`` section reads as
    all defaults, and a section or key the list does not name raises: the report cannot read past the list the
    gate below holds to the writer."""
    reader = SummaryReader(
        {"quantification": {"count_annotated": 5.0, "unlisted": 1}, "calibration": None}
    )
    q = reader["quantification"]
    assert q["count_annotated"] == 5.0
    assert q["count_nascent"] == SUMMARY_KEYS["quantification"]["count_nascent"]
    assert set(q) == set(SUMMARY_KEYS["quantification"])
    assert reader["calibration"] == SUMMARY_KEYS["calibration"]
    assert reader["strand_model.diagnostics"] == SUMMARY_KEYS["strand_model.diagnostics"]
    assert reader["fragment_length.*"] == {}
    assert reader.holds("quantification.count_annotated")
    assert not reader.holds("calibration.gdna_bulk_density")
    for unlisted in (
        lambda: reader["not_a_section"],
        lambda: q["unlisted"],
        lambda: reader.holds("calibration.unlisted"),
    ):
        with pytest.raises(KeyError):
            unlisted()


@pytest.fixture(scope="module")
def quant_output(tmp_path_factory) -> Path:
    """A real ``rigel quant`` output directory: a small stranded library with gDNA and nascent RNA, simulated,
    run through the pipeline and written by the writer ``rigel quant`` runs (``cli._write_quant_outputs``)."""
    work = tmp_path_factory.mktemp("quant")
    sc = Scenario("report", genome_length=12000, seed=42, work_dir=work / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(500, 900), (1500, 1900), (2500, 2900)], "abundance": 80},
            {"t_id": "t2", "exons": [(500, 900), (2500, 2900)], "abundance": 40},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t3", "exons": [(6000, 6400), (7000, 7400)], "abundance": 50}])
    sim = ReadSimConfig(
        frag_mean=200,
        frag_std=40,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.9,
        seed=42,
    )
    oracle = sc.build_oracle(
        n_fragments=2000, sim_config=sim, gdna_fraction=0.2, nrna_abundance=15.0
    )
    config = PipelineConfig(
        em=EMConfig(seed=42, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    result = run_pipeline(oracle.bam_path, oracle.index, config=config)
    out = work / "out"
    out.mkdir()
    args = argparse.Namespace(
        bam_file=str(oracle.bam_path),
        index_dir=str(oracle.index_dir),
        config=None,
        tsv=False,
        emit_locus_stats=False,
    )
    cli._write_quant_outputs(result, oracle.index, out, args)
    return out


def _held(summary: dict, path: str) -> bool:
    """Whether a summary holds a dotted key path, read off the raw JSON rather than through the report's reader."""
    node = summary
    for step in path.split("."):
        if not isinstance(node, dict) or step not in node:
            return False
        node = node[step]
    return True


def test_every_key_the_report_reads_is_one_the_writer_writes(quant_output):
    """Every ``SUMMARY_KEYS`` key, at every fragment-length category, is one the writer wrote; so is every
    ``gene_quant`` column the gene table shows; the fragment-length categories it lists are the ones the labels
    name (its histogram table holds the ones with observations); and the whole report builds from what it
    wrote."""
    summary = json.loads((quant_output / "summary.json").read_text())
    assert summary["schema_version"] == _SCHEMA, (
        "the writer has moved past the schema the report reads"
    )
    categories = list(summary["fragment_length"])
    paths = [
        f"{section}.{key}" if section else key
        for name, keys in SUMMARY_KEYS.items()
        for section in ([name.replace("*", c) for c in categories] if "*" in name else [name])
        for key in keys
    ]
    assert categories and paths
    missing = [p for p in paths if not _held(summary, p)]
    assert not missing, (
        f"the report reads {missing}, which the schema-{_SCHEMA} writer does not write"
    )

    columns = set(pd.read_feather(quant_output / "gene_quant.feather").columns)
    assert set(GENE_COLUMNS) <= columns, f"gene_quant has no {sorted(set(GENE_COLUMNS) - columns)}"
    assert set(categories) == set(FL_CATEGORY_LABELS)
    histograms = set(pd.read_feather(quant_output / "fragment_lengths.feather")["category"])
    assert histograms and histograms <= set(categories)
    assert build_report(quant_output).exists()


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
