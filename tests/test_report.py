"""Tests for the `rigel report` HTML builder (Phase 1).

Fast, pipeline-free: they synthesize a minimal-but-realistic report substrate
(v2 ``summary.json`` + companion feathers) in a temp dir and exercise the loader,
view model, and full HTML build. Vega-specific assertions are conditional on
``vl-convert-python`` being installed so the suite passes with the ``[dev]``
extra alone.
"""

import importlib.util
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import numpy as np

from rigel.calibration.diagnostics import CalibrationDiagnostics
from rigel.report.build import build_report
from rigel.report.capture import capture_kde_from_track
from rigel.report.model import build_view_model
from rigel.report.specs import build_charts, build_fl_specs, capture_kde_spec, genome_track_spec
from rigel.report.substrate import SubstrateError, load_substrate

_HAS_VEGA = importlib.util.find_spec("vl_convert") is not None


def _write_substrate(d: Path) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    summary = {
        "schema_version": 2,
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


def test_load_substrate_missing(tmp_path):
    with pytest.raises(SubstrateError):
        load_substrate(tmp_path / "does_not_exist")
    (tmp_path / "empty").mkdir()
    with pytest.raises(SubstrateError):
        load_substrate(tmp_path / "empty")  # no summary.json


def test_view_model_shape(tmp_path):
    d = _write_substrate(tmp_path / "run")
    sub = load_substrate(d)
    assert sub.schema_version == 2
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
            "gdna_mass": [1.0, 8.0, 2.0, 1.0, 3.0, 9.0],
            "rna_mass": [9.0, 2.0, 8.0, 9.0, 7.0, 1.0],
            "gdna_density": [0.01, 0.09, 0.02, 0.01, 0.03, 0.10],
            "gdna_frac": [0.1, 0.8, 0.2, 0.1, 0.3, 0.9],
        }
    )
    spec = genome_track_spec(track)
    assert spec is not None
    refs = {row["ref"] for row in spec["data"]["values"]}
    assert refs == {"chr1", "chr2"}
    # build_charts merges genome in when a track is present
    stub = SimpleNamespace(
        fragment_lengths=None,
        calibration_track=track,
        gdna_density_kde=None,
        gdna_density_nodes=None,
        summary={},
    )
    assert set(build_charts(stub)) == {"genome"}
    empty = SimpleNamespace(
        fragment_lengths=None,
        calibration_track=None,
        gdna_density_kde=None,
        gdna_density_nodes=None,
        summary={},
    )
    assert build_charts(empty) == {}


def test_capture_diagnostics_from_prior_labels_modes():
    # A bimodal density: depleted mode at x=-8, enriched at x=-2.
    x = np.linspace(-10, 0, 256)
    logp = np.log(
        np.exp(-0.5 * ((x + 8) / 0.4) ** 2) + 0.6 * np.exp(-0.5 * ((x + 2) / 0.5) ** 2) + 1e-9
    )
    prior = SimpleNamespace(
        x_grid=x,
        logP_grid=logp,
        bandwidth=0.4,
        n_eff=1234.0,
        train_x=np.array([-8.1, -7.9, -2.1, -1.9]),
        train_w=np.ones(4),
        train_kind=np.array([0, 1, 2, 3]),
        modes=((-8.0, 0.0), (-2.0, -0.5)),  # (x, logP) sorted by height desc
    )
    diag = CalibrationDiagnostics.from_prior(prior)
    assert diag.depleted_mode < diag.enriched_mode
    assert diag.separation_nats == pytest.approx(6.0, abs=1e-6)
    assert diag.enrichment_factor == pytest.approx(np.exp(6.0), rel=1e-6)
    assert diag.n_modes == 2


def test_capture_kde_from_track_mass_weighting_recovers_enrichment():
    # Many low-density (off-target) regions with tiny gDNA mass + a few
    # high-density (on-target) regions carrying large gDNA mass. Equal weight is
    # dominated by the low mode; mass weighting must surface the high mode.
    rng = np.random.default_rng(0)
    low_d = np.exp(rng.normal(-9.0, 0.4, 4000))
    high_d = np.exp(rng.normal(0.0, 0.4, 120))
    dens = np.concatenate([low_d, high_d])
    gmass = np.concatenate([np.full(4000, 0.01), np.full(120, 50.0)])
    track = pd.DataFrame(
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
    cap = capture_kde_from_track(track)
    assert cap is not None
    assert cap["enriched"] is True
    assert cap["enriched_mode_log_rho"] > cap["count_median_log_rho"]
    assert cap["separation_median_nats"] > 3.0
    assert cap["fold_vs_median"] > 3.0
    assert 0.0 < cap["mass_frac_ontarget"] <= 1.0
    # spec renders as an overlay (curve + mode rule + mode text)
    spec = capture_kde_spec(cap)
    assert spec is not None and len(spec["layer"]) == 3
    assert capture_kde_spec(None) is None


def test_capture_kde_from_track_unimodal_when_no_enrichment():
    rng = np.random.default_rng(1)
    dens = np.exp(rng.normal(-9.0, 0.5, 3000))
    track = pd.DataFrame(
        {
            "ref": pd.Categorical(["chr1"] * len(dens)),
            "start": np.arange(len(dens)),
            "end": np.arange(len(dens)) + 50,
            "gdna_mass": np.full(len(dens), 0.02),
            "rna_mass": np.ones(len(dens)),
            "gdna_density": dens,
            "gdna_frac": np.clip(dens, 0, 1),
        }
    )
    cap = capture_kde_from_track(track)
    assert cap is not None and cap["enriched"] is False


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
    m = re.search(r'<script id="rigel-data" type="application/json">(.*?)</script>', html, re.S)
    payload = json.loads(m.group(1).replace("<\\/", "</"))
    assert payload["model"]["meta"]["sample"] == "SampleX"
    assert set(payload["charts"]) == {"overlay", "small_multiples"}  # no track in this substrate


@pytest.mark.skipif(not _HAS_VEGA, reason="vl-convert-python not installed")
def test_build_report_inlines_vega_runtime(tmp_path):
    out = build_report(_write_substrate(tmp_path / "run"))
    html = out.read_text()
    assert "window.vegaEmbed" in html  # runtime bundled inline, offline-ready


def test_build_report_custom_output_path(tmp_path):
    d = _write_substrate(tmp_path / "run")
    dest = tmp_path / "reports" / "sample_x.html"
    out = build_report(d, out_path=dest, title="My QC")
    assert out == dest and dest.exists()
    assert "<title>My QC</title>" in dest.read_text()
