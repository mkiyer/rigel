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
from rigel.report.model import _reference_table, build_view_model
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
        gdna_density_kde=None,
        gdna_density_regions=None,
        summary={},
    )
    assert set(build_charts(stub)) == {"genome"}
    empty = SimpleNamespace(
        fragment_lengths=None,
        calibration_track=None,
        gdna_density_kde=None,
        gdna_density_regions=None,
        summary={},
    )
    assert build_charts(empty) == {}


def test_capture_diagnostics_from_abundance_landscape_labels_modes():
    """⭐ The QC panel is built from the `AbundanceLandscape`'s CENSUS, not from the top two maxima of
    a curve. On the bimodal fixture the depleted mode must be the one the pooled intergenic ANCHOR
    falls in — an independent measurement — and the separation must be the census's own mode ratio.

    ⛔ This replaces `from_prior`, which read a `DensityNPMLE`. The census is strictly more
    informative: basins rather than the two tallest peaks, a real `n_train`, and a REAL rug (the
    npmle carried no training points at all, so the report's rug was always empty)."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent / "calibration"))
    from test_abundance_landscape import bimodal_parts, parts

    from rigel.calibration.abundance_landscape import fit_abundance_landscape

    counts, lengths, sig, rho_lo, rho_hi = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    assert al is not None and al.enriched is not None

    diag = CalibrationDiagnostics.from_abundance_landscape(al)
    assert diag.depleted_mode < diag.enriched_mode
    # the census's own numbers, not a re-derivation
    assert diag.depleted_mode == pytest.approx(al.depleted.log_rho, rel=0, abs=0)
    assert diag.enriched_mode == pytest.approx(al.enriched.log_rho, rel=0, abs=0)
    assert diag.separation_nats == pytest.approx(np.log(al.span_R), rel=1e-12)
    assert diag.enrichment_factor == pytest.approx(al.span_R, rel=1e-12)
    assert diag.n_modes == len(al.modes)
    # the true separation of the fixture's two populations, recovered
    assert diag.separation_nats == pytest.approx(np.log(rho_hi / rho_lo), rel=0.25)
    # ⭐ a REAL rug — the thing the npmle could never supply
    assert diag.rug_log_rho.size > 0
    assert diag.rug_log_rho.size == diag.rug_kind.size
    assert set(np.unique(diag.rug_kind)) <= {0, 1, 2, 3}
    assert {0, 2} <= set(np.unique(diag.rug_kind))  # the fixture has intergenic and exon regions
    # n_eff is the training count, and the bandwidth is the smoothing ACTUALLY in force (the grid step)
    assert diag.n_eff == float(al.n_train)
    step = float(al.landscape.log_rho[1] - al.landscape.log_rho[0]) / np.log(10.0)
    assert diag.bandwidth == pytest.approx(step, rel=1e-12)


def test_capture_diagnostics_is_unimodal_safe():
    """A unimodal field has no enriched mode: the separation is 0 and the factor exactly 1, never
    None and never a crash — the report renders the panel with one labelled mode."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent / "calibration"))
    from test_abundance_landscape import parts, unimodal_parts

    from rigel.calibration.abundance_landscape import fit_abundance_landscape

    c, ln, sg = unimodal_parts()
    al = fit_abundance_landscape(*parts(c, ln, sg))
    assert al.enriched is None
    diag = CalibrationDiagnostics.from_abundance_landscape(al)
    assert diag.separation_nats == 0.0
    assert diag.enrichment_factor == 1.0
    assert diag.depleted_mode == diag.enriched_mode == al.depleted.log_rho


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
