"""Phase-0 report-substrate contract for ``rigel quant`` outputs.

Guards the v2 ``summary.json`` schema and the externalized
``fragment_lengths.feather`` companion:

* the fragment-length section is stats-only (no raw per-bin histograms bloat);
* the raw histograms live in a tidy ``(category, length, count)`` feather that
  reconciles with the summary observation counts;
* the splice-type breakdown (incl. implicit / artifact) and the strand
  contamination diagnostic are surfaced.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rigel import cli
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.frag_length_model import FragmentLengthModels
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType

SEED = 42


# ---------------------------------------------------------------------------
# Unit: the fragment-length split helper
# ---------------------------------------------------------------------------


def test_fragment_length_report_splits_lean_summary_from_histograms():
    flm = FragmentLengthModels(max_size=1000)
    lengths = np.array([100, 100, 150, 200, 200, 200, 300], dtype=np.intp)
    stypes = np.array(
        [
            SpliceType.UNSPLICED,
            SpliceType.UNSPLICED,
            SpliceType.SPLICED_ANNOT,
            SpliceType.SPLICED_ANNOT,
            SpliceType.SPLICED_ANNOT,
            SpliceType.SPLICED_IMPLICIT,
            SpliceType.SPLICE_ARTIFACT,
        ],
        dtype=np.intp,
    )
    flm.observe_batch(lengths, stypes)

    summary, hist_df = cli._fragment_length_report(flm, None)

    # summary is lean: stats only, never raw per-bin arrays
    assert "global" in summary
    for cat, s in summary.items():
        assert "histogram" not in s, f"{cat} still carries raw histogram"
        assert set(s) == {
            "n_observations",
            "mean",
            "std",
            "median",
            "mode",
            "max_size",
            "overflow_count",
            "overflow_fraction",
        }
    assert summary["global"]["n_observations"] == 7

    # tidy histogram table reconciles with the summary observation counts
    assert list(hist_df.columns) == ["category", "length", "count"]
    assert str(hist_df["length"].dtype) == "int32"
    for cat in hist_df["category"].unique():
        fsum = float(hist_df.loc[hist_df["category"] == cat, "count"].sum())
        assert abs(fsum - summary[cat]["n_observations"]) < 0.5


# ---------------------------------------------------------------------------
# Integration: the full quant summary.json + fragment_lengths.feather
# ---------------------------------------------------------------------------


def _build_scenario(work_dir: Path):
    sc = Scenario("summary_report", genome_length=12000, seed=SEED, work_dir=work_dir)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(500, 900), (1500, 1900), (2500, 2900)], "abundance": 80},
            {"t_id": "t2", "exons": [(500, 900), (2500, 2900)], "abundance": 40},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(6000, 6400), (7000, 7400)], "abundance": 50},
        ],
    )
    sim = ReadSimConfig(
        frag_mean=200,
        frag_std=40,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.9,
        seed=SEED,
    )
    return sc, sc.build_oracle(
        n_fragments=2000, sim_config=sim, gdna_fraction=0.2, nrna_abundance=15.0
    )


def test_summary_json_v2_schema_and_companion(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    sc, result = _build_scenario(tmp_path / "work")
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)

    args = argparse.Namespace(
        bam_file=str(result.bam_path),
        index_dir=str(result.index_dir),
        config=None,
        tsv=False,
        emit_locus_stats=False,
    )
    cli._write_quant_outputs(pr, result.index, out, args)

    summary = json.loads((out / "summary.json").read_text())

    # versioned schema
    assert summary["schema_version"] == cli.SUMMARY_SCHEMA_VERSION == 2

    # fragment_length is stats-only — the histogram bloat is gone
    fl = summary["fragment_length"]
    assert "global" in fl
    for d in fl.values():
        assert "histogram" not in d

    # splice breakdown surfaces implicit / artifact / blacklist
    splice = summary["fragment_stats"]["splice"]
    assert set(splice) == {
        "unspliced",
        "spliced_annotated",
        "spliced_unannotated",
        "spliced_implicit",
        "splice_artifact",
        "sj_blacklisted",
    }

    # strand contamination diagnostic present
    diag = summary["strand_model"]["diagnostics"]
    assert {
        "exonic_all_specificity",
        "exonic_all_p_r1_sense",
        "exonic_all_n_fragments",
        "contamination_gap",
    } <= set(diag)

    # companion histogram table exists and reconciles with the summary
    flf = pd.read_feather(out / "fragment_lengths.feather")
    assert list(flf.columns) == ["category", "length", "count"]
    for cat in flf["category"].unique():
        fsum = float(flf.loc[flf["category"] == cat, "count"].sum())
        assert abs(fsum - fl[cat]["n_observations"]) < 0.5

    # genome gDNA track: feather substrate + genome-browser bedGraph
    assert pr.calibration_track is not None
    track = pd.read_feather(out / "calibration_track.feather")
    assert list(track.columns) == [
        "ref",
        "start",
        "end",
        "gdna_mass",
        "rna_mass",
        "gdna_density",
        "gdna_frac",
    ]
    assert len(track) == pr.calibration.n_regions
    assert (track["gdna_frac"] >= 0).all() and (track["gdna_frac"] <= 1).all()
    bg = (out / "calibration_track.bedgraph").read_text().splitlines()
    assert bg[0].startswith("track type=bedGraph")
    assert len(bg) == len(track) + 1  # header + one line per region

    # Mass-weighted capture-enrichment summary — present when the track has
    # enough informative nodes (capture_summary returns non-None).
    if "capture" in summary["calibration"]:
        cap = summary["calibration"]["capture"]
        assert {
            "n_nodes",
            "enriched",
            "count_median_log_rho",
            "background_mode_log_rho",
            "enriched_mode_log_rho",
            "fold_peak_to_peak",
            "fold_vs_median",
            "mass_frac_ontarget",
            "kde_bandwidth_factor",
        } <= set(cap)

    # The prior's own equal-weight KDE is still persisted for provenance.
    if pr.calibration_diagnostics is not None:
        kde = pd.read_feather(out / "gdna_density_kde.feather")
        assert {"log_rho", "log_density", "density"} <= set(kde.columns)
        assert (out / "gdna_density_nodes.feather").exists()

    sc.cleanup()
