"""Regression tests for synthetic simulation analysis reporting."""

import json

import pandas as pd

from rigel.sim.analysis import (
    analyze_calibration,
    analyze_postfix_acceptance,
    discover_conditions,
    load_condition_truth,
    load_truth,
)


def _write_manifest(sim_base, condition: str, *, n_rna: int = 1_000_000, n_gdna: int = 100_000):
    manifest = {
        "simulation": {"frag_mean": 250, "frag_std": 50},
        "gdna": {"frag_mean": 350, "frag_std": 100},
        "conditions": [
            {
                "name": condition,
                "gdna_label": "low",
                "gdna_rate": n_gdna / n_rna,
                "strand_specificity": 0.99,
                "nrna_label": "none",
                "n_mrna": n_rna,
                "n_nrna": 0,
                "n_rna": n_rna,
                "n_gdna": n_gdna,
            }
        ],
    }
    (sim_base / "manifest.json").write_text(json.dumps(manifest))


def _write_summary(sim_base, condition: str):
    out_dir = sim_base / condition / "rigel_out"
    out_dir.mkdir(parents=True)
    summary = {
        # Acyclic-calibrator schema: one global gDNA density + strand-model scalars.
        "calibration": {
            "gdna_density_global": 0.1003,
            "rna_sense_frac": 0.99,
            "gdna_strand_overdispersion": 0.1,
            "rna_strand_overdispersion": 0.1,
            "n_regions": 51,
        },
        "fragment_length": {
            "rna": {"summary": {"mean": 257.0}},
            "gdna": {"summary": {"mean": 351.7}},
        },
        "quantification": {"gdna_fraction": 0.0913, "n_loci": 51},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary))
    return out_dir


def test_calibration_report_uses_manifest_fragment_lengths(tmp_path):
    condition = "gdna_low_ss_0.99_nrna_none"
    _write_manifest(tmp_path, condition)
    _write_summary(tmp_path, condition)

    report = analyze_calibration(tmp_path, [condition], pd.DataFrame())

    assert "gDNA=350bp" in report
    assert "gDNA=200bp" not in report
    assert "+0.005" in report
    assert "+0.758" not in report


def test_calibration_report_uses_condition_fl_truth_when_available(tmp_path):
    condition = "gdna_low_ss_0.99_nrna_none"
    _write_manifest(tmp_path, condition)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["conditions"][0]["truth_summary"] = f"{condition}/truth_summary.json"
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    (tmp_path / condition).mkdir()
    (tmp_path / condition / "truth_summary.json").write_text(json.dumps({
        "fragment_lengths": {
            "mrna": {"n": 10, "mean": 260.0},
            "gdna": {"n": 10, "mean": 300.0},
        }
    }))
    _write_summary(tmp_path, condition)

    report = analyze_calibration(tmp_path, [condition], pd.DataFrame())

    assert "+0.172" in report
    assert "+0.005" not in report


def test_postfix_acceptance_checks_pass_for_post_fix_profile(tmp_path):
    condition = "gdna_high_ss_0.99_nrna_none"
    _write_manifest(tmp_path, condition, n_gdna=2_000_000)
    _write_summary(tmp_path, condition)
    loci_path = tmp_path / condition / "rigel_out" / "loci.tsv"
    loci_path.write_text("nrna\n0\n")

    assignment_rows = [
        {"condition": condition, "total_gdna": 2_000_000, "gdna_as_rna": 9_430}
    ]

    report = analyze_postfix_acceptance(tmp_path, [condition], assignment_rows)

    assert "nRNA in nrna_none" in report
    assert "gDNA->RNA leak" in report
    assert "FAIL" not in report
    assert "PASS: all 2 evaluated acceptance checks passed." in report


def test_analysis_discovers_conditions_from_manifest(tmp_path):
    _write_manifest(tmp_path, "gdna_none_ss_1.00_nrna_none")
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["conditions"].append({
        "name": "gdna_none_ss_1.00_nrna_high",
        "gdna_label": "none",
        "gdna_rate": 0.0,
        "strand_specificity": 1.0,
        "nrna_label": "high",
        "n_rna": 150,
        "n_gdna": 0,
    })
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))

    assert discover_conditions(tmp_path) == [
        "gdna_none_ss_1.00_nrna_none",
        "gdna_none_ss_1.00_nrna_high",
    ]


def test_analysis_loads_condition_specific_truth(tmp_path):
    condition = "gdna_none_ss_1.00_nrna_high"
    _write_manifest(tmp_path, condition)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["truth_abundances"] = "truth_abundances_nrna_none.tsv"
    manifest["conditions"][0]["truth_abundances"] = "truth_abundances_nrna_high.tsv"
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))

    base_truth = "transcript_id\tmrna_abundance\nTX1\t10\n"
    high_truth = "transcript_id\tmrna_abundance\nTX1\t20\n"
    (tmp_path / "truth_abundances_nrna_none.tsv").write_text(base_truth)
    (tmp_path / "truth_abundances_nrna_high.tsv").write_text(high_truth)

    fallback = load_truth(tmp_path)
    loaded = load_condition_truth(
        tmp_path,
        condition,
        {condition: manifest["conditions"][0]},
        fallback,
    )

    assert fallback.loc[0, "mrna_abundance"] == 10
    assert loaded.loc[0, "mrna_abundance"] == 20