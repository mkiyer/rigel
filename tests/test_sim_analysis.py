"""Regression tests for synthetic simulation analysis reporting."""

import importlib.util
import json
from pathlib import Path

import pandas as pd


_SCRIPT_PATH = Path(__file__).parents[1] / "scripts/sim/run_rigel_analysis.py"
_SPEC = importlib.util.spec_from_file_location("run_rigel_analysis", _SCRIPT_PATH)
assert _SPEC is not None and _SPEC.loader is not None
run_rigel_analysis = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(run_rigel_analysis)

analyze_calibration = run_rigel_analysis.analyze_calibration
analyze_postfix_acceptance = run_rigel_analysis.analyze_postfix_acceptance


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
                "n_rna": n_rna,
                "n_gdna": n_gdna,
            }
        ],
    }
    (sim_base / "manifest.json").write_text(json.dumps(manifest))


def _write_summary(sim_base, condition: str, *, rho_ex: float = 0.0985):
    out_dir = sim_base / condition / "rigel_out"
    out_dir.mkdir(parents=True)
    summary = {
        "calibration": {
            "global_densities": {
                "INTERGENIC": {"rho": 0.1000},
                "INTRON": {"rho": 0.1005},
                "EXON-INTRON": {"rho": rho_ex},
            },
            "fl_models": {"rna_fl_mean": 257.0, "gdna_fl_mean": 351.7},
            "n_multi_loci": 51,
        },
        "quantification": {"gdna_fraction": 0.0913},
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


def test_postfix_acceptance_checks_pass_for_post_fix_profile(tmp_path):
    condition = "gdna_high_ss_0.99_nrna_none"
    _write_manifest(tmp_path, condition, n_gdna=2_000_000)
    _write_summary(tmp_path, condition, rho_ex=0.0985)
    loci_path = tmp_path / condition / "rigel_out" / "loci.tsv"
    loci_path.write_text("nrna\n0\n")

    assignment_rows = [
        {"condition": condition, "total_gdna": 2_000_000, "gdna_as_rna": 9_430}
    ]

    report = analyze_postfix_acceptance(tmp_path, [condition], assignment_rows)

    assert "rho_ex/rho_ig" in report
    assert "nRNA in nrna_none" in report
    assert "gDNA->RNA leak" in report
    assert "FAIL" not in report
    assert "PASS: all 3 evaluated acceptance checks passed." in report