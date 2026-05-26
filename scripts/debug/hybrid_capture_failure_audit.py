"""Summarize the hyb_capture_500kb failure mode from existing rigel outputs.

This is a read-only diagnostic script. It expects the synthetic suite produced by
``rigel.sim.analysis`` and prints compact tables used by the hybrid-capture
calibration diagnosis memo.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SIM = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb")

CONDS = [
    "gdna_none_ss_0.99_nrna_none_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_high_ss_0.99_nrna_none_capture_off",
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_off",
    "gdna_high_ss_0.50_nrna_none_capture_on",
]


def _condition_parts(cond: str) -> tuple[str, str, str]:
    gdna = "high" if "gdna_high" in cond else "none"
    ss = "0.99" if "ss_0.99" in cond else "0.50"
    capture = "on" if cond.endswith("capture_on") else "off"
    return gdna, ss, capture


def _summary(cond: str) -> dict:
    with (SIM / cond / "rigel_out" / "summary.json").open() as handle:
        return json.load(handle)


def _prior(priors: dict, name: str) -> dict:
    value = priors.get(name)
    return value if isinstance(value, dict) else {}


def summary_table() -> pd.DataFrame:
    rows = []
    for cond in CONDS:
        s = _summary(cond)
        gdna, ss, capture = _condition_parts(cond)
        cal = s.get("calibration", {})
        density = cal.get("density_evidence", {})
        priors = density.get("priors", {}) or {}
        diag = cal.get("diagnostics", {})
        fl = cal.get("fl_models", {})
        q = s.get("quantification", {})
        region_exposure = cal.get("region_exposure", {})
        prior_table = cal.get("prior_table", {})
        intron = _prior(priors, "INTRON")
        all_prior = _prior(priors, "ALL")
        rows.append(
            {
                "condition": cond,
                "gdna": gdna,
                "ss": ss,
                "capture": capture,
                "rho_ref": density.get("rho_ref", 0.0),
                "rho_intron": intron.get("mean_density", 0.0),
                "rho_all": all_prior.get("mean_density", 0.0),
                "intron_contained": diag.get("fl_pool_total", {}).get("INTRONIC_CONTAINED", 0.0),
                "intron_boundary": diag.get("fl_pool_total", {}).get("INTRONIC_BOUNDARY", 0.0),
                "exon_contained": diag.get("fl_pool_total", {}).get("EXONIC_CONTAINED", 0.0),
                "exon_boundary": diag.get("fl_pool_total", {}).get("EXONIC_BOUNDARY", 0.0),
                "intergenic_mass": diag.get("mass_by_coarse_class", {}).get("INTERGENIC", 0.0),
                "intron_mass": diag.get("mass_by_coarse_class", {}).get("INTRON", 0.0),
                "exon_mass": diag.get("mass_by_coarse_class", {}).get("EXON", 0.0),
                "gdna_fl_mean": fl.get("gdna_fl_mean", 0.0),
                "rna_fl_mean": fl.get("rna_fl_mean", 0.0),
                "gdna_fraction": q.get("gdna_fraction", 0.0),
                "A_mean": region_exposure.get("A_mean", 0.0),
                "A_max": region_exposure.get("A_max", 0.0),
                "prior_sum": prior_table.get("sum_gdna_prior_count_em", 0.0),
            }
        )
    return pd.DataFrame(rows)


def captured_truth_table() -> pd.DataFrame:
    truth = pd.read_csv(SIM / "truth_abundances_nrna_none.tsv", sep="\t")
    probes = pd.read_csv(SIM / "reference" / "capture_probes_on.tsv", sep="\t")
    probes["probe_len"] = probes["end"] - probes["start"]
    probe_cov = probes.groupby("transcript_id", as_index=False)["probe_len"].sum()
    truth = truth.merge(probe_cov, on="transcript_id", how="left").fillna({"probe_len": 0})
    truth["captured"] = truth["probe_len"] > 0
    truth["true_tpm"] = truth["mrna_abundance"] * 1_000_000.0 / truth["mrna_abundance"].sum()

    rows = []
    for cond in CONDS:
        q = pd.read_feather(SIM / cond / "rigel_out" / "quant.feather")
        merged = truth.merge(q[["transcript_id", "tpm"]], on="transcript_id", how="left")
        expressed = merged[merged["mrna_abundance"] > 0]
        for captured in (True, False):
            subset = expressed[expressed["captured"] == captured]
            rows.append(
                {
                    "condition": cond,
                    "captured": captured,
                    "true_tpm": subset["true_tpm"].sum(),
                    "rigel_tpm": subset["tpm"].sum(),
                }
            )
    return pd.DataFrame(rows)


def locus_comparison(ss: str = "0.99") -> pd.DataFrame:
    off = pd.read_feather(
        SIM / f"gdna_high_ss_{ss}_nrna_none_capture_off" / "rigel_out" / "loci.feather"
    )
    on = pd.read_feather(
        SIM / f"gdna_high_ss_{ss}_nrna_none_capture_on" / "rigel_out" / "loci.feather"
    )
    cols = [
        "locus_id",
        "n_em_fragments",
        "mrna",
        "gdna",
        "gdna_rate",
        "gdna_prior",
        "gdna_prior_count_em",
        "gdna_eff_len_weight_ratio",
    ]
    joined = off[cols].merge(on[cols], on="locus_id", suffixes=("_off", "_on"))
    joined["prior_fold_on_vs_off"] = joined["gdna_prior_on"] / joined["gdna_prior_off"]
    joined["gdna_fold_on_vs_off"] = (joined["gdna_on"] + 1.0) / (joined["gdna_off"] + 1.0)
    return joined.sort_values("gdna_on", ascending=False)


def main() -> None:
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.float_format", lambda value: f"{value:0.4g}")

    print("\nSUMMARY BY CONDITION")
    keep = [
        "condition",
        "rho_ref",
        "rho_intron",
        "intron_contained",
        "intron_boundary",
        "exon_contained",
        "exon_boundary",
        "gdna_fl_mean",
        "gdna_fraction",
        "A_mean",
        "A_max",
        "prior_sum",
    ]
    print(summary_table()[keep].to_string(index=False))

    print("\nCAPTURED VS UNCAPTURED EXPRESSED TPM")
    cap = captured_truth_table()
    pivot = cap.pivot(index="condition", columns="captured", values=["true_tpm", "rigel_tpm"])
    print(pivot.to_string())

    print("\nLOCUS COMPARISON, GDNA HIGH SS=0.99 CAPTURE OFF VS ON")
    print(locus_comparison("0.99").to_string(index=False))


if __name__ == "__main__":
    main()
