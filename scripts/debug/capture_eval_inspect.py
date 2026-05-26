"""Diagnostic inspector for the hyb_capture_500kb synthetic capture suite.

Joins truth abundances with probe coverage and per-condition rigel quant output
to surface where the tool is going wrong with hybrid capture data.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
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


def load_probe_cov() -> pd.DataFrame:
    probes = pd.read_csv(SIM / "reference/capture_probes_on.tsv", sep="\t")
    probes["len"] = probes["end"] - probes["start"]
    return probes.groupby("transcript_id")["len"].sum().rename("probe_cov").reset_index()


def main() -> None:
    truth = pd.read_csv(SIM / "truth_abundances_nrna_none.tsv", sep="\t")
    pcov = load_probe_cov()
    truth = truth.merge(pcov, on="transcript_id", how="left").fillna({"probe_cov": 0})
    truth["captured"] = truth["probe_cov"] > 0
    truth["cov_frac"] = truth["probe_cov"] / truth["spliced_length"]

    # Normalise truth abundance to TPM-like (sum to 1e6 across expressed)
    total_truth = truth["mrna_abundance"].sum()
    truth["true_tpm"] = truth["mrna_abundance"] * 1e6 / total_truth

    rows = []
    for cond in CONDS:
        qp = SIM / cond / "rigel_out" / "quant.feather"
        if not qp.exists():
            continue
        q = pd.read_feather(qp)
        m = truth.merge(
            q[["transcript_id", "count", "count_em", "tpm"]],
            on="transcript_id",
            how="left",
        ).fillna({"count": 0, "count_em": 0, "tpm": 0})
        # Renormalise rigel tpm to expressed-only sum for fair comparison.
        # (not strictly needed; tpm sums to 1e6 already)
        rows.append((cond, m.copy()))

    # 1) Show per-condition by-transcript table for expressed transcripts.
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.float_format", lambda v: f"{v:>10.2f}")
    for cond, m in rows:
        print("\n" + "=" * 100)
        print(f"  {cond}")
        print("=" * 100)
        exp = m[m["mrna_abundance"] > 0].sort_values("mrna_abundance", ascending=False)
        cols = ["transcript_id", "captured", "cov_frac", "true_tpm", "count", "count_em", "tpm"]
        print(exp[cols].to_string(index=False))
        # False positives (unexpressed with count > 1)
        fp = m[(m["mrna_abundance"] == 0) & (m["count"] > 1)]
        if not fp.empty:
            print("\n  FALSE POSITIVES (truth=0, rigel count>1):")
            print(fp[["transcript_id", "captured", "cov_frac", "count", "count_em", "tpm"]].to_string(index=False))

    # 2) Captured vs uncaptured aggregates (only relevant when capture on).
    print("\n" + "=" * 100)
    print("  CAPTURE EXPOSURE BIAS (captured vs uncaptured aggregates)")
    print("=" * 100)
    print(f"  {'condition':<48} {'cap_true_tpm':>12} {'cap_rigel':>10} {'uncap_true':>11} {'uncap_rigel':>12}")
    for cond, m in rows:
        e = m[m["mrna_abundance"] > 0]
        cap = e[e["captured"]]
        uncap = e[~e["captured"]]
        ct = cap["true_tpm"].sum()
        cr = cap["tpm"].sum()
        ut = uncap["true_tpm"].sum()
        ur = uncap["tpm"].sum()
        print(f"  {cond:<48} {ct:>12.0f} {cr:>10.0f} {ut:>11.0f} {ur:>12.0f}")

    # 3) Locus-level summary for each condition.
    print("\n" + "=" * 100)
    print("  LOCUS-LEVEL gDNA / MRNA")
    print("=" * 100)
    for cond in CONDS:
        lp = SIM / cond / "rigel_out" / "loci.feather"
        if not lp.exists():
            continue
        loci = pd.read_feather(lp)
        print(f"\n  --- {cond} ---")
        cols = [
            c
            for c in [
                "locus_id",
                "n_transcripts",
                "locus_span_bp",
                "count_unambig",
                "mrna",
                "nrna",
                "gdna",
                "gdna_rate",
                "gdna_prior",
            ]
            if c in loci.columns
        ]
        print(loci[cols].to_string(index=False))

    # 4) Summary calibration deltas: dump pi_gdna_pool, fl_models, density anchors.
    print("\n" + "=" * 100)
    print("  CALIBRATION DETAILS")
    print("=" * 100)
    for cond in CONDS:
        sp = SIM / cond / "rigel_out" / "summary.json"
        if not sp.exists():
            continue
        s = json.loads(sp.read_text())
        cal = s.get("calibration", {})
        pool = cal.get("pool", {}) or {}
        fl = cal.get("fl_models", {}) or {}
        density = cal.get("density_evidence", {}) or {}
        priors = density.get("priors", {}) or {}
        intron = priors.get("INTRON") or {}
        all_p = priors.get("ALL") or {}
        intergenic = priors.get("INTERGENIC") or {}
        q = s.get("quantification", {})
        print(f"\n  --- {cond} ---")
        print(f"    pool pi_gdna       : {pool.get('pi_gdna_pool')}")
        print(f"    fl rna mean        : {fl.get('rna_fl_mean'):.2f} (truth 250)")
        gdna_fl = fl.get("gdna_fl_mean")
        print(f"    fl gdna mean       : {gdna_fl:.2f} (truth 150)  quality={fl.get('gdna_quality')}")
        print(f"    intergenic prior   : {intergenic if intergenic else 'None'}")
        print(f"    intron mean_density: {intron.get('mean_density'):.6f} n_frags={intron.get('n_fragments')} n_reg={intron.get('n_regions')}" if intron else "    intron prior       : None")
        print(f"    ALL    mean_density: {all_p.get('mean_density'):.6f} n_frags={all_p.get('n_fragments')} n_reg={all_p.get('n_regions')}" if all_p else "    ALL prior          : None")
        print(f"    gdna_fraction est  : {q.get('gdna_fraction'):.4f}  intergenic_frac={q.get('intergenic_fraction'):.4f}")


if __name__ == "__main__":
    main()
