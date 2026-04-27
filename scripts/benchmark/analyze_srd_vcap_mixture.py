"""SRD v1 — VCaP mixture analysis.

Inputs: rigel quant outputs in
    /scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v1_vcap_mixture/<config>/<library>/

Reports per-library:
- gDNA mass fraction recovered vs nominal dna_reads / (rna+dna)
- mRNA + nRNA total mass stability across mixtures
- Per-transcript mRNA Pearson r vs the dna00m reference run
- SRD calibration summary (gdna_fl_quality, pi_pool, mixture_converged)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

LIB_TO_DNA_M = {
    "mctp_vcap_rna20m_dna00m": 0,
    "mctp_vcap_rna20m_dna01m": 1,
    "mctp_vcap_rna20m_dna02m": 2,
    "mctp_vcap_rna20m_dna05m": 5,
    "mctp_vcap_rna20m_dna10m": 10,
    "mctp_vcap_rna20m_dna20m": 20,
    "mctp_vcap_rna20m_dna40m": 40,
    "mctp_vcap_rna20m_dna80m": 80,
}
RNA_M = 20  # nominal RNA reads in millions


def load_run(rundir: Path) -> dict:
    summary = json.loads((rundir / "summary.json").read_text())
    loci = pd.read_feather(rundir / "loci.feather")
    quant = pd.read_feather(rundir / "quant.feather")
    return {"summary": summary, "loci": loci, "quant": quant}


def per_run_rollup(name: str, run: dict, dna_m: int) -> dict:
    loci = run["loci"]
    quant = run["quant"]
    summary = run["summary"]

    mrna_mass = float(loci["mrna"].sum())
    nrna_mass = float(loci["nrna"].sum())
    gdna_mass = float(loci["gdna"].sum())
    rna_mass = mrna_mass + nrna_mass
    total = mrna_mass + nrna_mass + gdna_mass

    nominal_dna_frac = dna_m / (RNA_M + dna_m) if (RNA_M + dna_m) > 0 else 0.0

    cal = summary.get("calibration", {})
    return {
        "library": name,
        "dna_M": dna_m,
        "nominal_gdna_frac": nominal_dna_frac,
        "n_loci": int(len(loci)),
        "total_fragments": total,
        "mrna_mass": mrna_mass,
        "nrna_mass": nrna_mass,
        "rna_mass": rna_mass,
        "gdna_mass": gdna_mass,
        "rigel_gdna_frac": gdna_mass / total if total > 0 else 0.0,
        "rigel_mrna_frac": mrna_mass / total if total > 0 else 0.0,
        "rigel_nrna_frac": nrna_mass / total if total > 0 else 0.0,
        "pi_pool": cal.get("pi_pool"),
        "n_pool": cal.get("n_pool"),
        "gdna_fl_quality": cal.get("gdna_fl_quality"),
        "mixture_converged": cal.get("mixture_converged"),
        "mixture_iterations": cal.get("mixture_iterations"),
        "strand_specificity": cal.get("strand_specificity"),
        "n_transcripts_quant": int(len(quant)),
        "n_transcripts_with_count": int((quant["count"] > 0).sum()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        default="/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v1_vcap_mixture/default",
    )
    ap.add_argument(
        "--out",
        default="/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v1_vcap_mixture/analysis",
    )
    args = ap.parse_args()
    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    libs = sorted(LIB_TO_DNA_M.keys(), key=lambda k: LIB_TO_DNA_M[k])
    available = []
    runs = {}
    for lib in libs:
        rundir = root / lib
        if (rundir / "summary.json").exists() and (rundir / "loci.feather").exists():
            runs[lib] = load_run(rundir)
            available.append(lib)
        else:
            print(f"[skip] {lib} (no output)")

    if not available:
        print("No completed runs found.")
        return

    rollups = [per_run_rollup(lib, runs[lib], LIB_TO_DNA_M[lib]) for lib in available]
    summary_df = pd.DataFrame(rollups)
    summary_df.to_csv(out / "summary.tsv", sep="\t", index=False)
    print("\n=== Per-run rollup ===")
    print(
        summary_df[
            [
                "library",
                "dna_M",
                "nominal_gdna_frac",
                "rigel_gdna_frac",
                "rigel_mrna_frac",
                "rigel_nrna_frac",
                "pi_pool",
                "gdna_fl_quality",
                "mixture_converged",
                "strand_specificity",
            ]
        ].to_string(index=False)
    )

    # Per-transcript mRNA stability vs dna00m reference
    ref_lib = "mctp_vcap_rna20m_dna00m"
    if ref_lib in runs:
        ref_q = runs[ref_lib]["quant"][["transcript_id", "count"]].rename(
            columns={"count": "count_ref"}
        )
        rows = []
        for lib in available:
            if lib == ref_lib:
                continue
            q = runs[lib]["quant"][["transcript_id", "count"]].rename(
                columns={"count": "count_lib"}
            )
            merged = ref_q.merge(q, on="transcript_id", how="inner")
            mask = (merged["count_ref"] > 0) | (merged["count_lib"] > 0)
            sub = merged[mask]
            r_log = np.corrcoef(
                np.log1p(sub["count_ref"]), np.log1p(sub["count_lib"])
            )[0, 1] if len(sub) > 1 else float("nan")
            r_lin = np.corrcoef(sub["count_ref"], sub["count_lib"])[0, 1] if len(sub) > 1 else float(
                "nan"
            )
            # spearman approximation via rank
            sub_rank_r = (
                sub["count_ref"].rank().corr(sub["count_lib"].rank())
                if len(sub) > 1
                else float("nan")
            )
            rows.append(
                {
                    "library": lib,
                    "dna_M": LIB_TO_DNA_M[lib],
                    "n_transcripts": int(len(sub)),
                    "pearson_log1p": r_log,
                    "pearson_linear": r_lin,
                    "spearman": sub_rank_r,
                    "ref_total": float(sub["count_ref"].sum()),
                    "lib_total": float(sub["count_lib"].sum()),
                }
            )
        stab_df = pd.DataFrame(rows)
        stab_df.to_csv(out / "transcript_stability.tsv", sep="\t", index=False)
        print("\n=== Per-transcript mRNA stability vs dna00m reference ===")
        print(stab_df.to_string(index=False))

    print(f"\nWrote: {out / 'summary.tsv'}")
    if (out / "transcript_stability.tsv").exists():
        print(f"Wrote: {out / 'transcript_stability.tsv'}")


if __name__ == "__main__":
    main()
