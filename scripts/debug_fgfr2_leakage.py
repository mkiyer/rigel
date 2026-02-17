#!/usr/bin/env python3
from pathlib import Path

import pandas as pd

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline
from benchmark_region_competition import parse_truth_from_fastq


def main() -> int:
    base = Path(
        "/Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs/seed_303/chr10_121476640_121601584"
    )
    index = HulkIndex.load(base / "hulkrna_index")
    bam = base / "align" / "reads.bam"
    truth_counts, _ = parse_truth_from_fastq(base / "reads" / "sim_R1.fq")

    pr = run_pipeline(
        bam,
        index,
        seed=42,
        sj_strand_tag="auto",
        gdna_threshold=0.0,
        em_pseudocount=0.01,
        em_iterations=1000,
    )
    df = pr.counter.get_counts_df(index)

    df["truth"] = df["transcript_id"].map(truth_counts).fillna(0.0)
    df["abs_err"] = (df["count"] - df["truth"]).abs()

    total_pred = float(df["count"].sum())
    total_truth = float(df["truth"].sum())
    zero_truth_mass = float(df.loc[df["truth"] == 0, "count"].sum())
    nonzero_truth_mass = float(df.loc[df["truth"] > 0, "count"].sum())

    print(f"total_truth={total_truth:.1f}")
    print(f"total_pred={total_pred:.1f}")
    print(f"pred_mass_on_zero_truth={zero_truth_mass:.1f} ({100.0*zero_truth_mass/max(total_pred,1.0):.2f}%)")
    print(f"pred_mass_on_nonzero_truth={nonzero_truth_mass:.1f}")

    print("\nTop 15 zero-truth transcripts by predicted count:")
    print(
        df.loc[df["truth"] == 0, ["transcript_id", "count"]]
        .sort_values("count", ascending=False)
        .head(15)
        .to_string(index=False)
    )

    print("\nTop 15 absolute errors (all transcripts):")
    print(
        df[["transcript_id", "truth", "count", "abs_err"]]
        .sort_values("abs_err", ascending=False)
        .head(15)
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
