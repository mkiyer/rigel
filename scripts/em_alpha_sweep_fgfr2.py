#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd

from hulkrna.index import HulkIndex
from hulkrna.pipeline import run_pipeline


def score(df: pd.DataFrame, truth: dict[str, int]) -> tuple[float, float, float, float]:
    observed = dict(zip(df["transcript_id"], df["count"]))
    tids = sorted(truth)
    t = np.array([truth[x] for x in tids], dtype=float)
    o = np.array([observed.get(x, 0.0) for x in tids], dtype=float)
    mae = float(np.mean(np.abs(o - t)))
    rmse = float(np.sqrt(np.mean((o - t) ** 2)))
    pear = float(np.corrcoef(t, o)[0, 1]) if np.std(t) > 0 and np.std(o) > 0 else 0.0
    return mae, rmse, pear, float(o.sum())


def main() -> int:
    base = Path(
        "/Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs/seed_303/chr10_121476640_121601584"
    )
    index = HulkIndex.load(base / "hulkrna_index")
    bam = base / "align" / "reads.bam"
    truth_df = pd.read_csv(base / "per_transcript_counts.csv")
    truth = dict(zip(truth_df["transcript_id"], truth_df["truth"]))

    print("alpha,mae,rmse,pearson,obs_total")
    for alpha in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        pr = run_pipeline(
            bam,
            index,
            seed=42,
            sj_strand_tag="auto",
            gdna_threshold=0.0,
            em_pseudocount=alpha,
            em_iterations=500,
        )
        df = pr.counter.get_counts_df(index)
        mae, rmse, pear, tot = score(df, truth)
        print(f"{alpha},{mae:.3f},{rmse:.3f},{pear:.5f},{tot:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
