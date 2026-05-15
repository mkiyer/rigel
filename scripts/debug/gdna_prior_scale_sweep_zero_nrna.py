#!/usr/bin/env python3
"""Counterfactual gDNA-prior scale sweep for the unstranded high-nRNA case."""

from __future__ import annotations

import argparse
import copy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from inspect_zero_nrna_candidate_sets import (
    DEFAULT_BASE,
    DEFAULT_CONDITION,
    rebuild_scored_fragments,
)
from rigel.config import EMConfig
from rigel.locus_partition import partition_and_free
from rigel.pipeline import _assign_locus_ids, _run_locus_em_partitioned
from rigel.scored_fragments import ScoredFragments


DEFAULT_OUT = Path("results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high")


def copy_scored_fragments(em_data: ScoredFragments) -> ScoredFragments:
    return replace(
        em_data,
        offsets=em_data.offsets.copy(),
        t_indices=em_data.t_indices.copy(),
        log_liks=em_data.log_liks.copy(),
        count_cols=em_data.count_cols.copy(),
        coverage_weights=em_data.coverage_weights.copy(),
        locus_t_indices=em_data.locus_t_indices.copy(),
        locus_count_cols=em_data.locus_count_cols.copy(),
        is_spliced=em_data.is_spliced.copy(),
        gdna_log_liks=em_data.gdna_log_liks.copy(),
        frag_ids=em_data.frag_ids.copy(),
        frag_class=em_data.frag_class.copy(),
        splice_type=em_data.splice_type.copy(),
    )


def run_scale(base_estimator, index, em_data, multi_loci, prior_table, scale: float) -> dict:
    em_config = EMConfig(mode="vbem", assignment_mode="fractional", n_threads=0)
    estimator = copy.deepcopy(base_estimator)
    estimator.em_config = em_config
    estimator.locus_results.clear()
    estimator._gdna_em_total = 0.0
    _assign_locus_ids(estimator, multi_loci)
    partitions = partition_and_free(copy_scored_fragments(em_data), multi_loci)
    gdna_prior = np.ascontiguousarray(prior_table.gdna_prior_count * float(scale), dtype=np.float64)
    _run_locus_em_partitioned(
        estimator,
        partitions,
        multi_loci,
        index,
        gdna_prior,
        em_config,
        enable_gdna=prior_table.enable_gdna,
        emit_locus_stats=False,
        annotations=None,
    )
    loci = estimator.get_loci_df(index)
    nrna = estimator.get_nrna_counts_df(index)
    total = float(loci["total"].sum())
    return {
        "scale": float(scale),
        "mrna": float(loci["mrna"].sum()),
        "nrna": float(loci["nrna"].sum()),
        "gdna": float(loci["gdna"].sum()),
        "total": total,
        "gdna_rate": float(loci["gdna"].sum() / total) if total > 0 else 0.0,
        "max_locus_nrna": float(loci["nrna"].max()) if len(loci) else 0.0,
        "max_entity_nrna": float(nrna["count"].max()) if len(nrna) else 0.0,
        "locus15_nrna": float(loci.loc[loci["locus_id"] == 15, "nrna"].sum()),
        "locus15_gdna": float(loci.loc[loci["locus_id"] == 15, "gdna"].sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--scales",
        type=float,
        nargs="+",
        default=[0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
    )
    args = parser.parse_args()

    index, base_estimator, em_data, multi_loci, prior_table, calibration = rebuild_scored_fragments(
        args.base.resolve(),
        args.condition,
    )
    del calibration
    rows = [run_scale(base_estimator, index, em_data, multi_loci, prior_table, s) for s in args.scales]
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "gdna_prior_scale_sweep.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()