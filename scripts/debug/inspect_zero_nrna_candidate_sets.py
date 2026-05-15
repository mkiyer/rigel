#!/usr/bin/env python3
"""Rebuild candidate sets for the zero-gDNA/unstranded/high-nRNA case.

This is a focused follow-up to ``interrogate_zero_nrna_high.py``. It rebuilds
Rigel's scored EM units for the target condition, joins them to the annotated
BAM fragment IDs, and verifies whether true nRNA fragments actually have their
synthetic nRNA parent in the candidate set.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from rigel.calibration import calibrate
from rigel.calibration.locus_prior import assemble_priors
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.index import TranscriptIndex
from rigel.locus import build_multi_loci
from rigel.pipeline import (
    _calibration_strand_summary,
    _native_detect_sj_tag,
    _score_fragments,
    _setup_geometry_and_estimator,
    scan_and_buffer,
)


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_CONDITION = "gdna_zero_ss_0.50_nrna_high"
DEFAULT_DIAG_DIR = Path("results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high")
DEFAULT_FRACTIONAL_DIR = DEFAULT_DIAG_DIR / "rigel_fractional_debug"


def safe_div(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else math.nan


def rebuild_scored_fragments(base: Path, condition: str):
    bam_path = base / condition / "sim_oracle.bam"
    index = TranscriptIndex.load(base / "rigel_index")
    sj_spec = _native_detect_sj_tag(str(bam_path))
    scan_cfg = BamScanConfig(sj_strand_tag=sj_spec)

    stats, strand_models, frag_length_models, buffer, calibration_payload = scan_and_buffer(
        str(bam_path),
        index,
        scan_cfg,
    )
    strand_models.finalize()
    calibration = calibrate(
        index=index,
        payload=calibration_payload,
        scan_trained=frag_length_models,
        fl_prior_ess=1000.0,
        pool_quality_good=5000,
        pool_quality_weak=200,
        strand_summary=_calibration_strand_summary(strand_models),
    )
    geometry, estimator = _setup_geometry_and_estimator(index, calibration.fl_models.rna, EMConfig())
    del geometry
    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        calibration.fl_models.rna,
        calibration.fl_models.gdna,
        stats,
        estimator,
        FragmentScoringConfig(),
        log_every=10_000_000,
        annotations=None,
    )
    multi_loci = build_multi_loci(em_data, index)
    prior_table = assemble_priors(
        multi_loci,
        em_data,
        index,
        calibration_payload,
        calibration.global_densities,
        gdna_fl=calibration.fl_models.gdna,
        splicing_anchor_tolerance=int(getattr(calibration_payload, "splicing_anchor_tolerance", 0)),
    )
    return index, estimator, em_data, multi_loci, prior_table, calibration


def build_unit_lookup(em_data) -> dict[int, int]:
    return {int(frag_id): unit_index for unit_index, frag_id in enumerate(em_data.frag_ids)}


def candidate_summary_for_fragment(
    row: pd.Series,
    em_data,
    unit_lookup: dict[int, int],
    t_df: pd.DataFrame,
) -> dict[str, Any]:
    frag_id = int(row.frag_id)
    unit_index = unit_lookup.get(frag_id, -1)
    if unit_index < 0:
        return {
            "frag_id": frag_id,
            "qname": row.qname,
            "origin_tx": row.origin_tx,
            "origin_nrna_id": row.origin_nrna_id,
            "unit_index": -1,
            "origin_nrna_candidate_present": False,
            "n_rna_candidates": 0,
            "n_synthetic_candidates": 0,
            "n_mrna_candidates": 0,
            "has_gdna_candidate": False,
        }

    start = int(em_data.offsets[unit_index])
    end = int(em_data.offsets[unit_index + 1])
    candidate_indices = em_data.t_indices[start:end].astype(int)
    candidates = t_df.iloc[candidate_indices]
    candidate_ids = candidates["t_id"].astype(str).tolist()
    is_synthetic = candidates["is_synthetic"].astype(bool).to_numpy()
    is_nrna = candidates["is_nrna"].astype(bool).to_numpy()
    has_gdna = (not bool(em_data.is_spliced[unit_index])) and np.isfinite(
        float(em_data.gdna_log_liks[unit_index])
    )
    return {
        "frag_id": frag_id,
        "qname": row.qname,
        "origin_tx": row.origin_tx,
        "origin_nrna_id": row.origin_nrna_id,
        "assigned_pool": row.assigned_pool,
        "region_class": row.region_class,
        "zc": row.zc,
        "winner_posterior": float(row.posterior),
        "unit_index": unit_index,
        "origin_nrna_candidate_present": str(row.origin_nrna_id) in set(candidate_ids),
        "n_rna_candidates": len(candidate_ids),
        "n_synthetic_candidates": int(is_synthetic.sum()),
        "n_nrna_candidates": int(is_nrna.sum()),
        "n_mrna_candidates": int((~is_synthetic).sum()),
        "has_gdna_candidate": bool(has_gdna),
        "gdna_log_lik": float(em_data.gdna_log_liks[unit_index]),
        "candidate_ids": ",".join(candidate_ids),
        "synthetic_candidate_ids": ",".join(candidates.loc[is_synthetic, "t_id"].astype(str)),
    }


def example_candidate_rows(
    examples: pd.DataFrame,
    em_data,
    unit_lookup: dict[int, int],
    t_df: pd.DataFrame,
    effective_lengths: np.ndarray,
) -> pd.DataFrame:
    rows = []
    for example in examples.itertuples(index=False):
        unit_index = unit_lookup.get(int(example.frag_id), -1)
        if unit_index < 0:
            continue
        start = int(em_data.offsets[unit_index])
        end = int(em_data.offsets[unit_index + 1])
        for flat_index in range(start, end):
            t_index = int(em_data.t_indices[flat_index])
            t_row = t_df.iloc[t_index]
            rows.append(
                {
                    "frag_id": int(example.frag_id),
                    "qname": example.qname,
                    "origin_tx": example.origin_tx,
                    "assigned_pool": example.assigned_pool,
                    "winner_posterior": float(example.posterior),
                    "candidate_t_index": t_index,
                    "candidate_id": str(t_row.t_id),
                    "candidate_gene": str(t_row.g_id),
                    "is_nrna": bool(t_row.is_nrna),
                    "is_synthetic": bool(t_row.is_synthetic),
                    "is_origin_nrna_parent": str(t_row.t_id) == str(example.origin_nrna_id),
                    "log_lik": float(em_data.log_liks[flat_index]),
                    "coverage_weight": float(em_data.coverage_weights[flat_index]),
                    "effective_length": float(effective_lengths[t_index]),
                }
            )
        rows.append(
            {
                "frag_id": int(example.frag_id),
                "qname": example.qname,
                "origin_tx": example.origin_tx,
                "assigned_pool": example.assigned_pool,
                "winner_posterior": float(example.posterior),
                "candidate_t_index": -2,
                "candidate_id": "GDNA",
                "candidate_gene": ".",
                "is_nrna": False,
                "is_synthetic": False,
                "is_origin_nrna_parent": False,
                "log_lik": float(em_data.gdna_log_liks[unit_index]),
                "coverage_weight": 1.0,
                "effective_length": 1.0,
            }
        )
    return pd.DataFrame(rows)


def load_fractional_alpha(
    fractional_dir: Path,
    t_df: pd.DataFrame,
    loci_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    quant = pd.read_csv(fractional_dir / "quant.tsv", sep="\t")
    nrna_quant = pd.read_csv(fractional_dir / "nrna_quant.tsv", sep="\t")

    transcript_counts = dict(zip(quant["transcript_id"].astype(str), quant["count"], strict=False))
    nrna_counts = dict(zip(nrna_quant["nrna_id"].astype(str), nrna_quant["count"], strict=False))

    alpha = np.empty(len(t_df), dtype=np.float64)
    for t_index, row in enumerate(t_df.itertuples(index=False)):
        t_id = str(row.t_id)
        if bool(row.is_synthetic):
            count = float(nrna_counts.get(t_id, 0.0))
        else:
            count = float(transcript_counts.get(t_id, 0.0))
        alpha[t_index] = count + 0.5

    max_locus = int(loci_df["locus_id"].max()) if len(loci_df) else -1
    gdna_alpha = np.zeros(max_locus + 1, dtype=np.float64)
    for row in loci_df.itertuples(index=False):
        gdna_alpha[int(row.locus_id)] = float(row.gdna) + float(row.gdna_prior_count) + 0.5
    return alpha, gdna_alpha


def reconstruct_posterior_distribution(
    fragment_rows: pd.DataFrame,
    em_data,
    unit_lookup: dict[int, int],
    t_df: pd.DataFrame,
    effective_lengths: np.ndarray,
    multi_loci,
    fractional_dir: Path,
    diag_dir: Path,
) -> None:
    loci_df = pd.read_csv(fractional_dir / "loci.tsv", sep="\t")
    alpha, gdna_alpha = load_fractional_alpha(fractional_dir, t_df, loci_df)

    unit_locus = np.full(em_data.n_units, -1, dtype=np.int32)
    for locus in multi_loci:
        unit_locus[locus.unit_indices.astype(np.int64)] = int(locus.multi_locus_id)

    t_ids = t_df["t_id"].astype(str).to_numpy()
    is_synthetic = t_df["is_synthetic"].astype(bool).to_numpy()

    rows = []
    true_nrna = fragment_rows[fragment_rows["origin"] == "nrna"]
    true_nrna = true_nrna[true_nrna["assigned_pool"].isin(["gdna", "mrna", "nrna"])]
    for row in true_nrna.itertuples(index=False):
        unit_index = unit_lookup.get(int(row.frag_id), -1)
        if unit_index < 0:
            continue
        start = int(em_data.offsets[unit_index])
        end = int(em_data.offsets[unit_index + 1])
        candidate_indices = em_data.t_indices[start:end].astype(np.int64)
        candidate_scores = (
            em_data.log_liks[start:end].astype(np.float64)
            + np.log(alpha[candidate_indices])
            - np.log(effective_lengths[candidate_indices])
        )
        score_parts = [candidate_scores]
        locus_id = int(unit_locus[unit_index])
        gdna_score = -np.inf
        has_gdna = (not bool(em_data.is_spliced[unit_index])) and np.isfinite(
            float(em_data.gdna_log_liks[unit_index])
        )
        if has_gdna and 0 <= locus_id < len(gdna_alpha):
            gdna_score = float(em_data.gdna_log_liks[unit_index]) + math.log(gdna_alpha[locus_id])
            score_parts.append(np.asarray([gdna_score], dtype=np.float64))
        all_scores = np.concatenate(score_parts)
        max_score = float(np.max(all_scores))
        denom = max_score + math.log(float(np.exp(all_scores - max_score).sum()))
        candidate_post = np.exp(candidate_scores - denom)
        gdna_post = math.exp(gdna_score - denom) if np.isfinite(gdna_score) else 0.0

        synthetic_mask = is_synthetic[candidate_indices]
        synthetic_posts = candidate_post[synthetic_mask]
        origin_mask = t_ids[candidate_indices] == str(row.origin_nrna_id)
        origin_posts = candidate_post[origin_mask]
        max_rna_index = int(np.argmax(candidate_post)) if len(candidate_post) else -1
        max_rna_post = float(candidate_post[max_rna_index]) if max_rna_index >= 0 else 0.0
        max_rna_is_synthetic = bool(synthetic_mask[max_rna_index]) if max_rna_index >= 0 else False
        max_component = "gdna" if gdna_post > max_rna_post else "synthetic_nrna" if max_rna_is_synthetic else "mrna"

        rows.append(
            {
                "frag_id": int(row.frag_id),
                "qname": row.qname,
                "origin_tx": row.origin_tx,
                "assigned_pool": row.assigned_pool,
                "region_class": row.region_class,
                "zc": row.zc,
                "winner_posterior_tag": float(row.posterior),
                "locus_id": locus_id,
                "n_candidates": int(len(candidate_indices)),
                "n_synthetic_candidates": int(synthetic_mask.sum()),
                "origin_nrna_posterior": float(origin_posts[0]) if len(origin_posts) else 0.0,
                "max_synthetic_posterior": float(synthetic_posts.max()) if len(synthetic_posts) else 0.0,
                "sum_synthetic_posterior": float(synthetic_posts.sum()) if len(synthetic_posts) else 0.0,
                "gdna_posterior": gdna_post,
                "max_rna_posterior": max_rna_post,
                "max_component": max_component,
            }
        )

    posterior_df = pd.DataFrame(rows)
    posterior_df.to_csv(diag_dir / "true_nrna_reconstructed_posteriors.tsv", sep="\t", index=False)

    summary = (
        posterior_df.groupby(["assigned_pool", "region_class"], dropna=False)
        .agg(
            n=("frag_id", "size"),
            max_origin_nrna_posterior=("origin_nrna_posterior", "max"),
            q99_origin_nrna_posterior=("origin_nrna_posterior", lambda values: values.quantile(0.99)),
            median_origin_nrna_posterior=("origin_nrna_posterior", "median"),
            max_synthetic_posterior=("max_synthetic_posterior", "max"),
            q99_sum_synthetic_posterior=("sum_synthetic_posterior", lambda values: values.quantile(0.99)),
            median_sum_synthetic_posterior=("sum_synthetic_posterior", "median"),
            max_sum_synthetic_posterior=("sum_synthetic_posterior", "max"),
            n_synthetic_map=("max_component", lambda values: int((values == "synthetic_nrna").sum())),
            n_gdna_map=("max_component", lambda values: int((values == "gdna").sum())),
            n_mrna_map=("max_component", lambda values: int((values == "mrna").sum())),
        )
        .reset_index()
    )
    summary.to_csv(diag_dir / "true_nrna_reconstructed_posterior_summary.tsv", sep="\t", index=False)

    top = posterior_df.sort_values("max_synthetic_posterior", ascending=False).head(50)
    top.to_csv(diag_dir / "true_nrna_top_synthetic_posteriors.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--diag-dir", type=Path, default=DEFAULT_DIAG_DIR)
    parser.add_argument("--fractional-dir", type=Path, default=DEFAULT_FRACTIONAL_DIR)
    args = parser.parse_args()

    base = args.base.resolve()
    diag_dir = args.diag_dir.resolve()
    fragment_rows = pd.read_csv(diag_dir / "fragment_rows.tsv", sep="\t")
    examples = pd.read_csv(diag_dir / "example_true_nrna_fragments.tsv", sep="\t")

    index, estimator, em_data, multi_loci, prior_table, calibration = rebuild_scored_fragments(
        base,
        args.condition,
    )
    del calibration
    unit_lookup = build_unit_lookup(em_data)
    t_df = index.t_df.reset_index(drop=True)

    true_nrna = fragment_rows[fragment_rows["origin"] == "nrna"].copy()
    true_nrna = true_nrna[true_nrna["assigned_pool"].isin(["gdna", "mrna", "nrna"])]
    candidate_rows = [
        candidate_summary_for_fragment(row, em_data, unit_lookup, t_df)
        for row in true_nrna.itertuples(index=False)
    ]
    candidate_df = pd.DataFrame(candidate_rows)
    candidate_df.to_csv(diag_dir / "true_nrna_candidate_presence.tsv", sep="\t", index=False)

    grouped = (
        candidate_df.groupby(["assigned_pool", "region_class"], dropna=False)
        .agg(
            n=("frag_id", "size"),
            origin_nrna_present=("origin_nrna_candidate_present", "sum"),
            any_synthetic=("n_synthetic_candidates", lambda values: int((values > 0).sum())),
            has_gdna=("has_gdna_candidate", "sum"),
            median_rna_candidates=("n_rna_candidates", "median"),
            median_synthetic_candidates=("n_synthetic_candidates", "median"),
            median_mrna_candidates=("n_mrna_candidates", "median"),
        )
        .reset_index()
    )
    grouped["origin_nrna_present_fraction"] = grouped.apply(
        lambda row: safe_div(float(row.origin_nrna_present), float(row.n)), axis=1
    )
    grouped["any_synthetic_fraction"] = grouped.apply(
        lambda row: safe_div(float(row.any_synthetic), float(row.n)), axis=1
    )
    grouped["has_gdna_fraction"] = grouped.apply(
        lambda row: safe_div(float(row.has_gdna), float(row.n)), axis=1
    )
    grouped.to_csv(diag_dir / "true_nrna_candidate_presence_summary.tsv", sep="\t", index=False)

    candidate_examples = example_candidate_rows(
        examples.head(18),
        em_data,
        unit_lookup,
        t_df,
        estimator.effective_lengths,
    )
    candidate_examples.to_csv(diag_dir / "example_candidate_sets.tsv", sep="\t", index=False)

    prior_df = pd.DataFrame(
        {
            "locus_id": np.arange(len(prior_table.gdna_prior_count)),
            "gdna_prior_count": prior_table.gdna_prior_count,
            "enable_gdna": prior_table.enable_gdna.astype(bool),
        }
    )
    prior_df.to_csv(diag_dir / "rebuilt_locus_prior_counts.tsv", sep="\t", index=False)

    fractional_dir = args.fractional_dir.resolve()
    if (fractional_dir / "loci.tsv").exists():
        reconstruct_posterior_distribution(
            fragment_rows,
            em_data,
            unit_lookup,
            t_df,
            estimator.effective_lengths,
            multi_loci,
            fractional_dir,
            diag_dir,
        )

    print(f"wrote candidate diagnostics to {diag_dir}")
    print(grouped.to_string(index=False))
    print("\nExample candidate sets:")
    print(candidate_examples.head(80).to_string(index=False))


if __name__ == "__main__":
    main()