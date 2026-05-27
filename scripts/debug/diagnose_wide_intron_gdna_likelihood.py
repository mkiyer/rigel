from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from rigel.calibration import calibrate
from rigel.calibration.prior import assemble_priors
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
from rigel.locus import build_multi_loci
from rigel.locus_partition import partition_and_free
from rigel.pipeline import (
    _assign_locus_ids,
    _calibration_strand_summary,
    _native_detect_sj_tag,
    _score_fragments,
    _setup_geometry_and_estimator,
    run_pipeline,
    scan_and_buffer,
)
from rigel.scoring import FragmentScorer
from rigel.sim import ReadSimConfig, Scenario, run_benchmark

sys.path.append(str(Path(__file__).resolve().parents[2]))
from scripts.debug.diagnose_single_exon_gdna_likelihood import (
    build_local_arrays,
    print_region_audit,
    row_for_unit,
    snapshot_buffer_features,
    summarize,
    transcript_index_by_id,
    transcript_label,
    unit_log_liks_with_global_rna_fl,
    vbem_fixed_point,
)

SIM_SEED = 42
PIPELINE_SEED = 42
WORK_DIR = Path("/tmp/rigel_diagnose_wide_intron_gdna_likelihood")
FRAGMENT_TSV = Path("/tmp/rigel_wide_intron_gdna_likelihood_fragments.tsv")


def sim_config() -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SIM_SEED,
    )


def build_scenario():
    sc = Scenario("wide_intron", genome_length=6000, seed=SIM_SEED, work_dir=WORK_DIR)
    sc.add_gene(
        "g1",
        "+",
        [{"t_id": "t1", "exons": [(1000, 2000), (3000, 4000)], "abundance": 100}],
    )
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(5000, 5300)], "abundance": 0}])
    return sc.build_oracle(
        n_fragments=500,
        sim_config=sim_config(),
        gdna_config=None,
        nrna_abundance=0,
    )


def main() -> None:
    result = build_scenario()
    index = result.index
    sample_config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED, assignment_mode="sample", rna_call_bias=0.5),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    sample_pr = run_pipeline(result.bam_path, index, config=sample_config)
    bench = run_benchmark(result, sample_pr, scenario_name="diagnose_wide_intron_sample")
    print("=== Reproduction with default sample assignment ===")
    print(bench.summary())
    print(
        sample_pr.estimator.get_loci_df(index)[
            [
                "locus_id",
                "n_em_fragments",
                "mrna",
                "nrna",
                "gdna",
                "alpha_gdna_add",
                "alpha_rna_add",
                "prior_rna_share_final",
                "enable_gdna",
            ]
        ].to_string(index=False)
    )

    em_config = EMConfig(seed=PIPELINE_SEED, assignment_mode="fractional", rna_call_bias=0.5)
    scan = BamScanConfig(sj_strand_tag="auto")
    if scan.sj_strand_tag == "auto":
        scan = replace(scan, sj_strand_tag=_native_detect_sj_tag(str(result.bam_path)))

    stats, strand_models, frag_length_models, buffer, payload = scan_and_buffer(
        str(result.bam_path), index, scan
    )
    strand_models.finalize()
    feature_by_frag_id = snapshot_buffer_features(buffer)
    calibration = calibrate(
        index=index,
        payload=payload,
        scan_trained=frag_length_models,
        strand_summary=_calibration_strand_summary(strand_models),
    )
    _geometry, estimator = _setup_geometry_and_estimator(index, calibration.fl_models.rna, em_config)
    scoring = FragmentScoringConfig()
    scorer = FragmentScorer.from_models(
        strand_models,
        calibration.fl_models.rna,
        calibration.fl_models.gdna,
        index,
        estimator,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
        pruning_min_posterior=scoring.pruning_min_posterior,
    )
    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        calibration.fl_models.rna,
        calibration.fl_models.gdna,
        stats,
        estimator,
        scoring,
        scan.log_every,
        None,
    )
    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)
    prior_table = assemble_priors(
        multi_loci=multi_loci,
        em_data=em_data,
        index=index,
        calibration=calibration,
        em_config=em_config,
    )
    t1_idx = transcript_index_by_id(index, "t1")
    t1_locus = next(locus for locus in multi_loci if t1_idx in set(map(int, locus.transcript_indices)))
    lid = int(t1_locus.multi_locus_id)
    print("\n=== t1 locus identity ===")
    print("t1_idx", t1_idx, "locus_id", lid)
    print(
        "locus transcript indices:",
        [(int(t), transcript_label(index, int(t))) for t in t1_locus.transcript_indices],
    )
    print(
        "prior:",
        {
            "alpha_gdna_add": float(prior_table.alpha_gdna_add[lid]),
            "alpha_rna_add": float(prior_table.alpha_rna_add[lid]),
            "prior_rna_share_final": float(prior_table.prior_rna_share_final[lid]),
            "gdna_eff_len": float(prior_table.gdna_eff_len[lid]),
            "gdna_eff_len_unweighted": float(prior_table.gdna_eff_len_unweighted[lid]),
            "gdna_em_exposure_weight": float(prior_table.gdna_em_exposure_weight[lid]),
            "enable_gdna": int(prior_table.enable_gdna[lid]),
        },
    )

    partitions = partition_and_free(em_data, multi_loci)
    part = partitions[lid]
    partition_tuple = (
        part.offsets,
        part.t_indices,
        part.log_liks,
        part.coverage_weights,
        part.count_cols,
        part.is_spliced,
        part.gdna_log_liks,
        part.locus_t_indices,
        part.locus_count_cols,
    )
    native_result = estimator.run_batch_locus_em_partitioned(
        [partition_tuple],
        [t1_locus.transcript_indices],
        np.array([prior_table.alpha_gdna_add[lid]], dtype=np.float64),
        index,
        rna_prior_count=np.array([prior_table.alpha_rna_add[lid]], dtype=np.float64),
        gdna_eff_len=np.array([prior_table.gdna_eff_len[lid]], dtype=np.float64),
        enable_gdna=np.array([prior_table.enable_gdna[lid]], dtype=np.uint8),
        em_iterations=em_config.iterations,
        em_convergence_delta=em_config.convergence_delta,
        emit_assignments=True,
    )
    print("\n=== Native fractional EM for t1 locus ===")
    print(
        {
            "native_locus_rna": float(native_result[1][0]),
            "native_locus_gdna": float(native_result[2][0]),
            "native_gdna_fraction": float(native_result[2][0] / part.n_units),
            "native_map_winner_gdna_count": int(np.count_nonzero(native_result[3] == -2)),
        }
    )

    (
        local_t,
        gdna_idx,
        unit_components,
        unit_log_liks,
        log_eff,
        unambig,
        alpha0,
        has_gdna,
    ) = build_local_arrays(
        part,
        t1_locus,
        estimator,
        prior_table.gdna_eff_len[lid],
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
    )
    theta, alpha, final_post, n_iter = vbem_fixed_point(
        unit_components,
        unit_log_liks,
        log_eff,
        unambig,
        alpha0,
        gdna_idx,
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
        has_gdna,
    )
    global_rna_log_liks = unit_log_liks_with_global_rna_fl(
        unit_components,
        unit_log_liks,
        local_t,
        gdna_idx,
        part,
        feature_by_frag_id,
        calibration,
    )
    theta_global, alpha_global, final_post_global, n_iter_global = vbem_fixed_point(
        unit_components,
        global_rna_log_liks,
        log_eff,
        unambig,
        alpha0,
        gdna_idx,
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
        has_gdna,
    )
    global_gdna_post_sum = float(
        sum(
            post[np.flatnonzero(comps == gdna_idx)[0]]
            for comps, post in zip(unit_components, final_post_global, strict=True)
            if np.any(comps == gdna_idx)
        )
    )
    rows = [
        row_for_unit(
            index,
            scorer,
            calibration,
            feature_by_frag_id,
            part,
            local_t,
            gdna_idx,
            unit_components,
            unit_log_liks,
            log_eff,
            theta,
            final_post,
            ui,
        )
        for ui in range(part.n_units)
    ]
    df = pd.DataFrame(rows)
    df.to_csv(FRAGMENT_TSV, sep="\t", index=False)

    print("\n=== Python VBEM check for t1 locus ===")
    print("iterations", n_iter)
    print("theta", {"gDNA": float(theta[gdna_idx]), **{transcript_label(index, t): float(theta[i]) for i, t in enumerate(local_t)}})
    print("alpha", {"gDNA": float(alpha[gdna_idx]), **{transcript_label(index, t): float(alpha[i]) for i, t in enumerate(local_t)}})
    print("sum p_gdna", float(df["gdna_post"].sum()), "mean", float(df["gdna_post"].mean()))
    print(
        "counterfactual RNA=global FL",
        {
            "iterations": int(n_iter_global),
            "theta_gdna": float(theta_global[gdna_idx]),
            "alpha_gdna": float(alpha_global[gdna_idx]),
            "sum_p_gdna": global_gdna_post_sum,
            "mean_p_gdna": global_gdna_post_sum / float(part.n_units),
        },
    )
    print("fragment table", str(FRAGMENT_TSV))

    print("\n=== Score deltas: gDNA minus best RNA ===")
    for col in (
        "raw_delta_gdna_minus_rna",
        "final_delta_with_theta_and_eff",
        "fl_delta",
        "fl_delta_if_rna_used_global_fl",
        "raw_delta_if_rna_used_global_fl",
        "strand_delta",
        "overhang_delta",
        "eff_delta",
        "theta_delta",
        "gdna_post",
    ):
        print(summarize(col, df[col].to_numpy()))
    print(
        "counts",
        {
            "n_units": int(len(df)),
            "raw_gdna_gt_rna": int((df["raw_delta_gdna_minus_rna"] > 0.0).sum()),
            "final_gdna_posterior_gt_0.5": int((df["gdna_post"] > 0.5).sum()),
            "overhang_positive": int((df["overhang"] > 0).sum()),
            "same_strand": int(df["same_strand"].sum()),
            "calc_rna_abs_resid_max": float(df["calc_resid_rna"].abs().max()),
            "calc_gdna_abs_resid_max": float(df["calc_resid_gdna"].abs().max()),
        },
    )
    show_cols = [
        "unit",
        "frag_id",
        "winner",
        "gdna_post",
        "raw_delta_gdna_minus_rna",
        "fl_delta",
        "strand_delta",
        "eff_delta",
        "theta_delta",
        "frag_length_rna",
        "genomic_footprint",
        "overhang",
        "same_strand",
        "splice_type",
        "genomic_start",
    ]
    print("\nTop 20 gDNA-favored fragments by raw score:")
    print(df.sort_values("raw_delta_gdna_minus_rna", ascending=False)[show_cols].head(20).to_string(index=False))
    print("\nFL model summary:")
    print(calibration.fl_models.to_summary_dict())
    target_regions = list(range(min(len(index.region_df), 6)))
    print_region_audit(index, payload, calibration, target_regions)


if __name__ == "__main__":
    main()
