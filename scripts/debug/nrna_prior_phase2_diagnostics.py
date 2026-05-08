"""Diagnose Phase 2 Bayesian-prior behavior in nRNA-heavy scenarios.

This script avoids heredoc/stdin execution so Copilot/VS Code terminal runs
remain reliable. It rebuilds the oracle nRNA double-counting scenario, runs the
pipeline under a few prior variants, and prints the quantities needed to decide
whether failures are caused by global gDNA prior mass, gDNA likelihood
competition, or something outside the prior path.
"""

from __future__ import annotations

import argparse
import dataclasses
import tempfile
from pathlib import Path

import numpy as np

import rigel.pipeline as pipeline_mod
from rigel.calibration.locus_prior import PriorTable
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

SIM_SEED = 42
PIPELINE_SEED = 42
N_FRAGMENTS = 2000


def _make_scenario(work_dir: Path) -> Scenario:
    scenario = Scenario(
        "nrna_double_count",
        genome_length=20000,
        seed=SIM_SEED,
        work_dir=work_dir / "nrna_double_count",
    )
    scenario.add_gene(
        "g1",
        "+",
        [
            {
                "t_id": "t1",
                "exons": [(2000, 4000), (8000, 10000)],
                "abundance": 100,
            }
        ],
    )
    scenario.add_gene(
        "g_ctrl",
        "-",
        [
            {
                "t_id": "t_ctrl",
                "exons": [(14000, 16000), (18000, 19000)],
                "abundance": 0,
            }
        ],
    )
    return scenario


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def _zero_gdna_prior(prior_table: PriorTable) -> PriorTable:
    zeros = np.zeros_like(prior_table.alpha_gdna, dtype=np.float64)
    multilocus_priors = tuple(
        dataclasses.replace(mlp, gdna_prior_count=0.0)
        for mlp in prior_table.multi_locus_priors
    )
    return dataclasses.replace(
        prior_table,
        multi_locus_priors=multilocus_priors,
        alpha_gdna=zeros.copy(),
        gdna_prior_count=zeros.copy(),
    )


def _disable_gdna(prior_table: PriorTable) -> PriorTable:
    disabled = np.zeros_like(prior_table.enable_gdna, dtype=np.uint8)
    return dataclasses.replace(prior_table, enable_gdna=disabled)


def _with_intergenic_only_density(args, kwargs):
    args = list(args)
    if "global_densities" in kwargs:
        global_densities = kwargs["global_densities"]
        replace_kwargs = kwargs.copy()
        replace_args = args
    else:
        global_densities = args[4]
        replace_kwargs = kwargs
        replace_args = args.copy()

    zero_intron = dataclasses.replace(
        global_densities.intron,
        rho=0.0,
        n_fragments=0,
    )
    zero_boundary = dataclasses.replace(
        global_densities.exon_intron,
        rho=0.0,
        n_fragments=0,
    )
    intergenic_only = dataclasses.replace(
        global_densities,
        intron=zero_intron,
        exon_intron=zero_boundary,
    )
    if "global_densities" in replace_kwargs:
        replace_kwargs["global_densities"] = intergenic_only
    else:
        replace_args[4] = intergenic_only
    return tuple(replace_args), replace_kwargs


def _make_prior_wrapper(original, variant: str):
    def wrapper(*args, **kwargs):
        if variant == "intergenic_only":
            args, kwargs = _with_intergenic_only_density(args, kwargs)
        prior_table = original(*args, **kwargs)
        if variant == "current":
            return prior_table
        if variant == "intergenic_only":
            return prior_table
        if variant == "alpha0":
            return _zero_gdna_prior(prior_table)
        if variant == "disabled":
            return _disable_gdna(prior_table)
        if variant == "alpha0_disabled":
            return _disable_gdna(_zero_gdna_prior(prior_table))
        raise ValueError(f"Unknown variant: {variant}")

    return wrapper


def _run_variant(result, variant: str, assignment_mode: str):
    original = pipeline_mod.assemble_priors
    pipeline_mod.assemble_priors = _make_prior_wrapper(original, variant)
    try:
        config = PipelineConfig(
            em=EMConfig(
                seed=PIPELINE_SEED,
                assignment_mode=assignment_mode,
            ),
            scan=BamScanConfig(sj_strand_tag="auto"),
            emit_locus_stats=True,
        )
        pipeline_result = pipeline_mod.run_pipeline(
            result.bam_path,
            result.index,
            config=config,
        )
    finally:
        pipeline_mod.assemble_priors = original

    benchmark = run_benchmark(result, pipeline_result, scenario_name=variant)
    total_rna_expected = benchmark.total_expected + benchmark.n_nrna_expected
    total_rna_observed = benchmark.total_rna_observed
    rna_rel_err = (
        abs(total_rna_observed - total_rna_expected) / total_rna_expected
        if total_rna_expected > 0
        else 0.0
    )
    prior_df = pipeline_result.calibration.multi_locus_prior_df
    cal_summary = pipeline_result.calibration.to_summary_dict()
    densities = cal_summary["global_densities"]

    return {
        "variant": variant,
        "assignment": assignment_mode,
        "truth_gdna": benchmark.n_gdna_expected,
        "truth_nrna": benchmark.n_nrna_expected,
        "pipe_gdna": benchmark.n_gdna_pipeline,
        "pipe_nrna": benchmark.n_nrna_pipeline,
        "total_rna_expected": total_rna_expected,
        "total_rna_observed": total_rna_observed,
        "rna_rel_err": rna_rel_err,
        "alpha_gdna_sum": float(pipeline_result.calibration.alpha_gdna.sum()),
        "alpha_rna_sum": float(pipeline_result.calibration.alpha_rna.sum()),
        "gdna_prior_count_sum": float(prior_df["gdna_prior_count"].sum())
        if len(prior_df) else 0.0,
        "enable_gdna_sum": int(pipeline_result.calibration.prior_table.enable_gdna.sum()),
        "n_loci": int(len(prior_df)),
        "rho_intergenic": float(densities["INTERGENIC"]["rho"]),
        "n_intergenic_density": int(densities["INTERGENIC"]["n_fragments"]),
        "rho_intron": float(densities["INTRON"]["rho"]),
        "n_intron_density": int(densities["INTRON"]["n_fragments"]),
        "rho_boundary": float(densities["EXON-INTRON"]["rho"]),
        "n_boundary_density": int(densities["EXON-INTRON"]["n_fragments"]),
        "gdna_em_count": float(pipeline_result.estimator.gdna_em_count),
        "nrna_em_count": float(pipeline_result.estimator.nrna_em_count),
        "n_intergenic_stats": int(pipeline_result.stats.n_intergenic),
        "strand_specificity": float(pipeline_result.strand_models.strand_specificity),
        "locus_results": pipeline_result.estimator.locus_results,
        "prior_df": prior_df,
        "per_locus_gdna_df": pipeline_result.calibration.per_locus_gdna_df,
    }


def _print_result(case_label: str, row: dict) -> None:
    print(
        "RESULT\t"
        f"case={case_label}\tvariant={row['variant']}\tassign={row['assignment']}\t"
        f"truth_gdna={row['truth_gdna']}\ttruth_nrna={row['truth_nrna']}\t"
        f"pipe_gdna={row['pipe_gdna']:.3f}\tpipe_nrna={row['pipe_nrna']:.3f}\t"
        f"rna_obs={row['total_rna_observed']:.3f}\trna_exp={row['total_rna_expected']}\t"
        f"rna_rel_err={row['rna_rel_err']:.4f}\t"
        f"alpha_gdna_sum={row['alpha_gdna_sum']:.6f}\t"
        f"alpha_rna_sum={row['alpha_rna_sum']:.6f}\t"
        f"enable_gdna_sum={row['enable_gdna_sum']}\t"
        f"rho_ig={row['rho_intergenic']:.8g}\tn_ig={row['n_intergenic_density']}\t"
        f"rho_in={row['rho_intron']:.8g}\tn_in={row['n_intron_density']}\t"
        f"rho_b={row['rho_boundary']:.8g}\tn_b={row['n_boundary_density']}\t"
        f"stats_intergenic={row['n_intergenic_stats']}\t"
        f"ss_est={row['strand_specificity']:.6f}"
    )

    for locus_result in row["locus_results"]:
        print(
            "LOCUS\t"
            f"case={case_label}\tvariant={row['variant']}\t"
            f"id={locus_result['locus_id']}\t"
            f"n_t={locus_result['n_transcripts']}\t"
            f"n_em={locus_result['n_em_fragments']}\t"
            f"rna_total={locus_result['rna_total']:.3f}\t"
            f"gdna={locus_result['gdna']:.3f}\t"
            f"alpha_gdna={locus_result['alpha_gdna']:.6f}\t"
            f"alpha_rna={locus_result['alpha_rna']:.6f}"
        )

    if row["variant"] == "current":
        print("PRIOR_DF")
        print(row["prior_df"].to_string(index=False))
        print("PER_LOCUS_GDNA_DF")
        print(row["per_locus_gdna_df"].to_string(index=False))


def _parse_case(text: str) -> tuple[int, int, float]:
    parts = text.split(",")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("cases must be gdna,nrna,ss")
    return int(parts[0]), int(parts[1]), float(parts[2])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        dest="cases",
        action="append",
        type=_parse_case,
        help="Case as gdna,nrna,ss. May be repeated.",
    )
    parser.add_argument(
        "--variant",
        dest="variants",
        action="append",
        choices=(
            "current",
            "alpha0",
            "intergenic_only",
            "disabled",
            "alpha0_disabled",
        ),
        help="Prior variant. May be repeated.",
    )
    parser.add_argument(
        "--assignment-mode",
        default="sample",
        choices=("sample", "fractional", "map"),
    )
    args = parser.parse_args()

    cases = args.cases or [(0, 30, 1.0), (0, 70, 1.0), (0, 30, 0.9), (20, 70, 0.65)]
    variants = args.variants or ["current", "alpha0", "intergenic_only", "disabled"]

    for gdna_abundance, nrna_abundance, strand_specificity in cases:
        case_label = f"g{gdna_abundance}_n{nrna_abundance}_s{int(strand_specificity * 100)}"
        with tempfile.TemporaryDirectory() as tmp:
            scenario = _make_scenario(Path(tmp))
            scenario_result = scenario.build_oracle(
                n_fragments=N_FRAGMENTS,
                sim_config=SimConfig(
                    frag_mean=200,
                    frag_std=30,
                    frag_min=80,
                    frag_max=450,
                    read_length=100,
                    strand_specificity=strand_specificity,
                    seed=SIM_SEED,
                ),
                gdna_config=_gdna_config(gdna_abundance),
                nrna_abundance=nrna_abundance,
            )
            try:
                print(f"CASE\t{case_label}")
                for variant in variants:
                    row = _run_variant(
                        scenario_result,
                        variant=variant,
                        assignment_mode=args.assignment_mode,
                    )
                    _print_result(case_label, row)
            finally:
                scenario.cleanup()


if __name__ == "__main__":
    main()
