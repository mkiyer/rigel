"""Inspect exposure diagnostics for PR06 known-xfail scenarios."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.density_model import fit_density_evidence
from rigel.calibration.density_observation import build_density_observation
from rigel.calibration.integration import fuse_density_and_strand
from rigel.calibration.region_count_ledger import build_region_count_ledger
from rigel.calibration.strand_deconv import (
    build_compartment_strand_counts,
    build_strand_region_counts,
)
from rigel.calibration.strand_summary import StrandSummary
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark


def _pipeline_config() -> PipelineConfig:
    return PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )


def _sim_config(strand_specificity: float) -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=42,
    )


def _print_summary(label: str, result, pipeline_result) -> None:
    bench = run_benchmark(result, pipeline_result, scenario_name=label)
    calibration = pipeline_result.calibration.region_calibration
    mass = calibration.region_unspliced_mass
    exposure = calibration.region_exposure
    quant = pipeline_result.estimator.get_counts_df(result.index)
    nrna = pipeline_result.estimator.get_nrna_counts_df(result.index)
    loci = pipeline_result.estimator.get_loci_df(result.index)
    region_arrays = RegionArrays.from_region_df(result.index.region_df, result.index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(pipeline_result.calibration_payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    observation = build_density_observation(
        region_arrays,
        ledger,
        pipeline_result.calibration.fl_models.gdna,
    )
    strand_summary = StrandSummary.from_model(pipeline_result.strand_models.exonic_spliced)
    strand_counts = build_strand_region_counts(
        region_arrays,
        payload_arrays,
        p_r1_sense=strand_summary.p_r1_sense,
    )
    density_evidence = fit_density_evidence(observation, confidence=0.95)
    fused = fuse_density_and_strand(
        region_arrays=region_arrays,
        ledger=ledger,
        density_observation=observation,
        density_evidence=density_evidence,
        strand_counts=strand_counts,
        strand_summary=strand_summary,
        kappa_d=pipeline_result.calibration.strand_channels.kappa_d,
        confidence=0.95,
    )
    compartment_counts = build_compartment_strand_counts(
        region_arrays,
        payload_arrays,
        p_r1_sense=pipeline_result.calibration.strand_channels.p_r1_sense,
    )
    strand_channels = pipeline_result.calibration.strand_channels
    region = result.index.region_df.iloc[region_arrays.order].reset_index(drop=True).copy()
    region["p_unexpressed"] = calibration.p_states[:, 0]
    region["gdna_mass"] = mass.gdna_mass
    region["rna_mass"] = mass.rna_mass
    region["total_mass"] = mass.total_mass
    region["method"] = mass.method
    region["omega"] = exposure.omega
    region["raw_ratio"] = exposure.raw_ratio
    region["shrink_weight"] = exposure.shrink_weight
    region["support_count"] = exposure.support_count
    region["contained_sense"] = compartment_counts.contained_sense
    region["contained_antisense"] = compartment_counts.contained_antisense
    region["contained_total"] = compartment_counts.contained_total
    region["contained_gdna_mean"] = strand_channels.contained_mean
    region["contained_rna_lower"] = strand_channels.contained_rna_lower
    region["contained_reliability"] = strand_channels.contained_reliability
    region["contained_logbf"] = strand_channels.contained_log_bayes_factor
    region["contained_precision"] = strand_channels.contained_precision
    region["left_sense"] = compartment_counts.boundary_left_sense
    region["left_antisense"] = compartment_counts.boundary_left_antisense
    region["left_total"] = compartment_counts.boundary_left_total
    region["left_gdna_mean"] = strand_channels.boundary_left_mean
    region["left_reliability"] = strand_channels.boundary_left_reliability
    region["right_sense"] = compartment_counts.boundary_right_sense
    region["right_antisense"] = compartment_counts.boundary_right_antisense
    region["right_total"] = compartment_counts.boundary_right_total
    region["right_gdna_mean"] = strand_channels.boundary_right_mean
    region["right_reliability"] = strand_channels.boundary_right_reliability

    gdna_mass = float(np.sum(mass.gdna_mass, dtype=np.float64))
    rna_mass = float(np.sum(mass.rna_mass, dtype=np.float64))
    total_mass = float(np.sum(mass.total_mass, dtype=np.float64))
    print(f"\n=== {label} ===")
    print(bench.summary())
    print(
        "mass: "
        f"gdna={gdna_mass:.3f} rna={rna_mass:.3f} total={total_mass:.3f} "
        f"gdna_share={gdna_mass / total_mass if total_mass > 0.0 else 0.0:.4f}"
    )
    print(
        "exposure: "
        f"tau2={exposure.tau2:.6g} tau2_hat={exposure.tau2_hat:.6g} "
        f"method={exposure.tau2_method} pool={exposure.tau2_pool_size} "
        f"omega_p50={np.quantile(exposure.omega, 0.50):.4f} "
        f"omega_p95={np.quantile(exposure.omega, 0.95):.4f} "
        f"omega_max={np.max(exposure.omega):.4f}"
    )
    strand_model = pipeline_result.strand_models.exonic_spliced
    print(
        "strand model: "
        f"n={strand_model.n_observations} same={strand_model.n_same} "
        f"opposite={strand_model.n_opposite} p_r1_sense={strand_model.p_r1_sense:.4f} "
        f"minor_alpha={strand_model.minor_rate_posterior_alpha:.1f} "
        f"minor_beta={strand_model.minor_rate_posterior_beta:.1f} "
        f"kappa_d={strand_channels.kappa_d:.4f}"
    )
    print(
        "density/strand fusion check: "
        f"rho_ref={density_evidence.rho_ref:.6g} source={density_evidence.rho_ref_source} "
        f"mean_sum={float(np.sum(fused.mean_count, dtype=np.float64)):.3f} "
        f"upper_sum={float(np.sum(fused.upper_count, dtype=np.float64)):.3f} "
        f"strand_used={int(np.count_nonzero(fused.strand_applicable))}"
    )
    print("quant:")
    print(
        quant[
            [
                "transcript_id",
                "effective_length",
                "em_effective_length",
                "em_exposure_factor",
                "count",
                "count_em",
            ]
        ].to_string(index=False)
    )
    if not nrna.empty:
        print("nrna:")
        print(
            nrna[
                [
                    "nrna_id",
                    "effective_length",
                    "em_effective_length",
                    "em_exposure_factor",
                    "count",
                ]
            ].to_string(index=False)
        )
    print("loci:")
    print(
        loci[
            [
                "locus_id",
                "mrna",
                "nrna",
                "gdna",
                "gdna_eff_len_unweighted",
                "gdna_exposure_factor",
                "gdna_eff_len_em",
            ]
        ].to_string(index=False)
    )
    print("regions with support:")
    with pd.option_context("display.max_rows", 80, "display.max_columns", 40):
        print(
            region.loc[
                region["support_count"] > 0,
                [
                    "ref_name",
                    "start",
                    "end",
                    "signature",
                    "p_unexpressed",
                    "gdna_mass",
                    "rna_mass",
                    "total_mass",
                    "method",
                    "omega",
                    "raw_ratio",
                    "shrink_weight",
                    "support_count",
                ],
            ].to_string(index=False)
        )
    print("strand deconvolution by region:")
    with pd.option_context("display.max_rows", 80, "display.max_columns", 80):
        print(
            region.loc[
                region["support_count"] > 0,
                [
                    "ref_name",
                    "start",
                    "end",
                    "signature",
                    "contained_sense",
                    "contained_antisense",
                    "contained_total",
                    "contained_gdna_mean",
                    "contained_rna_lower",
                    "contained_reliability",
                    "contained_logbf",
                    "contained_precision",
                    "left_sense",
                    "left_antisense",
                    "left_total",
                    "left_gdna_mean",
                    "left_reliability",
                    "right_sense",
                    "right_antisense",
                    "right_total",
                    "right_gdna_mean",
                    "right_reliability",
                ],
            ].to_string(index=False)
        )


def _antisense_intronic(work_dir: Path) -> tuple[object, object]:
    scenario = Scenario("anti_intron", genome_length=8000, seed=42, work_dir=work_dir)
    scenario.add_gene(
        "g1",
        "+",
        [
            {
                "t_id": "t1",
                "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
                "abundance": 100,
            }
        ],
    )
    scenario.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0}])
    scenario.add_gene(
        "g_ctrl",
        "+",
        [{"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0}],
    )
    result = scenario.build_oracle(
        n_fragments=2000,
        sim_config=_sim_config(0.65),
        nrna_abundance=50,
    )
    return result, run_pipeline(result.bam_path, result.index, config=_pipeline_config())


def _paralogs(work_dir: Path, gdna_abundance: int) -> tuple[object, object] | None:
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        print("Skipping paralog diagnostics: minimap2/samtools unavailable.")
        return None
    scenario = Scenario(f"paralogs_gdna_{gdna_abundance}", genome_length=12000, seed=42, work_dir=work_dir)
    scenario.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1000)], "abundance": 100}])
    scenario.add_gene("g2", "+", [{"t_id": "t2", "exons": [(5000, 5500)], "abundance": 100}])
    scenario.genome.edit(5000, scenario.genome[500:1000])
    scenario.add_gene(
        "g_helper",
        "+",
        [{"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50}],
    )
    scenario.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0}])
    gdna = GDNAConfig(abundance=gdna_abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
    result = scenario.build(
        n_fragments=3000,
        sim_config=_sim_config(1.0),
        gdna_config=gdna,
    )
    config = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
    )
    return result, run_pipeline(result.bam_path, result.index, config=config)


def main() -> None:
    root = Path(tempfile.mkdtemp(prefix="rigel_pr06_xfail_"))
    anti_result, anti_pipeline = _antisense_intronic(root / "anti")
    _print_summary("anti_intron_ss065_nrna50", anti_result, anti_pipeline)
    for gdna_abundance in (20, 100):
        paralog = _paralogs(root / f"paralog_{gdna_abundance}", gdna_abundance)
        if paralog is not None:
            result, pipeline_result = paralog
            _print_summary(f"paralogs_gdna_{gdna_abundance}", result, pipeline_result)


if __name__ == "__main__":
    main()