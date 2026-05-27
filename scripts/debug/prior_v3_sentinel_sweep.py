"""Sweep grouped prior edge/max settings on v3 sentinel scenarios."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

PIPELINE_SEED = 42
SIM_SEED = 42


@dataclass(frozen=True)
class ScenarioSpec:
    name: str
    n_fragments: int
    gdna_abundance: float = 0.0
    nrna_abundance: float = 0.0
    strand_specificity: float = 1.0


def sim_config(strand_specificity: float) -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=SIM_SEED,
    )


def gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def build_scenario(spec: ScenarioSpec, work_dir) -> Scenario:
    if spec.name.startswith("single_exon"):
        scenario = Scenario("single_exon", genome_length=8000, seed=SIM_SEED, work_dir=work_dir)
        scenario.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
        scenario.add_gene(
            "g_helper",
            "+",
            [{"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50}],
        )
        scenario.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0}])
        return scenario
    if spec.name.startswith("wide_intron"):
        scenario = Scenario("wide_intron", genome_length=6000, seed=SIM_SEED, work_dir=work_dir)
        scenario.add_gene("g1", "+", [{"t_id": "t1", "exons": [(1000, 2000), (3000, 4000)], "abundance": 100}])
        scenario.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(5000, 5300)], "abundance": 0}])
        return scenario
    if spec.name.startswith("gdna_diag"):
        scenario = Scenario("gdna_diag", genome_length=10_000, seed=SIM_SEED, work_dir=work_dir)
        scenario.add_gene(
            "G1",
            "+",
            [{"t_id": "T1", "exons": [(1000, 2000), (4000, 5000)], "abundance": 100.0}],
        )
        return scenario
    if spec.name.startswith("nonoverlap"):
        scenario = Scenario("non_overlapping", genome_length=10000, seed=SIM_SEED, work_dir=work_dir)
        scenario.add_gene("g1", "+", [{"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100}])
        scenario.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4400)], "abundance": 100}])
        scenario.add_gene("g_ctrl", "+", [{"t_id": "t_ctrl", "exons": [(8000, 8300)], "abundance": 0}])
        return scenario
    if spec.name.startswith("nrna_double"):
        scenario = Scenario("nrna_double_count", genome_length=20000, seed=SIM_SEED, work_dir=work_dir)
        scenario.add_gene(
            "g1",
            "+",
            [{
                "t_id": "t1",
                "exons": [
                    (2000, 2500),
                    (3000, 3500),
                    (4000, 4500),
                    (5000, 5500),
                    (6000, 6500),
                    (7000, 7500),
                    (8000, 8500),
                    (9500, 10000),
                ],
                "abundance": 100,
            }],
        )
        scenario.add_gene(
            "g_ctrl",
            "-",
            [{
                "t_id": "t_ctrl",
                "exons": [
                    (12000, 12500),
                    (13000, 13500),
                    (14000, 14500),
                    (15000, 15500),
                    (16000, 16500),
                    (17000, 17500),
                    (18000, 18500),
                    (19000, 19500),
                ],
                "abundance": 0,
            }],
        )
        return scenario
    raise ValueError(f"unknown spec {spec.name}")


def run_one(
    spec: ScenarioSpec,
    *,
    edge: float,
    cap: float,
    strength: float,
    bias: float,
) -> dict[str, float]:
    with tempfile.TemporaryDirectory(prefix=f"rigel_prior_{spec.name}_") as tmp:
        scenario = build_scenario(spec, tmp)
        result = scenario.build_oracle(
            n_fragments=spec.n_fragments,
            sim_config=sim_config(spec.strand_specificity),
            gdna_config=gdna_config(spec.gdna_abundance),
            nrna_abundance=spec.nrna_abundance,
        )
        config = PipelineConfig(
            em=EMConfig(
                seed=PIPELINE_SEED,
                aggregate_prior_edge_count=edge,
                aggregate_prior_max_count=cap,
                aggregate_prior_strength=strength,
                gdna_prior_logit_bias=bias,
            ),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)
        bench = run_benchmark(result, pr, scenario_name=spec.name)
        first = next(t for t in bench.transcripts if t.t_id not in {"t_ctrl", "t_helper"})
        loci = pr.estimator.get_loci_df(result.index)
        alpha_rna = float(loci["alpha_rna_add"].sum()) if "alpha_rna_add" in loci else 0.0
        alpha_gdna = float(loci["alpha_gdna_add"].sum()) if "alpha_gdna_add" in loci else 0.0
        gdna_eff = float(loci["gdna_eff_len"].sum()) if "gdna_eff_len" in loci else 0.0
        gdna_eff_unw = (
            float(loci["gdna_eff_len_unweighted"].sum())
            if "gdna_eff_len_unweighted" in loci
            else 0.0
        )
        gdna_weight = (
            float(loci["gdna_em_exposure_weight"].mean())
            if "gdna_em_exposure_weight" in loci and len(loci) > 0
            else 0.0
        )
        return {
            "observed": float(first.observed),
            "expected": float(first.expected),
            "diff": float(first.observed - first.expected),
            "gdna_pipeline": float(bench.n_gdna_pipeline),
            "gdna_expected": float(bench.n_gdna_expected),
            "nrna_pipeline": float(bench.n_nrna_pipeline),
            "alpha_rna_sum": alpha_rna,
            "alpha_gdna_sum": alpha_gdna,
            "gdna_eff_len": gdna_eff,
            "gdna_eff_len_unweighted": gdna_eff_unw,
            "gdna_weight": gdna_weight,
        }


def main() -> None:
    specs = [
        ScenarioSpec("single_exon_g0", n_fragments=500),
        ScenarioSpec("wide_intron_g0", n_fragments=500),
        ScenarioSpec("gdna_diag_g0", n_fragments=1000),
        ScenarioSpec("nonoverlap_g0", n_fragments=500),
        ScenarioSpec("nrna_double_g0_n0_s65", n_fragments=2000, strand_specificity=0.65),
        ScenarioSpec("nrna_double_g0_n0_s90", n_fragments=2000, strand_specificity=0.9),
        ScenarioSpec("nrna_double_g0_n0_s100", n_fragments=2000, strand_specificity=1.0),
        ScenarioSpec("single_exon_g20", n_fragments=500, gdna_abundance=20),
        ScenarioSpec("wide_intron_g20", n_fragments=500, gdna_abundance=20),
        ScenarioSpec("gdna_diag_g50", n_fragments=1000, gdna_abundance=50),
    ]
    configs = [
        (5.0, 10.0, 1.0, 0.0),
        (1000.0, 1000.0, 1.0, -2.0),
        (1000.0, 1000.0, 3.0, -2.0),
        (1000.0, 1000.0, 3.0, -4.0),
        (1000.0, 1000.0, 3.0, -6.0),
        (1000.0, 1000.0, 5.0, -4.0),
    ]
    print("edge\tcap\tstrength\tbias\tscenario\tobs\texp\tdiff\tgdna\tgdna_exp\tnrna\talpha_rna\talpha_gdna\tgdna_eff\tgdna_eff_unw\tgdna_weight")
    for edge, cap, strength, bias in configs:
        for spec in specs:
            metrics = run_one(spec, edge=edge, cap=cap, strength=strength, bias=bias)
            print(
                f"{edge:g}\t{cap:g}\t{strength:g}\t{bias:g}\t{spec.name}\t"
                f"{metrics['observed']:.1f}\t{metrics['expected']:.1f}\t{metrics['diff']:.1f}\t"
                f"{metrics['gdna_pipeline']:.1f}\t{metrics['gdna_expected']:.1f}\t"
                f"{metrics['nrna_pipeline']:.1f}\t{metrics['alpha_rna_sum']:.1f}\t{metrics['alpha_gdna_sum']:.1f}\t"
                f"{metrics['gdna_eff_len']:.1f}\t{metrics['gdna_eff_len_unweighted']:.1f}\t{metrics['gdna_weight']:.4g}"
            )


if __name__ == "__main__":
    main()
