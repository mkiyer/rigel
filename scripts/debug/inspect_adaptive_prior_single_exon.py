from __future__ import annotations

from pathlib import Path

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario, run_benchmark

SIM_SEED = 42
PIPELINE_SEED = 42


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

work_dir = Path("/tmp/rigel_inspect_adaptive_prior_single_exon")
sc = Scenario("single_exon", genome_length=8000, seed=SIM_SEED, work_dir=work_dir)
sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
sc.add_gene(
    "g_helper",
    "+",
    [{"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50}],
)
sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0}])
result = sc.build_oracle(n_fragments=500, sim_config=sim_config(), gdna_config=None, nrna_abundance=0)
for bias in (0.5, 0.9, 0.99, 0.9999):
    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED, rna_call_bias=bias),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=f"debug_single_exon_base_bias_{bias}")
    print("\n=== bias", bias, "===")
    print(bench.summary())
    print("\nLoci df:")
    loci_df = pr.estimator.get_loci_df(result.index)
    print(
        loci_df[
            [
                "locus_id",
                "n_em_fragments",
                "mrna",
                "gdna",
                "alpha_gdna_add",
                "alpha_rna_add",
                "prior_rna_share_final",
                "enable_gdna",
            ]
        ].to_string(index=False)
    )
    print("\nPrior policy:")
    print(pr.calibration.prior_table.to_summary_dict())
    print("\nRegion p_states and prior mass:")
    rc = pr.calibration.region_calibration
    print(rc.p_states)
    print(rc.region_unspliced_mass.total_mass)
print("region df:")
print(result.index.region_df[["region_id", "ref_name", "start", "end", "signature"]].to_string(index=False))
