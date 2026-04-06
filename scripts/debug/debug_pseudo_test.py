"""Debug script for TestPseudogene::test_parent_only failure."""
import logging
import pathlib
import tempfile

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, SimConfig, run_benchmark

SIM_SEED = 42
logging.basicConfig(level=logging.WARNING)

sim_config = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=1.0, seed=SIM_SEED,
)

with tempfile.TemporaryDirectory() as tmp:
    tmp = pathlib.Path(tmp)
    sc = Scenario("pseudogene_parent_only", genome_length=20000, seed=SIM_SEED,
                  work_dir=tmp / "pseudogene_parent_only")
    sc.add_gene("g_parent", "+", [
        {"t_id": "t_parent", "exons": [(1000, 1300), (2000, 2300)], "abundance": 100}
    ])
    sc.add_gene("g_pseudo", "+", [
        {"t_id": "t_pseudo", "exons": [(8000, 8600)], "abundance": 0}
    ])
    exon1 = sc.genome[1000:1300]
    exon2 = sc.genome[2000:2300]
    sc.genome.edit(8000, exon1 + exon2)
    sc.add_gene("g_helper", "+", [
        {"t_id": "t_helper", "exons": [(14000, 14300), (14700, 15000)], "abundance": 80}
    ])
    sc.add_gene("g_ctrl", "-", [
        {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0}
    ])

    result = sc.build(
        n_fragments=1000,
        sim_config=sim_config,
        gdna_config=None,
        nrna_abundance=0,
    )

    for em_mode in ("vbem", "map"):
        config = PipelineConfig(
            em=EMConfig(seed=42, mode=em_mode),
            scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)
        est = pr.estimator
        obs = est.t_counts.sum(axis=1)
        print(f"\n=== mode={em_mode} ===")
        for lr in est.locus_results:
            print(
                f"  locus {lr['locus_id']}: n_t={lr['n_transcripts']} "
                f"gdna_init={lr['gdna_init']:.4f} mrna={lr['mrna']:.1f} "
                f"gdna={lr['gdna']:.1f}"
            )
        print(f"t_parent={obs[0]:.1f} t_pseudo={obs[1]:.1f} "
              f"combined={obs[0]+obs[1]:.1f} n_gdna_em={est.gdna_em_count:.1f}")
