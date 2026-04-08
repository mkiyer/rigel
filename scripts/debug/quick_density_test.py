"""Quick test: density floor effect on g20_n0_s65."""
import logging
logging.basicConfig(level=logging.WARNING)

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.sim import Scenario, GDNAConfig, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline

SEED = 42

sc = Scenario("diag", genome_length=20000, seed=SEED, work_dir="/tmp/diag_density2")
sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(2000, 4000), (8000, 10000)], "abundance": 100}])
sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(14000, 16000), (18000, 19000)], "abundance": 0}])

sim = SimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                read_length=100, strand_specificity=0.65, seed=SEED)
gdna = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)

result = sc.build_oracle(n_fragments=2000, sim_config=sim, gdna_config=gdna, nrna_abundance=0)
cfg = PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="auto"))
pr = run_pipeline(result.bam_path, result.index, config=cfg)
bench = run_benchmark(result, pr, scenario_name="density_floor_2")

cal = pr.calibration
gamma = cal.region_e_gdna.sum() / max(cal.region_n_total.sum(), 1.0)
total_rna_exp = bench.total_expected + bench.n_nrna_expected
total_rna_obs = bench.total_rna_observed
err = abs(total_rna_obs - total_rna_exp) / max(total_rna_exp, 1)

print(f"SS={cal.strand_specificity:.3f} gamma={gamma:.4f} lambda={cal.lambda_gdna:.6f}")
print(f"gdna: exp={bench.n_gdna_expected} pipe={bench.n_gdna_pipeline:.0f}")
print(f"nrna: exp={bench.n_nrna_expected} pipe={bench.n_nrna_pipeline:.0f}")
print(f"Total RNA: exp={total_rna_exp} obs={total_rna_obs:.0f} err={err:.4f}")
print(f"Per-region e_gdna: {cal.region_e_gdna}")

sc.cleanup()
