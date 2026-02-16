"""Diagnostic: why does a pure-RNA single-exon scenario lose fragments to gDNA?"""
import tempfile
from pathlib import Path

from hulkrna.sim import Scenario, SimConfig, run_benchmark
from hulkrna.pipeline import run_pipeline

tmp = Path(tempfile.mkdtemp())
sc = Scenario("diag", genome_length=5000, seed=42, work_dir=tmp / "diag")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
])
sim_config = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80,
    frag_max=450, read_length=100, seed=42,
)
result = sc.build(n_fragments=100, sim_config=sim_config)
pr = run_pipeline(result.bam_path, result.index, sj_strand_tag="ts")
bench = run_benchmark(result, pr, scenario_name="diag")

print(bench.summary())
print()
print(f"n_intergenic={pr.stats.n_intergenic}")
print(f"shadow_init={pr.counter.shadow_init}")
print(f"gdna_total={pr.counter.gdna_total}")
print(f"gdna_em_count={pr.counter.gdna_em_count}")
print(f"gdna_unique_count={pr.counter.gdna_unique_count}")
print(f"unique_counts sum={pr.counter.unique_counts.sum()}")
print(f"em_counts sum={pr.counter.em_counts.sum()}")
print(f"t_counts sum={pr.counter.t_counts.sum()}")

sm = pr.strand_models
print(f"\nexonic p_r1_sense={sm.exonic.p_r1_sense:.4f}, n_obs={sm.exonic.n_observations}")
print(f"intergenic p_r1_sense={sm.intergenic.p_r1_sense:.4f}, n_obs={sm.intergenic.n_observations}")

rna_model = sm.model_for_category(0)
print(f"model_for_category → p_r1_sense={rna_model.p_r1_sense:.4f}")

if hasattr(pr.counter, '_converged_theta') and pr.counter._converged_theta is not None:
    theta = pr.counter._converged_theta
    print(f"\ntheta[transcript]={theta[0]:.6f}")
    print(f"theta[gDNA shadow]={theta[1]:.6f}")

sc.cleanup()
