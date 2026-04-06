"""Debug distinguishable paralogs test to understand fragment accountability."""
import tempfile
import logging
from pathlib import Path
import numpy as np

logging.basicConfig(level=logging.WARNING)

from rigel.sim import Scenario, SimConfig, run_benchmark
from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.calibration import calibrate_gdna, compute_region_stats, compute_sense_fraction

SIM_SEED = 42
tmpdir = Path(tempfile.mkdtemp())

sc = Scenario(
    "dist_paralogs_debug",
    genome_length=12000,
    seed=SIM_SEED,
    work_dir=tmpdir / "dist_paralogs_debug",
)
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(500, 800), (1200, 1500)], "abundance": 100},
])
sc.add_gene("g2", "+", [
    {"t_id": "t2", "exons": [(5000, 5300), (5700, 5900)], "abundance": 100},
])
sc.genome.edit(5000, sc.genome[500:800])
sc.add_gene("g_helper", "+", [
    {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50},
])
sc.add_gene("g_ctrl", "-", [
    {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
])

simcfg = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=1.0, seed=SIM_SEED,
)
result = sc.build(
    n_fragments=1000, sim_config=simcfg,
    gdna_config=None, nrna_abundance=0,
)

config = PipelineConfig(
    em=EMConfig(seed=42, mode="map"),
    scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
)
pr_map = run_pipeline(result.bam_path, result.index, config=config)

config_vbem = PipelineConfig(
    em=EMConfig(seed=42, mode="vbem"),
    scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
)
pr_vbem = run_pipeline(result.bam_path, result.index, config=config_vbem)

for label, pr in [("MAP", pr_map), ("VBEM", pr_vbem)]:
    print(f"\n=== {label} MODE ===")
    loci_df = pr.estimator.get_loci_df(result.index)
    print(loci_df[["locus_id", "n_em_fragments", "mrna", "nrna", "gdna",
                    "gdna_rate", "gdna_prior"]].to_string())

# Access calibration from VBEM run
cal = pr_vbem.calibration
print("\n=== CALIBRATION DIAGNOSTICS ===")
print(f"mixing_proportion (pi): {cal.mixing_proportion:.4f}")
print(f"gdna_density_global: {cal.gdna_density_global:.6f}")
print(f"kappa_strand: {cal.kappa_strand:.2f}")
print(f"n_iterations: {cal.n_iterations}")

# Region details
region_df = result.index.region_df
gammas = cal.region_posteriors
n_total = cal.region_n_total

print(f"\n=== REGION GAMMA VALUES ===")
print(f"{'idx':>4} {'ref':>6} {'start':>6} {'end':>6} {'len':>5} "
      f"{'tx_p':>4} {'tx_n':>4} {'n_tot':>6} {'gamma':>6}")
for i in range(len(region_df)):
    row = region_df.iloc[i]
    print(f"{i:4d} {str(row.get('ref','')):>6} "
          f"{row.get('start', '?'):>6} {row.get('end', '?'):>6} "
          f"{int(row['length']):>5} "
          f"{int(row['tx_pos']):>4} {int(row['tx_neg']):>4} "
          f"{n_total[i]:>6.0f} {gammas[i]:>6.3f}")

# Show only regions with gamma > 0 and data
gdna_regions = [(i, gammas[i], n_total[i]) for i in range(len(gammas))
                if gammas[i] > 0 and n_total[i] > 0]
print(f"\n=== REGIONS WITH gamma > 0 AND data ===")
for i, g, n in gdna_regions:
    print(f"  region {i}: gamma={g:.3f}, n_total={n:.0f}")

print(f"\nTotal fragments in gamma>0 regions: "
      f"{sum(g*n for _, g, n in gdna_regions):.0f}")

bench = run_benchmark(result, pr_vbem, scenario_name="debug_dist")
print(f"\n=== BENCHMARK ===")
print(f"total_rna_observed: {bench.total_rna_observed}")
print(f"total_expected: {bench.total_expected}")
print(f"n_gdna_pipeline: {bench.n_gdna_pipeline}")
for t in bench.transcripts:
    print(f"  {t.t_id}: expected={t.expected}, observed={t.observed:.1f}")

# Check locus-level details
print(f"\n=== LOCUS DETAILS ===")
loci_df = pr_vbem.estimator.get_loci_df(result.index)
print(loci_df.to_string())

sc.cleanup()
