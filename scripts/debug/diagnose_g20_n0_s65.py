"""Diagnose the g20_n0_s65 test failure in test_nrna_double_counting.

Reproduces the exact scenario, dumps calibration, prior, and EM details.
"""

import logging
import numpy as np

from rigel.config import EMConfig, PipelineConfig, BamScanConfig, CalibrationConfig
from rigel.sim import Scenario, GDNAConfig, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline
from rigel.calibration import calibrate_gdna
from rigel.locus import compute_locus_priors

logging.basicConfig(level=logging.WARNING)

# ── Reproduce the exact scenario ──
SEED = 42
N_FRAGMENTS = 2000

sc = Scenario(
    "diag_g20_n0_s65",
    genome_length=20000,
    seed=SEED,
    work_dir="/tmp/diag_g20_n0_s65",
)
sc.add_gene("g1", "+", [{
    "t_id": "t1",
    "exons": [(2000, 4000), (8000, 10000)],
    "abundance": 100,
}])
sc.add_gene("g_ctrl", "-", [{
    "t_id": "t_ctrl",
    "exons": [(14000, 16000), (18000, 19000)],
    "abundance": 0,
}])

sim_cfg = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.65, seed=SEED,
)
gdna_cfg = GDNAConfig(
    abundance=20, frag_mean=350, frag_std=100,
    frag_min=100, frag_max=1000,
)

result = sc.build_oracle(
    n_fragments=N_FRAGMENTS, sim_config=sim_cfg,
    gdna_config=gdna_cfg, nrna_abundance=0,
)

# ── Run pipeline with default VBEM ──
config_vbem = PipelineConfig(
    em=EMConfig(seed=SEED, mode="vbem"),
    scan=BamScanConfig(sj_strand_tag="auto"),
)
pr_vbem = run_pipeline(result.bam_path, result.index, config=config_vbem)
bench_vbem = run_benchmark(result, pr_vbem, scenario_name="g20_n0_s65_vbem")

# ── Run pipeline with MAP-EM ──
config_map = PipelineConfig(
    em=EMConfig(seed=SEED, mode="map"),
    scan=BamScanConfig(sj_strand_tag="auto"),
)
pr_map = run_pipeline(result.bam_path, result.index, config=config_map)
bench_map = run_benchmark(result, pr_map, scenario_name="g20_n0_s65_map")

# ── Run pipeline with higher c_base ──
config_high_c = PipelineConfig(
    em=EMConfig(seed=SEED, mode="vbem"),
    scan=BamScanConfig(sj_strand_tag="auto"),
    calibration=CalibrationConfig(gdna_prior_c_base=50.0),
)
pr_high_c = run_pipeline(result.bam_path, result.index, config=config_high_c)
bench_high_c = run_benchmark(result, pr_high_c, scenario_name="g20_n0_s65_vbem_c50")

# ── Use stats from VBEM run ──
stats = pr_vbem.stats
SS = pr_vbem.calibration.strand_specificity if pr_vbem.calibration else 0.5

print("=" * 70)
print("DIAGNOSIS: g20_n0_s65 (gDNA=20, nRNA=0, SS=0.65)")
print("=" * 70)

# ── Ground truth ──
print(f"\n--- Ground Truth ---")
print(f"  N fragments simulated: {result.n_simulated}")
print(f"  gDNA fragments expected: {bench_vbem.n_gdna_expected}")
print(f"  nRNA fragments expected: {bench_vbem.n_nrna_expected}")
for ta in bench_vbem.transcripts:
    print(f"  {ta.t_id}: expected={ta.expected}")
total_rna_expected = bench_vbem.total_expected + bench_vbem.n_nrna_expected
print(f"  Total RNA expected: {total_rna_expected}")
print(f"  Total mRNA expected: {bench_vbem.total_expected}")

# ── Strand model ──
print(f"\n--- Strand Model ---")
print(f"  Strand specificity: {SS:.4f}")
print(f"  Strand specificity (reported): {SS:.4f}")

# ── Calibration ──
cal = pr_vbem.calibration
if cal is not None:
    print(f"\n--- Calibration ---")
    print(f"  lambda_gdna: {cal.lambda_gdna:.6f}")
    print(f"  region_e_gdna (sum): {cal.region_e_gdna.sum():.2f}")
    print(f"  region_n_total (sum): {cal.region_n_total.sum():.2f}")
    n_regions = len(cal.region_e_gdna)
    gamma_global = cal.region_e_gdna.sum() / max(cal.region_n_total.sum(), 1.0)
    print(f"  Global γ: {gamma_global:.4f}")
    print(f"  Number of regions: {n_regions}")
    print(f"  Per-region e_gdna: {cal.region_e_gdna}")
    print(f"  Per-region n_total: {cal.region_n_total}")
    if cal.gdna_fl_model is not None:
        print(f"  gDNA FL model mean: {cal.gdna_fl_model.mean:.1f}")
        print(f"  gDNA FL model total_weight: {cal.gdna_fl_model.total_weight:.1f}")
else:
    print("  CALIBRATION IS NONE")

# ── Comparison ──
print(f"\n--- VBEM Results (c_base=5.0) ---")
print(bench_vbem.summary())
total_rna_obs_vbem = bench_vbem.total_rna_observed
rna_rel_err_vbem = abs(total_rna_obs_vbem - total_rna_expected) / total_rna_expected
print(f"  Total RNA observed: {total_rna_obs_vbem:.0f}")
print(f"  Total RNA relative error: {rna_rel_err_vbem:.4f}")

print(f"\n--- MAP-EM Results (c_base=5.0) ---")
print(bench_map.summary())
total_rna_obs_map = bench_map.total_rna_observed
rna_rel_err_map = abs(total_rna_obs_map - total_rna_expected) / total_rna_expected
print(f"  Total RNA observed: {total_rna_obs_map:.0f}")
print(f"  Total RNA relative error: {rna_rel_err_map:.4f}")

print(f"\n--- VBEM Results (c_base=50.0) ---")
print(bench_high_c.summary())
total_rna_obs_hc = bench_high_c.total_rna_observed
rna_rel_err_hc = abs(total_rna_obs_hc - total_rna_expected) / total_rna_expected
print(f"  Total RNA observed: {total_rna_obs_hc:.0f}")
print(f"  Total RNA relative error: {rna_rel_err_hc:.4f}")

# ── EM breakdown: where do fragments go? ──
print(f"\n--- EM Component Breakdown ---")
est_v = pr_vbem.estimator
est_m = pr_map.estimator
print(f"  VBEM: gdna_em={est_v.gdna_em_count:.0f}, nrna_em={est_v.nrna_em_count:.0f}, "
      f"mRNA_total={bench_vbem.total_observed:.0f}")
print(f"  MAP:  gdna_em={est_m.gdna_em_count:.0f}, nrna_em={est_m.nrna_em_count:.0f}, "
      f"mRNA_total={bench_map.total_observed:.0f}")

# ── How much gDNA is being misclassified as RNA? ──
gdna_expected = bench_vbem.n_gdna_expected
gdna_pipeline_vbem = bench_vbem.n_gdna_pipeline
gdna_pipeline_map = bench_map.n_gdna_pipeline
print(f"\n--- gDNA Separation ---")
print(f"  gDNA expected: {gdna_expected}")
print(f"  gDNA pipeline (VBEM): {gdna_pipeline_vbem:.0f} (error: {abs(gdna_pipeline_vbem - gdna_expected):.0f})")
print(f"  gDNA pipeline (MAP):  {gdna_pipeline_map:.0f} (error: {abs(gdna_pipeline_map - gdna_expected):.0f})")
print(f"  gDNA leaking to RNA (VBEM): {max(0, gdna_expected - gdna_pipeline_vbem):.0f}")
print(f"  gDNA leaking to RNA (MAP):  {max(0, gdna_expected - gdna_pipeline_map):.0f}")

# ── nRNA siphon check ──
print(f"\n--- nRNA Siphon Check ---")
print(f"  nRNA expected: {bench_vbem.n_nrna_expected}")
print(f"  nRNA pipeline (VBEM): {bench_vbem.n_nrna_pipeline:.0f}")
print(f"  nRNA pipeline (MAP):  {bench_map.n_nrna_pipeline:.0f}")

# ── Intergenic ──
print(f"\n--- Intergenic ---")
print(f"  Intergenic (VBEM): {bench_vbem.n_intergenic}")
print(f"  Intergenic (MAP):  {bench_map.n_intergenic}")

# ── c_base sweep ──
print(f"\n--- c_base Sweep ---")
for c_val in [0.5, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0]:
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, mode="vbem"),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=CalibrationConfig(gdna_prior_c_base=c_val),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    b = run_benchmark(result, pr, scenario_name=f"c_base={c_val}")
    total_rna = b.total_rna_observed
    err = abs(total_rna - total_rna_expected) / total_rna_expected
    print(f"  c_base={c_val:6.1f}: total_rna={total_rna:.0f}, "
          f"gdna_pipeline={b.n_gdna_pipeline:.0f}, "
          f"nrna_pipeline={b.n_nrna_pipeline:.0f}, "
          f"RNA_rel_err={err:.4f}, "
          f"gdna_err={abs(b.n_gdna_pipeline - gdna_expected):.0f}")

# ── Also try MAP-EM with c_base sweep ──
print(f"\n--- c_base Sweep (MAP-EM) ---")
for c_val in [0.5, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0]:
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, mode="map"),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=CalibrationConfig(gdna_prior_c_base=c_val),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    b = run_benchmark(result, pr, scenario_name=f"map_c_base={c_val}")
    total_rna = b.total_rna_observed
    err = abs(total_rna - total_rna_expected) / total_rna_expected
    print(f"  c_base={c_val:6.1f}: total_rna={total_rna:.0f}, "
          f"gdna_pipeline={b.n_gdna_pipeline:.0f}, "
          f"nrna_pipeline={b.n_nrna_pipeline:.0f}, "
          f"RNA_rel_err={err:.4f}, "
          f"gdna_err={abs(b.n_gdna_pipeline - gdna_expected):.0f}")

sc.cleanup()
