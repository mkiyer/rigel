"""Diagnose why SS=0.60 + FL-identical scenarios collapse with Gaussian density."""
import sys
sys.path.insert(0, "src")
sys.path.insert(0, "scripts/debug")

import numpy as np
import pandas as pd
from rigel.calibration import (
    _e_step, _m_step, _seed_initial_partition,
    _compute_density_llr_gaussian, _compute_strand_llr_binomial,
    compute_region_stats, compute_sense_fraction, compute_log_density,
    calibrate_gdna,
    _EPS,
)
from calibration_stress_test import ScenarioConfig, generate_scenario

# Run the actual failing scenario
cfg = ScenarioConfig(
    n_rna_regions=600,
    n_gdna_regions=400,
    n_nrna_regions=150,
    strand_specificity=0.60,
    rna_fl_mean=250.0,
    rna_fl_std=50.0,
    gdna_fl_mean=250.0,  # FL identical
    gdna_fl_std=50.0,
    gdna_count_mean=60.0,  # med_gdna
    nrna_count_mean=80.0,
    nrna_splice_rate=0.05,
    seed=42,
    label="ss0.60_med_gdna_fl_identical_mod_nrna",
)

data = generate_scenario(cfg)
region_counts = data["region_counts"]
region_df = data["region_df"]
fl_table = data["fl_table"]
truth = data["truth_labels"]  # 0=RNA, 1=gDNA, 2=nRNA
is_gdna = truth == 1

# Compute stats
stats = compute_region_stats(region_counts, region_df)
eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)
sense_frac = compute_sense_fraction(stats)
log_d, epsilon = compute_log_density(stats, eligible)
n_regions = len(stats["n_total"])

print(f"n_regions={n_regions}, n_eligible={eligible.sum()}")
print(f"Truth: {(truth==1).sum()} gDNA, {(truth==0).sum()} RNA, {(truth==2).sum()} nRNA")
print(f"n_spliced>0: {(stats['n_spliced']>0).sum()}")

has_splice = stats["n_spliced"] > 0
soft = eligible & ~has_splice
print(f"soft (eligible & unspliced): {soft.sum()}")
print(f"  of which gDNA: {(soft & is_gdna).sum()}")
print(f"  of which nRNA: {(soft & (truth==2)).sum()}")
print(f"  of which RNA: {(soft & (truth==0)).sum()}")
print()

# Log density by component
print("Log density by component:")
for label, mask in [("gDNA", truth==1), ("RNA", truth==0), ("nRNA", truth==2)]:
    ld = log_d[mask & eligible]
    print(f"  {label}: mean={ld.mean():.3f} std={ld.std():.3f} range=[{ld.min():.3f}, {ld.max():.3f}]")
print()

# Seed partition
gamma, pi_init, seed_diag = _seed_initial_partition(
    stats, log_d, sense_frac, 0.60, eligible
)
print(f"Seed: pi={pi_init:.4f}, n_expressed_seed={seed_diag['n_expressed_seed']}, "
      f"n_gdna_seed={seed_diag['n_gdna_seed']}")

# Init Gaussian params from seed
w_g = gamma[eligible]
w_e = 1.0 - gamma[eligible]
ld = log_d[eligible]
mu_overall = float(np.mean(ld))
var_overall = max(float(np.var(ld)), _EPS)

w_g_sum = max(float(w_g.sum()), _EPS)
mu_g = float(np.sum(w_g * ld) / w_g_sum)
var_g = max(float(np.sum(w_g * (ld - mu_g) ** 2) / w_g_sum), _EPS)
w_e_sum = max(float(w_e.sum()), _EPS)
mu_r = float(np.sum(w_e * ld) / w_e_sum)
var_r = max(float(np.sum(w_e * (ld - mu_r) ** 2) / w_e_sum), _EPS)

n_soft = int(soft.sum())
n_eligible_count = int(eligible.sum())
pi_soft = float(np.clip(pi_init * n_eligible_count / n_soft, _EPS, 1.0 - _EPS)) if n_soft > 0 else pi_init

print(f"pi_soft={pi_soft:.4f}")
print(f"Overall: mu={mu_overall:.3f}, var={var_overall:.3f}")
print(f"Init: mu_g={mu_g:.3f} var_g={var_g:.4f} | mu_r={mu_r:.3f} var_r={var_r:.4f}")
print()

# EM trace
pi = pi_soft
for it in range(20):
    gamma = _e_step(
        stats, pi, log_d, eligible, 0.60,
        mu_g, var_g, mu_r, var_r,
    )
    llr_d = _compute_density_llr_gaussian(log_d, mu_g, var_g, mu_r, var_r, soft, n_regions)
    llr_s = _compute_strand_llr_binomial(stats, 0.60, n_regions)

    gdna_gamma = gamma[is_gdna].mean()
    rna_gamma = gamma[truth==0].mean()
    nrna_gamma = gamma[truth==2].mean()

    new_pi, _, _, mu_g, var_g, mu_r, var_r = _m_step(stats, gamma, log_d, eligible)
    if n_soft > 0:
        pi = float(np.clip(gamma[soft].sum() / n_soft, _EPS, 1.0 - _EPS))
    else:
        pi = new_pi

    d_soft = llr_d[soft]
    s_soft = llr_s[soft]
    print(f"It {it:2d}: pi={pi:.4f} mu_g={mu_g:.3f} var_g={var_g:.4f} mu_r={mu_r:.3f} var_r={var_r:.4f}"
          f" | gDNA_g={gdna_gamma:.4f} RNA_g={rna_gamma:.4f} nRNA_g={nrna_gamma:.4f}"
          f" | density[{d_soft.min():.1f},{d_soft.max():.1f}] strand[{s_soft.min():.1f},{s_soft.max():.1f}]")

print("\n=== Now running actual calibrate_gdna ===")
result = calibrate_gdna(
    region_counts, fl_table, region_df, 0.60, diagnostics=True,
)
print(f"mixing_proportion={result.mixing_proportion:.6f}")
print(f"n_iterations={result.n_iterations}")
print(f"gDNA posterior mean={result.region_posteriors[is_gdna].mean():.4f}")

# Now verify: run EM WITH FL models to confirm FL is the cause
from rigel.calibration import build_gdna_fl_model, _compute_fl_llr
import math

gamma2, pi_init2, seed_diag2 = _seed_initial_partition(
    stats, log_d, sense_frac, 0.60, eligible
)
w_g2 = gamma2[eligible]
w_e2 = 1.0 - gamma2[eligible]
ld2 = log_d[eligible]
mu_o2 = float(np.mean(ld2))
var_o2 = max(float(np.var(ld2)), _EPS)
w_g_s2 = max(float(w_g2.sum()), _EPS)
mu_g2 = float(np.sum(w_g2 * ld2) / w_g_s2)
var_g2 = max(float(np.sum(w_g2 * (ld2 - mu_g2) ** 2) / w_g_s2), _EPS)
w_e_s2 = max(float(w_e2.sum()), _EPS)
mu_r2 = float(np.sum(w_e2 * ld2) / w_e_s2)
var_r2 = max(float(np.sum(w_e2 * (ld2 - mu_r2) ** 2) / w_e_s2), _EPS)

fl_rids = fl_table["region_id"].values.astype(np.intp)
fl_fls = fl_table["frag_len"].values.astype(np.intp)
gdna_fl = build_gdna_fl_model(fl_rids, fl_fls, gamma2)
rna_fl = build_gdna_fl_model(fl_rids, fl_fls, 1.0 - gamma2)

pi2 = float(np.clip(pi_init2 * n_eligible_count / n_soft, _EPS, 1.0 - _EPS)) if n_soft > 0 else pi_init2

print(f"\n=== EM with FL models (should reproduce collapse) ===")
print(f"Init FL weights: gDNA_fl_total={gdna_fl._total_weight:.0f}, RNA_fl_total={rna_fl._total_weight:.0f}")
print(f"FL normalization bias C = log((T_r+M)/(T_g+M)) = {math.log((rna_fl._total_weight + 1001) / (gdna_fl._total_weight + 1001)):.4f}")

for it in range(20):
    gamma2 = _e_step(
        stats, pi2, log_d, eligible, 0.60,
        mu_g2, var_g2, mu_r2, var_r2,
        fl_region_ids=fl_rids, fl_frag_lens=fl_fls,
        gdna_fl_model=gdna_fl, rna_fl_model=rna_fl,
    )

    llr_d2 = _compute_density_llr_gaussian(log_d, mu_g2, var_g2, mu_r2, var_r2, soft, n_regions)
    llr_s2 = _compute_strand_llr_binomial(stats, 0.60, n_regions)
    llr_f2 = _compute_fl_llr(fl_rids, fl_fls, gdna_fl, rna_fl, n_regions)

    gdna_g2 = gamma2[is_gdna].mean()

    new_pi2, _, _, mu_g2, var_g2, mu_r2, var_r2 = _m_step(stats, gamma2, log_d, eligible)
    gdna_fl = build_gdna_fl_model(fl_rids, fl_fls, gamma2)
    rna_fl = build_gdna_fl_model(fl_rids, fl_fls, 1.0 - gamma2)
    if n_soft > 0:
        pi2 = float(np.clip(gamma2[soft].sum() / n_soft, _EPS, 1.0 - _EPS))

    fl_gdna = llr_f2[is_gdna & soft]
    fl_rna = llr_f2[(truth != 1) & soft]
    fw_g = gdna_fl._total_weight
    fw_r = rna_fl._total_weight
    C = math.log((fw_r + 1001) / (fw_g + 1001)) if fw_g > 0 and fw_r > 0 else 0

    print(f"It {it:2d}: pi={pi2:.4f} gDNA_g={gdna_g2:.4f}"
          f" | FL_gDNA_mean={fl_gdna.mean():.1f} FL_rna_mean={fl_rna.mean():.1f}"
          f" | FL_T_g={fw_g:.0f} FL_T_r={fw_r:.0f} FL_bias_C={C:.3f}"
          f" | dens[{llr_d2[soft].min():.1f},{llr_d2[soft].max():.1f}]"
          f" | strand[{llr_s2[soft].min():.1f},{llr_s2[soft].max():.1f}]")
