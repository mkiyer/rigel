"""Diagnose VCaP v4.3 regression: per-class rho scales, kappa sensitivity."""
import json
import numpy as np

SUMMARY = "/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm/summary.json"

d = json.load(open(SUMMARY))
re = d["calibration"]["regional_exposure"]
pc = re["per_class"]

print(f"kappa_global = {re['kappa_global']}")
print(f"rho_global   = {re['rho_global']:.3e}")
print(f"rho_ref (Q95)= {re['rho_ref']:.3e}")
print(f"n_at_floor   = {re['n_at_floor']:,}")
print(f"observed_log_spread = {re['observed_log_spread']:.2f}")
print(f"null_log_spread     = {re['null_log_spread']:.2f}")
print()

print(
    f"{'class':<14}{'n_reg':>10}{'Y_tot':>14}{'E_tot':>14}"
    f"{'rho_avg':>12}{'rho_q50':>12}{'rho_q95':>12}{'q95/ref':>10}"
)
for c, p in pc.items():
    Y, E = p["evidence_count"], p["opportunity"]
    rho_avg = Y / E if E > 0 else 0
    print(
        f"{c:<14}{int(p['n_regions']):>10}{Y:>14.2e}{E:>14.2e}"
        f"{rho_avg:>12.3e}{p['rho_q50']:>12.3e}{p['rho_q95']:>12.3e}"
        f"{p['rho_q95'] / re['rho_ref']:>10.2f}"
    )

print()
print("=== Empty-region (Y=0) post-shrinkage rho and A under current kappa ===")
kappa = re["kappa_global"]
rg = re["rho_global"]
rref = re["rho_ref"]
for E in [1000, 5000, 10000, 50000, 1_300_000_000 / 33120]:
    rho_hat_empty = kappa * rg / (E + kappa)
    A = max(rho_hat_empty / rref, np.exp(-9.21))
    print(f"  E={E:>12.0f}  rho_hat(Y=0)={rho_hat_empty:.3e}  A={A:.3e}")

print()
print("=== Counterfactual: kappa = mean(E) per class ===")
for cls, mean_E in [
    ("INTRON", 1.5e9 / 296039),
    ("INTERGENIC", 1.3e9 / 33120),
    ("EXON", 180e6 / 318828),
]:
    k_alt = mean_E
    rho_hat_empty = k_alt * rg / (mean_E + k_alt)
    A = rho_hat_empty / rref
    print(
        f"  {cls:<12} mean_E={mean_E:>10.1f}  kappa={k_alt:>10.1f}  "
        f"rho_hat(Y=0)={rho_hat_empty:.3e}  A={A:.3f}"
    )

print()
print("=== Median A by component (from feathers) ===")
import pandas as pd
for fn, col in [
    ("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm/nrna_quant.feather", "em_exposure_weight"),
    ("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm/quant.feather", "em_exposure_weight"),
    ("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm/loci.feather", "gdna_em_exposure_weight"),
]:
    df = pd.read_feather(fn)
    print(f"  {fn.rsplit('/',1)[1]:<22} median={df[col].median():.3f}  mean={df[col].mean():.3f}")
