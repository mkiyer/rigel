"""Compute actual log-lik comparison with real calibration FL model."""
import numpy as np

# With the Dirichlet prior FL models:
gdna_fl_200 = -4.8574
rna_fl_200 = -5.0166

# Scoring equations (unspliced fragment, flen=200):
# RNA: log_strand + rna_fl + oh*pen + nm*pen
#   log_p_sense=-0.3185 (r1_antisense=False, SS=0.727)
#   log_p_antisense=-1.2993
# gDNA: gdna_fl + gdna_log_sp + LOG_HALF + nm*pen
#   gdna_log_sp=0.0 (unspliced)
#   LOG_HALF=-0.6931

log_half = -0.6931

rna_sense = -0.3185 + rna_fl_200
rna_anti = -1.2993 + rna_fl_200
gdna = gdna_fl_200 + 0.0 + log_half

print(f"Pre-bias-correction log-likelihoods (flen=200):")
print(f"  RNA (sense strand): {rna_sense:.4f}")
print(f"  RNA (anti strand):  {rna_anti:.4f}")
print(f"  gDNA:               {gdna:.4f}")
print()
print(f"  gDNA - RNA(sense): {gdna - rna_sense:.4f}")
print(f"  gDNA - RNA(anti):  {gdna - rna_anti:.4f}")

# Now add bias correction
# RNA t2 (length 8000, flen 200): eff_len = 8000-200+1 = 7801
# gDNA (span ~20000, flen 200): eff_len = 20000-200+1 = 19801
rna_bias = -np.log(7801)
gdna_bias = -np.log(19801)

rna_sense_bc = rna_sense + rna_bias
rna_anti_bc = rna_anti + rna_bias
gdna_bc = gdna + gdna_bias

print(f"\nAfter bias correction:")
print(f"  RNA (sense): {rna_sense_bc:.4f}")
print(f"  RNA (anti):  {rna_anti_bc:.4f}")
print(f"  gDNA:        {gdna_bc:.4f}")
print(f"  gDNA - RNA(sense): {gdna_bc - rna_sense_bc:.4f}")
print(f"  gDNA - RNA(anti):  {gdna_bc - rna_anti_bc:.4f}")

print(f"\nProbability ratio (exp of diff):")
print(f"  gDNA/RNA(sense): {np.exp(gdna_bc - rna_sense_bc):.4f}")
print(f"  gDNA/RNA(anti):  {np.exp(gdna_bc - rna_anti_bc):.4f}")
