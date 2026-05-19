"""Peek nrna_quant and quant schemas."""
import pandas as pd
nrna = pd.read_feather("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/nrna_quant.feather")
quant = pd.read_feather("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/quant.feather")
loci = pd.read_feather("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/loci.feather")
print("=== nrna_quant ===")
print(list(nrna.columns))
print(nrna.head(2).T)
print("\n=== quant ===")
print(list(quant.columns))
print(quant.head(2).T)
print("\n=== loci ===")
print(list(loci.columns))
print(loci.head(2).T)
