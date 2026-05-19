"""Inspect locus_stats.feather schema."""
import pandas as pd
df = pd.read_feather("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/locus_stats.feather")
print("rows:", len(df))
print("columns:", list(df.columns))
print(df.head(3).T)
