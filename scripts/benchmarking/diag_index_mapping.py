#!/usr/bin/env python3
"""Inspect the hulkrna index to check transcript-to-gene mapping."""
import sys
sys.path.insert(0, "src")

from hulkrna.index import HulkIndex, transcripts_to_dataframe

GTF = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/"
    "PVT1_MYC/region.gtf"
)

from hulkrna.transcript import Transcript

# Build transcripts and create index manually
transcripts = Transcript.read_gtf(GTF)
g_id_to_index = {}
for ti, t in enumerate(transcripts):
    t.t_index = ti
    if t.g_id not in g_id_to_index:
        g_id_to_index[t.g_id] = len(g_id_to_index)
    t.g_index = g_id_to_index[t.g_id]

import pandas as pd
from hulkrna.index import transcripts_to_dataframe

t_df = transcripts_to_dataframe(transcripts)
g_df = (
    t_df
    .groupby("g_index", sort=True)
    .agg(
        ref=("ref", "first"),
        start=("start", "min"),
        end=("end", "max"),
        strand=("strand", "first"),
        g_id=("g_id", "first"),
        g_name=("g_name", "first"),
    )
    .reset_index()
)

class FakeIndex:
    pass
index = FakeIndex()
index.t_df = t_df
index.g_df = g_df
index.t_to_g_arr = t_df["g_index"].values

# Find ENST00000624314.1
t_df = index.t_df
g_df = index.g_df

print("=== Gene table ===")
for i, row in g_df.iterrows():
    print(f"  g_idx={i:3d}  {str(row.get('g_id', 'N/A')):30s}  "
          f"{str(row.get('g_name', '')):15s}  strand={row.get('strand', '?')}")

print()
print("=== ENST00000624314.1 lookup ===")
mask = t_df["t_id"] == "ENST00000624314.1"
if mask.any():
    row = t_df[mask].iloc[0]
    idx = t_df.index[mask][0]
    print(f"  t_idx={idx}")
    print(f"  t_id={row.get('t_id')}")
    print(f"  g_id={row.get('g_id')}")
    print(f"  g_name={row.get('g_name')}")
    print(f"  strand={row.get('strand')}")
    print(f"  start={row.get('start')}")
    print(f"  end={row.get('end')}")

    # Check t_to_g mapping
    g_idx = int(index.t_to_g_arr[idx])
    print(f"  t_to_g_arr[{idx}] = {g_idx}")
    g_row = g_df.iloc[g_idx]
    print(f"  → gene: {g_row.get('g_id')} ({g_row.get('g_name')})")
else:
    print("  NOT FOUND in index!")

print()
print("=== t_to_g for all transcripts of ENSG00000280055.2 (LINC02912) ===")
for i, row in t_df.iterrows():
    if row.get("g_id") == "ENSG00000280055.2":
        g_idx = int(index.t_to_g_arr[i])
        g_row = g_df.iloc[g_idx]
        print(f"  t_idx={i} {row['t_id']} → g_idx={g_idx} ({g_row.get('g_id')})")

print()
print("=== How many genes total? ===")
print(f"  {len(g_df)} genes, {len(t_df)} transcripts")
