#!/usr/bin/env python3
"""Detailed GAPDH transcript breakdown across overhang alpha values.

Shows all GAPDH transcript TPMs for key alpha values to understand
the full effect of overhang penalty on isoform distribution.
"""
import sys
sys.path.insert(0, "src")

import math
import subprocess
import tempfile
import os
import pandas as pd
import pyarrow.feather as pf
from pathlib import Path

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"
ORACLE_DIR = f"{BASE}/rigel_oracle"

MANE_T_ID = "ENST00000229239.10"

# Load index to get GAPDH transcript IDs and lengths
from rigel.index import TranscriptIndex
print("Loading index...")
tx_index = TranscriptIndex.load(IDX)
tx_df = tx_index.t_df
gapdh_df = tx_df[(tx_df["g_name"] == "GAPDH") & (~tx_df["is_synthetic_nrna"])].copy()
print(f"  GAPDH mRNA transcripts: {len(gapdh_df)}")
gapdh_tids = set(gapdh_df["t_id"].tolist())

# Alpha values: key endpoints + default
ALPHA_VALUES = [0.01, 0.1, 0.3, 0.5, 1.0]


def load_quant(quant_dir: str) -> pd.DataFrame:
    tx_file = Path(quant_dir) / "quant.feather"
    return pf.read_feather(tx_file)


def run_quant(alpha: float, output_dir: str) -> None:
    cmd = [
        "conda", "run", "-n", "rigel",
        "rigel", "quant",
        "--bam", MM2_BAM, "--index", IDX,
        "-o", output_dir,
        "--overhang-alpha", str(alpha),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"quant failed:\n{result.stderr[-2000:]}")


# Load oracle
print("Loading oracle...")
oracle_all = load_quant(ORACLE_DIR)
tpm_col = next(c for c in oracle_all.columns if "tpm" in c.lower())
id_col = next(c for c in oracle_all.columns if c in ("t_id", "transcript_id", "name"))
oracle_tpm = dict(zip(oracle_all[id_col], oracle_all[tpm_col]))

# Filter to GAPDH
gapdh_oracle = gapdh_df[["t_id", "length"]].copy()
gapdh_oracle = gapdh_oracle.merge(
    oracle_all[[id_col, tpm_col]].rename(columns={id_col: "t_id", tpm_col: "oracle_tpm"}),
    on="t_id", how="left"
).fillna(0)
gapdh_oracle["is_mane"] = gapdh_df["is_mane"].values
gapdh_oracle = gapdh_oracle.sort_values("oracle_tpm", ascending=False)

print(f"\nOracle GAPDH transcript TPMs (oracle run):")
for _, r in gapdh_oracle.iterrows():
    marker = " <<< MANE" if r.get("is_mane") else ""
    print(f"  {r['t_id']:30s} len={r['length']:5.0f}  TPM={r['oracle_tpm']:10.1f}{marker}")

oracle_mane = oracle_tpm.get(MANE_T_ID, float("nan"))
print(f"\nOracle MANE TPM = {oracle_mane:.1f}")
print(f"\n{'='*80}")
print("SWEEP: Minimap2 BAM — GAPDH TPMs vs overhang alpha")
print(f"{'='*80}")

with tempfile.TemporaryDirectory(prefix="rigel_oh_") as tmpdir:
    for alpha in ALPHA_VALUES:
        oh_pen = math.log(alpha) if alpha > 0 else float("-inf")
        outdir = os.path.join(tmpdir, f"a{alpha:.3f}".replace(".", "p"))
        print(f"\n--- alpha={alpha:.2f}  oh_pen={oh_pen:+.4f} ---")
        run_quant(alpha, outdir)
        quant = load_quant(outdir)
        mm2_tpm = dict(zip(quant[id_col], quant[tpm_col]))

        total_gapdh = sum(mm2_tpm.get(tid, 0) for tid in gapdh_tids)
        for _, r in gapdh_oracle.iterrows():
            tid = r["t_id"]
            mm2 = mm2_tpm.get(tid, 0)
            ora = r['oracle_tpm']
            ratio = mm2 / ora if ora > 0 else float("nan")
            ratio_s = f"{ratio:.3f}" if not math.isnan(ratio) else " N/A "
            marker = " <<< MANE" if r.get("is_mane") else ""
            print(f"  {tid:30s} oracle={ora:8.1f}  mm2={mm2:8.1f}  ratio={ratio_s}{marker}")
        print(f"  {'TOTAL GAPDH':30s} oracle={sum(gapdh_oracle['oracle_tpm']):8.1f}  mm2={total_gapdh:8.1f}")

print("\nDONE")
