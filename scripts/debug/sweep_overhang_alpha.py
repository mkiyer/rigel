#!/usr/bin/env python3
"""Sweep overhang alpha and measure GAPDH MANE quantification accuracy.

For each alpha value, runs rigel quant on the GAPDH minimap2 BAM and
reports TPM for MANE (ENST00000229239.10) vs oracle.

Usage:
    conda activate rigel
    python scripts/debug/sweep_overhang_alpha.py
"""
import sys
sys.path.insert(0, "src")

import math
import subprocess
import tempfile
import os
import pandas as pd
from pathlib import Path

BASE = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
MM2_BAM = f"{BASE}/sim_minimap2.bam"
ORACLE_DIR = f"{BASE}/rigel_oracle"

MANE_T_ID = "ENST00000229239.10"
# t_indices of GAPDH mRNA transcripts (from previous analysis)
GAPDH_T_IDS = [
    "ENST00000229239.10",  # MANE
    "ENST00000396856.5",
    "ENST00000396859.5",
    "ENST00000396861.4",
    "ENST00000466525.1",
    "ENST00000619601.1",
]

# Alpha values to sweep (0.01 = default, 1.0 = no penalty)
ALPHA_VALUES = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]


def load_quant_tpm(quant_dir: str) -> dict:
    """Returns {t_id: tpm} from rigel quant output (feather format)."""
    import pyarrow.feather as pf
    tx_file = Path(quant_dir) / "quant.feather"
    if not tx_file.exists():
        raise FileNotFoundError(tx_file)
    df = pf.read_feather(tx_file)
    tpm_col = next(c for c in df.columns if "tpm" in c.lower())
    id_col = next(c for c in df.columns if c in ("t_id", "transcript_id", "name"))
    return dict(zip(df[id_col], df[tpm_col]))


def run_quant(alpha: float, output_dir: str) -> None:
    """Run rigel quant with given overhang alpha."""
    cmd = [
        "conda", "run", "-n", "rigel",
        "rigel", "quant",
        "--bam", MM2_BAM,
        "--index", IDX,
        "-o", output_dir,
        "--overhang-alpha", str(alpha),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR (alpha={alpha}):\n{result.stderr[-2000:]}")
        raise RuntimeError(f"quant failed for alpha={alpha}")


# Load oracle TPMs
print("Loading oracle quantification...")
try:
    oracle_tpm = load_quant_tpm(ORACLE_DIR)
    oracle_mane = oracle_tpm.get(MANE_T_ID, float("nan"))
    print(f"  Oracle MANE TPM: {oracle_mane:.1f}")
except FileNotFoundError as e:
    print(f"  WARNING: oracle not found ({e}), will skip oracle comparison")
    oracle_tpm = {}
    oracle_mane = float("nan")

# Header
print(f"\n{'alpha':>6}  {'oh_pen':>8}  {'MANE TPM':>10}  {'ratio/oracle':>13}  {'ENST466525 TPM':>15}")
print("-" * 60)

results = []
with tempfile.TemporaryDirectory(prefix="rigel_oh_sweep_") as tmpdir:
    for alpha in ALPHA_VALUES:
        oh_pen = math.log(alpha) if alpha > 0 else float("-inf")
        outdir = os.path.join(tmpdir, f"alpha_{alpha:.3f}".replace(".", "p"))
        try:
            run_quant(alpha, outdir)
            tpm = load_quant_tpm(outdir)
            mane_tpm = tpm.get(MANE_T_ID, float("nan"))
            comp_tpm = tpm.get("ENST00000466525.1", float("nan"))
            ratio = mane_tpm / oracle_mane if oracle_mane > 0 else float("nan")
            print(f"  {alpha:>6.3f}  {oh_pen:>8.4f}  {mane_tpm:>10.2f}  "
                  f"{'N/A' if math.isnan(ratio) else f'{ratio:.3f}':>13}  {comp_tpm:>15.2f}")
            results.append({
                "alpha": alpha, "oh_pen": oh_pen,
                "mane_tpm": mane_tpm,
                "oracle_mane_tpm": oracle_mane,
                "ratio": ratio,
                "comp_enst466525_tpm": comp_tpm,
            })
        except Exception as ex:
            print(f"  {alpha:>6.3f}  FAILED: {ex}")
            continue

# Summary
print(f"\nSummary: oracle MANE TPM = {oracle_mane:.1f}")
if results:
    best = max(results, key=lambda r: r["ratio"] if not math.isnan(r["ratio"]) else -1)
    print(f"Best alpha for MANE recovery: {best['alpha']:.3f} "
          f"(ratio={best['ratio']:.3f}, TPM={best['mane_tpm']:.1f})")

print("\nDONE")
