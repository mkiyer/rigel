#!/usr/bin/env python3
"""Investigate VBEM zero-forcing on specific loci.

Load the VBEM and MAP quant.feather files, identify the loci where VBEM
catastrophically fails, and analyze what's happening.
"""

import pandas as pd
import numpy as np
from pathlib import Path

BENCHMARK_DIR = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate")
CONDITION = "gdna_high_ss_0.90_nrna_none"
TRUTH_FILE = BENCHMARK_DIR / "truth_abundances_nrna_none.tsv"

def load_quant(tool_dir):
    """Load quant.feather and return DataFrame."""
    qf = tool_dir / "quant.feather"
    if qf.exists():
        return pd.read_feather(qf)
    return None

def load_nrna_quant(tool_dir):
    """Load nrna_quant.feather."""
    nf = tool_dir / "nrna_quant.feather"
    if nf.exists():
        return pd.read_feather(nf)
    return None

def load_loci(tool_dir):
    """Load loci.feather."""
    lf = tool_dir / "loci.feather"
    if lf.exists():
        return pd.read_feather(lf)
    return None

def main():
    # Load truth
    truth = pd.read_csv(TRUTH_FILE, sep="\t")
    print(f"Truth: {len(truth)} transcripts")
    print(f"Truth columns: {list(truth.columns)}")

    # Load VBEM and MAP results
    vbem_dir = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem"
    map_dir = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map"

    vbem_q = load_quant(vbem_dir)
    map_q = load_quant(map_dir)
    vbem_nrna = load_nrna_quant(vbem_dir)
    map_nrna = load_nrna_quant(map_dir)
    vbem_loci = load_loci(vbem_dir)
    map_loci = load_loci(map_dir)

    print(f"\nVBEM quant: {len(vbem_q)} transcripts")
    print(f"VBEM quant columns: {list(vbem_q.columns)}")
    print(f"\nMAP quant: {len(map_q)} transcripts")
    print(f"MAP quant columns: {list(map_q.columns)}")

    if vbem_nrna is not None:
        print(f"\nVBEM nRNA: {len(vbem_nrna)} entries")
        print(f"VBEM nRNA columns: {list(vbem_nrna.columns)}")

    if vbem_loci is not None:
        print(f"\nVBEM loci: {len(vbem_loci)} loci")
        print(f"VBEM loci columns: {list(vbem_loci.columns)}")

    # Merge VBEM and MAP quant with truth
    # Identify the top zero-forced transcripts
    zero_forced = [
        "ENST00000323345.11", "ENST00000222247.10", "ENST00000301821.11",
        "ENST00000396466.5", "ENST00000651669.1", "ENST00000009589.8",
        "ENST00000648437.1", "ENST00000329251.5", "ENST00000690968.2",
        "ENST00000479563.5", "ENST00000290902.10",
    ]

    # Transcripts where VBEM massively over-estimates
    over_estimated = [
        "ENST00000228306.8", "ENST00000468019.5", "ENST00000697730.1",
        "ENST00000344700.8", "ENST00000597648.5", "ENST00000424576.6",
        "ENST00000526409.5", "ENST00000406022.6", "ENST00000326092.8",
    ]

    target_txs = zero_forced + over_estimated

    print(f"\n{'='*80}")
    print("DETAILED ANALYSIS OF ZERO-FORCED AND OVER-ESTIMATED TRANSCRIPTS")
    print(f"{'='*80}")

    for tx_id in target_txs:
        vbem_row = vbem_q[vbem_q["transcript_id"] == tx_id]
        map_row = map_q[map_q["transcript_id"] == tx_id]
        truth_row = truth[truth["transcript_id"] == tx_id]

        if len(vbem_row) == 0 and len(map_row) == 0:
            print(f"\n{tx_id}: NOT FOUND in either result")
            continue

        print(f"\n--- {tx_id} ---")
        if len(truth_row) > 0:
            for col in truth_row.columns:
                print(f"  Truth.{col}: {truth_row[col].values[0]}")

        if len(vbem_row) > 0:
            print(f"  VBEM:")
            for col in vbem_row.columns:
                val = vbem_row[col].values[0]
                print(f"    {col}: {val}")

        if len(map_row) > 0:
            print(f"  MAP:")
            for col in map_row.columns:
                val = map_row[col].values[0]
                print(f"    {col}: {val}")

    # Now look at loci: which loci do these transcripts belong to?
    if vbem_loci is not None and "locus_id" in vbem_q.columns:
        print(f"\n{'='*80}")
        print("LOCI ANALYSIS FOR AFFECTED TRANSCRIPTS")
        print(f"{'='*80}")

        for tx_id in target_txs[:5]:  # just first 5
            vbem_row = vbem_q[vbem_q["transcript_id"] == tx_id]
            if len(vbem_row) > 0 and "locus_id" in vbem_row.columns:
                locus_id = vbem_row["locus_id"].values[0]
                locus_row = vbem_loci[vbem_loci["locus_id"] == locus_id]
                if len(locus_row) > 0:
                    print(f"\n--- Locus {locus_id} (contains {tx_id}) ---")
                    for col in locus_row.columns:
                        print(f"  {col}: {locus_row[col].values[0]}")

                    # All transcripts in this locus
                    locus_txs = vbem_q[vbem_q["locus_id"] == locus_id]
                    map_locus_txs = map_q[map_q.get("locus_id", pd.Series()) == locus_id] if "locus_id" in map_q.columns else pd.DataFrame()
                    print(f"\n  Transcripts in locus (VBEM): {len(locus_txs)}")
                    for _, r in locus_txs.iterrows():
                        t = truth[truth["transcript_id"] == r["transcript_id"]]
                        t_tpm = t["mrna_abundance"].values[0] if len(t) > 0 else "?"
                        print(f"    {r['transcript_id']}: vbem_tpm={r.get('tpm', r.get('est_count', '?'))}, truth={t_tpm}")

    # Analyze overall loci stats
    if vbem_loci is not None:
        print(f"\n{'='*80}")
        print("LOCI SUMMARY STATISTICS")
        print(f"{'='*80}")
        print(vbem_loci.describe().to_string())

    # Check quant columns for fragment counts
    print(f"\n{'='*80}")
    print("VBEM QUANT SAMPLE (first 5 rows)")
    print(f"{'='*80}")
    print(vbem_q.head().to_string())

    print(f"\n{'='*80}")
    print("VBEM nRNA SAMPLE (first 5 rows)")
    print(f"{'='*80}")
    if vbem_nrna is not None:
        print(vbem_nrna.head().to_string())

    # Where is the mass going?
    # Compare total est_count or tpm 
    print(f"\n{'='*80}")
    print("TOTAL MASS COMPARISON")
    print(f"{'='*80}")
    for col in vbem_q.columns:
        if vbem_q[col].dtype in [np.float64, np.float32, np.int64]:
            v_sum = vbem_q[col].sum()
            m_sum = map_q[col].sum() if col in map_q.columns else "N/A"
            if v_sum != 0 or m_sum != "N/A":
                print(f"  {col}: VBEM={v_sum:.2f}, MAP={m_sum if isinstance(m_sum, str) else f'{m_sum:.2f}'}")


if __name__ == "__main__":
    main()
