#!/usr/bin/env python3
"""Trace extreme FL outliers and compare gap correction approaches.

Focuses on SPLICED_ANNOT fragments with unusual FL values in Minimap2.
Also computes what gap correction would look like under strict containment
vs overlap by re-querying the SJ index.

Usage:
    conda activate rigel
    python scripts/debug/trace_fl_outliers.py
"""
import sys
import logging
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, "src")

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType

# --- Paths ---
BENCH = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3")
INDEX = BENCH / "rigel_index"
COND = BENCH / "gdna_none_ss_0.95_nrna_none"
ORACLE_BAM = COND / "align_oracle" / "reads_namesort.bam"
MM2_BAM = COND / "align_minimap2" / "reads_namesort.bam"

# --- Load index ---
logger.info("Loading index...")
idx = TranscriptIndex.load(str(INDEX))
tx_df = idx.t_df
is_nrna = tx_df["is_synthetic_nrna"].values
tx_names = tx_df["t_id"].values

def analyze_fl_details(buf, is_nrna, label, fl_threshold=400):
    """Deep analysis of fragments with unusual FL values."""
    print(f"\n{'='*70}")
    print(f"FL Outlier Analysis: {label}")
    print(f"{'='*70}")

    # Collect detailed info for fragments with FL anomalies
    outliers = []  # (frag_idx, splice_type, nh, t_inds, fls, footprint, read_len, nm)
    fl_hist_spliced_annot = defaultdict(int)  # FL value → count for SPLICED_ANNOT
    fl_hist_unspliced = defaultdict(int)  # FL value → count for UNSPLICED
    
    # Count how many fragments have FL disagreement within mRNA candidates
    mrna_fl_disagree = 0
    mrna_fl_agree = 0
    
    # Track FL per-transcript: what fraction of FL estimates come from gap correction?
    gap_correction_needed = 0  # fragments where mRNA candidates have different FLs
    total_spliced_annot = 0
    total_unspliced = 0
    
    for chunk in buf.iter_chunks():
        n = chunk.size
        for i in range(n):
            s = chunk.t_offsets[i]
            e = chunk.t_offsets[i + 1]
            t_inds = chunk.t_indices[s:e]
            fls = chunk.frag_lengths[s:e]
            st = chunk.splice_type[i]
            nh = chunk.num_hits[i]
            fp = chunk.genomic_footprint[i]
            rl = chunk.read_length[i]
            nm_val = chunk.nm[i]
            
            if len(t_inds) == 0:
                continue
            
            nrna_mask = is_nrna[t_inds]
            mrna_mask = ~nrna_mask
            mrna_fls = fls[mrna_mask]
            valid_mrna = mrna_fls[mrna_fls > 0]
            
            if st == SpliceType.SPLICED_ANNOT.value:
                total_spliced_annot += 1
                if len(valid_mrna) > 0:
                    fl_hist_spliced_annot[int(valid_mrna[0])] += 1
                if len(valid_mrna) > 1:
                    if len(set(valid_mrna.tolist())) > 1:
                        mrna_fl_disagree += 1
                    else:
                        mrna_fl_agree += 1
                        
            elif st == SpliceType.UNSPLICED.value:
                total_unspliced += 1
                if len(valid_mrna) > 0:
                    fl_hist_unspliced[int(valid_mrna[0])] += 1
            
            # Collect outliers (any FL > threshold or FL > 2 * read_length)
            valid_fls = fls[fls > 0]
            if len(valid_fls) > 0 and (valid_fls.max() > fl_threshold or 
                                        valid_fls.max() > 2 * rl):
                outlier = {
                    "splice_type": SpliceType(st).name,
                    "nh": int(nh),
                    "footprint": int(fp),
                    "read_length": int(rl),
                    "nm": int(nm_val),
                    "n_candidates": len(t_inds),
                    "fls": fls.tolist(),
                    "t_inds": t_inds.tolist(),
                    "is_nrna": nrna_mask.tolist(),
                    "max_fl": int(valid_fls.max()),
                    "min_fl": int(valid_fls.min()) if len(valid_fls) > 0 else 0,
                }
                outliers.append(outlier)
    
    print(f"\nTotal SPLICED_ANNOT: {total_spliced_annot}")
    print(f"Total UNSPLICED: {total_unspliced}")
    print(f"mRNA FL disagreement (SPLICED_ANNOT, >1 mRNA): agree={mrna_fl_agree}, "
          f"disagree={mrna_fl_disagree}")
    
    # Show FL distribution percentiles for SPLICED_ANNOT
    if fl_hist_spliced_annot:
        all_fls = []
        for fl, count in fl_hist_spliced_annot.items():
            all_fls.extend([fl] * count)
        arr = np.array(all_fls)
        pcts = [1, 5, 25, 50, 75, 95, 99, 99.9, 99.99, 100]
        print(f"\nSPLICED_ANNOT FL percentiles:")
        for p in pcts:
            print(f"  P{p:6.2f}: {np.percentile(arr, p):.0f}")
    
    if fl_hist_unspliced:
        all_fls = []
        for fl, count in fl_hist_unspliced.items():
            all_fls.extend([fl] * count)
        arr = np.array(all_fls)
        print(f"\nUNSPLICED FL percentiles:")
        for p in [1, 5, 25, 50, 75, 95, 99, 99.9, 99.99, 100]:
            print(f"  P{p:6.2f}: {np.percentile(arr, p):.0f}")

    # Show outlier details
    print(f"\nOutliers (FL > {fl_threshold}): {len(outliers)}")
    
    # Group by splice_type
    by_type = defaultdict(list)
    for o in outliers:
        by_type[o["splice_type"]].append(o)
    
    for st_name, outs in sorted(by_type.items()):
        print(f"\n  {st_name}: {len(outs)} outliers")
        # Sort by max_fl descending
        outs.sort(key=lambda x: x["max_fl"], reverse=True)
        for o in outs[:10]:  # top 10
            n_mrna = sum(1 for x in o["is_nrna"] if not x)
            n_nrna = sum(1 for x in o["is_nrna"] if x)
            mrna_fls = [o["fls"][j] for j, nr in enumerate(o["is_nrna"]) if not nr and o["fls"][j] > 0]
            nrna_fls = [o["fls"][j] for j, nr in enumerate(o["is_nrna"]) if nr and o["fls"][j] > 0]
            # Show transcript names for top outlier
            tx_list = [tx_names[t] for t in o["t_inds"][:5]]
            print(f"    FL={o['max_fl']:>6d} NH={o['nh']} "
                  f"fp={o['footprint']:>6d} rl={o['read_length']:>4d} nm={o['nm']} "
                  f"{n_mrna}m+{n_nrna}n "
                  f"mFL={mrna_fls[:3]} nFL={nrna_fls[:3]} "
                  f"tx={tx_list[:3]}")


# Run Oracle
logger.info("Scanning Oracle...")
scan_cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True)
orc = scan_and_buffer(str(ORACLE_BAM), idx, scan=scan_cfg)
analyze_fl_details(orc[3], is_nrna, "Oracle", fl_threshold=400)

# Run Minimap2
logger.info("Scanning Minimap2...")
mm2 = scan_and_buffer(str(MM2_BAM), idx, scan=scan_cfg)
analyze_fl_details(mm2[3], is_nrna, "Minimap2", fl_threshold=400)

print("\nDone.")
