#!/usr/bin/env python3
"""Check if any EC has nRNA ENST000624314 (comp 311) without mRNA (comp 118)."""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, "src")

from hulkrna.index import HulkIndex
from hulkrna.config import PipelineConfig
import hulkrna.estimator as estimator_mod

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")

TARGET_GLOBAL_T = 290  # ENST00000624314.1

original_run_locus_em = estimator_mod.AbundanceEstimator.run_locus_em

def instrumented_run_locus_em(self, locus_em, **kwargs):
    t_arr = locus_em.local_to_global_t
    n_t = locus_em.n_transcripts
    
    if TARGET_GLOBAL_T in t_arr:
        local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
        nrna_local = n_t + local_idx  # 311
        mrna_local = local_idx  # 118
        gdna_idx = 2 * n_t  # 386
        
        print(f"\n{'='*70}")
        print(f"EC ANALYSIS for nRNA comp={nrna_local}, mRNA comp={mrna_local}")
        
        offsets = locus_em.offsets
        t_indices = locus_em.t_indices
        n_units = len(offsets) - 1
        
        # Check each unit: does it have nrna_local? Does it also have mrna_local?
        units_with_nrna = 0
        units_nrna_only = 0  # has nrna but NOT mrna
        units_nrna_with_mrna = 0
        nrna_only_unit_sizes = []
        
        for u in range(n_units):
            s = int(offsets[u])
            e = int(offsets[u + 1])
            if s == e:
                continue
            seg = t_indices[s:e]
            has_nrna = nrna_local in seg
            has_mrna = mrna_local in seg
            if has_nrna:
                units_with_nrna += 1
                if has_mrna:
                    units_nrna_with_mrna += 1
                else:
                    units_nrna_only += 1
                    nrna_only_unit_sizes.append(e - s)
                    if units_nrna_only <= 5:
                        # Print details about this unit
                        seg_list = seg.tolist()
                        has_any_mrna = any(c < n_t for c in seg_list)
                        has_gdna = gdna_idx in seg_list
                        n_mrna = sum(1 for c in seg_list if c < n_t)
                        n_nrna = sum(1 for c in seg_list if n_t <= c < gdna_idx)
                        print(f"\n  nRNA-only unit #{units_nrna_only}: "
                              f"{e-s} candidates, {n_mrna} mRNA, {n_nrna} nRNA, "
                              f"has_gdna={has_gdna}, has_any_mrna={has_any_mrna}")
                        print(f"    comps: {seg_list[:20]}{'...' if len(seg_list)>20 else ''}")
                        # Check log_liks
                        lls = locus_em.log_liks[s:e]
                        nrna_pos = list(seg).index(nrna_local)
                        print(f"    ll[nRNA comp {nrna_local}] = {lls[nrna_pos]:.4f}")
                        if has_any_mrna:
                            for j, c in enumerate(seg_list):
                                if c < n_t:
                                    print(f"    ll[mRNA comp {c}] = {lls[j]:.4f}")
        
        print(f"\n  Summary:")
        print(f"    Units with nRNA comp {nrna_local}: {units_with_nrna}")
        print(f"    ... with mRNA comp {mrna_local}: {units_nrna_with_mrna}")
        print(f"    ... WITHOUT mRNA comp {mrna_local}: {units_nrna_only}")
        
        if nrna_only_unit_sizes:
            print(f"    nRNA-only unit sizes: {nrna_only_unit_sizes[:20]}")
        
        # Now check: units that have ONLY nRNA candidates (no mRNA at all)
        units_all_nrna = 0
        for u in range(n_units):
            s = int(offsets[u])
            e = int(offsets[u + 1])
            if s == e:
                continue
            seg = t_indices[s:e]
            has_nrna_comp = nrna_local in seg
            if not has_nrna_comp:
                continue
            has_any_mrna_comp = any(int(c) < n_t for c in seg)
            has_gdna_comp = gdna_idx in seg
            if not has_any_mrna_comp and not has_gdna_comp:
                units_all_nrna += 1
        
        print(f"    Units with nRNA {nrna_local} and NO mRNA/gDNA at all: {units_all_nrna}")
        
        # Also count: how many total units have NO mRNA candidates?
        total_no_mrna = 0
        total_no_mrna_has_nrna = 0
        for u in range(n_units):
            s = int(offsets[u])
            e = int(offsets[u + 1])
            if s == e:
                continue
            seg = t_indices[s:e]
            has_any_mrna_comp = any(int(c) < n_t for c in seg)
            if not has_any_mrna_comp:
                total_no_mrna += 1
                if nrna_local in seg:
                    total_no_mrna_has_nrna += 1
        
        print(f"    Total units with NO mRNA candidates: {total_no_mrna}")
        print(f"    ... of which have nRNA {nrna_local}: {total_no_mrna_has_nrna}")
        
        print(f"{'='*70}\n")
    
    return original_run_locus_em(self, locus_em, **kwargs)

estimator_mod.AbundanceEstimator.run_locus_em = instrumented_run_locus_em

print("Running pipeline to analyze EC structure...")
from hulkrna.pipeline import run_pipeline
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=cfg)
print(f"\nFinal: nrna_em_counts[290] = {pipe.estimator.nrna_em_counts[290]:.1f}")
