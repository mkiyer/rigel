#!/usr/bin/env python3
"""Deep diagnostic: trace EM candidates and posteriors for LINC02912 reads.

Runs the actual hulkrna pipeline and intercepts the locus EM data for the
locus containing ENST00000624314.1 to inspect candidate log-likelihoods,
nRNA priors, and strand scoring.
"""
import sys
sys.path.insert(0, "src")

import numpy as np
from hulkrna.index import HulkIndex

INDEX_DIR = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/"
    "PVT1_MYC/hulkrna_index"
)

index = HulkIndex.load(INDEX_DIR)
t_df = index.t_df
g_df = index.g_df

# Find ENST00000624314.1
t624_mask = t_df["t_id"] == "ENST00000624314.1"
t624_idx = int(t_df.index[t624_mask][0])
t624_gidx = int(index.t_to_g_arr[t624_idx])
t624_strand = int(index.t_to_strand_arr[t624_idx])
g624_strand = int(index.g_to_strand_arr[t624_gidx])

print(f"ENST00000624314.1: t_idx={t624_idx}")
print(f"  t_strand={t624_strand} (1=POS, 2=NEG)")
print(f"  gene={g_df.iloc[t624_gidx]['g_id']} ({g_df.iloc[t624_gidx]['g_name']})")
print(f"  g_strand={g624_strand}")
print(f"  span={t_df.iloc[t624_idx]['start']}-{t_df.iloc[t624_idx]['end']}")
print(f"  length={t_df.iloc[t624_idx]['length']}")

# Check PVT1 transcripts covering LINC02912 region
pvt1_mask = t_df["g_id"] == "ENSG00000249859.14"
pvt1_indices = t_df.index[pvt1_mask].values
pvt1_covering = []
for ti in pvt1_indices:
    ts = int(t_df.iloc[ti]["start"])
    te = int(t_df.iloc[ti]["end"])
    if ts <= 1503283 and te >= 1501119:
        pvt1_covering.append(ti)

pvt1_gidx = int(index.t_to_g_arr[pvt1_indices[0]])
pvt1_strand = int(index.g_to_strand_arr[pvt1_gidx])
print(f"\nPVT1: {len(pvt1_indices)} transcripts, g_strand={pvt1_strand}")
print(f"  {len(pvt1_covering)} transcripts span LINC02912 region")
for ti in pvt1_covering[:5]:
    row = t_df.iloc[ti]
    print(f"    t_idx={ti} {row['t_id']}  t_strand={int(index.t_to_strand_arr[ti])}  span={row['start']}-{row['end']}")

# Now run the pipeline and intercept the scan data
print("\n" + "="*60)
print("Running pipeline to check scan data...")
print("="*60)

from hulkrna.pipeline import run_pipeline
from hulkrna.config import PipelineConfig, EMConfig, ScanConfig, ScoringConfig

BAM = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam"
)

cfg = PipelineConfig(
    em=EMConfig(),
    scan=ScanConfig(),
    scoring=ScoringConfig(),
)
pipe = run_pipeline(BAM, index, config=cfg)

counter = pipe.counter
em_data = pipe.em_data
loci = pipe.loci

# Find the locus containing ENST00000624314.1
target_locus = None
for locus in loci:
    if t624_idx in locus.transcript_indices:
        target_locus = locus
        break

if target_locus is None:
    print("ERROR: No locus contains ENST00000624314.1!")
    sys.exit(1)

t_arr = target_locus.transcript_indices
n_t = len(t_arr)
print(f"\nLocus: {n_t} transcripts, {len(target_locus.unit_indices)} units")

# Find local index of ENST00000624314.1
local_624 = None
for i, gi in enumerate(t_arr):
    if gi == t624_idx:
        local_624 = i
        break
print(f"ENST00000624314.1 local index: {local_624}")

# Check nrna_init
print(f"\nnrna_init[ENST00000624314.1] = {counter.nrna_init[t624_idx]}")
print(f"transcript_intronic_sense[t624] = {counter.transcript_intronic_sense[t624_idx]}")
print(f"transcript_intronic_antisense[t624] = {counter.transcript_intronic_antisense[t624_idx]}")

# Check PVT1 nrna_init
pvt1_nrna_sum = counter.nrna_init[pvt1_indices].sum()
pvt1_nrna_nonzero = (counter.nrna_init[pvt1_indices] > 0).sum()
print(f"\nPVT1 nrna_init: sum={pvt1_nrna_sum:.2f}, nonzero={pvt1_nrna_nonzero}/{len(pvt1_indices)}")

# Check the EM input
from hulkrna.locus import build_locus_em_input

em_input = build_locus_em_input(
    target_locus, em_data, counter, index,
    mean_frag=200.0, gdna_init=0.0,
)

prior = em_input.prior
unique_totals = em_input.unique_totals
nrna_base_local = n_t

print(f"\nEM Input: {em_input.n_components} components")
print(f"  mRNA [0, {n_t})  nRNA [{n_t}, {2*n_t})  gDNA [{2*n_t}]")

# ENST00000624314.1 nRNA shadow
nrna_624 = nrna_base_local + local_624
print(f"\nENST00000624314.1 nRNA shadow (local component {nrna_624}):")
print(f"  prior = {prior[nrna_624]}")
print(f"  unique_totals = {unique_totals[nrna_624]}")

# mRNA of ENST00000624314.1
print(f"\nENST00000624314.1 mRNA (local component {local_624}):")
print(f"  prior = {prior[local_624]}")
print(f"  unique_totals = {unique_totals[local_624]}")

# Check how many EM units have candidates pointing to nRNA_624
offsets = em_input.offsets
t_indices = em_input.t_indices
log_liks = em_input.log_liks
n_units = len(offsets) - 1

units_with_nrna_624 = []
units_with_mrna_624 = []
for u in range(n_units):
    s = int(offsets[u])
    e = int(offsets[u + 1])
    seg_tidx = t_indices[s:e]
    seg_lls = log_liks[s:e]
    if nrna_624 in seg_tidx:
        idx_in_seg = int(np.where(seg_tidx == nrna_624)[0][0])
        nrna_ll = seg_lls[idx_in_seg]
        # Also find mRNA 624 if present
        mrna_ll = None
        if local_624 in seg_tidx:
            mrna_idx = int(np.where(seg_tidx == local_624)[0][0])
            mrna_ll = seg_lls[mrna_idx]
        units_with_nrna_624.append((u, nrna_ll, mrna_ll, len(seg_tidx)))
    if local_624 in seg_tidx:
        idx_in_seg = int(np.where(seg_tidx == local_624)[0][0])
        mrna_ll = seg_lls[idx_in_seg]
        units_with_mrna_624.append((u, mrna_ll, len(seg_tidx)))

print(f"\nUnits with nRNA candidate for ENST00000624314.1: {len(units_with_nrna_624)}")
print(f"Units with mRNA candidate for ENST00000624314.1: {len(units_with_mrna_624)}")

if units_with_nrna_624:
    print("\nSample units with nRNA ENST00000624314.1:")
    for u, nrna_ll, mrna_ll, n_cand in units_with_nrna_624[:5]:
        s = int(offsets[u])
        e = int(offsets[u + 1])
        seg_tidx = t_indices[s:e]
        seg_lls = log_liks[s:e]
        # Show all candidate types
        mrna_cands = seg_tidx[seg_tidx < n_t]
        nrna_cands = seg_tidx[(seg_tidx >= n_t) & (seg_tidx < 2*n_t)]
        gdna_cands = seg_tidx[seg_tidx == 2*n_t]
        print(f"  unit={u}: {n_cand} cands ({len(mrna_cands)} mRNA, {len(nrna_cands)} nRNA, {len(gdna_cands)} gDNA)")
        print(f"    nRNA ENST00000624314.1 LL = {nrna_ll:.4f}")
        if mrna_ll is not None:
            print(f"    mRNA ENST00000624314.1 LL = {mrna_ll:.4f}")
        # Show best mRNA and best nRNA candidates
        mrna_mask = seg_tidx < n_t
        nrna_mask = (seg_tidx >= n_t) & (seg_tidx < 2*n_t)
        if mrna_mask.any():
            best_mrna_idx = int(np.argmax(seg_lls[mrna_mask]))
            best_mrna_ll = float(seg_lls[mrna_mask][best_mrna_idx])
            best_mrna_t = int(seg_tidx[mrna_mask][best_mrna_idx])
            gt = int(t_arr[best_mrna_t])
            print(f"    best mRNA: local={best_mrna_t} ({t_df.iloc[gt]['t_id']}) LL={best_mrna_ll:.4f}")
        if nrna_mask.any():
            best_nrna_idx = int(np.argmax(seg_lls[nrna_mask]))
            best_nrna_ll = float(seg_lls[nrna_mask][best_nrna_idx])
            best_nrna_t = int(seg_tidx[nrna_mask][best_nrna_idx])
            gt = int(t_arr[best_nrna_t - n_t])
            print(f"    best nRNA: local={best_nrna_t} ({t_df.iloc[gt]['t_id']}) LL={best_nrna_ll:.4f}")
        # Show ALL nRNA candidates and their LLs for this unit
        if nrna_mask.any():
            nrna_ts = seg_tidx[nrna_mask]
            nrna_lls = seg_lls[nrna_mask]
            print(f"    ALL nRNA candidates ({len(nrna_ts)}):")
            sorted_idx = np.argsort(-nrna_lls)
            for j in sorted_idx[:10]:
                nt = int(nrna_ts[j])
                gt = int(t_arr[nt - n_t])
                gid = t_df.iloc[gt]["g_id"]
                gname = t_df.iloc[gt]["g_name"]
                ts = int(index.t_to_strand_arr[gt])
                print(f"      local={nt} ({t_df.iloc[gt]['t_id']}, {gname}) t_strand={ts} LL={float(nrna_lls[j]):.4f}")

if units_with_mrna_624:
    print(f"\nSample units with mRNA ENST00000624314.1:")
    for u, mrna_ll, n_cand in units_with_mrna_624[:5]:
        s = int(offsets[u])
        e = int(offsets[u + 1])
        seg_tidx = t_indices[s:e]
        seg_lls = log_liks[s:e]
        mrna_cands = seg_tidx[seg_tidx < n_t]
        nrna_cands = seg_tidx[(seg_tidx >= n_t) & (seg_tidx < 2*n_t)]
        print(f"  unit={u}: {n_cand} cands ({len(mrna_cands)} mRNA, {len(nrna_cands)} nRNA)")
        print(f"    mRNA ENST00000624314.1 LL = {mrna_ll:.4f}")

# Check PVT1 nRNA priors
pvt1_local_indices = []
for i, gi in enumerate(t_arr):
    if int(gi) in pvt1_indices:
        pvt1_local_indices.append(i)
pvt1_nrna_local = [n_t + i for i in pvt1_local_indices]
pvt1_nrna_priors = prior[pvt1_nrna_local]
pvt1_nrna_nonzero_prior = (pvt1_nrna_priors > 0).sum()
print(f"\nPVT1 nRNA priors: {pvt1_nrna_nonzero_prior}/{len(pvt1_local_indices)} nonzero")
if pvt1_nrna_nonzero_prior > 0:
    nz_idx = np.where(pvt1_nrna_priors > 0)[0]
    for idx in nz_idx[:5]:
        li = pvt1_local_indices[idx]
        gi = int(t_arr[li])
        print(f"  local={li} ({t_df.iloc[gi]['t_id']}) prior={pvt1_nrna_priors[idx]:.6f} nrna_init={counter.nrna_init[gi]:.4f}")

# Show all components with nonzero prior
nz_prior = np.where(prior > 0)[0]
print(f"\nComponents with nonzero prior: {len(nz_prior)} / {em_input.n_components}")
for c in nz_prior:
    if c < n_t:
        gt = int(t_arr[c])
        label = f"mRNA {t_df.iloc[gt]['t_id']} ({t_df.iloc[gt]['g_name']})"
    elif c < 2*n_t:
        gt = int(t_arr[c - n_t])
        label = f"nRNA {t_df.iloc[gt]['t_id']} ({t_df.iloc[gt]['g_name']})"
    else:
        label = "gDNA"
    if unique_totals[c] > 0.1:
        print(f"  comp={c:4d} {label:60s} prior={prior[c]:.8f} unique_totals={unique_totals[c]:.2f}")
