"""Analyze locus size distribution from a BAM+index to understand mega-locus bottleneck.

Usage:
    conda activate rigel
    python scripts/debug/analyze_locus_sizes.py
"""

import numpy as np
import time

from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _setup_geometry_and_estimator
from rigel.config import PipelineConfig, BamScanConfig
from rigel.scoring import FragmentScorer
from rigel.scan import FragmentRouter
from rigel.locus import build_loci
from rigel.stats import PipelineStats
from rigel.calibration import calibrate_gdna

BAM = "/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_salmon/gdna_high_ss_0.90_nrna_default/align_minimap2/reads_namesort.bam"
INDEX = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v5_phase3b/rigel_index"

print("Loading index...")
index = TranscriptIndex.load(INDEX)

scan_cfg = BamScanConfig(n_scan_threads=8)
print("Running scan...")
t0 = time.time()
stats, strand_models, frag_length_models, buffer, region_counts, fl_table = scan_and_buffer(
    BAM, index, scan_cfg
)
print(f"Scan: {time.time()-t0:.1f}s, {buffer.total_fragments:,} fragments")

# Finalize models (required before calibration/scoring)
strand_models.finalize()
frag_length_models.build_scoring_models()
frag_length_models.finalize()

# Calibration
print("Calibrating...")
t0 = time.time()
calibration = calibrate_gdna(
    region_counts,
    fl_table,
    index.region_df,
    strand_models.strand_specificity,
    intergenic_fl_model=frag_length_models.intergenic,
)
print(f"Calibration: {time.time()-t0:.1f}s")

# Apply calibrated gDNA FL model
frag_length_models.gdna_model = calibration.gdna_fl_model

# Set up geometry + estimator
cfg = PipelineConfig()
geometry, estimator = _setup_geometry_and_estimator(index, frag_length_models, cfg.em)

# Build scorer and route fragments
print("Building scorer + routing fragments...")
scorer = FragmentScorer.from_models(strand_models, frag_length_models, index, estimator)
router = FragmentRouter(scorer, estimator, stats, index, strand_models)
t0 = time.time()
em_data = router.scan(buffer, log_every=100_000_000)
print(f"Router scan: {time.time()-t0:.1f}s, {em_data.n_units:,} units")
buffer.release()

# Build loci
print("Building loci...")
t0 = time.time()
loci = build_loci(em_data, index)
print(f"build_loci: {time.time()-t0:.1f}s, {len(loci)} loci")

# ---- Size distribution ----
sizes = np.array([len(loc.transcript_indices) for loc in loci])
units_per = np.array([len(loc.unit_indices) for loc in loci])
print(f"\nLocus size distribution:")
print(f"  Total loci: {len(loci)}")
print(f"  Median transcripts: {np.median(sizes):.0f}")
print(f"  Mean transcripts: {np.mean(sizes):.0f}")
print(f"  Max transcripts: {np.max(sizes)}")
print(f"  Max units: {np.max(units_per)}")

order = np.argsort(sizes)[::-1]
print(f"\n  Top-10 loci by transcripts:")
for i in range(min(10, len(order))):
    li = order[i]
    print(
        f"    locus {loci[li].locus_id}: "
        f"{sizes[li]} tx, {units_per[li]:,} units, "
        f"gdna_span={loci[li].gdna_span:,}"
    )

print(f"\n  Size buckets:")
for lo, hi in [
    (1, 1), (2, 5), (6, 10), (11, 50), (51, 100),
    (101, 1000), (1001, 10000), (10001, 100000), (100001, 1000000),
]:
    mask = (sizes >= lo) & (sizes <= hi)
    n = np.sum(mask)
    if n > 0:
        total_units = units_per[mask].sum()
        print(f"    [{lo:>7}-{hi:>7}]: {n:>6} loci, {total_units:>12,} units")

# Estimate per-locus EM work
print("\n  Estimated EM time distribution (units × components heuristic):")
components = sizes + 1  # n_t + 1 gDNA
work = units_per.astype(np.float64) * components
total_work = work.sum()
for i in range(min(5, len(order))):
    li = order[i]
    pct = work[li] / total_work * 100
    print(
        f"    locus {loci[li].locus_id}: "
        f"{sizes[li]} tx × {units_per[li]:,} units = "
        f"{work[li]:,.0f} ops ({pct:.1f}%)"
    )
print(f"  Total estimated work: {total_work:,.0f}")
top1_pct = work[order[0]] / total_work * 100
print(f"  Top-1 locus share: {top1_pct:.1f}%")
top5_pct = work[order[:5]].sum() / total_work * 100
print(f"  Top-5 loci share: {top5_pct:.1f}%")
