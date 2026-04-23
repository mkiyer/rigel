#!/usr/bin/env python3
"""Direct profile of the calibration EM step on real data (dna01m buffer)."""
import cProfile, pstats, io, time, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/home/mkiyer/proj/rigel/src")
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.config import BamScanConfig
from rigel.calibration import calibrate_gdna
from rigel.native import detect_sj_strand_tag
from dataclasses import replace

BAM = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna01m/rigel/annotated.bam"
INDEX_DIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"

print("Loading index...", flush=True)
t0 = time.perf_counter()
index = TranscriptIndex.load(INDEX_DIR)
print(f"  {time.perf_counter()-t0:.1f}s", flush=True)

scan_cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True, n_scan_threads=10)
if scan_cfg.sj_strand_tag == "auto":
    scan_cfg = replace(scan_cfg, sj_strand_tag=detect_sj_strand_tag(BAM))

print("Scanning BAM...", flush=True)
t0 = time.perf_counter()
stats, sm, flm, buf, region_counts, fl_table = scan_and_buffer(BAM, index, scan_cfg)
sm.finalize()
flm.build_scoring_models()
flm.finalize(prior_ess=200.0)
print(f"  {time.perf_counter()-t0:.1f}s", flush=True)

print("\n=== Profiling calibrate_gdna ===", flush=True)
pr = cProfile.Profile()
t0 = time.perf_counter()
pr.enable()
cal = calibrate_gdna(
    region_counts, fl_table, index.region_df,
    sm.strand_specificity,
    mean_frag_len=flm.global_model.mean,
    intergenic_fl_model=flm.intergenic,
    fl_prior_ess=200.0,
)
pr.disable()
elapsed = time.perf_counter() - t0
print(f"TOTAL calibrate_gdna: {elapsed:.2f}s", flush=True)
print(f"  lambda_G={cal.lambda_gdna:.3e} kappa_G={cal.kappa_G:.2f} kappa_R={cal.kappa_R:.2f} iters={cal.em_n_iter}", flush=True)

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
ps.print_stats(40)
print(s.getvalue())

s2 = io.StringIO()
ps = pstats.Stats(pr, stream=s2).sort_stats("tottime")
ps.print_stats(25)
print(s2.getvalue())

pr.dump_stats("/tmp/rigel_cal_profile.prof")
