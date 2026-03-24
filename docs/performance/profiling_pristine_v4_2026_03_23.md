# Rigel C/C++ Performance Profiling Report (Pristine v4)

Date: 2026-03-23  
Author: GitHub Copilot

## 1. Goal

Profile current Rigel performance on a realistic simulated BAM, with focus on native/C++ runtime and memory behavior, then propose concrete optimization opportunities.

Input BAM:
- /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam

Index:
- /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index

## 2. Profiler status and fix applied

The existing stage profiler failed due to API drift:
- scripts/profiler.py used removed fields/args: TranscriptIndex.num_nrna and t_to_nrna
- scripts/profiler.py expected 4 return values from run_batch_locus_em, but current API returns 3

Fix applied in scripts/profiler.py:
- AbundanceEstimator constructor now uses is_synthetic_nrna=index.t_df["is_synthetic_nrna"].values
- profile_stages now unpacks run_batch_locus_em as:
  (total_gdna_em, locus_mrna_arr, locus_gdna_arr)

This restored stage profiling functionality.

## 3. Profiling runs executed

Commands run:

  conda activate rigel
  python scripts/profiler.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index \
    --stages --memory-interval 50 \
    --outdir /Users/mkiyer/proj/rigel/results/profile_pristine_v4_2026_03_23

  conda activate rigel
  python scripts/profiler.py ... --stages --threads 1 \
    --outdir /Users/mkiyer/proj/rigel/results/profile_pristine_v4_2026_03_23_t1

  conda activate rigel
  python scripts/profiler.py ... --stages --threads 8 \
    --outdir /Users/mkiyer/proj/rigel/results/profile_pristine_v4_2026_03_23_t8

Artifacts:
- results/profile_pristine_v4_2026_03_23/profile_summary.json
- results/profile_pristine_v4_2026_03_23/profile_report.txt
- results/profile_pristine_v4_2026_03_23/memory_timeline_default.csv
- results/profile_pristine_v4_2026_03_23_t1/profile_summary.json
- results/profile_pristine_v4_2026_03_23_t8/profile_summary.json

## 4. High-level results

Workload facts (all runs):
- 457,513 transcripts, 63,472 genes
- 31,242,867 total BAM records processed
- 9,999,650 fragment groups
- 17,685,478 buffered units (multimap expansion)
- 10,386 loci
- Largest locus: 91,456 transcripts, 4,573,447 units

### Runtime and memory summary

| Config | Wall time (s) | Throughput (frags/s) | Peak RSS (MB) | Scan (s) | Calibration (s) | Router scan (s) | Locus EM (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| threads=1 | 357.327 | 87,435 | 13,333 | 244.674 | 9.962 | 9.462 | 88.207 |
| default | 123.678 | 252,614 | 14,730 | 84.189 | 9.938 | 10.711 | 14.149 |
| threads=8 | 112.910 | 276,706 | 14,247 | 74.557 | 9.311 | 10.088 | 14.148 |

Scaling t1 -> t8:
- Wall speedup: 3.16x
- Scan speedup: 3.28x
- Locus EM speedup: 6.23x
- Calibration speedup: 1.07x (effectively serial)
- Router scan speedup: 0.94x (does not benefit from extra threads)

## 5. Deep analysis

## 5.1 Runtime bottlenecks on this dataset

Default run stage share:
- scan_and_buffer: 84.19s (68.1%)
- locus_em: 14.15s (11.4%)
- fragment_router_scan: 10.71s (8.7%)
- calibration: 9.94s (8.0%)
- everything else: ~3.8%

Interpretation:
1. C++ BAM scan/resolve dominates runtime for this workload.
2. Locus EM is still material, but much less dominant than in prior mega-locus-heavy runs.
3. Calibration is now non-trivial (8%), but not primary.
4. Router scan behaves as a serial stage in this setup (likely bound by single-thread native scoring path and Python orchestration overhead).

## 5.2 Mega-locus effect still visible

Even with better minimap2 settings and pristine scenario, largest locus is still very large:
- 91,456 transcripts
- 4.57M units

Locus EM cost remains meaningful (14.1s at default/8 threads). This locus scale will continue to be a tail-latency and memory-pressure risk in real data.

## 5.3 Memory behavior

Observed:
- RSS before: ~3.2 GB
- Peak RSS: 14.2-14.7 GB
- Peak reached near end of scan stage (~74-84s in multithreaded runs)
- RSS snapshots remain near peak after scan/buffer release/cleanup in this profiler run

Most important memory signal:
- Scanner spilled a single huge chunk: 17,685,478 units, 2,743.5 MB Arrow file.
- Yet process RSS rose by ~11 GB beyond baseline.

Likely causes:
1. Large in-memory native accumulation before spill/finalize (single giant accumulator lifecycle).
2. Arena/allocator retention on macOS after free (RSS not promptly returned).
3. Temporary duplication during finalize/serialization boundaries (native -> Python arrays -> Arrow).

## 5.4 Threading behavior

- Scan stage: good but sublinear scaling (3.28x at 8 threads), suggesting queue/decompression/memory bandwidth limits.
- Locus EM stage: strong scaling (6.23x), indicating C++ parallel EM is effective.
- Router stage: flat/slightly worse with more threads, suggesting it is effectively serial in current architecture.

## 6. Ideal detailed C/C++ profiling approach

The current stage profiler is good for pipeline-level attribution, but not enough for native function-level root cause. Recommended stack:

1. Build profileable native binaries (critical)
- Use RelWithDebInfo and preserve symbols for native modules.
- Avoid stripping debug symbols in profiling builds.

Suggested build mode for profiling sessions:

  conda activate rigel
  CMAKE_BUILD_TYPE=RelWithDebInfo pip install --no-build-isolation -e .

2. Native CPU profiling (macOS)
- Use Instruments Time Profiler (CLI via xctrace) on rigel quant command/workload.
- Capture separate traces for:
  - scan-heavy window
  - locus-EM-heavy window
- Export symbolized call trees and aggregate by self/cumulative time.

3. Native memory profiling
- Instruments Allocations + VM Tracker (macOS) to identify largest allocation sites and lifetime.
- Validate whether memory retention is true leaks, long-lived objects, or allocator caching.

4. Linux production-grade profiling (when moving to Linux)
- perf record/report + flamegraphs (CPU)
- heaptrack or massif (allocation/lifetime)
- Keep a stable benchmark BAM + index for A/B comparisons.

5. Add in-code stage/substage native timers
- For C++ scan: parse, resolve, multimap pairing, append/serialization.
- For C++ EM: build inputs, E-step, M-step, SQUAREM overhead, posterior assignment.
- Emit structured counters into profiler JSON.

## 7. Priority optimization opportunities

## 7.1 Highest impact: scan memory architecture

Problem:
- Peak RSS is dominated by scan accumulation lifecycle.

Opportunity:
- Stream from native scanner into bounded chunk buffers with immediate spill, instead of one giant end-of-scan finalize.

Expected benefit:
- Major peak RSS reduction (multi-GB; likely largest memory win).
- Potential runtime improvement from less allocator pressure and reduced large-copy overhead.

## 7.2 High impact: scan throughput

Problem:
- Scan is 68% of wall time and scales only 3.28x at 8 threads.

Opportunities:
1. Increase/de-couple decompression threading from worker count and tune dynamically.
2. Reduce per-group dynamic allocation churn in parse/group/build paths (preallocation/pool reuse).
3. Revisit queue/batching granularity to reduce producer/consumer overhead.
4. Consider SIMD/data-layout optimization in overlap resolution hot loops.

Expected benefit:
- 10-25% end-to-end wall-time reduction is plausible if scan improves materially.

## 7.3 Medium impact: router stage parallelization

Problem:
- fragment_router_scan is ~10s and effectively not scaling with threads.

Opportunities:
1. Parallelize fused_score_buffer across chunks/row blocks in native scoring path.
2. Reduce Python-side per-chunk loops for annotation handling by vectorized/native callbacks.
3. Evaluate lock-free or batched updates for deterministic-unambig accumulation.

Expected benefit:
- 3-8% end-to-end gain depending on achievable parallelization.

## 7.4 Medium impact: mega-locus handling

Problem:
- One mega-locus still has 91k transcripts and 4.57M units.

Opportunities:
1. Add optional pre-EM partition heuristics for weakly connected mega-components.
2. Candidate capping/pruning thresholds for ultra-large loci (careful with bias).
3. Incremental EM or block-EM strategy for giant loci.

Expected benefit:
- Better tail latency and memory stability; potential major wins on difficult real datasets.

## 7.5 Low-medium impact: calibration vectorization still worthwhile

Problem:
- Calibration is ~8% here; serial behavior persists.

Opportunity:
- Keep prior vectorization plan (gammaln/betaln/vectorized likelihood loops) for consistent gains across datasets.

Expected benefit:
- Moderate runtime reduction (larger on calibration-heavy datasets).

## 8. Proposed execution plan

Phase A (quick wins, 1-3 days):
1. Keep profiler.py API-compatible (done in this report).
2. Add native substage timers to scan and EM C++ code.
3. Tune decompression/thread settings and benchmark.

Phase B (major memory work, 3-7 days):
1. Implement streaming native-to-buffer chunk handoff.
2. Spill incrementally; avoid giant in-memory accumulator.
3. Re-profile peak RSS and runtime.

Phase C (algorithmic scaling, 1-2 weeks):
1. Mega-locus mitigation strategy.
2. Router path parallelization.
3. Regression checks on accuracy + benchmark suite.

## 9. Bottom line

For this target BAM, Rigel is already strongly parallel in locus EM, but total runtime and memory are dominated by C++ scan/accumulation behavior.

Best immediate ROI:
1. Redesign scan accumulation/spill lifecycle for bounded memory.
2. Further optimize scan-stage parallel efficiency.
3. Add native-level profiler visibility (Time Profiler + native substage counters) to drive next optimizations with function-level evidence.
