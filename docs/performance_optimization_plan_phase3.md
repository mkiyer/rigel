# Performance Optimization Plan — Phase 3

## Profiling Results (2026-03-06)

**Dataset**: 24M fragments (12M PE reads), 254K transcripts, 63K genes, simulated CCLE HeLa  
**Machine**: macOS 26.3 arm64 (Apple Silicon), 4 EM threads  
**Config**: production `run_pipeline()` path (batch C++ EM)

### Headline Numbers

| Metric               | Before Fused (Phase 2) | After Fused | Delta           |
|----------------------|------------------------|-------------|-----------------|
| **Peak RSS**         | 18,548 MB              | 13,883 MB   | **−4,665 MB (−25%)** |
| **Wall time**        | ~139 s                 | 131.5 s     | **−7.5 s (−5.4%)**   |
| **Throughput**       | 174 K frags/s          | 182 K frags/s | +5%           |
| `fragment_router_scan` | 6.1 s                | 4.1 s       | −2.0 s (−33%)  |
| `build_loci`         | 3.9 s                  | 2.1 s       | −1.8 s (−46%)  |
| `scan_and_buffer`    | 92.3 s                 | 91.5 s      | −0.8 s          |
| `locus_em`           | 31.4 s                 | 31.3 s      | similar         |

**Takeaway**: The fused two-pass C++ score buffer delivered a massive 4.7 GB memory
reduction (from 18.5 → 13.9 GB) by eliminating the triple-copy Python-C++ data path
and producing float32 arrays with zero-copy capsule return. Peak RSS is now **below
16 GB**, meeting the original memory target. The router scan and locus build both
sped up due to reduced allocation pressure and better cache behavior.

### Per-Stage Breakdown (--stages run)

| Stage                   | Time     | %     |
|-------------------------|----------|-------|
| `scan_and_buffer`       | 91.515 s | 69.6% |
| `locus_em`              | 31.320 s | 23.8% |
| `fragment_router_scan`  |  4.074 s |  3.1% |
| `build_loci`            |  2.124 s |  1.6% |
| `eb_gdna_priors`        |  1.711 s |  1.3% |
| `fragment_scorer`       |  0.717 s |  0.5% |
| `create_estimator`      |  0.075 s |  0.1% |
| **TOTAL**               |**131.5 s**| 100% |

### Memory Timeline

```
t=0s:    3,841 MB   (index loaded)
t=0s:    7,438 MB   (scan starts, initial buffer allocation)
t=91s:   9,649 MB   (end of scan, buffer at 11.2M frags → spilled 2.1 GB to disk)
t=92s:  11,827 MB   (chunk reloaded from disk + scoring begins)
t=93s:  13,883 MB   (peak: ScoredFragments allocated + EM starts)
t=132s: 13,883 MB   (flat through EM — macOS doesn't reclaim freed pages)
```

### cProfile Hotspots (production path)

| Func                         | Self Time | Calls    | Notes                          |
|------------------------------|-----------|----------|--------------------------------|
| `run_batch_locus_em`         | 0.36 s    | 1        | Python overhead only; C++ runs in threads (~31 s wall) |
| `_build_locus_meta`          | 0.18 s    | 29,266   | 1.68 s cumtime; pandas Arrow `__getitem__` dominates |
| `build_loci`                 | 0.60 s    | 1        | Python loop over 10.3 M units → `list.append` |
| `compute_eb_gdna_priors`     | 0.33 s    | 1        | EB hyperparameter estimation   |
| `scoring.py:180-181 genexpr` | 0.49 s    | 3.8 M    | exon data tuple comprehensions in `from_models` |
| `ndarray.copy`               | 0.23 s    | 40       | chunk reload copies            |
| `gc.collect`                 | 0.69 s    | 2        | explicit GC after del          |
| `list.append`                | 0.76 s    | 12.4 M   | locus.py component_map + unit assignment loops |

---

## Remaining Bottleneck Analysis

### 1. scan_and_buffer — 91.5 s (70%)

This is BAM I/O + C++ resolution + model training + buffering + disk spill.  
Breakdown (approximate from logs):
- **C++ native scan** (htslib BAM decode → overlap resolution → model
  observation collection → buffer fill): ~90 s
- **Spill to Arrow** (write 2.1 GB LZ4-compressed feather): ~0.3 s
- **Model finalize**: negligible

This stage is **I/O-bound** on BAM decode and **compute-bound** on interval
overlap resolution. The current design is single-threaded for BAM reading
(htslib sequential scan of name-sorted BAM). Multi-threaded BAM reading
would require coordinate-sorted input or a sharding strategy, which is a
major architectural change.

### 2. locus_em — 31.3 s (24%)

Batch C++ EM with 4 OpenMP threads, 29,266 loci, largest locus has 297
transcripts and 59,930 units. The C++ code performs per-locus
`extract_locus_sub_problem` + iterative EM.

### 3. Everything else — 8.7 s (6.6%)

`fragment_router_scan` (4.1 s), `build_loci` (2.1 s), `eb_gdna_priors`
(1.7 s), `fragment_scorer` (0.7 s), misc.

---

## Proposed Phase 3 Optimizations

### P3-A: Eliminate buffer spill round-trip (HIGH impact, MEDIUM effort)

**Problem**: The buffer spills 11.2 M fragments (2.1 GB) to Arrow on disk,
then reloads them immediately after the scan for scoring. This write+read
round-trip costs ~0.9 s in I/O and forces ~2 GB of temporary peak memory
for the deserialized copy.

**Proposal**: Add a **memory budget** option so that when sufficient RAM is
available, the buffer stays in memory (no spill). The spill path remains as
a fallback for memory-constrained machines.

- Add `max_memory_mb` parameter to `BamScanConfig` (default: auto-detect
  available memory).
- In `FragmentBuffer.finalize()`, skip spill if current RSS + buffer size
  < budget.
- Saves ~0.9 s wall time and ~0.3 s of feather write + 0.6 s of feather
  read.
- More importantly, avoids the temporary memory spike from having both the
  in-memory buffer and the reloaded copy during the load phase.

**Expected gain**: −0.9 s runtime, −500 MB peak memory (eliminate reload
allocation overhead).

### P3-B: Port `build_loci` to C++ (MEDIUM impact, LOW effort)

**Problem**: `build_loci` takes 2.1 s, spending most of its time in Python
loops with 12.4 M `list.append` calls to build `component_map` and
`root_to_units` dicts.

**Proposal**: Port the grouping logic after C++ `_cc_native` into C++:

1. C++ receives `labels[n_transcripts]` and `unit_labels[n_units]` arrays.
2. C++ sorts and groups into per-component arrays using `std::sort` +
   `std::partition`.
3. Returns packed CSR-style arrays: `locus_offsets_t`, `locus_t_indices`,
   `locus_offsets_u`, `locus_u_indices`.
4. Python wraps into `Locus` objects with zero slicing.

**Expected gain**: −1.5 s (2.1 → 0.6 s).

### P3-C: Port `_build_locus_meta` to C++ (LOW impact, LOW effort)

**Problem**: 29,266 calls to `_build_locus_meta` accumulate 1.68 s, dominated
by pandas Arrow `__getitem__` (0.94 s) for per-transcript chromosome lookups.

**Proposal**: Pre-extract `t_refs` as a numpy array (or integer-coded), then
compute primary_chrom and gene_set counts in C++ or vectorized numpy. Return
a structured array instead of building 29K Python dicts.

**Expected gain**: −1.3 s.

### P3-D: Port `FragmentScorer.from_models` exon iteration to C++ (LOW impact,
LOW effort)

**Problem**: `from_models` iterates over 254K transcripts in Python to build
`t_exon_data` dict with tuple comprehensions (lines 180-181), costing ~0.6 s
with 3.8 M generator calls.

**Proposal**: Build the exon coordinate arrays once in C++ during index load
(or lazily on first scorer creation) and store as packed CSR arrays. The
`FragmentScorer` already passes these to C++ scoring — eliminate the Python
materialization step.

**Expected gain**: −0.5 s.

### P3-E: Optimize EM thread utilization (HIGH impact, HIGH effort)

**Problem**: The batch C++ EM runs 29,266 loci across 4 threads in 31.3 s.
The largest locus has 59,930 units × 297 transcripts — a single large locus
can dominate wall time (long-pole effect). Current OpenMP `schedule(dynamic)`
might not balance well when a few loci are orders of magnitude larger.

**Proposal** (investigate first):
1. **Profile per-locus EM time distribution** — add timing output for top-N
   loci to confirm the long-pole hypothesis.
2. If confirmed, **sort loci largest-first** so large loci start early and
   small loci fill in the gaps.
3. Consider **sub-locus parallelism** for the handful of very large loci
   (the EM inner loop over candidates is embarrassingly parallel per unit).
4. **SIMD vectorization** of the EM inner loop (log-sum-exp, posterior
   computation) using NEON intrinsics.

**Expected gain**: −5–15 s (bring EM from 31 → 16–26 s), depending on
load-balancing improvement.

### P3-F: Port multimapper handling to C++ (MEDIUM impact, MEDIUM effort)

**Problem**: Multimapper scoring remains in Python (`_flush_mm_group`). While
this dataset has 0 multimappers, real datasets have 5-20% multimap rate. The
Python path uses `array.array` accumulators and Python loops.

**Proposal**: Extend `fused_score_buffer()` to handle multimapper molecules
natively. The C++ code already sees `num_hits > 1` fragments — instead of
skipping them to Python, score them in the same two-pass framework with
proper WTA gating.

**Expected gain**: 0 s on this dataset, but critical for real datasets with
multimappers (could be 1-3 s).

---

## Priority Ranking

| Priority | Item | Impact | Effort | Est. Gain |
|----------|------|--------|--------|-----------|
| 1 | **P3-E**: EM thread optimization | HIGH | HIGH | −5 to −15 s |
| 2 | **P3-B**: `build_loci` → C++ | MED | LOW | −1.5 s |
| 3 | **P3-A**: Skip buffer spill | MED | MED | −0.9 s, −500 MB |
| 4 | **P3-C**: `_build_locus_meta` → C++ | LOW | LOW | −1.3 s |
| 5 | **P3-D**: Exon iteration → C++ | LOW | LOW | −0.5 s |
| 6 | **P3-F**: MM → C++ | MED | MED | critical for real data |

**Combined best-case**: −9.7 s to −19.7 s (bringing wall time from 131 s to
112–121 s, throughput 198–214 K frags/s). EM is the single highest-leverage
target for runtime.

For **memory**, the 16 GB target is already met (13.9 GB peak). Further
reductions are possible via:
- Buffer spill avoidance (P3-A): −500 MB
- TranscriptIndex CSR packing (deferred from Phase 2): −150–200 MB
- Streaming scoring without full chunk reload: −2 GB

---

## Summary

The fused score buffer (Phase 2) delivered the flagship improvement: **−25%
peak memory** (18.5 → 13.9 GB). The 16 GB target is met. Phase 3 shifts
focus to **runtime**, where `scan_and_buffer` (70%) and `locus_em` (24%)
dominate. The scan is largely I/O-bound and hard to parallelize without
architectural changes. The EM is the best runtime target: investigate
load-balancing and SIMD optimization for an estimated −5 to −15 s improvement.
Low-hanging fruit (P3-B through P3-D) can collectively save another −3.3 s
by porting remaining Python loops to C++.
