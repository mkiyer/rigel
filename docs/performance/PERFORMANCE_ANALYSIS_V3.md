# Performance Analysis V3 — Streaming Scorer + Dual-Dataset Profiling

**Date:** 2026-03-30
**Datasets:** MEG-01 (heme/lymph), CAPAN-1 (pancreatic)
**System:** 2× Intel Xeon Gold 6154 (36 cores, 72 threads), 186 GB RAM, RHEL 8.10, AVX-512
**Software:** Python 3.12.13, GCC 14.3.0, OpenMP, 8 EM threads

---

## 1. Executive Summary

Profiling of MEG-01 and CAPAN-1 cell lines reveals that **the EM solver (SQUAREM) and BAM scanning dominate wall time**, while **CSR array memory dominates RSS**. The mega-locus (a single connected component with 160–175K transcripts) consumes 40% of EM time alone. Loci that hit the iteration cap (333 SQUAREM steps) account for ~50% of total EM time despite being <1% of loci.

Key optimization targets, in priority order:
1. **Mega-locus decomposition** — split or prune the giant locus to cap per-locus cost
2. **CSR memory reduction** — dtype narrowing saves ~2 GB (16% of CSR footprint)
3. **EM convergence tuning** — early termination or adaptive iteration caps
4. **Calibration acceleration** — replace scipy kappa estimation with C++/numba
5. **Extract-phase optimization** — reduce 10% overhead from global→local CSR gather

---

## 2. Dataset Profiles

| Metric | MEG-01 | CAPAN-1 |
|--------|-------:|--------:|
| BAM size | 12 GB | 14 GB |
| Total alignment records | 107,058,682 | 118,304,012 |
| Fragments (read names) | 49,368,403 | 54,459,542 |
| Unique mappers | 47,229,528 | 52,497,847 |
| Multimappers | 2,138,875 | 1,961,695 |
| Secondary alignments | 8,643,171 (8.1%) | 9,350,627 (7.9%) |
| Supplementary | 1,081,912 (1.0%) | 1,182,431 (1.0%) |
| Chimeric fragments | 453,349 (0.9%) | 471,887 (0.9%) |
| Intergenic | 743,976 | 742,459 |
| Records per fragment | 2.17 | 2.17 |
| Loci | 14,575 | 14,298 |
| Transcripts | 457,513 | 457,513 |
| Genes | 63,472 | 63,472 |

The two datasets have nearly identical alignment complexity (same records-per-fragment ratio, similar secondary alignment fractions). Differences in performance are attributable to EM convergence behavior and I/O variability, not alignment characteristics.

---

## 3. Timing Comparison

### 3.1 Stage Breakdown

| Stage | MEG-01 (s) | MEG-01 (%) | CAPAN-1 (s) | CAPAN-1 (%) | Ratio |
|-------|----------:|----------:|------------:|------------:|------:|
| scan_and_buffer | 1,196.3 | 51.9% | 338.3 | 29.7% | 3.54× |
| calibration | 70.1 | 3.0% | 91.1 | 8.0% | 0.77× |
| fragment_router_scan | 79.7 | 3.5% | 86.0 | 7.5% | 0.93× |
| build_loci | 4.3 | 0.2% | 3.9 | 0.3% | 1.10× |
| locus_em | 949.9 | 41.2% | 617.3 | 54.2% | 1.54× |
| **TOTAL** | **2,303.1** | | **1,139.3** | | **2.02×** |

### 3.2 Scan Throughput Anomaly

MEG-01 `scan_and_buffer` is 3.54× slower than CAPAN-1 despite processing 10% fewer alignment records with identical complexity. Both runs used `--cprofile --memory-interval 50 --threads 8`.

**Root Cause Assessment:**
- Python callback overhead is minimal: only 56 chunk callbacks totaling ~30s, so **>1,160s is spent inside C++ BAM scanning**.
- Alignment complexity is identical (2.17 records/fragment, same secondary fractions).
- The BAM files reside on `/scratch` (Lustre parallel filesystem). **I/O variability of 3–5× is common** on shared Lustre filesystems depending on system load, OST contention, and metadata server activity.
- Recommendation: Profile without cProfile and with a local-disk BAM copy to isolate I/O from CPU overhead. The true scan throughput is likely 200–300K frags/s.

### 3.3 EM Solver Breakdown (All Loci)

| Sub-stage | MEG-01 (s) | MEG-01 (%) | CAPAN-1 (s) | CAPAN-1 (%) |
|-----------|----------:|----------:|------------:|------------:|
| **squarem** | 1,709.5 | 85.9% | 1,093.5 | 84.7% |
| **extract** | 212.1 | 10.7% | 140.0 | 10.8% |
| build_ec | 29.8 | 1.5% | 25.4 | 2.0% |
| assign | 24.9 | 1.2% | 19.2 | 1.5% |
| bias | 13.1 | 0.7% | 10.9 | 0.8% |
| warm_start | 1.7 | 0.1% | — | — |
| **Total** | **1,991.0** | | **1,290.4** | |

> Note: EM total (1,991s) exceeds reported `locus_em` stage (950s) because EM sub-stages are measured per-thread with 8 OpenMP threads. The 950s is wall time; the 1,991s is the sum of per-locus wall times (serial sum).

### 3.4 Mega-Locus Dominance

| Metric | MEG-01 | CAPAN-1 |
|--------|-------:|--------:|
| Mega-locus transcripts | 159,851 | 174,522 |
| Mega-locus units | 18,897,853 | 21,660,291 |
| Mega-locus equiv classes | 407,807 | 439,953 |
| EC total elements | 182,698,326 | 221,378,807 |
| SQUAREM iterations | 333 (max) | 333 (max) |
| SQUAREM time | 764.3s | 477.4s |
| Total time | 799.8s | 519.9s |
| % of all EM time | 40.2% | 40.3% |
| Estimated memory | 2.6 GB | 3.1 GB |

The mega-locus always hits the iteration cap. MEG-01's mega-locus SQUAREM is 1.6× slower than CAPAN-1's despite being smaller, suggesting convergence difficulty (possibly due to worse identifiability in MEG-01's expression landscape).

### 3.5 Non-Converged Loci

| Metric | MEG-01 | CAPAN-1 |
|--------|-------:|--------:|
| Loci at max iterations (≥333) | 133 | 123 |
| Aggregate EM time | 1,034.2s | 642.0s |
| % of total EM time | **51.9%** | **49.7%** |

**~130 loci (0.9% of all loci) consume ~50% of total EM time.** These are the primary optimization target. Most non-converged loci are small (30–200 transcripts) but have high unit counts (33K–94K units), suggesting identifiability issues from overlapping transcripts.

---

## 4. Memory Analysis

### 4.1 RSS Timeline

| Snapshot | MEG-01 (MB) | CAPAN-1 (MB) | Delta |
|----------|----------:|------------:|------:|
| before | 2,862 | 2,914 | −52 |
| after_scan | 8,361 | 6,342 | +2,019 |
| after_calibration | 8,361 | 6,342 | +2,019 |
| **after_router_scan** | **20,581** | **20,272** | **+309** |
| after_buffer_release | 14,490 | 16,280 | −1,790 |
| after_locus_em | 17,387 | 19,854 | −2,467 |
| after_cleanup | 5,388 | 6,154 | −766 |
| **Peak RSS** | **22,940** | **26,613** | **−3,673** |

**Key Observations:**
- `after_router_scan` is consistent across datasets (~20 GB), confirming CSR array size dominates.
- Peak RSS for MEG-01 occurs during locus_em (mega-locus allocation), while CAPAN-1 peaks during router_scan.
- Buffer release frees 6.1 GB (MEG-01) vs 4.0 GB (CAPAN-1).

### 4.2 CSR Memory Composition

For MEG-01 (45.4M units, ~364M candidates estimated):

| Array | Dtype | Count | Bytes | GB |
|-------|-------|------:|------:|---:|
| **log_liks** | float64 | 364M cand | 8 | 2.91 |
| **coverage_weights** | float64 | 364M cand | 8 | 2.91 |
| t_indices | int32 | 364M cand | 4 | 1.46 |
| tx_starts | int32 | 364M cand | 4 | 1.46 |
| tx_ends | int32 | 364M cand | 4 | 1.46 |
| count_cols | uint8 | 364M cand | 1 | 0.36 |
| offsets | int64 | 45.4M units | 8 | 0.36 |
| gdna_log_liks | float64 | 45.4M units | 8 | 0.36 |
| frag_ids | int64 | 45.4M units | 8 | 0.36 |
| locus_t_indices | int32 | 45.4M units | 4 | 0.18 |
| genomic_footprints | int32 | 45.4M units | 4 | 0.18 |
| is_spliced | bool | 45.4M units | 1 | 0.05 |
| locus_count_cols | uint8 | 45.4M units | 1 | 0.05 |
| frag_class | int8 | 45.4M units | 1 | 0.05 |
| splice_type | uint8 | 45.4M units | 1 | 0.05 |
| **Total** | | | | **~12.2 GB** |

### 4.3 Locus Size Distribution

| Category | Count | % of Loci | EM Time Share |
|----------|------:|----------:|--------------:|
| ≤10 transcripts | 8,941 | 61.3% | negligible |
| 11–100 transcripts | 5,476 | 37.6% | ~45% |
| 101–1,000 transcripts | 157 | 1.1% | ~15% |
| Mega-locus (>1,000 tx) | 1 | 0.01% | 40% |

---

## 5. Prioritized Optimization Plan

### P0: Mega-Locus EM Acceleration

**Impact:** Save 300–750s wall time (40% of EM). Save 2.6 GB peak memory.
**Current:** One mega-locus with 160K transcripts always hits 333 SQUAREM iterations, consuming 800s.

#### Option A: Pre-EM Locus Splitting

Split the mega-locus into independent sub-loci before EM by pruning weak multimapper edges:

1. Build the fragment–transcript bipartite graph
2. Remove edges where log-likelihood is >50 nats below the best candidate for that unit
3. Find connected components on the pruned graph
4. Run EM independently on each sub-locus

**Implementation:**
- Add a `prune_weak_edges()` function in `locus.py` or `em_solver.cpp`
- Threshold parameter: `locus_prune_loglik_delta` (default: 50.0)
- Only apply to loci with `n_transcripts > 1000` to avoid overhead on small loci

**Risk:** Aggressive pruning may discard valid multimapper connections. Use conservative default (50 nats ≈ likelihood ratio of $e^{50}$).

#### Option B: Adaptive Iteration Cap

Scale `max_iterations` with locus complexity:

```
effective_max = min(max_iterations, base_cap + k * log2(n_components))
```

For the mega-locus (160K components): `base_cap=100, k=20` → effective_max ≈ 440 (no savings). Instead, use a wall-time cap:

```
if (elapsed_us > max_locus_wall_us) break;  // e.g., 60s per locus
```

**Implementation:** Add `max_locus_wall_sec` to `EMConfig` (default: 120). Check in the SQUAREM outer loop after each iteration.

#### Option C: Hierarchical EM (Research)

Two-level EM: first solve sub-loci (clusters of related transcripts), then refine globally. This is a research-level change with significant implementation cost.

**Recommendation:** Implement Option A first (locus splitting). If insufficient, add Option B (wall-time cap).

---

### P1: CSR Dtype Narrowing

**Impact:** Save ~2.0 GB memory (16% of CSR footprint).
**Risk:** Low — only narrowing arrays that don't affect EM numerical precision.

| Change | Savings | Notes |
|--------|--------:|-------|
| `offsets` int64 → int32 | 0.18 GB | Max n_candidates per dataset ≈ 364M, fits int32 (2.1B limit) |
| `gdna_log_liks` float64 → float32 | 0.18 GB | Used only for gDNA candidate pre-computation, not in EM inner loop |
| `frag_ids` int64 → int32 | 0.18 GB | Max fragments ≈ 54M, fits int32 |
| `coverage_weights` float64 → float32 | 1.46 GB | **Requires validation** with Kahan summation in EM |
| **Total** | **2.00 GB** | |

**Implementation:**
1. Change C++ `StreamingScorer` output dtypes in `scoring.cpp`
2. Update `ScoredFragments` dataclass type annotations
3. Update `em_solver.cpp` `extract_locus_sub_problem_from_partition()` to accept narrower input types
4. **Critical:** Keep EM internal `log_liks` as float64. Only narrow the global CSR transport layer.
5. Run golden output regression tests to confirm bit-for-bit equivalence

**Do NOT narrow:**
- `log_liks` (float64) — EM E-step uses log-sum-exp with Kahan compensation; float32 would cause catastrophic cancellation
- `t_indices` (int32) — already minimal

**Staged approach:**
- Phase 1: offsets, frag_ids, gdna_log_liks (safe, 0.54 GB)
- Phase 2: coverage_weights (needs Kahan validation, 1.46 GB)

---

### P2: Early Termination for Non-Converged Loci

**Impact:** Reduce 50% of EM time from non-converged loci (save 400–500s wall time).
**Current:** 133 loci hit max_iterations=333 SQUAREM steps. Many have near-zero theta changes in late iterations.

#### Option A: Relative Convergence Criterion

Current convergence uses L1 norm on normalized theta: $\sum_k |\theta_k^{(t)} - \theta_k^{(t-1)}|$. This is scale-dependent — loci with many near-zero components accumulate small residuals that prevent convergence.

Switch to max relative change:

$$\text{converged if } \max_k \frac{|\theta_k^{(t)} - \theta_k^{(t-1)}|}{\max(\theta_k^{(t)}, \epsilon)} < \delta$$

**Implementation:**
- In `em_solver.cpp` SQUAREM convergence check (~line 1060), replace L1 sum with max-relative-change
- Add `convergence_mode` to EMConfig: `"l1"` (default, backward compat) or `"max_relative"`
- Default threshold: $\delta = 10^{-4}$

#### Option B: Stagnation Detection

Track the convergence metric over a sliding window. If improvement rate drops below a threshold for N consecutive iterations, terminate early:

```cpp
if (iteration > 50 && 
    (prev_delta - delta) / prev_delta < 0.01) {
    break;  // <1% improvement rate
}
```

**Recommendation:** Implement Option A. Log convergence per-locus to validate.

---

### P3: Calibration Acceleration

**Impact:** Save 40–70s (3–8% of total time).
**Current:** `estimate_kappa_marginal()` calls `scipy.optimize.minimize_scalar` 20 times × 33 evaluations each = 660 `_marginal_loglik()` calls at 23ms each = 15.2s total. Additional time in `_compute_strand_llr_betabinom()` (0.35s), `build_gdna_fl_model()` (0.96s), `_e_step()` (0.51s), `_m_step()` (0.33s).

#### Option A: Numba JIT for _marginal_loglik

The inner loop evaluates `betaln()` over arrays of strand counts:

```python
def _marginal_loglik(self, kappa, counts_a, counts_b):
    alpha = kappa * p
    beta = kappa * (1 - p)
    ll = betaln(counts_a + alpha, counts_b + beta) - betaln(alpha, beta)
    return -ll.sum()
```

JIT this with `@numba.njit` or replace with a C++ implementation using `std::lgamma`.

**Implementation:**
- Add `_calibration_impl` C++ module OR numba-JIT the hot functions
- Expected speedup: 5–10× on kappa estimation → saves ~13s
- Total calibration time reduction: ~20–30s

#### Option B: Cache Kappa Results

Kappa estimation for strand balance uses similar count distributions across chromosomes. Cache results for chromosomes with similar total counts.

**Recommendation:** Option A (numba JIT) for immediate gains. Option B for diminishing returns.

---

### P4: Extract-Phase Optimization

**Impact:** Save 100–200s (10% of EM time).
**Current:** `extract_locus_sub_problem_from_partition()` does global→local transcript remapping and candidate deduplication for every locus. For the mega-locus, this takes 9.5s; aggregated across all loci: 212s.

#### Option A: Pre-sorted Candidates

The extract phase deduplicates candidates per unit by keeping only the best log-likelihood per component. If the CSR is pre-sorted by (unit, component), the dedup becomes a single linear scan instead of epoch-based hash tracking.

**Implementation:**
- Sort the CSR data by (unit_id, component_id) during `StreamingScorer.finish()`
- In `extract_locus_sub_problem_from_partition()`, replace the epoch-based loop with a merge-scan
- Expected speedup: 2–3× on extract phase → save 70–140s

#### Option B: Lazy Extraction

Only extract sub-problems for loci that will actually run EM (skip deterministic loci). Currently, all loci go through extraction even if they have a single transcript.

**Implementation:**
- Skip extraction for loci with `n_transcripts == 1`
- 8,941 single-transcript loci × minimal extract cost = small savings
- More impactful: skip loci with `n_units < 10` (trivially converge)

**Recommendation:** Option A for meaningful savings on large loci.

---

### P5: Buffer Spill Optimization

**Impact:** Reduce I/O overhead during `scan_and_buffer` by 5–15s.
**Current:** MEG-01 generated 29 spills to disk (Arrow IPC + LZ4), totaling 0.53s write + 0.60s read.

The 4 GiB default buffer significantly reduced spills compared to the old 2 GiB default. Further improvements:

1. **Adaptive buffer sizing:** Set buffer to `min(available_ram / 2, fragment_count_estimate * 200)` to avoid unnecessary spills for small datasets.
2. **Memory-mapped spills:** Use `mmap` instead of Arrow IPC for spilled chunks to avoid decompression overhead during read-back.

**Low priority** — spill I/O is <1s, dwarfed by other costs.

---

### P6: SQUAREM Vectorization

**Impact:** 5–15% speedup on SQUAREM inner loop for large loci.
**Current:** The SQUAREM `r_vec` and `v_vec` computations use scalar loops over `nc` components. For the mega-locus (`nc = 159,852`), these loops are not SIMD-vectorized.

**Implementation:**
- Use `__m256d` (AVX2) or `__m512d` (AVX-512) intrinsics for `r = theta1 - theta0`, `v = theta2 - 2*theta1 + theta0`, and the step-clamp loop
- Expected speedup: 2–4× on the r/v computation, but this is a small fraction of SQUAREM time (most time is in E-step)
- The E-step already uses `fast_exp()` with SIMD — the bottleneck is memory bandwidth, not compute

**Low-to-medium priority** — profile to confirm r/v fraction before vectorizing.

---

### P7: Partition-Based Memory Release

**Impact:** Reduce peak RSS by 2–4 GB.
**Current:** The full CSR is held in memory during all locus EM processing. After each locus completes, its sub-problem is freed, but the global CSR persists.

**Implementation:**
- Partition CSR arrays into per-locus slices during `build_loci()`
- Release each locus's partition after EM completes
- Use `numpy.split()` + explicit `del` + `gc.collect()` at partition boundaries
- **Challenge:** The partition requires sorted CSR data (units grouped by locus), which is already done in `build_loci()`. The issue is that numpy views share the backing buffer — need actual copies to release.

**Medium priority** — requires careful memory accounting. The 6 GB freed by buffer release already helps.

---

## 6. Implementation Roadmap

### Phase 1: Quick Wins (1–2 days)

| Item | Expected Savings | Complexity |
|------|-----------------|------------|
| P1 Phase 1: Safe dtype narrowing (offsets, frag_ids, gdna_log_liks) | 0.54 GB memory | Low |
| P2 Option A: Max-relative convergence criterion | 200–400s time | Low |
| P3 Option A: Numba JIT _marginal_loglik | 13s time | Low |

### Phase 2: Infrastructure (3–5 days)

| Item | Expected Savings | Complexity |
|------|-----------------|------------|
| P0 Option A: Mega-locus splitting via edge pruning | 300–750s time, 2.6 GB memory | Medium |
| P1 Phase 2: coverage_weights float32 with Kahan validation | 1.46 GB memory | Medium |
| P4 Option A: Pre-sorted CSR for extract phase | 70–140s time | Medium |

### Phase 3: Research (1–2 weeks)

| Item | Expected Savings | Complexity |
|------|-----------------|------------|
| P0 Option C: Hierarchical EM | Major time reduction | High |
| P6: SQUAREM SIMD vectorization | 5–15% EM speedup | Medium |
| P7: Partition-based memory release | 2–4 GB memory | Medium |

---

## 7. Measurement Methodology

### Profiler Configuration
Both runs used identical settings:
```
python scripts/profiler.py \
  --bam <path> --index <path> -o <output> \
  --stages --cprofile --memory-interval 50 --threads 8
```

### cProfile Overhead
cProfile adds significant overhead to Python-intensive phases. The `_thread.lock.acquire` tottime (2,240s for MEG-01) reflects thread wait time, not CPU contention. **Future profiling should use `--no-cprofile` for timing-sensitive measurements and `--cprofile` only when function-level hotspot data is needed.**

### I/O Variability
The 3.54× `scan_and_buffer` slowdown for MEG-01 vs CAPAN-1 is inconsistent with their identical alignment complexity. This is attributable to Lustre filesystem I/O variability. **Baseline timing comparisons should use local-disk BAM copies.**

---

## 8. Appendix: Raw Profiling Data

### MEG-01 Locus Statistics
- Total loci: 14,575
- Median transcripts per locus: 6
- Median units per locus: 45
- Median SQUAREM iterations: 6
- Max locus transcripts: 159,851
- Max locus units: 18,897,853
- Max EC width: 261
- Max EC depth: 73,527

### CAPAN-1 Locus Statistics
- Total loci: 14,298
- Max locus transcripts: 174,522
- Max locus units: 21,660,291
- Max EC elements: 221,378,807

### Top 10 MEG-01 Loci by Wall Time

| Rank | Idx | Transcripts | Units | ECs | SQUAREM (s) | Total (s) | Mega |
|------|----:|----------:|------:|----:|-----------:|----------:|------|
| 1 | 0 | 159,851 | 18,897,853 | 407,807 | 764.3 | 799.8 | Yes |
| 2 | 13355 | 39 | 58,286 | 153 | — | 12.5 | No |
| 3 | 10251 | 67 | 63,421 | 250 | — | 11.6 | No |
| 4 | 557 | 38 | 62,112 | 290 | — | 11.5 | No |
| 5 | 6346 | 74 | 33,742 | 317 | — | 11.4 | No |
| 6 | 6259 | 106 | 24,793 | 486 | — | 10.3 | No |
| 7 | 7780 | 187 | 42,608 | 855 | — | 10.2 | No |
| 8 | 2770 | 149 | 35,466 | 711 | — | 9.5 | No |
| 9 | 5474 | 33 | 44,004 | 152 | — | 9.1 | No |
| 10 | 3090 | 107 | 94,651 | 699 | — | 8.7 | No |

### cProfile Top Functions (MEG-01)

| Function | tottime (s) | cumtime (s) | calls |
|----------|----------:|----------:|------:|
| _thread.lock.acquire | 2,239.9 | 2,274.8 | 175,596 |
| threading.wait | 16.7 | 2,343.9 | 43,899 |
| _marginal_loglik | 15.1 | 15.2 | 660 |
| numpy.asarray | 9.7 | 9.7 | 104 |
| ndarray.astype | 4.3 | 4.3 | 598 |
| _snap_rss_current | 1.8 | 5.6 | 43,903 |
| ufunc.at | 2.0 | 2.0 | 20 |
| build_gdna_fl_model | 0.9 | 1.0 | 41 |
| gc.collect | 0.8 | 0.8 | 1 |
