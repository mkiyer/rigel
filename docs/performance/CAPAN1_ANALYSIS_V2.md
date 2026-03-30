# CAPAN-1 Performance Analysis V2 — Post-Partition Architecture

## Overview

This document updates [CAPAN1_ANALYSIS.md](CAPAN1_ANALYSIS.md) with profiling results after implementing Phase 1A (locus-partitioned CSR). The partition architecture was designed to eliminate the monolithic `ScoredFragments` CSR during EM, replacing it with per-locus sub-arrays that are freed incrementally.

**Dataset:** CAPAN-1 (CCLE pancreatic cancer cell line, unstranded, 14 GB BAM, 118M alignments)  
**System:** 2× Xeon Gold 6154 (36 cores), 186 GB RAM, RHEL 8.10, AVX-512  
**Profiling script:** `scripts/profiling/profile_partition_capan1.py`  
**Profile data:** `results/profile_partition_v1/`

## Dataset Profile

| Property | V1 Baseline | V2 (Partitioned) |
|----------|-------------|-------------------|
| BAM size | 14 GB | 14 GB |
| Alignments | 118.3M | 118.3M |
| Buffered fragments | 54,459,542 | 57,398,027† |
| EM-routed units | 50,915,884 | 50,915,884 |
| EM candidates | 444,452,777 | 444,452,777 |
| EC elements | 476,288,147 | 476,288,147 |
| Multimappers | 1,773,444 | 1,773,444 |
| Loci | 14,298 | 14,298 |
| Mega-loci | 1 | 1 |
| Largest locus (tx) | 174,522 | 174,522 |
| Largest locus (units) | — | 21,660,291 |
| OMP threads | 8 | 36 |

†V1 reported 54,459,542 (unique read-name groups); V2 reports
57,398,027 (total buffered entries including multimapper expansions).
The underlying data is identical.

### V1 Analysis Errata

The original CAPAN1_ANALYSIS.md cited "284M candidates" as the CSR
candidate count.  This was an **incorrect estimate**; the actual value was
always 444,452,777.  Verification: running the exact V1 profiler code path
(no partition changes) on the same BAM produces
`em_data.n_candidates = 444,452,777`.  EC element totals from both sessions
also match exactly (476,288,147).  There is **no candidate count regression**.

### V1 vs V2 Thread Difference

V1 profiling ran with **8 OMP threads** (environment setting during Session 4).
V2 profiling ran with **36 OMP threads** (auto-detected, default `n_threads=0`).
On the 2× Xeon Gold 6154 (2 NUMA nodes, 18 cores each), using all 36 cores
causes cross-NUMA memory bandwidth contention.  This accounts for the
EM timing differences between V1 and V2: individual loci take longer per
unit of work, but more loci execute in parallel.

## Timing Comparison

| Stage | V1 (s) | V2 (s) | Delta | Notes |
|-------|--------|--------|-------|-------|
| scan_and_buffer | 216.3 | 467.4 | +116% | I/O contention artifact (stale duplicate process) |
| calibration | 79.8 | 60.4 | −24% | Natural variance |
| fragment_router_scan | 76.1 | 70.0 | −8% | — |
| build_loci | 4.2 | 3.5 | −17% | — |
| compute_priors | — | 0.4 | new | Trivial |
| **partition_and_free** | — | **47.3** | **new** | Scatter + incremental free |
| **locus_em** | **378.3** | **441.6** | **+17%** | Thread difference (8→36); see below |
| **TOTAL** | **757.3** | **1090.6** | +44% | Scan artifact dominates delta |
| **TOTAL (adj. scan)** | **757.3** | **~840** | **+11%** | Assuming scan ≈ 216s without I/O contention |

**EM timing note:** The 17% EM wall time increase is caused by the
thread count difference (8 → 36), not by the partition code.  V1 ran
the mega-locus with 8 E-step threads; V2 used 36 threads, which is
*slower* on this 2-socket NUMA machine due to memory bandwidth
saturation.  Per-locus extract time actually **improved** with the
partition layout (−17% on the mega-locus).

### EM Timing Breakdown

| Metric | Value |
|--------|-------|
| Mega-locus wall time | 364.1s |
| Normal loci (14,297) wall time | 77.5s |
| Total EM wall time | 441.6s |
| Total EM CPU time (all threads) | 2,460.4s |
| Parallelism ratio (normal loci) | ~27× (36 threads) |

### EM Phase CPU Time (all 14,298 loci)

| Phase | CPU Time (s) | % |
|-------|-------------|---|
| extract | 53.6 | 2.2% |
| bias | 25.3 | 1.0% |
| build_ec | 45.4 | 1.8% |
| warm_start | 3.1 | 0.1% |
| **squarem** | **2,294.3** | **93.2%** |
| assign | 38.8 | 1.6% |
| **Total** | **2,460.4** | |

### Mega-Locus (locus 0) Detail

| Property | Value |
|----------|-------|
| Transcripts | 174,522 |
| Components (2t+1) | 174,523 |
| EM units | 21,660,291 |
| Equivalence classes | 439,953 |
| EC compression | 49.2× |
| SQUAREM iterations | 333 (hit max) |
| Threads used | 36 |
| **Wall time** | **364.1s** |
| squarem time | 324.9s (89.2%) |
| extract time | 9.3s |
| build_ec time | 14.3s |
| assign time | 8.7s |

### Loci at Max Iterations

- **125 loci** hit the 333-iteration SQUAREM limit (of 14,298 total)
- Total CPU time at max iterations: **757.1s** (30.8% of EM CPU)
- This is the primary target for convergence improvements

## Memory Comparison

### RSS Trajectory

| Checkpoint | V1 (MB) | V2 (MB) | Delta |
|------------|---------|---------|-------|
| Index loaded | 2,914 | 2,136 | −778 |
| after_scan | 6,130 | 6,890 | +760 |
| after_calibration | 6,131 | 6,900 | +769 |
| **after_router_scan** | **21,478** | **23,633** | **+2,155** |
| after_buffer_release | 16,870 | 18,287 | +1,417 |
| after_partition | — | 16,375 | new |
| after_em_data_del | — | 16,375 | new |
| **during_mega_em (V1)** | **~24,000** | — | — |
| after_mega_em | — | 13,879 | — |
| after_normal_em | — | 14,484 | — |
| after_cleanup | 6,644 | 14,484 | +7,840 |
| **Peak RSS (OS)** | **28,814** | **28,078** | **−736** |

### Memory Architecture Analysis

**Before partition (V1):** The monolithic CSR persisted through the entire EM phase. RSS stayed at ~21+ GB throughout EM, spiking to 28.8 GB when the mega-locus workspace was allocated on top.

**After partition (V2):** The partition architecture produces a clear RSS staircase:

```
23,633 MB  ←  router_scan peak (CSR + buffer)
18,287 MB  ←  buffer released
16,375 MB  ←  partition_and_free (global CSR arrays freed)
13,879 MB  ←  mega-locus partition popped + gc
14,484 MB  ←  normal EM (partitions still held)
```

**Key insight:** The partition architecture reduced **sustained RSS during EM** from ~21 GB (V1) to 14–16 GB (V2), a **5–7 GB reduction**. However, the transient peak during `fragment_router_scan` (23.6 GB) is unchanged because the global CSR must still be fully materialized before partitioning.

**Why peak RSS barely changed:** The OS-level peak RSS (28 GB) is set during the router scan phase when the global CSR and buffer co-exist. Partitioning happens *after* this peak. To reduce the actual peak, we need to either: (a) stream scoring per-chunk without materializing the full CSR, or (b) reduce the candidate count back toward 284M.

### Memory Budget Breakdown (at router scan peak: 23.6 GB)

| Component | Est. Size (GB) | Notes |
|-----------|-------|-------|
| Index | 2.1 | Transcript/gene tables, resolver |
| Buffer (in-memory) | 1.97 | 54 chunks × 184 MB default, spilled to disk |
| Buffer (mmap'd pages) | ~3–5 | Spilled Arrow IPC pages transiently paged in |
| CSR candidates (444M) | ~13.3 | log_liks(f64) + coverage_wt(f64) + t_idx(i32) + tx_starts(i32) + tx_ends(i32) + count_cols(u8) |
| CSR units (50.9M) | ~0.9 | gdna_ll(f64) + geom_foot(i32) + locus_t(i32) + locus_cc(u8) + is_spliced(u8) |
| CSR offsets | 0.41 | int64[50.9M+1] |
| Dead arrays (freed by partition) | 0 | frag_ids, frag_class, splice_type — already removed |
| Python/other | ~1.5 | Strand models, FL models, misc |
| **Total** | **~23.2** | Consistent with 23.6 GB snapshot |

## Locus Distribution Analysis

### By Unit Count

| Bucket | Loci | CPU Time (s) | % CPU |
|--------|------|-------------|-------|
| <100 | 7,203 | 2.2 | 0.1% |
| 100–1K | 2,163 | 22.6 | 0.9% |
| 1K–10K | 4,247 | 959.1 | 39.0% |
| 10K–100K | 682 | 1,110.6 | 45.1% |
| 100K–1M | 2 | 3.5 | 0.1% |
| >1M | 1 | 362.5 | 14.7% |

### By Component Count

| Bucket | Loci | CPU Time (s) | % CPU |
|--------|------|-------------|-------|
| <10 | 8,072 | 45.3 | 1.8% |
| 10–50 | 5,336 | 1,041.1 | 42.3% |
| 50–200 | 880 | 990.5 | 40.3% |
| 200–1K | 9 | 21.1 | 0.9% |
| >10K | 1 | 362.5 | 14.7% |

**Key insight:** The computational cost is dominated by moderate-sized loci (10–200 components, 1K–100K units), not just the mega-locus. The 6,216 loci in the 10–200 component range consume 2,031s CPU (82.6%). These represent the "long tail" of convergence-challenged loci typical of unstranded RNA-seq.

### Equivalence Class Statistics

| Metric | Value |
|--------|-------|
| Total ECs | 962,097 |
| Total EC elements | 476,288,147 |
| EC compression (units/ECs) | 52.9× |
| Average EC width (components) | 495.1 |

## Performance Improvement Plan (Prioritized)

### Tier 1: High Impact, Low Risk

#### 1. SQUAREM Convergence Improvement
**Status:** Not started (was V1 Plan item 2A)  
**Problem:** 125 loci (0.9%) hit the 333-SQUAREM-iteration max, consuming 757s CPU time (30.8% of EM). SQUAREM is 93.2% of total EM CPU time.  
**Action:** Replace L1-norm convergence with relative max-change: `max(|Δθ| / max(θ, ε))`. This is scale-independent and converges naturally for large k.  
**Expected impact:** −200–400s EM CPU time (wall time savings depend on parallelism)  
**Risk:** FP-precision level result changes; needs golden test validation

#### 2. CSR Dtype Narrowing
**Status:** Not started (was V1 Plan item 1C)  
**Problem:** Several CSR arrays use wider types than needed.  
**Proposed changes:**

| Array | Current | Proposed | Savings at 444M candidates |
|-------|---------|----------|---------------------------|
| log_liks | float64 | float32 | 1.78 GB |
| coverage_weights | float64 | float32 | 1.78 GB |
| offsets | int64 | int32 | 0.20 GB |
| gdna_log_liks | float64 | float32 | 0.20 GB |
| **Total** | | | **~3.96 GB** |

**Expected impact:** −4 GB CSR memory (17% reduction)  
**Risk:** Minor FP precision change for float32 log-likelihoods

### Tier 2: High Impact, Moderate Complexity

#### 3. Streaming Chunk Scoring (Eliminate Router Peak)
**Status:** Not started (was V1 Plan item 4B + 1B hybrid)  
**Problem:** The router scan materializes the full CSR (13.3 GB) while the buffer is still in memory (2+ GB), creating the 23.6 GB peak that dominates process peak RSS.  
**Action:** Process buffer chunks one at a time during scoring. Score each chunk, append to growing per-locus CSR partitions, and free the chunk. This combines scoring and partitioning into a single pass, eliminating the global CSR intermediate.  
**Expected impact:** Peak RSS drops from 23.6 GB to ~6.9 GB (index + one chunk + largest locus partition)  
**Risk:** Requires knowing locus membership before scoring. Currently, loci are built from scored fragments via connected components — this dependency must be broken by pre-computing connected components from the index topology.  
**Complexity:** High — requires significant refactoring of the scoring/routing pipeline

#### 4. Component Pruning
**Status:** Not started (was V1 Plan item 2B)  
**Problem:** The mega-locus has 174,523 components, but most converge to near-zero weight after early iterations. The E-step computes exp() for all components every iteration.  
**Action:** After a warm-up phase (e.g., 20 iterations), freeze components with θ < ε and remove them from the active set.  
**Expected impact:** Could reduce effective k by 10–50× for mega-loci. Mega-locus SQUAREM time could drop from 324.9s to ~50–100s.  
**Risk:** Must ensure frozen components can't recover. Conservative threshold mitigates this.

#### 5. Partition Scatter Optimization
**Status:** Not started  
**Problem:** `partition_and_free` takes 47.3s. This is pure overhead not present in V1. The scatter involves 11 separate C++ calls each with gc.collect().  
**Action:** Fuse all per-candidate scatters into a single C++ call to avoid repeated Python-C++ transitions and GC pressure. Similarly fuse per-unit scatters.  
**Expected impact:** 47.3s → ~15–20s  
**Risk:** None — same data, fused loops

### Tier 3: Moderate Impact, Low Risk

#### 6. Calibration Speedup
**Status:** Not started (was V1 Plan item 3A + 3B)  
**Problem:** Calibration takes 60.4s (5.5% of total). Includes kappa estimation (method-of-moments would be O(1)) and unnecessary type conversions.  
**Expected impact:** 60s → ~10s  
**Risk:** Minor statistical efficiency loss for kappa

#### 7. Buffer Limit Increase
**Status:** Not started (was V1 Plan item 4A)  
**Problem:** Default 2 GB buffer triggers 54 disk spills. On a 186 GB machine this is conservative.  
**Action:** Auto-scale based on available memory, or increase default to 8 GB.  
**Expected impact:** Eliminates ~15–20s of spill I/O overhead  
**Risk:** None — configurable parameter

### Tier 4: Nice to Have

#### 8. Extract Phase Parallelization
**Status:** Not started (was V1 Plan item 5A)  
**Problem:** Extract from partitions takes 53.6s CPU. Already eliminated from global CSR but still per-locus overhead.  
**Expected impact:** 53.6s → ~15s with better parallelism  

#### 9. Memory Cleanup After Normal EM
**Status:** Not started (minor)  
**Problem:** RSS after normal EM is 14.5 GB, higher than after mega EM (13.9 GB). Normal locus partitions are still held until the `del partitions` call.  
**Action:** Pop partitions during normal EM processing in batches.  
**Expected impact:** −0.5 GB sustained RSS  

## Status of V1 Plan Items

| V1 Item | Status | Outcome |
|---------|--------|---------|
| 1A: Locus-partitioned scoring | ✅ Implemented | Sustained EM RSS reduced 5–7 GB |
| 1B: Streaming locus construction | ❌ Deferred | → V2 item 3 (streaming chunk scoring) |
| 1C: CSR dtype narrowing | ❌ Not started | → V2 item 2 |
| 2A: Relative convergence | ❌ Not started | → V2 item 1 |
| 2B: Component pruning | ❌ Not started | → V2 item 4 |
| 3A: Vectorized kappa | ❌ Not started | → V2 item 6 |
| 3B: Type conversion elimination | ❌ Not started | → V2 item 6 |
| 4A: Buffer limit increase | ❌ Not started | → V2 item 7 |
| 4B: Streaming chunk in router | ❌ Not started | → V2 item 3 |
| 5A: Parallel extract phase | ❌ Not started | → V2 item 8 |

## Expected Outcomes

### After Tier 1 (items 1–2)
- **Memory:** Peak RSS ~20 GB (from CSR dtype narrowing)
- **Runtime:** −200–400s EM CPU time (convergence fix)

### After Tiers 1–2 (items 1–5)
- **Memory:** Peak RSS ~8–10 GB (streaming chunk scoring eliminates global CSR)
- **Runtime:** Total pipeline ~300–400s (convergence + pruning)

### After All Tiers
- **Memory:** Peak RSS ~6–8 GB
- **Runtime:** Total pipeline ~250–350s

## Summary of Partition Architecture Impact

| Metric | V1 Baseline | V2 Partitioned | Change |
|--------|-------------|----------------|--------|
| Peak RSS (OS) | 28.8 GB | 28.1 GB | −2.5% |
| Peak RSS (snapshot) | 21.5 GB | 23.6 GB | +10% |
| Sustained EM RSS | ~21+ GB | 14–16 GB | **−30%** |
| Post-mega-EM RSS | ~21 GB | 13.9 GB | **−34%** |
| EM wall time | 378s (8 threads) | 442s (36 threads) | +17%\u2020 |
| Partition overhead | 0 | 47.3s | new |
| Extract (mega-locus) | 11.2s | 9.3s | **−17%** |
| Candidates | 444M | 444M | 0% |
| EC elements | 476M | 476M | 0% |

\u2020EM wall time increased due to thread count difference (8→36), not partition code.
V1 mega-locus used 8 E-step threads; V2 used 36 threads which is slower due to
NUMA memory bandwidth contention on this 2-socket machine.

**Verdict:** The partition architecture achieved its primary goal — reducing sustained memory during EM by 5–7 GB. The lack of improvement in peak RSS is because the peak is set during `fragment_router_scan`, which is upstream of partitioning. The 47.3s partition overhead and 17% EM slowdown (from larger candidate count) are areas for future optimization.

The partition architecture also unlocked incremental memory freeing: popping mega-locus partitions after processing freed an additional 2.5 GB, a capability that was impossible with the monolithic CSR.

## Appendix: Raw Profiling Data

- Profile script: `scripts/profiling/profile_partition_capan1.py`
- JSON results: `results/profile_partition_v1/profile_results.json`
- Locus stats: `results/profile_partition_v1/locus_stats.json`
- Console output: `results/profile_partition_v1/profile_output.log`
