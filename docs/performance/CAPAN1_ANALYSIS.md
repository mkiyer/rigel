# CAPAN-1 Performance Analysis & Improvement Plan

## Dataset Profile

| Property | Value |
|----------|-------|
| Sample | CCLE CAPAN-1 (pancreatic cancer cell line) |
| Library | Unstranded RNA-seq, high RNA purity, low gDNA |
| BAM file | 14 GB, 118.3M alignments |
| Buffered fragments | 54,459,542 |
| EM-routed units | 50,915,884 |
| Deterministic-unambig | 2,329,312 (4.4% of exonic) |
| Multimappers | 1,961,695 (3.6%) |
| Strand-ambiguous EM units | 48,747,533 (95.7% of EM units) |

## Timing Summary

| Stage | Time | % |
|-------|------|---|
| scan_and_buffer | 216.3s | 28.6% |
| calibration | 79.8s | 10.5% |
| fragment_router_scan | 76.1s | 10.1% |
| build_loci | 4.2s | 0.6% |
| **locus_em** | **378.3s** | **50.0%** |
| **TOTAL** | **757.3s** | |

## Memory Summary

| Checkpoint | RSS (MB) | Delta |
|------------|----------|-------|
| Index loaded | 2,914 | — |
| after_scan | 6,130 | +3,216 |
| after_calibration | 6,131 | +1 |
| **after_router_scan** | **21,478** | **+15,348** |
| after_buffer_release | 16,870 | −4,609 |
| after_locus_em | 20,441 | +3,572 |
| after_cleanup | 6,644 | −13,797 |
| **Peak RSS** | **28,814** | |

## Cross-Dataset Comparison

|  | Small (LBX0338) | Large (MI_1190) | **XL (CAPAN-1)** |
|--|-----------------|-----------------|-------------------|
| BAM size | 203 MB | 1.2 GB | **14 GB** |
| Buffered frags | 729K | 9.25M | **54.5M** |
| EM-routed units | ~400K | 7.4M | **50.9M** |
| Strand type | Stranded | Stranded | **Unstranded** |
| Det-unambig % | ~70% | ~55% | **4.4%** |
| Peak RSS | 2.9 GB | 5.8 GB | **28.8 GB** |
| Router memory | 328 MB | 1,612 MB | **15,348 MB** |
| Total time | 43s | 778s | **757s** |
| Loci at max iters | ~3 | 4 | **126** |

## Root Cause Analysis

### Why This Dataset Is Different

1. **Unstranded → strand ambiguity explosion.** In stranded RNA-seq, ~55-70% of fragments are deterministic-unambig (single spliced candidate on the correct strand). In unstranded data, only 4.4% qualify — fragments overlap transcripts on both strands, creating ≥2× more candidates per unit. This inflates the ScoredFragments CSR from 1.6 GB (stranded) to **15.3 GB** (unstranded).

2. **High depth → massive unit count.** 54.5M fragments (6× more than MI_1190) with 50.9M EM units generates enormous CSR arrays that must persist in memory throughout the entire EM phase.

3. **Convergence breakdown is widespread.** 126 loci hit the 333-iteration max (vs 4 for MI_1190). The unstranded ambiguity creates more complex loci where VBEM struggles to converge. These 126 loci consume 408.9s (48.5% of EM time).

### Memory Architecture Problem

The fundamental issue is the **monolithic ScoredFragments** structure. The current architecture:

```
Buffer (54.5M frags, 3.2 GB)
    ↓ fragment_router_scan (76s)
ScoredFragments (50.9M units, ~15 GB) ← stays in memory
    ↓ build_loci (4s)
Loci[14,298] (references into ScoredFragments)
    ↓ locus_em (378s, processes loci one-at-a-time)
    Per-locus: extract → build_ec → squarem → assign
```

The entire ScoredFragments must stay in memory because each locus references arbitrary ranges within the global CSR. During the mega-locus EM, the system holds:
- ScoredFragments: 15.3 GB
- Mega-locus extract copy: ~3 GB
- Mega-locus EC data: ~2.8 GB
- Base (index): ~2.9 GB
- **Peak: 28.8 GB**

## Performance Improvement Plan

### Phase 1: Memory Architecture (Critical — Enables Scaling)

#### 1A. Locus-Partitioned Scoring [HIGH IMPACT]

**Problem:** ScoredFragments is monolithic (15 GB for CAPAN-1).
**Solution:** After `build_loci` determines which units belong to which locus, partition the global CSR into per-locus sub-arrays and free the global structure.

**Approach:**
1. Run the router scan as today (produces global CSR)
2. Run `build_loci` to determine locus assignments
3. For each locus, gather its unit data into a compact per-locus structure
4. Free the global ScoredFragments
5. Process loci from their local copies

**Impact:** Peak memory drops from ~29 GB to ~6 GB + largest locus data (~3 GB) = ~9 GB. The partitioning itself is O(n_units) and would take ~5-10s.

**Risk:** None — data is the same, just reorganized.

#### 1B. Streaming Locus Construction [HIGHER IMPACT, MORE COMPLEX]

**Problem:** The global CSR is an intermediate — we only need per-locus data.
**Solution:** Merge the router_scan and build_loci steps: score fragments and immediately assign them to loci, building per-locus CSR incrementally. Never materialize the full global CSR.

**Approach:**
1. Pre-compute locus membership (from the transcript-to-locus mapping in the index)
2. During fragment scoring, route each fragment directly to its locus's CSR buffer
3. Each locus accumulates its own CSR arrays independently

**Impact:** Eliminates the 15 GB peak entirely. Memory = buffer + largest locus. Requires refactoring the scoring/routing pipeline.

**Risk:** Medium — requires ensuring locus membership can be determined before scoring (currently, loci are built FROM scored fragments via connected components).

#### 1C. CSR Dtype Optimization [MODERATE IMPACT, LOW RISK]

**Problem:** Several CSR arrays use wider types than necessary.
**Solution:** Use narrower types where range permits:

| Array | Current | Proposed | Savings (50.9M units / 444M cands†) |
|-------|---------|----------|-------------------------------------|
| log_liks | float64 | float32 | 1.14 GB |
| coverage_weights | float64 | float32 | 1.14 GB |
| offsets | int64 | int32 | 0.20 GB |
| frag_ids | int64 | int32 | 0.20 GB |
| gdna_log_liks | float64 | float32 | 0.20 GB |

**Total savings: ~2.9 GB** (19% reduction in CSR memory).

For `offsets`, int32 supports up to 2.1 billion candidates — sufficient since CAPAN-1 has 444M (†the 284M originally cited here was an incorrect estimate). For `frag_ids`, int32 supports 2.1B fragments — sufficient since CAPAN-1 has 54.5M. The `log_liks` and `coverage_weights` arrays are consumed in the extract phase where they're copied into the EM EC structures (which remain float64).

**Risk:** FP-precision level changes for log_liks/coverage_weights. None for integer narrowing.

### Phase 2: Convergence (High Impact on Wall Time)

#### 2A. Relative Convergence Criterion

**Problem:** 126 loci hit the 333-iteration max, consuming 409s (48.5% of EM time). The L1 convergence norm scales with the number of components, making convergence harder for larger loci.

**Solution:** Replace L1 sum with max relative change: `max(|Δθ| / max(θ, ε))`. This is scale-independent and converges more naturally for large k.

**Impact:** Most of those 126 loci would converge in 50-100 iterations instead of 333. Expected savings: ~250-300s.

**Risk:** FP-precision level result changes. Needs ground truth validation.

#### 2B. Component Pruning

**Problem:** The mega-locus has 174,523 components but most have negligible weight after early iterations. Computing exp() for all components each iteration wastes cycles. The early-zero skip helps per-row but doesn't reduce column-sum accumulation.

**Solution:** After a warm-up phase (e.g., 20 iterations), freeze components with θ < ε (e.g., 1e-10) and remove them from the active set. This reduces effective k, improving both compute and convergence.

**Impact:** Could reduce effective k by 10-50× for the mega-locus, dramatically speeding convergence.

**Risk:** Must ensure frozen components can be "thawed" if their weight would increase. Or use a conservative threshold.

### Phase 3: Calibration Optimization (10.5% of Total)

#### 3A. Vectorized Kappa Estimation

**Problem:** `estimate_kappa_marginal` consumes 28.3s (726 calls to `_marginal_loglik`, each calling `scipy.special.betaln`). 

**Solution:** Pre-compute betaln lookup tables or use analytical moment-based κ estimator. The method-of-moments approximation `κ̂ = p̄(1-p̄)/s² - 1` is O(1) vs O(n_regions × n_eval) for marginal likelihood.

**Impact:** 28.3s → <1s.

**Risk:** Some loss of statistical efficiency. Validate on simulation data.

#### 3B. Type Conversion Elimination

**Problem:** `np.asarray` (9.6s, 110 calls) and `astype` (7.1s, 812 calls) are pure overhead from dtype mismatches, consuming 16.7s total.

**Solution:** Ensure arrays originate with correct dtypes. Trace each call site and fix the source. Many may be from `int64 → intp`, `float32 → float64`, or contiguity checks that could be avoided.

**Impact:** 16.7s → ~1s.

**Risk:** None — pure optimization.

### Phase 4: Buffer Spilling Optimization

#### 4A. Increase Buffer Memory Limit

**Problem:** The 2 GB buffer limit triggers 47 spills + 48 reloads, adding ~14s during scan and ~3s during reload. The spill files also consume temporary disk space.

**Solution:** Allow a higher default buffer limit (e.g., 6-8 GB) or make it auto-scale based on available system memory. On the 186 GB system, spilling at 2 GB is unnecessarily conservative.

**Impact:** Eliminates spill I/O overhead (~17s). Increases scan-phase RSS by ~4 GB but still well within limits.

**Risk:** None — configurable parameter.

#### 4B. Streaming Chunk Processing in Router

**Problem:** `scan.py` loads ALL chunks into `chunk_arrays` before calling `fused_score_buffer`, causing buffer + chunk_arrays coexistence (additional ~2.7 GB).

**Solution:** Process chunks one at a time in `fused_score_buffer` — score each chunk, append to growing CSR, free the chunk.

**Impact:** −2.7 GB peak memory during router phase.

**Risk:** None if the C++ API supports incremental accumulation.

### Phase 5: Extract Phase Optimization

#### 5A. Vectorized Extract for Small Loci

**Problem:** The extract phase takes 65.4s across all loci, with 54.2s on 14,297 non-mega loci. These are small loci (mostly <10K units) processed sequentially in C++. The 1K-100K unit bucket (4,929 loci) dominates at 52s.

**Solution:** Batch small loci together and process their extractions in parallel. CurrentlyOpenMP parallelism is used for the E-step kernel but not for extract/build_ec.

**Impact:** 65.4s → ~15-25s with 8-thread parallelism.

**Risk:** None — read-only extraction.

## Priority Summary

| # | Optimization | Time Saved | Memory Saved | Risk | Complexity |
|---|-------------|-----------|-------------|------|-----------|
| 1A | Locus-partitioned scoring | — | **−12 GB** | None | Medium |
| 1C | CSR dtype narrowing | — | **−2.9 GB** | FP-precision | Low |
| 2A | Relative convergence | **~300s** | — | FP-precision | Low |
| 3B | Type conversion elimination | **~16s** | — | None | Low |
| 3A | Vectorized kappa estimation | **~27s** | — | FP-precision | Medium |
| 4A | Increase buffer limit | **~17s** | −4 GB (tradeoff) | None | Trivial |
| 5A | Parallel extract phase | **~40s** | — | None | Medium |
| 2B | Component pruning | **~100-200s** | **−1 GB** | FP-precision | Medium-High |
| 1B | Streaming locus construction | — | **−15 GB** | None | High |
| 4B | Streaming chunk in router | — | **−2.7 GB** | None | Medium |

### Recommended Execution Order

**Immediate (low-risk, high-impact):**
1. 1A — Locus-partitioned scoring (fixes OOM)
2. 1C — CSR dtype narrowing (easy memory win)
3. 4A — Increase buffer limit (trivial config change)
4. 3B — Type conversion elimination (pure cleanup)

**Next (FP-precision changes, needs validation):**
5. 2A — Relative convergence criterion
6. 3A — Vectorized kappa estimation

**Later (medium complexity):**
7. 5A — Parallel extract phase
8. 2B — Component pruning
9. 4B — Streaming chunk in router
10. 1B — Streaming locus construction (if 1A isn't sufficient)

## Expected Outcomes

**After Phase 1 (memory):** Peak RSS drops from 28.8 GB to ~12-15 GB. OOM eliminated for CAPAN-1 class datasets on 32 GB machines.

**After Phases 1-3:** Wall time drops from 757s to ~350-400s. Memory stays at ~12-15 GB.

**After all phases:** Wall time ~250-350s, peak RSS ~8-12 GB.

## Appendix: Profiling Data

- Profile output: `results/v2_profile_capan1/`
- Locus stats: `results/v2_profile_capan1/locus_stats_default.json` (14,298 loci)
- System: 2× Xeon Gold 6154, 186 GB RAM, 8 threads used
- Python 3.12.13, GCC 8.5.0, AVX-512 enabled
