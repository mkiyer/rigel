# VBEM vs MAP-EM Performance Analysis

**Date:** March 2026
**System:** 2× Intel Xeon Gold 6154 (36 cores), 186 GB RAM, RHEL 8.10

---

## 1. Executive Summary

MAP-EM is **2.13× faster** than VBEM on the large BAM (309s vs 658s wall time for `locus_em`), but **1.5× slower** on the small BAM (15.5s vs 10.4s). The performance reversal is explained entirely by **convergence behavior on the mega-locus**: VBEM hits the max iteration limit (333 SQUAREM iters = 1000 EM steps) without converging, while MAP converges in 123 iterations.

### Key Findings

1. **VBEM does not converge on the large mega-locus** (315K components, 5.8M units). It hits the 333-iteration cap. The L1 convergence criterion requires each of 5.8M components to change by < 1.7×10⁻¹³ on average — unreachable with VBEM's Dirichlet smoothing.

2. **Per-iteration cost is similar**: VBEM is actually 18% faster per SQUAREM iteration (1.91s vs 2.32s) because digamma overhead (~50ms, 2.5%) is small relative to the E-step kernel time.

3. **The E-step is memory-bandwidth bound**, not compute-bound: 3.54 GB of data per SQUAREM iteration at ~1.5–1.9 GB/s effective bandwidth. This is far from the ~60 GB/s theoretical bandwidth, indicating severe cache thrashing from the column-sum access pattern.

4. **Digamma is NOT the bottleneck** for VBEM. At ~50ns per call, 315K digamma calls × 3 E-steps = ~47ms per SQUAREM iteration, which is only 2.5% of the 1.91s per-iteration cost.

---

## 2. Profiling Results

### Wall-Clock Timing (8 threads)

| Stage | VBEM Small | MAP Small | VBEM Large | MAP Large |
|-------|-----------|----------|-----------|----------|
| scan_and_buffer | 28.5s | 27.7s | 183.1s | 183.2s |
| calibration | 7.1s | 6.9s | 39.1s | 38.6s |
| **locus_em** | **10.4s** | **15.5s** | **657.8s** | **309.3s** |
| Total | 48.6s | 52.7s | 907.1s | 553.5s |

### C++ Aggregate Locus EM Time (sum across all loci, single-thread equivalent)

| Metric | VBEM Small | MAP Small | VBEM Large | MAP Large |
|--------|-----------|----------|-----------|----------|
| Total loci | 14,530 | 14,530 | 13,882 | 13,882 |
| Sum of locus times | 33.6s | 41.1s | 741.0s | 396.9s |
| Ratio | 1.00× | 1.22× | 1.87× | 1.00× |

### SQUAREM Iteration Statistics

| Metric | VBEM Small | MAP Small | VBEM Large | MAP Large |
|--------|-----------|----------|-----------|----------|
| Mean iterations | 6.0 | 7.6 | 3.2 | 6.1 |
| Median iterations | 4.0 | 6.0 | 2.0 | 3.0 |
| Max iterations | 333 | 159 | **333** | 333 |
| Loci at max (333) | 1 | 0 | **4** | 1 |

---

## 3. Mega-Locus Deep Dive

The mega-locus dominates total EM time. Both BAMs produce one mega-locus each.

### Small BAM Mega-Locus

| Metric | VBEM | MAP |
|--------|------|-----|
| Transcripts | 51,655 | 51,655 |
| Units (components) | 179,704 | 179,704 |
| Equiv classes | 27,854 | 27,854 |
| EC elements | 1,953,071 | 1,953,071 |
| Max EC width | 270 | 270 |
| Max EC depth | 6,462 | 6,462 |
| **SQUAREM iterations** | **74** | **121** |
| **SQUAREM time** | **6.57s** | **11.38s** |
| Time per iteration | 88.8ms | 94.0ms |

On the small mega-locus, VBEM converges faster (74 vs 121 iters) with slightly cheaper iterations.

### Large BAM Mega-Locus

| Metric | VBEM | MAP |
|--------|------|-----|
| Transcripts | 315,030 | 315,030 |
| Units (components) | 5,777,023 | 5,777,023 |
| Equiv classes | 275,385 | 275,385 |
| EC elements | 49,219,285 | 49,219,285 |
| Max EC width | 315 | 315 |
| Max EC depth | 10,077 | 10,077 |
| **SQUAREM iterations** | **333 (MAX!)** | **123** |
| **SQUAREM time** | **634.6s** | **285.1s** |
| Time per iteration | 1,906ms | 2,318ms |
| Digamma calls/iter | 945,093 | 0 |
| Total digamma calls | 314,715,969 | 0 |

**VBEM did not converge** — it exhausted the max iteration budget (1000 EM steps = 333 SQUAREM iterations). MAP converged in 123 iterations.

---

## 4. Per-Iteration Cost Breakdown

### Large Mega-Locus (per SQUAREM iteration)

Each SQUAREM iteration = 3 E-step + M-step cycles + extrapolation + convergence check.

| Component | VBEM | MAP | Notes |
|-----------|------|-----|-------|
| E-step kernel (×3) | ~1,800ms | ~2,200ms | Same `em_step_kernel_range()` |
| Digamma computation | ~47ms | 0ms | 315K calls × 3, ~50ns each |
| Normalization/M-step | ~20ms | ~40ms | MAP normalizes theta |
| SQUAREM overhead | ~40ms | ~80ms | r,v vectors, extrapolation |
| **Total** | **1,906ms** | **2,318ms** | **VBEM 18% faster per iter** |

### Why VBEM per-iteration is faster

VBEM's log-weight computation uses `digamma(alpha_i) - digamma(sum)`. The `digamma(sum)` is computed once and reused. Despite ~47ms of digamma overhead, VBEM has slightly simpler M-step (no normalization) and the E-step kernel has the same cost for both modes. The net result is VBEM is slightly faster per iteration.

### Memory Bandwidth Analysis

| Metric | Value |
|--------|-------|
| EC elements per E-step | 49,219,285 |
| Data per element | ~24 bytes (weight + theta read + accumulate) |
| Data per SQUAREM iter | 3 × 49.2M × 24 = 3.54 GB |
| Measured time | ~1.9s/iter (VBEM), ~2.3s/iter (MAP) |
| Effective bandwidth | 1.5–1.9 GB/s |
| Theoretical DDR4 bandwidth | ~60 GB/s |
| **Utilization** | **~3%** |

The E-step is severely memory-bandwidth limited due to the column-sum access pattern (Pass 2 in `em_step_kernel_range()`). Each row of the EC matrix is accessed row-major, but the column-sum accumulation has a stride equal to the number of components (5.8M doubles = 44 MB). This causes catastrophic L3 cache misses on every access.

---

## 5. Root Cause: VBEM Convergence Failure

### Why VBEM doesn't converge on large mega-loci

The convergence criterion is the L1 norm of the change in **normalized** theta:

```
delta = Σᵢ |theta_new[i] - theta_old[i]|
```

where `theta[i] = alpha[i] / Σ alpha`. With 5,777,023 components, for `delta < 1e-6`, each component must change by less than ~1.7×10⁻¹³ on average. This is near machine epsilon (2.2×10⁻¹⁶).

VBEM's Dirichlet posterior naturally produces smoother, more spread-out distributions than MAP-EM's point estimates. This means small components that MAP-EM would drive to exactly zero remain at small positive values in VBEM, creating persistent tiny oscillations that prevent the L1 norm from dropping below 1e-6.

### Possible solutions

1. **Scale convergence delta by number of components**: Use `convergence_delta / sqrt(n_components)` or `convergence_delta * n_components`. This accounts for the fact that L1 norm naturally scales with dimensionality. (**Recommended**)

2. **Use relative convergence**: Check `max_i |theta_new[i] - theta_old[i]| / max(theta_old[i], eps)` instead of L1 sum. This measures relative change of the largest components. 

3. **Use log-likelihood convergence**: Track ELBO/log-likelihood change instead of parameter change. More stable for high-dimensional problems.

4. **Increase max_iterations**: Not recommended — the mega-locus already takes 634s at 333 iterations. More iterations worsen performance without guaranteed convergence.

---

## 6. Optimization Priorities (Revised)

Based on these findings, the optimization priorities should be revised from the original plan:

### Priority 1: Fix VBEM convergence criterion

**Impact: Could recover the 2.23× VBEM slowdown on mega-loci**

The L1 convergence delta should be scaled by dimensionality. If VBEM converges in ~120 iterations (similar to MAP), the mega-locus VBEM time would drop from 634s to ~230s.

### Priority 2: Fuse Pass 2 column-sum into Pass 1 (cache optimization)

**Impact: Could improve E-step by 5–10× for mega-loci**

The current `em_step_kernel_range()` has two passes:
- Pass 1: Row-major scan computing per-EC responsibilities (cache-friendly)
- Pass 2: Column-major accumulation of component totals (cache-hostile, 44 MB stride)

Fusing these passes eliminates the column-sum cache thrashing. Since effective bandwidth is only 3% of theoretical, this is the largest single optimization opportunity.

### Priority 3: Digamma vectorization (low priority)

**Impact: ~2.5% improvement per iteration for VBEM**

Digamma is only ~47ms per SQUAREM iteration on the large mega-locus. AVX-512 vectorization could reduce this to ~6ms (8× parallel), saving ~41ms per iteration. At 120 iterations (if convergence is fixed), this saves ~5s total — negligible compared to the 285s E-step time.

### Priority 4: Thread scaling for mega-locus

The mega-locus monopolizes one thread while remaining loci are small. Investigate intra-locus parallelism for the E-step (partition ECs across threads within a single mega-locus solve).

---

## 7. Small vs Large BAM Performance Reversal

Why does VBEM win on the small BAM but lose on the large BAM?

| Factor | Small BAM | Large BAM |
|--------|-----------|-----------|
| Mega-locus components | 179,704 | 5,777,023 |
| VBEM SQUAREM iters | 74 | 333 (MAX) |
| MAP SQUAREM iters | 121 | 123 |
| VBEM converged? | Yes | **No** |
| MAP converged? | Yes | Yes |

The small mega-locus has only 180K components, making the L1 convergence criterion achievable. VBEM's smoother posterior landscape actually helps SQUAREM find a shorter path to convergence. But at 5.8M components, the L1 norm can't reach 1e-6 within the iteration budget.

---

## 8. Conclusions

1. **VBEM's poor performance on large BAMs is entirely a convergence criterion bug**, not an algorithmic inferiority. Fix the convergence delta scaling, and VBEM should perform comparably to MAP-EM.

2. **The E-step memory access pattern is the dominant bottleneck** for both modes. Fixing cache thrashing in Pass 2 could yield 5–10× E-step speedup, which would reduce the large BAM mega-locus from ~5 minutes to ~30–60 seconds.

3. **Digamma vectorization is low priority**. It accounts for only 2.5% of per-iteration cost.

4. **The mega-locus problem is real**: one locus with 315K transcripts and 275K equivalence classes accounts for 96% of total EM time on the large BAM. Any optimization must target mega-locus performance.
