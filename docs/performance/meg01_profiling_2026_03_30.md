# MEG-01 Profiling: Extract Phase Optimization Results

**Date**: 2026-03-30  
**Dataset**: CCLE MEG-01 (12 GB BAM, 107M reads, 49.4M fragments, 457K transcripts)  
**System**: 2× Intel Xeon Gold 6154 (36 cores), 186 GB RAM, RHEL 8.10, Lustre filesystem

## Changes Tested

Extract phase simplification (see `extract_phase_simplification.md`):
- Deleted ~870 lines of dead code (non-partitioned EM path)
- Rewrote `extract_locus_sub_problem_from_partition()`:
  - Removed epoch-based dedup machinery (`BestCandidate`, `seen_epoch[nc]`, `best_buf[nc]`, `dirty_comps`)
  - Pre-allocated output arrays via `resize(max_out)` + write cursor (no `push_back()`)
  - Reusable `std::vector<LocalCandidate>` sort buffer across units (zero allocs after first unit)
- All 977 tests pass, golden outputs bit-for-bit identical

## Results

### Wall-Time Comparison

| Stage | V1 (baseline) | V2 (optimized) | V2 (normalized†) |
|-------|--------------|----------------|-------------------|
| scan_and_buffer | 1196s (51.9%) | 722s (30.6%) | — (I/O variance) |
| locus_em | 950s (41.2%) | 1457s (61.7%) | ~950s |
| fragment_router_scan | 80s (3.5%) | 92s (3.9%) | ~80s |
| calibration | 70s (3.0%) | 82s (3.5%) | ~70s |
| **TOTAL** | **2303s** | **2361s** | **~2120s** |

†Normalized by system load factor of 2.6× (measured via untouched `build_ec` code path).
V2 ran under heavy cluster contention; V1 ran on a quieter node.

### C++ EM Sub-Stage Breakdown

| Sub-stage | V1 | V2 (normalized) | Speedup |
|-----------|-----|-----------------|---------|
| **extract** | **212.1s (10.7%)** | **31.3s (1.6%)** | **6.78×** |
| squarem | 1709.5s (85.9%) | ~1855s (93.4%) | ~1.0× (noise) |
| build_ec | 29.8s (1.5%) | ~same | — |
| assign | 24.9s (1.2%) | ~same | — |
| bias | 13.1s (0.7%) | ~same | — |
| **TOTAL** | **1991.0s** | **~1810s** | **1.10×** |

Estimated wall-time savings: **~180s (8% reduction in C++ EM time)**.

### Correctness Verification

- All 14,575 loci produce **identical equivalence classes** (matched by locus index)
- Zero differences in `n_equiv_classes` or `ec_total_elements` for any locus
- 977 tests pass, golden output bit-for-bit identical

## Current Bottleneck Analysis

### Time Composition (V1 baseline, 2303s wall)

```
scan_and_buffer     1196s  ██████████████████████████  51.9%  (I/O bound, Lustre)
locus_em             950s  █████████████████████       41.2%  (CPU bound)
  └─ squarem         1710s  ══════════════════         85.9% of C++ EM
  └─ extract          212s  ══                         10.7% of C++ EM → NOW 31s
  └─ build_ec          30s                              1.5% of C++ EM
  └─ assign            25s                              1.2% of C++ EM
  └─ bias              13s                              0.7% of C++ EM
frag_router_scan      80s  ███                          3.5%
calibration           70s  ███                          3.0%
everything else        7s                               0.3%
```

### Non-Converged Loci: The Dominant Problem

133 loci (0.9% of all loci) hit the 333-iteration cap and consume **57.6% of all SQUAREM time**:

| Category | Count | SQUAREM Time | % of SQUAREM |
|----------|-------|-------------|--------------|
| Mega (159K tx) | 1 | 764s | 44.7% |
| Large (100-10K tx) | 16 | 51s | 3.0% |
| Medium (10-100 tx) | 109 | 165s | 9.7% |
| Small (≤10 tx) | 7 | 4s | 0.2% |
| **Total non-converged** | **133** | **984s** | **57.6%** |
| Converged | 14,442 | 725s | 42.4% |

### Mega-Locus Deep Dive

The single mega-locus dominates all EM time:

| Metric | Value |
|--------|-------|
| Transcripts | 159,851 |
| Units (fragments) | 18,897,853 |
| Equivalence classes | 407,807 |
| EC total elements | 182,698,326 |
| Max EC width | 261 |
| Max EC depth | 73,527 |
| SQUAREM iterations | 333 (maxed out) |
| Per-iteration time | 2.295s |
| Working set | 2.2 GB |
| Effective bandwidth | 1.0 GB/s (8 threads) |
| Theoretical bandwidth | ~100 GB/s per socket |
| **Bandwidth utilization** | **~1%** |

## Optimization Plan (Priority Order)

### P0: SQUAREM Early Termination / Stall Detection
**Target savings: 200-400s wall time (10-20%)**

133 non-converged loci run all 333 iterations even when making negligible progress.
Options:
- **Stale detection**: track `max|θ_new - θ_old|` history; stop if plateau detected for K consecutive iterations
- **Damped SQUAREM**: reduce acceleration factor near convergence to prevent oscillation that defeats the convergence criterion
- **Adaptive max iterations**: scale max_iterations based on locus complexity

### P1: Mega-Locus Component Pruning Mid-EM
**Target savings: 200-600s wall time (10-25%)**

The mega-locus has 159K components, but most converge to near-zero within ~50 iterations. Currently every E-step computes posteriors for all 159K components even when only ~1K are non-negligible.
- After N warmup iterations, identify components with θ < ε and freeze them at zero
- Remove frozen columns from EC CSR (compress in-place or rebuild)
- Reduces E-step scan volume proportionally to active fraction

### P2: E-Step Cache Locality / EC Reordering
**Target savings: 100-300s wall time (5-15%)**

The mega-locus achieves only 1% of theoretical memory bandwidth. The CSR traversal has poor spatial locality because EC component indices are scattered across the 159K-element θ vector.
- **Sort ECs** to group those sharing component indices → better θ cache reuse
- **Blocked E-step**: partition ECs by component range, process each block while θ slice is cache-hot
- **Prefetch hints** for θ lookups during CSR scan

### P3: Router-Scorer Fusion
**Target savings: 50-80s wall time (2-4%)**

`fragment_router_scan` (80s) does a second pass over scored fragments to build CSR data for EM. Could be fused with the scoring pass to eliminate the separate stage entirely.

### P4: GC / Allocation Pressure
**Target savings: 10-15s wall time (~1%)**

V2 triggered 16 GC cycles (15.4s) vs V1's 1 cycle (0.8s). The partition path creates more temporary Python objects. Investigate buffer reuse in `partition_and_free`.

## Key Takeaway

The extract optimization achieved its goal (6.78× speedup, extract dropped from 10.7% to 1.6% of EM time). The dominant remaining bottleneck is SQUAREM convergence: **133 non-converged loci (0.9%) consume 57.6% of SQUAREM time, with a single 159K-transcript mega-locus responsible for 44.7%**. The next optimization should target convergence behavior (P0/P1), not per-iteration throughput.
