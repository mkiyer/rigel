# Performance Analysis v2 — Post-Optimization

**Date:** 2026-03-28
**System:** 2× Intel Xeon Gold 6154 (36 cores, 3.0 GHz), 186 GB RAM, DDR4-2666
**Software:** Python 3.12.13, GCC 8.5.0, Rigel with fused column-sum + stack-local row buffer optimizations

## Benchmark Configuration

| Property | Small BAM | Large BAM |
|----------|-----------|-----------|
| BAM file | mctp_LBX0338 (203 MB) | mctp_MI_1190 (1.2 GB) |
| Fragments | 4,145,225 | 24,149,525 |
| Buffered | 728,785 | 9,250,691 |
| Unique | 562,699 | 8,691,392 |
| Multimapping | 166,086 | 559,299 |
| Intergenic (gDNA) | 21,006 | 1,731,595 |
| Index | Human GENCODE, 457,513 transcripts, 63,472 genes |
| EM mode | VBEM | VBEM |
| Threads | 8 | 8 |

## Current Timing Summary

### Stage Breakdown

| Stage | Small BAM | % | Large BAM | % |
|-------|----------|---|-----------|---|
| scan_and_buffer | 25.1s | 58.3% | 174.1s | 22.4% |
| calibration | 6.6s | 15.2% | 37.3s | 4.8% |
| fragment_scorer | 2.0s | 4.8% | 2.0s | 0.3% |
| fragment_router_scan | 0.9s | 2.0% | 8.9s | 1.1% |
| build_loci | 0.7s | 1.6% | 1.3s | 0.2% |
| eb_gdna_priors | 0.6s | 1.4% | 0.7s | 0.1% |
| **locus_em** | **7.1s** | **16.6%** | **553.4s** | **71.2%** |
| **TOTAL** | **43.0s** | | **777.7s** | |

### Throughput

| BAM | Throughput |
|-----|-----------|
| Small | 96,383 frags/s |
| Large | 31,053 frags/s |

### Improvements vs. Pre-Optimization Baseline

| Stage | Original (Large) | Current (Large) | Speedup |
|-------|-----------------|-----------------|---------|
| locus_em | 657.8s | 553.4s | 1.19× |

*Note: The 1.53× speedup measured in isolation was from the optimized E-step kernel. The wall-clock difference is smaller because locus_em includes extract, build_ec, and assign phases that were not optimized.*

## Memory Profile

### RSS Snapshots (MB)

| Checkpoint | Small BAM | Large BAM | Delta (Large) |
|------------|-----------|-----------|---------------|
| before (index loaded) | 2,860 | 2,861 | — |
| after_scan | 2,398 | 3,885 | +1,024 |
| after_calibration | 2,404 | 3,885 | +0 |
| after_router_scan | 2,404 | 5,497 | +1,612 |
| after_buffer_release | 2,076 | 3,775 | −1,722 |
| after_locus_em | 2,098 | 4,510 | +735 |
| after_cleanup | 2,098 | 2,899 | −1,612 |
| **Peak RSS** | **2,860** | **5,783** | |

### Memory Breakdown (Large BAM)

| Component | Size | Notes |
|-----------|------|-------|
| Base (index) | ~2,861 MB | Loaded once, persistent |
| Fragment buffer | ~1,024 MB | 9.25M frags × ~110 bytes |
| Router CSR arrays | ~1,612 MB | ScoredFragments output |
| EM working set | ~735 MB | Mega-locus dominated |
| Peak additional | ~2,922 MB | Buffer + router overlap |

**Peak RSS** occurs between router_scan and buffer_release when both the fragment buffer and the CSR arrays coexist. Releasing the buffer before the EM runs recovers ~1.7 GB.

### Mega-Locus EM Memory (Large BAM)

The single mega-locus (315K components, 5.8M units, 275K equiv classes) dominates:

| Array | Size |
|-------|------|
| ll_flat (log-likelihoods) | 393.8 MB |
| comp_idx (component indices) | 196.9 MB |
| wt_flat (weights) | 46.2 MB |
| col_acc (Kahan accumulators) | 7.6 MB |
| theta + em_totals | 5.0 MB |
| **Total** | **~650 MB** |

*Previously, `ec.scratch` added ~394 MB on top of this. That allocation was eliminated in the stack-local row buffer optimization.*

## Locus-Level EM Analysis

### Time Distribution

| Property | Small BAM | Large BAM |
|----------|-----------|-----------|
| Total loci | 14,530 | 13,882 |
| Mega-loci | 1 | 1 |
| Top-1 locus time | 5.5s (32%) | 547.3s (92.3%) |
| Top-5 loci time | 6.7s (39%) | 552.0s (93.1%) |
| Median SQUAREM iters | 2 | 2 |
| Mean SQUAREM iters | 3.0 | 3.2 |
| P95 SQUAREM iters | 8 | 9 |
| P99 SQUAREM iters | 18 | 22 |
| Max SQUAREM iters | 333 | 333 |
| Loci hitting max iters | ~3 | 4 |

### Mega-Locus Detail (Large BAM)

| Property | Value |
|----------|-------|
| Transcripts | 315,030 |
| Units (fragments) | 5,777,023 |
| Components (2n_t + 1) | 315,031 |
| Equiv classes | 275,385 |
| EC total elements | 49,219,285 |
| Max EC width | 315 |
| Max EC depth | 10,077 |
| SQUAREM iterations | 333 (max) |

### Sub-Phase Breakdown (Large BAM, All Loci)

| Sub-phase | Time | % of EM |
|-----------|------|---------|
| extract | 16.3s | 2.7% |
| bias | 2.0s | 0.3% |
| build_ec | 5.8s | 1.0% |
| warm_start | 0.3s | 0.0% |
| **squarem** | **564.5s** | **95.2%** |
| assign | 4.1s | 0.7% |

### Per-Iteration Cost (Mega-Locus)

| Metric | Value |
|--------|-------|
| Total SQUAREM time | 535.9s |
| SQUAREM iterations | 333 |
| Per-iteration | 1,609 ms |
| Per E-step (~3 per iter) | ~536 ms |
| Data read per E-step | 0.64 GB |
| Effective bandwidth | 1.2 GB/s |
| Theoretical DDR4-2666 (6-channel) | ~120 GB/s |
| Bandwidth utilization | ~1% |

Low bandwidth utilization indicates the bottleneck is **not** memory bandwidth but rather **compute-bound work** within the E-step (log-sum-exp, SIMD exp, Kahan accumulation, col-sum reduction) and/or **poor cache locality** across the large equivalence class matrices.

## Calibration Analysis

### cProfile Hotspots (Large BAM, 37.3s total)

| Function | Self-time | Calls | Cost/call |
|----------|-----------|-------|-----------|
| `_marginal_loglik` (betaln) | 13.0s | 315 | 41 ms |
| `estimate_kappa_marginal` | 0.3s | 18 | 16 ms |
| `build_gdna_fl_model` (bincount) | 1.3s | 37 | 36 ms |
| `_e_step` | 0.3s | 18 | 15 ms |
| `_compute_fl_llr` (add.at scatter) | 0.4s | 18 | 21 ms |
| `_m_step` | 0.4s | 17 | 25 ms |
| `_compute_strand_llr_betabinom` | 0.6s | 18 | 36 ms |
| `_compute_density_llr_gaussian` | 0.2s | 18 | 10 ms |
| `np.add.at` (scatter-add) | 0.8s | 18 | 47 ms |
| `astype` (int→float conversions) | 1.9s | 178 | 10 ms |
| `np.asarray` | 1.6s | 98 | 17 ms |

**Key insight:** `_marginal_loglik` dominates calibration (35% of calibration time). It calls `scipy.special.betaln` vectorized over all valid regions, 315 times (18 regions × ~17.5 Brent iterations). The betaln function involves two `lgamma` calls per element.

### Calibration EM Structure

- 18 EM iterations (large BAM) / 41 iterations (small BAM)
- Per iteration: E-step → M-step → FL model rebuild → convergence check
- `estimate_kappa_marginal` called once per EM iteration — its inner Brent optimization calls `_marginal_loglik` ~17 times per call

## Scan & Buffer Analysis

### scan_and_buffer Timing

| BAM | Time | Throughput |
|-----|------|-----------|
| Small | 25.1s | 165K reads/s |
| Large | 174.1s | 139K reads/s |

*Note: scan_and_buffer timing is highly cache-dependent. Cold-cache runs can be 2× slower (174s vs ~92s for the large BAM).*

### Architecture

- 1 reader thread (htslib BAM I/O with BGZF decompression)
- N worker threads (fragment resolution, model training)
- GIL cycling for Python callbacks (chunk finalization)
- `n_decomp_threads = 4` for BGZF decompression (default)

The reader thread is the likely bottleneck — single-threaded htslib I/O with BGZF decompression fans out to N workers.

### fragment_router_scan (8.9s on Large BAM)

Fused C++ scoring path: per-fragment candidate enumeration, strand/FL/splice scoring, CSR assembly. This is compute-bound in C++ with SIMD optimizations already applied.

---

## Optimization Priorities (v2)

### Priority 1: VBEM Convergence on Mega-Loci [CRITICAL]

**Impact:** Eliminates 535.9s → expects ~10-50s (10-50× improvement on mega-locus)
**Risk:** Results will change (FP-precision level)

The single mega-locus accounts for 92.3% of all EM time on the large BAM. It runs 333 SQUAREM iterations (the hard cap) without converging. The convergence criterion is L1 norm across 5.8M component weights, which doesn't scale well to mega-loci.

**Potential approaches:**
1. **Relative convergence criterion**: Use max relative change `max(|θ_new - θ_old| / max(θ_old, ε))` instead of L1 sum. This scales independently of component count.
2. **Convergence by log-likelihood**: Monitor ELBO or bound change between iterations instead of parameter L1 norm.
3. **Hierarchical decomposition**: Break connected component into sub-loci based on shared fragment density (effectively "cutting" weak bridges in the multimapper graph).
4. **Pruning near-zero components**: After initial convergence, freeze components with θ < ε to reduce the effective k.

**Status:** Deferred pending ground truth validation (see `docs/performance/VBEM_CONVERGENCE_ISSUE.md`).

### Priority 2: Calibration Speedup [HIGH]

**Impact:** 37.3s → ~5-10s (4-7× improvement)
**Risk:** FP-precision-level changes acceptable; gross changes unlikely

The calibration EM (`calibrate_gdna`) consumes 37.3s on the large BAM (4.8% of total). The dominant cost is `_marginal_loglik` calling scipy's `betaln` 315 times.

**Potential approaches:**

#### 2a. Vectorize kappa estimation across regions
Currently `estimate_kappa_marginal` is called once per EM iteration with a Brent scalar optimizer. Each call evaluates `_marginal_loglik` ~17 times, each computing `betaln` over all valid regions. Pre-computing the (k, n) pairs and evaluating all kappa candidates in a single vectorized batch could reduce Python overhead.

#### 2b. Analytical kappa approximation
The beta-binomial κ estimation can be approximated using method-of-moments instead of marginal likelihood maximization:
$$\hat{\kappa} = \frac{\bar{p}(1-\bar{p})}{s^2} - 1$$
This eliminates the entire Brent + betaln cost at the expense of statistical efficiency.

#### 2c. Cache betaln intermediate values
When Brent evaluates nearby kappa values, the (k, n) data is the same — only the (a, b) parameters change. Pre-computing lgamma lookup tables for common (k, n) combinations could reduce redundant computation.

#### 2d. Pre-compute FL histogram matrix (attempted, reverted)
Previously attempted: pre-compute per-region FL histogram matrix (n_regions × 1001), replace per-iteration `np.bincount` with matrix-vector multiply. Reverted because FP summation order changes affected results. Could re-attempt with user's permission for FP-precision-level changes.

#### 2e. Reduce type conversion overhead
`astype` conversions consume 1.9s and `np.asarray` 1.6s. These are scattered across calibration hot paths and could potentially be eliminated by ensuring arrays are already the correct dtype when created.

### Priority 3: Scan Throughput [MEDIUM]

**Impact:** 174s → ~90-120s (1.5-2× improvement)
**Risk:** No result changes

The BAM scan is I/O-bound (single reader thread + BGZF decompression). Variable performance depends on OS page cache state.

**Potential approaches:**

#### 3a. Increase BGZF decompression threads
Default is 4. On a 36-core system, increasing to 8-12 may improve decompression throughput. This is a simple config change.

#### 3b. Pre-load BAM into page cache
Use `vmtouch` or `cat > /dev/null` before profiling to ensure consistent measurements. Not a code change but important for benchmarking.

#### 3c. Concurrent read + process
The current architecture (reader → worker queue → callback) is already pipelined. The bottleneck may be the Python GIL cycling during chunk callbacks. Profiling the C++ reader thread in isolation could reveal if pure I/O throughput is the limit.

### Priority 4: Router Memory Optimization [MEDIUM]

**Impact:** Peak RSS 5,783 → ~4,200 MB (−1.6 GB)
**Risk:** No result changes

The peak RSS occurs when both the fragment buffer (~1 GB) and router CSR arrays (~1.6 GB) coexist in memory. The buffer is released only after the router completes.

**Potential approaches:**

#### 4a. Stream-and-release buffer chunks during router scan
Instead of loading all chunks into `chunk_arrays` before calling `fused_score_buffer`, process chunks one at a time and release each before loading the next. This would reduce peak overlap to one chunk (~65 MB) instead of the full buffer.

#### 4b. In-place router output
Write CSR arrays directly from the buffer scan rather than building a temporary copy. This requires changes to the C++ `fused_score_buffer` API.

### Priority 5: EM Data Size Reduction [MEDIUM-LOW]

**Impact:** −200-400 MB working set, potential cache-locality improvement
**Risk:** FP-precision-level changes

The mega-locus EC data is 650 MB. The ll_flat array (394 MB) uses float64 for log-likelihoods.

**Potential approaches:**

#### 5a. Single-precision log-likelihoods
Use float32 for `ll_flat` storage (half the memory: 394 → 197 MB). The E-step kernel would upcast to float64 for accumulation. This may affect numerical stability in the Kahan summation but the impact would be at the FP-precision level.

#### 5b. Sparse ll_flat representation
Many EC rows have only a few non-zero entries (median EC width is small). A CSR-within-CSR scheme for ll_flat could significantly reduce memory for sparse ECs while maintaining dense access for wide ECs.

#### 5c. Component index compression
`comp_idx` uses int32 (197 MB). For loci with k < 65536, int16 would halve this. For the mega-locus (k=315031), int32 is required but a mixed scheme could help smaller loci.

### Priority 6: build_exon_csr Python Loop [LOW]

**Impact:** 0.3-0.4s → ~0.01s (30× improvement)
**Risk:** No result changes

`index.build_exon_csr()` iterates over 228K multi-exon transcripts in Python, calling `np.cumsum` per-transcript. This is called once per run and costs ~0.4s.

**Potential approaches:**
- Vectorize the entire CSR construction using a single `np.cumsum` over the flattened lengths array with offset-based slicing. This eliminates 228K Python-level function calls.
- Or move to C++ alongside the exon interval data.

### Priority 7: Post-EM Assignment [LOW]

**Impact:** 4.1s → ~1s (4× improvement)
**Risk:** No result changes (deterministic mode) or seed-dependent changes (sample mode)

`assign_posteriors` takes 4.1s across all loci. For the mega-locus it's 2.6s. The main cost is log-sum-exp per unit and scatter-add to accumulator arrays.

**Potential approaches:**
- Batch similar-width units for SIMD-friendly processing
- Pre-allocate accumulator arrays at maximum size and zero-fill instead of reallocating

---

## Summary Table

| # | Optimization | Est. Impact | Risk Level | Result Change |
|---|-------------|-------------|------------|---------------|
| 1 | VBEM convergence criterion | 500s+ | High | FP-precision |
| 2 | Calibration speedup | 25-30s | Medium | FP-precision |
| 3 | Scan throughput | 50-80s | Low | None |
| 4 | Router memory reduction | −1.6 GB | Low | None |
| 5 | EM data size reduction | −200 MB + cache | Medium | FP-precision |
| 6 | build_exon_csr vectorize | 0.3s | Very low | None |
| 7 | Post-EM assignment | 3s | Low | None |

## Comparison with v1 Analysis

| Change | v1 Status | v2 Status |
|--------|-----------|-----------|
| Fused column-sum | Planned | ✅ Done (1.15×) |
| Stack-local row buffer | Planned | ✅ Done (+0.33× additional) |
| Scratch elimination | Planned | ✅ Done (−394 MB) |
| AVX-512 fast_exp | Planned | ✅ Done (minimal impact due to early-zero skip) |
| VBEM convergence | Documented | → Priority 1 |
| Calibration FL matrix | Attempted | Reverted (FP changes) → Priority 2 |
| Scan throughput | Not analyzed | → Priority 3 |
| Memory optimization | Not analyzed | → Priority 4 |

## Appendix: Raw Profiling Data

- Small BAM: `results/v2_profile_small/profile_summary.json`, `locus_stats_default.json`
- Large BAM: `results/profile_vbem_large/profile_summary.json`, `locus_stats_vbem_8t.json`
