# Phase 5 — Scaling Bottleneck Analysis & Improvement Plan

**Date**: 2026-03-24
**Status**: Planned
**Prerequisite**: Phase 4 (quantification-stage cleanup) complete

---

## Motivation

Phase 4 optimized the Python quantification stages, bringing the real BAM
(1.6M fragments) to 17.0s. However, profiling on a **large simulated dataset**
(15M fragments, 50% gDNA contamination, GENCODE v46 human transcriptome)
reveals fundamentally different bottlenecks that only emerge at scale.

### Dataset: sim_ccle_hela_salmon / gdna_high_ss_0.90_nrna_default

| Metric | Value |
|--------|-------|
| Total BAM records | 42,036,167 |
| Buffered fragments | 14,732,713 |
| Unique alignments | 13,145,960 |
| Multimapping | 1,586,753 |
| Intergenic (gDNA) | 1,734,280 |
| Loci | 12,402 |
| EM units | 12,894,819 |
| Transcripts | 457,513 |
| Peak RSS | 15,810 MB |

### Profile: 205.6s wall time

| Stage | Time | % | Notes |
|-------|------|---|-------|
| **scan_and_buffer** | **93.1s** | **45.3%** | C++ htslib + resolution |
| **locus_em** | **90.0s** | **43.8%** | C++ EM solver |
| calibration | 14.1s | 6.9% | Scipy optimizer loop |
| fragment_router_scan | 6.4s | 3.1% | C++ scoring |
| fragment_scorer | 1.0s | 0.5% | Phase 4b optimized |
| build_loci | 0.8s | 0.4% | Phase 4a optimized |
| eb_gdna_priors | 0.3s | 0.2% | Phase 4a optimized |

### The Mega-Locus Problem

The dominant EM bottleneck is a **single mega-locus**:

| Metric | Value |
|--------|-------|
| Transcripts | 376,144 (82% of total) |
| EM units | 12,393,179 (96.1% of all units) |
| Components | 376,145 (n_t + 1 gDNA) |
| gdna_span | ~1.7 Gbp |
| Est. work share | **100.0%** of total EM work |

This mega-locus is a connected component formed when multimapping gDNA
fragments bridge otherwise-independent gene loci. With 50% gDNA
contamination, intergenic reads resolve to multiple genomic regions,
linking them. This is a worst-case scenario but realistic for heavily
contaminated samples.

**EM solver behavior on the mega-locus**:
- SQUAREM runs 3 E-steps per iteration × ~333 iterations = ~1000 E-steps
- Each E-step iterates over 12.4M units × variable k candidates
- Equivalence class building: hash 12.4M unit keys → merging/sorting
- The `parallel_estep` function parallelizes across 8 threads, but the
  per-iteration overhead of thread synchronization and Kahan reduction
  across 376K components limits scaling
- Deterministic sorting of EC rows (for reproducibility) is O(N log N)
  on the 12.4M-unit dataset

### cProfile Evidence (Python-visible)

| Function | Self Time | Calls | Notes |
|----------|-----------|-------|-------|
| `_thread.lock.acquire` | 188.6s | 7,608 | GIL wait → all time in C++ |
| `calibration._marginal_loglik` | 7.6s | 283 | Scipy optimizer inner loop |
| `calibration.build_gdna_fl_model` | 1.3s | 33 | Fragment-length deconv |
| `gc.collect` | 0.7s | 1 | Still present (post-EM) |
| `numpy.asarray` | 1.3s | 92 | Array conversion overhead |

**Key insight**: 188.6s in lock acquire = almost all time is in C++ code
(scan and EM). Python overhead is negligible at this scale. The remaining
optimizations must focus on C++ algorithmic improvements.

---

## Comparison: Real BAM vs Simulated BAM

| Metric | Real BAM (1.6M) | Sim BAM (14.7M) | Ratio |
|--------|-----------------|------------------|-------|
| Fragments | 1,130,994 | 14,732,713 | 13.0× |
| scan_and_buffer | 9.0s | 93.1s | 10.3× |
| locus_em | 3.6s | 90.0s | 24.7× |
| calibration | 2.4s | 14.1s | 5.9× |
| Max locus tx | 162,574 | 376,144 | 2.3× |
| Max locus units | 510,690 | 12,393,179 | 24.3× |
| Peak RSS | 3,910 MB | 15,810 MB | 4.0× |

**The EM scales super-linearly (24.7× for 13× more data)** because:
1. The mega-locus grows quadratically: more gDNA → more bridging → bigger component
2. SQUAREM convergence slows with more components (376K vs 162K)
3. Equiv class building is O(N log N) on the mega-locus

---

## Design Principles

1. **Measure before optimizing.** py-spy native profiling is unavailable on
   macOS ARM64. Use C++ `std::chrono` instrumentation in the EM solver to
   decompose the 90s into sub-stages (EC build, sorting, SQUAREM iterations,
   assignment).
2. **Tackle the biggest bottleneck first.** The mega-locus EM (90s) and
   BAM scan (93s) dominate — focus there.
3. **Algorithmic improvements over micro-optimization.** The EC sort is
   O(N log N) on 12.4M rows — no amount of SIMD will fix that. Need to
   reduce the problem size or eliminate redundant work.
4. **Preserve correctness and reproducibility.** Golden output tests must
   still pass. Deterministic output is required.

---

## Phase 5a — EM Solver Instrumentation

**Target**: Identify which sub-steps of `process_locus` dominate for the
mega-locus.

### Approach

Add optional `std::chrono` timing to the `process_locus` lambda and
accumulate across all loci. Report timing breakdown back to Python.

#### Sub-stages to measure:
1. `extract_locus_sub_problem` — data extraction and global→local mapping
2. `apply_bias_correction_uniform` — bias correction
3. `build_equiv_classes` — hash map construction + EC building
4. EC sorting (deterministic ordering) — sort ECs + sort rows within ECs
5. `compute_ovr_prior_and_warm_start` — prior computation
6. `run_squarem` — SQUAREM iterations (E-steps + convergence)
7. `assign_posteriors` — posterior assignment

### Implementation

Add a `TimingAccumulator` struct with atomic accumulators for each
sub-stage, passed through to `process_locus`. Return the timings in the
Python result dict.

### Files Modified

- `src/rigel/native/em_solver.cpp` — Add timing instrumentation
- `src/rigel/estimator.py` — Report timings if present

### Verification

Run the sim BAM profile and report per-stage breakdown within the 90s.

---

## Phase 5b — Reduce Equivalence Class Sorting Cost

**Target**: EC deterministic sort overhead in mega-loci.

### Problem

The `build_equiv_classes` function sorts 12.4M rows by log-likelihood
fingerprint for determinism. This is O(N × k × log N) where k is the
number of candidates per EC. For the mega-locus with N=12.4M, this
is extremely expensive.

### Option 1: Tiered determinism

For mega-loci (>100K units), use a faster deterministic ordering that
doesn't require full comparison-based sorting:

- Hash each row's log-likelihood vector to a 64-bit fingerprint
- Sort by fingerprint (single comparison instead of k comparisons)
- Collisions are rare and produce ULP-level differences (acceptable)

This reduces per-comparison cost from O(k) to O(1).

### Option 2: Batch row insertion order

Instead of sorting rows after the fact, insert rows into equivalence
classes in deterministic order by pre-sorting the unit indices before
the hash map pass. Since unit indices come from CSR offsets which are
already sorted, this might be achievable with zero sorting.

### Files Modified

- `src/rigel/native/em_solver.cpp` — Optimize EC sort strategy

---

## Phase 5c — Reduce Mega-Locus EM Convergence Cost

**Target**: SQUAREM iteration count on the mega-locus.

### Problem

With 376K components, SQUAREM typically converges slowly because:
1. Many near-zero components create numerical noise
2. The gDNA component competes with 376K transcript components
3. The convergence check (normalized L1 delta) accumulates error
   across all 376K components

### Option 1: Component pruning during EM

After N initial iterations, prune components with theta < epsilon.
This reduces the effective n_components and speeds up subsequent
iterations. The pruned components get zero abundance.

### Option 2: Adaptive convergence threshold

For mega-loci, relax the convergence delta proportionally to n_components.
The per-component accuracy is already limited by the number of fragments
per component.

### Option 3: Early termination for gDNA-dominated loci

If the gDNA component posterior exceeds 0.99 after initial iterations,
skip remaining EM iterations and assign all fragments to gDNA.

### Files Modified

- `src/rigel/native/em_solver.cpp` — EM convergence improvements

---

## Phase 5d — BAM Scan I/O Optimization

**Target**: `scan_and_buffer` (93.1s)

### Problem

The BAM scan reads 42M records through htslib. The three-thread
architecture (reader, workers, main) has potential bottlenecks:

1. **htslib I/O** — sequential BAM decompression + parsing dominates
   the reader thread. With `n_decomp_threads=4`, BGZF decompression
   is parallel, but `sam_read1()` itself is sequential.

2. **Fragment resolution** — each qname group goes through cgranges
   interval overlap queries. With heavy gDNA contamination, intergenic
   reads resolve to multiple overlapping regions.

3. **Main thread callback** — GIL acquisition for each chunk callback
   introduces latency.

### Option 1: Increase chunk size

The default `chunk_size=1M` means ~15 chunk callbacks for 14.7M fragments.
Each callback requires GIL acquisition. Increasing to 2M or 4M reduces
callback overhead.

### Option 2: Profile htslib vs resolution

Without py-spy native profiling on macOS, instrument the C++ scanner
with `std::chrono` timers to decompose:
- BAM record reading time
- qname group resolution time
- Model training time
- Chunk finalization time

### Option 3: Reduce intergenic resolution overhead

For intergenic reads (gDNA), the fragment resolver does cgranges overlap
queries that return multiple hits. If gDNA reads are identifiable early
(no transcript overlap), they could take a fast path that skips full
resolution.

### Files Modified

- `src/rigel/native/bam_scanner.cpp` — Add timing instrumentation
- `src/rigel/config.py` — Potentially adjust default chunk_size

---

## Phase 5e — Calibration Optimization

**Target**: `calibration` (14.1s, 6.9%)

### Problem

The calibration EM runs 16 regions × ~18 scipy optimizer iterations =
283 `_marginal_loglik` calls at 27ms each. This is 8× slower on the
sim data because:
1. More regions with evidence (more fragments)
2. Each `_marginal_loglik` call evaluates the beta-binomial PMF over
   larger count vectors
3. The M-step `build_gdna_fl_model` rebuilds the gDNA FL model 33 times

### Option 1: Analytical kappa gradient

The current `estimate_kappa_marginal` uses scipy.optimize.minimize_scalar
with bracket search. Providing an analytical gradient would reduce the
number of function evaluations.

### Option 2: Reduce calibration iterations

The calibration EM uses `max_iterations=50` but often converges in ~15.
Add convergence checking to exit early.

### Option 3: Vectorize beta-binomial computation

The `_marginal_loglik` function computes `betaln` element-wise over
region count arrays. Batching with numpy could reduce overhead.

### Files Modified

- `src/rigel/calibration.py` — Convergence checking, vectorization

---

## Phase 5f — Memory Optimization

**Target**: Peak RSS 15.8 GB

### Problem

The 15.8 GB peak RSS is dominated by:
1. `FragmentBuffer` in-memory chunks: 14.7M fragments × ~100 bytes = ~1.5 GB
2. `ScoredFragments` CSR: 12.4M units × ~5 candidate arrays = ~2 GB
3. EM solver working set: 376K components × multiple temp vectors = ~500 MB
4. htslib BAM buffer + decompression: ~500 MB

RSS snapshots from the profile:
- after_scan: 9,394 MB (buffer + scan working set)
- after_router_scan: 13,422 MB (+ scored fragments CSR)
- after_locus_em: 15,810 MB (+ EM working set)

### Option 1: Release buffer after scoring

Currently `buffer.release()` is called after `_score_fragments`. Verify
this is happening correctly and the RSS actually drops.

### Option 2: Streaming EM (process loci during scoring)

Instead of building the full CSR (ScoredFragments) and then iterating
loci, interleave scoring and EM. This would avoid holding the full CSR
in memory simultaneously.

### Option 3: Disk-backed buffer spilling

The buffer already supports spilling to disk (`spill_dir` config) but
defaults to in-memory. For large datasets, enable spilling automatically
when RSS exceeds a threshold.

---

## Summary & Prioritization

| Phase | Target | Est. Impact | Difficulty | Priority |
|-------|--------|-------------|------------|----------|
| **5a** | EM instrumentation | 0s (diagnostic) | Low | **1st** |
| **5b** | EC sort optimization | 5-20s | Medium | **2nd** |
| **5c** | EM convergence | 10-30s | Medium | **3rd** |
| **5d** | BAM scan optimization | 5-15s | Medium | **4th** |
| **5e** | Calibration optimization | 5-10s | Low | **5th** |
| **5f** | Memory optimization | RSS only | High | **6th** |

### Execution Plan

```
5a → 5b → 5c → 5d → 5e → 5f
```

Phase 5a (instrumentation) must come first — it gives us precise data
to guide the remaining phases. Phases 5b and 5c are the highest-impact
EM optimizations. Phase 5d addresses the other major bottleneck (scan).
Phase 5e is lower priority but straightforward. Phase 5f is a
longer-term architectural improvement.

### What's NOT in Phase 5

- **Multithreaded Python** — GIL prevents parallelism; already not a bottleneck
- **GPU acceleration** — Overkill for this problem size; adds significant complexity
- **Alternative EM algorithms** — Switching away from SQUAREM would require
  extensive validation. SQUAREM is well-understood and robust.
- **BAM format alternatives** — CRAM support would help I/O but doesn't
  address the EM bottleneck.

---

## Appendix: Locus Size Distribution (Measured)

From the profiled run (12,402 loci, 12,894,819 total EM units):

| Size Bucket (tx) | Loci | Total Units | % Units |
|-------------------|------|-------------|---------|
| 1 | 5,553 | 19,395 | 0.15% |
| 2-5 | 4,512 | 83,430 | 0.65% |
| 6-10 | 1,057 | 84,024 | 0.65% |
| 11-50 | 1,199 | 271,876 | 2.11% |
| 51-100 | 74 | 40,706 | 0.32% |
| 101-1,000 | 6 | 2,209 | 0.02% |
| 100,001+ | **1** | **12,393,179** | **96.1%** |

### The Mega-Locus Domination

The single mega-locus (locus 0) is extraordinary:

| Metric | Value |
|--------|-------|
| Transcripts | 376,144 (82% of 457K total) |
| EM units | 12,393,179 (96.1% of all units) |
| gdna_span | 1,705,490,533 bp (~1.7 Gbp) |
| Estimated work (n_t × n_u) | **4.66 trillion** (100.0% of total) |

**The mega-locus contains 100.0% of the EM computational work** (by the
units × components heuristic). The remaining 12,401 loci combined contribute
<0.01% of EM work.

**Top-10 loci by transcript count**:

| Locus | Tx | Units | gdna_span |
|-------|------|-------|-----------|
| 0 | 376,144 | 12,393,179 | 1,705,490,533 |
| 9792 | 122 | 580 | 134,499 |
| 1732 | 116 | 1,018 | 16,123 |
| 9334 | 113 | 240 | 23,736 |
| 5754 | 109 | 57 | 35,716 |
| 3491 | 108 | 280 | 10,242 |
| 2061 | 105 | 34 | 3,976 |
| 8833 | 99 | 722 | 10,725 |
| 4065 | 97 | 108 | 63,827 |
| 489 | 93 | 127 | 28,961 |

The gap between locus 0 (376K tx) and locus 9792 (122 tx) is **3,082:1**.
This is a fundamentally bimodal distribution — one monster component
formed by multimapping gDNA fragments bridging transcripts across the
genome, plus 12,401 small well-behaved loci.

### Implication for Phase 5

Optimizing the EM solver for small/medium loci is irrelevant — they
collectively take <0.01% of EM time. **All EM optimization effort
must target the mega-locus.** Strategies that reduce the effective
problem size (component pruning, sparse E-step, subgraph decomposition)
will have orders-of-magnitude more impact than algorithmic speedups
to the core E-step loop.
