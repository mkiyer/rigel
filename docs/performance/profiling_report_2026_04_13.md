# Rigel Performance Profiling Report

**Date**: 2026-04-13
**Platform**: Linux 4.18.0 (RHEL 8), x86-64, 8 threads
**Python**: 3.12.13 (conda-forge)
**BAM**: VCaP simulation, STAR-aligned, 19.4M fragments (1.8 GB)
**Index**: Human genome, 457,513 transcripts (254K real + 203K synthetic nRNA), 63,472 genes

---

## Executive Summary

A critical O(N×M) performance regression was found in `calibration.py:_build_gdna_fl_model()`, introduced in commit `26001ab`. A per-region Python for-loop iterated 542K regions, each scanning 13.3M fragment-length observations — total ~7.2 trillion element comparisons. This transformed a sub-second operation into a **~45-minute wall-clock hang**.

**Fix**: Replaced the loop with three vectorized NumPy broadcasts. Calibration dropped from ~45 min → **0.58 seconds** (4,600× speedup). Full pipeline now completes in **120 seconds** with no test regressions (1,013/1,013 passed).

---

## Stage Timing Breakdown (Post-Fix)

| Stage | Time | % | Notes |
|-------|------|---|-------|
| `scan_and_buffer` | 80.7s | 67.1% | C++ BAM scan, resolve, model training, buffering |
| `locus_em` | 23.4s | 19.5% | Per-locus SQUAREM EM (29,458 loci, 8 threads) |
| `partition` | 13.2s | 11.0% | Connected components + scatter |
| `fragment_router_scan` | 11.1s | 9.2% | Scoring + CSR construction |
| `fragment_scorer` | 2.3s | 1.9% | FragmentScorer initialization |
| `build_loci` | 1.6s | 1.3% | Locus construction |
| `calibration` | 0.6s | 0.5% | gDNA deconvolution (post-fix) |
| `eb_gdna_priors` | 0.5s | 0.4% | Empirical Bayes priors |
| Other | <0.1s | <0.1% | |
| **Total** | **120.2s** | | **347K fragments/sec** |

## Memory Profile

| Checkpoint | RSS (MB) | Delta |
|-----------|----------|-------|
| After index load | 2,835 | baseline |
| After BAM scan | 5,344 | +2,509 MB (fragment buffer) |
| After scoring/routing | 8,549 | +3,205 MB (CSR arrays) |
| After buffer release | 5,499 | −3,050 MB freed |
| After EM | 5,222 | −277 MB |
| **Peak** | **8,549** | |

**Key finding**: Peak RSS occurs when both the fragment buffer (2.5 GB) and CSR arrays (3.2 GB) coexist in memory before the buffer is released. This is the single largest memory optimization opportunity.

## Locus EM Analysis

| Metric | Value |
|--------|-------|
| Total loci | 29,458 |
| Mega-loci (>50K transcripts) | 1 |
| Total EM time | 72.8s (wall), 23.4s (parallel) |
| Top 10 loci time | 17.7s (24% of EM) |

### Mega-Locus #0 (98,608 transcripts, 4.06M units, 140K ECs)

| Sub-stage | Time | % |
|-----------|------|---|
| SQUAREM iterations (59 iters) | 9.69s | 69.6% |
| Build equivalence classes | 1.43s | 10.3% |
| Assign posteriors | 1.06s | 7.6% |
| Extract fragment data | 0.95s | 6.8% |
| Bias correction | 0.69s | 4.9% |
| Warm start | 0.11s | 0.8% |

### Deterministic Fast Path

Only **3.5%** of scoring units take the deterministic-unambiguous fast path (680K / 19.4M). This is far below the typical 60-80% for mRNA-only indexes, and is a direct consequence of the 203K synthetic nRNA transcripts doubling the candidate set per fragment.

---

## Root Cause: Calibration Regression

### The Bug

In `_build_gdna_fl_model()` ([calibration.py](../../src/rigel/calibration.py#L365)), introduced in commit `26001ab`:

```python
# BEFORE (O(N×M) — ~45 minutes)
for rid_val in range(len(gene_strand)):          # N = 542,510 regions
    if not eligible[rid_val] or not has_strand[rid_val]:
        continue
    gs = gene_strand[rid_val]
    rid_match = rids == rid_val                   # M = 13,266,376 comparisons
    if gs == 1:
        strand_mask |= rid_match & (fstrands == STRAND_POS)
    elif gs == -1:
        strand_mask |= rid_match & (fstrands == STRAND_NEG)
```

### The Fix

```python
# AFTER (O(M) — 0.58 seconds)
obs_eligible = eligible[rids] & has_strand[rids]  # fancy index: O(M)
obs_gs = gene_strand[rids]                         # fancy index: O(M)
strand_mask |= obs_eligible & (obs_gs == 1) & (fstrands == STRAND_POS)
strand_mask |= obs_eligible & (obs_gs == -1) & (fstrands == STRAND_NEG)
```

Replaced N per-region iterations (each scanning all M observations) with three O(M) vectorized broadcasts via NumPy fancy indexing. Semantically identical; validated by 106 calibration tests + 1,013 full suite tests.

---

## Proposed Optimizations

### Priority 1: scan_and_buffer (80.7s, 67.1%)

The BAM scan dominates runtime. It is already parallelized in C++ with OpenMP and htslib multi-threaded I/O. Further optimization opportunities:

1. **StreamingScorer vector pre-allocation**: All 17 `std::vector` members in `StreamingScorer` ([scoring.cpp](../../src/rigel/native/scoring.cpp)) grow via `push_back()` with no `reserve()` calls. For 19.4M fragments at ~2.5 candidates each, this causes ~48.5M push_backs with 15-20 reallocations per vector. **Expected impact: 5-10% of scan time (~4-8s).**

2. **Fragment buffer memory**: 2,509 MB for 19.4M fragments (129 bytes/fragment). The columnar CSR layout is efficient, but `to_scoring_arrays()` upcasts `t_offsets` and `frag_id` from int32→int64 on every chunk ([buffer.py](../../src/rigel/buffer.py#L253)). Accepting int32 in the C++ kernel would save ~155 MB of temporary copies.

### Priority 2: locus_em (23.4s, 19.5%)

1. **Mega-locus dominance**: A single mega-locus (98K transcripts) consumes 13.9s (19% of all EM time). This is intrinsic to the human genome's complex loci and nRNA architecture.

2. **`MAX_K_STACK` overflow**: The EM solver uses a stack-allocated buffer of 512 elements. Loci with >512 components (common with nRNA) fall back to heap allocation. Increasing to 1024 would cover 99%+ of cases.

3. **EC hash cost**: `VecHash` computes O(k) hashes per equivalence class. For mega-loci with k=278 max EC width, this adds up. Caching intermediate hashes could help.

### Priority 3: Memory — Peak RSS Reduction (8.5 GB)

1. **Streaming CSR construction**: Currently, the full fragment buffer (2.5 GB) and full CSR (3.2 GB) coexist at peak. A streaming approach that builds CSR chunks while releasing buffer chunks could cap peak at ~5.5 GB.

2. **nRNA index bloat**: 203K synthetic nRNA transcripts nearly double the index size. This inflates candidate sets, CSR arrays, and EM component counts. Coalescing overlapping nRNA spans or increasing merge tolerance (currently 20bp) could reduce nRNA count by 30-50%.

### Priority 4: partition (13.2s, 11.0%)

The partition stage includes C++ union-find (fast) plus Python scatter operations. The 13.2s is spent building per-locus arrays from global CSR data. Optimizing the numpy scatter/indexing could yield 2-5s savings.

### Not Recommended

- **cProfile**: Adds 15-30× overhead; unsuitable for production profiling. Use `--stages` mode instead.
- **Scoring `-ffast-math` changes**: Already enabled; further SIMD optimization shows diminishing returns given scan I/O dominance.
- **EM convergence tuning**: SQUAREM already provides 3-5× speedup over vanilla EM; further convergence changes risk accuracy.

---

## Summary of Actions Taken

| Action | Status | Impact |
|--------|--------|--------|
| Vectorized `_build_gdna_fl_model()` loop | **Done** | 45 min → 0.58s (4,600× speedup) |
| Full profiling with stage decomposition | **Done** | Baseline established |
| Full test suite regression check | **Done** | 1,013/1,013 passed |

---

## Appendix: Pipeline Statistics

| Metric | Value |
|--------|-------|
| Total reads in BAM | 41,757,162 |
| QC fail / unmapped / supplementary | 1,166,628 / 4,540 |
| Secondary alignments | 1,753,375 |
| Name-sorted fragments | 19,416,686 |
| Proper pairs | 19,416,273 |
| With exonic overlap | 15,803,639 |
| With annotated splice junction | 5,096,688 |
| Same-strand fragments | 14,487,295 |
| Ambig-strand fragments | 1,316,344 |
| Intergenic / gDNA | 3,611,669 |
| Multimapper groups | 215,611 |
| Strand specificity (spliced) | 0.9999 |
| Fragment length mean (global) | 297.5 |
| Fragment length mean (spliced) | 251.2 |
| Calibrated gDNA density (λ_G) | 3.34e-3 |
| Expected gDNA fragments | 8,653,168 / 15,074,279 unspliced |
