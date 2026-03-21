# Phase 3 Performance Plan: Early Pipeline Optimization

**Date:** March 21, 2026  
**Status:** Proposal  
**Prerequisite:** Phases 1–2 complete (calibration + medium-impact optimizations)

---

## Current State After Phase 1 + 2

### Wall-Time Summary

| Stage | Baseline | Phase 2 | % of Total | Speedup |
|-------|----------|---------|------------|---------|
| index_load | 8.2s | 8.2s | — (not timed) | — |
| scan_and_buffer | 108.1s | 109.8s | 43.2% | — |
| calibration | 149.6s | 11.7s | 4.6% | **12.8×** |
| fragment_scorer | 1.5s | 1.4s | 0.5% | — |
| fragment_router_scan | 18.6s | 18.1s | 7.1% | — |
| build_loci | 2.4s | 2.4s | 0.9% | — |
| eb_gdna_priors | 1.5s | 1.5s | 0.6% | — |
| locus_em | 111.0s | 109.3s | 43.0% | — |
| **Total** | **392.8s** | **254.3s** | **100%** | **1.54×** |

### Memory Summary (Simulated BAM, 14.7M fragments)

| Checkpoint | RSS (MB) |
|------------|----------|
| After index load | 3,981 |
| After scan | 17,282 |
| After locus_em (peak) | 20,788 |

### Throughput

| Dataset | Baseline | Phase 2 | Improvement |
|---------|----------|---------|-------------|
| **Simulated (14.7M frags)** | 98,500 frags/s | 152,100 frags/s | 1.54× |
| **Real (309K frags)** | 149,000 frags/s | 290,200 frags/s | 1.95× |

### cProfile Top Functions (Phase 2, self-time)

After removing cProfile artifacts (`_thread.lock.acquire`, `threading.wait`):

| Function | Self Time | Calls | Notes |
|----------|-----------|-------|-------|
| `_marginal_loglik` | 6.5s | 251 | Kappa estimation (Brent's method, vectorized) |
| `numpy.asarray` | 1.4s | 19K | Index/scoring setup overhead |
| `pyarrow.compute.take` | 1.0s | 19K | Pandas ArrowStringArray lookups in `build_loci` |
| `build_gdna_fl_model` | 1.0s | 27 | Still shows 1s (see note below) |
| `pandas.ArrowArray.__getitem__` | 0.7s | 511K | Locus gdna_prior computation |
| `ndarray.astype` | 0.5s | 119 | Type casting in calibration + buffer |
| `_compute_strand_llr_betabinom` | 0.4s | 13 | Calibration strand LLR |
| `numpy.ufunc.at` | 0.4s | 14 | scatter-add in `_compute_fl_llr` |
| `gc.collect` | 0.4s | 1 | Full GC cycle after pipeline stages |
| `_compute_fl_llr` | 0.3s | 13 | Vectorized FL LLR (was 12.6s baseline) |
| `_merged_intervals` | 0.15s | 80K | Per-locus genomic span calculation |

### What's Left to Optimize: Budget Breakdown

**Python-visible time: ~15s** (calibration 11.7s + scorer 1.4s + build_loci 2.4s)

**C++-opaque time: ~239s** (scan 110s + router 18s + locus_em 109s)

Python-level profiling has reached diminishing returns. The remaining 94% of
wall time is in C++ code: BAM scanning (43%), locus EM (43%), and fragment
routing (7%). Optimizing further requires C++ profiling and architectural
changes.

---

## Optimization Targets

### Target 1: Index Loading (8.2s — startup cost)

**Current architecture:**
- `TranscriptIndex.load()` reads 4–5 feather files via `pyarrow.feather.read_table()`
- Constructs Pandas DataFrames (`t_df`, `region_df`, `nrna_df`, `sj_df`)
- Builds Python dicts: `_t_exon_intervals` (254K entries), splice junction maps
- Constructs C++ `FragmentResolver` (copies interval data to C++)
- Builds `cgranges` indexes for overlap queries

**Breakdown (estimated from read_table dominance):**
- Feather read + DataFrame: ~3s (4–5 files, 254K–1.2M rows)
- C++ FragmentResolver construction: ~2s (copies index data to C++)
- cgranges build: ~1s (interval insertion + indexing)
- Dict construction (exon intervals, SJ maps): ~2s

**Optimization opportunities:**

| # | Change | Expected Savings | Effort | Notes |
|---|--------|-----------------|--------|-------|
| I1 | **Lazy-load region_df and nrna_df** | 0.5–1s | Low | Only needed during calibration and locus EM, respectively; defer read until first access |
| I2 | **Cache C++ FragmentResolver to disk** | 2–3s | Medium | Serialize the built C++ resolver state to a binary cache file during `rigel index`; reload via mmap at `rigel quant` time instead of rebuilding from DataFrames |
| I3 | **Memory-map feather files** | 0.5–1s | Low | Use `memory_map=True` in `pf.read_feather()` to avoid copying data into process heap; columns become zero-copy views |
| I4 | **Consolidate index into single file** | 1–2s | Medium | Replace 4–5 separate feather files with one Arrow IPC file or a custom binary format; one file open + one mmap vs five |
| I5 | **Pre-build cgranges binary** | 0.5–1s | Medium | Serialize cgranges index during `rigel index`; reload at quant time |

**Estimated Phase 3 index savings: 3–5s → index load drops to 3–5s.**

### Target 2: BAM Scanning (109.8s — 43% of wall time)

**Current architecture:**
- C++ multi-threaded: 1 reader thread (htslib, 2 BGZF decomp threads) + N worker threads (resolve + accumulate)
- Reader: sequential `sam_read1()` → group by qname → push to bounded queue
- Workers: pop qname groups → resolve fragments via cgranges → append to per-worker `FragmentAccumulator`
- After scan: merge all worker accumulators → finalize (compute ambig_strand) → return raw bytes to Python
- Python: `_FinalizedChunk.from_raw(raw)` → `buffer.inject_chunk(chunk)` → optional spill to disk

**Data volume:**
- 38.7M BAM records → 14.7M fragments → ~19M candidate entries (CSR)
- Per-fragment storage: ~37 bytes fixed + ~50 bytes CSR = ~87 bytes
- Total accumulator: ~1.3 GB data + std::vector overhead

**Bottleneck analysis:**

The 110s scan breaks down roughly as:
1. **BAM I/O + BGZF decompression** (~25–35s): htslib `sam_read1()` with 2 decompression threads on a 1.7 GB BAM
2. **Record parsing + qname grouping** (~15–20s): `parse_bam_record()` extracts tags, cigar, coordinates
3. **Fragment resolution** (~30–40s): cgranges interval queries, exon block assembly, chimera detection
4. **Accumulator append + merge** (~15–20s): Vector appends with occasional reallocation; final merge of N worker accumulators

(Exact breakdown requires Instruments/perf profiling — cProfile only sees the C++ call as one opaque block.)

**Optimization opportunities:**

| # | Change | Expected Savings | Effort | Notes |
|---|--------|-----------------|--------|-------|
| S1 | **Increase BGZF decompression threads** | 5–15s | Trivial | Current: `n_decomp_threads=2` (hardcoded default). On Apple M-series with 10 cores, 4 decomp threads should improve I/O throughput. Expose as config parameter. |
| S2 | **Pre-allocate accumulator vectors with `reserve()`** | 2–5s | Low | Currently vectors grow with amortized doubling. For 14.7M fragments, ~24 doublings per vector × 17 vectors = significant realloc overhead. After first chunk/scan, we know approximate fragment count → `reserve()` based on BAM record count estimate. |
| S3 | **Eliminate worker accumulator merge** | 5–10s | Medium | Instead of N workers → merge → one accumulator, have each worker finalize independently → yield N independent chunks to Python. Avoids the O(total_data) copy during merge. |
| S4 | **Avoid double-copy in finalize + from_raw** | 3–5s | Medium | Currently: C++ vectors → `nb::bytes` (copy 1) → `np.frombuffer().copy()` (copy 2). Use nanobind's buffer protocol or capsule to transfer ownership of C++ vectors directly to NumPy arrays without copying. |
| S5 | **Streaming chunk output** | 10–20s + memory | Hard | Instead of accumulating ALL fragments in C++ then returning to Python, have workers produce _FinalizedChunks periodically (e.g., every 1M fragments). Python can buffer/spill while scan continues. Eliminates the accumulator merge, reduces peak memory, and enables pipelining with calibration. |
| S6 | **Pre-fetch BAM index for random access** | 0s (arch change) | Hard | Not directly applicable (name-sorted BAM), but if we supported coordinate-sorted BAMs, per-chromosome parallelism would be possible. |

**Estimated Phase 3 scan savings: 15–35s → scan drops to 75–95s.**

### Target 3: Fragment Buffering & Spill (included in scan time)

**Current architecture:**
- Python `FragmentBuffer` holds chunks in memory until `max_memory_bytes` (default 5 GB) is exceeded
- Spill: Arrow IPC (Feather v2) with LZ4 compression to temp dir
- For 14.7M fragments: 1 in-memory chunk (~1 GB) + 0 spills (below threshold)
- BUT: the C++ accumulator holds ALL data before finalization (~1.3 GB + overhead)
- So effective peak = C++ accumulator + Python chunk = ~2.5 GB for fragment data alone

**The real memory problem is the C++ accumulator lifecycle:**

```
┌──────────────────────────────────┐
│ scan() running                   │  C++ accumulator growing: 0 → 1.3 GB
│                                  │  Per-worker accumulators: N × (1.3GB / N)
│ ← all 38.7M reads processed →   │
│                                  │
│ merge workers → main acc         │  Peak: main_acc(1.3GB) + worker_accs(1.3GB) = 2.6 GB
│ worker_accs freed                │  Main acc: 1.3 GB
│ scan() returns                   │
├──────────────────────────────────┤
│ finalize_accumulator()           │  nb::bytes created: + 1.3 GB (copy)
│                                  │  Peak: 2.6 GB (acc + bytes)
├──────────────────────────────────┤
│ _FinalizedChunk.from_raw()       │  NumPy arrays: + 1.3 GB (copy)
│                                  │  Peak: 3.9 GB (acc + bytes + numpy)
├──────────────────────────────────┤
│ buffer.inject_chunk(chunk)       │  Buffer holds chunk reference
│ GC frees raw bytes + acc         │  Drops to ~1.3 GB
└──────────────────────────────────┘
```

The ~13.3 GB delta observed is the accumulator data (1.3 GB) plus heap fragmentation,
`std::vector` over-allocation (up to 2× capacity), per-worker duplicates, and the
fact that macOS doesn't return freed pages to the OS (RSS never shrinks). Much of
the 13.3 GB is "virtual" overhead, not useful data.

**Optimization opportunities:**

| # | Change | Expected Memory Savings | Effort | Notes |
|---|--------|------------------------|--------|-------|
| B1 | **Streaming chunks from C++** (= S5 above) | 5–10 GB peak | Hard | Workers produce chunks of 1M fragments; Python spills to disk as they arrive. C++ never holds more than N_workers × 1M fragments. Peak accumulator: ~N×87 MB instead of 1.3 GB. |
| B2 | **Zero-copy finalization** (= S4 above) | Eliminates 1.3 GB copy | Medium | Transfer C++ vector ownership to Python via capsules instead of bytes+copy. |
| B3 | **Lower max_memory_bytes default** | Reduces Python buffer peak | Trivial | Currently 5 GB; lowering to 1–2 GB triggers earlier spills but fragment data is mostly from C++ accumulator, not Python buffer. Only helps when *multiple* chunks are in-memory. |
| B4 | **Use int16 for exon_bp/intron_bp/unambig_intron_bp** | ~25% of CSR data | Medium | These columns rarely exceed 32K. int16 halves their storage: saves ~30 MB per 1M fragments. |
| B5 | **Compress CSR in-memory** | 20–40% of CSR | Hard | Delta-encode t_indices within each fragment (sorted), run-length encode frag_lengths. Adds decode overhead but reduces memory. |

**Estimated memory reduction: 5–10 GB peak (primarily from streaming + zero-copy).**

### Target 4: Fragment Router / Scoring (18.1s — 7% of wall time)

**Current architecture:**
- C++ `fused_score_buffer()` performs two passes over all chunks:
  - **Pass 1 (count):** Count total EM units + candidates (no output)
  - **Pass 2 (fill):** Pre-allocate output arrays, score each fragment, fill CSR
- Per-fragment scoring: strand log-prob, fragment-length log-prob, overhang penalty, NM mismatch penalty, coverage weight

**Optimization opportunities:**

| # | Change | Expected Savings | Effort | Notes |
|---|--------|-----------------|--------|-------|
| R1 | **Single-pass scoring** | 5–9s | Medium | The two-pass design exists to pre-allocate exact output sizes. Alternative: use growing vectors (like accumulator), then transfer to NumPy. Saves one full pass over 14.7M fragments. Trade: ~10% vector over-allocation. |
| R2 | **Pipeline scoring with scan** | 5–10s | Hard | If streaming chunks (S5/B1), scoring can start while scan is still running. Score-as-you-go on each chunk, then merge CSR outputs. |
| R3 | **Skip re-scoring deterministic fragments** | 1–3s | Medium | Deterministic unambiguous fragments (79K units) are scored but don't enter EM. Their scores can be computed more cheaply or skipped entirely if not needed for output. |

**Estimated Phase 3 router savings: 5–9s → drops to 9–13s.**

---

## Proposed Phase 3 Implementation Plan

### Phase 3a: Quick Wins (low effort, immediate impact)

| # | Item | Target | Savings | Effort |
|---|------|--------|---------|--------|
| 3a.1 | Expose `n_decomp_threads` in BamScanConfig | scan | 5–15s wall time | Trivial: add config field + pass to `scanner.scan()` |
| 3a.2 | Pre-allocate accumulator vectors with `reserve()` | scan | 2–5s wall time | Low: estimate from BAM record count |
| 3a.3 | Lazy-load `region_df` and `nrna_df` | index | 0.5–1s startup | Low: `@cached_property` patterns |
| 3a.4 | Memory-map feather files | index | 0.5–1s startup | Low: pass `memory_map=True` |

**Expected Phase 3a impact: scan drops ~10–20s, index load drops ~1–2s.**

### Phase 3b: Zero-Copy & Merge Elimination (medium effort)

| # | Item | Target | Savings | Effort |
|---|------|--------|---------|--------|
| 3b.1 | Eliminate double-copy in finalize + from_raw | scan + buffer | 3–5s + 1.3 GB memory | Medium: nanobind capsule-based ownership transfer |
| 3b.2 | Eliminate worker accumulator merge | scan | 5–10s wall time | Medium: each worker finalizes independently |
| 3b.3 | Single-pass fused_score_buffer | router | 5–9s wall time | Medium: switch to growing vectors |

**Expected Phase 3b impact: scan drops another 8–15s, router drops 5–9s. Memory -1.3 GB.**

### Phase 3c: Streaming Architecture (high effort, high reward)

| # | Item | Target | Savings | Effort |
|---|------|--------|---------|--------|
| 3c.1 | Streaming chunks from C++ scanner to Python | scan + buffer + memory | 10–20s + 5–10 GB memory | Hard: requires callback or coroutine interface between C++ and Python |
| 3c.2 | Score-while-scanning pipeline | router + scan | 5–10s (overlap) | Hard: requires chunk-level CSR merge |
| 3c.3 | Cache C++ FragmentResolver to binary | index | 2–3s startup | Medium: serialize during `rigel index` |
| 3c.4 | Consolidate index to single Arrow IPC file | index | 1–2s startup | Medium: changes `rigel index` output format |

**Expected Phase 3c impact: total pipeline time drops to ~180–200s (30–40% from current Phase 2). Peak memory drops to ~12–15 GB.**

---

## Projected Timeline

| Phase | Target Time | Peak Memory | Key Changes |
|-------|-------------|-------------|-------------|
| **Baseline** | 392.8s | 21.3 GB | — |
| **After Phase 1+2** | 254.3s | 20.8 GB | Calibration scipy + Brent + vectorize |
| **After Phase 3a** | ~235–244s | 20.8 GB | BGZF threads, reserve(), lazy index |
| **After Phase 3b** | ~215–225s | ~19.5 GB | Zero-copy, merge elimination, 1-pass scorer |
| **After Phase 3c** | ~180–200s | ~12–15 GB | Streaming scan + pipeline scoring |

### Beyond Phase 3

The remaining ~180–200s is dominated by:
- **Locus EM: 109s** — Algorithmic challenge (mega-locus with 226K transcripts). Requires locus decomposition (graph partitioning) or transcript pruning.
- **BAM I/O: ~25–35s** — htslib sequential read on name-sorted BAM. Fundamentally limited by BGZF decompression throughput.
- **Overlap resolution: ~30–40s** — cgranges queries per fragment are O(log n + k). Potential for batch query optimization.

These would be Phase 4+ targets requiring deeper algorithmic changes.

---

## Appendix: Phase 2 Results Detail

### Phase 2 Changes Applied

| # | Change | Before (Phase 1) | After (Phase 2) | Improvement |
|---|--------|-------------------|------------------|-------------|
| 5 | Brent optimizer replacing golden-section | `_golden_section_max`: ~35 evals/call, 0.574s cumtime (13 calls) | `_bounded_max` via `minimize_scalar`: ~19 evals/call, 7.6s cumtime (13 calls) for kappa estimation (improvement in eval count, total kappa time dominated by per-eval cost of `_marginal_loglik` at 26ms/eval) | 34% fewer evaluations (**251 vs ~455**); calibration: 17.8s → 11.7s |
| 6 | Vectorize exon tuple construction | `tuple(int(x) for x in ...)`: 1.9M generator calls, ~1.5s | `tuple(.tolist())` + `np.cumsum()`: 254K cumsum calls (0.26s) + 763K tolist calls (0.13s) | fragment_scorer 1.5s → 1.4s (7% faster) |
| 7 | Cache `fragment_classes` in `_FinalizedChunk` | `fragment_classes` property recomputed on every access | Computed once, cached via `object.__setattr__` in slots dataclass | Negligible in practice (called only once per chunk) |

### Simulated BAM (large, 14.7M fragments)

| Stage | Baseline | Phase 1 | Phase 2 | Cumulative Speedup |
|-------|----------|---------|---------|-------------------|
| scan_and_buffer | 108.1s | 109.4s | 109.8s | — |
| calibration | 149.6s | 17.8s | 11.7s | **12.8×** |
| fragment_scorer | 1.5s | 1.5s | 1.4s | 1.1× |
| fragment_router_scan | 18.6s | 18.3s | 18.1s | — |
| locus_em | 111.0s | 110.6s | 109.3s | — |
| **Total** | **392.8s** | **261.5s** | **254.3s** | **1.54×** |
| Peak RSS | 21,339 MB | 20,956 MB | 20,788 MB | -551 MB |

### Real BAM (small, 309K fragments)

| Stage | Baseline | Phase 1 | Phase 2 | Cumulative Speedup |
|-------|----------|---------|---------|-------------------|
| scan_and_buffer | 2.5s | 2.5s | 2.5s | — |
| calibration | 8.9s | 1.4s | 1.2s | **7.4×** |
| fragment_scorer | 1.6s | 1.5s | 1.4s | 1.1× |
| locus_em | 1.7s | 1.4s | 1.5s | — |
| **Total** | **16.4s** | **8.7s** | **8.4s** | **1.95×** |

### Function Call Count Reduction

| Metric | Baseline | Phase 2 | Reduction |
|--------|----------|---------|-----------|
| Total function calls | ~160M | 9.25M | **17×** |
| `builtins.min` calls | 145M | 0 | Eliminated |
| `_vec_lgamma` calls | 2,808 | 0 | Replaced by scipy |
| `_marginal_loglik` evals | ~455 | **251** | 1.8× fewer (Brent) |
| Generator calls (exon tuples) | 3.8M | 0 | Replaced by tolist/cumsum |

---

## Appendix: Platform

- macOS-26.3-arm64 (Apple Silicon, M-series)
- Python 3.12.13
- Clang 19.1.7 (C++ compiler via conda-forge)
- htslib (via conda-forge)
- condaforge + bioconda dependencies
