# Performance Optimization Plan

## Profiling Summary

**Workload**: 1.2 GB simulated oracle BAM (24M reads → 12M paired-end fragments),
254,461 transcripts, 63,472 genes, HeLa cervix gDNA-high (ss=0.90).
**Platform**: macOS ARM64, Python 3.13, 4 threads.

### Stage Timing Breakdown (85.2s total)

| Stage                   | Time (s) |     % | Notes                                                   |
|-------------------------|----------|-------|---------------------------------------------------------|
| `scan_and_buffer`       |    48.19 | 56.6% | C++ BAM scan (4 workers + 2 htslib decomp threads)     |
| `locus_em`              |    23.68 | 27.8% | C++ batch EM (4 threads, 29 266 loci, 10.3M ambig frags)|
| `fragment_router_scan`  |     4.46 |  5.2% | C++ fused two-pass scoring + CSR build                  |
| `build_loci`            |     3.84 |  4.5% | C++ union-find + Python post-processing                 |
| `eb_gdna_priors`        |     3.03 |  3.6% | Python/numpy hierarchical EB                            |
| `fragment_scorer`       |     1.90 |  2.2% | Python `from_models()` — per-transcript exon tuple build|
| `create_estimator`      |     0.08 |  0.1% | TSS groups + estimator init                             |
| Other                   |    <0.05 |  0.0% | geometry, nRNA init, finalize_models                    |

**Index load**: 7.4s (separate from pipeline wall time).

### Memory Profile

| Phase               | RSS (MB) | Delta (MB) |
|---------------------|----------|------------|
| After index load    |    3,827 |          — |
| After scan (t=46s)  |   10,197 |     +6,370 |
| Peak (steady)       |   14,980 |     +4,783 |

- Buffer: 12M fragments, 1 chunk spilled to disk (2.1 GB Arrow/LZ4).
- RSS never reclaimed (Python memory fragmentation — `malloc` doesn't return pages to OS).
- 15 GB peak is ~12× the BAM file size.

### cProfile Hotspots (Python-visible, 39.5s CPU)

| Function                          | Self (s) | Cum (s) | Calls    | Notes                          |
|-----------------------------------|----------|---------|----------|--------------------------------|
| `_thread.lock.acquire`            |    33.04 |   37.05 | 1,280    | thread sync (scan + EM)        |
| `threading.wait`                  |     0.80 |   43.44 | 320      | memory sampling thread          |
| `ArrowArray.__getitem__`          |     0.96 |    1.44 | 742,064  | pandas Arrow column access      |
| `list.append`                     |     0.79 |    0.79 | 12.4M    | CSR array.array growth          |
| `numpy.asarray`                   |     0.89 |    0.89 | 453      | array conversion                |
| `gc.collect`                      |     0.35 |    0.35 | 1        | manual GC                       |
| `scoring.py genexpr` ×2           |     0.49 |    0.49 | 3.85M    | from_models() tuple build       |
| `_snap_rss_current`               |     0.16 |    1.05 | 320      | fixed (was 47s pre-fix)        |
| `_compute_intergenic_density`     |     0.07 |    0.12 | 1        | interval merging                |

---

## Optimization Targets

### P1 — BAM Scan I/O: Eliminate `bam_dup1` Copies (HIGH)

**Impact**: scan_and_buffer = 48.2s (56.6% of total)

**Problem**: The producer thread deep-copies every BAM record via `bam_dup1()`
before pushing to the bounded work queue. For 24M records at ~200–500 bytes
each, this creates 4.8–12 GB of transient allocations + copies.

**Solution**: Pool-based `bam1_t` recycling.

1. Pre-allocate a pool of `bam1_t` objects (size = queue capacity × avg records/group).
2. Producer grabs from pool, reads into pre-allocated record, pushes group.
3. Worker processes group, returns records to pool when done.
4. Eliminates malloc/free per record.

**Alternative**: Batch mode — producer reads all records for a qname group into
a thread-local buffer, then moves the entire buffer to the queue. Workers parse
from the raw buffer (zero-copy within the group).

**Estimated gain**: 5–15s (10–30% of scan time) — depends on malloc cost share.

### P2 — Locus EM Convergence Tuning (HIGH)

**Impact**: locus_em = 23.7s (27.8%)

**Problem**: The C++ batch EM processes 29,266 loci with SQUAREM acceleration.
With 4 threads and dynamic chunk dispatch (16 loci/chunk), the current EM is
already well-parallelized. However:

- Default convergence delta and iteration cap may allow over-iteration on easy loci.
- Equivalence class construction uses `std::unordered_map` per locus.
- `digamma()` is called O(n_transcripts × iterations) per locus.

**Optimizations**:

1. **Early termination for trivial loci**: Loci with ≤2 transcripts or ≤10
   ambiguous units can bypass full EM and use direct assignment.
2. **Tighter convergence guard**: Reduce max iterations for small loci where
   convergence is typically fast (<5 iterations).
3. **Tabulated digamma**: Pre-compute digamma for common α ranges
   (0.01–1000) in a LUT with interpolation, avoiding the asymptotic series.
4. **Pool equivalence class storage**: Pre-allocate `EmEquivClass` vectors
   per-thread instead of per-locus hash map construction.

**Estimated gain**: 3–8s (15–35% of EM time).

### P3 — `build_loci` Python Post-Processing (MEDIUM)

**Impact**: build_loci = 3.8s (4.5%)

**Problem**: After C++ connected_components (fast), the Python post-processing
builds per-component transcript/unit lists using dicts, `.setdefault()`, and
Python loops over 10.3M units. The wall-clock (3.8s) far exceeds cProfile's
Python-visible time (0.09s), suggesting memory pressure (page faults during
large array copies on the 14+ GB working set).

**Optimizations**:

1. **Return locus assignments from C++**: Extend `connected_components_native()`
   to return pre-sorted transcript lists and unit assignments per component,
   eliminating the Python dict-building and sorting loops.
2. **Avoid `t_indices.copy()`**: The 50M-element array copy triggers page
   faults at high RSS — operate on the original array in-place or in C++.

**Estimated gain**: 2–3s.

### P4 — Hierarchical EB gDNA Priors (MEDIUM)

**Impact**: eb_gdna_priors = 3.0s (3.6%)

**Problem**: `compute_eb_gdna_priors()` calls `_compute_chrom_gdna_rates()`
(0.042s cProfile) and `_compute_per_locus_gdna_rates()` (0.039s cProfile), but
wall-clock is 3.0s — again memory-pressure dominated.

**Optimizations**:

1. **Vectorize per-locus accumulation**: Replace Python `defaultdict(float)`
   loops with `np.add.at()` grouped scatter.
2. **Fuse with locus EM setup**: The gDNA prior computation could be folded
   into the batch EM preparation, avoiding a separate pass.
3. **Reduce intermediate allocations**: Reuse arrays from earlier stages.

**Estimated gain**: 1–2s.

### P5 — `fragment_scorer.from_models()` Tuple Build (LOW)

**Impact**: fragment_scorer = 1.9s (2.2%)

**Problem**: `FragmentScorer.from_models()` iterates all 254,461 transcripts
in a Python loop, calling `index.get_exon_intervals(t_idx)` and building
`(starts, ends, cumsum)` tuples. The two genexprs account for 0.5s of Python
overhead (3.85M calls).

**Optimizations**:

1. **Pre-compute exon data in C++**: Build the `_t_exon_data` dict during
   index loading or in the `NativeFragmentScorer` constructor.
2. **Vectorized construction**: Use the existing CSR exon arrays
   (`exon_offsets`, `exon_starts`, `exon_ends`) directly in numpy slicing
   instead of per-transcript Python iteration.

**Estimated gain**: 1–1.5s.

### P6 — `fragment_router_scan` CSR Construction (LOW)

**Impact**: fragment_router_scan = 4.5s (5.2%)

**Problem**: The fused C++ scoring pass is already well-optimized. The 4.5s
includes loading 1 spilled chunk from disk (0.09s) and the fused scoring
of 11.2M fragments. The `_scan_native()` self-time is only 0.003s in cProfile
— almost all time is in C++ (invisible to cProfile).

**Optimizations**:

1. **Avoid chunk spilling altogether**: Increase `max_memory_bytes` or
   implement streaming scoring during the initial BAM scan (eliminate the
   two-pass architecture for the scoring stage).
2. **Reduce `list.append` overhead**: The CSR builder uses `array.array`
   for dynamic growth (12.4M appends = 0.79s). Pre-size arrays with the
   known fragment count from the scan pass.

**Estimated gain**: 0.5–1.5s.

### P7 — Index Loading (LOW)

**Impact**: index_load = 7.4s

**Problem**: Index loading reads multiple files (intervals, splice junctions,
transcript metadata) and builds in-memory data structures. This is a one-time
startup cost.

**Optimizations**:

1. **Memory-mapped index**: Store the index in a binary format that can be
   directly mmap'd, avoiding deserialization.
2. **Lazy loading**: Load only the structures needed for the requested
   operation (e.g., skip splice junction cgranges if not needed).

**Estimated gain**: 3–5s startup (not pipeline throughput).

---

## Memory Optimization Targets

### M1 — Reduce Peak RSS from 15 GB (HIGH)

**Target**: Reduce peak from 15 GB to <8 GB for the 1.2 GB BAM workload.

**Problem breakdown**:
- Index: 3.8 GB (254K transcripts, cgranges, splice junction maps)
- Scan buffer: ~2.1 GB in-memory (11.2M fragments, columnar)
  → spilled to disk, but RSS not reclaimed
- CSR arrays (array.array accumulators): ~1.3 GB during scoring
- ScoredFragments CSR: ~2–3 GB (10.3M units × candidates)
- EM workspace: ~0.5–1 GB (per-thread sub-problems)
- Python overhead + fragmentation: ~2–3 GB

**Strategies**:

1. **Streaming architecture** (eliminate buffer-then-score): Score fragments
   during the BAM scan pass instead of buffering all fragments and rescoring.
   This eliminates the 2.1 GB buffer + 1.3 GB array.array CSR accumulators.

2. **Explicit memory return to OS**: Call `malloc_trim(0)` (Linux) or
   equivalent after freeing large buffers. On macOS, consider using
   `madvise(MADV_FREE)` on freed pages.

3. **Compact index representation**: Use int32 for coordinates instead of
   int64 where genomic positions allow. Use packed structs for transcript
   metadata.

4. **Reduce ScoredFragments memory**: The CSR stores per-candidate
   log-likelihoods, coverage weights, tx_starts, tx_ends. Some of these
   could be recomputed during EM instead of stored.

5. **Memory-mapped spill chunks**: Use mmap for spilled Arrow chunks
   instead of reading them entirely into memory.

### M2 — Buffer Spill Overhead (LOW)

**Current**: 1 chunk spilled (2.1 GB Arrow/LZ4), loaded back in 0.09s.

This is already efficient. Arrow IPC + LZ4 is fast. No optimization needed
unless workloads produce many spilled chunks.

---

## Priority Ranking

| ID | Target                              | Impact   | Effort | Est. Gain  |
|----|-------------------------------------|----------|--------|------------|
| P1 | BAM scan `bam_dup1` pool            | HIGH     | Medium | 5–15s      |
| P2 | EM convergence + trivial-loci skip  | HIGH     | Medium | 3–8s       |
| M1 | Reduce peak RSS to <8 GB            | HIGH     | Large  | 7+ GB      |
| P3 | `build_loci` C++ return             | MEDIUM   | Small  | 2–3s       |
| P4 | EB gDNA priors vectorization        | MEDIUM   | Small  | 1–2s       |
| P5 | `from_models()` vectorization       | LOW      | Small  | 1–1.5s     |
| P6 | CSR pre-sizing                      | LOW      | Small  | 0.5–1.5s   |
| P7 | Index loading                       | LOW      | Large  | 3–5s start |
| M2 | Buffer spill                        | LOW      | —      | Already OK |

---

## Profiler Bug Fix (Done)

The `_snap_rss_current()` function in `scripts/profiler.py` was re-importing
`ctypes`, calling `ctypes.util.find_library("c")` (6 `posix.stat` calls each),
and redefining a ctypes Structure class **on every invocation**. With 1-second
memory sampling, this added 47s of pure overhead to the original cProfile run.

**Fix applied**: Hoisted ctypes initialization to module level. The mach
task_info handle, libc reference, and structure class are now created once at
import time and reused for all subsequent `_snap_rss_current()` calls.

Result: `_snap_rss_current` dropped from 47.3s to 1.05s cumulative (320 calls).

---

## Recommended Implementation Order

1. **P2** — EM trivial-loci early exit (quickest HIGH-impact win, C++ only)
2. **P3** — C++ locus builder return (small effort, clear gain)
3. **P1** — bam_dup1 pool (largest absolute gain, moderate effort)
4. **P5** — from_models vectorization (quick Python fix)
5. **P4** — EB priors vectorization (Python/numpy)
6. **M1** — Streaming architecture (largest structural change, highest memory gain)
