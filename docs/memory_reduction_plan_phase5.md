# Memory Reduction Plan — Phase 5

## 1. Motivation

RNA-seq libraries range enormously in size:

| Tier | Library Size | PE Reads | Frequency |
|---|---|---|---|
| Small | 1–5M fragments | 0.5–2.5M PE | ~10% |
| Typical | 20–50M fragments | 10–25M PE | ~60% |
| Large | 50–100M fragments | 25–50M PE | ~20% |
| Very large | 100–200M fragments | 50–100M PE | ~9% |
| Extreme | 200–400M+ fragments | 100–200M+ PE | ~1% |

The tool must **comfortably** handle up to 200M fragments (100M PE
reads, the 95th percentile).  Typical usage (20–50M fragments) should
run on a machine with 16 GB RAM.  The current implementation, profiled
at 24M fragments, peaked at 20.3 GB RSS — already too much for 16 GB
machines and linearly worse at higher library sizes.  After Phase
5A+5B (implemented), peak is 18.5 GB — still over budget.


## 2. Current Memory Profile (24M Fragments, 254K Transcripts)

### Measured Memory Timeline

```
Time(s)   RSS(MB)   Delta    Phase
───────   ───────   ─────    ──────────────────────────
  0       3,868       —      Index loaded
 91      11,867    +8,000    Buffer chunks finalized
 98      16,846    +4,979    ScoredFragments CSR allocated
105      20,212    +3,366    Estimator + batch EM started
139      20,212        0     End — nothing freed
```

**Peak RSS: 20.2 GB** (16.3 GB delta above 3.9 GB baseline).

### Per-Structure Memory Accounting

#### FragmentBuffer (`_FinalizedChunk`)

Per-fragment fixed fields: **38 bytes**

| Field | Type | Bytes |
|---|---|---|
| `splice_type` | uint8 | 1 |
| `exon_strand` | uint8 | 1 |
| `sj_strand` | uint8 | 1 |
| `num_hits` | uint16 | 2 |
| `merge_criteria` | uint8 | 1 |
| `chimera_type` | uint8 | 1 |
| `ambig_strand` | uint8 | 1 |
| `frag_id` | int64 | 8 |
| `read_length` | uint32 | 4 |
| `genomic_footprint` | int32 | 4 |
| `genomic_start` | int32 | 4 |
| `nm` | uint16 | 2 |
| `t_offsets` | int64 | 8 |

Per-candidate variable fields: **20 bytes**

| Field | Type | Bytes |
|---|---|---|
| `t_indices` | int32 | 4 |
| `frag_lengths` | int32 | 4 |
| `exon_bp` | int32 | 4 |
| `intron_bp` | int32 | 4 |
| `unambig_intron_bp` | int32 | 4 |

**Total per fragment:** ~38 + 20 × avg_candidates ≈ **98 bytes** (at
~3 candidates/fragment average).

With the default 2 GiB in-memory spill budget, the buffer holds ~21M
fragments before spilling to disk.

#### ScoredFragments (CSR)

Per EM-unit: **36 bytes**

| Field | Type | Bytes |
|---|---|---|
| `offsets` | int64 | 8 |
| `locus_t_indices` | int32 | 4 |
| `locus_count_cols` | uint8 | 1 |
| `is_spliced` | bool (1B) | 1 |
| `gdna_log_liks` | float64 | 8 |
| `genomic_footprints` | int32 | 4 |
| `frag_ids` | int64 | 8 |
| `frag_class` | int8 | 1 |
| `splice_type` | uint8 | 1 |

Per candidate: **29 bytes**

| Field | Type | Bytes |
|---|---|---|
| `t_indices` | int32 | 4 |
| `log_liks` | float64 | 8 |
| `count_cols` | uint8 | 1 |
| `coverage_weights` | float64 | 8 |
| `tx_starts` | int32 | 4 |
| `tx_ends` | int32 | 4 |

**ScoredFragments has no spill mechanism — it is entirely in RAM.**

#### AbundanceEstimator

~22 arrays sized by `num_transcripts` (254K).  Total: **~100 MB**.
Negligible relative to buffer and ScoredFragments.

#### Scanner Accumulator (`FragmentRouter`)

During `builder.scan(buffer)`, the `FragmentRouter` accumulates EM
candidate data into Python `array.array` objects (one per field).
These grow dynamically to the same total size as ScoredFragments.

When `_to_np()` converts `array.array → numpy` (via `.copy()`), **both**
the source `array.array` and the destination `numpy` array exist
simultaneously.  This doubles the peak memory during ScoredFragments
construction.


## 3. Memory Scaling Projections

Assumptions: 25% of fragments enter EM, average 5 candidates per EM
unit, average 3 candidates per buffer fragment.

| Frags (M) | Buffer | SF (CSR) | Scanner Peak | Peak w/ Spill |
|---:|---:|---:|---:|---:|
| 1 | 0.1 GB | 0.04 GB | 0.1 GB | **0.2 GB** |
| 5 | 0.5 GB | 0.21 GB | 0.4 GB | **0.9 GB** |
| 10 | 0.9 GB | 0.42 GB | 0.8 GB | **1.8 GB** |
| **24** | **2.2 GB** | **1.0 GB** | **2.0 GB** | **4.0 GB** |
| 50 | 4.6 GB | 2.1 GB | 4.2 GB | **6.2 GB** |
| 100 | 9.1 GB | 4.2 GB | 8.4 GB | **10.4 GB** |
| **200** | **18.3 GB** | **8.4 GB** | **16.9 GB** | **18.9 GB** |
| 400 | 36.5 GB | 16.9 GB | 33.7 GB | **35.7 GB** |

**"Peak w/ Spill"** = min(Buffer, 2 GB) + 2 × ScoredFragments (the
`array.array` + numpy duplication during CSR construction).  Baseline
RSS (~4 GB for index) should be added.

**At 200M fragments, the tool would need ~23 GB** (18.9 + 4 GB
baseline) — far beyond a typical 16 GB machine.  The dominant cost
is ScoredFragments + the scanner accumulator duplication.

### Observation

The profiled 24M-fragment dataset shows 20.2 GB peak RSS vs the model's
prediction of ~8 GB (4 GB baseline + 4 GB working).  The ~12 GB gap
is attributable to:

1. The profiled dataset has a **much higher candidate-to-unit ratio**
   (~10-15 candidates per EM unit rather than the assumed 5), inflating
   ScoredFragments from 1 GB to ~4-5 GB.
2. C++ htslib internal buffers and Python interpreter overhead (~2-3 GB).
3. The scanner accumulator + numpy copy duplication (~4-5 GB peak).
4. Pandas DataFrames in TranscriptIndex (~1-2 GB for 254K transcripts).

The scaling model is **conservative** — real-world memory may exceed
these projections by 1.5–2× depending on the candidate ratio, which
varies with genome annotation complexity (more overlapping transcripts
= more candidates per fragment).


## 4. Memory Lifecycle and Waste

```
Phase                    Buffer    Scanner    SF (CSR)    Estimator
───────────────────────  ──────    ───────    ────────    ─────────
scan_and_buffer          ALLOC       —          —          —
FragmentRouter.scan()    READ      ALLOC      ALLOC*       —
  └ _to_np()             READ      PEAK**     ALLOC        —
  └ scanner GC'd         READ       FREE      LIVE         —
build_loci               —           —         READ        —
EB gDNA priors           —           —         READ        —
nrna_frac priors         —           —          —         READ/WRITE
batch_locus_em           —           —         READ       WRITE
_build_locus_meta        —           —          —         READ
return estimator         —           —        WASTE†     RETURN

* CSR allocated at end of scan, from scanner accumulator
** Both array.array AND numpy copies exist simultaneously
† ScoredFragments is no longer needed but persists in scope
```

**Key waste points:**

1. **Buffer** persists through the entire EM phase but is never read
   after `FragmentRouter.scan()` completes.  Holding 2 GB (or disk
   files on an SSD) for no reason.

2. **ScoredFragments** persists after `batch_locus_em()` returns but
   is never read again.  Holding GB of CSR data for no reason.

3. **Scanner accumulator → numpy duplication**: During `_to_np()`,
   the Python `array.array` and the new numpy array coexist.  For
   large datasets, this doubles ScoredFragments memory at peak.

4. **Temp directory**: Buffer spill files default to the system temp
   dir (e.g. `/var/folders/.../rigel_buf_XXXXX/chunk_NNNN.arrow`),
   which may be on a slow filesystem or have limited space.  Production
   environments need control over this location (fast SSD mount).


## 5. Proposed Changes

### Phase 5A: Free Buffer After Scan (Low Risk, Easy)

**Target: −2 GB at 24M frags, −2 GB at any size (capped by spill budget)**

In `quant_from_buffer()`, immediately after `em_data = builder.scan(buffer)`:

```python
em_data = builder.scan(buffer, log_every)
# Buffer is no longer needed — free in-memory chunks + disk files
buffer.cleanup()            # removes spilled Arrow files
buffer._chunks.clear()      # releases in-memory chunk references
buffer._memory_bytes = 0
gc.collect()                # prompt CPython to reclaim
```

This must be gated on the **annotated BAM** case — when
`annotations is not None`, the pipeline does a second BAM pass
*after* `quant_from_buffer` returns, but that second pass reads
the BAM file directly, not the buffer.  So early buffer cleanup is
safe in all cases.

The `run_pipeline()` `finally: buffer.cleanup()` call becomes a
no-op but should remain for safety.

**Expected savings:** In-memory chunks (~2 GiB budget) + disk space
from spilled chunks.

### Phase 5B: Free ScoredFragments After EM (Low Risk, Easy)

**Target: −1 to −5 GB depending on library size**

In `quant_from_buffer()`, immediately after the EM section ends
(after both the C++ batch path and the Python fallback path):

```python
# ScoredFragments CSR is no longer needed — free all arrays
del em_data
gc.collect()
```

The `em_data` variable is local to `quant_from_buffer` and is not
returned or stored on the estimator.

**Expected savings at various library sizes:**

| Frags (M) | ScoredFragments freed |
|---:|---:|
| 24 | ~1–5 GB |
| 50 | ~2 GB |
| 100 | ~4 GB |
| 200 | ~8 GB |

### Phase 5C: `--temp-dir` CLI Parameter (Low Risk, Easy)

**Target: User control over buffer spill location**

Currently `spill_dir` exists in `BamScanConfig` and `FragmentBuffer`
but is not exposed to the CLI.  The default is `None` → Python
`tempfile.mkdtemp()` → system temp (often `/var/folders/...` on macOS
or `/tmp` on Linux).

**Add `--temp-dir` to`rigel quant`:**

```
--temp-dir PATH   Directory for temporary buffer files (default: system temp).
                  Use a fast SSD mount in production environments.
```

Wire through:

1. `cli.py`: `parser.add_argument("--temp-dir", ...)`
2. `cli.py`: `scan=BamScanConfig(spill_dir=args.temp_dir, ...)`
3. `profiler.py`: Same parameter support
4. `BamScanConfig.spill_dir` already exists — just needs CLI wiring

Also add an info-level log message when spilling occurs, showing
the temp directory path, to aid debugging.

### Phase 5D: Eliminate Scanner Accumulator Duplication (Moderate Risk)

**Target: Halve peak memory during ScoredFragments construction**

The current flow:

```
array.array.append() × N_candidates   →   ~X GB in array.array
_to_np: np.frombuffer(arr, dtype).copy()  →   ~X GB in numpy
                                          Total: ~2X GB peak
```

**Option D1: Two-pass approach** (recommended)

1. **Counting pass**: Iterate buffer chunks once, counting candidates
   per unit (cheap — only needs `t_offsets` per chunk).  Compute
   `n_units` and `n_candidates`.
2. **Allocation**: Pre-allocate all numpy arrays at exact final size.
3. **Fill pass**: Iterate chunks again, writing directly into the
   pre-allocated numpy arrays at the correct offsets.

This eliminates `array.array` entirely.  Peak memory = 1× SF instead
of 2×.

**Cost**: Two iterations over buffer chunks instead of one.  For
spilled chunks this doubles I/O, but each chunk is only ~100 MB
compressed (LZ4) and sequential read is fast.  The second pass can
be optimized by computing candidate counts during the first pass
and pre-sizing the output.

Actually the two passes can be fused with the existing scan phases:
- Phase 1 (existing): `FragmentRouter.scan()` iterates chunks once
  but currently accumulates into `array.array`.  Instead, make this
  a **counting pass** that only computes `n_units`, per-unit candidate
  counts, and per-unit metadata (splice_type, is_spliced, etc.).
- Phase 2 (new): Pre-allocate numpy arrays.  Re-iterate chunks to
  fill per-candidate arrays (log_liks, coverage_weights, etc.)

**Option D2: In-place array.array → numpy conversion**

Replace `np.frombuffer(arr, dtype).copy()` with `np.frombuffer(arr,
dtype)` (zero-copy), then immediately `del arr`.  This avoids the
copy but leaves the numpy array backed by the `array.array`'s buffer.

**Risk**: `array.array` uses over-allocated buffers (2× growth factor),
so the zero-copy numpy array wastes up to 2× the actual data size.
Also, `np.frombuffer` creates a non-owning array that becomes invalid
if the source is garbage collected — need to transfer ownership
carefully.

**Option D3: C++ accumulator → numpy directly**

Replace the Python `array.array` accumulators with a C++ accumulator
(in `_scoring_impl` or a new module) that owns contiguous memory and
exposes it as numpy arrays via nanobind capsules.  This is the cleanest
solution but requires non-trivial C++ work.

**Recommendation**: Start with Option D1 (two-pass).  It's pure
Python, low risk, and eliminates the entire duplication problem.  If
I/O cost of the second chunk-read is unacceptable, upgrade to D3.

### Phase 5E: Type Narrowing (Moderate Risk, Moderate Impact)

**Target: −25% to −30% reduction in ScoredFragments and buffer size**

Several fields use wider types than needed:

#### ScoredFragments per-candidate

| Field | Current | Proposed | Savings/candidate | Risk |
|---|---|---|---|---|
| `log_liks` | float64 | float32 | 4 bytes | EM uses exp(log_lik); float32 has 7 decimal digits — likely sufficient for likelihood ratios |
| `coverage_weights` | float64 | float32 | 4 bytes | Weights are ≥1.0 and typically <10.0; float32 is more than adequate |

**Savings**: 8 bytes/candidate.  At 5 candidates/unit average:
- 24M fragments: ~300 MB saved (from 1.5 GB nominal)
- 200M fragments: ~2 GB saved (from 8.4 GB)

#### ScoredFragments per-unit

| Field | Current | Proposed | Savings/unit | Risk |
|---|---|---|---|---|
| `frag_ids` | int64 | int32 | 4 bytes | Max 2^31 ≈ 2.1B fragments before overflow; 200M is safe. Could use uint32 for 4.3B. |
| `gdna_log_liks` | float64 | float32 | 4 bytes | Same reasoning as `log_liks` |
| `offsets` | int64 | int32/uint32 | 4 bytes | Only safe if n_candidates < 2^31 (~2.1B). At 200M frags × 5 = 1B candidates, int32 works. At 400M+ frags with high candidate ratios, need int64. **Keep int64 for safety.** |

**Savings**: 8 bytes/unit (excluding offsets).
- 24M fragments: ~48 MB saved
- 200M fragments: ~400 MB saved

#### FragmentBuffer per-fragment

| Field | Current | Proposed | Savings/frag | Risk |
|---|---|---|---|---|
| `frag_id` | int64 | int32 | 4 bytes | Safe up to 2.1B fragments |
| `t_offsets` | int64 | int32 | 4 bytes | Safe if max CSR offset fits int32 (per chunk only — each chunk has ≤1M fragments × ~3 candidates = ~3M, well within int32) |

**Savings**: 8 bytes/fragment.
- 24M fragments: ~192 MB saved in buffer
- 200M fragments: ~1.6 GB saved in buffer

**Total type narrowing savings estimate (200M fragments):**
- ScoredFragments: ~2.4 GB
- Buffer: ~1.6 GB
- **Total: ~4 GB savings**

**Important**: The C++ `batch_locus_em()` function signatures must be
updated to accept the narrower types.  All callers in Python and C++
must be coordinated.  The EM solver internally should continue to use
float64 for numerical stability — only the storage types narrow.

### Phase 5F: Spill ScoredFragments to Memory-Mapped File (High Risk, Future)

**Target: Support 200M+ fragment libraries on 16 GB machines**

For very large libraries, even with all the above optimizations,
ScoredFragments may exceed available RAM.  A memory-mapped approach
would allow the OS to page CSR data in/out:

1. Write CSR arrays to a temporary file (or set of files) in the
   temp directory.
2. `mmap` the file and pass the pointer to C++ `batch_locus_em`.
3. The OS handles paging — frequently-accessed loci stay in RAM,
   cold regions get evicted to disk.
4. Delete the temporary file after EM completes.

**Advantages:**
- Memory usage bounded by working set, not total data size
- The C++ locus EM accesses data per-locus (good locality)
- No code-level streaming/batching complexity

**Disadvantages:**
- Random I/O on locus access (mitigated by OS page cache)
- Sorting loci by CSR offset would optimize sequential access
- File I/O adds ~10-30% wall time
- Platform-specific mmap behavior (Linux vs macOS)

**Recommendation**: Defer to Phase 6.  Only needed if users actually
encounter >16 GB machines running 200M+ fragment libraries.  The
`--temp-dir` parameter (Phase 5C) provides the prerequisite
infrastructure.


## 6. Implementation Priority

| Phase | Change | Savings | Risk | Effort |
|---|---|---|---|---|
| **5A** | Free buffer after scan | −2 GB (fixed) | Low | ~10 lines |
| **5B** | Free ScoredFragments after EM | −1 to −8 GB | Low | ~5 lines |
| **5C** | `--temp-dir` CLI parameter | Usability | Low | ~20 lines |
| **5D** | Two-pass scanner (eliminate duplication) | Halve peak | Moderate | ~200 lines |
| **5E** | Type narrowing | −25-30% CSR | Moderate | ~100 lines + C++ |
| **5F** | Memory-mapped ScoredFragments | Extreme scale | High | ~500 lines |

**Recommended execution order: 5A → 5B → 5C → 5E → 5D → 5F**

- 5A+5B are trivial, high-impact, and should be done immediately.
- 5C is a usability improvement needed for production deployment.
- 5E is straightforward type changes with significant savings.
- 5D is the most impactful for scaling but requires refactoring the
  scanner loop.
- 5F is only needed if scaling to extreme library sizes on constrained
  memory.

## 7. Projected Memory After Optimizations

### After 5A + 5B (Easy Wins)

| Frags (M) | Before | After | Savings |
|---:|---:|---:|---:|
| 24 | 20.2 GB | ~14.5 GB | ~5.7 GB |
| 50 | ~27 GB | ~21 GB | ~6 GB |
| 100 | ~37 GB | ~29 GB | ~8 GB |
| 200 | ~53 GB | ~39 GB | ~14 GB |

*Based on freeing ~2 GB buffer + ScoredFragments size.*

### After 5A + 5B + 5D (Halve Scanner Peak)

| Frags (M) | Before | After 5A+5B | After 5A+5B+5D | Total Savings |
|---:|---:|---:|---:|---:|
| 24 | 20.2 GB | ~14.5 GB | ~12 GB | ~8 GB |
| 50 | ~27 GB | ~21 GB | ~17 GB | ~10 GB |
| 100 | ~37 GB | ~29 GB | ~22 GB | ~15 GB |
| 200 | ~53 GB | ~39 GB | ~28 GB | ~25 GB |

*Eliminating the array.array duplication halves the CSR construction peak.*

### After 5A + 5B + 5D + 5E (Full Plan)

| Frags (M) | Original | Full Plan | Savings | Feasible on |
|---:|---:|---:|---:|---|
| 24 | 20.2 GB | ~9 GB | ~11 GB | 16 GB machine |
| 50 | ~27 GB | ~13 GB | ~14 GB | 16 GB machine |
| 100 | ~37 GB | ~19 GB | ~18 GB | 32 GB machine |
| 200 | ~53 GB | ~27 GB | ~26 GB | 32 GB machine |

*Type narrowing reduces CSR by ~25-30% on top of duplication elimination.*

## 8. Temp Directory Architecture

### Current State

```
spill_dir: None  →  tempfile.mkdtemp(prefix="rigel_buf_")
                →  /var/folders/ph/_lsypp750wj0q9_44ksdh3_00000gq/T/rigel_buf_6ysdvofr/
```

- Not user-controllable from CLI
- System temp dir may be slow, small, or shared
- No logging of the temp directory path
- No cleanup on crash (only `__del__` or context manager)

### Proposed Architecture

```
--temp-dir /fast_ssd/rigel_tmp
    └── rigel_buf_<random>/
        ├── chunk_0000.arrow   (LZ4 compressed, ~100 MB each)
        ├── chunk_0001.arrow
        └── ...
```

Requirements:
1. `--temp-dir` CLI parameter (default: system temp)
2. Info-level log when temp directory is created
3. Info-level log when each chunk is spilled (size, path)
4. `atexit` handler as safety net for cleanup on unclean exit
5. Validate directory exists and is writable at startup
6. Future: ScoredFragments mmap files go here too (Phase 5F)


## 9. Risk Analysis

| Risk | Mitigation |
|---|---|
| Freeing buffer too early when annotations need it | Annotations use a second BAM pass, not the buffer; verify code path |
| Freeing ScoredFragments while still referenced | `em_data` is local to `quant_from_buffer`; verify no references escape |
| float32 precision loss in EM | EM solver internally uses float64; only storage narrows; verify tolerance in tests |
| Two-pass scanner doubles chunk I/O | LZ4-compressed chunks are ~100 MB; sequential read is fast; re-verify with profiling |
| int32 frag_id overflow at >2.1B fragments | Use uint32 (4.3B) or guard with runtime check; 200M target is well within int32 |
| Memory-mapped CSR performance on macOS | macOS mmap has different page cache behavior than Linux; benchmark before committing |

## 10. Testing Strategy

1. **Add memory tracking to test suite**: A few integration tests that
   check RSS stays below expected thresholds for known-size inputs.
2. **Type narrowing parity tests**: Run existing 823 tests with
   narrowed types; verify numerical outputs within `rtol=1e-5`.
3. **Large-dataset smoke test**: Generate a 50M-fragment synthetic
   BAM and verify the pipeline completes under 16 GB RSS.
4. **Temp directory tests**: Verify `--temp-dir` creates files in the
   right place, cleanup happens on success, cleanup happens on error.
5. **Re-profile after each phase**: Confirm memory reduction matches
   projections.
