# Performance Optimization Plan — Real Data Profiling (2026-03-01)

## Dataset

| Property | Value |
|----------|-------|
| BAM | `mctp_LBX0069_SI_42153_HFFFMDRX7/bam/star.collate.markdup.bam` |
| Aligner | STAR |
| BAM size | 1.7 GB |
| Reference | GRCh38 + controls (286 contigs) |
| GTF | GENCODE genes + controls (254,461 transcripts, 63,472 genes) |
| Total BAM records | 21,595,608 |
| Unique fragments | 462,410 |
| Duplicates | 19,180,177 (88.8%) |
| Secondary alignments | 1,302,186 |
| Supplementary | 67,176 |
| Multimapping | 217,404 |
| Chimeric | 4,902 |
| EM loci | 18,136 |

## Baseline Performance

| Metric | Value |
|--------|-------|
| **Wall time** | **133.1s** |
| Index load time | 50.1s |
| Throughput | 162,197 BAM records/sec |
| Total function calls | 294,272,680 |

## Phase Breakdown

| Phase | Self-time | Cumul-time | % of Wall |
|-------|-----------|------------|-----------|
| BAM parsing (`parse_bam_file`) | 34.1s | 65.4s | 49.1% |
| EM solver (`_em_step`) | 14.0s | 22.6s | 17.0% |
| Stats proxy (`__setitem__` + `__getitem__`) | 20.0s | 25.7s | 19.3% |
| Locus construction (`build_locus_em_data`) | 1.6s | 12.9s | 9.7% |
| BiasProfile / `np.arange` | — | 7.5s + 1.8s | 7.0% |
| Fragment resolution (`resolve_fragment`) | 3.4s | 3.5s | 2.6% |
| Fragment construction (`from_reads`) | 2.1s | 5.7s | 4.3% |
| Scan (`_add_single_fragment` + `scan`) | 2.6s | 8.0s | 6.0% |
| Buffer (`__getitem__`) | 1.1s | 1.2s | 0.9% |
| Locus EM driver (`run_locus_em`) | 1.4s | 25.7s | 19.3% |

## Top 10 Hotspots by Self-Time

| Rank | Function | Self (s) | Calls | Notes |
|------|----------|----------|-------|-------|
| 1 | `bam.parse_bam_file` | 34.1 | 462K | groupby + per-record filter loop |
| 2 | `estimator._em_step` | 14.0 | 278K | EM inner loop (numpy) |
| 3 | `stats.__setitem__` | 14.0 | 44M | `_BamStatsProxy` via setattr per BAM record |
| 4 | `{numpy.arange}` | 7.5 | 361K | `BiasProfile.uniform()` per transcript per locus |
| 5 | `{ufunc.reduce}` | 6.7 | 9.2M | numpy sum/max/any inside EM |
| 6 | `stats.__getitem__` | 6.0 | 44M | `_BamStatsProxy` via getattr per BAM record |
| 7 | `bam.<lambda>` (groupby key) | 3.5 | 21.6M | `lambda r: r.query_name` |
| 8 | `resolution.resolve_fragment` | 3.4 | 462K | C++ dispatch (already native) |
| 9 | `{builtins.setattr}` | 3.1 | 44M | called by stats proxy |
| 10 | `{builtins.getattr}` | 2.8 | 45.5M | called by stats proxy |

## Key Findings

### 1. Stats Proxy is the #1 Surprise (~20s, 15% of wall)

The `_BamStatsProxy` class wraps a `PipelineStats` dataclass to give
`parse_bam_file()` a dict-like interface. Every `stats['total'] += 1`
triggers **both** `__getitem__` (→ `getattr`) and `__setitem__` (→
`setattr`). With 21.6M BAM records, this amounts to 44M proxy calls
costing ~20s total (`__setitem__` 14.0s + `__getitem__` 6.0s +
`setattr` 3.1s + `getattr` 2.8s).

**Fix**: Replace the proxy with direct attribute access on the dataclass
inside `parse_bam_file`. Use local variables for hot counters
(`total`, `duplicate`, `secondary`, `supplementary`, `qc_fail`,
`unmapped`) and flush to the dataclass once at the end.

**Expected savings**: ~20s → ~0.5s (local int increments).

### 2. BAM Parsing Inner Loop (34.1s self, 49% of wall)

`parse_bam_file` iterates 21.6M pysam records in Python with:
- `groupby(bam_iter, key=lambda r: r.query_name)` — the lambda alone
  costs 3.5s (21.6M calls)
- per-record flag checks (`is_qcfail`, `is_unmapped`, `is_duplicate`,
  `is_secondary`, `is_supplementary`, `has_tag('NH')`, `get_tag('NH')`)
- Python list.append per usable record

The 34.1s self-time in `parse_bam_file` is dominated by pysam attribute
access overhead — each flag check crosses the C ↔ Python boundary.
The duplicate rate is 88.8% which means most records are touched
only to be discarded.

**Fix (Phase F)**: Move the BAM record filtering and grouping loop into
C++ using htslib directly. A C++ function would stream BAM records,
filter, group by query name, and return only usable read groups. This
eliminates 21.6M Python-level pysam attribute accesses.

**Expected savings**: 34s + 3.5s lambda → ~5s in C++ (~32s saved).

### 3. BiasProfile.uniform() Allocation (7.5s + 1.8s, 7%)

`BiasProfile.uniform(length)` is called ~307K times (once per
transcript-component per locus), each creating a new
`np.arange(length + 1, dtype=float64)` array. These are immediately
consumed by `_apply_bias_correction` which already detects the uniform
fast-path.

**Fix**: Since **all** profiles are currently uniform, skip the
`BiasProfile` object creation entirely. Instead, pass a flat
`np.ndarray` of transcript lengths to `_apply_bias_correction_uniform`
directly, avoiding 307K numpy array allocations.

**Expected savings**: ~9s → ~0.1s.

### 4. EM Solver (_em_step) (14.0s self, 278K calls)

Each `_em_step` call performs numpy operations on small per-locus
matrices. With 18,136 loci and SQUAREM acceleration (3 EM steps per
SQUAREM iteration), this creates 278K calls. The per-call overhead of
numpy on small arrays (often < 100 elements) dominates.

**Fix (Phase E)**: Move `_em_step` / `_vbem_step` into C++ operating
on pre-packed per-locus data structures. Eliminates Python loop and
numpy per-call overhead for small arrays.

**Expected savings**: 14s + 6.7s numpy reduce → ~2s in C++ (~18s saved).

### 5. Locus Construction (1.6s self + subcalls → 12.9s cumulative)

`build_locus_em_data` takes 12.9s cumulative across 18,136 loci.
The 7.5s `np.arange` cost (item 3 above) is the largest subcall.
After fixing the BiasProfile allocation, the remaining ~4s is
vectorized numpy work (fancy indexing, dedup, sorting) that could
move to C++ alongside the EM solver.

### 6. Fragment Resolution (3.4s, already native)

`resolve_fragment` at 3.4s for 462K calls is already C++ (via
`_resolve_impl`). At 7.3 µs/fragment, this is acceptable.

### 7. Fragment Construction (2.1s) + BAM grouping (1.0s)

`Fragment.from_reads` at 2.1s and `_group_records_by_hit` at 1.0s
are Python-side fragment construction from pysam records. These
would be subsumed by a C++ BAM parsing phase (Phase F).

## Optimization Phases (Prioritized by ROI)

### Phase G: Eliminate Stats Proxy Overhead
- **Target**: `_BamStatsProxy.__setitem__/__getitem__` (20s)
- **Approach**: Replace proxy with local counter variables in
  `parse_bam_file`; flush to `PipelineStats` dataclass once at end
- **Complexity**: Low (1-2 hours)
- **Expected savings**: ~20s
- **New wall time**: ~113s

### Phase H: Eliminate BiasProfile Allocation
- **Target**: `BiasProfile.uniform()` + `np.arange` (9.3s)
- **Approach**: Pass transcript-length arrays directly to
  `_apply_bias_correction_uniform` instead of creating BiasProfile
  objects; skip `all(p.is_uniform ...)` check
- **Complexity**: Low-medium (2-4 hours)
- **Expected savings**: ~9s
- **New wall time**: ~104s

### Phase E: C++ EM Kernel
- **Target**: `_em_step` (14.0s) + numpy reduce overhead (6.7s)
- **Approach**: C++ implementation of `_em_step` / `_vbem_step`
  operating on packed per-locus data. Include `_build_equiv_classes`
  and bias correction in the C++ path.
- **Complexity**: High (1-2 days)
- **Expected savings**: ~18s
- **New wall time**: ~86s

### Phase F: C++ BAM Parsing
- **Target**: `parse_bam_file` (34.1s) + lambda (3.5s) +
  `from_reads` (2.1s) + `_group_records_by_hit` (1.0s)
- **Approach**: C++ htslib-based BAM reader that filters, groups by
  query name, and returns structured fragment data directly. Bypasses
  21.6M pysam Python attribute accesses.
- **Complexity**: Very high (2-4 days)
- **Expected savings**: ~35s
- **New wall time**: ~51s

### Phase I: C++ Locus Construction
- **Target**: `build_locus_em_data` residual (4s after Phase H)
- **Approach**: Fuse with Phase E C++ EM kernel — extract per-locus
  subproblem, build equiv classes, and run EM in one C++ call
- **Complexity**: Medium (absorbed into Phase E)
- **Expected savings**: ~3s
- **New wall time**: ~48s

## Projected Performance Timeline

| Phase | Wall Time | Savings | Cumulative Speedup |
|-------|-----------|---------|-------------------|
| Baseline (current) | 133s | — | 1.0× |
| Phase G (stats proxy) | 113s | 20s | 1.2× |
| Phase H (bias alloc) | 104s | 9s | 1.3× |
| Phase E (C++ EM) | 86s | 18s | 1.5× |
| Phase F (C++ BAM) | 51s | 35s | 2.6× |
| Phase I (C++ locus) | 48s | 3s | 2.8× |

## Recommendations

1. **Start with Phase G** — pure Python fix, 20s savings, < 2 hours
   of work. Highest ROI.
2. **Phase H next** — another easy Python-level fix, 9s savings.
3. **Phase E** for the EM solver is the highest-impact C++ task.
4. **Phase F** (BAM parsing) is the largest single opportunity (35s)
   but also the most complex, requiring direct htslib integration.

After all phases, the projected wall time is ~48s for a 1.7GB
BAM file with 21.6M records — a 2.8× speedup from current baseline.

## Notes

- This profiling was run on the full BAM (no fragment limit), unlike
  the previous 200K-fragment simulated data profiles.
- The 88.8% duplicate rate means the BAM parsing loop processes
  ~19M records that are immediately discarded. C++ BAM parsing
  (Phase F) would short-circuit these at the htslib level.
- Index loading takes 50.1s (not included in pipeline wall time).
  This is dominated by cgranges interval tree construction from
  4M genomic intervals and is a separate optimization target.
- The real-data profile reveals that `_BamStatsProxy` overhead
  was invisible in simulated data (200K fragments vs 21.6M records)
  because the proxy cost scales with BAM records, not fragment count.
