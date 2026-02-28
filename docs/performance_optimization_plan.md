# hulkrna Performance Optimization Plan

**Date:** 2026-02-28
**Benchmark:** PVT1 locus, 100k reads, oracle BAM, 379 transcripts, 35 genes
**Conditions:** gDNA 10% + nRNA 10% (realistic mixed-pool scenario)

---

## Current Profiling Summary (after P1 + P2 + collapsed index + P3)

**Total pipeline time: 16.1s** (profiler stage timings)
**cProfile total: 31.3s** (includes profiler overhead from double scan_and_buffer)

| Stage | Time | % | Description |
|---|---|---|---|
| `scan_and_buffer` | 7.0s | **43.4%** | BAM parse → resolve → buffer |
| `EmDataBuilder.scan` | 7.7s | **48.1%** | Buffer → scored CSR candidates |
| Per-locus EM | 1.0s | **6.5%** | SQUAREM iterations + posterior assignment |
| Everything else | 0.3s | 2.0% | Index load, geometry, loci, priors |

### scan_and_buffer internals (7.0s)

| Function | Time | % of stage |
|---|---|---|
| `resolve_fragment` | 4.6s | **73.0%** |
| `parse_bam_file` | 0.7s | 11.3% |
| `Fragment.from_reads` | 0.6s | 9.7% |
| `buffer.append` | 0.4s | 6.0% |

### EmDataBuilder.scan internals (7.7s)

| Function | Self time | Calls | Notes |
|---|---|---|---|
| `_score_wta_nrna` | 3.36s | 80,167 | Dominant due to nRNA 10%; absent in nRNA=0 runs |
| `compute_coverage_weight_batch` | 1.60s | 89,763 | Called by both mRNA + nRNA |
| `_emit_nrna` | 0.93s | 80,167 | CSR list appends for nRNA candidates |
| `_score_wta_mrna` | 0.40s | 80,323 | Already vectorized |
| `_finalize_unit_metadata` | 0.10s | 80,323 | |

### Per-locus EM internals (1.0s)

| Locus | Transcripts | Units | Components | Build | EM | Assign | Total |
|---|---|---|---|---|---|---|---|
| Locus 8 | 195 | 20,740 | 391 | 0.166s | 0.499s | 0.024s | **0.689s** |
| Locus 4 | 141 | 42,199 | 283 | 0.046s | 0.224s | 0.007s | **0.278s** |
| Other 10 | ≤20 | ≤7,173 | ≤41 | ~0.01s | ~0.06s | ~0.001s | 0.07s |

### cProfile Top-12 by self time (31.3s total)

| # | Function | Self time | Calls | Description |
|---|---|---|---|---|
| 1 | `resolve_fragment` | 5.63s | 100k | Single-pass query + overlap |
| 2 | `_score_wta_nrna` | 3.36s | 80k | Vectorized nRNA scoring |
| 3 | `buffer.append` | 1.87s | 80k | Per-fragment list accumulation |
| 4 | `compute_coverage_weight_batch` | 1.60s | 90k | Trapezoid model, already vectorized |
| 5 | `list.append` | 1.31s | 19M | CSR list accumulation |
| 6 | `array.array.append` | 1.19s | 13M | CSR array accumulation |
| 7 | `dict.get` | 0.99s | 13.7M | From resolve_fragment inner loops |
| 8 | `_emit_nrna` | 0.93s | 80k | nRNA CSR list accumulation |
| 9 | `builtins.max` | 0.67s | 6.1M | Scalar max in inner loops |
| 10 | `numpy.ufunc.reduce` | 0.59s | 651k | NumPy internal reductions |
| 11 | `dict.setdefault` | 0.59s | 6.4M | resolve_fragment inner loops |
| 12 | `_em_step` | 0.51s | 216 | Per-EC dense matrix operations |

### Total function calls: 100M (down from 113M pre-P3, 323M baseline)
### Total pipeline (profiler): 16.1s (down from 23.6s pre-P3, 52.7s baseline)

---

## Profiling History

| Stage | Baseline | After P1 | After P1+P2 | +Collapsed idx | **+P3** |
|---|---|---|---|---|---|
| `scan_and_buffer` | 22.5s | 17.7s | 16.9s | 6.7s | **7.0s** |
| `EmDataBuilder.scan` | 15.8s | 15.6s | 12.8s | 7.3s | **7.7s** |
| Per-locus EM | 14.1s | 13.0s | 13.5s | 9.3s | **1.0s** |
| **Total** | **52.7s** | **46.8s** | **43.4s** | **23.6s** | **16.1s** |

| Metric | Baseline | Current | Reduction |
|---|---|---|---|
| Pipeline time | 52.7s | **16.1s** | **69.4%** |
| Function calls | 323M | **100M** | **69.0%** |
| Speedup | 1.0× | **3.3×** | — |

> **P3 impact:** Per-locus EM collapsed from 9.3s → 1.0s (89% reduction).
> `BiasProfile.effective_length` (4.84s, 183k calls) and
> `fragment_weight` (0.99s, 2.14M calls) were completely eliminated
> from the hot path by the vectorized uniform fast-path. The
> `_apply_bias_correction_uniform` function does not even appear in
> the top-60 cProfile entries — it executes in ~1ms total.

---

## Completed Optimizations

### P1: Eliminate double cgranges query ✅

**Merged** two cgranges queries per exon block (transcript set + overlap)
into a single `query_exon_with_coords` call.

- `resolve_fragment`: 19.4s → 14.7s
- Pipeline total: 52.7s → **46.8s** (5.9s saved, 11.2% speedup)
- All 654 tests pass, MAE unchanged

### P2: Vectorize per-candidate scoring ✅

**Replaced** Python per-candidate loops with batch NumPy operations in
`_score_wta_mrna`, `_score_wta_nrna`, `compute_coverage_weight_batch`,
and `genomic_to_transcript_pos_bisect`.

- `EmDataBuilder.scan`: 15.6s → 12.8s
- Function calls: 309M → **227M** (82M fewer)
- Pipeline total: 46.8s → **43.4s** (3.4s saved)
- All 654 tests pass, MAE unchanged

### Collapsed interval index ✅

**Merged** EXON+INTRON cgranges into a single collapsed index with one
`query()` call per exon block returning `(start, end, itype, t_set)`.

- `resolve_fragment`: 14.4s → 4.4s
- Function calls: 227M → **113M** (50% reduction)
- Pipeline total: 43.4s → **23.6s** (19.8s saved)
- All 654 tests pass, MAE unchanged

### P3: Uniform bias fast-path ✅

**Added** `is_uniform` flag to `BiasProfile` and a vectorized
`_apply_bias_correction_uniform` that computes
`−log(max(L − frag_len + 1, 1))` in 3 lines of NumPy, completely
bypassing the per-pair LUT loop, `effective_length`, and
`fragment_weight` calls.

- Per-locus EM: 9.3s → **1.0s** (89% reduction)
  - Locus 8: 6.74s → 0.69s
  - Locus 4: 2.32s → 0.28s
- `BiasProfile.effective_length`: 4.84s → **0s** (eliminated)
- `BiasProfile.fragment_weight`: 0.99s → **0s** (eliminated)
- `numpy.ufunc.reduce`: 3.39s → 0.59s (83% reduction)
- Function calls: 113M → **100M** (13M fewer)
- Pipeline total: 23.6s → **16.1s** (7.5s saved, 32% reduction)
- All 663 tests pass (9 new tests added), MAE unchanged

---

## Remaining Optimization Opportunities

### P4: Reduce nRNA scoring overhead ★★★ HIGH IMPACT (now #1)

**Problem:** `_score_wta_nrna` is now the **#2 self-time function** at
**3.36s** (7.31s cumulative including `_emit_nrna` + coverage weight).
With per-locus EM collapsed to 1.0s, nRNA scoring is the dominant
remaining bottleneck in the EM-builder path.

It's called for every buffered fragment (80k), scoring nRNA candidates
even when:
- The transcript is single-exon (always filtered out → wasted computation)
- The locus has zero intronic evidence (nRNA prior will be zeroed)

**Fix A — Short-circuit single-exon transcripts before vectorization:**

Currently the single-exon filter happens inside `_score_wta_nrna` after
building arrays. Check it before allocation:
```python
multi_exon_mask = self.ctx.multi_exon_mask[t_inds]
if not np.any(multi_exon_mask):
    return {}
```

**Fix B — Skip nRNA scoring for genes with zero intronic evidence:**

Pre-compute a per-gene "has intronic evidence" flag during the scan pass.
Skip `_score_wta_nrna` for fragments in genes where intronic evidence is
zero. Many genes in the PVT1 benchmark have zero nRNA evidence.

**Fix C — Stop emitting nRNA candidates for zero-prior components:**

nRNA components with zero prior are zeroed in the EM anyway. Don't emit
them in `_emit_nrna`. This reduces CSR size and EM candidate count.

**Expected savings:** ~2-4s (15-25% of total)

---

### P5: Reduce resolve_fragment inner-loop overhead ★★ MEDIUM-HIGH IMPACT

**Problem:** `resolve_fragment` is the **#1 self-time function** at
**5.63s**. The inner loop processes each cgranges hit, clipping intervals
and accumulating overlap BP per transcript. Key overhead sources:

- `dict.get`: 13.7M calls, 0.99s (overlap_bp lookups)
- `builtins.max`: 6.1M calls, 0.67s (interval clipping)
- `dict.setdefault`: 6.4M calls, 0.59s
- `_subtract_intervals`: 3.3M calls, 0.41s
- `set.update`: 3.0M calls, 0.33s

**Fix A — Use arrays instead of dicts for overlap accumulation:**

Transcript indices are small integers (0-378). Use fixed-size NumPy
arrays indexed by `t_idx` instead of `defaultdict(int)`, eliminating
`dict.setdefault`, `dict.get` overhead (~1.6s).

**Fix B — Vectorize interval clipping:**

Batch `max()`/`min()` calls into NumPy operations when `index.query()`
returns multiple hits (~0.7s).

**Fix C — Defer unambig_intron_bp computation:**

`_subtract_intervals` (3.3M calls, 0.41s) is only needed for nRNA init
on unique-gene fragments. Skip for multi-gene fragments.

**Expected savings:** ~2-3s (12-19% of total)

---

### P6: Reduce buffer/CSR append overhead ★★ MEDIUM IMPACT

**Problem:** Buffer and CSR list appends consume substantial time:

- `buffer.append` (AccumulatorChunk): 1.87s, 80k calls
- `list.append`: 19M calls, 1.31s
- `array.array.append`: 13M calls, 1.19s
- `_emit_nrna`: 0.93s, 80k calls (Python loop over winners)

**Fix A — Pre-allocate NumPy arrays for buffer accumulation:**

Replace Python lists with pre-allocated NumPy arrays and a cursor index.

**Fix B — Use `extend()` for batch CSR emission:**

Collect winners into temporary lists and `extend()` once per fragment.

**Expected savings:** ~1-2s

---

### P7: Coverage weight optimization ★ LOW-MEDIUM IMPACT

**Problem:** `compute_coverage_weight_batch` takes 1.60s across 90k calls.
Small-batch NumPy overhead dominates (mean ~25 candidates per call).

**Fix — Fuse operations and scalar fallback for small batches:**

For batches < 32, use a scalar loop to avoid 6× array allocation overhead.

**Expected savings:** ~0.5-1s

---

### P8: EM inner loop ★ LOW IMPACT (collapsed by P3)

With P3 implemented, EM itself takes only 1.0s total. The remaining cost
is dominated by `_build_equiv_classes` (0.14s) and `_em_step` (0.71s).
Further optimization here has diminishing returns.

**Fix — Vectorize equivalence class construction:**

Replace Python dict grouping with `np.unique` on compound keys.

**Expected savings:** ~0.1-0.3s

---

### P9: Memory optimizations ★ LOW IMPACT (aids scalability)

**9a — BiasProfile uniform case:**
Don't allocate `np.arange(L+1, float64)` prefix-sum arrays for uniform.
Store only integer length. Saves ~16 MB per locus (PVT1), ~100+ MB (chr22).

**9b — ScanData CSR float32 for coverage_weights.**

**9c — Buffer field consolidation.**

---

## Recommended Execution Order

| Phase | Optimizations | Status | Savings |
|---|---|---|---|
| **Phase 1a** | P1 (merge cgranges query) | ✅ Done | **5.9s (52.7→46.8s)** |
| **Phase 1b** | P2 (vectorize scoring) | ✅ Done | **3.4s (46.8→43.4s)** |
| **Phase 1c** | Collapsed interval index | ✅ Done | **19.8s (43.4→23.6s)** |
| **Phase 2** | P3 (uniform bias fast-path) | ✅ Done | **7.5s (23.6→16.1s)** |
| **Phase 3** | P4 (nRNA gating) + P5 (resolve inner loops) | 🔲 Pending | **~4-7s (16→~10s)** |
| **Phase 4** | P6 (buffer/CSR) + P7 (cov weight) | 🔲 Pending | **~2-3s (10→~7s)** |
| **Phase 5** | P8 (EM) + P9 (memory) | 🔲 Pending | ~0.5s + memory |

### Projected timeline

| Milestone | Total time | Speedup vs original |
|---|---|---|
| Baseline | 52.7s | 1.0× |
| After P1+P2 | 43.4s | 1.2× |
| After collapsed index | 23.6s | 2.2× |
| **After P3 (current)** | **16.1s** | **3.3×** |
| After P4+P5 | ~**10s** | **5.3×** |
| After all (P4-P7) | ~**7s** | **7.5×** |

**Next target: P4 (nRNA scoring gating) + P5 (resolve_fragment inner loops).**

After P3, the pipeline is dominated by two stages of roughly equal cost:
`scan_and_buffer` (7.0s) and `EmDataBuilder.scan` (7.7s). P4 targets the
scan stage (3.36s in nRNA scoring), P5 targets the resolve stage (5.63s).
Together they could bring the total under 10s.

---

## cProfile comparison: Before vs After P3

| Function | Before P3 | After P3 | Change |
|---|---|---|---|
| `BiasProfile.effective_length` | 4.84s (183k calls) | **0s** | **eliminated** |
| `numpy.ufunc.reduce` | 3.39s (834k calls) | 0.59s (651k calls) | **-83%** |
| `BiasProfile.fragment_weight` | 0.99s (2.14M calls) | **0s** | **eliminated** |
| `_apply_bias_correction` (self) | 1.06s | **~0s** | **eliminated** |
| `_em_step` | 0.50s | 0.51s | unchanged |
| `_build_equiv_classes` | 0.10s | 0.10s | unchanged |
| Per-locus EM (wall time) | 9.30s | **1.04s** | **-89%** |
| **Pipeline total** | **23.6s** | **16.1s** | **-32%** |
