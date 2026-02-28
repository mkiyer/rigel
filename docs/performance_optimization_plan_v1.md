# hulkrna Performance Optimization Plan

**Date:** 2026-02-27  
**Benchmark:** PVT1 locus, 100k reads, oracle BAM, 379 transcripts, 35 genes  

---

## Current Profiling Summary (after P1 + P2)

**Total pipeline time: 43.4s** (scan_and_buffer + count_from_buffer)

| Stage | Time | % | Description |
|---|---|---|---|
| `scan_and_buffer` | 16.9s | **38.9%** | BAM parse → resolve → buffer |
| `EmDataBuilder.scan` | 12.8s | **29.5%** | Buffer → scored CSR candidates |
| Per-locus EM | 13.5s | **31.1%** | SQUAREM iterations + posterior assignment |
| Everything else | 0.2s | 0.5% | Index load, geometry, loci, priors |

### scan_and_buffer internals (16.9s)

| Function | Time | % of stage |
|---|---|---|
| `resolve_fragment` | 14.4s | **86.9%** |
| `Fragment.from_reads` | 0.8s | 4.9% |
| `parse_bam_file` | 0.8s | 4.7% |
| `buffer.append` | 0.6s | 3.5% |

### Per-locus EM internals (13.5s)

| Locus | Transcripts | Units | Components | EM time | Total |
|---|---|---|---|---|---|
| Locus 6 | 193 | 43,937 | 387 | 10.0s | 10.3s |
| Locus 2 | 136 | 41,331 | 273 | 2.9s | 3.0s |
| Locus 5 | 20 | 5,143 | 41 | 0.1s | 0.1s |

### cProfile Top-10 by self time (after P1 + P2)

| # | Function | Self time | Calls | Why expensive |
|---|---|---|---|---|
| 1 | `resolve_fragment` | 12.1s | 100k | Single-pass query + overlap computation |
| 2 | `query_exon_with_coords` | 5.9s | 187k | cgranges query per exon block |
| 3 | `_score_wta_mrna` | 5.0s | 91k | Vectorized mRNA scoring (batch NumPy) |
| 4 | `BiasProfile.effective_length` | 4.8s | 188k | Vectorized prefix-sum diff + .sum() |
| 5 | `numpy.ufunc.reduce` | 3.5s | 1.5M | NumPy internal reductions |
| 6 | `list.append` | 3.2s | 46.8M | CSR list accumulation |
| 7 | `builtins.max` | 3.0s | 31.4M | resolve_fragment inner loops |
| 8 | `_em_step` | 3.0s | 405 | Per-EC dense matrix operations |
| 9 | `_score_wta_nrna` | 3.0s | 51k | Vectorized nRNA scoring (batch NumPy) |
| 10 | `buffer.append` | 3.0s | 100k | Per-fragment list accumulation |

### Profiling History

**Original total pipeline time: 52.7s**

| Stage | Baseline | After P1 | After P1+P2 | Total savings |
|---|---|---|---|---|
| `scan_and_buffer` | 22.5s | 17.7s | 16.9s | **5.6s (25%)** |
| `EmDataBuilder.scan` | 15.8s | 15.6s | 12.8s | **3.0s (19%)** |
| Per-locus EM | 14.1s | 13.0s | 13.5s | 0.6s |
| **Total** | **52.7s** | **46.8s** | **43.4s** | **9.3s (17.6%)** |

cProfile totals: 323M → 309M → **227M** function calls (96M fewer, 30% reduction)

---

## Optimization Opportunities

### P1: Eliminate double cgranges query ✅ IMPLEMENTED

**Problem:** `resolve_fragment` queried cgranges TWICE per exon block:
1. `query_exon(block)` — returns (t_idx, g_idx, itype) for transcript set construction
2. `query_exon_with_coords(block)` via `compute_overlap_profile` — returns (t_idx, g_idx, itype, start, end) for overlap computation

That was **374k cgranges queries** (187k × 2) for 100k fragments, burning **11.7s** combined.

**Fix applied:** Merged into a SINGLE `query_exon_with_coords` per block inside `resolve_fragment`. Transcript sets and overlap profiles are now built in one pass. The standalone `compute_overlap_profile` call was removed from `resolve_fragment` (function retained for test compatibility).

**Measured results:**
- `resolve_fragment`: 19.4s → 14.7s (**4.7s saved, 24% reduction**)
- `scan_and_buffer`: 22.5s → 17.7s (**4.8s saved, 21% reduction**)
- Pipeline total: 52.7s → 46.8s (**5.9s saved, 11.2% speedup**)
- cgranges queries: 374k → 187k (50% reduction)
- Function calls: 323M → 309M (14.4M fewer)
- `query_exon` eliminated, `compute_overlap_profile` no longer called from hot path
- All 654 tests pass, benchmark MAE unchanged (15.86)

---

### P2: Replace per-candidate Python loops with vectorized ops ✅ IMPLEMENTED

**Problem:** `_score_wta_mrna` (3.9s), `_score_wta_nrna` (4.7s), `compute_coverage_weight` (6.4s), and `genomic_to_transcript_pos` (2.2s) were pure Python per-candidate loops called millions of times. The 73.5M `max()` and 55.6M `min()` calls were builtin-level overhead from scalar arithmetic in tight loops.

**Fixes applied:**

**Fix A — `compute_coverage_weight_batch`:** Rewrote the trapezoid coverage-weight model as a single vectorized NumPy function operating on arrays of (start, end, length). Eliminates 4.2M individual `compute_coverage_weight` calls and their 73.5M `max()`/55.6M `min()` scalar calls. Uses `np.clip`, `np.minimum`, `np.maximum`, and vectorized arithmetic.

**Fix B — `genomic_to_transcript_pos_bisect`:** Replaced the Python for-loop over exon intervals with `bisect.bisect_right` on pre-computed `_t_exon_data` (per-transcript exon starts/ends/cumulative offsets stored in `ScoringContext`). Reduces 1.8M calls from O(n_exons) Python loop to O(log n_exons) binary search. Total: 1.84M calls → 0.5s (was 2.2s).

**Fix C — Batch `_score_wta_mrna` / `_score_wta_nrna`:** Rewrote both functions to collect all candidates into NumPy arrays, then call `compute_coverage_weight_batch`, `frag_len_log_lik_batch`, and `genomic_to_transcript_pos_bisect` over the full batch. Scalar per-candidate Python loops replaced with vectorized NumPy operations.

Also added `frag_len_log_lik_batch` — vectorized fragment-length log-likelihood for arrays of fragment lengths.

**Measured results:**
- `EmDataBuilder.scan`: 15.6s → 12.8s (**2.8s saved, 18% reduction**)
- `compute_coverage_weight`: 6.4s/4.2M calls → 2.6s/142k calls (`compute_coverage_weight_batch`)
- `genomic_to_transcript_pos`: 2.2s/1.8M calls → 0.5s/1.8M calls (`genomic_to_transcript_pos_bisect`)
- `_score_wta_mrna`: 3.9s → 5.0s self (but cumulative 13.6→10.6s, includes batch callees)
- `_score_wta_nrna`: 4.7s → 3.0s self (cumulative 13.3→5.6s)
- `builtins.max`: 73.5M calls → 31.4M calls (57% reduction)
- `builtins.min`: 55.6M calls → 18.0M calls (68% reduction)
- Function calls: 309M → **227M** (82M fewer, 27% reduction)
- Pipeline total: 46.8s → **43.4s** (3.4s additional savings)
- All 654 tests pass, benchmark MAE unchanged (15.86)

---

### P3: Fast-path for uniform bias profiles (HIGH IMPACT)

**Problem:** `_apply_bias_correction` takes 11.3s cumulative. Under uniform bias (the current default):
- `effective_length(frag_len) = max(L - frag_len + 1, 1)` — a single subtraction — but the code does vectorized prefix-sum slicing + `.sum()` (188k calls, 4.2s)
- `fragment_weight(start, end) = 1.0` always — but the code clips, does prefix-sum lookup, and divides (3.9M calls, 1.7s + 2.9s numpy `.sum()`)

**Fix:** Add a uniform fast-path that uses closed-form arithmetic:
```python
def _apply_bias_correction_uniform(log_liks, t_indices, tx_starts, tx_ends,
                                    profile_lengths):
    """Optimized: uniform profiles → pure vectorized NumPy."""
    frag_lens = (tx_ends - tx_starts).astype(np.int64)
    eff_lens = np.maximum(profile_lengths[t_indices] - frag_lens + 1, 1)
    log_liks -= np.log(eff_lens.astype(np.float64))
    # fragment_weight is always 1.0 → skip entirely
```

Also: don't allocate `np.arange(L+1, float64)` prefix-sum arrays per transcript for the uniform case. Store only the length as an integer. This saves ~16 MB per locus.

**Estimated savings:** ~10s (11.3s → ~1s)

---

### P4: Optimize EM equivalence class operations (MEDIUM IMPACT)

**Problem:** `_em_step` (405 calls, 4.1s cumulative) iterates over equivalence classes doing per-class dense matrix operations. The two large loci (193t/387 components and 136t/273 components) dominate at 13.4s combined.

**Fix A — Profile and reduce equivalence class count:**
Currently groups by exact candidate tuple. With 387 components, many classes may have size 1 (overhead-dominated). Measure actual class distribution.

**Fix B — Fuse log-sum-exp:**
The inner loop does 4 separate NumPy calls per class: `np.add`, `.max(axis=1)`, `np.exp`, `.sum()`. Could combine into a single fused operation or use `scipy.special.logsumexp` which is optimized for this exact pattern.

**Fix C — Pre-compute `log_weights[comp_idx]` slicing:**
Currently indexes `log_weights` by `comp_idx` array every iteration. Pre-computing this as a contiguous subarray per class avoids repeated fancy indexing.

**Estimated savings:** 2-3s

---

### P5: Reduce buffer append overhead (LOW-MEDIUM IMPACT)

**Problem:** `buffer.append` (2.9s self, 5.1s cumulative) does 18 `list.append` per fragment + per-candidate appends. 60.5M `list.append` calls burn 4.0s globally.

**Fix A — Direct NumPy accumulation:**
Instead of Python list → convert to NumPy, write directly to pre-allocated NumPy arrays using an index cursor. The total fragment count is known (or can be bounded from BAM header).

**Fix B — Reduce per-candidate field count:**
`overlap_bp` stores 3 values per candidate (exon_bp, intron_bp, unambig_intron_bp). `unambig_intron_bp` is only used for nRNA init, which only needs unique-gene fragments. Could compute it lazily.

**Estimated savings:** 2-3s

---

### P6: Memory optimizations

**6a — Buffer memory (65.6 MB for 100k fragments = 656 bytes/frag):**
- `frag_lengths` (int32 per candidate) is only used during EmDataBuilder.scan scoring. Could compute on-the-fly from genomic coords + exon intervals instead of storing in buffer.
- `unambig_intron_bp` is only needed for unique-gene nRNA init accumulation. Could be accumulated during resolve_fragment and not stored per-candidate.

**6b — ScanData CSR (3.86M candidates × ~50 bytes = ~193 MB):**
- `count_cols` (uint8) encodes `splice_type * 2 + anti` — reconstructable from per-unit splice_type + per-candidate strand. Saves 3.86M bytes.
- `coverage_weights` could be float32 instead of float64, halving 31 MB.

**6c — BiasProfile uniform case:**
Currently allocates `np.arange(L+1, float64)` per transcript. For 379 transcripts totaling ~2M bases ≈ 16 MB of prefix sums that are literally `[0, 1, 2, ..., L]`. Store only the length integer for uniform profiles.

**6d — LocusEMInput duplication:**
`build_locus_em_data` re-extracts CSR arrays via fancy indexing (copies). The two large loci copy ~3.4M entries. Could use contiguous index ranges where possible to enable views.

---

### P7: Algorithmic simplifications

**7a — Lazy `compute_frag_lengths`:**
For unspliced fragments (no introns), all transcripts share the same genomic footprint. Skip the SJ query loop and return a single value. Currently 0.6s for 100k fragments.

**7b — Skip chimera detection when unnecessary:**
`_detect_intrachromosomal_chimera` (0.4s) runs on every single-ref fragment. Gate behind a config flag for oracle BAMs, or skip when all exon blocks share at least one transcript (quick frozenset intersection check).

**7c — Short-circuit nRNA scoring for zero-nRNA-evidence genes:**
`_score_wta_nrna` (4.7s, 51k calls) scores nRNA for every unspliced fragment. For genes with zero intronic evidence (common in this dataset where nRNA=0), the nRNA prior will be zeroed later anyway. Could gate on a pre-computed "gene has any intronic reads" flag.

**7d — Deduplicate gDNA log-likelihood computation:**
`_finalize_unit_metadata` recomputes strand/splice checks already done in `_score_wta_mrna`. Return gDNA ll as a byproduct of mRNA scoring.

---

## Recommended Execution Order

| Phase | Optimizations | Status | Measured/Expected savings |
|---|---|---|---|
| **Phase 1a** | P1 (merge double query) | ✅ Done | **5.9s saved (52.7→46.8s)** |
| **Phase 1b** | P2 (vectorize scoring) | ✅ Done | **3.4s saved (46.8→43.4s)** |
| **Phase 2** | P3 (uniform bias fast-path) | Pending | **~10s expected** |
| **Phase 3** | P4 (EM inner loop) + P5 (buffer append) | Pending | **~4-6s expected** |
| **Phase 4** | P6 (memory) + P7a/b/c/d (algorithmic) | Pending | Memory + ~2s |

### Projected timeline (revised after P1 + P2)

- **After P1:** 52.7s → **46.8s** (1.13× speedup) ✅ Measured
- **After P2:** 46.8s → **43.4s** (1.21× total speedup) ✅ Measured
- **After P3:** 43.4s → **~33s** (1.6× total speedup) — uniform bias fast-path
- **After P4+P5:** 33s → **~27s** (2.0× total speedup) — EM + buffer tuning
- **After P6+P7:** 27s → **~24s** (2.2× total speedup) + significant memory reduction

**Target: 52.7s → ~24-27s (2.0-2.2× speedup) with pure Python/NumPy, zero C/Cython.**

P3 (uniform bias fast-path) is the next highest-impact optimization,
targeting the 11.9s spent in `_apply_bias_correction` (4.8s in
`effective_length` + 2.9s in `fragment_weight` + numpy overhead)
that is pure overhead under the default uniform bias profiles.
