# hulkrna Performance Optimization Plan

## Executive Summary

hulkrna is **14× slower** than salmon/kallisto (mean 7.2s vs 0.5s per 50k-fragment
region).  Detailed profiling reveals two dominant bottlenecks that account for
**85%** of wall time:

| Bottleneck | Time (chr17) | % of Total | Root Cause |
|---|---|---|---|
| Per-locus EM | 14.1s | 63% | 48K units processed individually instead of 169 equivalence classes |
| `resolve_fragment` | 5.1s | 22% | Pure-Python per-exon-block overlap loop |
| `_scan_and_build_em_data` | 1.9s | 8% | Python list.append CSR building (22.8M appends) |
| `buffer.append` | 1.9s | 8% | Per-fragment attribute copying |

The EM bottleneck alone contributes **63%** of total time on the slowest region
(chr17).  The fix in Priority 1 below is expected to deliver a **~20–50×** EM
speedup, cutting chr17 from 22.6s to **~5–7s** with just algorithmic changes
in pure Python + numpy.

---

## Profiling Data

### Baseline Timing (50k fragments/region, oracle aligner, pristine)

| Region | hulkrna | salmon | kallisto | hulkrna/salmon |
|---|---|---|---|---|
| chr17_43044295_43170245 | 22.55s | 1.17s | 0.55s | 19.3× |
| chr10_121476640_121601584 | 15.02s | 1.23s | 0.52s | 12.2× |
| chr11_35138870_35232402 | 7.06s | 0.40s | 0.67s | 17.7× |
| chr7_55019017_55211628 | 4.89s | 0.34s | 0.40s | 14.4× |
| **Mean** | **7.22s** | **0.52s** | **0.49s** | **13.9×** |

### chr17 Stage Breakdown (worst case: 22.6s)

```
  scan_and_buffer                6.381s   28.3%
    parse_bam_file:                0.408s   6.7%
    Fragment.from_reads:           0.425s   6.9%
    resolve_fragment:              5.051s  82.5%   ← Bottleneck #2
    buffer.append:                 0.242s   3.9%
  scan_and_build_em_data         1.862s    8.2%
  Per-locus EM                  14.253s   63.2%   ← Bottleneck #1
    Build locus EM data:           0.128s
    Run locus EM:                 14.087s          ← 93% of EM stage
    Assign posteriors:             0.038s
```

### chr11_5225000 (fast region: 2.4s, 24 transcripts, 7 loci)

```
  scan_and_buffer                1.808s   76.3%
    resolve_fragment:              0.878s  55.2%
  scan_and_build_em_data         0.440s   18.6%
  Per-locus EM                   0.065s    2.8%
```

### cProfile Hot Functions (chr17, sorted by tottime)

| Function | tottime | cumtime | Calls | Notes |
|---|---|---|---|---|
| `run_locus_em` | 6.51s | 10.86s | 1 | EM inner loop |
| `compute_overlap_profile` | 2.18s | 5.85s | 50K | Per-fragment overlap |
| `np.ufunc.at` (add.at) | 1.92s | 1.92s | 1007 | EM M-step scatter |
| `list.append` | 1.56s | 1.56s | 22.8M | CSR building |
| `index.query_exon` | 1.49s | 1.72s | 112K | cgranges query |
| `index.query_exon_with_coords` | 1.45s | 1.67s | 112K | With coords |
| `ndarray.repeat` | 1.35s | 1.35s | 2009 | EM E-step expand |
| `_add_transcript_candidates` | 1.24s | 1.55s | 48K | mRNA scoring |
| `buffer._AccumulatorChunk.append` | 1.11s | 1.94s | 50K | Columnar append |
| `ufunc.reduceat` | 1.08s | 1.08s | 2004 | EM E-step reduce |
| `_add_nrna_candidates` | 0.91s | 1.05s | 26K | nRNA scoring |

---

## Key Findings

### Finding 1: EM Processes Individual Units Instead of Equivalence Classes

The chr17 mega-locus has:
- **48,006 units** with 1,679,800 candidate entries (35 cand/unit)
- But only **169 unique candidate-set patterns** — a **284× compression ratio**
- Top pattern: 15,205 units all share the same 48 candidate set

The current EM iterates over all 1.68M entries per iteration using CSR
`reduceat`/`repeat`/`add.at` operations.  With equivalence-class grouping,
each class becomes a single dense matrix operation, reducing overhead from
48K individual units to 169 batched operations.

### Finding 2: EM Fails to Converge

The EM hits the 1000-iteration cap with delta=1.98e-5 (threshold=1e-6).
Convergence profile:
- iter 50: delta = 4.0e-3
- iter 100: delta = 2.1e-3
- iter 200: delta = 3.7e-4
- iter 500: delta = 3.9e-5
- iter 1000: delta = 2.0e-5 (still above threshold)

The linear convergence rate is slow.  At the current rate, convergence to 1e-6
would require ~3000+ iterations.

### Finding 3: resolve_fragment is a Pure-Python Per-Overlap Loop

`compute_overlap_profile` iterates over each fragment exon block (typically 2
for paired-end), queries cgranges for overlapping intervals, and accumulates
per-transcript overlap in Python dicts via `defaultdict(int)`.  For chr17 with
52 transcripts, each query returns many hits, producing millions of Python dict
operations.

### Finding 4: CSR Building Uses 22.8M list.append Calls

`_scan_and_build_em_data` builds CSR arrays via Python `list.append()` for
`t_indices_list`, `log_liks_list`, and `count_cols_list`.  The 22.8M calls
consume 1.56s.

---

## Optimization Plan

### Priority 1: Equivalence-Class EM (Expected: 20–50× EM speedup)

**Impact:** Reduce EM from 14.1s → ~0.3–0.7s on chr17  
**Effort:** Medium (pure Python/numpy, ~200 lines)  
**Risk:** Low (mathematically equivalent)

**Strategy:**  After building per-locus EM data, group units by their
candidate index set.  For each equivalence class:

```python
# For class c with n_c units and candidates S_c = {c1, ..., ck}:
# log_lik_matrix: (n_c, k_c) — each row is one unit's log-likelihoods
# log_weights: (k_c,) = log(theta[S_c]) - log(eff_len[S_c])

log_post = log_lik_matrix + log_weights  # broadcast (n_c, k_c)
max_per_row = log_post.max(axis=1, keepdims=True)
post = np.exp(log_post - max_per_row)
post /= post.sum(axis=1, keepdims=True)
em_totals[S_c] += post.sum(axis=0)  # column sums
```

**Concrete steps:**
1. Add `_build_equiv_classes()` helper that groups local EM units by sorted
   candidate tuple → returns list of `(candidate_indices, log_lik_matrix)`.
2. Modify `run_locus_em()` to accept equiv-class data structure and iterate
   over classes instead of the flat CSR.
3. Each class processes a dense `(n_units_in_class × n_candidates)` matrix — 
   numpy can vectorize this efficiently.

**Estimated compression:** 169 classes from 48K units (chr17).  Each iteration
touches the same 1.68M candidate entries but in 169 cache-friendly blocks
instead of 48K tiny segments.  The `repeat`/`reduceat`/`add.at` overhead is
eliminated.

### Priority 2: EM Convergence Acceleration (Expected: 2–5× fewer iterations)

**Impact:** Reduce iterations from 1000 → 200–500  
**Effort:** Low (~30 lines)  
**Risk:** Low

**Strategy A — Relax tolerance:**  Raise `_EM_CONVERGENCE_DELTA` from 1e-6 to
1e-4.  For count estimation, the difference between theta at delta=1e-4 vs 1e-6
is negligible.  Would converge at ~200 iterations for chr17.

**Strategy B — SQUAREM acceleration:**  Replace vanilla EM with SQUAREM
(Varadhan & Roland 2008), a 3-step accelerator that achieves superlinear
convergence with minimal code changes.  Salmon uses a similar acceleration.
Expected: 50–100 iterations to converge to 1e-6.

**Strategy C — Combined:**  Relax tolerance to 1e-5 + SQUAREM = ~30–50 iterations.

### Priority 3: Vectorize resolve_fragment (Expected: 3–5× scan speedup)

**Impact:** Reduce scan_and_buffer from 6.4s → ~1.5–2.5s on chr17  
**Effort:** High (~300 lines, careful correctness)  
**Risk:** Medium

**Strategy A — Batch overlap computation:**  Instead of calling
`query_exon_with_coords` per fragment-exon-block (Python loop), collect all
fragment blocks first, batch-query cgranges, and compute overlaps in numpy:

1. Collect all exon blocks: `(chrom, start, end, frag_id, block_idx)`.
2. Batch cgranges query → array of `(frag_id, t_idx, interval_type, overlap_bp)`.
3. `groupby(frag_id, t_idx)` to aggregate exon/intron overlap via numpy.
4. Compute `read_length`, `exon_bp`, `intron_bp`, `unambig_intron_bp` vectorized.

**Strategy B — Cython extension for compute_overlap_profile:**  The inner loop
is simple arithmetic (clamp, subtract, accumulate).  A Cython version accessing
cgranges data directly could achieve 10–20× speedup.

**Strategy C — Pre-compute interval tree in array form:**  Replace cgranges
Python callback with a sorted-array binary-search approach that returns raw
numpy arrays.

### Priority 4: Pre-allocate CSR Arrays (Expected: 1.5s savings)

**Impact:** Reduce `_scan_and_build_em_data` by ~1.5s  
**Effort:** Low (~50 lines)  
**Risk:** Low

**Strategy:**  Replace `list.append()` CSR building with pre-allocated numpy
arrays.  Estimate total candidates from buffer metadata (n_buffered ×
avg_candidates), allocate, and fill with an index pointer:

```python
# Pre-allocate with generous estimate
max_candidates = buffer.total_fragments * avg_cand_estimate
t_indices_arr = np.empty(max_candidates, dtype=np.int32)
log_liks_arr = np.empty(max_candidates, dtype=np.float64)
pos = 0
# In inner loop:
t_indices_arr[pos] = t_idx
log_liks_arr[pos] = log_lik
pos += 1
# After:
t_indices_arr = t_indices_arr[:pos]
```

Eliminates 22.8M `list.append` calls and the final `np.array()` conversion.

### Priority 5: Buffer Append Optimization (Expected: 1s savings)

**Impact:** Reduce buffer.append from 1.9s → ~0.5s  
**Effort:** Medium (~100 lines)  
**Risk:** Low

**Strategy:**  Batch-append resolved fragments in groups of ~100–1000 instead
of one-at-a-time.  Accumulate a list of `ResolvedFragment` objects, then
copy all fields into columnar arrays in one bulk operation per field.

### Priority 6: Cython/C Extension for Hot Loops (Expected: 5–10× overall)

**Impact:** Approach salmon/kallisto speed  
**Effort:** Very high  
**Risk:** Medium (build complexity, platform portability)

**Long-term strategy:**  Port the following to Cython or C:
- `compute_overlap_profile` inner loop
- `_add_transcript_candidates` / `_add_nrna_candidates` scoring
- Fragment.from_reads + parse_read (pysam interaction)
- EM inner loop (if equiv-class approach isn't sufficient)

---

## Expected Impact Summary

| Optimization | Est. chr17 Savings | Cumulative Time |
|---|---|---|
| Baseline | — | 22.6s |
| P1: Equiv-class EM | −13.5s | **9.1s** |
| P2: EM convergence accel | −0.3s | **8.8s** |
| P3: Vectorize resolve | −3.5s | **5.3s** |
| P4: Pre-alloc CSR | −1.5s | **3.8s** |
| P5: Buffer batch append | −1.0s | **2.8s** |

With Priorities 1–5 (all pure Python/numpy), chr17 would drop from 22.6s →
~2.8s.  The **mean** across 10 regions would drop from 7.2s → ~1.5–2.0s,
putting hulkrna within **3–4×** of salmon/kallisto.

Priority 6 (Cython) could close the remaining gap to **1–2×** of salmon/kallisto.

---

## Recommended Execution Order

1. **P1 + P2** (EM optimizations) — biggest impact, medium effort, low risk.
   Implement equiv-class EM and relax convergence tolerance together.
   **Expected outcome: 22.6s → ~9s.**

2. **P4** (Pre-alloc CSR) — quick win, low effort.
   **Expected outcome: 9s → ~7.5s.**

3. **P3** (Vectorize resolve) — biggest remaining bottleneck.
   **Expected outcome: 7.5s → ~4s.**

4. **P5** (Buffer batch) — diminishing returns, do if time permits.
   **Expected outcome: 4s → ~3s.**

5. **P6** (Cython) — only if ≤2× of salmon/kallisto is required.

All optimizations preserve mathematical equivalence and should not affect
accuracy.  Each should be validated by confirming the existing 723 tests pass
and benchmark MAE values remain unchanged.
