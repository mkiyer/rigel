# hulkrna Performance Profiling Report & Optimization Plan

## Profiling Setup

- **Regions tested:** EGFR (14 tx, 2 genes), chr19_dense (64 tx, 13 genes), chr1_dense (55 tx, 12 genes)
- **Fragments per region:** 100,000
- **Simulation:** gDNA + nRNA injected, strand_specificity=0.95

## Throughput Summary

| Region | Tx | Genes | scan_and_buffer | count_from_buffer | Total (s) | Frags/sec |
|--------|---:|------:|----------------:|------------------:|----------:|----------:|
| EGFR | 14 | 2 | 3.88 | 8.93 | 12.81 | 7,807 |
| chr19_dense | 64 | 13 | 4.65 | 9.85 | 14.49 | 6,899 |
| chr1_dense | 55 | 12 | 4.69 | 8.71 | 13.40 | 7,462 |

**Current throughput: ~7,400 frags/sec** (100k fragments → ~13s).
For a typical 30M-fragment RNA-seq BAM this translates to ~68 minutes.
Salmon processes the same in ~1–2 minutes.

## Bottleneck Breakdown (EGFR, 100k fragments, 27.4s w/ cProfile overhead)

### Stage 1: `count_from_buffer` — 18.5s (67.6%)

| Function | Calls | Cum. Time | % Pipeline |
|----------|------:|----------:|-----------:|
| `_scan_and_build_em_data` | 1 | 11.35s | 41.5% |
| → `_add_transcript_candidates` | 90,425 | 5.61s | 20.5% |
| → → `_score_candidate` | 688,419 | 4.68s | 17.1% |
| → `_add_nrna_candidates` | 81,421 | 4.04s | 14.8% |
| → → `_score_nrna_candidate` | 637,613 | 3.07s | 11.2% |
| `run_locus_em` | 1 | 2.97s | 10.9% |
| `_build_locus_em_data` | 1 | 1.71s | 6.3% |
| `assign_locus_ambiguous` | 1 | 1.43s | 5.2% |
| `_build_loci` (union-find) | 1 | 0.98s | 3.6% |

### Stage 2: `scan_and_buffer` — 8.95s (32.7%)

| Function | Calls | Cum. Time | % Pipeline |
|----------|------:|----------:|-----------:|
| `resolve_fragment` | 100,000 | 4.51s | 16.5% |
| → `compute_overlap_profile` | 100,000 | 1.74s | 6.4% |
| → `query_exon` / `query_exon_with_coords` | 142,050 + 142,050 | 1.23s | 4.5% |
| `parse_bam_file` | 100,001 | 1.30s | 4.8% |
| `fragment.from_reads` | 100,000 | 0.97s | 3.5% |
| `buffer.append` | 100,000 | 0.92s | 3.4% |

### Callee Hotspots (Leaf Functions)

| Function | Calls | Own Time | Notes |
|----------|------:|----------:|-------|
| `insert_model.log_likelihood` | 1,403,914 | 0.94s | Called per (fragment × candidate); 2 np.log per call |
| `insert_model.total_weight` | 1,403,923 | 0.28s | `counts.sum()` per call (1001-element array) |
| `strand_model.strand_likelihood` | 1,517,025 | 0.35s | 3 property lookups per call (α, β, divisions) |
| `strand_model.p_r1_sense/antisense` | 1.5M+1.4M | 1.96s | Property chain: `alpha → n_same + prior`, `beta → …` |
| `strand_model.alpha/beta` | 3M+1.5M | 0.67s+0.35s | Python @property overhead at this call volume |
| `list.append` | 15.8M | 1.13s | CSR accumulation via Python lists |
| `enum.__call__` | 2.45M | 0.78s | `Strand(int)`, `SpliceType(int)` conversions |
| `builtins.max`/`min` | 5.7M+2.7M | 0.99s | Clamping in scoring functions |

## Root Cause Analysis

### Problem 1: Per-candidate scoring in pure Python (49% of pipeline)

`_score_candidate` is called **688k times** (7.6 per fragment) and calls:
- `strand_model.strand_likelihood` → `p_r1_sense` → `alpha` property → `n_same + prior_alpha`
- `insert_model.log_likelihood` → `total_weight` → `counts.sum()` (sums a 1001-element numpy array!)
- `np.log(max(p_strand, floor))` + `np.log(count + 1)` − `np.log(total + size)` → 3 scalar np.log calls
- `overhang_bp * overhang_log_penalty` → scalar multiply

`_score_nrna_candidate` has the same pattern, called 637k times.

These are **scalar Python loops calling numpy scalar operations** — the worst-case performance pattern. Each call does:
- 1 numpy `.sum()` on a 1001-element array (via `total_weight` property, called every time!)
- 2–3 scalar `np.log()` calls
- Multiple Python property lookups with arithmetic

### Problem 2: `total_weight` recomputed per call (28% of insert_model time)

`insert_model.total_weight` calls `self.counts.sum()` on a 1001-element array
**1.4 million times**. The histogram doesn't change during scoring — this is
computed once during model training but the property recomputes it every call.

### Problem 3: Strand model property chain overhead (7% of pipeline)

`strand_likelihood` → `p_r1_sense` → `self.alpha / (self.alpha + self.beta)`.
Each of `alpha`, `beta` is a @property computing `self.prior_alpha + self.n_same`.
With 1.5M+ calls, the Python property dispatch overhead alone is ~2s.

### Problem 4: Enum construction overhead (2.8%)

`Strand(int(…))` and `SpliceType(int(…))` are called millions of times
in scoring. Python enum construction is slow (~300ns/call).

### Problem 5: `compute_overlap_profile` per-fragment (6.4%)

Loops over fragment exon blocks, queries cgranges, accumulates per-transcript
overlap using Python dicts. Called 100k times with multiple index queries each.

### Problem 6: CSR accumulation via Python lists (4.1%)

15.8M `list.append` calls build the CSR candidate matrix. This creates
enormous Python list objects before conversion to numpy arrays.

---

## Optimization Plan

### Tier 1: Low-Hanging Fruit (Expected: 3–5× speedup)

#### 1.1 Cache `total_weight` in InsertSizeModel
**Impact:** Eliminates 1.4M × `arr.sum()` → instant lookup
**Effort:** Trivial (add `self._total_weight` updated in `add()`)

```python
# insert_model.py
def add(self, insert_size, weight=1.0):
    idx = min(max(insert_size, 0), self.max_size)
    self.counts[idx] += weight
    self.n_observations += 1
    self._total_weight += weight  # ← cache

@property
def total_weight(self):
    return self._total_weight
```

#### 1.2 Cache strand model probabilities
**Impact:** Eliminates 3M+ property lookups → instant float read
**Effort:** Trivial (compute once after training, store as plain float)

```python
# strand_model.py — after training is complete:
model._cached_p_sense = model.alpha / (model.alpha + model.beta)
model._cached_p_antisense = 1.0 - model._cached_p_sense
```

#### 1.3 Pre-compute `log_likelihood` lookup table
**Impact:** Replace 1.4M × (2 np.log + array lookup) with 1.4M × array[idx]
**Effort:** Small (pre-compute log-prob array after training)

```python
# insert_model.py — after training:
total = self._total_weight
self._log_prob = np.log(self.counts + 1.0) - np.log(total + self.max_size + 1)

def log_likelihood(self, insert_size):
    return self._log_prob[min(max(insert_size, 0), self.max_size)]
```

#### 1.4 Replace enum construction with int comparisons
**Impact:** Eliminates 2.5M enum.__call__s
**Effort:** Small (use int constants instead of Strand(int_val))

### Tier 2: Vectorized Scoring (Expected: 5–10× additional speedup)

#### 2.1 Batch `_score_candidate` for all candidates of a fragment
**Impact:** Replace per-candidate Python loop with vectorized numpy
**Effort:** Medium

Instead of calling `_score_candidate` per candidate, collect all candidates'
strand/insert/overhang data into arrays and score with one vectorized call:

```python
def _score_candidates_batch(t_strands, splice_type, insert_size,
                            exon_strand, overhang_bps,
                            p_sense, p_antisense, log_prob_table, ovh_lp):
    # All operations vectorized over N candidates
    same_strand = (exon_strand == t_strands)  # bool array
    p_strand = np.where(same_strand, p_sense, p_antisense)
    log_strand = np.log(np.maximum(p_strand, 1e-30))
    log_insert = log_prob_table[min(max(insert_size, 0), max_size)]  # scalar, broadcast
    log_lik = log_strand + log_insert + overhang_bps * ovh_lp
    return log_lik
```

#### 2.2 Batch `resolve_fragment` overlap computation
**Impact:** Replace per-exon Python dict accumulation with array ops
**Effort:** Medium-High

#### 2.3 Vectorize `_build_locus_em_data`

Currently loops over all units building per-locus CSR data.
This can be done with numpy fancy indexing and `np.unique` + `np.add.at`.

### Tier 3: Structural Optimizations (Expected: 2–3× additional)

#### 3.1 Replace Python-list CSR accumulation with pre-allocated arrays
**Impact:** Eliminate 15.8M list.append calls
**Effort:** Medium (estimate max candidates per fragment, pre-allocate)

#### 3.2 Store buffer data as columnar numpy arrays
The buffer already uses columnar storage, but the per-fragment iteration
re-creates Python objects. A fully vectorized path could operate directly
on the columnar arrays without per-fragment object creation.

#### 3.3 Cython/C extension for the inner scoring loop
**Impact:** 10–50× for the scoring inner loop
**Effort:** High
If Tiers 1+2 are insufficient, the inner `_scan_and_build_em_data` loop
could be compiled with Cython. The logic is simple arithmetic — ideal
for Cython.

### Tier 4: Algorithmic (Expected: variable)

#### 4.1 Reduce candidate multiplicity
Many fragments have 7+ transcript candidates in overlapping loci.
Pre-filtering candidates more aggressively (e.g., stricter overlap
thresholds, early exon-only filtering) reduces the inner loop count.

#### 4.2 EM warm-starting
`run_locus_em` takes 3s for EGFR's single large locus (90k units ×
14 transcripts). Warm-starting from the high-confidence unique counts
could reduce iteration count.

---

## Recommended Execution Order

1. **Tier 1.1 + 1.2 + 1.3 + 1.4** — Cache total_weight, strand probs, log-likelihood LUT, remove enums. ~1 hour of work, expected **2–3× speedup** on the scoring path.

2. **Tier 2.1** — Vectorize `_score_candidate` batch. ~2–4 hours, expected **3–5× additional** on `_scan_and_build_em_data` (the 41.5% bottleneck).

3. **Tier 3.1** — Pre-allocated arrays for CSR. ~1–2 hours, eliminates list.append overhead.

4. **Re-profile** after Tiers 1+2+3 to see if EM or `resolve_fragment` becomes the new dominant cost.

5. **Tier 2.2** — Vectorize overlap computation if scan_and_buffer becomes dominant.

**Target:** 10× overall improvement → ~70k frags/sec → ~7 min for 30M-fragment BAM.

---

## Profiling Artifacts

- Full cProfile data: `profiling/<region>/profile_stats.prof`
- Text summaries: `profiling/<region>/profile_stats.txt`
- Caller info: `profiling/<region>/profile_callers.txt`
- Stage breakdowns: `profiling/<region>/pipeline_breakdown.md`
- Combined JSON: `profiling/profile_summary.json`
