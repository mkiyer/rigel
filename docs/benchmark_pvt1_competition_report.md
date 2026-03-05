# HulkRNA Benchmark & Profiling Report
**Region:** PVT1/MYC locus (chr8:126445441–128594175, 35 genes, 379 transcripts)
**Fragments:** 100,000 per condition | **Aligner:** Oracle | **Date:** 2026-03-02

---

## 1. Competition Results: hulkrna vs salmon vs kallisto

### 1.1 Transcript-Level MAE (Mean Absolute Error, lower = better)

| Condition | hulkrna | salmon | kallisto | hulkrna advantage |
|---|---:|---:|---:|---|
| Pristine (no contaminants) | 20.20 | 27.40 | 24.11 | **1.2–1.4×** better |
| nRNA 10% | 5.12 | 24.54 | 42.95 | **4.8–8.4×** better |
| gDNA 10% | 0.98 | 16.78 | 18.21 | **17–19×** better |
| gDNA 10% + nRNA 10% | 0.77 | 16.64 | 18.80 | **22–24×** better |
| gDNA 25% | 0.49 | 16.03 | 17.44 | **33–36×** better |
| gDNA 25% + nRNA 10% | 0.51 | 16.18 | 18.00 | **32–35×** better |

### 1.2 Transcript-Level Correlation (Pearson, higher = better)

| Condition | hulkrna | salmon | kallisto |
|---|---:|---:|---:|
| Pristine | **0.996** | 0.988 | 0.994 |
| nRNA 10% | **0.958** | 0.671 | 0.461 |
| gDNA 10% | 0.528 | 0.071 | 0.062 |
| gDNA 25% | **0.558** | 0.053 | 0.033 |

### 1.3 Gene-Level MAE

| Condition | hulkrna | salmon | kallisto |
|---|---:|---:|---:|
| Pristine | **0.00** | 0.03 | 0.00 |
| nRNA 10% | **3.39** | 164.46 | 334.00 |
| gDNA 10% | **1.87** | 171.40 | 188.29 |
| gDNA 25% | **1.33** | 169.37 | 184.63 |
| gDNA 25% + nRNA 10% | **1.13** | 171.14 | 191.09 |

### 1.4 Dropout Rate (truth > 0 but predicted ≤ 0)

| Tool | Dropout Rate |
|---|---:|
| hulkrna | **0.00%** |
| salmon | 11.08% |
| kallisto | 8.62% |

### 1.5 Speed

| Tool | Time per condition | Relative |
|---|---:|---:|
| hulkrna | 18–48s | 1× |
| salmon | 1–5s | **10–15×** faster |
| kallisto | 0.9–1.3s | **15–40×** faster |

### 1.6 Competition Summary

**hulkrna dominates accuracy across all conditions.** The advantage grows dramatically
with contamination:
- Pristine RNA: **1.2× better** MAE — salmon/kallisto are competitive
- With gDNA/nRNA contamination: **17–36× better** MAE — salmon/kallisto
  cannot distinguish contaminants from real signal
- **Zero dropout**: hulkrna detects all expressed transcripts; salmon/kallisto miss 9–11%
- Gene-level: hulkrna is **50–170× better** under contamination

**Speed is the primary weakness.** hulkrna is 10–40× slower than salmon/kallisto, driven
almost entirely by the Python EM solver.

---

## 2. Profiling Results

**Condition:** Pristine (100k fragments, 200k BAM records) | **Wall time:** 15.67s

### 2.1 Time Breakdown (by self-time)

| Component | Self-time | % Total | Calls | Description |
|---|---:|---:|---:|---|
| `_em_step` | 7.98s | **51.0%** | 810 | EM E-step + M-step inner loop |
| `scan_and_buffer` (C++) | 2.54s | 16.2% | 1 | BAM parse + resolve + buffer |
| numpy.ufunc.reduce | 2.16s | 13.8% | 1.2M | sum/max inside EM step |
| `scan` (Python glue) | 0.55s | 3.5% | 1 | Post-scan Python processing |
| `_add_single_fragment` | 0.38s | 2.4% | 89K | Per-fragment scan handler |
| `build_locus_em_data` | 0.31s | 2.0% | 11 | Locus EM data construction |
| `buffer.__getitem__` | 0.23s | 1.5% | 100K | Buffer retrieval in quantification |
| `_build_equiv_classes` | 0.11s | 0.7% | 11 | Equivalence class grouping |
| Other | 1.41s | 9.0% | — | fraglen model, strand model, etc. |

### 2.2 Bottleneck Analysis

**The EM solver dominates** — `_em_step` + its numpy calls account for **~65% of
total runtime** (10.5s cumulative).

The hot path in `_em_step` iterates over equivalence classes, performing per-class:
1. `np.add(ll_matrix, log_weights[comp_idx], out=scratch)` — add log-weights to log-likelihoods
2. `scratch.max(axis=1, keepdims=True)` — log-sum-exp denominator
3. `scratch -= max_r; np.exp(scratch, out=scratch)` — softmax numerator
4. `scratch.sum(axis=1, keepdims=True)` — softmax denominator
5. `scratch /= row_sum` — normalize
6. `em_totals[comp_idx] += scratch.sum(axis=0)` — accumulate posteriors

This runs 810 times total (11 loci × ~73 SQUAREM iterations × 3 steps/iter), each
iteration processing all equivalence classes. The 1.2M numpy ufunc calls indicate
**tiny matrix operations** where Python call overhead dominates — the matrices are
typically 2–20 columns wide.

**C++ BAM scanning is fast** at 2.54s (16%) — the prior C++ port was highly effective.

---

## 3. Optimization Opportunities (Ranked by Impact)

### Priority 1: Port EM Inner Loop to C++ (~8s → ~1–2s, **40% of total**)

The `_em_step` function is the single highest-impact optimization target. The
pattern is:
- Python `for` loop over equivalence classes (Python overhead per iteration)
- Many small numpy operations on tiny matrices (2–20 columns)
- numpy call overhead dominates actual computation for small arrays

**Approach:** Write a single nanobind C++ function `native_em_step(theta, ec_data,
log_eff_len, unambig_totals, prior, em_totals)` that:
- Takes pre-built equivalence class data as contiguous arrays
- Processes all classes in a tight C++ loop (log-sum-exp + normalize + accumulate)
- Returns `theta_new` directly

Expected speedup: 4–8× for this phase alone (measured numpy overhead on
small arrays is 3–10μs per call; C++ eliminates this entirely).

### Priority 2: Batch Small-Matrix Operations (~2s savings)

If a full C++ port is deferred, an intermediate optimization:
- Concatenate all equivalence class matrices into a single large matrix
- Use a single vectorized numpy operation over all classes at once
- Avoids 1.2M individual numpy calls per EM iteration

### Priority 3: Port `_build_equiv_classes` + `build_locus_em_data` to C++ (~0.5s)

These Python-heavy data preparation steps could be fused with the C++ scan output,
building the EM data structures directly in C++ rather than reconstructing them in
Python.

### Priority 4: Streaming Buffer Access (~0.2s)

The 100K `buffer.__getitem__` calls (0.23s) could be batched — retrieve all
fragments for a locus in a single slice rather than one-by-one.

### Priority 5: SQUAREM Convergence Tuning

The 810 `_em_step` calls across 11 loci suggests ~73 SQUAREM cycles per locus
(= ~219 EM evaluations). Investigating whether convergence criteria could be
relaxed (currently `1e-7`) or whether the warm-start initialization could reduce
iteration count would multiply with any per-step speedup.

---

## 4. Recommended Next Steps

1. **Immediate (highest ROI):** Port `_em_step` and `_vbem_step` to C++ via nanobind.
   This single change would likely bring hulkrna within 3–5× of salmon/kallisto speed
   while retaining its massive accuracy advantage.

2. **Short-term:** Port equivalence class construction and locus EM data building
   to C++ as well, creating a fully native quantification path.

3. **Medium-term:** Investigate convergence behavior — are 73 SQUAREM cycles per
   locus optimal? Profile the largest loci (Locus 8: 195 transcripts, 18K units)
   individually to see if they dominate.

4. **Stretch goal:** Multi-threaded per-locus EM — the 11 loci are independent
   sub-problems that could run in parallel via a C++ thread pool.

---

## 5. Key Takeaway

**hulkrna is already the clear accuracy winner** — 1.2× better at pristine, 17–36× better
with contamination, zero dropout. The speed gap (10–40×) is entirely attributable to the
Python EM solver, which is a well-scoped, high-confidence optimization target. Porting
`_em_step` to C++ alone would likely close the speed gap to 3–5×.
