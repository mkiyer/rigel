# Performance Improvement Plan — Post C++ EM Porting

**Date:** 2026-03-03
**Baseline:** PVT1/MYC locus (200k fragments, 379 transcripts, 35 genes)
**Environment:** macOS ARM64, Apple M-series, Python 3.13, C++17/nanobind

---

## 1 Current Profile Summary

### 1.1 Wall Clock Times

| Condition | Total | EM Solver | Scan+Buffer | Score/Build | Other |
|---|---:|---:|---:|---:|---:|
| **Pristine** | 14.84s | 5.80s (39%) | 5.23s (35%) | 2.76s (19%) | 1.05s (7%) |
| **nRNA low** | 67.62s | 54.92s (81%) | 4.60s (7%) | 6.05s (9%) | 2.05s (3%) |
| **gDNA+nRNA low** | 21.24s | 15.59s (73%) | 2.10s (10%) | 2.85s (13%) | 0.70s (3%) |

### 1.2 Key Observations

1. **C++ EM solver is still the dominant bottleneck** in all conditions.
   - Pristine: 10 locus calls, 0.58s avg → 5.80s total
   - nRNA: 9 locus calls, 6.10s avg → 54.92s total
   - The nRNA condition creates **1–2 giant merged loci** (intronic reads
     link many genes into one connected component with hundreds of
     transcripts × 3 components = 700+ EM components)

2. **scan_and_buffer** (C++ BAM reader) takes 4.6–5.2s consistently.
   Good as-is but could be optimized further.

3. **scan.scan()** (Python fragment scoring loop) takes 2.8–6.0s:
   - `_add_single_fragment`: 0.7–1.0s (170k–197k calls)
   - `buffer.__getitem__`: 0.4–0.5s (200k calls, creates Python objects)
   - `_finalize_unit_metadata`: 0.1–0.4s
   - `_gdna_log_lik`: 0.06–0.13s (90k–189k calls)

4. **build_locus_em_data** (NumPy vectorized): 0.3–1.3s

5. **Enum overhead**: 520k `enum.Enum.__call__` invocations cost 0.14s
   (integer strand/splice type lookups going through Python enum machinery).

---

## 2 Optimization Roadmap

### Phase A: C++ EM Solver Optimizations (Target: 3–10× speedup on EM)

These require no Python changes — pure C++ improvements inside
`em_solver.cpp`.

#### A1. SIMD-Vectorized Inner Loop (Priority: HIGH, Impact: 2–4×)

The `em_step_kernel` inner loop iterates over `n × k` elements doing
`exp`, multiply, divide operations that are trivially SIMD-parallelizable.

**Current code** (scalar):
```cpp
for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
        scratch[offset + j] = ll[offset + j] + log_weights[cidx[j]];
    }
    // log-sum-exp, normalize, accumulate...
}
```

**Optimization:**
- Use NEON intrinsics (ARM) or SSE/AVX (x86) for the inner `k`-loop
- For common small `k` (1–5), specialize with unrolled templates
- For large `k` (>8), use `vfmaq_f64` / `_mm256_fmadd_pd` for
  add + exp + reduce pipeline
- Column accumulation can use horizontal SIMD reduction

**Estimated impact:** 2–4× on the EM kernel, or 20–60% total for nRNA
condition.

#### A2. Sorted Equivalence Classes + Early Termination (Priority: HIGH, Impact: 1.5–3×)

Sort equivalence classes by size (largest first). For the nRNA mega-locus,
a few classes may contain >80% of units. Iterate large classes first and
check convergence between classes within a SQUAREM step.

Also: track per-iteration `delta` and skip the stabilisation SQUAREM step
when convergence is already fast (falling back to plain EM). For
well-conditioned loci, SQUAREM's 3× cost per iteration is wasted.

#### A3. Sparse Component Representation (Priority: MEDIUM, Impact: 2–5× on large loci)

For a mega-locus with 700+ components but where each equivalence class
touches only 2–5 components, the `log_weights[cidx[j]]` gather is already
sparse. But the M-step still sums over all `n_components` to normalize.

**Optimization:**
- Track active component set per iteration (components receiving any
  posterior mass)
- Skip zero-weight components in M-step normalization
- Use a bitset to track active components (fits in a few cache lines
  for 1024 components)

#### A4. Cache-Friendly Memory Layout (Priority: MEDIUM, Impact: 1.3–2×)

Current layout stores `ll_flat`, `wt_flat`, `scratch` as separate
contiguous buffers. For each equivalence class access pattern (`ll` read +
`scratch` write per row), interleaving them would improve locality.

**Optimization:**
- Structure-of-arrays → array-of-structures for the row-major data
- Align rows to cache line boundaries (64B)
- Pre-sort equivalence classes by `k` to improve branch prediction

#### A5. OpenMP Parallelism for Large Loci (Priority: LOW, Impact: 2–8×)

The E-step (iterate over equivalence classes) is embarrassingly parallel —
each class's posterior computation is independent. The M-step accumulation
needs a reduction.

```cpp
#pragma omp parallel for reduction(+:em_totals[:n_components]) schedule(dynamic)
for (size_t ec_idx = 0; ec_idx < ec_data.size(); ++ec_idx) {
    em_step_kernel(ec_data[ec_idx], log_weights, em_totals);
}
```

This would mainly help the nRNA mega-locus (54.9s → potentially 7–15s
on 8 cores). The overhead for small loci would be negligible if
conditional on locus size (>1000 units).

### Phase B: Fused Scan Pipeline (Target: 2–3× on scan stage)

The scan stage currently has a Python-level per-fragment loop that calls
into C++ for scoring but does Python bookkeeping between calls.

#### B1. Fused C++ Scan Loop (Priority: HIGH, Impact: 2–3×)

Move the entire `FragmentRouter.scan()` loop into C++:
- Currently: Python iterates chunks → Python `__getitem__` per fragment →
  Python calls `score_emit_fragment` → Python appends to `array.array`
- Target: Single C++ function takes the entire buffer and produces the
  final ScoredFragments arrays in one pass

This eliminates:
- 200k `BufferedFragment.__getitem__` calls (0.5s)
- 200k Python-level append operations (0.15s)
- 170k `_finalize_unit_metadata` calls (0.08–0.4s)
- 90k `_gdna_log_lik` calls (0.06–0.13s)
- Python loop overhead

**Estimated impact:** Save 1.0–2.0s (20–40% of scan time).

#### B2. Eliminate Enum Overhead (Priority: MEDIUM, Impact: 0.15s)

Replace `enum.Enum.__call__` for strand/splice type lookups with plain
integer comparisons. The enum machinery costs 520k calls × 0.27µs = 0.14s.

Either use `IntEnum` (faster) or raw integer constants throughout the hot
path.

### Phase C: build_locus_em_data Optimization (Target: 2×)

#### C1. Move to C++ (Priority: MEDIUM, Impact: 0.3–1.3s saved)

`build_locus_em_data` does NumPy vectorized operations (fancy indexing,
repeat, lexsort, dedup) that are fast but have per-call overhead of
60–130ms. For the nRNA condition with only 9 loci, this is 1.3s total.

A C++ implementation would:
- Avoid Python GIL overhead per locus
- Use cache-friendly single-pass CSR extraction
- Combine with EM solver: `build_locus_em_data` + `run_locus_em` in one
  C++ call, eliminating intermediate numpy array creation

### Phase D: Algorithmic Improvements (Target: reduce nRNA mega-locus)

#### D1. Locus Decomposition for nRNA (Priority: HIGH, Impact: 10–50×)

The nRNA condition creates giant loci because intronic reads link genes
across overlapping genomic intervals. A hierarchical decomposition:

1. Build initial loci using only mRNA (spliced) evidence
2. Assign nRNA candidates to their nearest mRNA locus
3. Only merge loci if they share strong ambiguous mRNA evidence

This would break the 700-component mega-locus into many smaller loci
(~20–50 components each), where EM converges in milliseconds instead of
seconds.

#### D2. nRNA Pre-assignment (Priority: MEDIUM, Impact: reduce unit count)

Fragments with >80% intronic overlap and matching strand pattern could be
pre-assigned to nRNA with high confidence, removing them from the EM
problem. This shrinks the unit count for the mega-locus.

#### D3. Convergence-Adaptive Iteration Budget (Priority: LOW, Impact: 1.5×)

Currently `max_iterations` is fixed at 1000 (333 SQUAREM iterations).
Track convergence rate and dynamically reduce budget for fast-converging
loci. Most small loci converge in <10 SQUAREM iterations.

---

## 3 Estimated Impact Summary

| Phase | Change | Time Saved (nRNA) | Time Saved (Pristine) |
|---|---|---:|---:|
| A1 | SIMD inner loop | 15–25s | 1.5–3.0s |
| A2 | Sorted ECs + early term | 5–15s | 0.5–1.5s |
| A3 | Sparse component repr | 5–10s | 0.2–0.5s |
| A5 | OpenMP parallelism | 20–40s | 1.0–3.0s |
| B1 | Fused C++ scan loop | 1.5–2.5s | 1.0–2.0s |
| C1 | C++ build_locus_em_data | 0.8–1.3s | 0.3–0.6s |
| D1 | Locus decomposition | 40–50s | — |

### Projected Performance After All Optimizations

| Condition | Current | Projected | Speedup |
|---|---:|---:|---:|
| Pristine | 14.84s | 3–6s | 2.5–5× |
| nRNA low | 67.62s | 4–10s | 7–17× |
| gDNA+nRNA low | 21.24s | 3–6s | 3.5–7× |

### Comparison with Competitors

| Tool | Current (Pristine) | Projected (Pristine) | Current (nRNA) | Projected (nRNA) |
|---|---:|---:|---:|---:|
| hulkrna | 14.84s | 3–6s | 67.62s | 4–10s |
| salmon | 1.4s | — | 1.1s | — |
| kallisto | 1.1s | — | 1.3s | — |

---

## 4 Recommended Implementation Order

1. **D1 — Locus decomposition** (highest ROI: fixes the 67s nRNA case)
2. **A5 — OpenMP parallelism** (easy win for large loci, 2–8×)
3. **B1 — Fused C++ scan loop** (eliminates Python per-fragment overhead)
4. **A1 — SIMD inner loop** (purely mechanical, portable with intrinsics)
5. **A2 — Sorted ECs + early termination** (algorithmic, modest code change)
6. **C1 — C++ build_locus_em_data** (combine with EM for less overhead)
7. **A3 — Sparse components** (only matters for very large loci)
8. **B2 — Enum elimination** (quick, small impact)
9. **D2 — nRNA pre-assignment** (improves both accuracy and speed)
10. **A4 — Cache layout** (micro-optimization, measure before committing)
