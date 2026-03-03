# Plan: Port EM Solver to C++ via nanobind

## 1. Motivation

Profiling on the PVT1/MYC benchmark (100k fragments, pristine reads) shows:

| Function              | Time (s) | % Total |
|-----------------------|---------:|--------:|
| `_em_step`            |     7.98 |     51% |
| `numpy.ufunc.reduce`  |     2.16 |     14% |
| `scan_and_buffer`     |     2.54 |     16% |
| Other                 |     2.99 |     19% |
| **Total**             | **15.67** | **100%** |

`_em_step` is a Python `for`-loop over equivalence classes, each calling
5-6 tiny numpy operations on 2–20 column matrices. With ~810 calls per
run (11 loci × ~73 SQUAREM iters × 3 evals), Python call overhead and
numpy dispatch dominate actual computation. The 2.16 s `ufunc.reduce` is
*inside* `_em_step`.

A single C++ function replacing the entire per-locus SQUAREM loop
(eq-class building + EM iterations + convergence check) eliminates:

- ~810 Python→numpy roundtrips per run
- ~1.2M numpy ufunc dispatches
- Python `dict` hashing in `_build_equiv_classes`
- Per-iteration `ndarray` allocation

Target: **5–15× speedup** on the EM solver, bringing total pipeline time
from ~16 s to ~4–6 s for the benchmark workload (putting hulkrna in the
same performance ballpark as salmon/kallisto while retaining its accuracy
advantage).

---

## 2. Scope — What Moves to C++

### In scope (new `_em_impl` module)

| Python function | Why |
|---|---|
| `_em_step()` | THE bottleneck — 51% of runtime |
| `_vbem_step()` | Identical inner loop, digamma variant |
| `_build_equiv_classes()` | Builds the data fed to EM; keeping it in C++ avoids materializing Python lists of tuples |
| `_apply_bias_correction_uniform()` | Pre-EM log-lik adjustment; fusing it avoids a separate numpy pass over the same arrays |
| SQUAREM acceleration loop | Currently in `run_locus_em()`; porting it avoids 810 Python↔C++ calls per run |
| Coverage-weighted warm start + OVR prior | Part of `run_locus_em()`'s pre-loop setup; pure arithmetic on the eq-class data |
| Post-prune re-EM loop | 10 plain EM iterations after pruning; trivial once the above is in C++ |

### Out of scope (stays in Python)

| Component | Why |
|---|---|
| `AbundanceEstimator` class (state management) | Accumulator arrays, config, output DataFrames — not hot |
| `assign_locus_ambiguous()` | Post-EM posterior scatter; runs once per locus, vectorized numpy, ~0% of runtime |
| `ScanData` / `LocusEMInput` dataclasses | Data containers; only constructed once per locus |
| `build_locus_em_data()` (locus.py) | Locus graph construction; measured at 0.31 s (2%) |
| `build_loci()` (locus.py) | scipy connected_components; not hot |
| `compute_nrna_init()`, `compute_eb_gdna_priors()` | Pre-EM initialization; negligible time |
| Output: `get_counts_df()`, `get_gene_counts_df()` | DataFrame construction; not hot |

---

## 3. C++ API Design

### 3.1 Entry point: `run_locus_em_native()`

A single function that takes CSR data + config and returns converged
`(theta, alpha, em_totals)`:

```cpp
// hulkrna._em_impl

std::tuple<nb::ndarray<double>, nb::ndarray<double>, nb::ndarray<double>>
run_locus_em_native(
    // --- CSR per-locus data (from LocusEMInput) ---
    nb::ndarray<const int64_t,  nb::ndim<1>, nb::c_contig> offsets,       // [n_units+1]
    nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig> t_indices,     // [n_candidates]
    nb::ndarray<const double,   nb::ndim<1>, nb::c_contig> log_liks,      // [n_candidates]
    nb::ndarray<const double,   nb::ndim<1>, nb::c_contig> coverage_wts,  // [n_candidates]
    nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig> tx_starts,     // [n_candidates]
    nb::ndarray<const int32_t,  nb::ndim<1>, nb::c_contig> tx_ends,       // [n_candidates]
    nb::ndarray<const int64_t,  nb::ndim<1>, nb::c_contig> bias_profiles, // [n_components]

    // --- Per-component vectors ---
    nb::ndarray<const double, nb::ndim<1>, nb::c_contig> unique_totals,   // [n_components]
    nb::ndarray<const double, nb::ndim<1>, nb::c_contig> effective_lens,  // [n_components]
    nb::ndarray<const double, nb::ndim<1>, nb::c_contig> prior_eligible,  // [n_components] bool-ish

    // --- Scalar config ---
    int n_components,
    double prior_alpha,
    double prior_gamma,
    int max_iterations,
    double convergence_delta,
    bool use_vbem,                  // true = VBEM, false = MAP-EM
    double prune_threshold,         // < 0 means disabled
    int post_prune_iters
);
// Returns: (theta[n_components], alpha[n_components], em_totals[n_components])
```

### 3.2 Python call site (in `run_locus_em`)

```python
from hulkrna._em_impl import run_locus_em_native

theta, alpha, em_totals = run_locus_em_native(
    locus_em.offsets,
    locus_em.t_indices,
    locus_em.log_liks,      # mutated in-place by bias correction inside C++
    locus_em.coverage_weights,
    locus_em.tx_starts,
    locus_em.tx_ends,
    locus_em.bias_profiles,
    locus_em.unique_totals,
    locus_em.effective_lengths,
    (locus_em.prior > 0).astype(np.float64),  # prior eligibility mask
    locus_em.n_components,
    self.em_config.prior_alpha,
    self.em_config.prior_gamma,
    em_iterations,
    em_convergence_delta,
    self.em_config.mode == "vbem",
    self.em_config.prune_threshold if self.em_config.prune_threshold is not None else -1.0,
    _POST_PRUNE_ITERS,
)
```

This replaces the entire body of `run_locus_em()` from the
`_apply_bias_correction()` call through the post-prune block (lines
~800-1040 of `estimator.py`), leaving only the early-exit check and the
return statement.

---

## 4. Internal C++ Architecture

### 4.1 File layout

```
src/hulkrna/native/
    em_solver.cpp          # All EM logic + nanobind module definition
    em_solver.h            # (optional) shared declarations if needed
```

Single file preferred for now — the entire solver is ~400-500 lines of
C++. Can split later if headers are needed by other modules.

### 4.2 Key internal functions

```cpp
namespace hulk {
namespace em {

// --- Bias correction (uniform fast-path) ---
void apply_bias_correction_uniform(
    double* log_liks,           // mutated in-place
    const int32_t* t_indices,
    const int32_t* tx_starts,
    const int32_t* tx_ends,
    const int64_t* profile_lengths,
    size_t n_candidates,
    int n_components
);

// --- Equivalence class builder ---
struct EquivClass {
    std::vector<int32_t> comp_idx;   // k component indices
    std::vector<double>  ll_flat;    // n*k log-likelihoods (row-major)
    std::vector<double>  wt_flat;    // n*k coverage weights (row-major)
    std::vector<double>  scratch;    // n*k workspace
    int n;  // number of units in this class
    int k;  // number of components
};
std::vector<EquivClass> build_equiv_classes(
    const int64_t* offsets,
    const int32_t* t_indices,
    const double*  log_liks,
    const double*  coverage_wts,
    int n_units
);

// --- Single EM step (MAP-EM) ---
// Writes em_totals in-place, returns theta_new as vector.
void em_step(
    const double* theta,
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    const double* prior,
    double* em_totals,    // output: zeroed then accumulated
    double* theta_new,    // output: normalized
    int n_components
);

// --- Single VBEM step ---
void vbem_step(
    const double* alpha,
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    const double* prior,
    double* em_totals,
    double* alpha_new,
    int n_components
);

// --- digamma (ψ) function ---
double digamma(double x);

} // namespace em
} // namespace hulk
```

### 4.3 Equivalence class hashing

Replace Python `dict[tuple, list[int]]` with a C++ hash map.
The key is the sorted candidate component set. Since each unit's
candidates are already sorted in the CSR, we can hash the raw
`int32_t*` slice:

```cpp
struct SliceHash {
    size_t operator()(const std::pair<const int32_t*, size_t>& s) const {
        // FNV-1a or xxhash on the raw bytes
        size_t h = 14695981039346656037ULL;
        for (size_t i = 0; i < s.second; ++i) {
            h ^= static_cast<size_t>(s.first[i]);
            h *= 1099511628211ULL;
        }
        return h;
    }
};
```

Actually, since we need stable keys (the CSR data doesn't move), we
can use `std::string_view` over the raw bytes, or an
`std::unordered_map<std::vector<int32_t>, ...>` for simplicity given
the small number of unique classes per locus (typically < 100).

### 4.4 EM inner loop — the hot kernel

The critical loop iterates over equivalence classes. For each class
with `n` units and `k` candidates:

```cpp
void em_step_kernel(
    const EquivClass& ec,
    const double* log_weights,   // log(theta) - log(eff_len), [n_components]
    double* em_totals            // accumulated posterior sums [n_components]
) {
    const int n = ec.n;
    const int k = ec.k;
    const double* ll = ec.ll_flat.data();
    double* scratch = ec.scratch.data();
    const int32_t* cidx = ec.comp_idx.data();

    // 1. Add log_weights to log-likelihoods
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            scratch[i * k + j] = ll[i * k + j] + log_weights[cidx[j]];
        }
    }

    // 2. Log-sum-exp normalize (row-wise)
    for (int i = 0; i < n; ++i) {
        double* row = scratch + i * k;
        double max_val = row[0];
        for (int j = 1; j < k; ++j)
            max_val = std::max(max_val, row[j]);

        double sum = 0.0;
        for (int j = 0; j < k; ++j) {
            row[j] = std::exp(row[j] - max_val);
            sum += row[j];
        }

        if (sum > 0.0 && std::isfinite(sum)) {
            double inv_sum = 1.0 / sum;
            for (int j = 0; j < k; ++j)
                row[j] *= inv_sum;
        } else {
            for (int j = 0; j < k; ++j)
                row[j] = 0.0;
        }
    }

    // 3. Accumulate posterior column sums → em_totals
    for (int j = 0; j < k; ++j) {
        double col_sum = 0.0;
        for (int i = 0; i < n; ++i)
            col_sum += scratch[i * k + j];
        em_totals[cidx[j]] += col_sum;
    }
}
```

This tight loop has no Python calls, no numpy dispatch, no memory
allocation. The compiler can auto-vectorize the inner `j` loop with
SSE/AVX. For large equivalence classes (n > 1000, k > 10), the
row-major layout gives good cache behavior.

### 4.5 digamma implementation

A self-contained C++ digamma using the asymptotic series (same
accuracy as scipy, ~12 significant digits):

```cpp
double digamma(double x) {
    // Shift x up so asymptotic series is accurate
    double result = 0.0;
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    // Asymptotic: ψ(x) ~ ln(x) - 1/(2x) - Σ B_{2k}/(2k·x^{2k})
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;
    result += std::log(x) - 0.5 * inv_x
        - inv_x2 * (1.0/12.0
        - inv_x2 * (1.0/120.0
        - inv_x2 * (1.0/252.0
        - inv_x2 * (1.0/240.0
        - inv_x2 * (1.0/132.0)))));
    return result;
}
```

No external dependency (no Boost, no scipy at C++ level).

### 4.6 SQUAREM wrapper

The entire SQUAREM loop (both MAP-EM and VBEM variants) moves to C++.
The high-level structure:

```cpp
struct EMResult {
    std::vector<double> theta;      // converged normalized theta
    std::vector<double> alpha;      // un-normalized Dirichlet params
    std::vector<double> em_totals;  // final E-step accumulation
};

EMResult run_squarem(
    const std::vector<EquivClass>& ec_data,
    const double* log_eff_len,
    const double* unique_totals,
    double* prior,          // may be modified by pruning
    int n_components,
    int max_sq_iters,
    double convergence_delta,
    bool use_vbem,
    double prune_threshold, // < 0 = disabled
    int post_prune_iters
);
```

This is a direct translation of lines 870–1040 of `estimator.py`.
Since both MAP-EM and VBEM share the same SQUAREM structure (just
different step functions), we use a function pointer or template
parameter:

```cpp
using StepFn = void(*)(
    const double* state,
    const std::vector<EquivClass>& ec,
    const double* log_eff_len,
    const double* unique, const double* prior,
    double* em_totals, double* state_new,
    int n_comp
);
```

### 4.7 Coverage-weighted warm start + OVR prior

Currently in `run_locus_em()` lines 840-895. Moves to C++ since it
iterates over the same equivalence classes:

```cpp
void compute_ovr_prior(
    const std::vector<EquivClass>& ec_data,
    const bool* eligible,        // [n_components]
    const double* unique_totals, // [n_components]
    double alpha,                // flat Dirichlet pseudocount
    double gamma,                // OVR scale factor
    double* prior_out,           // [n_components] written
    double* theta_init_out       // [n_components] written
);
```

---

## 5. Build Integration

### 5.1 CMakeLists.txt addition

```cmake
# --- Build the _em_impl extension -------------------------------------------

nanobind_add_module(
  _em_impl
  NOMINSIZE
  STABLE_ABI
  src/hulkrna/native/em_solver.cpp
)

install(TARGETS _em_impl LIBRARY DESTINATION hulkrna)
```

No external dependencies (no htslib, no cgranges). Pure C++17 + nanobind + numpy.

### 5.2 Header includes

```
<nanobind/nanobind.h>
<nanobind/ndarray.h>        // for nb::ndarray<>
<nanobind/stl/tuple.h>      // for returning std::tuple
<cmath>
<cstdint>
<vector>
<unordered_map>
<algorithm>
<numeric>
```

---

## 6. Python-side Changes

### 6.1 `estimator.py` modifications

1. **Add import** at top:
   ```python
   from hulkrna._em_impl import run_locus_em_native
   ```

2. **Replace `run_locus_em()` body** — the section from `_apply_bias_correction()` call through post-prune
   (roughly lines 800-1040) becomes a single call to `run_locus_em_native()`.

3. **Delete** `_em_step()`, `_vbem_step()`, `_build_equiv_classes()`,
   `_apply_bias_correction_uniform()` functions entirely.

4. **Keep** `_apply_bias_correction()` dispatcher IF non-uniform
   `BiasProfile` support is needed. Currently only the uniform fast-path
   is used, so this can likely also be deleted.

### 6.2 Other files — no changes

- `locus.py` — unchanged (builds `LocusEMInput`, passes to `run_locus_em`)
- `pipeline.py` — unchanged (calls `run_locus_em` via `AbundanceEstimator`)
- `config.py` — unchanged (EMConfig dataclass)

---

## 7. Test Strategy

Per the project's direction: **discard existing unit tests for the
ported functions and write new ones from scratch.**

### 7.1 What to delete

- All tests that directly call `_em_step`, `_vbem_step`,
  `_build_equiv_classes`, or `_apply_bias_correction_uniform`.

### 7.2 New tests for `_em_impl`

Write a new `test_em_impl.py` that tests the C++ module directly:

**Unit tests for `run_locus_em_native()`:**

1. **Single-transcript locus** — 1 mRNA + 1 nRNA + 1 gDNA component.
   All units are unique. Theta should converge to the unique-count
   distribution.

2. **Two-transcript locus** — known ground truth where transcript A
   has 2× the evidence of B. Verify theta_A ≈ 2 × theta_B.

3. **MAP-EM vs VBEM** — same input, verify both modes converge.
   VBEM should suppress low-evidence components more aggressively.

4. **Convergence** — verify iteration count is reasonable (not hitting
   max_iterations for well-conditioned inputs).

5. **Prior eligibility** — zero-prior components should get theta → 0.

6. **Pruning** — components with no unique evidence and low evidence
   ratio should be pruned to zero.

7. **Empty locus** — zero units → returns uniform theta.

8. **Numerical stability** — log-likelihoods near -1000, verify no
   NaN/Inf in output.

9. **Bias correction** — verify that effective-length correction
   produces the expected log-lik adjustment (compare a few values
   against Python reference).

**Integration tests:**

10. **End-to-end pipeline** — run the full pipeline on a small test BAM
    and compare transcript counts against a saved reference. This is
    the ultimate correctness check.

11. **Benchmark regression** — run the PVT1/MYC benchmark and verify
    MAE/correlation are within tolerance of the Python baseline.

### 7.3 Validation approach

Before deleting the Python functions, run both implementations
side-by-side on the test suite and verify identical (within floating-point
tolerance) theta/alpha outputs. This can be a temporary validation
script, not a permanent test.

---

## 8. Estimated Complexity

| Task | Lines of C++ | Effort |
|---|---:|---|
| `em_solver.cpp` core (step functions, eq classes, bias) | ~250 | Medium |
| SQUAREM + warm start + prune | ~150 | Medium |
| nanobind module definition + array handling | ~50 | Low |
| digamma implementation | ~25 | Low |
| CMakeLists.txt | ~10 | Trivial |
| Python-side cleanup (estimator.py) | ~-200 (deletion) | Low |
| New test file | ~300 | Medium |
| **Total new C++** | **~475** | |
| **Net Python change** | **~-200** | |

Estimated implementation time: **1 session** (the EM algorithm is
well-understood, the C++ patterns are established from the 4 existing
extensions, and the inner loop is straightforward numerics).

---

## 9. Risks and Mitigations

| Risk | Mitigation |
|---|---|
| Floating-point divergence between Python and C++ | Validate side-by-side before cutting over; use identical constants (`_EM_LOG_EPSILON = 1e-300`) |
| digamma accuracy | Test against `scipy.special.digamma` for a range of inputs; asymptotic series with shift ≥ 6 gives ~15 digits |
| Large loci (thousands of units) causing memory pressure | EC scratch arrays are pre-allocated once; total memory is bounded by `n_candidates × sizeof(double)` which is already materialized in Python |
| SQUAREM extrapolation producing negative theta | Same clamp logic as Python: `std::max(theta, 0.0)` + renormalize |
| nanobind ndarray ABI issues | Stable ABI already proven with 4 existing extensions |
| Compiler auto-vectorization insufficient | Inner loop structure (contiguous row-major, no branches) is vectorization-friendly; can add `#pragma` hints or explicit SIMD later if needed |

---

## 10. Success Criteria

1. **Correctness**: Full test suite passes (671+ tests), PVT1/MYC
   benchmark MAE within 1% of Python baseline.

2. **Performance**: EM solver ≥ 5× faster (from ~8 s to < 1.6 s on the
   benchmark workload). Total pipeline time ≤ 8 s (from ~16 s).

3. **Code cleanliness**: `_em_step`, `_vbem_step`,
   `_build_equiv_classes`, `_apply_bias_correction_uniform` deleted from
   Python. No dual code paths. No fallback imports.

4. **Build**: `pip install .` succeeds on macOS ARM64 and Linux x86_64.
   No new external dependencies.
