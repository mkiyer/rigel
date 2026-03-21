# Post-EM Discrete Fragment Assignment — Implementation Plan

**Date:** 2026-03-20  
**Status:** Implemented

## Motivation

The EM algorithm produces optimal *fractional* transcript abundances, but each sequenced molecule originates from exactly one biological source. A post-EM discretization step assigns each fragment to a single transcript/nRNA/gDNA component, yielding integer counts that are:

- **Interpretable**: downstream count-based tools (DESeq2, edgeR) assume integer inputs
- **Coherent**: total fragment count is preserved exactly — no fractional fragments
- **Faithful**: each molecule is attributed to exactly one source, consistent with the underlying biology

The EM itself remains fractional — only the final output step changes.

## Design

### Assignment modes

Three modes, selected by `--assignment-mode`:

| Mode | Behavior | Deterministic? | Default? |
|------|----------|---------------|----------|
| `fractional` | Current behavior. Scatter posterior probabilities to accumulators (no discretization). | Yes | No |
| `map` | Maximum A Posteriori. Each fragment assigned to its highest-posterior component. | Yes | No |
| `sample` | Multinomial sampling. Each fragment drawn from its posterior distribution using a seeded PRNG. | Yes (given seed) | **Yes** |

**Default: `sample`** — unbiased in expectation (averages to the fractional EM counts over many seeds), while producing meaningful integer counts per run.

### Posterior thresholding

Before assignment (in `map` or `sample` modes), candidates with posterior below `--assignment-min-posterior` (default 0.01) are zeroed and the remaining posteriors are renormalized. This prevents pathological assignments to negligible-probability components.

Not applied in `fractional` mode (backward-compatible).

### Parameters

| CLI flag | Config field | Type | Default | Description |
|----------|-------------|------|---------|-------------|
| `--assignment-mode` | `em.assignment_mode` | str | `"sample"` | `"fractional"`, `"map"`, or `"sample"` |
| `--assignment-min-posterior` | `em.assignment_min_posterior` | float | `0.01` | Minimum posterior to be eligible for assignment (map/sample modes only) |

The existing `--seed` parameter is used for the `sample` mode PRNG. When seed is `None`, sample mode uses the current timestamp (existing behavior for seed).

## Integration points

### 1. C++ `assign_posteriors()` — em_solver.cpp (core change)

The `assign_posteriors` function (line 1423) already iterates per-unit, computes the full posterior vector via log-sum-exp, and scatters to accumulators. The modification is surgical:

**New parameters added to `assign_posteriors`:**
```cpp
int    assignment_mode,       // 0=fractional, 1=map, 2=sample
double min_posterior,         // threshold (map/sample only)
uint64_t rng_seed             // base seed for sample mode
```

**Logic after posterior normalization (per-unit):**

```
if mode == fractional:
    scatter p[j] to accumulators (current behavior)
else:
    1. zero candidates with p[j] < min_posterior, renormalize
    2. if mode == map:  chosen_j = argmax(p)
       if mode == sample: chosen_j = draw from categorical(p) using PRNG
    3. scatter 1.0 to chosen_j's accumulator, 0.0 to all others
```

**PRNG strategy:** SplitMix64 — fast, stateless, excellent statistical properties. Seed = `rng_seed XOR unit_global_index` per fragment for determinism independent of thread scheduling and locus processing order. Each fragment gets a unique deterministic RNG state derived from the global seed and its position.

**High-confidence semantics under discrete modes:**
- `fractional`: unchanged (max mRNA posterior ≥ confidence_threshold gates all mRNA posteriors)
- `map`/`sample`: if the winning component is mRNA and its posterior ≥ confidence_threshold, the unit contributes 1.0 to `em_high_conf_counts`

**Posterior tracking semantics:**
- `fractional`: `posterior_sum[t] += p²`, `n_assigned[t] += p` (current)
- `map`/`sample`: `posterior_sum[t] += p_winner`, `n_assigned[t] += 1.0` (posterior of the assigned component)

### 2. C++ `batch_locus_em()` — signature and plumbing

Add three new parameters to `batch_locus_em`:
```cpp
int      assignment_mode,          // 0=fractional, 1=map, 2=sample
double   assignment_min_posterior,  // min posterior threshold
uint64_t rng_seed                   // base seed for sampling
```

These are threaded through to `process_locus` → `assign_posteriors`. The `rng_seed` is passed directly; each unit derives its own seed deterministically.

In the `process_locus` lambda, pass the new parameters to `assign_posteriors`:
```cpp
assign_posteriors(
    sub, result.theta.data(), confidence_threshold,
    assignment_mode, min_posterior, rng_seed,  // NEW
    em_out, hc_out, nrna_out, gdna_out,
    psum_out, nass_out,
    locus_mrna, locus_nrna, locus_gdna,
    N_T, N_NRNA, N_COLS);
```

### 3. C++ nanobind bindings

Add three new `nb::arg()` entries to `batch_locus_em` binding:
```cpp
nb::arg("assignment_mode"),
nb::arg("assignment_min_posterior"),
nb::arg("rng_seed"),
```

`run_locus_em_native` is **not modified** — it does not call `assign_posteriors` (used only in unit tests that verify theta/alpha/em_totals, not assignment).

### 4. Python `EMConfig` — config.py

Add two new fields:
```python
assignment_mode: str = "sample"
assignment_min_posterior: float = 0.01
```

Update `__post_init__` validation:
```python
if self.assignment_mode not in ("fractional", "map", "sample"):
    raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")
```

### 5. Python `AbundanceEstimator.run_batch_locus_em()` — estimator.py

Pass the new config values to the C++ call:
```python
# Assignment config
_ASSIGNMENT_MODE_MAP = {"fractional": 0, "map": 1, "sample": 2}
...
_ASSIGNMENT_MODE_MAP[self.em_config.assignment_mode],
self.em_config.assignment_min_posterior,
np.uint64(effective_seed),
```

The `effective_seed` is derived from `self.em_config.seed` (or timestamp if None).

### 6. CLI — cli.py

Add to `_PARAM_SPECS`:
```python
_ParamSpec("assignment_mode", "em.assignment_mode"),
_ParamSpec("assignment_min_posterior", "em.assignment_min_posterior"),
```

Add CLI arguments to `quant_parser` (common options group):
```python
quant_parser.add_argument(
    "--assignment-mode", dest="assignment_mode",
    choices=["fractional", "map", "sample"], default=None,
    help="Post-EM fragment assignment mode. 'fractional' preserves EM "
         "posterior weights (traditional). 'map' assigns each fragment to "
         "its highest-posterior component. 'sample' draws from the posterior "
         "distribution (default: sample).",
)
```

And in advanced options:
```python
adv.add_argument(
    "--assignment-min-posterior", dest="assignment_min_posterior",
    type=float, default=None,
    help="Minimum posterior for a component to be eligible for discrete "
         "assignment (map/sample modes). Default: 0.01.",
)
```

## SplitMix64 RNG

Inline implementation in em_solver.cpp (public domain algorithm):

```cpp
struct SplitMix64 {
    uint64_t state;
    explicit SplitMix64(uint64_t seed) : state(seed) {}
    uint64_t next() {
        uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    // Uniform double in [0, 1)
    double uniform() {
        return (next() >> 11) * 0x1.0p-53;
    }
};
```

Per-unit seeding: `SplitMix64 rng(rng_seed ^ unit_global_offset)` where `unit_global_offset` is derived from the unit's position in the global CSR. This ensures:
- Determinism given a fixed seed (independent of thread count or locus scheduling)
- No correlation between nearby fragments

## Categorical draw from posterior

```cpp
int categorical_draw(const double* posteriors, int n, SplitMix64& rng) {
    double u = rng.uniform();
    double cumsum = 0.0;
    for (int j = 0; j < n - 1; ++j) {
        cumsum += posteriors[j];
        if (u < cumsum) return j;
    }
    return n - 1;  // numerical safety
}
```

## Test plan

1. **Unit test — MAP mode**: 2-component locus with known posteriors. Verify MAP always picks argmax. Integer-valued output.
2. **Unit test — sample mode**: Run sample assignment many times with different seeds. Verify empirical distribution matches posteriors (chi-squared test).
3. **Unit test — min_posterior threshold**: Verify components below threshold are never assigned.
4. **Unit test — fractional mode**: Verify backward compatibility (identical to current behavior).
5. **Integration test**: Pipeline smoke test with each mode.
6. **Golden output**: Regenerate with default mode (`sample`).

## Files modified

| File | Change |
|------|--------|
| `src/rigel/native/em_solver.cpp` | SplitMix64 struct, modify `assign_posteriors`, `batch_locus_em` signature, nanobind bindings |
| `src/rigel/config.py` | Add `assignment_mode`, `assignment_min_posterior` to `EMConfig` |
| `src/rigel/estimator.py` | Pass new params to `batch_locus_em` |
| `src/rigel/cli.py` | Add `_ParamSpec` entries and CLI arguments |
| `tests/test_em_impl.py` | New `TestDiscreteAssignment` class |
| `docs/METHODS.md` | Document discrete assignment step |
| `docs/MANUAL.md` | Document new CLI options |
| `docs/parameters.md` | Add parameter entries |

## Implementation order

1. `EMConfig` + CLI (Python-only, no compilation)
2. C++ `assign_posteriors` modification (SplitMix64, mode dispatch, thresholding)
3. C++ `batch_locus_em` signature + nanobind bindings
4. Python `estimator.py` plumbing
5. Tests
6. Recompile + full test suite
7. Regenerate golden outputs
8. Documentation updates
