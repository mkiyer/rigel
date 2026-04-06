# Calibration V3: Implementation Plan

## Design Summary

Replace the iterative calibration EM (~1400 lines) with an analytical gDNA–RNA deconvolution (~100 lines). Calibration's sole job is to separate gDNA from Total RNA. It does not distinguish mature from nascent RNA.

**Core outputs:**
- Per-region $\mathbb{E}[N_{\text{gDNA},r}]$ — continuous physical fragment counts
- gDNA fragment-length model
- Global gDNA background density $\lambda_G$

**Core mechanism:** Two signal pathways, blended by strand specificity:

1. **Strand pathway:** Analytical decomposition using strand balance. Dominant when $SS \to 1.0$.
2. **Density pathway:** Background density estimation from low-density regions. Dominant when $SS \to 0.5$.

Blending weight: $w = (2 \cdot SS - 1)^2$. Applied uniformly to three aspects:
- Per-region $\mathbb{E}[N_{\text{gDNA}}]$
- Global $\lambda_G$ estimation
- gDNA FL model source fragment selection

**EM solver changes:** Delete the binary `if (gamma == 0)` gate and the `VBEM_CLAMP_FLOOR = 0.1`. Pass discounted physical counts ($\alpha_{\text{gDNA}} = \beta \cdot \mathbb{E}[N_{\text{gDNA}}]$) as Dirichlet priors directly.

---

## Global $\lambda_G$ Estimation (Blended)

The global gDNA background density $\lambda_G$ is itself a blend of two estimates:

$$\lambda_G = w \cdot \lambda_G^{\text{strand}} + (1 - w) \cdot \lambda_G^{\text{density}}$$

**Strand-derived:** Aggregate the analytical decomposition across all regions with unambiguous gene strand:

$$\lambda_G^{\text{strand}} = \frac{\sum_{r:\, \text{gene\_strand} \neq 0} \mathbb{E}[N_{\text{gDNA},r}^{\text{strand}}]}{\sum_{r:\, \text{gene\_strand} \neq 0} L_r}$$

**Density-derived:** The Kth percentile of unspliced density across eligible regions:

$$\lambda_G^{\text{density}} = \text{percentile}_K\!\left(\left\{\frac{N_{\text{unspliced},r}}{L_r} : r \in \text{eligible}\right\}\right)$$

Length-weighted percentile (large regions anchor the estimate, small boundary regions have negligible influence).

This blending ensures:
- At $SS = 1.0$: $w = 1$, $\lambda_G$ is purely strand-derived (no arbitrary percentile)
- At $SS = 0.5$: $w = 0$, $\lambda_G$ is purely density-derived (only option)
- At $SS = 0.75$: $w = 0.25$, mostly density with a strand correction

The density percentile $K$ is exposed as a user parameter (`density_percentile`, default 10).

---

## Phases

### Phase 1: C++ — Add Fragment Strand to FL Table

**Goal:** The FL table currently stores `(region_id, frag_len)` per unique unspliced single-region fragment. Add `frag_strand` so calibration can select antisense fragments for the gDNA FL model.

**Files:**

1. **`src/rigel/native/resolve_context.h`** — `RegionAccumulator`

   Add a new vector alongside `fl_frag_lens`:
   ```cpp
   std::vector<int32_t> fl_frag_strands;  // STRAND_POS or STRAND_NEG
   ```

   In `accumulate()`, where `fl_frag_lens.push_back(frag_len)` is called (line ~1416), also push:
   ```cpp
   fl_frag_strands.push_back(frag_strand);
   ```

   In the `clear()` / `reset()` method, also clear `fl_frag_strands`.

2. **`src/rigel/native/bam_scanner.cpp`** — Expose `fl_frag_strands` in the result dictionary

   Where `fl_region_ids` and `fl_frag_lens` are exposed to Python (search for `"fl_frag_lens"`), add:
   ```cpp
   re["fl_frag_strands"] = ...;  // same pattern as fl_frag_lens
   ```

3. **`src/rigel/pipeline.py`** — Extract `fl_frag_strands` from scan result

   In `scan_and_buffer()` (line ~315), where `fl_table` is constructed:
   ```python
   fl_frag_strands = np.asarray(re["fl_frag_strands"], dtype=np.int32)
   fl_table = pd.DataFrame({
       "region_id": fl_region_ids,
       "frag_len": fl_frag_lens,
       "frag_strand": fl_frag_strands,  # NEW
   })
   ```

**Tests:** Verify that `fl_table` now has a `frag_strand` column. Run existing tests to ensure no regressions.

**Build:** Requires `pip install --no-build-isolation -e .` after C++ changes.

---

### Phase 2: Python — New Calibration Implementation

**Goal:** Replace the entire calibration EM with the analytical V3 algorithm.

**File: `src/rigel/calibration.py`**

#### 2a. New Output Dataclass

Replace `GDNACalibration` with:

```python
@dataclass(frozen=True)
class CalibrationResult:
    """Continuous gDNA deconvolution results.

    For every genomic region, provides the expected number of gDNA
    fragments. No binary classification, no probabilities — physical
    fragment count estimates.
    """
    region_e_gdna: np.ndarray      # (n_regions,) float64 — E[N_gDNA] per region
    region_n_total: np.ndarray     # (n_regions,) float64 — N_total per region
    gdna_fl_model: FragmentLengthModel | None  # gDNA fragment-length distribution
    lambda_gdna: float             # Global gDNA background density (frags/bp)
    strand_specificity: float      # Library SS (echoed for reference)

    def to_summary_dict(self) -> dict:
        total_e_gdna = float(self.region_e_gdna.sum())
        total_n = float(self.region_n_total.sum())
        return {
            "lambda_gdna": round(self.lambda_gdna, 8),
            "total_expected_gdna": round(total_e_gdna, 1),
            "gdna_fraction": round(total_e_gdna / max(total_n, 1.0), 4),
            "strand_specificity": round(self.strand_specificity, 4),
            "gdna_fl_observations": (
                self.gdna_fl_model.n_observations
                if self.gdna_fl_model else 0
            ),
        }
```

#### 2b. Core Calibration Function

New function `calibrate_gdna_v3()`:

```
def calibrate_gdna_v3(
    region_counts: pd.DataFrame,
    fl_table: pd.DataFrame,
    region_df: pd.DataFrame,
    strand_specificity: float,
    density_percentile: float = 10.0,
    intergenic_fl_model: FragmentLengthModel | None = None,
) -> CalibrationResult:
```

**Algorithm steps (no iteration):**

1. **`compute_region_stats()`** — unchanged, produces per-region summary arrays.

2. **Eligibility mask:** `(n_total > 0) & (region_length > 0)`. Same as current.

3. **Strand decomposition (all regions with `gene_strand != 0`):**
   ```
   N_anti_r = n_pos[gene_strand == +1]  or  n_neg[gene_strand == -1]
   e_gdna_strand_r = max(0, (N_anti_r - N_unspliced_r * (1 - SS)) / (SS - 0.5))
   ```
   Set `e_gdna_strand = 0` for regions with `gene_strand == 0` or ineligible.

4. **Global λ_G (strand-derived):**
   ```
   λ_G_strand = Σ e_gdna_strand[eligible & has_strand] / Σ L_r[eligible & has_strand]
   ```

5. **Global λ_G (density-derived):**
   ```
   d_unspliced = N_unspliced / L_r  (eligible regions)
   λ_G_density = length-weighted percentile(d_unspliced, K)
   ```

6. **Blended λ_G:**
   ```
   w = (2 * SS - 1)²
   λ_G = w * λ_G_strand + (1 - w) * λ_G_density
   ```

7. **Per-region density pathway:**
   ```
   e_gdna_density_r = min(N_unspliced_r, λ_G * L_r)
   ```

8. **Per-region blend:**
   ```
   w_r = w  if gene_strand != 0  else  0.0
   e_gdna_r = w_r * e_gdna_strand_r + (1 - w_r) * e_gdna_density_r
   ```

9. **gDNA FL model** — from `build_gdna_fl_v3()` (see 2c).

10. **Return `CalibrationResult`.**

#### 2c. gDNA FL Model Builder

New function `build_gdna_fl_v3()`:

```
def build_gdna_fl_v3(
    fl_table: pd.DataFrame,         # region_id, frag_len, frag_strand
    stats: dict[str, np.ndarray],
    eligible: np.ndarray,
    strand_specificity: float,
    intergenic_fl_model: FragmentLengthModel | None,
    max_fl: int = 1000,
    min_ess: int = 50,
) -> FragmentLengthModel | None:
```

Two source pathways, blended by the same $w = (2 \cdot SS - 1)^2$:

**Stranded source:** From `fl_table`, select unique unspliced fragments that are:
- In regions with `gene_strand != 0`
- Antisense (strand opposite to gene strand)

These are ≥(2·SS−1) pure gDNA. Build a FL histogram from these.

**Density source:** From `fl_table`, select fragments in the bottom Kth percentile of unspliced density (the same regions used for `λ_G_density`).

If $w > 0$ and enough antisense fragments exist, use the stranded source.
If $w < 1$ or insufficient stranded observations, supplement with density source fragments.
If neither source yields enough (ESS < `min_ess`), fall back to `intergenic_fl_model`.

#### 2d. Deletions

Remove the following functions and constants from `calibration.py`:
- `GDNA_INIT_DENSITY_PERCENTILE`, `GDNA_INIT_MIN_REGIONS`
- `_compute_strand_llr_binomial()`
- `_compute_strand_llr_betabinomial()`
- `estimate_kappa_marginal()`
- `_seed_initial_partition()`
- `_e_step()`
- `_m_step()`
- The old `calibrate_gdna()` function (keep `compute_region_stats()` and `compute_sense_fraction()`)
- `build_gdna_fl_model()` (replaced by `build_gdna_fl_v3()`)
- `compute_log_density()` (no longer needed)

Retain:
- `compute_region_stats()` — used by V3
- `compute_sense_fraction()` — may be useful for diagnostics
- `FragmentLengthModel` import and usage

#### 2e. Length-Weighted Percentile Helper

New utility function:

```python
def _length_weighted_percentile(
    values: np.ndarray,
    weights: np.ndarray,
    percentile: float,
) -> float:
    """Weighted percentile using linear interpolation."""
    ...
```

This ensures large regions (which provide more reliable density estimates) anchor the percentile computation. A 50,000bp region should count for more than a 50bp boundary fragment.

---

### Phase 3: Config and Pipeline Wiring

**Goal:** Update `CalibrationConfig`, the pipeline orchestration, and the locus prior computation.

#### 3a. Config (`src/rigel/config.py`)

Replace `CalibrationConfig`:

```python
@dataclass(frozen=True)
class CalibrationConfig:
    """Configuration for gDNA calibration.

    The calibration phase analytically deconvolves gDNA from Total RNA
    using strand balance and density signals. Outputs continuous expected
    gDNA fragment counts per region.
    """
    beta: float = 0.05
    """Discount factor for converting calibration estimates to EM priors.
    Prior α_gDNA = beta × E[N_gDNA]. 0.05 means calibration evidence
    carries 5% the weight of per-fragment evidence."""

    density_percentile: float = 10.0
    """Percentile of unspliced density distribution used to estimate
    the gDNA background density λ_G. Used directly for unstranded
    libraries (SS=0.5); blended with the strand-derived estimate
    for stranded libraries. Lower values are more conservative
    (underestimate gDNA). Range: 1–50."""

    min_fl_ess: int = 50
    """Minimum effective sample size for the gDNA FL model.
    Below this threshold, fall back to the intergenic FL model."""
```

Deleted fields: `max_iterations`, `convergence_tol`, `min_gdna_regions` (no longer relevant — no iteration, no minimum region count).

#### 3b. Pipeline Calibration Call (`src/rigel/pipeline.py`)

In `quant()`, replace the `calibrate_gdna()` call (~line 844) with:

```python
from .calibration import calibrate_gdna_v3

cal_cfg = config.calibration
calibration = calibrate_gdna_v3(
    region_counts,
    fl_table,
    index.region_df,
    strand_models.strand_specificity,
    density_percentile=cal_cfg.density_percentile,
    intergenic_fl_model=frag_length_models.intergenic,
)
```

Update the log line to reflect V3 output:

```python
logger.info(
    f"[CAL] gDNA calibration: λ_G={calibration.lambda_gdna:.2e}, "
    f"E[gDNA]={calibration.region_e_gdna.sum():.0f}, "
    f"SS={calibration.strand_specificity:.3f}"
)
```

#### 3c. Locus Prior Computation (`src/rigel/locus.py`)

Replace `compute_gdna_locus_gammas()` with `compute_locus_gdna_priors()`:

```python
def compute_locus_gdna_priors(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: CalibrationResult,
    beta: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-locus Dirichlet priors from calibration.

    Returns
    -------
    alpha_gdna : np.ndarray, shape (n_loci,), float64
        Discounted expected gDNA fragments per locus.
    alpha_rna : np.ndarray, shape (n_loci,), float64
        Discounted expected RNA fragments per locus.
    """
    n_loci = len(loci)
    alpha_gdna = np.empty(n_loci, dtype=np.float64)
    alpha_rna = np.empty(n_loci, dtype=np.float64)

    e_gdna = calibration.region_e_gdna
    n_total = calibration.region_n_total

    region_cr = getattr(index, "region_cr", None)
    if region_cr is None or e_gdna is None:
        alpha_gdna[:] = 0.0
        alpha_rna[:] = beta  # minimal RNA prior
        return alpha_gdna, alpha_rna

    for li, locus in enumerate(loci):
        e_gdna_sum = 0.0
        n_total_sum = 0.0
        for ref, start, end in locus.merged_intervals:
            for _s, _e, rid in region_cr.overlap(ref, start, end):
                e_gdna_sum += float(e_gdna[rid])
                n_total_sum += float(n_total[rid])

        alpha_gdna[li] = beta * e_gdna_sum
        alpha_rna[li] = beta * max(n_total_sum - e_gdna_sum, 0.0)

    return alpha_gdna, alpha_rna
```

Update `_compute_priors()` in `pipeline.py` to call the new function and return `(alpha_gdna, alpha_rna)` instead of `locus_gammas`.

#### 3d. Pipeline EM Dispatch

In `_run_locus_em_partitioned()`, change:
- `batch_gammas` → `(batch_alpha_gdna, batch_alpha_rna)` two arrays
- `_call_batch_em()` → pass both alpha arrays to the C++ solver
- `_build_locus_meta()` → log `alpha_gdna` instead of `locus_gamma`

The `estimator.run_batch_locus_em_partitioned()` signature must be updated to accept `(alpha_gdna, alpha_rna)` instead of `locus_gammas`.

---

### Phase 4: C++ EM Solver Changes

**Goal:** Accept physical-count priors. Delete the binary gate and VBEM clamp floor.

**File: `src/rigel/native/em_solver.cpp`**

#### 4a. Delete `VBEM_CLAMP_FLOOR`

Replace `static constexpr double VBEM_CLAMP_FLOOR = 0.1;` (line 59) with deletion (or set to `EM_LOG_EPSILON`).

At both clamp sites (lines ~946 and ~962):
```cpp
// OLD:
double floor_i = std::max(prior[i], VBEM_CLAMP_FLOOR);
// NEW:
double floor_i = std::max(prior[i], EM_LOG_EPSILON);
```

This prevents SQUAREM from creating negative α (invalid for Dirichlet) while allowing arbitrarily small positive values. At α = 1e-300, digamma returns ≈ −1e300. After max-subtraction in the E-step kernel, exp(−1e300) = 0.0 (IEEE underflow). Component receives zero responsibility and stays dead.

#### 4b. Modify `batch_locus_em_partitioned()` Signature

Replace `locus_gammas` with two arrays:

```cpp
batch_locus_em_partitioned(
    ...
    f64_1d   locus_alpha_gdna,   // was: locus_gammas
    f64_1d   locus_alpha_rna,    // NEW
    i64_1d   gdna_spans,
    ...
)
```

Remove `total_pseudocount` parameter — the prior budget is now the sum of `alpha_gdna + alpha_rna`, not a separate scalar.

#### 4c. Modify `compute_ovr_prior_and_warm_start()`

Change signature from `(locus_gamma, total_pseudocount, gdna_idx)` to `(alpha_gdna, alpha_rna_total, gdna_idx)`:

```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        alpha_gdna,        // was: locus_gamma
    double        alpha_rna_total,   // was: total_pseudocount
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components)
```

**Prior distribution:**
```cpp
// gDNA prior
if (gdna_idx >= 0 && gdna_idx < n_components && eligible[gdna_idx] > 0.0) {
    prior_out[gdna_idx] = std::max(alpha_gdna, EM_LOG_EPSILON);
} else if (gdna_idx >= 0) {
    prior_out[gdna_idx] = 0.0;
}

// RNA budget distributed proportionally to coverage
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] <= 0.0 || i == gdna_idx) continue;
    if (total_rna_coverage > 0.0) {
        prior_out[i] = std::max(
            alpha_rna_total * coverage_totals[i] / total_rna_coverage,
            EM_LOG_EPSILON);
    } else if (n_rna_eligible > 0) {
        prior_out[i] = std::max(
            alpha_rna_total / n_rna_eligible, EM_LOG_EPSILON);
    }
}
```

**Warm start:**
```cpp
// gDNA warm-start: proportional to alpha ratio
if (gdna_idx >= 0 && eligible[gdna_idx] > 0.0 && alpha_gdna > 0.0) {
    double total_alpha = alpha_gdna + alpha_rna_total;
    double others = total_theta - theta_init_out[gdna_idx];
    if (others > 0.0 && total_alpha > 0.0) {
        theta_init_out[gdna_idx] = (alpha_gdna / alpha_rna_total) * others;
    }
}
```

#### 4d. Delete the Binary Gate

In `extract_locus_sub_problem_from_partition()` (line ~1879):

```cpp
// DELETE:
if (locus_gamma == 0.0) {
    sub.prior[sub.gdna_idx] = 0.0;
}
```

This block is no longer needed. When `alpha_gdna = 0`, the prior computation in 4c sets `prior[gdna_idx] = EM_LOG_EPSILON ≈ 0`. The VBEM digamma naturally suppresses the component to zero responsibility. For MAP-EM, the mode formula `max(N + α − 1, 0)` with α ≈ 0 requires at least 1 fragment to be nonzero.

However, `extract_locus_sub_problem_from_partition()` must be updated to receive `alpha_gdna` and `alpha_rna_total` instead of `locus_gamma`, and pass them through to `compute_ovr_prior_and_warm_start()`.

#### 4e. Update Python Binding

In the `NB_MODULE` section, update the function signature for `batch_locus_em_partitioned()` to match the new parameters.

#### 4f. Update `src/rigel/estimator.py`

Change `run_batch_locus_em_partitioned()` to accept `(alpha_gdna_array, alpha_rna_array)` instead of `locus_gammas`. Pass them as contiguous float64 arrays to the C++ function. Remove `prior_pseudocount` from the C++ call.

---

### Phase 5: Test Updates

**Goal:** All tests pass with the new calibration and EM interface.

#### 5a. Calibration Tests

Files affected:
- `tests/test_calibration.py` — primary calibration test suite
- `tests/test_calibration_integration.py` — integration tests
- `tests/test_calibration_stress.py` — stress tests
- `tests/test_calibrated_density.py` — density-specific tests

Update approach:
- Replace `GDNACalibration` references with `CalibrationResult`
- Remove assertions on deleted fields (`kappa_strand`, `mixing_proportion`, `gdna_density_global`, `expressed_density`, `n_iterations`, `converged`)
- Add assertions on new fields (`region_e_gdna`, `lambda_gdna`)
- Tests that verify EM convergence properties are deleted (no iteration)
- Tests that verify binary γ classification → rewrite as continuous E[N_gDNA] assertions

#### 5b. EM Tests

Files affected:
- `tests/test_em_impl.py` — EM solver tests

Update approach:
- Tests calling the C++ EM with `locus_gammas` → use `(alpha_gdna, alpha_rna)` interface
- Add specific test: small α_gDNA (0.001) → VBEM naturally suppresses to ~0
- Add specific test: zero α_gDNA → component disabled, no crash
- Add specific test: large α_gDNA → component persists with appropriate share

#### 5c. Pipeline / Integration Tests

Files affected:
- `tests/test_pipeline_smoke.py`
- `tests/test_pipeline_routing.py`
- `tests/test_golden_output.py`
- `tests/scenarios_aligned/*.py`

Update approach:
- Update mock/fixture code that constructs `GDNACalibration` → `CalibrationResult`
- Update any code that references `locus_gammas` → `(alpha_gdna, alpha_rna)`
- Regenerate golden outputs: `pytest tests/ --update-golden`

#### 5d. Config Tests

Any test that constructs `CalibrationConfig` with old fields must be updated.

#### 5e. New Test Cases

Add to `tests/test_calibration.py`:

1. **Pure RNA scenario:** All reads are strand-biased → `e_gdna ≈ 0` everywhere, `lambda_gdna ≈ 0`
2. **Pure gDNA scenario:** All reads are strand-symmetric → `e_gdna ≈ n_total` everywhere
3. **Mixed scenario:** Known strand composition → verify analytical formula produces correct estimates
4. **Unstranded (SS=0.5):** Strand pathway weight = 0, density pathway dominates
5. **Intermediate SS (0.75):** Both pathways contribute
6. **gDNA FL from antisense:** Verify FL model is built from antisense unspliced fragments
7. **No antisense fragments:** FL falls back to intergenic model

---

## Implementation Order

```
Phase 1  (C++: fl_table strand)     ─── independent
Phase 4a (C++: delete clamp floor)  ─── independent of Phase 1
                                         can be done in parallel
         ↓
Phase 2  (Python: calibration.py)   ─── depends on Phase 1 (fl_table has strand)
Phase 3a (Python: config.py)        ─── independent, can parallel with Phase 2
         ↓
Phase 3b-d (Python: pipeline + locus wiring) ─── depends on Phases 2 + 3a
Phase 4b-f (C++: EM interface change)        ─── depends on Phase 3
         ↓
Phase 5  (Tests)                    ─── depends on all above
```

**Minimum compilation cycles:** Two.
1. After Phase 1 + 4a (both C++ changes, can be compiled together)
2. After Phase 4b–f (EM interface change)

Phases 2 and 3a are pure Python and need no compilation.

---

## Files Changed

| File | Phase | Change Summary |
|------|-------|----------------|
| `src/rigel/native/resolve_context.h` | 1 | Add `fl_frag_strands` vector to `RegionAccumulator` |
| `src/rigel/native/bam_scanner.cpp` | 1 | Expose `fl_frag_strands` in result dict |
| `src/rigel/pipeline.py` | 1, 3b-d | Extract `frag_strand` col; rewire calibration call and locus priors |
| `src/rigel/calibration.py` | 2 | New `CalibrationResult`, `calibrate_gdna_v3()`, `build_gdna_fl_v3()`. Delete old EM. |
| `src/rigel/config.py` | 3a | Replace `CalibrationConfig` fields |
| `src/rigel/locus.py` | 3c | Replace `compute_gdna_locus_gammas()` with `compute_locus_gdna_priors()` |
| `src/rigel/estimator.py` | 3d, 4f | Accept `(alpha_gdna, alpha_rna)` instead of `locus_gammas` |
| `src/rigel/native/em_solver.cpp` | 4 | Delete clamp floor, binary gate; new prior interface |
| `tests/test_calibration.py` | 5 | Rewrite for V3 API and assertions |
| `tests/test_calibration_integration.py` | 5 | Update for V3 |
| `tests/test_calibration_stress.py` | 5 | Update for V3 |
| `tests/test_calibrated_density.py` | 5 | Update for V3 |
| `tests/test_em_impl.py` | 5 | Update EM tests for new prior interface |
| `tests/test_pipeline_smoke.py` | 5 | Update for V3 |
| `tests/test_golden_output.py` | 5 | Regenerate golden outputs |
| `tests/scenarios_aligned/*.py` | 5 | Update fixtures for V3 |

---

## Parameters

| Parameter | Location | Default | Purpose |
|-----------|----------|---------|---------|
| `beta` | `CalibrationConfig` | 0.05 | Discount factor: α = β × E[N]. Controls prior strength relative to data. |
| `density_percentile` | `CalibrationConfig` | 10.0 | Percentile of unspliced density for λ_G estimation. Lower = more conservative. |
| `min_fl_ess` | `CalibrationConfig` | 50 | Minimum observations for gDNA FL model before falling back to intergenic. |

---

## Verification Checklist

- [ ] All existing tests pass (with updated assertions)
- [ ] Golden outputs regenerated
- [ ] `ruff check src/ tests/` clean
- [ ] Pure RNA scenario: α_gDNA ≈ 0 for all loci
- [ ] Pure gDNA scenario: α_gDNA dominates
- [ ] Unstranded (SS=0.5): density pathway produces reasonable estimates
- [ ] VBEM with tiny α: component suppressed naturally (no clamp)
- [ ] VBEM with zero α: component disabled (no gate, no crash)
- [ ] MAP-EM with same priors: consistent behavior
- [ ] Stranded FL model: built from antisense unspliced fragments
- [ ] Unstranded FL model: built from low-density region fragments
- [ ] Intermediate SS (0.75): both pathways contribute to λ_G and E[N_gDNA]
