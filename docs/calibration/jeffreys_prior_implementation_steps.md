# Mode-Aware Dirichlet Prior — Detailed Implementation Plan

**Date:** 2026-04-07
**Status:** Approved — ready to implement
**Theory:** See [jeffreys_corrected_prior_plan.md](jeffreys_corrected_prior_plan.md)

---

## Overview

Replace the κ×N_locus prior scaling with a mode-aware Dirichlet baseline that
eliminates VBEM digamma sparsification while remaining bias-free for MAP-EM.

**Core formula:**

```
prior[j] = baseline + γ_j × C_base

where baseline = 0.5 (VBEM) or 0.0 (MAP)
```

**Phases:**
1. Global strand aggregation (calibration fix)
2. Mode-aware Dirichlet prior (config → Python → C++)
3. Test updates and golden output regeneration
4. Validation

---

## Phase 1: Global Strand Aggregation

**File:** `src/rigel/calibration.py`

### Change 1.1: Aggregate before rectification

**Current code** (lines ~490–515):
```python
e_gdna_strand = np.maximum(
    0.0,
    (n_anti - n_unspliced * (1.0 - SS)) / denom,       # per-region max(0,...)
)
e_gdna_strand[~has_strand | ~eligible] = 0.0
e_gdna_strand = np.minimum(e_gdna_strand, n_unspliced)
```

**New code:**
```python
# Global aggregation before rectification (eliminates rectifier bias)
strand_mask = has_strand & eligible
total_anti = float(n_anti[strand_mask].sum())
total_unspliced = float(n_unspliced[strand_mask].sum())

global_e_gdna = max(0.0, (total_anti - total_unspliced * (1.0 - SS)) / denom)
global_e_gdna = min(global_e_gdna, total_unspliced)

# Distribute global estimate proportionally to per-region unspliced count
e_gdna_strand = np.zeros_like(n_unspliced)
if total_unspliced > 0 and global_e_gdna > 0:
    e_gdna_strand[strand_mask] = (
        global_e_gdna * n_unspliced[strand_mask] / total_unspliced
    )
```

**Logic:** Sum numerator/denominator across all eligible stranded regions
*before* applying max(0, ...). This eliminates the rectifier bias that
creates false γ > 0 for pure RNA at moderate SS.

### Downstream code unchanged

The subsequent `strand_e_total`, `strand_l_total`, `lambda_g_strand`
computations and the density/blending pathways remain as-is. They now
consume the globally-aggregated per-region e_gdna values.

---

## Phase 2: Mode-Aware Dirichlet Prior

### Change 2.1: CalibrationConfig (`src/rigel/config.py`)

**Remove** fields (lines 170–181):
```python
gdna_prior_kappa_min: float = 0.001
gdna_prior_kappa_max: float = 0.50
```

**Add** field:
```python
#: Calibration evidence strength for the Dirichlet prior.
#: Combined with a mode-aware baseline (+0.5 for VBEM, 0.0 for MAP)
#: per component.  The equilibrium is insensitive to this value
#: for degenerate likelihoods — any value in [1, 50] gives identical
#: results.  Default 5.0.
gdna_prior_c_base: float = 5.0
```

**Update** docstring accordingly.

### Change 2.2: compute_locus_priors() (`src/rigel/locus.py`)

**Replace** entire function (lines 386–477).

**New signature:**
```python
def compute_locus_priors(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    *,
    c_base: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
```

**New body:**
```python
    """Compute per-locus Dirichlet priors from calibration.

    For each locus, finds overlapping calibration regions via
    ``index.region_cr`` (cgranges), computes a local gDNA mixing
    fraction γ, and sets:

        α_gDNA = γ × c_base
        α_RNA  = (1 − γ) × c_base

    The C++ EM solver adds a mode-aware baseline per eligible
    component (+0.5 for VBEM, 0.0 for MAP) before running EM.
    This compensates for VBEM's digamma sparsification bias
    while remaining bias-free for MAP.

    Falls back to global γ when no regions overlap a locus.

    Returns
    -------
    alpha_gdna : np.ndarray, shape (n_loci,), float64 ≥ 0
    alpha_rna : np.ndarray, shape (n_loci,), float64 ≥ 0
    """
    n_loci = len(loci)
    alpha_gdna = np.empty(n_loci, dtype=np.float64)
    alpha_rna = np.empty(n_loci, dtype=np.float64)

    region_e_gdna = calibration.region_e_gdna
    region_n = calibration.region_n_total

    region_cr = getattr(index, "region_cr", None)
    if region_cr is None or region_e_gdna is None or region_n is None:
        # Global fallback
        total_e = float(region_e_gdna.sum()) if region_e_gdna is not None else 0.0
        total_n = float(region_n.sum()) if region_n is not None else 0.0
        gamma = total_e / max(total_n, 1.0)
        for li in range(n_loci):
            alpha_gdna[li] = gamma * c_base
            alpha_rna[li] = (1.0 - gamma) * c_base
        return alpha_gdna, alpha_rna

    # Global fallback γ for loci with no overlapping regions
    total_e = float(region_e_gdna.sum())
    total_n = float(region_n.sum())
    fallback_gamma = total_e / max(total_n, 1.0)

    for li, locus in enumerate(loci):
        e_sum = 0.0
        n_sum = 0.0
        has_overlap = False
        for ref, start, end in locus.merged_intervals:
            for _s, _e, rid in region_cr.overlap(ref, start, end):
                e_sum += float(region_e_gdna[rid])
                n_sum += float(region_n[rid])
                has_overlap = True

        if has_overlap and n_sum > 0:
            gamma = e_sum / n_sum
        else:
            gamma = fallback_gamma

        alpha_gdna[li] = gamma * c_base
        alpha_rna[li] = (1.0 - gamma) * c_base

    return alpha_gdna, alpha_rna
```

**Key removals:**
- `kappa_min`, `kappa_max` parameters
- SS/w/kappa computation block
- `N_locus = len(locus.unit_indices)` and `C = kappa * N_locus`

### Change 2.3: _compute_priors() wrapper (`src/rigel/pipeline.py`)

**Current** (lines 420–443):
```python
def _compute_priors(
    estimator: AbundanceEstimator,
    loci: list,
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    *,
    kappa_min: float = 0.001,
    kappa_max: float = 0.50,
) -> tuple[np.ndarray, np.ndarray]:
    ...
    return compute_locus_priors(
        loci, index, calibration, kappa_min=kappa_min, kappa_max=kappa_max,
    )
```

**New:**
```python
def _compute_priors(
    estimator: AbundanceEstimator,
    loci: list,
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    *,
    c_base: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    ...
    return compute_locus_priors(
        loci, index, calibration, c_base=c_base,
    )
```

### Change 2.4: quant_from_buffer() signature (`src/rigel/pipeline.py`)

**Remove** parameters (lines 687–688):
```python
gdna_prior_kappa_min: float = 0.001,
gdna_prior_kappa_max: float = 0.50,
```

**Add** parameter:
```python
gdna_prior_c_base: float = 5.0,
```

### Change 2.5: quant_from_buffer() call to _compute_priors() (`src/rigel/pipeline.py`)

**Current** (lines 780–786):
```python
alpha_gdna, alpha_rna = _compute_priors(
    estimator,
    loci,
    index,
    calibration=calibration,
    kappa_min=gdna_prior_kappa_min,
    kappa_max=gdna_prior_kappa_max,
)
```

**New:**
```python
alpha_gdna, alpha_rna = _compute_priors(
    estimator,
    loci,
    index,
    calibration=calibration,
    c_base=gdna_prior_c_base,
)
```

### Change 2.6: run_pipeline() call to quant_from_buffer() (`src/rigel/pipeline.py`)

**Current** (lines 936–937 area):
```python
gdna_prior_kappa_min=config.calibration.gdna_prior_kappa_min,
gdna_prior_kappa_max=config.calibration.gdna_prior_kappa_max,
```

**New:**
```python
gdna_prior_c_base=config.calibration.gdna_prior_c_base,
```

### Change 2.7: C++ compute_ovr_prior_and_warm_start() (`src/rigel/native/em_solver.cpp`)

**Add** `use_vbem` parameter to the function signature. Apply mode-aware
baseline per eligible component.

**Current signature** (line 772):
```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        alpha_gdna,
    double        alpha_rna,
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components)
```

**New signature:**
```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        alpha_gdna,
    double        alpha_rna,
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components,
    bool          use_vbem)
```

**Replace** the prior construction block (lines 828–847):

**Current:**
```cpp
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = std::max(alpha_gdna, EM_LOG_EPSILON);
        } else if (total_rna_coverage > 0.0) {
            prior_out[i] = std::max(
                alpha_rna * coverage_totals[i] / total_rna_coverage,
                EM_LOG_EPSILON);
        } else if (n_rna_eligible > 0) {
            prior_out[i] = std::max(alpha_rna / n_rna_eligible,
                                    EM_LOG_EPSILON);
        } else {
            prior_out[i] = EM_LOG_EPSILON;
        }
    }
```

**New:**
```cpp
    // Mode-aware baseline: VBEM needs +0.5 (Jeffreys) to cancel
    // digamma sparsification bias; MAP has no bias, needs 0.0.
    const double baseline = use_vbem ? 0.5 : 0.0;

    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = baseline + std::max(alpha_gdna, EM_LOG_EPSILON);
        } else if (total_rna_coverage > 0.0) {
            prior_out[i] = baseline + std::max(
                alpha_rna * coverage_totals[i] / total_rna_coverage,
                EM_LOG_EPSILON);
        } else if (n_rna_eligible > 0) {
            prior_out[i] = baseline + std::max(
                alpha_rna / n_rna_eligible,
                EM_LOG_EPSILON);
        } else {
            prior_out[i] = baseline;
        }
    }
```

### Change 2.8: gDNA gate override (`src/rigel/native/em_solver.cpp`)

**Current** (lines ~1871–1873):
```cpp
    if (alpha_gdna <= 0.0) {
        sub.prior[sub.gdna_idx] = 0.0;
    }
```

This code runs **before** `compute_ovr_prior_and_warm_start`, setting
`sub.prior[gdna] = 0.0` which feeds into `eligible[gdna] = 0.0`. The
`eligible` array then prevents the baseline from being added in the new
code (the `eligible[i] <= 0.0` branch). **No change needed** — the gate
already works correctly because ineligible components get `prior = 0.0`
before `compute_ovr_prior_and_warm_start` is called.

Wait — re-reading: `sub.prior` is initialized to `EM_PRIOR_EPSILON` for
all components, then gDNA is set to 0.0 if alpha_gdna <= 0. Then
`sub.eligible[c] = (sub.prior[c] > 0.0) ? 1.0 : 0.0`. So effectively
the gate disables gDNA via the eligible array, and
`compute_ovr_prior_and_warm_start` sees `eligible[gdna] = 0.0` and sets
`prior_out[gdna] = 0.0`. **No change needed.**

### Change 2.9: Update call sites of compute_ovr_prior_and_warm_start()

**Call site 1** (line 1201, `run_locus_em_native`):
```cpp
    compute_ovr_prior_and_warm_start(
        ec_data, unambig_totals.data(), pe_ptr,
        0.0, alpha_rna, -1,
        prior.data(), theta_init.data(), n_components);
```

**New** — add `use_vbem` argument. This function already receives `use_vbem`
at line 1105:
```cpp
    compute_ovr_prior_and_warm_start(
        ec_data, unambig_totals.data(), pe_ptr,
        0.0, alpha_rna, -1,
        prior.data(), theta_init.data(), n_components,
        use_vbem);
```

**Call site 2** (line 2118, `batch_locus_em_partitioned`):
```cpp
            compute_ovr_prior_and_warm_start(
                ec_data,
                sub.unambig_totals.data(),
                sub.eligible.data(),
                ag_ptr[li], ar_ptr[li], sub.gdna_idx,
                prior.data(), theta_init.data(), nc);
```

**New** — add `use_vbem` (already in scope at line 2020):
```cpp
            compute_ovr_prior_and_warm_start(
                ec_data,
                sub.unambig_totals.data(),
                sub.eligible.data(),
                ag_ptr[li], ar_ptr[li], sub.gdna_idx,
                prior.data(), theta_init.data(), nc,
                use_vbem);
```

### Change 2.10: Recompile

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

---

## Phase 3: Test Updates

### Change 3.1: test_locus_priors.py

**Full rewrite.** Remove all κ-related tests. New test structure:

```python
class TestCBaseComputation:
    """Verify α_gDNA = γ × c_base, α_RNA = (1−γ) × c_base."""

    def test_basic_split(self):
        """α_gDNA + α_RNA = c_base."""

    def test_pure_rna(self):
        """γ=0 → α_gDNA=0, α_RNA=c_base."""

    def test_pure_gdna(self):
        """γ=1 → α_gDNA=c_base, α_RNA=0."""

    def test_c_base_scales_linearly(self):
        """c_base=10 gives 2× priors vs c_base=5."""

    def test_no_dependence_on_N_locus(self):
        """Priors independent of locus size (no κ×N)."""


class TestRegionOverlap:
    """Verify per-locus γ from cgranges overlap."""

    def test_single_region(self): ...
    def test_multiple_regions(self): ...
    def test_no_overlap_fallback(self): ...
    def test_no_cgranges_fallback(self): ...


class TestEdgeCases:
    def test_empty_locus(self): ...
    def test_zero_total_fragments(self): ...
```

### Change 3.2: test_calibration_integration.py

Replace:
```python
assert cfg.gdna_prior_kappa_min == 0.001
assert cfg.gdna_prior_kappa_max == 0.50
```

With:
```python
assert cfg.gdna_prior_c_base == 5.0
```

Remove any `kappa_min`/`kappa_max` references. Add test for `c_base`:
```python
def test_custom_c_base(self):
    cfg = CalibrationConfig(gdna_prior_c_base=10.0)
    assert cfg.gdna_prior_c_base == 10.0
```

### Change 3.3: Regenerate golden outputs

```bash
pytest tests/ --update-golden
```

### Change 3.4: Run full test suite

```bash
pytest tests/ -v
```

Expected: all tests pass, including the 4 previously-failing scenario tests:
- `test_nrna_double_counting: g0_n0_s65, g0_n0_s90`
- `test_antisense_intronic: ss_0.65, ss_0.9`

---

## Phase 4: Validation

### Run synthetic benchmarks

```bash
for cfg in scripts/benchmark/configs/*.yaml; do
  name=$(basename "$cfg" .yaml)
  python scripts/synthetic_sim_sweep.py -c "$cfg" -o "scripts/benchmark/golden/$name"
done
python scripts/benchmark/analyze_golden.py
python scripts/benchmark/analyze_deep.py
```

### Check for regressions

Compare mRNA/nRNA/gDNA relative error against previous golden results.
Focus on:
- nRNA siphon: should be reduced (smaller priors)
- gDNA leakage: should be unchanged (gDNA gate still active for γ=0)
- Pure RNA scenarios: should now correctly have γ=0

---

## File Change Summary

| File | Action |
|------|--------|
| `src/rigel/calibration.py` | Modify strand aggregation (~10 lines) |
| `src/rigel/config.py` | Remove 2 fields, add 1 field |
| `src/rigel/locus.py` | Rewrite `compute_locus_priors()` (simpler: remove κ/SS/w/N logic) |
| `src/rigel/pipeline.py` | s/kappa_min,kappa_max/c_base/ in 4 locations |
| `src/rigel/native/em_solver.cpp` | Add `use_vbem` param + baseline in `compute_ovr_prior_and_warm_start()`, update 2 call sites |
| `tests/test_locus_priors.py` | Full rewrite (simpler tests) |
| `tests/test_calibration_integration.py` | Update config field assertions |
| `tests/golden/` | Regenerate |

**Estimated C++ lines changed:** ~15
**Estimated Python lines changed:** ~100
**Estimated test lines changed:** ~200
