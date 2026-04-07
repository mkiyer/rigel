# Implementation Plan: Anchored Empirical Bayes Prior + Honest gDNA Flanking

**Status**: Approved for implementation  
**Date**: 2026-04-07

## Overview

Two coordinated changes to resolve gDNA/nRNA degeneracy at low strand specificity:

1. **Phase 1 (Prior)**: Dynamic pseudocount budget $C = \kappa \cdot N_{\text{locus}}$ where $\kappa$ is inversely proportional to strand specificity
2. **Phase 2 (Likelihood)**: Extend gDNA effective length by one mean gDNA fragment length (honest flank)

No C++ changes required. All modifications are Python-only.

---

## Phase 1: Dynamic Prior Strength

### Mathematical Specification

The strand reliability weight (already used in calibration blending):
$$w = (2 \cdot SS - 1)^2$$

Prior strength interpolation:
$$\kappa = \kappa_{\min} + (1 - w) \cdot (\kappa_{\max} - \kappa_{\min})$$

Per-locus Dirichlet pseudocount budget:
$$C = \kappa \cdot N_{\text{locus}}$$

Per-locus priors (γ from calibration, unchanged):
$$\alpha_{\text{gDNA}} = \gamma \cdot C, \quad \alpha_{\text{RNA}} = (1 - \gamma) \cdot C$$

**Behavioral table** (N_locus = 2000, γ = 0.05):

| SS | w | κ | C | α_gDNA | α_RNA | Prior weight |
|----|---|---|---|--------|-------|--------------|
| 1.00 | 1.00 | 0.001 | 2.0 | 0.10 | 1.90 | 0.1% |
| 0.95 | 0.81 | 0.10 | 190 | 9.5 | 181 | 9% |
| 0.90 | 0.64 | 0.18 | 360 | 18 | 342 | 15% |
| 0.65 | 0.09 | 0.46 | 910 | 46 | 865 | 31% |
| 0.50 | 0.00 | 0.50 | 1000 | 50 | 950 | 33% |

At SS=0.50 (degenerate), the EM converges to θ_gDNA = γ = 0.05, allocating exactly 100 fragments to gDNA. At SS=1.0, the prior is comparable to current C=1.0 and the strand signal dominates.

### File: `src/rigel/config.py`

**Remove** `total_pseudocount: float = 1.0` from `CalibrationConfig`.

**Add** two new fields:

```python
@dataclass(frozen=True)
class CalibrationConfig:
    ...

    #: Prior strength at SS=1.0 (perfectly stranded).  Fraction of locus
    #: fragment count used as Dirichlet pseudocount budget.  The strand
    #: signal in the likelihood resolves gDNA locally, so the prior is
    #: a light touch.  Value of 0.001 gives C ≈ 2 for a 2000-fragment
    #: locus, comparable to the previous fixed C=1.0.
    prior_kappa_min: float = 0.001

    #: Prior strength at SS=0.50 (unstranded).  Fraction of locus
    #: fragment count.  When the likelihood is degenerate (gDNA ≈ nRNA),
    #: the prior anchors gDNA allocation to the calibrated mixing
    #: fraction γ.  Value of 0.5 gives the prior ~33% influence for a
    #: 2000-fragment locus, enough to stably resolve the degeneracy.
    prior_kappa_max: float = 0.50
```

### File: `src/rigel/locus.py` — `compute_locus_priors()`

**Replace** the current function signature and implementation.

Current signature:
```python
def compute_locus_priors(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    total_pseudocount: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
```

New signature:
```python
def compute_locus_priors(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: "CalibrationResult",
    *,
    kappa_min: float = 0.001,
    kappa_max: float = 0.50,
) -> tuple[np.ndarray, np.ndarray]:
```

New docstring:
```
Compute per-locus Dirichlet priors from calibration.

For each locus, computes a local gDNA mixing fraction γ from
overlapping calibration regions, then sets:

    C = κ × N_locus
    α_gDNA = γ × C
    α_RNA  = (1 − γ) × C

where κ interpolates between κ_min (at SS=1.0, strand signal
strong) and κ_max (at SS=0.50, likelihood degenerate):

    w = (2·SS − 1)²
    κ = κ_min + (1 − w) · (κ_max − κ_min)

When the likelihood is degenerate (identical gDNA/nRNA
per-fragment scores), the EM converges to θ_gDNA = γ
regardless of κ.  The κ value controls noise dampening
and convergence speed in the partially-degenerate case.
```

New core logic (inside the `for li, locus in enumerate(loci)` loop):
```python
    N_locus = len(locus.unit_indices)
    C = kappa * N_locus

    if has_overlap and n_sum > 0:
        gamma = e_sum / n_sum
    else:
        gamma = fallback_gamma

    alpha_gdna[li] = gamma * C
    alpha_rna[li] = (1.0 - gamma) * C
```

Compute `kappa` once before the loop:
```python
    SS = calibration.strand_specificity
    w = (2.0 * SS - 1.0) ** 2
    kappa = kappa_min + (1.0 - w) * (kappa_max - kappa_min)
```

### File: `src/rigel/pipeline.py` — `_compute_priors()`

**Change** signature from `total_pseudocount` to `kappa_min`/`kappa_max`:

Current:
```python
def _compute_priors(
    estimator, loci, index, calibration,
    total_pseudocount: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    ...
    return compute_locus_priors(
        loci, index, calibration, total_pseudocount=total_pseudocount,
    )
```

New:
```python
def _compute_priors(
    estimator, loci, index, calibration,
    kappa_min: float = 0.001,
    kappa_max: float = 0.50,
) -> tuple[np.ndarray, np.ndarray]:
    ...
    return compute_locus_priors(
        loci, index, calibration,
        kappa_min=kappa_min,
        kappa_max=kappa_max,
    )
```

### File: `src/rigel/pipeline.py` — `quant_from_buffer()`

**Remove** `calibration_total_pseudocount: float = 1.0` from signature.

**Add** `prior_kappa_min` and `prior_kappa_max` keyword arguments:

```python
def quant_from_buffer(
    ...
    calibration: "CalibrationResult" = None,
    emit_locus_stats: bool = False,
    prior_kappa_min: float = 0.001,
    prior_kappa_max: float = 0.50,
    fl_prior_ess: float | None = None,
) -> AbundanceEstimator:
```

**Change** the call site inside `quant_from_buffer`:
```python
    alpha_gdna, alpha_rna = _compute_priors(
        estimator, loci, index,
        calibration=calibration,
        kappa_min=prior_kappa_min,
        kappa_max=prior_kappa_max,
    )
```

### File: `src/rigel/pipeline.py` — `run_pipeline()`

**Change** the `quant_from_buffer()` invocation:

Current:
```python
    estimator = quant_from_buffer(
        ...
        calibration_total_pseudocount=config.calibration.total_pseudocount,
        fl_prior_ess=config.calibration.fl_prior_ess,
    )
```

New:
```python
    estimator = quant_from_buffer(
        ...
        prior_kappa_min=config.calibration.prior_kappa_min,
        prior_kappa_max=config.calibration.prior_kappa_max,
        fl_prior_ess=config.calibration.fl_prior_ess,
    )
```

### Preservation of gDNA Gate

The existing gate in `build_locus_em_data()` (line ~352):
```python
if gdna_init == 0.0:
    prior[gdna_idx] = 0.0
```

And the C++ equivalent (em_solver.cpp line ~1870):
```cpp
if (alpha_gdna <= 0.0) {
    sub.prior[sub.gdna_idx] = 0.0;
}
```

These remain unchanged. When α_gDNA = γ × κ × N = 0 (because γ = 0, i.e., calibration found zero gDNA), the gate fires as before. The dynamic C doesn't create gDNA from nothing.

---

## Phase 2: Honest gDNA Flanking

### Mathematical Specification

The gDNA effective length extends beyond the locus boundary by up to one mean gDNA fragment length:

$$L_{\text{gDNA}} = L_{\text{locus}} + F_{\text{mean}}$$

where $F_{\text{mean}}$ is the calibrated gDNA fragment length mean (typically ~200bp).

The flanking is applied to `bias_profiles[gdna_idx]` which determines the denominator
in the bias correction: `eff_len = max(L_gDNA - footprint + 1, 1)`.

**Per-fragment impact** (200bp fragments, 200bp flank):

| Locus span | nRNA eff_len | gDNA eff_len | Difference (log) |
|-----------|-------------|-------------|-----------------|
| 1,000 bp | 801 | 1,001 | 0.22 / fragment |
| 5,000 bp | 4,801 | 5,001 | 0.04 / fragment |
| 10,000 bp | 9,801 | 10,001 | 0.02 / fragment |
| 50,000 bp | 49,801 | 50,001 | 0.004 / fragment |

Most impactful for short loci. For large loci, the prior (Phase 1) does the heavy lifting.

### Biological Justification

A gDNA fragment centered at position P has a footprint of ~F_mean bp. If P is within F_mean/2 of the locus boundary, the fragment overlaps the locus and is captured by the BAM scanner. The furthest a captured fragment can start outside the locus is F_mean − 1 bp. Using one mean fragment length as the flank is conservative (the expected overhang is F_mean/2) and avoids:

- Arbitrary extension choices (no magic numbers — the flank IS the data)
- Hybrid capture boundary violations (probes cover at least the transcript body ± fragment length)
- Locus merging (a 200bp flank won't bridge to the next gene)

### Implementation

Two options for where to apply the flank:

**Option A: In pipeline.py where spans are assembled.**

This is the simplest. The `gdna_span` values are assembled into `np.int64` arrays right before the C++ batch EM call. Adding the flank here requires no changes to `locus.py` or C++.

In `_run_locus_em_partitioned()`, where spans are built:

```python
# Current (line ~591, mega-loci):
np.array([locus.gdna_span], dtype=np.int64)

# Current (line ~627, normal loci):
normal_spans = np.array([loc.gdna_span for loc in normal_loci], dtype=np.int64)
```

Change to:
```python
# Mega-loci:
np.array([locus.gdna_span + gdna_flank], dtype=np.int64)

# Normal loci:
normal_spans = np.array(
    [loc.gdna_span + gdna_flank for loc in normal_loci], dtype=np.int64
)
```

Where `gdna_flank` is passed as a parameter to `_run_locus_em_partitioned()`.

**Option B: Store on Locus.gdna_span itself.**

Modify `build_loci()` to incorporate the flank when computing `gdna_span`. This is conceptually cleaner but changes the Locus dataclass semantics.

**Recommendation**: Option A. It's minimal, reversible, and keeps the Locus dataclass reporting the true genomic span (useful for diagnostics).

### Data Flow for gdna_flank

The gDNA FL mean is available from `calibration.gdna_fl_model.mean`. Compute the flank in `quant_from_buffer()`:

```python
# After calibration injection, before scoring:
gdna_flank = int(calibration.gdna_fl_model.mean) if (
    calibration and calibration.gdna_fl_model
) else 0
```

Pass `gdna_flank` as a parameter through `_run_locus_em_partitioned()`:

Current signature (line ~501):
```python
def _run_locus_em_partitioned(
    estimator, partitions, loci, index,
    alpha_gdna, alpha_rna, em_config,
    *, emit_locus_stats=False, annotations=None,
):
```

New signature:
```python
def _run_locus_em_partitioned(
    estimator, partitions, loci, index,
    alpha_gdna, alpha_rna, em_config,
    *, emit_locus_stats=False, annotations=None,
    gdna_flank: int = 0,
):
```

---

## Test Changes

### `tests/test_locus_priors.py`

This file has extensive tests for `compute_locus_priors()`. All tests currently pass `total_pseudocount=X`. These must be rewritten to use the new `kappa_min`/`kappa_max` API.

**Key changes:**

1. **`test_budget_sums_to_C_custom`**: Currently checks `alpha_g + alpha_r == 5.0` with `total_pseudocount=5.0`. Under the new API, the budget depends on `kappa × N_locus`. The mock locus has `unit_indices=np.array([0])` → N=1. Rewrite to test `alpha_g + alpha_r == kappa × 1` with explicit kappa values.

2. **`TestBudgetSensitivity::test_doubling_C_doubles_both_alphas`**: Currently compares `total_pseudocount=1.0` vs `2.0`. Rewrite to compare `kappa_max=0.5` vs `kappa_max=1.0` (with SS=0.50 so kappa=kappa_max).

3. **All tests**: The mock `_make_locus()` helper creates loci with `unit_indices=np.array([0])` (N=1). This means C = kappa × 1 = kappa. Tests should either:
   - Set predictable kappa values and verify math, OR
   - Expand `_make_locus()` to accept a configurable `n_units` parameter

4. **New test**: `test_kappa_ss_interpolation` — verify that kappa varies correctly with SS:
   - At SS=1.0: kappa = kappa_min
   - At SS=0.50: kappa = kappa_max
   - Intermediate SS produces intermediate kappa
   - Verify C = kappa × N_locus

5. **New test**: `test_degenerate_convergence` — verify that with identical per-fragment likelihoods, the EM converges to γ:
   - Create a locus with N unspliced fragments all having identical gDNA LL and RNA LL
   - Run EM with dynamic prior
   - Verify θ_gDNA ≈ γ

### `tests/test_calibration_integration.py`

Line 96: `assert cfg.total_pseudocount == 1.0` → Replace with:
```python
assert cfg.prior_kappa_min == 0.001
assert cfg.prior_kappa_max == 0.50
```

Line 108: `cfg.total_pseudocount = 0.5` (frozen test) → Replace with:
```python
cfg.prior_kappa_min = 0.01
```

### `tests/scenarios/test_nrna_double_counting.py`

The scenario tests at SS=0.65 should improve (stronger prior aids gDNA resolution). The SS=1.0 tests should see minimal change. Golden outputs will need regeneration.

**New scenario dimensions**: Consider adding SS=0.50 to the sweep grid to verify unstranded behavior. This is a new `SS_LEVELS = [0.50, 0.65, 0.9, 1.0]` entry.

### `scripts/debug/trace_em_likelihoods.py`

Line 106: `total_pseudocount=config.calibration.total_pseudocount` → update to new parameter names.

### Golden Outputs

All golden outputs (`tests/golden/`) will change because the prior strength is different. Regenerate with `pytest tests/ --update-golden` after all changes pass.

---

## Migration Checklist

### Source Files

| File | Change | Phase |
|------|--------|-------|
| `src/rigel/config.py` | Remove `total_pseudocount`, add `prior_kappa_min`, `prior_kappa_max` | 1 |
| `src/rigel/locus.py` | Rewrite `compute_locus_priors()` signature + logic | 1 |
| `src/rigel/pipeline.py` | Update `_compute_priors()`, `quant_from_buffer()`, `run_pipeline()` | 1 |
| `src/rigel/pipeline.py` | Add `gdna_flank` param to `_run_locus_em_partitioned()`, apply to spans | 2 |

### Test Files

| File | Change | Phase |
|------|--------|-------|
| `tests/test_locus_priors.py` | Rewrite for new API, add kappa/SS tests | 1 |
| `tests/test_calibration_integration.py` | Update config default assertions | 1 |
| `tests/scenarios/test_nrna_double_counting.py` | Adjust tolerances, optionally add SS=0.50 | 1 |
| `tests/golden/*` | Regenerate all golden outputs | 1+2 |

### Script Files

| File | Change | Phase |
|------|--------|-------|
| `scripts/debug/trace_em_likelihoods.py` | Update parameter name | 1 |

### No Changes Required

| File | Reason |
|------|--------|
| `src/rigel/native/em_solver.cpp` | Receives alpha_gdna/alpha_rna as before; C++ code unchanged |
| `src/rigel/native/scoring.cpp` | Per-fragment likelihoods unchanged |
| `src/rigel/calibration.py` | Calibration logic unchanged; γ computation unchanged |
| `src/rigel/frag_length_model.py` | Dirichlet FL prior unchanged |
| `src/rigel/strand_model.py` | Strand model unchanged |

---

## Verification Strategy

1. **Unit tests**: All `test_locus_priors.py` tests pass with new API
2. **Scenario tests**: All 36 nRNA double-counting tests pass (including SS=0.65)
3. **Full suite**: All 1009 tests pass
4. **Diagnostic**: Run `trace_fl_models.py` / `trace_em_likelihoods.py` at SS=0.50 to verify EM converges to γ
5. **Regression**: Golden output regeneration, confirm diffs are consistent with stronger prior

---

## Future Work (Deferred)

- **Locoregional density prior**: Use region partition data within a sliding window (e.g., 1Mb) to compute local γ instead of global. Captures CNV-driven variation.
- **Intergenic fragment utilization**: Route intergenic fragments near locus boundaries to the locus EM as gDNA-only candidates. Provides direct gDNA evidence.
- **Intra-RNA effective-length weighting**: Switch from coverage-based to effective-length-based α_RNA splitting among transcripts.
- **SS=0.50 scenario tests**: Add unstranded dimension to the nRNA double-counting sweep grid.
