# EM Configuration Before Simplification

> **Date**: 2025-01  
> **Purpose**: Record of all EM parameters and design decisions prior to the
> simplification refactor.  Retained for rollback reference if the simplified
> EM proves insufficient on future datasets.

## Overview

The original EM used a **hybrid VBEM** design with model re-estimation:

1. **Standard EM warm-up** (5 iterations) using `log(theta)` to seed α.
2. **VBEM convergence** using `digamma(alpha)` for sparsity-inducing
   regularization.
3. **Model re-estimation** — after convergence, re-estimate strand models
   from converged posteriors, then re-converge (2 passes total).
4. **Hybrid E-step** — `log(theta)` for transcript components,
   `digamma(alpha)` for gDNA shadow components (to preserve isoform
   estimation while inducing gDNA sparsity).

### Why it was simplified

A systematic sweep of 23 scenarios × 33 EM configurations
(`scripts/em_baseline_and_sweep.py`) showed:

- **Warmup**: All warmup values (0–20) produced identical results (MAE=30.45,
  MaxRelErr=10.55).  Warmup was effectively dead code.
- **Shadow prior**: All values (0.001–0.5) produced identical results
  (MAE=30.45 ±0.001).  A separate shadow prior was unnecessary.
- **Re-estimation**: 0 passes gave MAE=30.43 with *better* MaxRelErr=7.1
  vs baseline 10.6.  Re-estimation actively hurt.
- **Standard EM vs Hybrid/VBEM**: Pure standard EM (MAE=30.42) matched the
  hybrid baseline (MAE=30.45).  VBEM complexity was unnecessary.
- **Pure VBEM**: MAE=47.5+, caused superset collapse on two-isoform scenarios.
  VBEM was actively harmful for isoform quantification.
- **Alpha (pseudocount)**: Only ~0.2 MAE variation across 0.01–1.0 in standard
  EM family.  Not critical; 0.5 is a safe default.

---

## Module-Level Constants (`counter.py`)

| Constant | Value | Purpose |
|---|---|---|
| `_SHADOW_PRIOR` | `0.05` | Separate Dirichlet prior for gDNA shadow components (non-sparse). |
| `_VBEM_WARMUP_ITERS` | `5` | Number of standard EM iterations before switching to VBEM. |
| `_REEST_PASSES` | `2` | Number of VBEM re-estimation passes (strand model updates). |
| `_EM_LOG_EPSILON` | `1e-300` | Epsilon added to theta before log to prevent log(0). |
| `_VBEM_CONVERGENCE_DELTA` | `1e-6` | Relative convergence tolerance for alpha updates. |
| `_STRAND_PRIOR_N` | `10.0` | Pseudo-observations for strand model re-estimation. |
| `_INSERT_PRIOR_N` | `50.0` | Pseudo-observations for insert model re-estimation. |
| `_REEST_LOG_FLOOR` | `1e-10` | Numeric floor for strand log-probability clamping. |

## CLI Parameters (`cli.py`)

| Flag | Default | Purpose |
|---|---|---|
| `--em-pseudocount` | `0.5` | Dirichlet prior α₀ for VBEM. |
| `--em-iterations` | `1000` | Maximum VBEM iterations. |
| `--gdna-threshold` | `0.5` | RNA/gDNA classification threshold. |
| `--gdna-splice-penalty-unannot` | `0.01` | gDNA splice penalty for unannotated spliced reads. |
| `--confidence-threshold` | `0.95` | High-confidence posterior threshold. |

## EMData Fields (now removed)

The following `EMData` fields supported model re-estimation and have been
removed in the simplification:

```python
cand_strand_dir: np.ndarray | None = None
    # int8[n_candidates] — strand direction per candidate.
    # +1 = sense, −1 = antisense, 0 = ambiguous.

cand_is_shadow: np.ndarray | None = None
    # bool[n_candidates] — True for gDNA shadow candidates.

cand_insert_sizes: np.ndarray | None = None
    # int32[n_candidates] — insert size per candidate.

cand_splice_ll: np.ndarray | None = None
    # float64[n_candidates] — splice penalty log-likelihood per candidate.

init_p_rna_sense: float = 0.5
    # Initial RNA strand specificity (from exonic_spliced model).

init_p_gdna_sense: float = 0.5
    # Initial gDNA strand balance (from intergenic model).

init_rna_insert_ll: np.ndarray | None = None
    # float64[max_insert+1] — initial RNA insert log-probabilities.

init_gdna_insert_ll: np.ndarray | None = None
    # float64[max_insert+1] — initial gDNA insert log-probabilities.

insert_size_max: int = 1000
    # Maximum insert size for histogram dimensioning.
```

## Key Algorithm Details

### Hybrid E-Step (removed)

```python
# Transcript components: standard EM (log theta)
log_weight_t = log(theta_t)

# gDNA shadow components: VBEM (digamma alpha)
log_weight_s = digamma(alpha_s)

# Both normalized by effective length:
log_posterior ∝ log_lik + log_weight - log(eff_len)
```

### Model Re-Estimation (`_reestimate_and_rescore`, removed)

After VBEM convergence:
1. Decompose fragments into RNA/gDNA using converged posteriors.
2. Re-estimate RNA strand specificity (p_rna_sense) using Beta prior
   with `_STRAND_PRIOR_N` pseudo-observations.
3. Re-estimate insert size histograms for RNA/gDNA (monitoring only;
   not used for scoring due to self-reinforcing asymmetry).
4. Recompute `log_liks` with updated strand parameters.
5. Re-converge VBEM with updated log-likelihoods.
6. Repeat for `_REEST_PASSES` total passes.

### Per-Component Prior (removed)

```python
prior = np.full(n_total, self.vbem_prior)  # 0.5 for transcripts
prior[gdna_base:] = _SHADOW_PRIOR           # 0.05 for shadows
```

### Assign Ambiguous (changed)

Previously used `digamma(alpha)` for final posterior computation:
```python
log_weights = digamma(alpha) - np.log(eff_len)
```

Now uses `log(theta)`:
```python
log_weights = np.log(theta + _EM_LOG_EPSILON) - np.log(eff_len)
```

## Dependencies Removed

- `scipy.special.digamma` — no longer imported in `counter.py`.

## Pipeline Changes (`pipeline.py`)

The following data collection was removed from `_scan_and_build_em_data()`:

- Per-candidate strand direction tracking (`cand_strand_dir_list`)
- Per-candidate shadow flag tracking (`cand_is_shadow_list`)
- Per-unit insert size / count category expansion to per-candidate arrays
- Splice penalty log-likelihood computation per candidate
- Initial strand model parameter extraction (p_rna_sense, p_gdna_sense)
- Initial insert log-probability array computation (rna_insert_ll, gdna_insert_ll)
