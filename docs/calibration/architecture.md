# gDNA Calibration Architecture

## Overview

The gDNA calibration system (`calibrate_gdna()` in `calibration.py`) deconvolves
genomic DNA contamination from RNA signal using a two-component EM mixture model
operating on genomic regions. Three independent likelihood-ratio channels—density,
strand balance, and fragment length—are fused into per-region posterior
probabilities γ_r ∈ [0, 1] representing P(not expressed | data).

## Standard Pathway (5 Phases)

### Phase 1: Summary Statistics
`compute_region_stats()` aggregates per-region fragment tallies:
- `n_unspliced`, `n_spliced`, `n_total`, `strand_ratio`, `splice_rate`, `density`
- `gene_strand`: +1 (sense only), −1 (antisense only), 0 (ambiguous)

### Phase 2: Preprocessing
- `sense_frac`: gene-strand-normalized sense fractions (NaN when unavailable)
- `log_d`: log-density with region-adaptive pseudocount ε = 1/median(L)
- **Eligible mask**: `(n_total > 0) & (L > 0)`

### Phase 3: Initialization (Two-Phase Seed Partition)

**Expressed Seed (3A):**
- Hard constraint: any region with `n_spliced > 0` → γ = 0 (expressed)
- Fallback if < 2 spliced regions: median eligible log-density threshold, top 50th percentile as expressed seed

**gDNA Seed (3B):**
- Compute eCDF of expressed-region log-densities
- Collect unspliced-only regions below `density_percentile` (default 0.10)
- Minimum guarantee: at least `min_gdna_regions` (default 100) taken by sorting
- Optional strand filter: symmetric regions (|sense_frac − 0.5| < 0.1)

**π initialization:** `π_init = max(n_gdna,1) / n_eligible`

### Phase 4: EM Iteration (up to 50 iterations)

**E-Step:**
- Hard constraints applied first:
  - `n_spliced > 0` → γ = 0 (definitively expressed)
  - not eligible → γ = π (prior)
- Three independent LLR channels for eligible unspliced regions:

  1. **Density LLR (Gaussian)**: log N(x|μ_G,σ²_G) − log N(x|μ_R,σ²_R)
  2. **Strand LLR (Beta-Binomial or Binomial fallback)**:
     - BetaBin: gDNA ~ BetaBin(k|n,κ/2,κ/2), RNA ~ BetaBin(k|n,κ·SS,κ·(1−SS))
     - When SS = 0.5 → LLR ≡ 0 automatically
     - Binomial fallback when κ estimation fails
  3. **Fragment-Length LLR** (optional): shape-normalized density log-ratio

- Combined posterior: `γ_r = sigmoid(log[π/(1−π)] + LLR_density + LLR_strand + LLR_fl)`

**M-Step:**
- π: fragment-count weighted, clipped to [ε, 1−ε]
- Gaussian density parameters: region-level γ-weighted stats
- κ: marginal likelihood maximization via golden-section search on [0.01, 500]
- FL models: γ-weighted histogram accumulation with ESS threshold

**Convergence:** |Δπ_soft| < 1e-4 (default)

### Phase 5: Final Outputs
- Final E-step with converged parameters
- FL model with min_fl_ess check (default 50); fallback to intergenic FL model
- Returns `GDNACalibration` frozen dataclass

## Fallback Pathways

| # | Trigger | Action | EM Runs? |
|---|---------|--------|----------|
| 1 | **No eligible regions** (`not eligible.any()`) | All γ=1.0, π=1.0, density=0 | No |
| 2 | **Too few eligible** (`n_eligible < min_gdna_regions`, default 100) | Algebraic fallback: spliced→γ=0, unspliced→γ=1, π=mean(γ), κ=2.0 | No |
| 3 | **κ estimation returns None** (< 3 valid strand regions) | E-step uses Binomial LLR fallback | Yes |
| 4 | **Extreme posterior imbalance** (sum(γ)<1 or sum(1−γ)<1) | Reset γ:=0.5, fit symmetric Beta-Binomial | Yes |
| 5 | **FL model ESS too low** (< 50) | Use intergenic FL model fallback | Yes |
| 6 | **Non-convergence** (50 iterations reached) | Return current parameters with `converged=False` | Yes |

## Key Constants

| Constant | Default | Purpose |
|----------|---------|---------|
| `_EPS` | 1e-12 | Numerical floor |
| `GDNA_INIT_DENSITY_PERCENTILE` | 0.10 | eCDF threshold for gDNA seed |
| `GDNA_INIT_MIN_REGIONS` | 100 | Minimum gDNA seed regions / bailout threshold |
| `max_iterations` | 50 | Max EM iterations |
| `convergence_tol` | 1e-4 | π_soft convergence tolerance |
| `min_fl_ess` | 50 | Minimum effective sample size for gDNA FL model |

## GDNACalibration Dataclass Fields

| Field | Type | Description |
|-------|------|-------------|
| `region_posteriors` | ndarray float64 | γ_r per region |
| `gdna_density_global` | float | λ_G = gDNA fragments/bp |
| `kappa_strand` | float | Beta-Binomial concentration |
| `gdna_fl_model` | FragmentLengthModel or None | gDNA fragment-length shape |
| `mixing_proportion` | float | π = P(gDNA) globally |
| `expressed_density` | float | λ_E = RNA fragments/bp |
| `n_iterations` | int | Iterations run |
| `converged` | bool | Convergence flag |
| `region_n_total` | ndarray or None | Per-region fragment counts |
| `region_stats` | dict or None | Full region statistics (diagnostics) |
| `iteration_history` | list or None | Convergence trace (diagnostics) |

## Locus-Level Integration

`compute_gdna_locus_gammas()` aggregates calibration posteriors per-locus:
1. Compute merged genomic intervals per chromosome via `_merged_intervals()`
2. Query cgranges spatial index for overlapping calibration regions
3. Fragment-weighted average: γ_locus = Σ(γ_r × n_r) / Σ(n_r)
4. Fallback: global `mixing_proportion` when no regions overlap
