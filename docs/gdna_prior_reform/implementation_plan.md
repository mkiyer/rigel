# gDNA Prior Reform: Unified OVR Prior + Calibration Warm-Start

## Revision History

- **v1** (2026-03-18): Flat prior `λ_g × L` + warm-start.  Fails for hybrid
  capture / CNV — uniform background assumption is wrong.
- **v2** (2026-03-18): Proportional prior `τ × γ × N`.  Fixes CNV/capture
  robustness but the prior vanishes at low N, exactly where regularization
  is most needed.
- **v3** (2026-03-18): Fixed fractional pseudocount `γ × C` for gDNA only.
  Correct Bayesian scaling but the prior for nRNA (flat α=0.9) remains
  ad-hoc and disconnected from the evidence.
- **v4** (2026-03-18): Unified OVR prior.  Total budget C = 1 virtual
  fragment, partitioned γ×C to gDNA and (1-γ)×C to RNA, with coverage-
  weighted distribution across ALL RNA candidates (mRNA + nRNA).
  Eliminates 5 parameters, replaces 3 separate prior mechanisms with one.

## Motivation

The current EM solver uses three independent, ad-hoc prior mechanisms:

| Component | Mechanism | Source |
|-----------|-----------|--------|
| mRNA | `prior_alpha` (0.01) + OVR coverage-weighted bonus | Data-informed but limited to mRNA |
| nRNA | `nrna_sparsity_alpha` (0.9) | Arbitrary constant, ignores evidence |
| gDNA | `gdna_prior_scale × gdna_init` | Conflated prior + initialization, ~10-1000 pseudo-counts |

These were developed independently and have no unified theoretical
framework.  The nRNA sparsity prior is particularly troublesome: it
ignores evidence from the data, applies the same α=0.9 whether the locus
has strong or zero intronic coverage, and has no principled justification
beyond "we want nRNA to be sparse."

The v4 reform recognizes that the OVR concept — distributing a fixed
prior mass proportional to evidence — can be applied to ALL components:

1. **gDNA prior** = γ × C (from calibration, not coverage)
2. **RNA prior** (mRNA + nRNA) = (1-γ) × C, distributed by coverage weights
3. **Total** = C = 1.0 (one virtual fragment)

The key insight: nRNA is just another RNA isoform.  It represents the
full genomic span from transcript start to end, competing with annotated
mRNA isoforms for fragment assignment.  It should participate in the same
evidence-based prior distribution.

## Theoretical Framework

### The Two Inputs to the C++ EM

In the SQUAREM implementation, the initial state is:

```
state0[i] = theta_init[i] + prior[i]     // VBEM: initial Dirichlet α
state0[i] = normalize(theta_init + prior) // MAP-EM: initial proportions
```

And the M-step is:

```
alpha_new[i] = unambig_totals[i] + em_totals[i] + prior[i]
```

So **`theta_init`** is a one-time starting signal, and **`prior`** is additive
pseudo-counts that persist every M-step iteration.

### A. The Warm-Start: Calibration-Informed γ Proportion

The calibration EM produces per-region posteriors γ_r ∈ [0,1] encoding the
probability that each region is pure gDNA, fusing density + strand + FL
signals across the entire genome.

For each locus, we aggregate these into a single `locus_gamma`:

$$\bar\gamma_{\text{locus}} = \frac{\sum_{r \in R} \gamma_r \cdot n_r}{\sum_{r \in R} n_r}$$

(fragment-weighted to respect hybrid-capture enrichment differences).

In C++, after the coverage-weighted warm-start builds `theta_init` from
unambiguous + coverage-weighted fragment assignments (which assigns ~0 to
gDNA since gDNA has no unambiguous reads), we **override** the gDNA
component to achieve the correct starting proportion:

$$\frac{\theta_{\text{gDNA}}}{\sum_i \theta_i} = \bar\gamma_{\text{locus}}$$

This is **inherently locally adaptive**:
- If the locus has 10× CNV amplification, N is 10× larger → θ_gDNA scales.
- If enriched by hybrid capture (1000×), θ_gDNA scales proportionally.
- If off-target and depleted, θ_gDNA is small but correctly proportioned.

### B. The Prior: Unified Budget-Constrained OVR

The v3 plan proposed changing only gDNA's prior while leaving mRNA (OVR)
and nRNA (sparsity α=0.9) as-is.  The v4 insight is that all three should
share a single, unified prior mechanism.

**The total prior budget is C virtual fragments** (default: 1.0).

The budget is partitioned in two stages:

#### Stage 1: gDNA vs RNA split (calibration-informed)

$$\alpha_{\text{gDNA}} = \bar\gamma_{\text{locus}} \times C$$
$$\text{RNA budget} = (1 - \bar\gamma_{\text{locus}}) \times C$$

The gDNA fraction comes from calibration — genome-level information
fusing density, strand symmetry, and fragment-length signals.

#### Stage 2: RNA distribution (coverage-weighted OVR)

The RNA budget is distributed across ALL RNA candidates (mRNA + nRNA)
proportional to their coverage weights from ambiguous fragment evidence:

$$\alpha_i = \text{RNA budget} \times \frac{\text{coverage}_i}{\sum_{j \in \text{RNA}} \text{coverage}_j}$$

where $\text{coverage}_i$ is the likelihood-weighted coverage total
for RNA component $i$, already computed for the warm-start.

**Why this works:**

1. **One virtual fragment**: The total Σα = C = 1.0.  This is a weak,
   minimally informative prior that lets the data dominate for any locus
   with more than a handful of fragments.

2. **Evidence-proportional**: RNA components with more coverage evidence
   get more prior mass.  This is the same OVR principle that works for
   mRNA — now extended to nRNA.

3. **Natural sparsification**: With C = 1.0 distributed across ~20 RNA
   components, the typical per-component α is ~0.05.  This is sub-unity,
   so the VBEM digamma naturally drives low-evidence components toward
   zero.  The sparsification is now data-driven, not hard-coded.

4. **nRNA gets earned prior mass**: An nRNA span with strong intronic
   evidence (reads in introns) gets substantial OVR prior mass.  An nRNA
   span with no evidence gets ε (effectively zero) → sparsified away.
   Compare to the current flat α=0.9 which gives equal prior mass to all
   nRNA components regardless of evidence.

5. **Budget-constrained**: Unlike the current OVR (which adds an
   unconstrained bonus), the total prior sums to C.  This gives the prior
   a clear interpretation and prevents unexpected scaling with locus size.

### C. Why This Is Better Than Separate Mechanisms

The current system uses three disconnected mechanisms:

| Parameter | Default | Role |
|-----------|---------|------|
| `prior_alpha` | 0.01 | mRNA base prior (minimum floor) |
| `prior_gamma` | 1.0 | mRNA OVR scaling (virtual fragments for coverage bonus) |
| `nrna_sparsity_alpha` | 0.9 | nRNA flat prior (sub-unity for sparsification) |
| `gdna_prior_scale` | 1.0 | gDNA prior multiplier |
| `gdna_kappa_shrink` | None | EB shrinkage strength |

The unified OVR replaces all five with a single parameter:

| Parameter | Default | Role |
|-----------|---------|------|
| C (`prior_pseudocount`) | 1.0 | Total prior budget (virtual fragments) |

The relative influence of the prior scales naturally:

$$\text{Influence} = \frac{C}{N + C} = \frac{1}{N + 1}$$

| N | Influence |
|---|-----------|
| 3 | 25% — strong regularization for sparse loci |
| 20 | 5% — moderate stabilization |
| 100 | 1% — mild |
| 10000 | 0.01% — negligible, EM trusts data |

### D. Quantitative Analysis of gDNA Prior vs v3

The gDNA prior α_gDNA = γ × C is unchanged from v3.  The key numbers:

| Scenario | γ | N | α = γ × 1.0 | Influence |
|----------|---|---|-------------|-----------|
| **Enriched target, contaminated** | 0.3 | 10000 | 0.3 | 0.003% |
| **Off-target, mostly gDNA** | 0.8 | 20 | 0.8 | 3.8% |
| **Sparse locus, contaminated** | 0.5 | 3 | 0.5 | 14% |
| **Standard WTS locus** | 0.1 | 500 | 0.1 | 0.02% |

### E. Quantitative Analysis of nRNA Prior: Unified OVR vs Flat α=0.9

Consider a locus with 2 mRNA isoforms and 1 nRNA span:

| Scenario | Current (α=0.9) | Unified OVR (C=1.0, γ=0.1) |
|----------|------------------|-----------------------------|
| Strong intronic evidence (50% RNA coverage) | 0.9 | ~0.45 |
| Weak intronic evidence (5% RNA coverage) | 0.9 | ~0.045 |
| Zero intronic evidence (0% RNA coverage) | 0.9 | ε ≈ 0 |

The current flat 0.9 gives the same prior regardless of evidence.
The unified OVR adapts: strong nRNA evidence → strong prior, no
evidence → prior collapses.  This is strictly more principled.

The flat α=0.9 was designed to create sparsity.  The unified OVR
achieves sparsity more effectively: zero-coverage components get ε
prior, which under VBEM (digamma(ε) → -∞) ensures instant collapse.

## What Gets Deleted

| Component | Reason |
|-----------|--------|
| `compute_gdna_density_hybrid()` | Noisy per-locus strand algebra |
| `_compute_per_locus_gdna_densities()` | Iteration wrapper for the above |
| `compute_eb_gdna_priors()` | Entire EB shrinkage pipeline |
| `EMConfig.prior_alpha` | Replaced by unified OVR budget |
| `EMConfig.prior_gamma` | Replaced by unified OVR budget |
| `EMConfig.nrna_sparsity_alpha` | Replaced by coverage-weighted OVR |
| `EMConfig.gdna_prior_scale` | Replaced by γ × C |
| `EMConfig.gdna_kappa_shrink` | Shrinkage no longer needed |
| CLI `--prior-alpha` | Config deleted |
| CLI `--prior-gamma` | Config deleted |
| CLI `--nrna-sparsity-alpha` | Config deleted |
| CLI `--gdna-prior-scale` | Config deleted |
| CLI `--gdna-kappa-shrink` | Config deleted |
| `estimator.transcript_unspliced_sense/antisense` | Only consumed by deleted code |
| C++ sense/antisense accumulation in scanner | Only feeds deleted code |
| Tripartite `base_prior` construction in `batch_locus_em()` | Prior computed inside unified OVR |

(Note: `transcript_intronic_sense/antisense` were dead code and have been removed.)

## What Changes

### 1. `GDNACalibration` — Store Per-Region Fragment Counts (calibration.py)

Add field:

```python
region_n_total: np.ndarray  # (n_regions,) float64
```

Populated from `stats["n_total"]` inside `calibrate_gdna()`.  Cost:
~8 KB for 1000 regions.  Required for fragment-weighted γ aggregation.

### 2. New Function: `compute_gdna_locus_gammas()` (locus.py)

```python
def compute_gdna_locus_gammas(
    loci: list[Locus],
    index: TranscriptIndex,
    calibration: GDNACalibration,
) -> np.ndarray:
    """Aggregate calibration region posteriors per locus.

    For each locus, finds overlapping calibration regions via
    ``index.region_cr`` (cgranges) and computes a fragment-weighted
    aggregate gDNA posterior:

        γ_locus = Σ(γ_r × n_r) / Σ(n_r)

    Fallback: ``calibration.mixing_proportion`` (global π) when no
    regions overlap a locus.

    Returns
    -------
    locus_gammas : np.ndarray, shape (n_loci,), float64 ∈ [0, 1]
    """
```

Implementation: one cgranges overlap query per locus (~5000 queries,
O(log N) each, negligible).

### 3. C++ EM Solver — Unified OVR (em_solver.cpp)

The current `compute_ovr_prior_and_warm_start()` is restructured.
The function signature changes:

```cpp
// OLD signature:
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    const double* base_prior,   // pre-built per-component prior
    int           n_t,          // number of mRNA (OVR scope)
    double        gamma,        // OVR scaling factor
    double*       prior_out,
    double*       theta_init_out,
    int           n_components);

// NEW signature:
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        locus_gamma,       // calibration gDNA fraction ∈ [0,1]
    double        total_pseudocount, // C (total prior budget, default 1.0)
    int           gdna_idx,          // index of gDNA component (-1 if none)
    double*       prior_out,
    double*       theta_init_out,
    int           n_components);
```

Parameters removed: `base_prior`, `n_t`, `gamma`.
Parameters added: `locus_gamma`, `total_pseudocount`, `gdna_idx`.

The implementation:

```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        locus_gamma,
    double        total_pseudocount,
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components)
{
    // Initialize theta_init from unambig_totals
    std::copy(unambig_totals, unambig_totals + n_components, theta_init_out);

    // Accumulate coverage totals from all equivalence classes
    std::vector<double> coverage_totals(n_components, 0.0);

    for (const auto& ec : ec_data) {
        const int n = ec.n;
        const int k = ec.k;
        const int32_t* cidx = ec.comp_idx.data();
        const double* wt = ec.wt_flat.data();

        for (int i = 0; i < n; ++i) {
            double row_sum = 0.0;
            for (int j = 0; j < k; ++j) {
                double w = wt[i * k + j] * eligible[cidx[j]];
                row_sum += w;
            }
            if (row_sum == 0.0) row_sum = 1.0;
            double inv_row_sum = 1.0 / row_sum;

            for (int j = 0; j < k; ++j) {
                double share = (wt[i * k + j] * eligible[cidx[j]])
                             * inv_row_sum;
                theta_init_out[cidx[j]] += share;
                coverage_totals[cidx[j]] += share;
            }
        }
    }

    // --- Unified OVR prior: budget-constrained distribution ---

    // Compute total RNA coverage (excluding gDNA and ineligible)
    double total_rna_coverage = 0.0;
    int n_rna_eligible = 0;
    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] > 0.0 && i != gdna_idx) {
            total_rna_coverage += coverage_totals[i];
            ++n_rna_eligible;
        }
    }

    double rna_budget = (1.0 - locus_gamma) * total_pseudocount;
    double gdna_alpha = std::max(locus_gamma * total_pseudocount,
                                 EM_LOG_EPSILON);

    for (int i = 0; i < n_components; ++i) {
        if (eligible[i] <= 0.0) {
            prior_out[i] = 0.0;
        } else if (i == gdna_idx) {
            prior_out[i] = gdna_alpha;
        } else if (total_rna_coverage > 0.0) {
            prior_out[i] = std::max(
                rna_budget * coverage_totals[i] / total_rna_coverage,
                EM_LOG_EPSILON);
        } else if (n_rna_eligible > 0) {
            // No coverage data — distribute uniformly
            prior_out[i] = rna_budget / n_rna_eligible;
        } else {
            prior_out[i] = EM_LOG_EPSILON;
        }
    }

    // --- gDNA warm-start override ---
    if (gdna_idx >= 0 && gdna_idx < n_components
        && eligible[gdna_idx] > 0.0)
    {
        double total_theta = 0.0;
        for (int i = 0; i < n_components; ++i)
            total_theta += theta_init_out[i];

        double others = total_theta - theta_init_out[gdna_idx];
        if (others > 0.0 && locus_gamma < 1.0 - 1e-10) {
            theta_init_out[gdna_idx] =
                (locus_gamma / std::max(1.0 - locus_gamma, 1e-10))
                * others;
        } else {
            theta_init_out[gdna_idx] =
                std::max(locus_gamma * total_theta, EM_LOG_EPSILON);
        }
    }
}
```

**Why this is correct:**
- theta_init for RNA components is built from unambig_totals + coverage
  shares (unchanged from current code).
- theta_init for gDNA is overridden to achieve θ_gDNA / Σθ = γ_locus.
- The prior for ALL components comes from a single budget C = 1.0:
  gDNA gets γ × C, RNA gets (1-γ) × C distributed by coverage.
- The `base_prior` array and tripartite construction are eliminated.

### 4. `batch_locus_em()` Per-Locus Loop (em_solver.cpp)

The tripartite base_prior setup is eliminated.  The per-locus code
simplifies from:

```cpp
// OLD: tripartite base_prior + OVR
int gdna_idx = sub.gdna_idx;
std::vector<double> base_prior(nc, prior_alpha);
for (int i = sub.n_t; i < gdna_idx; ++i)
    base_prior[i] = nrna_sparsity_alpha;
if (gdna_idx >= 0 && gdna_idx < nc)
    base_prior[gdna_idx] = gdna_prior_scale * gdna_inits.data()[li];

std::vector<double> prior(nc);
std::vector<double> theta_init(nc);
compute_ovr_prior_and_warm_start(
    ec_data, sub.unambig_totals.data(), sub.eligible.data(),
    base_prior.data(), sub.n_t, prior_gamma,
    prior.data(), theta_init.data(), nc);
```

To:

```cpp
// NEW: unified OVR prior + warm-start
std::vector<double> prior(nc);
std::vector<double> theta_init(nc);
compute_ovr_prior_and_warm_start(
    ec_data, sub.unambig_totals.data(), sub.eligible.data(),
    locus_gammas.data()[li], total_pseudocount, sub.gdna_idx,
    prior.data(), theta_init.data(), nc);
```

### 5. Nanobind Interface (em_solver.cpp binding)

```cpp
// Parameters REMOVED from batch_locus_em:
f64_1d gdna_inits,           // [n_loci]
double prior_alpha,
double prior_gamma,
double nrna_sparsity_alpha,
double gdna_prior_scale,

// Parameters ADDED:
f64_1d locus_gammas,          // [n_loci] — calibration gDNA fraction ∈ [0,1]
double total_pseudocount,     // C — total prior budget (default 1.0)
```

Net: 5 parameters removed, 2 added.

### 6. Estimator (estimator.py)

```python
def run_batch_locus_em(
    self, loci, em_data, index,
    locus_gammas: np.ndarray,   # replaces gdna_inits
    ...
)
```

Delete references to `prior_alpha`, `prior_gamma`, `nrna_sparsity_alpha`,
`gdna_prior_scale`.  Pass `total_pseudocount` from config.
Delete `transcript_unspliced_sense/antisense` fields.

### 7. Pipeline (pipeline.py)

`_compute_priors()` is replaced by a simpler call:

```python
from .locus import compute_gdna_locus_gammas

locus_gammas = compute_gdna_locus_gammas(loci, index, calibration)
estimator.run_batch_locus_em(loci, em_data, index, locus_gammas, ...)
```

### 8. Config / CLI Cleanup

Delete from `EMConfig`:
- `prior_alpha` (0.01)
- `prior_gamma` (1.0)
- `nrna_sparsity_alpha` (0.9)
- `gdna_prior_scale` (1.0)
- `gdna_kappa_shrink` (None)

Add to `EMConfig`:
- `prior_pseudocount: float = 1.0` — total OVR budget C

Delete from CLI (5 arguments):
- `--prior-alpha`
- `--prior-gamma`
- `--nrna-sparsity-alpha`
- `--gdna-prior-scale`
- `--gdna-kappa-shrink`

Add to CLI (1 argument, optional):
- `--prior-pseudocount` — total OVR budget

Net: 5 config fields deleted, 1 added.  5 CLI args deleted, 1 added.

### 9. C++ Scanner — Delete Unspliced Sense/Antisense Accumulation

The `transcript_unspliced_sense/antisense` arrays are only consumed by
deleted code.  Remove their accumulation in `scan.py` and their
declaration in `estimator.py`.

## Implementation Order

1. Add `region_n_total` to `GDNACalibration` + populate in `calibrate_gdna()`
2. Add `compute_gdna_locus_gammas()` to `locus.py`
3. Restructure C++ `compute_ovr_prior_and_warm_start()` — new signature
   with `locus_gamma`, `total_pseudocount`, `gdna_idx`
4. Simplify `batch_locus_em()` — remove tripartite base_prior, remove
   5 old parameters, add `locus_gammas` + `total_pseudocount`
5. Update nanobind binding
6. Update `estimator.py` to pass `locus_gammas` + `total_pseudocount`
7. Update `pipeline.py` to call `compute_gdna_locus_gammas()`
8. Delete old code: `compute_eb_gdna_priors`, `compute_gdna_density_hybrid`,
   `_compute_per_locus_gdna_densities`, strand accumulators
9. Delete config: 5 fields from EMConfig + CLI
10. Add: `prior_pseudocount` to EMConfig + CLI
11. Update tests, regenerate goldens
12. Build, full test suite, lint

## Data Flow

```
calibration.py              locus.py                    C++ EM solver
┌──────────────┐     ┌─────────────────────┐     ┌──────────────────────────┐
│              │     │                     │     │ compute_ovr_prior_and_   │
│ γ_r (region  │────>│ compute_gdna_       │     │  warm_start():           │
│  posteriors) │     │  locus_gammas()     │     │                          │
│              │     │                     │     │ 1. theta_init from       │
│ n_r (region  │────>│  γ_locus[i] =      │────>│    unambig + coverage    │
│  frag counts)│     │   Σ γ_r n_r / Σ n_r│     │                          │
│              │     │                     │     │ 2. coverage_totals for   │
│ region_df    │     └─────────────────────┘     │    all components        │
│  (intervals) │            ↑                    │                          │
└──────────────┘      index.region_cr             │ 3. prior:               │
                      (cgranges overlap)          │    gDNA = γ × C         │
                                                  │    RNA_i = (1-γ)×C ×    │
                                                  │      cov_i / Σ cov_RNA  │
                                                  │                          │
                                                  │ 4. θ[gDNA] warm-start:  │
                                                  │    γ/(1-γ) × others     │
                                                  │                          │
                                                  │ 5. Run SQUAREM           │
                                                  └──────────────────────────┘
```

## Edge Cases

| Scenario | γ_locus | N | gDNA prior | RNA prior | Behavior |
|----------|---------|---|------------|-----------|----------|
| No overlapping regions | = π (global) | any | π × C | (1-π) × C by coverage | Conservative fallback |
| Mega-locus (500kb) | 0.3 | 50000 | 0.3 | 0.7 distributed | Prior 0.002% of data — EM drives |
| Clean RNA | 0.01 | 100 | 0.01 | 0.99 distributed | gDNA prior negligible → collapses |
| Sparse off-target | 0.8 | 5 | 0.8 | 0.2 distributed | Prior 17% of data — stabilizes |
| Enriched target | 0.3 | 10000 | 0.3 | 0.7 distributed | Prior 0.01% — EM uses likelihoods |
| Cancer CNV 10× | 0.3 | 5000 | 0.3 | 0.7 distributed | Warm-start scales with N |
| Pure gDNA locus | 0.99 | 200 | 0.99 | 0.01 distributed | gDNA dominates, converges fast |
| No ambig fragments | any | any | γ × C | uniform | Prior irrelevant — all unambig |
| nRNA with evidence | 0.1 | 500 | 0.1 | nRNA gets coverage share | nRNA prior earned from data |
| nRNA without evidence | 0.1 | 500 | 0.1 | nRNA gets ε → 0 | Natural VBEM sparsification |

## Risks & Mitigations

1. **region_cr not available** (test fixtures without regions):
   Fall back to `γ_locus = calibration.mixing_proportion` (global π).

2. **Warm-start inaccuracy**: If calibration γ is wrong for a specific
   locus, the EM self-corrects — fragment-level strand/FL likelihoods
   override the starting point within a few iterations.  The weak prior
   (C = 1.0) does not resist this correction.

3. **Zero-coverage RNA components**: These get prior ε ≈ 0.  Under
   VBEM, digamma(ε) → -∞, ensuring instant collapse.  This is
   BETTER than the current flat α=0.9 which keeps dead components alive.

4. **All RNA coverage zero** (no ambiguous fragments): Fall back to
   uniform distribution of RNA budget.  But since all fragments are
   unambiguous, the EM solution is trivially the unambiguous counts
   and the prior is irrelevant.

5. **Removal of prior_alpha floor**: Currently mRNA components get
   at least 0.01 base_prior even with zero coverage.  Under the
   unified OVR, zero-coverage mRNA gets ε.  This is correct: if a
   transcript has no likelihood for any fragment in the locus, it
   should not receive prior mass.  (The eligible[] mask already
   handles truly absent components.)

6. **Circular reasoning (empirical Bayes)**: The prior uses γ (from
   calibration) and the EM processes fragments from the same data.
   This is standard empirical Bayes: calibration operates on region-level
   aggregate summary statistics, while the locus EM evaluates individual
   fragment likelihoods.  The information channels are distinct.

7. **C = 1.0 sensitivity**: One virtual fragment is the minimally
   informative proper prior.  If tuning is needed, the single
   `--prior-pseudocount` parameter can be adjusted.  The Bayesian
   self-scaling property means the exact value matters little for N >> C.

## Tests

### Unit: `compute_gdna_locus_gammas()`
- High-γ regions → γ_locus ≈ 1.0
- Low-γ regions → γ_locus ≈ 0.0
- Mixed regions, fragment-weighted → dominated by high-count region
- No overlapping regions → falls back to π
- Single-region locus → γ_locus = γ_r exactly
- γ_locus always ∈ [0, 1]

### Unit: C++ unified OVR prior
- gDNA prior = γ × C
- RNA prior sums to (1-γ) × C (within float tolerance)
- RNA prior proportional to coverage_totals
- Zero-coverage RNA component gets ε
- mRNA/nRNA both participate in coverage distribution
- gDNA warm-start: θ_gDNA / Σθ ≈ γ_locus
- With no ambiguous fragments: uniform RNA distribution

### Unit: VBEM sparsification
- Zero-coverage nRNA component collapses to ε after 1 iteration
- Strong-coverage nRNA component persists

### Integration
- `GDNACalibration.region_n_total` populated with correct shape
- End-to-end pipeline produces valid, non-negative counts
- gDNA counts reflect calibration predictions (high γ → high gDNA)

### Regression
- Golden output regeneration
- Scenario test tolerances (may need widening)
