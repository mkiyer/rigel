# Mode-Aware Dirichlet Prior for EM Quantification

**Date:** 2026-04-07
**Status:** Proposed — awaiting review before implementation
**Supersedes:** `anchored_empirical_bayes_plan.md`, `anchored_empirical_bayes_implementation.md`

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Root Cause: VBEM Digamma Sparsification](#2-root-cause-vbem-digamma-sparsification)
3. [Theory: Jeffreys-Corrected Prior](#3-theory-jeffreys-corrected-prior)
4. [MAP-EM: No Correction Needed](#4-map-em-no-correction-needed)
5. [Unified Design: Mode-Aware Baseline](#5-unified-design-mode-aware-baseline)
6. [Simulation Results](#6-simulation-results)
7. [Design Specification](#7-design-specification)
8. [Cleanup: What to Remove](#8-cleanup-what-to-remove)
9. [Implementation Plan](#9-implementation-plan)
10. [Testing Strategy](#10-testing-strategy)
11. [Risk Analysis](#11-risk-analysis)

---

## 1. Problem Statement

The V3 calibration system computes a per-locus gDNA mixing fraction γ from
strand decomposition and density baselines. This γ must be translated into a
Dirichlet prior `(α_gDNA, α_RNA)` for the VBEM solver. The current
implementation uses:

```
C = κ × N_locus
α_gDNA = γ × C
α_RNA  = (1 − γ) × C
```

where κ interpolates between κ_min (0.001) and κ_max (0.50) based on strand
specificity via `w = (2·SS − 1)²`.

### Failures

This design has two coupled failure modes:

**A. α_RNA amplification (siphoning).** At SS=0.65 with N=2000 fragments:
w = 0.09, κ = 0.455, C = 910. The RNA pseudocount α_RNA ≈ 819 is distributed
proportionally across mRNA and nRNA components by coverage. This massive prior
dominates the EM for small-to-medium loci, causing nRNA components to siphon
fragments from mRNA.

**B. gDNA collapse.** When C is made small (e.g., fixed C=2), α_gDNA = 0.20
for γ=0.10. With only 0.20 pseudocounts, VBEM systematically drains the gDNA
component to zero over ~150 iterations. The gDNA estimate collapses even though
the warm-start initialized it correctly.

These failures create a **Goldilocks dilemma**: C must be large enough to
prevent gDNA collapse but small enough to avoid RNA siphoning. The κ×N scaling
makes this worse, not better, because it ties these two requirements together.

### Failed approaches

| Approach | Problem |
|----------|---------|
| Fixed C=2 | gDNA collapses for N ≤ 2000 (α_gDNA < 0.5) |
| Fixed C=10 | Prevents collapse but causes siphoning (α_RNA = 9.0) |
| Fixed C=100 | Heavy siphoning, warps all component ratios |
| κ×N (current) | Couples C to N; simultaneously too large (siphon) and too small (collapse) |
| α floor (max(0.5, γ·C)) | Distorts the prior ratio; γ=0.001 gets 20% gDNA prior |

---

## 2. Root Cause: VBEM Digamma Sparsification

The root cause is a **systematic bias in the VBEM E-step** arising from the
digamma function's asymptotic behavior.

### Derivation

The VBEM E-step computes responsibilities:

$$r_j \propto \ell_j \cdot \exp\!\bigl(\psi(\theta_j + \alpha_j)\bigr)$$

where ψ is the digamma function, θ_j is the current count for component j, and
α_j is the Dirichlet pseudocount.

For the **degenerate case** (likelihoods ℓ_j identical across components, which
occurs when gDNA reads have no distinguishing features — e.g., at SS=0.50):

$$r_j = \frac{\exp(\psi(\theta_j + \alpha_j))}{\sum_k \exp(\psi(\theta_k + \alpha_k))}$$

Using the digamma asymptotic expansion:

$$\psi(x) \approx \ln(x) - \frac{1}{2x} - \frac{1}{12x^2} + \cdots$$

we get:

$$\exp(\psi(x)) \approx x \cdot \exp\!\left(-\frac{1}{2x}\right) \approx x - \frac{1}{2} + O(x^{-1})$$

Therefore the VBEM update effectively computes:

$$\theta_j^{(\text{new})} = N \cdot \frac{\theta_j + \alpha_j - \tfrac{1}{2}}{\sum_k (\theta_k + \alpha_k - \tfrac{1}{2})} = N \cdot \frac{\theta_j + \alpha_j - \tfrac{1}{2}}{N + C - K/2}$$

where K is the number of components and C = Σ α_k.

### Per-iteration drift

The per-iteration change in θ_j is:

$$\Delta\theta_j = \theta_j^{(\text{new})} - \theta_j = \frac{N \cdot (\alpha_j - \tfrac{1}{2}) - \theta_j \cdot (C - K/2)}{N + C - K/2}$$

At equilibrium (Δθ_j = 0):

$$\theta_j^* = \frac{N \cdot (\alpha_j - \frac{1}{2})}{C - K/2}$$

**Key insight:** When α_j < 0.5, the effective count (α_j − 0.5) is **negative**,
driving θ_j toward zero. This is the sparsification bias.

### Numerical validation

For a 2-component system with θ_gDNA=200, θ_RNA=1800 (N=2000):

| α_gDNA | drift/iter | iterations to collapse |
|:------:|:----------:|:---------------------:|
| 0.01   | −1.38      | 145                   |
| 0.10   | −1.29      | 155                   |
| 0.50   | −0.90      | 223                   |
| 1.00   | −0.40      | 502                   |
| 2.00   | +0.60      | (grows instead)       |

Even α = 0.5 shows downward drift. The neutral point is at α ≈ 1.5 for this
specific configuration (it depends on the other components' α values). This is
why no simple floor or scaling law fixes the problem.

---

## 3. Theory: Jeffreys-Corrected Prior

### The correction

The Jeffreys prior for a K-component Multinomial is Dir(½, ½, …, ½). This is
the well-established non-informative reference prior that naturally compensates
for the digamma bias.

We combine the Jeffreys baseline with calibration information:

$$\alpha_j = \frac{1}{2} + \gamma_j \cdot C_{\text{base}}$$

where:
- **½** = Jeffreys baseline (compensates the digamma −½ bias)
- **γ_j** = calibrated mixing fraction for component j (from V3 calibration)
- **C_base** = calibration evidence parameter (a fixed constant, e.g. 5.0)

Total concentration: C = K/2 + C_base.

### Proof of correctness (degenerate case)

Substituting into the VBEM equilibrium:

$$\theta_j^* = \frac{N \cdot (\alpha_j - \frac{1}{2})}{C - K/2} = \frac{N \cdot \gamma_j \cdot C_{\text{base}}}{C_{\text{base}}} = N \cdot \gamma_j$$

The C_base **cancels exactly**. The equilibrium is θ_j/N = γ_j regardless of:
- The value of C_base
- The number of fragments N
- The number of components K

This is why Part 3 of our simulation shows zero sensitivity to C_base.

### Intuition

The standard Dirichlet prior Dir(γ₁·C, γ₂·C, …) works correctly for rigel's
MAP-EM, which uses additive pseudocounts and has no algorithmic bias (see §4).
VBEM, however, uses the **mean of the log-posterior** via the digamma function,
which shifts the effective count by −½. The Jeffreys baseline pre-compensates
for this shift, restoring the correct equilibrium.

### Interaction with warm-start

The Jeffreys correction fixes the **equilibrium** but not the **basin of
attraction**. For large N (≥ 100k) with degenerate likelihoods, the VBEM
has two stable equilibria: one near γ and one near 0.5 (uniform). The
warm-start θ_gDNA = γ × N ensures convergence to the correct equilibrium.

For small N (≤ 2000), the basins merge and the correction alone suffices —
warm-start merely accelerates convergence.

The existing warm-start logic in `compute_ovr_prior_and_warm_start()` already
implements γ-proportional initialization via `θ[gdna] = (α_gdna/α_rna) × others`,
which with the Jeffreys prior gives:

$$\theta_{gDNA}^{(0)} = \frac{0.5 + \gamma \cdot C_{\text{base}}}{0.5 + (1-\gamma) \cdot C_{\text{base}}} \cdot \theta_{RNA} \approx \frac{\gamma}{1-\gamma} \cdot \theta_{RNA}$$

This is the correct initialization for any C_base > 0.

### Interaction with strand likelihoods

When strand specificity provides per-fragment log-likelihood ratios between
gDNA and RNA (LLR < 0 for gDNA), the **likelihood dominates** the prior for
moderate-to-high SS. At SS=0.65, the average per-fragment LLR is −0.046 nats,
giving a total LLR of −91 nats for N=2000. This overwhelms C_base=5.0,
correctly suppressing gDNA when the strand signal says so.

Simulation confirms: at SS≥0.55 with N≥500, gDNA fraction drops below 0.001
regardless of prior (Jeffreys or standard), because the likelihood dominates.

---

## 4. MAP-EM: No Correction Needed

### Rigel's MAP-EM is not the textbook Dirichlet mode formula

A natural question: does the Jeffreys +0.5 baseline also apply to MAP-EM mode?
The answer is **no** — and applying it would introduce bias.

The textbook MAP estimate for a Dirichlet-Multinomial uses the posterior
**mode**:

$$\theta_j^{MAP} = \frac{N_j + \alpha_j - 1}{\sum_k (N_k + \alpha_k - 1)}$$

If rigel used this formula, a +1.0 (Laplace) baseline would be needed to cancel
the −1, analogous to how +0.5 cancels the digamma −0.5 in VBEM.

**However, rigel's MAP-EM does not use the mode formula.** The actual M-step
in `em_solver.cpp` (line 704–711) is:

```cpp
// Plain MAP-EM: theta_new = (unambig + em + prior)
for (int i = 0; i < n_components; ++i) {
    theta_new[i] = unambig_totals[i] + em_totals[i] + prior[i];
    total += theta_new[i];
}
// ... normalize by total ...
```

This computes the posterior **mean** (additive pseudocounts), not the mode:

$$\theta_j^{\text{new}} = \frac{N_j + \text{prior}[j]}{\sum_k (N_k + \text{prior}[k])} = \frac{N_j + \text{prior}[j]}{N + C}$$

Compare with the VBEM M-step (line 766–769):
```cpp
// M-step: alpha_new = unambig_totals + em_totals + prior (unnormalized)
alpha_new[i] = unambig_totals[i] + em_totals[i] + prior[i];
```

Both modes add `prior[j]` identically. The **only difference** is the
E-step weight computation:

| Mode | E-step log-weights | Bias |
|------|-------------------|------|
| MAP  | `log(θ_j + ε)`   | None |
| VBEM | `ψ(α_j) − ψ(Σα)` | −0.5 per component (digamma) |

### MAP-EM equilibrium (degenerate likelihoods)

For degenerate likelihoods (all components identical), the MAP-EM E-step gives:

$$r_j = \frac{\theta_j}{\sum_k \theta_k}$$

And the M-step:

$$\theta_j^{\text{new}} = \frac{N \cdot \theta_j / \sum\theta + \text{prior}[j]}{N + C}$$

At equilibrium ($\theta_j = \theta_j^{\text{new}}$ up to normalization):

$$\theta_j \cdot C = \text{prior}[j]$$

$$\theta_j = \frac{\text{prior}[j]}{C} = \gamma_j$$

This is **exact** with no baseline correction, for any N and any C_base > 0.
There is no digamma, no −0.5, no mode shift. The prior acts as pure additive
pseudocounts and the equilibrium is simply the prior's own mean.

Furthermore, with a warm-start at θ = [γ, 1−γ], MAP-EM converges in **one
iteration** for degenerate likelihoods (since the E-step returns the same
ratios as the initial θ, and the M-step preserves them).

### Why Laplace (+1.0) is wrong for rigel's MAP-EM

If we add a +1.0 Laplace baseline, the prior becomes:

$$\text{prior}[j] = 1.0 + \gamma_j \cdot C_{\text{base}}$$

The implied equilibrium is:

$$\theta_j = \frac{1.0 + \gamma_j \cdot C_{\text{base}}}{K + C_{\text{base}}}$$

For γ=0.10, C_base=5.0, K=9: θ_gDNA = (1.0 + 0.5) / (9 + 5) = 0.107.
This is upward-biased. At small N the bias is worse because the +1.0 per
component dominates: for K=9, the +1.0 baseline adds 9.0 symmetric
pseudocounts, pushing everything toward uniform.

Simulation confirms: with Laplace +1.0 at N=100, K=9, θ_gDNA = 0.107 instead
of the correct 0.100, and the mRNA/nRNA split distorts from 0.585/0.315 to
0.495/0.398. At N=2000 the bias is similar (0.107 vs 0.100).

### Simulation validation

From `scripts/debug/analyze_map_vs_vbem_baseline.py`:

MAP-EM with no baseline, γ=0.10, 9 components:

| N | C_base | gDNA | mRNA | nRNA |
|--:|:------:|:----:|:----:|:----:|
| 100 | 1 | 0.1000 | 0.5850 | 0.3150 |
| 100 | 5 | 0.1000 | 0.5850 | 0.3150 |
| 2,000 | 1 | 0.1000 | 0.5850 | 0.3150 |
| 2,000 | 5 | 0.1000 | 0.5850 | 0.3150 |
| 100,000 | 5 | 0.1000 | 0.5850 | 0.3150 |

**Perfect convergence** at every N and C_base, with no baseline. The additive
pseudocount formulation needs no correction.

---

## 5. Unified Design: Mode-Aware Baseline

The two EM modes have different algebraic biases requiring different baselines:

| EM Mode | E-step bias | Required baseline | Prior formula |
|---------|:-----------:|:-----------------:|:-------------|
| VBEM    | −0.5 (digamma) | +0.5 (Jeffreys) | `prior[j] = 0.5 + γ_j × C_base` |
| MAP     | 0 (log)        | 0 (none)         | `prior[j] = γ_j × C_base` |

The unified rule:

$$\text{prior}[j] = b + \gamma_j \cdot C_{\text{base}}$$

where $b$ is the **mode-aware baseline**:

$$b = \begin{cases} 0.5 & \text{if mode = VBEM (Jeffreys)} \\ 0.0 & \text{if mode = MAP} \end{cases}$$

### Why this is the right abstraction

Both baselines serve the same purpose: **cancel the inherent bias of the
algorithm so that the degenerate equilibrium is exactly γ_j**. VBEM's digamma
function subtracts ≈0.5 per component; the Jeffreys +0.5 compensates. MAP's
log function has no such shift; no compensation needed.

The rest of the design — C_base, warm-start, gDNA gate, global strand
calibration — is identical for both modes.

### Implementation: where the baseline lives

The baseline is applied in the C++ `compute_ovr_prior_and_warm_start()`
function, which already receives a `use_vbem` flag. This is the natural place
to select the baseline:

```cpp
const double baseline = use_vbem ? 0.5 : 0.0;

for (int i = 0; i < n_components; ++i) {
    if (eligible[i] > 0.0) {
        prior_out[i] = baseline;
        if (i == gdna_idx) {
            prior_out[i] += alpha_gdna;
        } else if (total_rna_coverage > 0.0) {
            prior_out[i] += alpha_rna * coverage[i] / total_rna_coverage;
        }
    }
}
```

Python passes `alpha_gdna = γ × c_base` and `alpha_rna = (1−γ) × c_base`,
independent of EM mode. The C++ side adds the mode-aware baseline.

---

## 6. Simulation Results

All results from `scripts/debug/analyze_jeffreys_correction.py` (VBEM) and
`scripts/debug/analyze_map_vs_vbem_baseline.py` (MAP vs VBEM comparison).

### 6.1 Two-component degenerate (VBEM)

γ = 0.10, warm-start at γ, C_base = 5.0:

| N | Jeffreys θ_gDNA | Standard θ_gDNA (C=6.0) |
|--:|:---------------:|:-----------------------:|
| 50 | 0.1013 | 0.0255 |
| 100 | 0.1007 | 0.0233 |
| 500 | 0.1000 | 0.0208 |
| 2,000 | 0.1000 | 0.0207 |
| 10,000 | 0.1000 | 0.0494 |
| 100,000 | 0.1000 | 0.0924 |
| 1,000,000 | 0.1000 | 0.1000 |

Jeffreys converges to γ=0.10 for **every** N. Standard prior only converges
at N=1M.

### 6.2 Nine-component (VBEM)

γ = 0.10, warm-start: 65% mRNA / 35% nRNA within RNA:

| N | C_base | gDNA | mRNA | nRNA |
|--:|:------:|:----:|:----:|:----:|
| 100 | 1 | 0.1005 | 0.5761 | 0.3234 |
| 100 | 10 | 0.1000 | 0.5841 | 0.3158 |
| 2,000 | 1 | 0.1000 | 0.5850 | 0.3150 |
| 2,000 | 10 | 0.1000 | 0.5850 | 0.3150 |
| 100,000 | 1 | 0.1000 | 0.5850 | 0.3150 |
| 100,000 | 10 | 0.1000 | 0.5850 | 0.3150 |

All configurations converge to the correct fractions. No siphoning.

### 6.3 C_base sensitivity (VBEM)

Nine-component, N=2000, γ=0.10:

| C_base | gDNA | mRNA | nRNA |
|:------:|:----:|:----:|:----:|
| 0.5 | 0.10000 | 0.58500 | 0.31500 |
| 1.0 | 0.10000 | 0.58500 | 0.31500 |
| 5.0 | 0.10000 | 0.58500 | 0.31500 |
| 20.0 | 0.10000 | 0.58500 | 0.31500 |
| 50.0 | 0.10000 | 0.58500 | 0.31500 |

**Zero sensitivity.** C_base is no longer a fragile hyperparameter.

### 6.4 Small γ behavior (VBEM)

C_base = 5.0, N = 2000:

| γ | Jeffreys θ_gDNA | Standard θ_gDNA |
|:---:|:-------:|:-------:|
| 0.001 | 0.0026 | 0.0000 |
| 0.005 | 0.0057 | 0.0000 |
| 0.01 | 0.0104 | 0.0000 |
| 0.05 | 0.0500 | 0.0000 |
| 0.10 | 0.1000 | 0.0207 |

The Jeffreys correction prevents collapse at all γ values. At very small γ
(0.001) with small N, there is slight **upward** bias (~0.15% absolute) because
the ½ baseline is large relative to γ·C_base. This is cosmetic: the absolute
error is ~3 fragments out of 2000, and it vanishes at larger N.

---

## 7. Design Specification

### 7.1 Prior construction

For a locus with K components (1 gDNA + K−1 RNA), calibrated γ, and
per-component coverage fractions `cov_j` (summing to 1 over RNA components):

```
b = 0.5 if mode == VBEM else 0.0    (mode-aware baseline)

α_gDNA  = b + γ × C_base
α_RNA_j = b + (1 − γ) × cov_j × C_base    for each RNA component j
```

Total: C = K × b + C_base

**C_base** is a single fixed constant. Default: **5.0** (value is insensitive;
anything in [1, 50] gives identical results in the degenerate case).

The baseline `b` is selected by the EM mode and applied in C++. Python passes
only the calibration portion: `alpha_gdna = γ × c_base`, `alpha_rna = (1−γ) × c_base`.

### 7.2 gDNA gate

When calibration returns γ = 0 (pure RNA), the gDNA component is **disabled**:
- α_gDNA = 0.0 (overriding any baseline)
- θ_gDNA^(0) = 0.0

This preserves the existing gate behavior. The baseline is only meaningful
when gDNA is eligible. Applies identically to both VBEM and MAP modes.

### 7.3 Warm-start

No change to the warm-start mechanism. The existing logic in
`compute_ovr_prior_and_warm_start()` uses `θ[gdna] = (α_gdna/α_rna) × others`,
which correctly initializes at:

$$\frac{\theta_{gDNA}}{\theta_{total}} = \frac{\alpha_{gDNA}}{\alpha_{gDNA} + \alpha_{RNA}} = \frac{0.5 + \gamma \cdot C_{\text{base}}}{K/2 + C_{\text{base}}} \approx \gamma$$

(exact when K=2 or when C_base >> K/2)

### 7.4 Global strand calibration (rectifier bias fix)

**Separate from the Jeffreys correction**, but needed for the 4 failing tests.

The current per-region strand formula applies `max(0, ...)` to each region
independently, then sums. Because the rectifier max(0, x) has E[max(0, x)] > 0
even when E[x] = 0 (for symmetric noise), this creates a systematic upward
bias in the gDNA estimate for pure RNA samples.

**Fix:** Aggregate strand counts globally (sum numerator and denominator across
regions) **before** applying max(0, ...). This eliminates the rectifier bias
because the aggregated estimate has much lower variance.

```python
# BEFORE (per-region, then sum — biased):
e_gdna = max(0, (n_anti_1 - n_unspliced_1·(1-SS)) / (SS-0.5))
       + max(0, (n_anti_2 - n_unspliced_2·(1-SS)) / (SS-0.5))
       + ...

# AFTER (sum, then rectify — unbiased):
total_anti = sum(n_anti_i)
total_uns  = sum(n_unspliced_i)
e_gdna_global = max(0, (total_anti - total_uns·(1-SS)) / (SS-0.5))
```

The per-region γ values (needed for `compute_locus_priors` to look up
per-locus γ) are then computed by distributing the global estimate
proportionally: `e_gdna_i = e_gdna_global × (n_unspliced_i / total_uns)`.

Note: the user raised that a global aggregator does not capture region-to-region
overdispersion (beta-binomial, not binomial). This is a known limitation.
However, for the purpose of the prior (which is now insensitive to C_base),
the precise local γ matters less — what matters is that γ > 0 vs γ = 0 is
correct, which the global aggregator handles well.

### 7.5 Configuration changes

**Remove:**
- `CalibrationConfig.gdna_prior_kappa_min`
- `CalibrationConfig.gdna_prior_kappa_max`

**Add:**
- `CalibrationConfig.gdna_prior_c_base: float = 5.0` — Calibration evidence
  strength for the Dirichlet prior. Combined with a mode-aware baseline
  (+0.5 for VBEM, +0.0 for MAP) per component.

**No change:**
- `EMConfig.mode` — already supports `"vbem"` and `"map"`. The mode selection
  flows through to the C++ solver via the existing `use_vbem` flag.

---

## 8. Cleanup: What to Remove

The mode-aware baseline eliminates all κ-based scaling machinery. Here is the
complete list of code to remove or modify.

### 8.1 Config (`src/rigel/config.py`)

| Change | Details |
|--------|---------|
| Remove field | `gdna_prior_kappa_min: float = 0.001` |
| Remove field | `gdna_prior_kappa_max: float = 0.50` |
| Add field | `gdna_prior_c_base: float = 5.0` |

### 8.2 Locus priors (`src/rigel/locus.py`)

| Change | Details |
|--------|---------|
| Remove parameter | `kappa_min` keyword arg from `compute_locus_priors()` |
| Remove parameter | `kappa_max` keyword arg from `compute_locus_priors()` |
| Add parameter | `c_base: float = 5.0` keyword arg |
| Remove logic | κ interpolation block: `w = (2·SS − 1)²`, `kappa = kappa_min + ...` |
| Remove logic | `C = kappa * N_locus` |
| Rewrite | Prior computation to: `α_gDNA = 0.5 + γ × c_base`, `α_RNA = 0.5 × n_rna_components + (1−γ) × c_base` |
| Note | The function currently returns `(alpha_gdna, alpha_rna)` per locus. The C++ side distributes `alpha_rna` across components proportionally to coverage. With the Jeffreys correction, **the C++ side must add +0.5 per RNA component** (see §6.4). This means `alpha_rna` from Python should represent only the calibration portion: `(1−γ) × c_base`. The +0.5 per component is added in C++. |

### 8.3 Pipeline (`src/rigel/pipeline.py`)

| Change | Details |
|--------|---------|
| Remove | Extraction of `gdna_prior_kappa_min`, `gdna_prior_kappa_max` from config |
| Add | Extraction of `gdna_prior_c_base` from config |
| Update | Call to `compute_locus_priors()` / `_compute_priors()` wrapper |

### 8.4 C++ EM solver (`src/rigel/native/em_solver.cpp`)

| Change | Details |
|--------|---------|
| Modify `compute_ovr_prior_and_warm_start()` | For each RNA component: `prior[j] = 0.5 + alpha_rna × coverage[j] / total_coverage` instead of `prior[j] = alpha_rna × coverage[j] / total_coverage` |
| Modify gDNA prior | Keep as-is: `prior[gdna] = alpha_gdna` (Python passes 0.5 + γ·c_base) |
| Modify gDNA gate | Keep as-is: `if (alpha_gdna <= 0.0) { prior[gdna] = 0.0; }` |
| No change | Warm-start logic unchanged (already uses α_gdna/α_rna ratio) |

**Alternative (simpler):** Python passes the full α values including the +0.5
baseline, and C++ uses them directly. This avoids C++ knowing about the Jeffreys
correction. To do this:
- Python computes per-locus arrays: `alpha_gdna_full[li] = 0.5 + γ × c_base`
  and `alpha_rna_full[li] = n_rna_components × 0.5 + (1−γ) × c_base`
- C++ distributes `alpha_rna_full` across RNA components, adding 0.5 to each
  before distribution... **No, this doesn't work** because C++ distributes
  proportionally to coverage and doesn't know how many components get +0.5.

**Preferred approach:** see §8.4 for the unified C++ implementation with
`const double baseline = use_vbem ? 0.5 : 0.0`.

### 8.5 Calibration (`src/rigel/calibration.py`)

| Change | Details |
|--------|---------|
| Add | Global strand aggregation: sum numerator/denominator before max(0,...) |
| Keep | Per-region γ distribution (proportional to unspliced counts) for per-locus lookup |
| Keep | Density pathway, w blending (these are about **estimating γ**, not about the prior) |
| Note | The w=(2·SS−1)² blending remains for γ estimation. It was previously also used for κ interpolation — that usage is removed. |

### 8.6 Tests

| File | Change |
|------|--------|
| `tests/test_locus_priors.py` | Rewrite `TestKappaInterpolation` → `TestJeffreysCorrection`. Remove κ-specific tests. Add tests for +0.5 baseline, C_base insensitivity, gDNA gate. |
| `tests/test_calibration_integration.py` | Update config construction (remove kappa, add c_base). |
| `tests/test_golden_output.py` | Regenerate golden outputs with `--update-golden`. |
| `tests/scenarios/test_nrna_double_counting.py` | Expect 4 currently-failing tests to pass. |
| `tests/scenarios/test_antisense_intronic.py` | Expect 4 currently-failing tests to pass. |

### 8.7 Debug scripts (no cleanup needed, informational)

The following scripts were created during investigation and can remain in
`scripts/debug/`:
- `analyze_warmstart_vs_prior.py` — warm-start vs prior interaction
- `analyze_digamma_sparsification.py` — digamma bias quantification
- `analyze_jeffreys_correction.py` — Jeffreys correction validation
- `analyze_rectifier_bias.py` — rectifier bias in strand pathway
- `analyze_vbem_c_sensitivity.py` — C sensitivity for standard prior
- `diagnose_pure_rna_calibration.py` — calibration dump for test scenarios
- `analyze_map_vs_vbem_baseline.py` — MAP vs VBEM baseline comparison

---

## 9. Implementation Plan

### Phase 1: Global strand aggregation (calibration fix)

**Goal:** Eliminate the rectifier bias that causes false γ > 0 for pure RNA.

1. Modify `calibrate_gdna()` in `src/rigel/calibration.py`:
   - Aggregate strand numerator/denominator across all eligible regions
   - Apply max(0, ...) to the aggregate
   - Distribute global e_gdna proportionally to per-region unspliced counts
2. Run the 4 currently-failing scenario tests (these may pass with this fix
   alone even before the prior changes, since the root cause is false γ)
3. Run full test suite

### Phase 2: Mode-aware Dirichlet prior

**Goal:** Replace κ×N scaling with mode-aware baseline + C_base.

1. Update `CalibrationConfig` in `config.py`: remove κ fields, add `gdna_prior_c_base`
2. Rewrite `compute_locus_priors()` in `locus.py`: remove κ interpolation,
   output `alpha_gdna = γ × c_base`, `alpha_rna = (1−γ) × c_base`
3. Update `pipeline.py`: wire new config field through to `compute_locus_priors()`
4. Modify `compute_ovr_prior_and_warm_start()` in `em_solver.cpp`:
   add `const double baseline = use_vbem ? 0.5 : 0.0` per eligible component
5. Update gDNA gate in `em_solver.cpp` to override the baseline when disabling
6. Recompile: `pip install --no-build-isolation -e .`

### Phase 3: Test updates

1. Rewrite `tests/test_locus_priors.py` for mode-aware API (test both VBEM
   and MAP baselines)
2. Update `tests/test_calibration_integration.py` for new config fields
3. Add MAP-EM specific tests: verify no baseline for MAP, +0.5 for VBEM,
   both converge to γ
4. Run full test suite, fix any integration issues
5. Regenerate golden outputs: `pytest tests/ --update-golden`
6. Verify all tests pass including the 4 previously-failing scenario tests

### Phase 4: Validation

1. Run synthetic benchmark sweep (`scripts/benchmark/`) to verify no regressions
2. Compare gDNA/mRNA/nRNA relative error across the grid
3. If available: run full-scale benchmarks on Armis2

---

## 10. Testing Strategy

### Unit tests for the mode-aware prior

```python
def test_vbem_jeffreys_baseline():
    """VBEM: every component gets +0.5 pseudocount."""
    # Setup: 3-component locus (gDNA + 2 RNA), γ=0.10, c_base=5.0
    # Expected: α_gDNA = 0.5 + 0.50 = 1.0
    #           α_RNA_1 = 0.5 + 0.45 * cov_1/cov_total * 5.0
    #           α_RNA_2 = 0.5 + 0.45 * cov_2/cov_total * 5.0

def test_map_no_baseline():
    """MAP: components get NO baseline (just calibration share)."""
    # Setup: same locus, γ=0.10, c_base=5.0
    # Expected: α_gDNA = 0.50
    #           α_RNA_1 = 0.45 * cov_1/cov_total * 5.0
    #           α_RNA_2 = 0.45 * cov_2/cov_total * 5.0

def test_c_base_insensitivity_vbem():
    """VBEM: result identical for c_base ∈ {1, 5, 20, 50}."""

def test_c_base_insensitivity_map():
    """MAP: result identical for c_base ∈ {1, 5, 20, 50}."""

def test_gdna_gate_overrides_baseline():
    """When γ=0, gDNA is fully disabled (both VBEM and MAP)."""
    # α_gdna passed as 0 → gate sets prior[gdna] = 0, theta[gdna] = 0

def test_small_locus_both_modes():
    """N=50 fragments converges to γ in both VBEM and MAP."""

def test_mega_locus():
    """N=1M fragments with K=1000 components converges correctly."""

def test_pure_rna_no_false_gdna():
    """With global strand aggregation, pure RNA gets γ=0."""
    # SS=0.65, pure RNA scenario → γ=0 → gDNA disabled
```

### Regression tests

- Golden output tests regenerated
- Scenario tests pass at SS=0.65 and SS=0.90

---

## 11. Risk Analysis

### Low risk: C_base value

C_base has **zero sensitivity** in the degenerate case and rapidly decreasing
sensitivity as strand signal increases. Default of 5.0 is well within the
insensitive range [0.5, 50].

### Low risk: Warm-start unchanged

The warm-start mechanism is not modified. It continues to use the α_gdna/α_rna
ratio, which correctly initializes at γ.

### Medium risk: K-dependent total concentration (VBEM only)

The total concentration C = K/2 + C_base grows linearly with K in VBEM mode
(each component gets +0.5). For a mega-locus with K=1000 components, C ≈ 505.
This is still small relative to N (mega-loci have millions of fragments), so
the prior remains non-informative. In MAP mode, C = C_base (no K scaling),
so this concern does not apply.

### Medium risk: Global strand aggregation

The global aggregation eliminates per-region overdispersion information (the
beta-binomial signal). For most samples this is fine because the calibration
value is used as a prior (not a hard constraint) and the likelihood dominates.
For edge cases with high region-to-region strand variation, the global estimate
may under- or over-estimate γ for specific loci. Monitoring via benchmarks.

### Not a risk: Small γ upward bias (VBEM only)

At γ=0.001, N=2000: VBEM with Jeffreys gives θ_gDNA ≈ 0.003 (absolute error:
4 fragments). This bias arises because the +0.5 baseline is large relative to
γ·C_base = 0.005. MAP mode does not have this issue (no baseline). The bias is
within noise for any practical purpose and vanishes at larger N.

### Not a risk: MAP-EM correctness

Rigel's MAP-EM uses additive pseudocounts (posterior mean), not the textbook
Dirichlet mode formula (α−1). This was verified by inspection of `em_solver.cpp`
line 704–711 and confirmed by simulation (`scripts/debug/analyze_map_vs_vbem_baseline.py`).
MAP-EM converges to prior[j]/C = γ_j with zero baseline for any N ≥ 1. If the
MAP-EM implementation is ever changed to use the true mode formula, the baseline
would need to change to +1.0 (Laplace) — but this is a hypothetical scenario.
