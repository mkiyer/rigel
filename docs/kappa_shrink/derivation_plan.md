# κ_shrink Derivation Plan

## Problem Statement

The EB gDNA prior (`compute_eb_gdna_priors`) uses a shrinkage pseudo-count
κ_shrink to blend per-locus gDNA density toward the globally calibrated density:

```
w = n_locus / (n_locus + κ_shrink)
shrunk_density = w × density_locus + (1 − w) × density_global
```

Currently, `κ_shrink` defaults to `calibration.kappa_strand` — the
Beta-Binomial strand concentration parameter — which is conceptually
unrelated.  It happens to land in a usable range (5–200), but the
conflation is fragile and theoretically unjustified.

**Goal**: Estimate κ_shrink from the calibration data itself, so that
the shrinkage strength automatically adapts to the observed between-region
variance in gDNA density.

## Mathematical Framework

### The Shrinkage Model

The current shrinkage estimator is a standard normal-means model.
For locus $i$ with evidence $n_i$ unspliced reads spread over $L_i$ bp:

$$\hat\lambda_i = \frac{n_i}{L_i}$$

is the raw per-locus density estimate, and the shrinkage weight is:

$$w_i = \frac{n_i}{n_i + \kappa}$$

so the shrunk estimate is:

$$\tilde\lambda_i = w_i \hat\lambda_i + (1 - w_i) \lambda_G$$

This is the posterior mean of a conjugate normal-normal (or
Gamma-Gamma) hierarchical model where κ controls the "prior sample
size" — the number of pseudo-observations of the global density
to admix with each locus's local estimate.

### Interpretation of κ_shrink

In a Gamma-Poisson framework:

- Each locus has a true gDNA rate $\lambda_i \sim \text{Gamma}(\alpha, \beta)$
- The observed count is $x_i | \lambda_i \sim \text{Poisson}(\lambda_i \cdot L_i)$
- The posterior mean is $\frac{x_i + \alpha}{L_i + \beta}$

The shrinkage weight $w_i = L_i / (L_i + \beta)$ is analogous to
our formula with $n_i \approx L_i$ (using evidence, not raw reads)
and $\kappa \approx \beta$.

The key quantity is $\beta$ — the inverse-rate of the Gamma prior.
It controls the strength of the prior:

- **Small β (small κ_shrink)**: Trust the local data; minimal pooling.
  Appropriate when gDNA density truly varies across loci.
- **Large β (large κ_shrink)**: Strong pooling toward the global
  mean.  Appropriate when gDNA density is uniform genome-wide
  (clean library, random shearing).

## Estimation Strategy: Method of Moments from Calibration Posteriors

### Available Data

After calibration EM converges, we have:
- **Per-region posteriors** $\gamma_r = P(\text{not expressed} | \text{data})$
- **Per-region summary stats**: $n_{\text{total},r}$, $L_r$, etc.
- **Global gDNA density** $\lambda_G$
- **Expressed density** $\lambda_E$
- **Mixing proportion** $\pi$

### Method of Moments Derivation

Define the weighted sample of log-densities for gDNA-classified regions:

$$d_r = \log\left(\frac{n_r}{L_r} + \epsilon\right)$$

Already computed as `log_d` in the calibration code.

For regions with high $\gamma_r$ (gDNA-dominant), the variance of $d_r$
reflects the **between-region** variability of gDNA density.

In the hierarchical model:

$$\text{Var}(d_r) = \underbrace{\sigma^2_{\text{between}}}_{\text{true variability}} + \underbrace{\sigma^2_{\text{within}}}_{\text{sampling noise}}$$

The within-region sampling variance for a Poisson count $n_r \sim
\text{Poisson}(\lambda_r L_r)$ on the log scale is approximately
$1/n_r$ (delta method).

So:

$$\sigma^2_{\text{between}} = \text{Var}_\gamma(d_r) - \overline{1/n_r}$$

where $\text{Var}_\gamma$ is the γ-weighted variance (exactly the
`var_g` computed in the calibration M-step), and $\overline{1/n_r}$
is the γ-weighted mean of $1/n_r$.

### Converting to κ_shrink

On the log-density scale, the shrinkage weight translates to:

$$\kappa_{\text{shrink}} = \frac{\overline{n} \cdot \sigma^2_{\text{within}}}{\sigma^2_{\text{between}}}$$

Intuitively: when between-region variance is small relative to
sampling noise, κ is large (strong pooling).  When between-region
variance is large (true heterogeneity), κ is small (trust local
estimates).

A simpler alternative: define κ_shrink directly as the effective
number of "prior pseudo-reads" worth of global density data.
Using the coefficient of variation of gDNA densities:

$$\text{CV}^2 = \frac{\sigma^2_{\text{between}}}{\mu_G^2}$$

then:

$$\kappa_{\text{shrink}} = \frac{1}{\text{CV}^2}$$

This gives "pseudo-sample-size equivalent" of the prior. With
CV = 0.1 (densities within 10% of global), κ = 100 (strong
pooling).  With CV = 1.0, κ = 1 (minimal pooling).

### Practical Formula

```python
def estimate_kappa_shrink(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    log_d: np.ndarray,
    eligible: np.ndarray,
) -> float:
    """Estimate shrinkage κ from calibration region posteriors.

    Uses Method of Moments on γ-weighted log-density variance:
      σ²_between = Var_γ(log_d) - E_γ[1/n]
      κ_shrink   = 1 / max(σ²_between, ε)
    """
    g = gamma[eligible]
    d = log_d[eligible]
    n = stats["n_total"][eligible].astype(np.float64)

    # γ-weighted mean and variance of log-density
    w_sum = max(g.sum(), 1e-10)
    mu = (g * d).sum() / w_sum
    var_total = (g * (d - mu) ** 2).sum() / w_sum

    # γ-weighted mean of inverse counts (within-region noise)
    safe_n = np.maximum(n, 1.0)
    var_within = (g / safe_n).sum() / w_sum

    # Between-region variance (floor at small positive)
    var_between = max(var_total - var_within, 0.01)

    # κ as inverse between-region variance  
    kappa = 1.0 / var_between

    # Bound to reasonable range
    return float(np.clip(kappa, 1.0, 500.0))
```

## Implementation Plan

### Step 1: Add `estimate_kappa_shrink()` to `calibration.py`

Place it alongside `estimate_kappa_marginal()`.  Uses the same
calibration arrays (`stats`, `gamma`, `log_d`, `eligible`).

### Step 2: Store in `GDNACalibration`

Add a new field:
```python
kappa_shrink: float
```

This field is independent of `kappa_strand`.

### Step 3: Wire through to `compute_eb_gdna_priors()`

In `locus.py`, the fallback changes from:
```python
k_shrink = calibration.kappa_strand if kappa_shrink is None else kappa_shrink
```
to:
```python
k_shrink = calibration.kappa_shrink if kappa_shrink is None else kappa_shrink
```

### Step 4: Compute during calibration EM

At the end of `calibrate_gdna()`, after the final E-step:
```python
kappa_shrink = estimate_kappa_shrink(stats, final_gamma, log_d, eligible)
```

Pass it into the `GDNACalibration` constructor.  The two bailout
paths should use a sensible default (e.g., 20.0 — moderate pooling).

### Step 5: Update `EMConfig.gdna_kappa_shrink`

Keep as a manual override.  When `None` (default), the calibrated
value is used.  This allows users to force a specific shrinkage
strength via `--gdna-kappa-shrink`.

### Step 6: Update logging

Log `κ_shrink` alongside `κ_strand` in the calibration summary line
and the EB gDNA initialization log.

### Step 7: Tests

- Unit test `estimate_kappa_shrink()` with synthetic scenarios:
  - Uniform density → large κ (strong pooling)
  - Variable density → small κ (weak pooling)
  - Degenerate cases: all γ=0, single region, etc.
- Integration test: verify that `GDNACalibration.kappa_shrink`
  is populated and positive after calibration.
- Regression: golden output regen (κ_shrink changes initialization).

### Step 8: Remove the kappa_strand fallback

Once `kappa_shrink` is estimated independently, remove the
`calibration.kappa_strand` fallback from `compute_eb_gdna_priors()`.
The code should fail loudly if `kappa_shrink` is missing rather
than silently reusing a wrong parameter.

## Risks and Mitigations

1. **Variance estimate instability with few gDNA regions**: When
   very few regions have high γ (e.g., clean RNA-seq), the MoM
   estimate is noisy.  Mitigation: floor κ_shrink at 1.0 and cap
   at 500.0.  Alternatively, use a prior on κ_shrink itself
   (hyper-shrinkage), but this adds complexity.

2. **Log-density variance is not on the right scale**: The MoM
   derivation assumes normal approximation on the log scale.  For
   very low-count regions, the Poisson approximation breaks down.
   Mitigation: filter to regions with $n \geq 5$ for the MoM
   estimate.

3. **Correlation with κ_strand**: If the library has extreme strand
   bias, the hybrid density estimator already accounts for this
   through `compute_gdna_density_hybrid()`.  The between-region
   variance in the raw density (not the strand-corrected density)
   could be inflated by strand effects.  Mitigation: use the
   calibration's posterior-weighted density (which already
   integrates strand, density, and FL signals).

## Alternative Approaches (Considered)

### REML / Marginal Likelihood
Maximize the restricted marginal likelihood of a Gamma-Poisson
hierarchical model.  More principled but significantly more complex
to implement and harder to debug.  The MoM approach is adequate
for a pseudo-count and matches the informality of the shrinkage
formula.

### Cross-validation
Hold out regions and minimize prediction error.  Principled but
expensive and unnecessary given that we have hundreds–thousands of
regions.

### Fixed default
Use a fixed κ_shrink (e.g., 20 or 50).  Simpler, but doesn't
adapt to the data.  The current kappa_strand stand-in is
effectively this — accidentally data-adaptive through completely
wrong reasoning.

## Recommendation

Implement the MoM approach (Steps 1–8).  It is:
- **Simple**: ~15 lines of new code
- **Principled**: derives from the hierarchical model
- **Adaptive**: automatically adjusts pooling strength
- **Testable**: clear expected behavior for uniform vs variable densities
- **Backwards-compatible**: the `--gdna-kappa-shrink` CLI override remains
