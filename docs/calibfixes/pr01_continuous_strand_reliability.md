# PR 1: Continuous Strand Reliability

## Goal

Replace hard strand-source decisions with a continuous, exact-at-low-depth reliability weight. Strand-derived gDNA mass should be downweighted when the observed strand imbalance is compatible with all RNA, and trusted smoothly when the mixed RNA/gDNA model is much more likely.

This PR repairs source mass only. It does not change capture exposure, ESS policy, transcript floors, or unstranded capture-on source modeling.

## Non-goals

- Do not add any transcript-level RNA prior floor.
- Do not tune `kappa_d` downward to force overdispersion.
- Do not derive RNA/gDNA source mass from latent-state probabilities.
- Do not introduce capture exposure caps or pseudocount constants.
- Do not attempt to solve unstranded capture-on source splitting.

## Files to edit

| File | Purpose of edit |
|---|---|
| `src/rigel/strand_model.py` | Expose same/opposite counts and minor-rate posterior parameters from the trained strand model. |
| `src/rigel/calibration/strand_summary.py` | Carry the posterior parameters through the public calibration summary. |
| `src/rigel/calibration/strand_deconv.py` | Add exact marginal-likelihood strand reliability and attach it to compartment deconvolution outputs. |
| `src/rigel/calibration/calibration_iteration.py` | Weight strand-derived source mass before building `PriorMassDeconvolution`. |
| `src/rigel/calibration/_result.py` | Add compact reliability diagnostics to `CalibrationResult.to_summary_dict()`. |
| `tests/test_strand_deconv.py` | Unit tests for exact predictive probabilities, continuity, and edge cases. |
| `tests/test_calibration_iteration.py` | Integration tests for reliability-weighted source handoff. |
| Optional: `tests/test_strand_summary.py` | Small tests for summary fields if the existing test layout is clearer with a new file. |

## Public contract after this PR

`PriorMassDeconvolution.gdna_unspliced_mean` remains mass-conserving and explicit. When strand evidence is available, it should represent source-reliable gDNA mass, not raw posterior gDNA mass. RNA mass is still derived as:

```text
rna_unspliced_mean = unspliced_total - gdna_unspliced_mean
```

The reliability weight must be stored for diagnostics. A downstream reviewer should be able to distinguish:

- no strand information;
- low depth with broad predictive uncertainty;
- all-RNA-compatible strand imbalance;
- strong mixed-source strand imbalance.

### Model-symmetry decision

This PR must remove the statistical asymmetry between the all-RNA null and the mixed RNA/gDNA slab. The all-RNA model integrates uncertainty in the RNA minor-orientation rate with a beta-binomial predictive distribution. The mixed model must use the same RNA beta-binomial component, not a fixed point-estimate binomial RNA component.

The PR therefore has two linked responsibilities:

1. compute `P(data | H1)` with beta-binomial RNA uncertainty and symmetric beta-binomial gDNA strand balance;
2. compute `E[D | data, H1]` under that same slab likelihood before multiplying by `P(H1 | data)`.

Using `sigmoid(log_bf) * old_binomial_slab_mean` is not acceptable, because the posterior model for the slab mean would not match the Bayes factor. If backward compatibility requires keeping a legacy binomial helper, give the new beta-binomial slab helper a distinct name and make `PriorMassDeconvolution` consume the new slab mean.

Audit these existing code paths while implementing the PR:

- `strand_log_likelihood_d_grid()`;
- `_exact_posterior_R()`;
- `_summarize_normal()`;
- `deconvolve_regions_by_strand()` and `deconvolve_compartments_by_strand()`;
- `calibration/integration.py` call sites that combine density and strand likelihoods.

Any call site that contributes source mass or posterior gDNA summaries should either pass the RNA minor-rate posterior or explicitly document why it remains a legacy point-binomial approximation.

## Step 1: expose minor-rate posterior from `StrandModel`

Add properties to `StrandModel`:

```python
@property
def n_minor(self) -> int:
    return min(self.n_same, self.n_opposite)

@property
def n_major(self) -> int:
    return max(self.n_same, self.n_opposite)

@property
def minor_rate_posterior_alpha(self) -> float:
    return float(self.n_minor + 1)

@property
def minor_rate_posterior_beta(self) -> float:
    return float(self.n_major + 1)
```

Use the same Laplace/Beta(1, 1) convention already used by `strand_specificity_ci_epsilon()`. Do not change the existing MLE behavior of `p_r1_sense`, `p_r1_antisense`, or `strand_specificity`.

If `StrandModel.to_dict()` includes derived strand fields, add these diagnostic keys without removing old keys:

```text
n_same
n_opposite
n_minor
n_major
minor_rate_posterior_alpha
minor_rate_posterior_beta
```

## Step 2: extend `StrandSummary`

Add frozen dataclass fields with backward-compatible defaults:

```python
n_same: int = 0
n_opposite: int = 0
minor_rate_alpha: float = 1.0
minor_rate_beta: float = 1.0
```

Validation rules:

- counts are non-negative integers;
- alpha and beta are finite and positive;
- when counts are present, `minor_rate_alpha == min(n_same, n_opposite) + 1` and `minor_rate_beta == max(n_same, n_opposite) + 1` unless the summary was explicitly built as uninformative;
- `uninformative()` returns `p_r1_sense=0.5`, `n_observations=0`, `n_same=0`, `n_opposite=0`, `minor_rate_alpha=1.0`, `minor_rate_beta=1.0`.

Update `StrandSummary.from_model()` to copy `model.n_same`, `model.n_opposite`, and the new posterior alpha/beta fields.

Add convenience properties:

```python
@property
def minor_rate_mean(self) -> float:
    return self.minor_rate_alpha / (self.minor_rate_alpha + self.minor_rate_beta)

@property
def minor_rate_concentration(self) -> float:
    return self.minor_rate_alpha + self.minor_rate_beta
```

## Step 3: add reliability dataclasses and flags in `strand_deconv.py`

Add imports for vectorized stable probability calculations:

```python
from scipy.special import expit, log_ndtr, logsumexp
```

Use `expit()` for array sigmoid conversion. Use `log_ndtr()` inside the large-count Normal approximation so deep-tail interval probabilities do not underflow to `-inf`.

Add new flags near the existing strand-deconvolution flags:

```python
FLAG_LOW_STRAND_RELIABILITY: int = 1 << 5
FLAG_RELIABILITY_APPROX: int = 1 << 6
```

If `uint8` flags become too tight, promote the per-compartment flags to `uint16` in this PR and update tests accordingly.

Add a dataclass:

```python
@dataclass(frozen=True, slots=True)
class StrandReliabilityEstimate:
    reliability: np.ndarray          # float32[R], in [0, 1]
    log_bayes_factor: np.ndarray     # float32[R]
    log_p_all_rna: np.ndarray        # float32[R]
    log_p_mixed: np.ndarray          # float32[R]
    slab_gdna_mean: np.ndarray       # float32[R], E[D | data, H1]
    k_minor: np.ndarray              # float32[R]
    n_total: np.ndarray              # float32[R]
    flags: np.ndarray                # uint16[R]
```

Extend `RegionGdnaChannelEstimate` with compartment reliability arrays:

```python
contained_reliability: np.ndarray
contained_log_bayes_factor: np.ndarray
boundary_left_reliability: np.ndarray
boundary_left_log_bayes_factor: np.ndarray
boundary_right_reliability: np.ndarray
boundary_right_log_bayes_factor: np.ndarray
```

Either update the existing `contained_mean` / `boundary_*_mean` fields to be slab means under the beta-binomial RNA mixed model, or add explicit fields:

```python
contained_slab_gdna_mean: np.ndarray
boundary_left_slab_gdna_mean: np.ndarray
boundary_right_slab_gdna_mean: np.ndarray
```

Prefer updating the existing mean fields if all current consumers expect strand-derived gDNA posterior means. Add explicit fields only if another call site still needs the old point-binomial posterior for diagnostics.

Keep these arrays shape `(R,)`, finite, and clipped to valid ranges in `__post_init__` if a post-init exists. If no post-init exists for `RegionGdnaChannelEstimate`, add a small private validator used by `deconvolve_compartments_by_strand()` before construction.

## Step 4: implement exact all-RNA predictive likelihood

Use the existing `_log_beta_binom_pmf()` helper. For each region/compartment:

```text
n = round(k_sense + k_antisense)
k_minor = k_antisense if p_r1_sense >= 0.5 else k_sense
log_p0 = log BetaBinomial(k_minor | n, minor_rate_alpha, minor_rate_beta)
```

Important details:

- Clamp rounded counts into `[0, n]` exactly as `_exact_posterior_R()` does.
- If `n == 0`, set reliability to `0`, `log_p0 = 0`, `log_p_mixed = 0`, and flag ineligible.
- If strand contrast is below `STRAND_CONTRAST_NUMERICAL_FLOOR`, set reliability to `0` and keep source split inactive.
- Do not use a Normal approximation for the all-RNA predictive distribution at low depth.

## Step 5: implement mixed-model marginal likelihood

Preferred production target:

```text
p1 = P(k_minor | n, mixed RNA/gDNA)
   = average_D P(k_minor | n, D, alpha_q, beta_q, kappa_d)
```

Use the same uniform prior over `R`/`D` that `_exact_posterior_R()` uses. Fold the observation into minor-orientation coordinates before computing both `H0` and `H1`:

```text
k_minor = k_antisense if p_r1_sense >= 0.5 else k_sense
q ~ Beta(minor_rate_alpha, minor_rate_beta)
```

For a fixed slab value `D=d` and `R=n-d`:

```text
J_rna ~ BetaBinomial(R, minor_rate_alpha, minor_rate_beta)
J_dna ~ BetaBinomial(D, kappa_d / 2, kappa_d / 2)
K_minor = J_rna + J_dna
```

So:

```text
log L(D=d) = logsumexp_j [
    log BetaBinomial(j | R, minor_rate_alpha, minor_rate_beta)
  + log BetaBinomial(k_minor - j | D, kappa_d / 2, kappa_d / 2)
]
```

This replaces the current binomial RNA component in `strand_log_likelihood_d_grid()` for PR01 source/reliability work. The gDNA component remains the existing symmetric beta-binomial with mean `0.5`.

Exact path for `n <= MAX_EXACT_POSTERIOR_N`:

```python
d_grid = np.arange(n + 1, dtype=np.int64)
log_l = strand_log_likelihood_d_grid_minor_beta_binom(
    k_minor,
    n,
    d_grid,
    kappa_d=kappa_d,
    minor_rate_alpha=minor_rate_alpha,
    minor_rate_beta=minor_rate_beta,
)
log_p1 = logsumexp(log_l) - math.log(n + 1)
slab_post = np.exp(log_l - logsumexp(log_l))
slab_gdna_mean = float(np.sum(d_grid * slab_post))
```

Implementation choices:

- Option A: add `minor_rate_alpha` and `minor_rate_beta` keyword-only parameters to `strand_log_likelihood_d_grid()` and preserve the old binomial behavior only when they are omitted.
- Option B: add a new helper such as `strand_log_likelihood_d_grid_minor_beta_binom()` and use it for PR01 reliability/source mass.

Prefer Option B if it makes the audit safer. Whichever option is chosen, exact slab marginal likelihood and exact slab posterior mean must use the same likelihood array.

Approximate path for large `n`:

1. Use a fixed deterministic grid, for example `np.linspace(0, n, 257)`, rounded and uniqued.
2. Give each grid point an interval-width weight so the grid approximates the uniform prior over `D`.
3. Approximate `K_minor | D` by a continuity-corrected Normal with:

```text
R = n - D
q_mean = minor_rate_alpha / (minor_rate_alpha + minor_rate_beta)
q_conc = minor_rate_alpha + minor_rate_beta
mean = R * q_mean + 0.5 * D
var_rna = R * q_mean * (1 - q_mean) * (q_conc + R) / (q_conc + 1)
var_dna = 0.25 * D * (D + kappa_d) / (1 + kappa_d)
var = max(var_rna + var_dna, tiny)
```

4. Compute `log P(k - 0.5 <= K <= k + 0.5)` with a `log_ndtr()`-based helper, not `log(norm.cdf(high) - norm.cdf(low))`.
5. `log_p1 = logsumexp(log_grid_prob + log_grid_weight)`.
6. Compute approximate `slab_gdna_mean` from the normalized grid posterior weights.
7. Set `FLAG_RELIABILITY_APPROX`.

Add a private helper for stable Normal interval probabilities:

```python
def _log_normal_interval_prob(low: np.ndarray, high: np.ndarray) -> np.ndarray:
    """Return log(P(low <= Z <= high)) for standard Normal Z."""
```

Implementation notes:

- For negative-tail intervals, use `log_ndtr(high)` and `log_ndtr(low)` with a `logdiffexp` helper.
- For positive-tail intervals, use survival symmetry: `P(low <= Z <= high) = Phi(-low) - Phi(-high)` and `log_ndtr()` on the negated bounds.
- For intervals crossing zero, either branch to the more stable side or still use `logdiffexp` with `log_ndtr()` if the difference is well-conditioned.
- Clamp only as a final numerical guard, for example to `np.finfo(float).tiny`, and flag non-finite results.

Keep the exact path unit-tested. The approximate path only needs sanity tests: finite, monotone in strong evidence examples, and no discontinuity at `MAX_EXACT_POSTERIOR_N` large enough to affect source mass qualitatively.

## Step 6: convert likelihoods into a smooth weight

Compute:

```text
log_bf = log_p_mixed - log_p_all_rna
w_strand = sigmoid(log_bf)
```

Use the vectorized stable sigmoid from SciPy:

```python
w = expit(log_bf)
```

Do not write a scalar `if log_bf >= 0` loop over regions. `expit()` implements the same stability guard in C and works directly on `(R,)` arrays.

Then apply conservative gates only for structural absence of evidence:

- ineligible region -> `w = 0`;
- near-unstranded -> `w = 0`;
- non-finite likelihood -> `w = 0` and flag;
- do not threshold valid finite weights.

Add `FLAG_LOW_STRAND_RELIABILITY` when `0 < w < 0.5`. This flag is diagnostic only; it must not change the value.

## Step 7: attach reliability to compartment deconvolution

Update `_single_compartment_estimate()` or `deconvolve_compartments_by_strand()` so each compartment gets a matching `StrandReliabilityEstimate` using the same `k_sense`, `k_antisense`, `n_total`, `eligible`, `p_r1_sense`, and `kappa_d` values used for deconvolution.

Update the function signatures so the RNA minor-rate posterior is available where the slab is computed:

```python
def deconvolve_regions_by_strand(
    counts: StrandRegionCounts,
    *,
    kappa_d: float,
    strand_summary: StrandSummary | None = None,
    ...
) -> RegionGdnaEstimate: ...

def deconvolve_compartments_by_strand(
    counts: CompartmentStrandCounts,
    *,
    kappa_d: float,
    strand_summary: StrandSummary | None = None,
) -> RegionGdnaChannelEstimate: ...
```

If `strand_summary is None`, synthesize a concentrated beta posterior from `counts.p_r1_sense` only for backward-compatible legacy callers, and set a diagnostic flag indicating point-estimate fallback. Production calibration in `_orchestrator.py` must pass the real `StrandSummary`.

For `deconvolve_regions_by_strand()`, either:

- return reliability in `RegionGdnaEstimate`; or
- keep `RegionGdnaEstimate` unchanged for public compatibility but ensure the compartment wrapper computes both reliability and slab means from the new beta-binomial slab helper.

Prefer the second option for a narrow PR if it avoids churn. The hard requirement is that `RegionGdnaChannelEstimate` fields consumed by `build_prior_mass_deconvolution()` use the beta-binomial RNA slab mean, not the old point-binomial slab mean.

## Step 8: weight source mass in `build_prior_mass_deconvolution()`

Change only the strand-channel branch.

Current behavior:

```python
gdna = contained + left + right
```

New behavior:

```python
contained_gdna = strand_channels.contained_reliability * strand_channels.contained_mean
left_gdna = strand_channels.boundary_left_reliability * strand_channels.boundary_left_mean
right_gdna = strand_channels.boundary_right_reliability * strand_channels.boundary_right_mean
gdna = contained_gdna + left_gdna + right_gdna
```

Here `contained_mean` and boundary means must be `E[D | data, H1]` under the beta-binomial RNA mixed model from Step 5. If explicit `*_slab_gdna_mean` fields were added, use those fields instead of the legacy means.

Precision should reflect both posterior sharpness and source reliability:

```python
precision = np.maximum.reduce([
    contained_precision * contained_reliability,
    boundary_left_precision * boundary_left_reliability,
    boundary_right_precision * boundary_right_reliability,
])
```

Keep the final clipping to `[0, unspliced_total]` and derive RNA from float32-converted totals exactly as the current function does.

## Step 9: diagnostics

Update `CalibrationResult.to_summary_dict()` so the strand-deconvolution block reports compact reliability stats:

```text
contained_reliability: min/p25/p50/p75/p90/p99/max/mean
boundary_left_reliability: same
boundary_right_reliability: same
n_regions_low_reliability
n_regions_approx_reliability
n_regions_near_unstranded
```

If summary helpers already exist in `_result.py`, reuse them. Do not add large per-region arrays to `summary.json`.

Optional debug table for benchmark analysis:

```text
region_id
signature
unspliced_total
contained_mean
contained_reliability
contained_log_bayes_factor
prior_gdna_unspliced_mean
prior_rna_unspliced_mean
flags
```

Keep this as an opt-in diagnostic script, not default quant output.

## Tests

### `tests/test_strand_deconv.py`

Add tests for:

1. Exact all-RNA predictive probability matches a direct beta-binomial calculation for `n=10`, `k_minor=0..3`, `alpha=1`, `beta=100`.
2. Exact mixed likelihood matches a direct convolution of RNA beta-binomial and gDNA beta-binomial components for several `D` values.
3. Exact slab posterior mean uses the same likelihood array as `log_p_mixed`; verify a hand-computed small example.
4. Reliability is continuous: for fixed `n=50`, increasing `k_minor` by one changes `w_strand` smoothly and never jumps from exactly `0` to exactly `1`.
5. Pure RNA example: `n=100`, minor count close to expected under `q=0.01` gives low reliability.
6. Strong mixed example: `n=100`, strand balance near 50/50 in a highly stranded protocol gives reliability near `1`.
7. Small strand-training set: `minor_rate_alpha=1`, `minor_rate_beta=3` produces lower reliability than `minor_rate_alpha=1`, `minor_rate_beta=100` for the same region without causing a slab/null variance mismatch.
8. Near-unstranded protocol: `p_r1_sense=0.5001` gives reliability `0` and sets the near-unstranded/inactive flag.
9. Large-count approximate path returns finite values, uses beta-binomial RNA variance, and sets `FLAG_RELIABILITY_APPROX`.
10. Deep-tail Normal interval probabilities remain finite by exercising the `log_ndtr()` helper with extreme z-scores.
11. Vectorized `expit()` path accepts and returns `(R,)` arrays without scalar branching.

### `tests/test_calibration_iteration.py`

Add a focused test for `build_prior_mass_deconvolution()`:

- build a fake `RegionGdnaChannelEstimate` with `contained_mean=[80]`, `contained_reliability=[0.25]`, no boundary mass, and `unspliced_total=[100]`;
- assert `gdna_unspliced_mean == 20` and `rna_unspliced_mean == 80` within float32 tolerance;
- assert mass conservation remains exact after float32 conversion.

Add a second test with reliability `0` proving strand-derived gDNA does not leak into priors when evidence is all-RNA-compatible.

### `tests/test_strand_summary.py` or equivalent

- `StrandSummary.from_model()` copies same/opposite counts.
- `StrandSummary.uninformative()` has Beta(1, 1) minor-rate posterior.
- invalid alpha/beta values raise `ValueError`.

## Validation commands

```bash
conda activate rigel && ruff check src/rigel/strand_model.py src/rigel/calibration/strand_summary.py src/rigel/calibration/strand_deconv.py src/rigel/calibration/calibration_iteration.py tests/test_strand_deconv.py tests/test_calibration_iteration.py
conda activate rigel && pytest tests/test_strand_deconv.py tests/test_calibration_iteration.py -v
```

Run broader calibration tests before merge:

```bash
conda activate rigel && pytest tests/test_calibrate.py tests/test_calibration_result.py tests/test_per_locus_gdna_mass.py tests/test_pipeline_wiring.py -v
```

If `calibration/integration.py` is changed to pass the RNA minor-rate posterior into density/strand fusion, include its focused tests too:

```bash
conda activate rigel && pytest tests/test_calibration_integration.py tests/test_calibration_prior.py -v
```

## Benchmark gate

After PR 1, before PR 2, run the existing eight-condition suite only as a diagnostic checkpoint. Expected behavior:

- stranded capture-off improves or maintains source split;
- stranded capture-on source mass becomes less noisy, but `A_r` is still old and may remain wrong;
- unstranded scenarios should be unchanged except for diagnostics;
- no-gDNA stranded capture-on should not show large false `prior_gdna_unspliced_mean` from strand noise.

Do not judge final capture-on performance until PR 2 lands.

## Review checklist

- [ ] `PriorMassDeconvolution` still conserves mass per region.
- [ ] Reliability is continuous; no p-value threshold controls source mass.
- [ ] Low-depth exact path uses beta-binomial probabilities for both all-RNA `H0` and mixed-model RNA under `H1`.
- [ ] The slab mean `E[D | data, H1]` and `log_p_mixed` are computed from the same likelihood array.
- [ ] Large-count approximation uses beta-binomial RNA variance, not point-binomial RNA variance.
- [ ] Sigmoid conversion uses vectorized `scipy.special.expit()`.
- [ ] Deep-tail Normal interval probabilities use `scipy.special.log_ndtr()` or an equivalent log-space helper.
- [ ] Near-unstranded data turns reliability off without creating false gDNA.
- [ ] No transcript-level RNA floor was added.
- [ ] No latent-state probability is used as an RNA/gDNA source label.
