# Negative-Binomial density channel: implementation plan

**Date:** 2026-04-22
**Status:** Design proposal. Implements §4.2 B.1 of [channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md).
**Scope:** Replace the per-region Poisson kernel inside the calibration density LLR with a Negative Binomial that absorbs un-modeled spatial variance in the per-bp gDNA capture rate. Borrows estimation strategy from DESeq2.

---

## 1. Conceptual model

### 1.1 The current model and its asymmetry

The current density channel ([`_count_llr_poisson_ln`](src/rigel/calibration/_em.py#L134)) is a two-component mixture:

$$P(k_i | G) = \text{Poisson}(k_i;\;\lambda_G E_i)$$

$$P(k_i | R) = \int_0^\infty \text{Poisson}(k_i;\;(\lambda_G + \mu) E_i)\,\text{LogNormal}(\mu;\;\mu_R, \sigma_R^2)\,d\mu$$

Note that the **R class is already a compound (overdispersed) distribution** — it is Poisson-LogNormal, which has variance ≫ mean. The R class can therefore accommodate arbitrary cross-region heterogeneity in expression rate μ. This is why the R class is rarely the source of mis-specification.

The **G class is a rigid Poisson** with a single global rate λ_G. It assumes that *every* gDNA-only region in the genome produces fragments at the same per-bp rate, with only Poisson sampling noise around it. This is the broken assumption under hybrid capture: λ_G varies 100–140× between on-target and off-target regions.

So the density channel suffers from an **asymmetry**: R accommodates heterogeneity, G does not. The fix is to give G the same kind of compound structure.

### 1.2 The proposed model

Replace the G class with a Gamma-Poisson (= Negative Binomial) compound:

$$\lambda_{G,i} \sim \text{Gamma}(\text{shape}=1/\phi_G,\;\text{scale}=\phi_G \cdot \lambda_G) \quad\text{(}\mathbb{E}[\lambda_{G,i}] = \lambda_G,\;\text{Var}[\lambda_{G,i}] = \phi_G \lambda_G^2\text{)}$$

$$k_i \mid \lambda_{G,i} \sim \text{Poisson}(\lambda_{G,i} \cdot E_i)$$

Marginalizing over λ_{G,i} gives:

$$P(k_i | G) = \text{NegBin}(k_i;\;\text{mean}=\lambda_G E_i,\;\text{dispersion}=\phi_G)$$

with the DESeq2 parameterization:

$$\text{Var}(k_i | G) = \lambda_G E_i + \phi_G (\lambda_G E_i)^2$$

A single dispersion parameter φ_G ≥ 0 controls cross-region heterogeneity in the gDNA capture rate. As φ_G → 0, this collapses back to Poisson — so WGS-like data and synthetic benchmarks degrade gracefully.

### 1.3 Does NB apply to both R and G?

**Yes — both classes get an NB kernel, but for different reasons.**

The R class becomes:

$$P(k_i | R) = \int_0^\infty \text{NegBin}(k_i;\;(\lambda_G + \mu) E_i,\;\phi_G)\,\text{LogNormal}(\mu;\;\mu_R, \sigma_R^2)\,d\mu$$

This is **NB with the same φ_G**, integrated over a LogNormal expression rate μ. Two reasons for keeping a single shared φ_G:

1. **Conceptual:** φ_G models per-region variation in the gDNA capture/mappability process. That process operates on every region regardless of whether it expresses RNA. So the R class inherits the same Gamma-distributed background rate.
2. **Identifiability:** the LogNormal in R already absorbs cross-region expression heterogeneity. Adding a separate φ_R would be confounded with σ_R² (both inflate the R-class variance). One free overdispersion parameter per class is the maximum that remains identifiable from a single library.

The R class is already a doubly-compound distribution (Gamma over λ_G, LogNormal over μ, then Poisson on top). The NB-LogNormal compound has no closed form but is straightforward to evaluate via the same 21-node Gauss-HermiteE quadrature already used today — we just swap the Poisson PMF inside the integrand for an NB PMF.

### 1.4 Are we doing mixture deconvolution of two NB distributions?

**Effectively yes**, with two technical caveats:

- The G class is a single NB. The R class is an *NB-LogNormal compound* (NB integrated over the LogNormal prior on μ). So strictly: we are deconvolving NB versus NB-LogNormal. In practice both are evaluated as discrete probabilities at each observed `k_i`.
- The deconvolution is performed per region inside the EM E-step (LLR contribution), with parameters (λ_G, μ_R, σ_R, φ_G) updated jointly in the M-step. This is NOT a separate two-NB-mixture EM problem; it's a drop-in replacement for the Poisson kernel inside the existing two-component mixture EM in `_em.py`.

---

## 2. What we borrow from DESeq2

DESeq2 [Love et al. 2014, Genome Biology] is the canonical NB-based count model in genomics. Its estimation strategy is well-tested for low-count, high-overdispersion data and applies directly here. Key elements we adopt:

### 2.1 Parameterization

- **Mean–dispersion form:** `k ~ NB(μ, α)`, `Var(k) = μ + α·μ²`. The dispersion α is on a relative-variance scale (CV² of the underlying Gamma rate), which is interpretable and stable across orders-of-magnitude variation in μ.
- We use **`φ_G`** for the symbol to avoid confusion with the Bayesian prior π. Identical mathematically to DESeq2's `α`.

### 2.2 Cox-Reid adjusted profile likelihood

When estimating dispersion jointly with the mean(s), the marginal MLE of φ is biased downward in finite samples. DESeq2 corrects this with the **Cox-Reid bias adjustment**:

$$\ell_{\text{CR}}(\phi) = \ell(\phi, \hat\mu(\phi)) - \tfrac{1}{2}\log\det(I_{\mu\mu})$$

where `I_{μμ}` is the Fisher information of the mean parameters at the profile maximum. For our setup with a single shared mean λ_G across N eligible regions, `I_{μμ}` has a simple analytic form and the correction is one extra term per dispersion update.

This matters: without Cox-Reid, the MLE for φ_G can collapse to zero whenever the data happens to look "Poisson-enough" by chance, which would defeat the purpose of the model.

### 2.3 Single shared dispersion (no per-region φ)

DESeq2's full pipeline fits one φ per gene, then shrinks toward a fitted trend `φ(μ)`. **We do not.** We fit a single shared `φ_G` across all eligible regions because:

- We have one library (one observation per region), not multiple replicates. Per-region φ is unidentifiable.
- The mis-specification we want to absorb is genome-scale (capture-targeting heterogeneity), so a shared scalar captures the bulk of it.
- One scalar parameter is the minimum that fixes the failure mode in [§3.1 of channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md#31-the-density-llr-is-systematically-wrong) without adding identifiability risk.

**Future option:** B.2 of the parent document (`φ_G(μ)` trend or per-window φ) is the natural next refinement once B.1 is shown to work.

### 2.4 What we do NOT borrow

- **Empirical Bayes shrinkage of φ.** Not applicable — we have one shared φ already.
- **Wald / LRT testing.** We are not testing differential expression; we are fitting parameters of a generative model.
- **Independent filtering / outlier replacement.** Our eligibility filter (mappable_bp ≥ mean_frag_len) is upstream. We do not down-weight outlier regions inside the EM.
- **GLM design matrices.** We fit a single intercept (λ_G), not a regression.

---

## 3. Implementation

### 3.1 Files to modify

| File | Change |
|---|---|
| [src/rigel/calibration/_em.py](src/rigel/calibration/_em.py) | Replace `_count_llr_poisson_ln` with `_count_llr_nb_ln`. Add `_estimate_phi_G_marginal` (golden-section MLE, Cox-Reid adjusted). Update `_m_step` to call φ_G estimator. Update `run_em` state to carry `phi_G` across iterations. |
| [src/rigel/calibration/_calibrate.py](src/rigel/calibration/_calibrate.py) | Pass `phi_G` through to the orchestrator. Surface it in `CalibrationResult` for QC reporting. |
| [src/rigel/calibration/_types.py](src/rigel/calibration/_types.py) | Add `phi_G: float` field to `CalibrationResult` and `CalibrationFit`. |
| [src/rigel/calibration/_config.py](src/rigel/calibration/_config.py) | Add `nb_density: bool = True` flag (default on; set False to fall back to Poisson exactly for regression testing). Add `phi_G_init: float = 0.0` and `phi_G_max: float = 50.0` bounds. |
| [src/rigel/calibration/_summary.py](src/rigel/calibration/_summary.py) | Emit `phi_G` and `nb_density` in `summary.json`. |
| [tests/test_calibration*.py](tests/) | Add NB-specific unit tests (see §5). Update golden outputs once. |
| [docs/METHODS.md](docs/METHODS.md) | Update §5 once empirically validated (see §6). |

### 3.2 New core kernels

```python
# src/rigel/calibration/_em.py

def _nb_logpmf(k: np.ndarray, mu: np.ndarray, phi: float) -> np.ndarray:
    """Negative Binomial log PMF, DESeq2 mean-dispersion parameterization.

    Var(k) = mu + phi * mu**2.
    As phi -> 0, this converges to Poisson(mu).

    Parameters mu and k may broadcast (e.g. (K,N) for quadrature).
    """
    if phi <= _EPS:
        # Exact Poisson limit; avoids the lgamma(1/phi) singularity.
        return k * np.log(np.maximum(mu, _EPS)) - mu - _vlgamma(k + 1.0)
    r = 1.0 / phi  # NB "size" parameter
    # log NB(k | mu, phi) = lgamma(k+r) - lgamma(r) - lgamma(k+1)
    #                       + r*log(r/(r+mu)) + k*log(mu/(r+mu))
    log_r_over_rmu = np.log(r) - np.log(r + mu)
    log_mu_over_rmu = np.log(np.maximum(mu, _EPS)) - np.log(r + mu)
    return (
        _vlgamma(k + r) - _vlgamma(r) - _vlgamma(k + 1.0)
        + r * log_r_over_rmu + k * log_mu_over_rmu
    )


def _count_llr_nb_ln(
    k_u: np.ndarray,
    E: np.ndarray,
    lam_G: float,
    mu_R: float,
    sigma_R: float,
    phi_G: float,
) -> np.ndarray:
    """Per-region count LLR with NB density: log P(k|G) − log P(k|R).

    P(k|G) = NegBin(k | λ_G·E, φ_G)
    P(k|R) = ∫ NegBin(k | (λ_G+μ)·E, φ_G) · LogNormal(μ; μ_R, σ_R²) dμ
             evaluated by 21-node Gauss-HermiteE quadrature on log μ.

    Reduces to _count_llr_poisson_ln when φ_G → 0.
    """
    k = np.asarray(k_u, dtype=np.float64)
    E = np.asarray(E, dtype=np.float64)

    # G class
    rate_G = lam_G * E
    ll_G = _nb_logpmf(k, rate_G, phi_G)

    # R class: quadrature over log μ
    mu_j = np.exp(mu_R + sigma_R * _GH_NODES)            # (K,)
    rate_R = (lam_G + mu_j[:, None]) * E[None, :]        # (K, N)
    log_p_R_nodes = _nb_logpmf(k[None, :], rate_R, phi_G)  # (K, N)
    log_p_R = _logsumexp(log_p_R_nodes + _GH_LOGW[:, None], axis=0)

    return ll_G - log_p_R
```

### 3.3 Dispersion estimator

We fit `φ_G` by a 1-D golden-section search on the **Cox-Reid-adjusted profile log-likelihood**, profiling out (λ_G, μ_R, σ_R) at each candidate φ. Because λ_G already has a closed-form weighted MLE (`Σγk / Σγ E`) and (μ_R, σ_R) are fit on the strand-anchored hard subset, the profile is cheap.

```python
# src/rigel/calibration/_em.py

def _estimate_phi_G_marginal(
    k_u: np.ndarray,
    E: np.ndarray,
    gamma: np.ndarray,
    lam_G: float,
    mu_R: float,
    sigma_R: float,
    phi_lo: float = 0.0,
    phi_hi: float = 50.0,
    n_iter: int = 60,
) -> float:
    """Golden-section MLE for shared NB dispersion φ_G.

    Maximizes the γ-weighted marginal log-likelihood of the two-component
    mixture (G = NB, R = NB-LogNormal compound), with Cox-Reid bias correction
    on the profiled-out λ_G mean parameter.

    Cox-Reid correction (single-mean case):
        ℓ_CR(φ) = ℓ(φ, λ̂_G(φ)) − ½ log Σ_i [γ_i E_i / (1 + φ·λ̂_G·E_i)]

    The denominator is the Fisher information of λ_G under NB at the MLE.
    For Poisson (φ=0) this reduces to ½ log Σ γ_i E_i, the standard correction.
    """
    # ... bracket-and-shrink golden-section over [phi_lo, phi_hi] on
    # log-likelihood evaluated by _count_llr_nb_ln components,
    # γ-weighted across regions.
```

Estimator runs once per M-step. Convergence target: |φ_new − φ_prev| < 1e-3 over the EM outer loop, in addition to the existing `|π_soft − prev_π_soft| < 1e-4` criterion.

### 3.4 EM loop changes

In [`run_em`](src/rigel/calibration/_em.py):

```python
# Initialize
phi_G = config.phi_G_init  # 0.0 → starts at Poisson, EM warms up dispersion

for it in range(config.max_iter):
    # E-step
    llr_count = _count_llr_nb_ln(k_u, E, lam_G, mu_R, sigma_R, phi_G)
    llr_strand = ...   # unchanged
    llr_fl = ...       # unchanged (FL channel has its own histogram model)
    log_odds = log_prior_odds + llr_count + llr_strand + llr_fl
    gamma = sigmoid(np.clip(log_odds, -500, 500))

    # M-step
    lam_G = weighted_poisson_mle(gamma, k_u, E)              # unchanged form
    mu_R, sigma_R = fit_lognormal_anchors(...)               # unchanged
    phi_G = _estimate_phi_G_marginal(                        # NEW
        k_u, E, gamma, lam_G, mu_R, sigma_R,
        phi_hi=config.phi_G_max,
    )
    pi, pi_soft = update_priors(gamma)
```

`lam_G` MLE is unchanged: under the NB-Gamma-Poisson hierarchy, the MLE of the mean given φ has the same `Σγk / ΣγE` form (the Gamma compound preserves the conditional mean).

### 3.5 Numerical considerations

- **`lgamma(k + 1/φ)` cost:** at typical k ≤ a few thousand and 1/φ on the order of 0.02–10, the lgamma function is well-behaved. We use the existing `_vlgamma` (vectorized scipy.gammaln).
- **Poisson limit:** branch on `φ ≤ _EPS` to use the closed-form Poisson PMF. The NB formula is numerically unstable as φ → 0 (1/φ blows up).
- **Large φ:** clamp at `phi_G_max=50` (corresponds to CV² of background rate = 50, which is already higher than expected even for severe hybrid capture). Keeps the search bounded; values pinning at the upper bound will be flagged in QC.
- **No autodiff needed:** golden section is gradient-free; the cost of one φ optimization is O(60 evaluations × 21 quadrature nodes × N regions) ≈ 1–2 s per iteration on the 700k-region human index. Acceptable.

---

## 4. Expected behavior and theoretical predictions

| Scenario | Expected φ_G | Expected effect on density LLR |
|---|---|---|
| Synthetic uniform λ_G (current tests) | ≈ 0 | Identical to Poisson; no regression |
| WGS-like genomic input | 0.1–1 | Small per-region down-weighting; tighter γ posteriors |
| MCTP hybrid capture (dna10m+) | 10–50 | Per-region density LLR shrinks 5–20×, letting strand channel dominate inside expressed loci |
| Ultra-clean dna00m pure-RNA | ≈ 0 (no gDNA to fit) | No density signal either way; γ collapses to π_soft. No change vs current behavior. |

The **decisive theoretical prediction**: across the VCaP titration, learned φ_G should correlate strongly with the empirical CV² of `k_u_i / E_i` over the high-eligibility intergenic subset. If the model is right, this correlation will be near 1.

---

## 5. Validation plan

### 5.1 Unit tests (new)

- `test_nb_logpmf_poisson_limit`: assert `_nb_logpmf(k, μ, 1e-12) ≈ poisson.logpmf(k, μ)` to 1e-9 across a grid of (k, μ).
- `test_nb_logpmf_against_scipy`: assert agreement with `scipy.stats.nbinom.logpmf` for φ ∈ {0.1, 1, 10}, k ∈ {0, 5, 50, 500}.
- `test_count_llr_nb_reduces_to_poisson`: assert `_count_llr_nb_ln(..., phi=0)` ≈ `_count_llr_poisson_ln(...)` to 1e-9.
- `test_estimate_phi_G_recovers_truth`: synthesize NB(μ, φ_true) data with φ_true ∈ {0.1, 1, 10}, fit, assert recovery within 10% RSE.
- `test_estimate_phi_G_cox_reid_unbiased`: simulate small N (50 regions) Poisson data, confirm CR-corrected estimator concentrates near 0 while uncorrected MLE biases positive.

### 5.2 Integration tests (existing, regenerate goldens)

- `test_golden_output.py`: with `nb_density=False`, assert bit-identical output to current code (regression guard).
- With `nb_density=True` (default), regenerate goldens once; confirm `phi_G` field appears in `summary.json` and is small for synthetic fixtures.

### 5.3 VCaP titration validation

Run calibration on all 8 libraries `mctp_vcap_rna20m_dna{00,01,02,05,10,20,40,80}m`, both with `nb_density=False` (current) and `nb_density=True`. Per the parent doc's §5.1 plan:

1. Plot fitted `φ_G` vs nominal gDNA spike (M reads). Expect monotone increase from ~0 (dna00m) to 10–50 (dna80m).
2. Compare per-region `|llr_count|` distributions. Expect 5–10× shrinkage at dna10m+, near-zero change at dna00m.
3. Compare final per-region γ to the strand-only baseline `γ_strand` derived from per-region Binomial posteriors. Pearson correlation should improve from current ~0.5 to >0.9.
4. Compare `quant.gdna_total / strand_inferred_gdna_total`. Expect drop from 1.1–1.5× to 1.0–1.1×.

### 5.4 No-regression on synthetic benchmarks

`scripts/benchmark/configs/locus_simple_*.yaml` (uniform-λ_G synthetic data). Expect:
- `φ_G` fitted ≈ 0 across all sweeps.
- mRNA / nRNA / gDNA relative errors within 1% of current.

---

## 6. Rollout

| Phase | Action | Gate |
|---|---|---|
| 1 | Implement `_nb_logpmf`, `_count_llr_nb_ln`, `_estimate_phi_G_marginal`. Pure additions, default `nb_density=False`. | All existing tests pass unchanged. New unit tests pass. |
| 2 | Run VCaP titration validation (§5.3). Document φ_G trajectory. | If §5.3 predictions hold, proceed. If not, root-cause before flipping default. |
| 3 | Flip `nb_density=True` as default. Regenerate goldens. | Synthetic benchmarks (§5.4) within 1% of baseline. |
| 4 | Update [docs/METHODS.md](docs/METHODS.md) §5: remove the historical `W = (2·SS−1)²` weighting description, document the NB density model. | Only after §5.3 confirms the NB model resolves the channel-fusion failure on hybrid capture. |
| 5 (later) | Add `--targets BED` mode for unstranded hybrid-capture support (separate plan, follows this work). | — |

---

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| `φ_G` over-inflates and **down-weights density even when density is correct**, causing strand-poor (low-SS) libraries to lose useful information. | Validate on dUTP synthetic and unstranded synthetic. If observed, add `--max_phi_G` user override. The golden-section bound + Cox-Reid correction should already prevent runaway. |
| LogNormal-NB compound integral becomes numerically unstable for very small `μ_j` or very large `k`. | Reuse the existing Gauss-HermiteE node grid which is already validated for the Poisson-LogNormal version. The NB PMF in log-space remains well-conditioned for k up to ~10⁵, more than enough for per-region counts. |
| `lgamma(k + 1/φ)` cost dominates EM time at very small φ_G (1/φ → ∞). | Branch to Poisson PMF when φ ≤ 1e-6. |
| Cox-Reid correction itself biased when N is small. | Not a concern — we have ~700k eligible regions. CR correction is asymptotically negligible at this N but harmless. |
| φ_G absorbs mappability/coverage-bias mis-specification not just hybrid capture, possibly hiding real bias problems. | Surface φ_G in `summary.json` as a QC metric. Document expected ranges per library type. |

---

## 8. References

- Love, M.I., Huber, W., Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550.
- Cox, D.R., Reid, N. (1987). *Parameter orthogonality and approximate conditional inference*. JRSS-B 49(1):1–39.
- McCarthy, D.J., Chen, Y., Smyth, G.K. (2012). *Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation*. Nucleic Acids Research 40(10):4288–4297. (Cox-Reid in edgeR.)
- Parent document: [channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md).
- Code: [src/rigel/calibration/_em.py](src/rigel/calibration/_em.py).
