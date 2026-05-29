# 03 — Inference

This document specifies the inference algorithm. Inputs are
**per-region (and per-boundary) sufficient statistics**. There is no
per-fragment iteration anywhere in the calibrator.

Every step is either:
- A weighted sum over $r$ (E-step soft allocation; every M-step except $\phi$, $\rho_d^{\text{BB}}$, $\rho_r^{\text{BB}}$),
- A 1-D Newton on a single global scalar ($\phi$, $\rho_d^{\text{BB}}$, $\rho_r^{\text{BB}}$),
- A closed-form Gamma posterior per region (exposure).

> **FL is not an inference channel.** Earlier drafts evaluated a
> per-region FL log-Bayes-factor and per-bin FL responsibility over
> a `(R, FL_BINS)` grid. Those terms have been removed — see
> [`00_overview.md`](00_overview.md) §1 and
> [`01_generative_model.md`](01_generative_model.md) §2. The
> downstream EM scorer continues to use FL per-fragment at full
> resolution.

## 1. Algorithm at a glance

```
Initialize: rho_0, phi, rho_d_BB, rho_r_BB, eps_s  (§5)
            kappa_d = 0.5 (fixed); kappa_rna from StrandModel (fixed)
            pi_g[r] = 0.5;  omega[r] = 1.0                         (per-region)

for iter in 1..max_outer:
    # ---- E-step: per-region soft allocation (§3) ----
    for each region/boundary r:                                    # vectorized over |S|
        compute strand log-Bayes-factor (BB-gDNA vs BB-RNA)
        compute count log-Bayes-factor (NB-gDNA-only vs NB-mixture)
        update pi_g[r] = sigmoid( logit(pi_g_prior[r]) + LLR_count[r] + LLR_strand[r] )
        M_g_unspl[r] = n_unspl[r] * pi_g[r]
        M_d_unspl[r] = n_unspl[r] - M_g_unspl[r]
        # Spliced is deterministic RNA (with eps_s failsafe correction)
        M_g_cont[r] = M_g_unspl[r] + eps_s * n_spliced[r]
        M_d_cont[r] = M_d_unspl[r] + (1 - eps_s) * n_spliced[r]

    # ---- Exposure update: closed-form Gamma posterior (§4) ----
    M_g_tot = M_g_cont + 0.5 * (M_g_L + M_g_R)
    alpha_post = 1.0/phi + M_g_tot
    beta_post  = 1.0/phi + rho_0 * L_eff
    omega          = alpha_post / beta_post                        # posterior mean
    log_omega_var  = 1.0 / alpha_post                              # delta-method log-variance

    # ---- M-step (§5) ----
    rho_0       <- sum(M_g_tot) / sum(omega * L_eff)
    eps_s       <- update via Beta(1,1) + soft-allocated artifacts
    phi         <- 1-D Newton on NB log-likelihood gradient
    rho_d_BB    <- 1-D Newton on BB log-likelihood gradient (kappa_d = 0.5 fixed)
    rho_r_BB    <- 1-D Newton on BB log-likelihood gradient (kappa_rna fixed)
    pi_g_prior  <- omega * rho_0 * L_eff / (omega * rho_0 * L_eff + max(N_r - omega * rho_0 * L_eff, eps))

    # ---- Convergence ----
    delta = max(|M_g_tot - M_g_tot_prev| / (M_g_tot_prev + 1))
    if delta < mass_rel_tol: break
```

Total LoC estimate: ~300 numpy + ~75 scipy for the three 1-D Newton routines.

## 2. Inputs

Vectorized substrate arrays. All shapes named relative to the
"region-or-boundary set" $S$. The calibrator runs **three times** —
once with $S$ = contained regions, once with $S$ = left boundaries,
once with $S$ = right boundaries — sharing the global hyperparameters
through a single outer-loop schedule. Concretely:

| Array | Shape | Dtype | Meaning |
|---|---|---|---|
| `n_unspl[s]` | `(|S|,)` | `int64` | unspliced count |
| `n_spliced[s]` | `(|S|,)` | `int64` | spliced count |
| `k_plus[s]` | `(|S|,)` | `int64` | sense-strand unspliced count |
| `L_eff[s]` | `(|S|,)` | `float64` | effective length |
| `kappa_rna[s]` | `(|S|,)` | `float64` | per-region RNA sense-strand mean |

The three pre-allocated index sets share the underlying region
identity so that boundary outputs can be split back into the
neighbouring regions for exposure aggregation (§4).

## 3. E-step (per region, vectorized over `|S|`)

### 3.1 Channel log-Bayes-factor: count

Compare the gDNA-only mean $\mu_r^{(g)} = \omega_r \rho_0 L_r^{\text{eff}}$
against the mixture mean $\mu_r^{(g)} + \mu_r^{(d)}$ under the NB
likelihood at the observed $n_r^{\text{u}}$:

```python
mu_g = omega * rho_0 * L_eff
mu_mix = mu_g + M_d_unspl_prev                # M_d from previous outer iter (start: 0)
LL_count_g = nbinom.logpmf(n_unspl, n=1.0/phi, p=(1.0/phi)/(1.0/phi + mu_g))
LL_count_mix = nbinom.logpmf(n_unspl, n=1.0/phi, p=(1.0/phi)/(1.0/phi + mu_mix))
LLR_count = LL_count_g - LL_count_mix
```

Under the gDNA-only hypothesis the observed count is consistent with
the exposure-driven mean; under the mixture hypothesis the count is
explained partly by inferred RNA mass. The log-ratio is positive when
the count matches the gDNA-only expectation and negative when the
count exceeds it by a margin the gDNA model cannot explain.

### 3.2 Channel log-Bayes-factor: strand (BB vs BB)

Both gDNA and RNA strand counts are Beta-Binomial. The gDNA mean is
fixed at $\kappa_d = 0.5$ (biology); the RNA mean $\kappa_r^{\text{rna}}$
is pre-computed by `StrandModel` from spliced uniques during the BAM
scan. Only the two dispersion parameters $\rho_d^{\text{BB}}$ and
$\rho_r^{\text{BB}}$ are fitted by the calibrator.

For each region $s$, the strand log-evidence ratio against the
"all-gDNA" vs "all-RNA" hypothesis on the unspliced count:

```python
# Under gDNA: K^+ ~ BetaBinom(n_unspl, kappa_d=0.5, rho_d_BB)
alpha_d = 0.5 * (1.0 - rho_d_BB) / rho_d_BB
beta_d  = alpha_d                                                  # symmetric about 0.5
LL_strand_g = betabinom.logpmf(k_plus, n_unspl, alpha_d, beta_d)
# Under RNA: K^+ ~ BetaBinom(n_unspl, kappa_rna, rho_r_BB)
alpha_r = kappa_rna * (1.0 - rho_r_BB) / rho_r_BB
beta_r  = (1.0 - kappa_rna) * (1.0 - rho_r_BB) / rho_r_BB
LL_strand_d = betabinom.logpmf(k_plus, n_unspl, alpha_r, beta_r)
LLR_strand = LL_strand_g - LL_strand_d
```

For unstranded libraries `kappa_rna ≈ 0.5`, so the BB-RNA and
BB-gDNA likelihoods are similar (both centred at 0.5) and
`LLR_strand ≈ 0` — the strand channel contributes nothing to the
discrimination, by construction, no special case. For perfectly
stranded libraries `kappa_rna ≈ 1.0` (or 0.0), and the BB-RNA
likelihood is sharply asymmetric, making strand a strong discriminator
— again automatic. The 126/26 paralog asymmetry produces near-zero
`LLR_strand` because BB-gDNA at $\rho_d^{\text{BB}} > 0$ has enough
variance mass at moderate asymmetry to absorb it.

### 3.3 Combining channels into per-region $\pi_r^{(g)}$

```python
logit_pi = logit(pi_g_prior) + LLR_count + LLR_strand
pi_g     = expit(logit_pi)
```

Stable in log-space via `scipy.special.expit`.

### 3.4 Soft-allocated masses (unspliced)

```python
M_g_unspl = n_unspl * pi_g                                         # (|S|,)
M_d_unspl = n_unspl - M_g_unspl
```

### 3.5 Soft-allocated strand-plus counts (for $\rho_r^{\text{BB}}$ M-step)

The per-region observed sense count $k_r^{(+)}$ is the sum of a gDNA
component (mean $0.5 \cdot M_r^{(g, \text{unspl})}$) and an RNA
component (mean $\kappa_r^{\text{rna}} \cdot M_r^{(d, \text{unspl})}$).
For the RNA-side $\rho_r^{\text{BB}}$ M-step we form the
soft-allocated RNA sense count

```python
k_plus_d_hat = np.maximum(k_plus - 0.5 * M_g_unspl, 0.0)           # RNA-attributable +
k_plus_g_hat = np.maximum(k_plus - kappa_rna * M_d_unspl, 0.0)     # gDNA-attributable +
```

The two soft-allocated counts feed the two BB Newton M-steps (§5.2
and 5.3) directly. No method-of-moments mean estimation is needed
because both means ($\kappa_d = 0.5$, $\kappa_r^{\text{rna}}$) are
fixed inputs.

### 3.6 Total masses (G2 outputs)

```python
M_g_cont = M_g_unspl + eps_s * n_spliced
M_d_cont = M_d_unspl + (1.0 - eps_s) * n_spliced
```

Same machinery on the boundary substrate gives `M_g_L, M_d_L, M_g_R, M_d_R` (G3 outputs).

## 4. Per-region exposure: closed-form Gamma posterior (G4)

For each region $r$:

```python
M_g_tot       = M_g_cont + 0.5 * (M_g_L + M_g_R)                   # half-split boundary attribution
alpha_post    = 1.0 / phi + M_g_tot
beta_post     = 1.0 / phi + rho_0 * L_eff
omega         = alpha_post / beta_post                              # E[omega | data]
log_omega_var = 1.0 / alpha_post                                    # Var(log omega | data) by delta method
```

This is `O(R)` and fully vectorized. No iteration.

## 5. M-step (library hyperparameters)

### 5.1 $\rho_0$ — closed form

```python
rho_0 = M_g_tot.sum() / (omega * L_eff).sum()
```

### 5.2 $\rho_d^{\text{BB}}$ — 1-D Newton (gDNA strand dispersion)

$\kappa_d = 0.5$ is fixed (not estimated); only $\rho_d^{\text{BB}}$
is fitted. The BetaBinom log-likelihood over all regions, as a
function of $\rho_d^{\text{BB}}$ alone:

```python
def neg_log_lik_bb_gdna(rho_bb):
    alpha = 0.5 * (1.0 - rho_bb) / rho_bb
    beta  = alpha                                                  # symmetric about 0.5
    return -betabinom.logpmf(
        np.round(k_plus_g_hat).astype(int),
        np.round(M_g_unspl).astype(int),
        alpha, beta,
    ).sum()
```

Single global scalar; `scipy.optimize.minimize_scalar` on the bounded
interval $(10^{-6}, 1 - 10^{-6})$. Moment estimator warm start:

$$
\hat{\rho}_d^{\text{BB}, \text{init}} = \mathrm{clip}\!\left(\frac{\widehat{\mathrm{Var}}\bigl(k_r^{(g,+)} / M_r^{(g, \text{unspl})}\bigr) - \overline{0.25 / M_r^{(g, \text{unspl})}}}{0.25 - \overline{0.25 / M_r^{(g, \text{unspl})}}}, \;10^{-6},\; 1 - 10^{-6}\right).
$$

### 5.3 $\rho_r^{\text{BB}}$ — 1-D Newton (RNA strand dispersion)

$\kappa_r^{\text{rna}}$ is the **pre-computed** per-region RNA strand
mean from `StrandModel` (not estimated by the calibrator). Only
$\rho_r^{\text{BB}}$ is fitted, again as a single global scalar:

```python
def neg_log_lik_bb_rna(rho_bb):
    alpha = kappa_rna * (1.0 - rho_bb) / rho_bb
    beta  = (1.0 - kappa_rna) * (1.0 - rho_bb) / rho_bb
    return -betabinom.logpmf(
        np.round(k_plus_d_hat).astype(int),
        np.round(M_d_unspl).astype(int),
        alpha, beta,
    ).sum()
```

Same `scipy.optimize.minimize_scalar` pattern, same bounds. Moment
estimator warm start uses per-region $\kappa_r^{\text{rna}}(1 - \kappa_r^{\text{rna}})$ in place of $0.25$.

### 5.4 $\phi$ — 1-D Newton

NB log-likelihood gradient on the unspliced counts, treating the
mixture mean $\mu_r^{(g)} + \mu_r^{(d)}$ as known (plug in current
$\omega_r \rho_0 L_r^{\text{eff}}$ for $\mu_r^{(g)}$ and $M_d^{\text{unspl}}_r$ for $\mu_r^{(d)}$):

```python
def neg_log_lik_nb(phi):
    mu = omega * rho_0 * L_eff + M_d_unspl
    return -nbinom.logpmf(n_unspl, n=1.0/phi, p=(1.0/phi)/(1.0/phi + mu)).sum()
```

Same approach as §5.3 with `scipy.optimize.minimize_scalar` bounded on
$(10^{-6}, 10^{2})$. Moment estimator warm start: $\hat{\phi}_{\text{init}} = \max(\widehat{\mathrm{Var}(n)}/\bar{\mu}^2 - 1/\bar{\mu}, 10^{-6})$.

### 5.5 $\epsilon_s$ — closed form (failsafe)

The upstream BAM scanner already removes most splice artifacts before
they reach the substrate. $\epsilon_s$ is a failsafe that quantifies
any residual mis-classifications rather than letting them silently
degrade the model.

```python
eps_s_num = 1.0 + (M_g_unspl / (M_g_unspl + M_d_unspl) * n_spliced).sum()
eps_s_den = 1.0 + n_spliced.sum()
eps_s = eps_s_num / eps_s_den                                       # Beta(1,1) prior
```

In practice $\epsilon_s$ converges very close to its Beta(1,1) prior
mean unless the upstream artifact filter has a real gap, in which case
the value rises above $10^{-3}$ and surfaces in the logged hyperparameters.

### 5.6 $\pi_r^{(g, \text{prior})}$ — next iteration

```python
mu_g    = omega * rho_0 * L_eff
mu_d    = np.maximum(n_unspl - mu_g, _PI_FLOOR)
pi_g_prior_next = mu_g / (mu_g + mu_d)
pi_g_prior_next = np.clip(pi_g_prior_next, _PI_CLIP, 1.0 - _PI_CLIP)
```

The `np.maximum(..., eps)` is the only soft-guard in the algorithm.
With $\mu_g \to N_r$ the prior tends to 1 (region is all-gDNA); with
$\mu_g \to 0$ it tends to 0 (region is all-RNA). No qualitative cliff.

## 6. Boundaries: same code, parallel substrate

The boundary E-step is the same code with boundary sufficient
statistics. After all three E-steps (contained, left, right) have run
within an outer iteration, the boundary masses contribute to
$M_r^{(g, \text{tot})}$ via the half-split in §4. The boundary
populations have their own per-boundary $\pi_b^{(g, \text{prior})}$
that converges independently of the region prior.

## 7. Initialization

All data-driven, no magic:

| Param | Initial value |
|---|---|
| $\pi_r^{(g)}, \pi_b^{(g)}$ | 0.5 |
| $\omega_r$ | 1.0 |
| $\rho_0$ | $0.5 \cdot \sum_r n_r^{\text{u}} / \sum_r L_r^{\text{eff}}$ (half-and-half assumption) |
| $\phi$ | $\max\bigl(\widehat{\mathrm{Var}}(n^{\text{u}})/\bar{n^{\text{u}}}^2 - 1/\bar{n^{\text{u}}},\; 10^{-3}\bigr)$ |
| $\kappa_d$ | **0.5 (fixed; biology)** |
| $\rho_d^{\text{BB}}$ | 0.01 |
| $\rho_r^{\text{BB}}$ | 0.01 |
| $\epsilon_s$ | $10^{-3}$ |
| $\kappa_r^{\text{rna}}$ | **from `StrandModel.p_r1_sense` (fixed; pre-computed during BAM scan)** |

## 8. Numerical constants (exhaustive list)

Eight numerics total. Each is a Bayesian prior of unit strength, a
machine-precision floor, or a documented iteration cap.

| Name | Value | Role |
|---|---|---|
| `_PHI_FLOOR` | `1e-6` | Floor on $\phi$ Newton (prevents $1/\phi$ blowup). |
| `_BB_FLOOR` | `1e-6` | Floor on $\rho_d^{\text{BB}}$ and $\rho_r^{\text{BB}}$ Newton. |
| `_PI_FLOOR` | `1e-9` | Floor inside `np.maximum` for $\pi_r^{(g)}$ update. |
| `_PI_CLIP` | `1e-6` | Keeps $\pi_r^{(g)} \in (\epsilon, 1-\epsilon)$. |
| `_TOL_MASS` | `1e-4` | Relative max-region mass-change convergence. |
| `_MAX_OUTER` | `25` | Outer cap; converges typically in 5–10. |
| `_LOG_EPS` | `np.finfo(np.float64).tiny` | Additive floor inside `np.log()`. |
| `_K_RNA_FALLBACK` | `0.5` | Used when `StrandModel` has no estimate (no spliced uniques on region's strand). |

Count: 8. Each derivable; none flip qualitative behavior.

## 9. Convergence

Standard EM theory: each E-step and each M-step weakly increases the
observed-data log-likelihood. Convergence is guaranteed; in practice
3–10 outer iterations.

The only failure mode that could violate this is an implementation
bug in the gradient calculations for the two 1-D Newton steps. The
mass-change diagnostic in §1 is the runtime sentinel: if it
*increases* between iterations, raise `CalibrationConvergenceError`
immediately rather than silently continuing.

## 10. Complexity

Per outer iteration (one substrate set, e.g. contained regions):
- E-step soft allocation: `O(|S|)` per channel.
- Count + strand log-evidence: `O(|S|)` scipy calls.
- M-step: closed forms `O(|S|)` + three 1-D Newton (constant time).

For 3 substrate sets × 10 outer iterations × $|S| = 2 \times 10^5$:
about 3 × 10 × (~30 ms) ≈ 1 s wall-clock. Negligible relative to BAM
scan.

## 11. Robustness properties (by construction)

1. **Bounded posteriors.** $\pi_r^{(g)} \in (0, 1)$ via `expit`. No clipping needed.
2. **Bounded masses.** $M_r^{(g, \text{cont})}, M_r^{(d, \text{cont})} \in [0, n_r^{\text{u}} + n_r^{\text{s}}]$.
3. **Bounded exposure.** $\hat{\omega}_r$ strictly positive and finite for all $M_r^{(g, \text{tot})} \geq 0$. The $1/\phi$ floor in numerator and denominator prevents $0/0$ in empty-region cases.
4. **Monotone EM theory.** Mass-change diagnostic must decrease; runtime sentinel catches violations.
5. **Identifiability on degenerate cases.** Region with zero fragments: $M_r = 0$, $\hat{\omega}_r = 1$ (prior mean), $\log\sigma_r^2 = \phi$ (prior variance). Downstream consumer down-weights via inverse variance.
6. **No FL channel.** The previously-cited identifiability hazard of coupling FL pmf estimation with mass deconvolution is structurally absent.

## 12. Comparison to previous drafts

| Aspect | Previous draft (per-fragment Bayes) | This draft (per-region sufficient stats) |
|---|---|---|
| Inference level | Per-fragment posterior | Per-region sufficient-statistic scalar algebra |
| Streaming pass over fragments | Required (additional) | Not required; native already aggregates |
| Count overdispersion model | Poisson (with $\tau^2$ exposure variance separately) | NB unified with Gamma exposure prior via single $\phi$ |
| Strand overdispersion | Binomial | Beta-Binomial on **both** gDNA ($\kappa_d = 0.5$ fixed) and RNA ($\kappa_r^{\text{rna}}$ pre-computed) |
| Splice channel | $\bar{q}$ rate hyperparameter | Removed; spliced fragments are deterministic RNA (with $\epsilon_s$ failsafe) |
| Per-region E-step | $O(F_r)$ per region | $O(\text{FL\_BINS})$ per region |
| Total time complexity | $O(F + R)$ per iter | $O(R \cdot \text{FL\_BINS})$ per iter |
| Memory | per-fragment table ~450 MB | per-region FL histogram ~75 MB |
| Newton solvers | None per region; 1 global ($\phi$) | None per region; 3 global ($\phi$, $\rho_d^{\text{BB}}$, $\rho_r^{\text{BB}}$) |
| Total LoC | ~250 | ~300 |
