# Strand calibration redesign: implementation plan

**Date:** 2026-04-23
**Status:** Implementation plan. Ready to code.
**Parent:** [strand_channel_theory_and_redesign.md](strand_channel_theory_and_redesign.md)
**Scope:** Replace the strand likelihood in
[`src/rigel/calibration/_em.py`](../../src/rigel/calibration/_em.py)
with the $\rho$-integrated Beta-Binomial model, separate
$\kappa_G$ / $\kappa_R$, remove the library-level z-gate, and
surface a per-region strand-only posterior.

---

## 1. Scope and non-goals

**In scope:**

- New strand LLR kernel `_strand_llr_bb_rho` (BB under G, $\rho$-integrated BB under R).
- New dispersion estimators `_estimate_kappa_G_soft`, `_estimate_kappa_R_soft` (γ-weighted soft MLE per class).
- Integrate into the existing `run_em` EM loop — strand LLR always on, no z-gate.
- Expose per-region $\tilde\gamma^{\text{strand}}_i$ on `EMFit` / `CalibrationResult` as a diagnostic field.
- Simplify: delete the z-gate code path, the shared-κ estimator, the Binomial `_strand_llr`, and the `strand_llr_mode` toggle. One model, cleanly replacing the old ones.

**Out of scope (deferred):**

- Density channel redesign ([density_nb_model_plan.md](density_nb_model_plan.md)).
- Beta-distributed $\rho$ prior (future upgrade; start with Uniform).
- Fragment-length channel changes.
- Any fusion-layer veto / rail / gate.

**Deliberately removed:**

- `_strand_llr` (Binomial).
- `_strand_llr_betabinom` (shared-κ BB without ρ integration).
- `_estimate_kappa_marginal` (shared-κ mixture-marginal estimator).
- `_aggregate_strand_z` (library-level z-gate).
- `strand_used`, `strand_z` fields on `EMFit` / `CalibrationResult` (always true, always 0 — meaningless).
- `strand_llr_mode` field on `CalibrationConfig`. One model. No mode.
- `strand_z_threshold` argument on `run_em`.

---

## 2. Mathematical reference (what the code computes)

### 2.1 Per-region variables (from `_stats.compute_region_stats`)

- $n_i$ = `n_unspliced[i]` = number of unspliced fragments in region $i$.
- $k_i$ = sense-strand count, derived from `n_unspliced`, `n_pos`, `tx_strand` (R1-antisense convention, existing `_compute_k_sense_valid` helper).
- Valid region: $n_i \ge 1$ and $\texttt{tx\_strand}[i] \ne 0$.
  (We relax from $n_i \ge 2$ to $n_i \ge 1$: a single oriented read contributes a real, bounded LLR under BB, and discarding it is a legacy artifact.)

### 2.2 G-class log-PMF

$$\log P(k \mid n, G, \kappa_G) \;=\; \log \binom{n}{k} + \log B(k + \tfrac{\kappa_G}{2},\, n - k + \tfrac{\kappa_G}{2}) - \log B(\tfrac{\kappa_G}{2}, \tfrac{\kappa_G}{2})$$

with $\log B(a, b) = \Gamma(a) + \Gamma(b) - \Gamma(a + b)$ (use `scipy.special.gammaln`). The $\log\binom{n}{k}$ cancels in the LLR but is kept in the per-region log-likelihoods for the M-step objective.

### 2.3 R-class log-PMF (ρ-integrated)

With $p(\rho) = \rho \cdot \mathrm{SS} + (1 - \rho)/2$, $\rho \sim \mathrm{Uniform}(0, 1)$:

$$P(k \mid n, R, \kappa_R) \;=\; \int_0^1 \mathrm{BB}\!\left(k; n, \kappa_R p(\rho), \kappa_R (1 - p(\rho))\right) d\rho$$

Evaluated by 17-node Gauss–Legendre quadrature on $[0, 1]$:

$$\log P(k \mid n, R, \kappa_R) \;\approx\; \mathrm{logsumexp}_{j=1}^{17}\!\left[\,\log w_j + \log \mathrm{BB}(k; n, \kappa_R p(\rho_j), \kappa_R (1 - p(\rho_j)))\,\right]$$

where $\{\rho_j, w_j\}$ come from `np.polynomial.legendre.leggauss(17)` mapped from $[-1, 1]$ to $[0, 1]$ (so weights are pre-multiplied by $1/2$; sum of weights = 1).

Nodes and log-weights are module-level constants, computed once at import.

### 2.4 Strand LLR

$$\mathrm{LLR}^{\text{strand}}_i \;=\; \log P(k_i \mid n_i, G, \kappa_G) - \log P(k_i \mid n_i, R, \kappa_R)$$

For invalid regions ($n_i = 0$ or $\texttt{tx\_strand}_i = 0$): LLR = 0.

### 2.5 Strand-only posterior (diagnostic)

$$\tilde\gamma^{\text{strand}}_i \;=\; \mathrm{sigmoid}\!\left(\log\tfrac{\pi_0}{1-\pi_0} + \mathrm{LLR}^{\text{strand}}_i\right)$$

with $\pi_0 = 0.5$ (uninformative). Computed once per E-step on the current LLR vector; essentially free.

### 2.6 κ M-step

$\kappa_G$: maximise $\sum_i \gamma_i \cdot \log P(k_i \mid n_i, G, \kappa)$ over $\kappa \in [\kappa_\min, \kappa_\max]$.
$\kappa_R$: maximise $\sum_i (1 - \gamma_i) \cdot \log P(k_i \mid n_i, R, \kappa)$ over the same bracket.

Default bracket: $[\kappa_\min, \kappa_\max] = [0.1, 1000]$. Chosen so that:

- Lower bound $\kappa = 0.1$: near-U-shaped Beta(0.05, 0.05), extreme overdispersion, rate-floor protection against degenerate data.
- Upper bound $\kappa = 1000$: effectively Binomial (the BB and Binomial variances differ by < 0.1% at $\kappa = 1000$). No reason to go higher.

Each fit uses the existing `_golden_section_max` helper. γ-weighted log-likelihood evaluations dominate cost: each $\kappa$ candidate touches all $N$ regions; R-class touches $17 N$ BB evaluations per candidate. Golden section converges in ~30 candidates on $[0.1, 1000]$ at $10^{-3}$ relative tolerance, so one κ_R M-step = $\sim 510 N$ BB PMF evaluations. At $N = 700\text{k}$ this is $\sim 3.6 \times 10^8$ `gammaln` calls per κ_R M-step. Acceptable — use vectorized `scipy.special.gammaln`, avoid Python loops over regions.

---

## 3. Files to modify

| File | Change |
|---|---|
| [`src/rigel/calibration/_em.py`](../../src/rigel/calibration/_em.py) | Delete `_strand_llr`, `_strand_llr_betabinom`, `_estimate_kappa_marginal`, `_aggregate_strand_z`, `_betaln_scalar`, `_betaln_vec` (inline them or keep `_betaln_vec` as helper). Add `_bb_logpmf`, `_strand_llr_bb_rho`, `_estimate_kappa_G_soft`, `_estimate_kappa_R_soft`, `_strand_posterior`. Refactor `_e_step` / `_m_step` / `run_em`. |
| [`src/rigel/calibration/_calibrate.py`](../../src/rigel/calibration/_calibrate.py) | Drop `strand_llr_mode` argument; drop `strand_used` / `strand_z` propagation; add `region_gamma_strand` propagation; add `kappa_G`, `kappa_R`. |
| [`src/rigel/calibration/_result.py`](../../src/rigel/calibration/_result.py) | Drop `strand_used`, `strand_z`, `strand_llr_mode`, `kappa` fields. Add `kappa_G`, `kappa_R`, `region_gamma_strand` (np.ndarray or None). Update `to_summary_dict`. |
| [`src/rigel/config.py`](../../src/rigel/config.py) | Drop `strand_llr_mode` field on `CalibrationConfig`. Field removal is a breaking change to the YAML config; acceptable. |
| [`scripts/profiling/profiler.py`](../../scripts/profiling/profiler.py) | Drop `strand_llr_mode=` kwarg forwarding. |
| [`scripts/benchmark/betabinom_pilot.py`](../../scripts/benchmark/betabinom_pilot.py) | This script's reason for existence (Binomial vs BetaBinomial comparison) is obsolete. Either delete it outright or replace with a new strand-model validation script (§6 below). |
| [`tests/test_betabinom_strand.py`](../../tests/test_betabinom_strand.py) | Rewrite: tests the new BB+ρ model. See §5. |
| [`tests/test_calibration*.py`](../../tests/) | Audit all existing calibration tests. Many will change numerics (new model, new LLR magnitudes). Update goldens deliberately. |
| [`tests/test_golden_output.py`](../../tests/test_golden_output.py) | Regenerate goldens once after the model is validated. |
| [`docs/METHODS.md`](../../docs/METHODS.md) | Update §5 (calibration) to describe the new strand model. Deferred until validation completes. |

---

## 4. Code sketch

### 4.1 Module-level constants

```python
# src/rigel/calibration/_em.py, near existing _GH_NODES setup

# 17-node Gauss–Legendre on [0, 1] for the ρ integral.
# np.polynomial.legendre.leggauss returns nodes x_j ∈ [-1, 1] with weights
# summing to 2 under the unit-weight measure.  We remap to ρ ∈ [0, 1]:
#   ρ_j = (x_j + 1) / 2,   w'_j = w_j / 2,   Σ w'_j = 1.
_RHO_NODES_X, _RHO_WEIGHTS_X = np.polynomial.legendre.leggauss(17)
_RHO_NODES = 0.5 * (_RHO_NODES_X + 1.0)          # shape (17,)
_RHO_LOGW = np.log(0.5 * _RHO_WEIGHTS_X)         # shape (17,)
```

### 4.2 Beta-Binomial log-PMF (numerically robust)

```python
def _bb_logpmf(
    k: np.ndarray,
    n: np.ndarray,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    """Log PMF of BetaBinomial(k; n, α, β).

    log BB = [lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)]   # C(n,k)
           + [lgamma(k+α) + lgamma(n-k+β) - lgamma(n+α+β)] # log B(k+α, n-k+β)
           - [lgamma(α)   + lgamma(β)     - lgamma(α+β)]   # log B(α, β)

    All inputs broadcastable; returns ndarray of broadcast shape.
    Handles k=0 and k=n cleanly (lgamma(0) is not invoked; lgamma(α) with
    α>0 is always finite).  Caller must ensure α > 0 and β > 0.
    """
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    log_binom = _vlgamma(n + 1.0) - _vlgamma(k + 1.0) - _vlgamma(n - k + 1.0)
    log_B_num = _vlgamma(k + alpha) + _vlgamma(n - k + beta) - _vlgamma(n + alpha + beta)
    log_B_den = _vlgamma(alpha) + _vlgamma(beta) - _vlgamma(alpha + beta)
    return log_binom + log_B_num - log_B_den
```

Shapes: for the G-class call, `alpha = beta = κ_G/2` (scalars) broadcast against `k, n` of shape `(N,)`. For the R-class call, `alpha, beta` have shape `(17, 1)` (one per ρ node) and `k, n` have shape `(1, N)`, so the broadcast result is `(17, N)`.

### 4.3 G and R log-likelihoods

```python
def _bb_g_loglik(k: np.ndarray, n: np.ndarray, kappa_G: float) -> np.ndarray:
    """Per-region log P(k | n, G, κ_G).  Shape (N,)."""
    alpha = beta = 0.5 * kappa_G
    return _bb_logpmf(k, n, alpha, beta)


def _bb_r_loglik(
    k: np.ndarray,
    n: np.ndarray,
    kappa_R: float,
    strand_specificity: float,
) -> np.ndarray:
    """Per-region log P(k | n, R, κ_R) under ρ ~ Uniform(0,1).  Shape (N,)."""
    ss = float(strand_specificity)
    p_rho = _RHO_NODES * ss + (1.0 - _RHO_NODES) * 0.5  # shape (17,)
    alpha = (kappa_R * p_rho)[:, None]                  # shape (17, 1)
    beta = (kappa_R * (1.0 - p_rho))[:, None]           # shape (17, 1)
    # Evaluate BB log-PMF at all 17 ρ nodes simultaneously.
    log_p_nodes = _bb_logpmf(k[None, :], n[None, :], alpha, beta)  # (17, N)
    # Integrate via logsumexp with pre-log Gauss-Legendre weights.
    return _logsumexp(log_p_nodes + _RHO_LOGW[:, None], axis=0)
```

### 4.4 Strand LLR and strand-only posterior

```python
def _strand_llr_bb_rho(
    stats: dict[str, np.ndarray],
    strand_specificity: float,
    kappa_G: float,
    kappa_R: float,
) -> np.ndarray:
    """Per-region strand LLR: log P(k|n,G,κ_G) − log P(k|n,R,κ_R).

    Zero for invalid regions (n < 1 or tx_strand == 0).  No library-level
    gate; SS near 0.5 causes the R-class integrand to collapse onto the
    G-class and the LLR vanishes naturally.
    """
    n_regions = stats["n_unspliced"].shape[0]
    valid, k_valid, n_valid = _compute_k_sense_valid(stats)  # existing helper,
                                                             # relaxed to n ≥ 1
    llr = np.zeros(n_regions, dtype=np.float64)
    if k_valid.size == 0:
        return llr
    ll_g = _bb_g_loglik(k_valid, n_valid, kappa_G)
    ll_r = _bb_r_loglik(k_valid, n_valid, kappa_R, strand_specificity)
    llr[valid] = ll_g - ll_r
    return llr


def _strand_posterior(llr: np.ndarray, pi_0: float = 0.5) -> np.ndarray:
    """Strand-only γ̃_i = sigmoid(logit(π_0) + LLR_strand_i).  Diagnostic."""
    pi_safe = float(np.clip(pi_0, _EPS, 1.0 - _EPS))
    logit = math.log(pi_safe / (1.0 - pi_safe))
    return 1.0 / (1.0 + np.exp(-np.clip(logit + llr, -500.0, 500.0)))
```

### 4.5 κ estimators (γ-weighted, soft, per class)

```python
def _estimate_kappa_G_soft(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    kappa_lo: float = 0.1,
    kappa_hi: float = 1000.0,
) -> float:
    """γ-weighted soft MLE for κ_G on the G-class BB.

    Maximises Σ_i γ_i · log BB(k_i; n_i, κ/2, κ/2) over κ ∈ [lo, hi].

    Returns the current midpoint if the valid region count is < 2 (not
    enough data to fit; caller uses the previous iteration's κ).
    """
    valid, k, n = _compute_k_sense_valid(stats)
    if k.size < 2:
        return math.sqrt(kappa_lo * kappa_hi)  # geometric midpoint; harmless
    w = gamma[valid].astype(np.float64)
    if float(w.sum()) <= 0.0:
        return math.sqrt(kappa_lo * kappa_hi)

    def obj(kappa: float) -> float:
        ll = _bb_g_loglik(k, n, kappa)
        return float(np.sum(w * ll))

    return _golden_section_max(obj, kappa_lo, kappa_hi, tol=1e-3)


def _estimate_kappa_R_soft(
    stats: dict[str, np.ndarray],
    gamma: np.ndarray,
    strand_specificity: float,
    kappa_lo: float = 0.1,
    kappa_hi: float = 1000.0,
) -> float:
    """(1−γ)-weighted soft MLE for κ_R on the ρ-integrated R-class BB."""
    valid, k, n = _compute_k_sense_valid(stats)
    if k.size < 2:
        return math.sqrt(kappa_lo * kappa_hi)
    w = (1.0 - gamma[valid]).astype(np.float64)
    if float(w.sum()) <= 0.0:
        return math.sqrt(kappa_lo * kappa_hi)

    def obj(kappa: float) -> float:
        ll = _bb_r_loglik(k, n, kappa, strand_specificity)
        return float(np.sum(w * ll))

    return _golden_section_max(obj, kappa_lo, kappa_hi, tol=1e-3)
```

### 4.6 E-step simplification

```python
def _e_step(
    stats,
    eligible,
    lam_G,
    mu_R,
    sigma_R,
    pi_soft,
    strand_specificity,
    kappa_G,
    kappa_R,
    llr_fl,
    *,
    return_strand_posterior=False,
):
    """Compute γ_i from fused LLR.  No strand_enabled flag, no mode."""
    n = stats["n_unspliced"].shape[0]
    n_s = stats["n_spliced"]
    gamma = np.full(n, pi_soft, dtype=np.float64)

    hard_expr = (n_s > 0) & eligible
    gamma[hard_expr] = 0.0

    soft = eligible & ~hard_expr
    llr_strand_full = _strand_llr_bb_rho(stats, strand_specificity, kappa_G, kappa_R)

    if not soft.any():
        if return_strand_posterior:
            return gamma, _strand_posterior(llr_strand_full)
        return gamma

    pi_safe = float(np.clip(pi_soft, _EPS, 1.0 - _EPS))
    log_prior_odds = math.log(pi_safe / (1.0 - pi_safe))

    llr_count = _count_llr_poisson_ln(
        stats["n_unspliced"][soft],
        stats["mappable_bp"][soft],
        lam_G, mu_R, sigma_R,
    )
    llr_strand = llr_strand_full[soft]
    llr_fl_soft = llr_fl[soft] if llr_fl is not None else 0.0

    log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft
    log_odds = np.clip(log_odds, -500.0, 500.0)
    gamma[soft] = 1.0 / (1.0 + np.exp(-log_odds))

    if return_strand_posterior:
        return gamma, _strand_posterior(llr_strand_full)
    return gamma
```

### 4.7 `run_em` changes

- Remove `strand_llr_mode`, `strand_z_threshold` parameters.
- Remove the `_aggregate_strand_z` block and the `strand_used` flag.
- Initialize `kappa_G = kappa_R = 20.0` before the first E-step.
  (Rationale: $\kappa = 20$ gives asymptotic $\mathrm{SD}(\hat p) \approx 0.11$ at $p = 0.5$ — wide enough to be uninformative about dispersion but not so wide as to collapse the LLR on iter 0. Any $\kappa_0 \in [10, 50]$ converges to the same fit within 3–5 iterations in tests.)
- After each E-step (from iter 0 γ-seed onwards), call `_estimate_kappa_G_soft` and `_estimate_kappa_R_soft` in the M-step.
- Record `kappa_G`, `kappa_R` in the history trace (replace `kappa`).
- Compute `region_gamma_strand` once at the end using `_strand_posterior(llr_final)`.
- Add `kappa_G`, `kappa_R`, `region_gamma_strand` fields to `EMFit`; drop `strand_used`, `strand_z`, `strand_llr_mode`, `kappa`.

Convergence criterion unchanged: `|π_soft − prev_π_soft| < 1e-4`. (The κ values converge on the same timescale as π_soft in practice; no separate convergence check needed.)

### 4.8 `CalibrationResult` and summary

- Add: `kappa_G: float`, `kappa_R: float`, `region_gamma_strand: np.ndarray | None`.
- Drop: `strand_used`, `strand_z`, `strand_llr_mode`, `kappa`.
- `to_summary_dict`: emit `kappa_G`, `kappa_R`, plus (optionally) summary stats of `region_gamma_strand` (e.g. fraction with γ̃ > 0.95, γ̃ < 0.05, for QC).

### 4.9 `CalibrationConfig`

Remove the `strand_llr_mode: str = "binomial"` field outright. Not
deprecation — removal. Any user YAML still containing this field will
fail with a clear `dataclass` error, which is the right signal.

---

## 5. Tests

### 5.1 New unit tests (`tests/test_strand_bb_rho.py`)

1. **`test_bb_logpmf_matches_scipy`**: assert `_bb_logpmf` agrees with `scipy.stats.betabinom.logpmf` across a grid of $(k, n, \alpha, \beta)$ to 1e-10.
2. **`test_bb_logpmf_binomial_limit`**: as $\kappa \to \infty$ (e.g. 1e6), BB collapses to Binomial. Check $|_bb\_logpmf(k, n, \kappa/2, \kappa/2) - \mathrm{Binom}\log\mathrm{pmf}(k; n, 0.5)| < 10^{-4}$.
3. **`test_bb_logpmf_edge_k_0_and_n`**: $k = 0$ and $k = n$ return finite, correct values.
4. **`test_rho_integration_normalizes`**: $\int_0^1 1 \,d\rho = 1$, i.e. the 17-node scheme applied to $f(\rho) = 1$ gives `sum(exp(_RHO_LOGW)) ≈ 1` to 1e-14.
5. **`test_rho_integration_linear`**: apply to $f(\rho) = \rho$; expect $1/2$ to 1e-14.
6. **`test_ss_half_gives_zero_llr`**: $\mathrm{SS} = 0.5$ → G and R are identical → LLR = 0 everywhere.
7. **`test_49_51_small_n`**: 49/51 at $n = 100$, $\mathrm{SS} = 0.81$, $\kappa_G = \kappa_R = 20$ → LLR $\in [0.5, 1.5]$ (weak G).
8. **`test_4912_5123_large_n`**: as above but $n = 10035$ → LLR $\in [2.5, 5.5]$ (moderate G).
9. **`test_mid_p_gdna_fluctuation`**: 60/140 at $n = 200$, $\mathrm{SS} = 0.81$, $\kappa_G = 20$ → LLR positive (current Binomial gives ≈ −30; regression of the actual win we're seeking).
10. **`test_strand_posterior_sanity`**: strong-G region → γ̃ > 0.9; strong-R region → γ̃ < 0.1.

### 5.2 κ estimator tests

11. **`test_kappa_G_recovers_from_synthetic`**: generate $N = 500$ synthetic regions drawn from $\mathrm{BB}(n_i, \kappa_\mathrm{true}/2, \kappa_\mathrm{true}/2)$ with $\kappa_\mathrm{true} \in \{5, 20, 100\}$, γ = 1 everywhere; assert `_estimate_kappa_G_soft` recovers $\kappa_\mathrm{true}$ within 15% relative error.
12. **`test_kappa_R_recovers_from_synthetic`**: same for R-class.
13. **`test_kappa_G_tiny_data`**: $N = 2$ regions → estimator returns a finite value in the bracket (no crash).
14. **`test_kappa_soft_weights_respected`**: synthesize 50/50 mixture of $\kappa_G = 5$ and $\kappa_R = 100$ regions; set $\gamma_i = \mathbb{1}[z_i = G]$; assert the two estimators recover the correct per-class κ.

### 5.3 End-to-end EM tests

15. **`test_run_em_no_z_gate`**: assert that EM is always-on for strand (no `strand_enabled` conditional remains in the public API).
16. **`test_run_em_neutral_init_escapes`**: initialize at all-γ=0.5; assert EM separates into G/R as expected (e.g. for a synthetic dataset with 50% G-class and 50% R-class regions, final γ bimodally distributes).
17. **`test_run_em_tiny_test_case`**: 10 regions, no crashes; κ fits bounded.
18. **`test_calibration_result_has_new_fields`**: `region_gamma_strand`, `kappa_G`, `kappa_R` present; old fields (`strand_used`, `kappa`, etc.) absent.

### 5.4 Regression / deletion of obsolete tests

- Delete `tests/test_betabinom_strand.py` — its premise (comparing binomial vs betabinom mode) is obsolete.
- Audit `tests/test_calibration*.py` for references to `strand_used`, `strand_z`, `strand_llr_mode`, `kappa`. Replace with new field accesses or remove.
- Tests that assert specific LLR numeric values will change. Update goldens once after the new model is validated; do not "fix" them pre-emptively.
- `tests/test_golden_output.py`: regenerate goldens via `pytest tests/ --update-golden` **after** manual validation on VCaP titration (§6) confirms the new model is correct.

---

## 6. Validation

After the implementation compiles and unit tests pass, validate on `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna{00,01,02,05,10,20,40,80}m/`.

### 6.1 Per-library calibration comparison

For each of the 8 libraries, capture:
- Fitted $\kappa_G$, $\kappa_R$.
- Mean / median / tails of $\mathrm{LLR}^{\text{strand}}_i$.
- Distribution of $\tilde\gamma^{\text{strand}}_i$ (strand-only posterior).
- Distribution of final $\gamma_i$ (fused).
- `quant.gdna_total` and ratio vs strand-inferred gDNA total.

Compare against current code (same libraries, current shipping model) on the same metrics.

### 6.2 Expected outcomes

Per the theory doc §5:

- $\kappa_G$: expect 10–50 for hybrid capture, possibly higher (near-Binomial) for WGS-style input. Should not pin at $\kappa_\min$ or $\kappa_\max$.
- $\kappa_R$: expect higher than $\kappa_G$ (RNA regions cluster tighter around SS than gDNA around 0.5), likely 30–200.
- $\tilde\gamma^{\text{strand}}$ on intronic hotspots in expressed genes: expect > 0.8 (strand says "this is gDNA"). Under current Binomial model many of these flip to $\tilde\gamma < 0.5$ due to $\hat p$ noise.
- Fused $\gamma$ should move *toward* strand on high-confidence strand regions. If density still dominates incorrectly on hybrid-capture hotspots, that is expected and diagnostic — the NB density plan will fix it.

### 6.3 Validation script

Create [`scripts/calibration/strand_rho_validation.py`](../../scripts/calibration/strand_rho_validation.py) (new). This is the functional replacement for the deleted `betabinom_pilot.py`. It:

1. Runs calibration on each VCaP library with the new model.
2. Dumps per-library: κ_G, κ_R, fraction of regions with γ̃ > 0.95, fraction with γ̃ < 0.05, quant.gdna_total, strand-inferred gDNA total (from γ̃ alone), their ratio.
3. Plots (if matplotlib available): $\tilde\gamma^{\text{strand}}$ vs current-model $\gamma$ per region, colored by $n_i$; $\hat p$ vs $n_i$ colored by γ̃.

### 6.4 Synthetic benchmark regression

Run `scripts/benchmark/configs/locus_simple_baseline.yaml` and confirm synthetic mRNA/nRNA/gDNA relative errors are within 1% of current. If not, diagnose before proceeding.

---

## 7. Rollout sequence

1. **Implement** §4.1–4.6 kernels + §4.7 EM driver changes. No config plumbing yet — hardwire hooks in `_em.py` first.
2. **Run new unit tests** (§5.1–5.3). All pass.
3. **Wire through** `CalibrationResult` (§4.8) and delete `CalibrationConfig.strand_llr_mode` (§4.9). Update `_calibrate.py`. Fix test compilation.
4. **Audit existing tests** (§5.4). Fix test failures that are legitimate numerical changes (do not change tolerances — update goldens). Do not update `tests/golden/` yet.
5. **Run validation** (§6) on VCaP titration. Document results inline in this doc's §6 (results section to be added).
6. **If validation passes**, regenerate goldens (`pytest tests/ --update-golden`).
7. **Update `docs/METHODS.md`** §5 to describe the new strand model.
8. **Commit** as a single atomic change (the whole strand rewrite is interdependent; no halfway state is meaningful).
9. **Return to NB density plan**.

---

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| `_estimate_kappa_R_soft` becomes an EM-time bottleneck on 700k-region human index. | Profile first. If >1s per M-step iter, cache per-ρ-node intermediates that don't depend on κ (the `lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)` term), or switch to a Laplace approximation around the current κ. |
| Numerical instability at boundary $\kappa \to 0^+$ (Beta(0, 0) limit). | `kappa_lo = 0.1` bracket avoids it. Add a single assertion `κ > 0` in `_bb_logpmf`. |
| γ-weighted MLE collapses when all γ ≈ 0 or all γ ≈ 1. | Handled in both κ estimators: check `w.sum() > 0`, return bracket midpoint otherwise. Also: such states only occur in pathological test cases; real data has γ distributed across [0, 1]. |
| $\rho$-integral under-resolved for very large $n_i$ (integrand becomes peaked). | 17-node GL is good enough for $n_i \le 10^4$. For bigger $n_i$ the integrand is peaked near $\hat p$; 17 nodes cover the whole domain so the peak is captured. If we see issues in practice, upgrade to an adaptive scheme. Defer. |
| Existing goldens change dramatically → regression worry. | Validation (§6) first, goldens regenerated only after manual confirmation on VCaP. Golden regeneration is a deliberate sign-off step, not an automatic byproduct. |
| Removing `strand_llr_mode` from `CalibrationConfig` breaks user YAML files. | Acceptable. This is a major internal change; breakage is expected and the error message is clear. Document in CHANGELOG. |

---

## 9. Estimated LOC footprint

- `_em.py`: **~+200 / −250** net (new kernels add ~250 lines; deleted Binomial + shared-κ BB + z-gate code removes ~300 lines; simplified `run_em` is net shorter).
- `_calibrate.py`: **~−10** (fewer arguments, simpler result construction).
- `_result.py`: **~+5 / −10**.
- `config.py`: **~−15**.
- Tests: **~+400** new, **~−350** deleted/refactored.
- Scripts: `betabinom_pilot.py` deleted (~250 lines); `strand_rho_validation.py` new (~200 lines).

Net: simpler codebase, more tests.

---

## 10. Sign-off checklist before coding

- [ ] User confirms: delete `strand_llr_mode` outright (no deprecation cycle).
- [ ] User confirms: delete `tests/test_betabinom_strand.py` outright.
- [ ] User confirms: $\kappa_0 = 20$ initial value (alternative: compute once from pooled γ=0.5 data before EM starts; more rigorous but adds code).
- [ ] User confirms: $\kappa$ bracket $[0.1, 1000]$.
- [ ] User confirms: 17-node Gauss–Legendre for ρ-integration (alternative: 33 nodes, trivially more expensive).
- [ ] User confirms: relax valid-region criterion from $n_i \ge 2$ to $n_i \ge 1$.

Once these are green, begin implementation at §7 step 1.
