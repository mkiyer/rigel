# Phase III - Four-State Tensor Solver Implementation Plan

Date: 2026-05-26
Status: initial additive implementation complete
Scope: Implement the four-state log Bayes factor tensor, annotation-derived priors, calibration E/M iteration, and `RegionCalibration` output.

---

## Overview

Phase III builds the core four-state solver that combines evidence from:

1. **Expression evidence**: spliced mass, strand RNA lower bounds, contained gDNA vs. expected
2. **Capture evidence**: boundary mixture posteriors, captured density enrichment
3. **gDNA density evidence**: Negative Binomial predictive under off-target vs. enriched density
4. **Strand evidence**: existing exact/approx strand posterior machinery, shared `kappa_d`

The solver runs an E/M iteration loop that refits `rho_off`, `kappa_d`, and optionally `capture_enrichment_target` while maintaining state posteriors.

Output is `RegionCalibration`, which contains per-region state probabilities, gDNA estimates, `A_r` exposure multipliers, `gamma_r` captured-density multipliers, and diagnostic flags.

---

## 1. State Definitions And Priors

### 1.1 Four biological states

```python
STATE_BACKGROUND = 0              # not expressed, not captured
STATE_GDNA_ONLY_CAPTURE = 1       # not expressed, captured
STATE_EXPRESSED_CAPTURE = 2       # expressed, captured
STATE_EXPRESSED_OFFTARGET = 3     # expressed, not captured
N_STATES = 4

STATE_IS_EXPRESSED = np.array([False, False, True, True], dtype=bool)
STATE_IS_CAPTURED = np.array([False, True, True, False], dtype=bool)
STATE_NAMES = ("background", "gdna_only_capture", "expressed_capture", "expressed_offtarget")
```

Derived marginals:

```text
p_expressed(r) = p_states[r, 2] + p_states[r, 3]
p_captured(r) = p_states[r, 1] + p_states[r, 2]
```

### 1.2 Annotation-derived state log priors

Initialize `state_log_prior[R, 4]` using region signature and background seed mask:

**Intergenic regions** (no annotation):

```text
log_prior[background] = 1.0       (soft preference)
log_prior[gdna_only_capture] = 0.0  (possible, but not annotation-favored)
log_prior[expressed_capture] = -2.0  (penalize expression)
log_prior[expressed_offtarget] = -2.0  (penalize expression)
```

**Intronic regions** (annotation exists but not exonic):

```text
log_prior[background] = 0.75      (soft preference)
log_prior[gdna_only_capture] = 0.0  (possible, but not annotation-favored)
log_prior[expressed_capture] = -1.0  (weakly penalize expression)
log_prior[expressed_offtarget] = -1.0  (weakly penalize expression)
```

Rationale: the annotation prior should encode only what annotation tells us before evidence. Intergenic and intron-only regions are most likely background by default. Capture enrichment should be inferred from boundary/density likelihoods, not granted by annotation alone. The captured state remains close enough to background that strong boundary evidence can overcome the prior, while expressed states start lower unless spliced or strand-RNA evidence supports them.

**Exonic regions**:

```text
log_prior[background] = 0.0       (neutral)
log_prior[gdna_only_capture] = 0.0  (neutral)
log_prior[expressed_capture] = 0.0  (neutral)
log_prior[expressed_offtarget] = 0.0  (neutral)
```

**Background seed regions** (from Phase I):

For the first 1-2 passes, add soft anchor boost to `background` state:

```text
if background.seed_mask[r]:
    log_prior[r, background] += 1.0  (soft prior boost; can be overcome by data)
```

**No hard clamping**: regions are not impossible. Even intergenic regions can be expressed (novel transcription); exonic regions can be background (covered by gDNA). The prior is soft guidance only.

### 1.3 Construction method

```python
def build_state_log_prior(
    region_arrays: RegionArrays,
    background: BackgroundModel,
    *,
    pass_index: int = 0,
    background_boost: float = 1.0,
) -> np.ndarray:
    """Return [R, 4] state log priors.
    
    Args:
        region_arrays: RegionArrays with signature and ts_class fields
        background: BackgroundModel with seed_mask
        pass_index: calibration pass (0-indexed); apply boost on pass < 2
        background_boost: magnitude of soft background anchor
    
    Returns:
        state_log_prior: float64[R, 4]
    """
```

---

## 2. Log Bayes Factor Terms

The tensor assembly uses orthogonal likelihood terms, avoiding double-counting of priors.

### 2.1 Expression term: `logBF_expression[R, 2]`

Columns: `[not_expressed, expressed]`

Evidence:

- Spliced mass strongly penalizes not-expressed states.
- Strand RNA lower bound (when available) penalizes not-expressed.
- Contained gDNA excess over expected (when available) supports expressed.

Implementation note: current `DensityObservation` labels `spliced_count` as diagnostic-only for the v5 density model. Phase III intentionally changes that contract for the v6 tensor path only: `spliced_count` becomes expression evidence, while the old density model remains unchanged until Phase IV wiring.

Implementation:

```python
region_count = int(observation.contained_count.shape[0])
logBF_expression = np.zeros((region_count, 2), dtype=np.float64)

# Spliced mass term
spliced_log_penalty = -10.0  # strong penalty; tune if needed
logBF_expression[np.asarray(observation.spliced_count) > 0, 0] += spliced_log_penalty

# Strand RNA lower bound term (when stranded)
if strand_channels is not None:
    rna_lb = np.asarray(strand_channels.contained_rna_lower)
    logBF_expression[rna_lb > 0.1, 0] += -5.0  # penalize not-expressed

# Contained gDNA vs. expected gDNA term
# If contained count exceeds upper_gdna, that supports expressed state
# This will be added in the density term instead to avoid double-counting
```

### 2.2 Capture term: `logBF_capture[R, 2]`

Columns: `[not_captured, captured]`

Evidence:

- Boundary mixture posterior predicts captured probability.
- Shrunk local enrichment gamma supports captured density.

Implementation:

```python
region_count = int(observation.contained_count.shape[0])
logBF_capture = np.zeros((region_count, 2), dtype=np.float64)

# Boundary imputation provides posterior on captured vs. offtarget
# p_captured_local = mu_boundary_captured / (mu_boundary_captured + mu_boundary_offtarget)
off_target_mu = background.rho_off_mean * observation.contained_leff
excess_mu = np.maximum(sweep.mu_sweep - off_target_mu, 0.0)
p_cap = np.clip(excess_mu / (excess_mu + off_target_mu + eps), eps, 1.0 - eps)
logBF_capture[:, 1] = np.log(p_cap)
logBF_capture[:, 0] = np.log(1.0 - p_cap)
```

Implementation note: the first implementation uses excess over the off-target background, not raw `sweep.mu_sweep`, for capture probability. This keeps no-excess regions from starting at `p_captured ~= 0.5` simply because `sweep.mu_sweep` includes the background Gamma prior.

### 2.3 gDNA density term: `logBF_gdna_density[R, 2]`

Columns: `[off_target_density, captured_density]`

Evidence:

- Negative Binomial predictive of contained gDNA under `rho_off` vs. enriched density.
- Strand-deconvolved contained gDNA estimate (when available).

Implementation uses the same NB predictive math from Phase II boundary model:

```python
def logBF_gdna_density(
    observation: DensityObservation,
    background: BackgroundModel,  # rho_off_alpha, rho_off_beta
    sweep: BoundarySweepResult,   # mu_sweep, upper_sweep, alpha_excess, beta_excess
    strand_channels: RegionGdnaChannelEstimate | None = None,
) -> np.ndarray:  # [R, 2]
    """Log Bayes factor for gDNA density: off-target vs. enriched.
    
    Uses boundary sweep excess evidence to predict enrichment level.
    Negative Binomial predictive likelihood over observed contained counts.
    """
    contained_count = np.asarray(observation.contained_count, dtype=np.float64)
    contained_leff = np.asarray(observation.contained_leff, dtype=np.float64)
    region_count = int(contained_count.shape[0])
    logBF = np.zeros((region_count, 2), dtype=np.float64)
    
    if strand_channels is not None:
        # Use strand-deconvolved gDNA counts
        observed_gdna = np.asarray(strand_channels.contained_mean, dtype=np.float64)
    else:
        # Use unspliced contained counts as gDNA-compatible proxy
        observed_gdna = np.asarray(contained_count, dtype=np.float64)
    
    leff = np.asarray(contained_leff, dtype=np.float64)
    leff = np.maximum(leff, 1.0e-12)
    
    # Off-target density posterior (background Gamma only)
    alpha0 = float(background.rho_off_alpha)
    beta0 = float(background.rho_off_beta)
    p_nb_offtarget = beta0 / (beta0 + leff)
    p_nb_offtarget = np.clip(p_nb_offtarget, 1.0e-10, 1.0 - 1.0e-10)
    
    # Enriched density posterior (background + boundary excess)
    alpha_enriched = alpha0 + np.maximum(sweep.alpha_excess, 0.0)
    beta_enriched = beta0 + np.maximum(sweep.beta_excess, 0.0)
    p_nb_enriched = beta_enriched / (beta_enriched + leff)
    p_nb_enriched = np.clip(p_nb_enriched, 1.0e-10, 1.0 - 1.0e-10)
    
    # NB logpmf for observed count. Region counts can be fractional from
    # fractional routing, so use rounded counts for the NB likelihood term.
    from scipy.stats import nbinom
    observed_gdna_int = np.rint(np.maximum(observed_gdna, 0.0)).astype(np.int64)
    logBF[:, 0] = nbinom.logpmf(observed_gdna_int, alpha0, p_nb_offtarget)
    logBF[:, 1] = nbinom.logpmf(observed_gdna_int, alpha_enriched, p_nb_enriched)
    
    return np.nan_to_num(logBF, nan=-20.0, posinf=0.0, neginf=-20.0)
```

### 2.4 Strand term: `logBF_strand[R, 4]`

Columns: `[background, gdna_only, expressed_capture, expressed_offtarget]`

When strand channels are available, the first implementation uses the compartment strand summary (`contained_rna_lower` and `contained_precision`) as a conservative state-specific term. Exact small-count strand tensor wiring is deferred until raw folded strand counts are passed into this layer.

When not available, return zeros (no evidence to distinguish states).

```python
def logBF_strand(
    strand_channels: RegionGdnaChannelEstimate | None,
    region_count: int,
) -> np.ndarray:  # [R, 4]
    """Log Bayes factor for strand evidence.
    
    Uses exact strand posteriors for gDNA vs. RNA separation.
    Reuses existing helpers from strand_deconv.py.
    """
    if strand_channels is None:
        return np.zeros((int(region_count), 4), dtype=np.float64)
    
    actual_region_count = int(strand_channels.flags.shape[0])
    if actual_region_count != int(region_count):
        raise ValueError(
            "logBF_strand: strand_channels region count does not match observation; "
            f"got {actual_region_count} and {region_count}."
        )
    logBF = np.zeros((actual_region_count, 4), dtype=np.float64)
    
    # For not-expressed states, add a conservative penalty when RNA lower-bound
    # evidence exists. Exact `strand_log_likelihood_d_grid()` wiring is deferred
    # because RegionGdnaChannelEstimate does not carry raw folded counts.
    
    return logBF
```

---

## 3. Tensor Assembly And Normalization

```python
def build_state_tensor(
    state_log_prior: np.ndarray,           # [R, 4]
    logBF_expression: np.ndarray,          # [R, 2]
    logBF_capture: np.ndarray,             # [R, 2]
    logBF_gdna_density: np.ndarray,        # [R, 2]
    logBF_strand: np.ndarray,              # [R, 4]
) -> np.ndarray:  # [R, 4]
    """Assemble log Bayes factor tensor and normalize to probabilities.
    
    STATE_IS_EXPRESSED = [False, False, True, True]
    STATE_IS_CAPTURED = [False, True, True, False]
    """
    from scipy.special import logsumexp
    
    log_tensor = state_log_prior.copy()
    log_tensor[:, STATE_IS_EXPRESSED] += logBF_expression[:, 1:2]
    log_tensor[:, ~STATE_IS_EXPRESSED] += logBF_expression[:, 0:1]
    log_tensor[:, STATE_IS_CAPTURED] += logBF_capture[:, 1:2]
    log_tensor[:, ~STATE_IS_CAPTURED] += logBF_capture[:, 0:1]
    log_tensor[:, STATE_IS_CAPTURED] += logBF_gdna_density[:, 1:2]
    log_tensor[:, ~STATE_IS_CAPTURED] += logBF_gdna_density[:, 0:1]
    log_tensor += logBF_strand
    
    # Normalize
    log_norm = logsumexp(log_tensor, axis=1, keepdims=True)
    log_p_states = log_tensor - log_norm
    p_states = np.exp(log_p_states).astype(np.float32)
    p_states = np.clip(p_states, 0.0, 1.0)  # safety clip for numerical errors
    
    return p_states
```

---

## 4. E-Step: Evidence-Based State Update

Given current model parameters (`rho_off`, `kappa_d`, `capture_enrichment_target`):

```python
def calibration_e_step(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    boundaries: BoundaryTable,
    background: BackgroundModel,
    local_posterior: BoundaryLocalPosterior | None = None,
    strand_channels: RegionGdnaChannelEstimate | None = None,
    *,
    pass_index: int = 0,
) -> CalibrationStepResult:
    """E-step: compute state posteriors, gDNA/RNA bounds, and diagnostics.
    
    Returns:
        p_states: [R, 4] state probabilities
        mu_gdna: [R] expected contained gDNA
        upper_gdna: [R] upper credible bound on contained gDNA
        rna_lower: [R] lower bound on RNA (protected from gDNA siphon)
        A_r: [R] exposure multiplier
        gamma_r: [R] captured-density multiplier
    """
    
    # 1. Local boundary posterior and boundary sweep
    if local_posterior is None:
        local_posterior = build_boundary_local_posterior(
            observation,
            boundaries,
            background,
            strand_channels=strand_channels,
            confidence=0.95,
        )

    sweep = run_boundary_sweep(
        local_posterior, boundaries, observation, background,
        strand_channels=strand_channels,
        confidence=0.95
    )
    
    # 2. Build state priors
    state_prior = build_state_log_prior(
        region_arrays, background,
        pass_index=pass_index,
        background_boost=1.0
    )
    
    # 3. Build likelihood terms
    logBF_expr = build_logBF_expression(observation, strand_channels)
    logBF_capt = build_logBF_capture(sweep, observation, background)
    logBF_density = build_logBF_gdna_density(observation, background, sweep, strand_channels)
    logBF_str = build_logBF_strand(strand_channels, observation.contained_leff.shape[0])
    
    # 4. Assemble tensor and normalize
    log_tensor = build_state_log_tensor(
        state_prior, logBF_expr, logBF_capt, logBF_density, logBF_str
    )
    p_states = normalize_state_log_tensor(log_tensor)
    sum_log_evidence = float(np.sum(logsumexp(log_tensor, axis=1), dtype=np.float64))
    
    # 5. Derive gDNA/RNA bounds from sweep and state weights
    p_captured = p_states[:, STATE_IS_CAPTURED].sum(axis=1)
    
    mu_offtarget = background.rho_off_mean * observation.contained_leff
    upper_offtarget = predict_contained_gdna_from_excess(
        background,
        observation.contained_leff,
        np.zeros_like(observation.contained_leff),
        np.zeros_like(observation.contained_leff),
        confidence=0.95,
    )[1]
    mu_captured = np.asarray(sweep.mu_sweep, dtype=np.float64)
    upper_captured = np.asarray(sweep.upper_sweep, dtype=np.float64)

    mu_gdna = (p_captured * mu_captured + (1.0 - p_captured) * mu_offtarget).astype(
        np.float32
    )
    upper_gdna = (
        p_captured * upper_captured + (1.0 - p_captured) * upper_offtarget
    ).astype(np.float32)
    upper_gdna = np.maximum(upper_gdna, mu_gdna)
    
    # RNA lower bound: conservative estimate assuming observed gDNA is wrong
    if strand_channels is None:
        rna_lower = np.zeros_like(observation.contained_leff, dtype=np.float32)
    else:
        rna_lower = np.maximum(
            np.asarray(strand_channels.contained_rna_lower, dtype=np.float64),
            0.0,
        ).astype(np.float32)
    
    # Exposure multipliers
    denom = background.rho_off_mean * observation.contained_leff + 1.0e-12
    A_r = mu_gdna / denom
    A_r = np.clip(A_r, 0.0, np.inf).astype(np.float32)
    gamma_r = mu_captured / denom
    gamma_r = np.clip(gamma_r, 0.0, np.inf).astype(np.float32)
    flags = derive_region_flags(p_states, local_posterior, strand_channels)
    
    return CalibrationStepResult(
        p_states=p_states,
        mu_gdna=mu_gdna,
        upper_gdna=upper_gdna,
        rna_lower=rna_lower,
        A_r=A_r,
        gamma_r=gamma_r,
        flags=flags,
        local_posterior=local_posterior,
        sweep=sweep,
        log_tensor=log_tensor,
        sum_log_evidence=sum_log_evidence,
    )
```

---

## 5. M-Step: Model Parameter Refit

After E-step, use posterior state weights to refit model components:

### 5.1 Background weight and refit

```text
w_background[r] = p_states[r, STATE_BACKGROUND]
```

Refit `rho_off` using background-weighted gDNA evidence from Phase I background training regions. This is a posterior update: do not retrain candidate regions; only reweight existing candidates.

### 5.2 Strand deconvolution refit (when stranded)

```text
w_notexpr[r] = p_states[r, STATE_BACKGROUND] + p_states[r, STATE_GDNA_ONLY_CAPTURE]
```

Refit `kappa_d` using high-confidence not-expressed regions (w_notexpr > 0.8). Use existing `estimate_kappa_d()` helper, but gate it on the weight threshold.

### 5.3 Capture enrichment shrinkage (optional for v6)

```text
w_captured[r] = p_states[r, STATE_GDNA_ONLY_CAPTURE] + p_states[r, STATE_EXPRESSED_CAPTURE]
```

Estimate enriched density from captured-weighted boundary gDNA. Apply damping toward 1.0 to avoid overcommitting to noisy local ratios.

---

## 6. Convergence And Iteration Control

### 6.1 Stop rule

Stop when both conditions are met:

```text
max_abs(p_states_new - p_states_old) < p_tol  (e.g., 0.01)
relative_change(rho_off) < rho_tol  (e.g., 0.02)
```

Also stop after `max_calibration_passes` (default 5).

### 6.2 Diagnostics reported per pass

```text
pass_index
max_state_shift: max absolute change in any state probability
rho_off: current off-target density estimate
kappa_d: current strand specificity parameter
capture_enrichment_target: (optional) current enrichment shrinkage target
converged: boolean, whether both thresholds met
n_regions_expressed: count of regions with p_expressed > 0.5
n_regions_captured: count of regions with p_captured > 0.5
n_regions_background: count of regions with p_background > 0.5
```

---

## 7. RegionCalibration Output Schema

```python
@dataclass(frozen=True, slots=True)
class RegionCalibration:
    p_states: np.ndarray              # float32[R, 4]
    mu_gdna: np.ndarray               # float32[R]
    upper_gdna: np.ndarray            # float32[R]
    rna_lower: np.ndarray             # float32[R]
    A_r: np.ndarray                   # float32[R]
    gamma_r: np.ndarray               # float32[R]
    rho_off: float
    kappa_d: float | None
    capture_enrichment_target: float  # scalar internal stabilizer, reported diagnostically
    n_passes: int
    converged: bool
    flags: np.ndarray                 # uint16[R]
    pass_diagnostics: tuple[dict[str, object], ...]  # per-pass summary
    
    @property
    def p_background(self) -> np.ndarray:
        return self.p_states[:, 0]
    
    @property
    def p_gdna_only_capture(self) -> np.ndarray:
        return self.p_states[:, 1]
    
    @property
    def p_expressed_capture(self) -> np.ndarray:
        return self.p_states[:, 2]
    
    @property
    def p_expressed_offtarget(self) -> np.ndarray:
        return self.p_states[:, 3]
    
    @property
    def p_expressed(self) -> np.ndarray:
        return self.p_states[:, [2, 3]].sum(axis=1)
    
    @property
    def p_captured(self) -> np.ndarray:
        return self.p_states[:, [1, 2]].sum(axis=1)
```

### 7.1 Flags (uint16 per region)

```python
FLAG_STATE_AMBIGUOUS = 1 << 0      # max state prob < 0.6
FLAG_STRAND_UNINFORMATIVE = 1 << 1  # when stranded, kappa_d could not be fit
FLAG_BOUNDARY_SPARSE = 1 << 2      # minimal boundary evidence
FLAG_EXPRESSED_UNCERTAIN = 1 << 3  # p_expressed in (0.3, 0.7)
FLAG_EXACT_STRAND_POSTERIOR = 1 << 4  # used exact enumeration
```

---

## 8. Implementation File Structure

Primary files to create:

- `src/rigel/calibration/latent_states.py`: state constants, priors, likelihood assembly
- `src/rigel/calibration/calibration_iteration.py`: E/M loop, convergence, RegionCalibration builder
- Tests:
  - `tests/test_latent_states.py`: state prior construction, likelihood term correctness
  - `tests/test_calibration_iteration.py`: E/M convergence, synthetic fixtures

Modified files:

- `src/rigel/calibration/__init__.py`: export new Phase III API
- `src/rigel/calibration/_result.py`: add `RegionCalibration` to result schema (Phase IV)

---

## 9. Verification / Test Strategy

### 9.1 Unit tests (latent_states.py)

- State prior construction produces correct log priors for intergenic/intronic/exonic regions
- State prior boost applies and can be overcome by data
- Likelihood term sign/magnitude: spliced mass penalizes not-expressed, boundary-enriched captures supports captured
- Tensor assembly: rows sum to one, no NaN/Inf output

### 9.2 Integration tests (calibration_iteration.py)

- Four archetype synthetic fixture (pure background, pure gDNA capture, pure expression, mixed) converges to correct dominant state
- State posteriors sum to 1.0 ± 1e-5
- Convergence stop rule activates within 5 passes
- E-step produces valid gDNA bounds: mu_gdna <= upper_gdna, both >= 0
- E-step produces finite non-negative `A_r` and `gamma_r` values
- M-step refits rho_off and kappa_d without crashing

### 9.3 Synthetic fixture: four-state convergence

Create synthetic observation where:

- Region 1: only background (no fragments)
- Region 2: only gDNA capture (high unspliced, no strand evidence)
- Region 3: only expression off-target (high spliced, no boundary)
- Region 4: expression + capture (spliced + boundary, mixed strand evidence)

Run iteration, verify:
- Region 1 converges to background
- Region 2 converges to gdna_only_capture
- Region 3 converges to expressed_offtarget
- Region 4 shows high probability on expressed_capture

---

## 10. Design Decisions And Alternatives

### 10.1 Why `state_log_prior` vs. `state_posterior`?

**Decision**: Use log priors to avoid double-counting priors when building the tensor.

Alternative considered: Hard clamping of states (e.g., intergenic can't be expressed). Rejected because novel/unannotated expression is real and we want to detect it.

### 10.2 Why background boost only in early passes?

**Decision**: Background seed mask provides useful anchor in pass 0-1, but if it persists after pass 2, it can lock regions into background incorrectly.

Tune: `background_boost = 1.0` if `pass_index < 2` else `0.0`.

### 10.3 Why NB predictive for gDNA density, not exact Poisson?

**Decision**: The Negative Binomial predictive marginalizes over unobserved effective length uncertainty (via beta parameter), which Poisson does not. This gives wider uncertainty bands when boundary evidence is sparse, which is honest.

### 10.4 Why not fit `capture_enrichment_target` in v6?

**Decision**: Per-region enrichment is very noisy in unstranded or low-capture data. A global shrinkage target is useful, but identifying it requires strong consistent captured-state posterior across many regions. Defer to Phase V benchmarks.

For v6, keep it fixed at 1.0 or use a simple damped moving average.

### 10.5 Why separate `strand_channels` parameter instead of merging with `observation`?

**Decision**: `DensityObservation` records raw counts; `RegionGdnaChannelEstimate` records strand-deconvolved estimates. Keeping them separate makes the dependency graph clear and makes unstranded mode trivial (pass None).

---

## 11. Known Residual Issues

### R1. Exact strand posterior in tensor is not yet implemented

The initial Phase III implementation uses `RegionGdnaChannelEstimate` summaries to add a conservative strand term when RNA lower-bound evidence exists. Exact integration with `strand_log_likelihood_d_grid()` remains deferred because the tensor layer does not yet receive raw folded strand counts (`k_sense`, `n_total`) per compartment.

### R2. Capture enrichment shrinkage is deferred

v6 uses fixed `capture_enrichment_target = 1.0` or simple moving average. Gamma-Poisson EB shrinkage is a future refinement.

### R3. Convergence diagnostics do not include ELBO

We report max state shift and rho_off shift, but not a formal ELBO. The iteration is empirical convergence, not formal variational EM.

### R4. Region-level pass-by-pass debugging is not yet captured

Optional per-region debug output can be added in Phase IV if needed.

### R5. `spliced_count` changes role in the v6 tensor path

Current `DensityObservation` treats `spliced_count` as density-model diagnostics only. Phase III intentionally consumes it as expression evidence, but only inside `latent_states.py`; this should be called out in code comments/tests so the old density model contract is not accidentally changed during Phase III.

### R6. Scalar M-step is partial in the additive implementation

The initial implementation refits `rho_off` from background-weighted seed evidence and updates `capture_enrichment_target` diagnostically. It carries `kappa_d` from existing strand deconvolution output rather than refitting it because `run_calibration_iteration()` does not yet receive `PayloadArrays`/folded strand counts. This should be completed during Phase IV wiring if benchmarks show the scalar loop needs live kappa updates.

---

## 12. Next Steps

After the initial additive implementation:

1. Review exact strand tensor requirements and decide whether to pass raw folded strand counts into `latent_states.py`.
2. Wire `RegionCalibration` into `CalibrationResult` behind the Phase IV migration flag.
3. Update exposure and prior assembly to consume `A_r`, `mu_gdna`, and `upper_gdna` when v6 is enabled.
4. Run v6-on smoke tests before any golden updates.
5. Move to Phase V benchmarks before making v6 the default.
