# New Calibration Plan v6

Date: 2026-05-26
Status: implementation-ready plan after outside review
Supersedes: `new_cal_plan_v5.md`, `new_cal_plan_v4.md`, `new_cal_plan_v3.md`, `new_cal_plan_v2.md`, `new_cal_plan_v1.md`
Primary target: strand-specific hybrid capture. Secondary target: non-capture and unstranded data through the same region/boundary architecture.

v6 keeps the v5 direction: four expression x capture state probabilities, region/boundary terminology, boundary-to-contained gDNA imputation, and a vectorized state tensor. It incorporates the outside review only where the argument improves correctness or reduces implementation risk. It avoids speculative machinery that benchmarks have not yet justified.

---

## 0. Outside Review Disposition

| Concern | Decision | v6 action |
|---|---|---|
| C1 cumulative-product sweep is numerically unsafe | Accept | Drop cumulative-product vectorization. Use a simple per-reference sequential scan over excess Gamma evidence. Optimize later only if profiling demands it. |
| C2 `gamma_global` is not directly identifiable | Accept with simplification | Rename to `capture_enrichment_target`. Default to 1.0 and only move it with strong posterior captured evidence. No bait BED in v6. |
| C3 four-state model needs anchors | Accept | Add annotation-derived state priors and carry `background_seed_mask` into early passes as a soft prior boost, not a hard clamp. |
| C4 M-step couplings need damping | Accept lightly | Add fixed internal damping and report progress. Do not claim formal variational monotonicity. |
| C5 transfer weight conflates purity and information | Accept | Define boundary transfer weight as evidence reliability, not gDNA fraction. Keep it inside `boundary_sweep.py`; no separate module yet. |
| C6 tensor should use log Bayes factors, not posteriors | Accept | Rename likelihood terms to `logBF_*` and keep `state_log_prior` as the only prior term. |
| C7 `BoundaryTable` needs exactness contract | Accept | Specify exact region-boundary mapping and add `validate_boundary_table()` for tests. |
| C8 `A_r` units need pinning | Accept | Define `A_r` as unitless density multiplier with invariants. |
| C9 smooth gDNA vs boolean `enable_gdna` tension | Accept in principle | v6 path should keep gDNA smooth and enabled by default. Keep compatibility shims while migrating. |
| C10 `kappa_d` should be one library scalar | Accept | Compartment posteriors share one `kappa_d`. No per-compartment kappa in v6. |
| C11 Phase IV migration risk | Accept | Add v6 behind a feature flag initially and derive compatibility fields until benchmarks pass. |
| C12 minor items | Mixed | Accept opt-in region debug and strand helper reuse. Defer mappability hook and extra aliases unless touched by implementation. |

The guiding rule: fix known mathematical hazards now; do not build broad frameworks, alternative shrinkers, optional BED paths, or native boundary storage until benchmark evidence forces the issue.

---

## 1. Core Contract

For each region `r`, estimate posterior probabilities over four biological states:

```text
state 0: background             = P(not expressed, not captured | data)
state 1: gdna_only_capture      = P(not expressed, captured     | data)
state 2: expressed_capture      = P(expressed,     captured     | data)
state 3: expressed_offtarget    = P(expressed,     not captured | data)
```

The state order is fixed:

```python
STATE_BACKGROUND = 0
STATE_GDNA_ONLY_CAPTURE = 1
STATE_EXPRESSED_CAPTURE = 2
STATE_EXPRESSED_OFFTARGET = 3
N_STATES = 4

STATE_IS_EXPRESSED = np.array([False, False, True, True])
STATE_IS_CAPTURED = np.array([False, True, True, False])
STATE_NAMES = (
    "background",
    "gdna_only_capture",
    "expressed_capture",
    "expressed_offtarget",
)
```

Derived marginals:

```text
p_expressed(r) = p_states[r, expressed_capture] + p_states[r, expressed_offtarget]
p_captured(r) = p_states[r, gdna_only_capture] + p_states[r, expressed_capture]
```

Downstream EM receives:

```text
mu_gdna(r)       expected regional gDNA count
upper_gdna(r)    conservative upper credible bound for regional gDNA count
A_r              regional gDNA exposure multiplier for the denominator
rna_lower(r)     conservative RNA lower bound for protecting RNA from gDNA siphon
```

`A_r` is unitless:

```text
A_r = E[lambda_gdna(r) | data] / E[rho_off | data]
```

`A_r == 1` means the region's gDNA density matches the off-target background density. `A_r > 1` means the region has enriched gDNA opportunity, usually because of capture. `A_r` is not a raw count ratio and not the old `rho_post / rho_ref` pathway.

Required invariants:

```text
A_r is finite
A_r >= 0
A_r dtype is float32 on RegionCalibration
mu_gdna >= 0
upper_gdna >= mu_gdna when both are finite
p_states rows sum to 1 within tolerance
```

---

## 2. Regions And Boundaries

### 2.1 Regions

A region is a non-overlapping contiguous genomic interval. It is one sorted row in `RegionArrays` / `DensityObservation`:

```text
region_id = row index
fields    = ref_id, start, end, signature, ts_class, length
```

Regions keep contained evidence:

```text
contained_unspliced_pos
contained_unspliced_neg
contained_spliced_pos
contained_spliced_neg
contained_leff_gdna
```

### 2.2 Boundaries

For a reference with `N` local regions, there are `N + 1` local boundary slots:

```text
boundary 0       terminal boundary before region 0
boundary i       internal boundary between region i - 1 and region i
boundary N       terminal boundary after region N - 1
```

Region-to-boundary indexing:

```text
left_boundary(region i)  = boundary i
right_boundary(region i) = boundary i + 1
```

Global indexing:

```text
local_i = region_id - ref_region_offsets[ref_id]
left_b  = ref_boundary_offsets[ref_id] + local_i
right_b = ref_boundary_offsets[ref_id] + local_i + 1
```

The boundary table does not store `left_region` or `right_region` arrays. Those relationships are derivable from offsets. Helper methods may compute them for tests or debugging, but they are not primary storage.

### 2.3 `BoundaryTable`

```python
@dataclass(frozen=True, slots=True)
class BoundaryTable:
    boundary_pos: np.ndarray                  # int64[B]
    ref_id: np.ndarray                        # int32[B]
    ref_region_offsets: np.ndarray            # int64[n_ref + 1]
    ref_boundary_offsets: np.ndarray          # int64[n_ref + 1]
    is_terminal: np.ndarray                   # bool[B]

    # Mass landing in the region to the left of this boundary.
    left_region_unspliced_pos: np.ndarray     # float32[B]
    left_region_unspliced_neg: np.ndarray
    left_region_spliced_pos: np.ndarray
    left_region_spliced_neg: np.ndarray

    # Mass landing in the region to the right of this boundary.
    right_region_unspliced_pos: np.ndarray    # float32[B]
    right_region_unspliced_neg: np.ndarray
    right_region_spliced_pos: np.ndarray
    right_region_spliced_neg: np.ndarray

    left_region_boundary_leff: np.ndarray     # float64[B]
    right_region_boundary_leff: np.ndarray    # float64[B]
```

Exact current-ledger mapping for an internal boundary `b` between local regions `i - 1` and `i`:

```text
left_region_*[b]  = region i - 1 boundary_right_* channels
right_region_*[b] = region i     boundary_left_* channels
```

Terminal boundaries are zero-filled and flagged. Consumers must either skip `is_terminal` or expect zero crossing mass.

Add test-only validation:

```python
def validate_boundary_table(
    boundaries: BoundaryTable,
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
) -> None:
    """Recompute boundary table relationships and assert equality with ledger slots."""
```

First implementation is pure Python in `src/rigel/calibration/boundaries.py`. No scanner changes in v6.

---

## 3. Evidence Products

### 3.1 Compartment-aware strand deconvolution

The existing aggregate strand deconvolution remains useful, especially for the current path and for estimating `kappa_d`. v6 adds compartment-aware posteriors for calibration:

```python
@dataclass(frozen=True, slots=True)
class RegionGdnaChannelEstimate:
    contained_mean: np.ndarray
    contained_upper: np.ndarray
    contained_rna_lower: np.ndarray
    boundary_left_mean: np.ndarray
    boundary_left_upper: np.ndarray
    boundary_left_rna_lower: np.ndarray
    boundary_right_mean: np.ndarray
    boundary_right_upper: np.ndarray
    boundary_right_rna_lower: np.ndarray
    kappa_d: float
    p_r1_sense: float
    confidence: float
    flags: np.ndarray
```

Decision: `kappa_d` is one scalar per library. It is fit jointly from gDNA-like regions and shared by contained, boundary-left, and boundary-right posteriors. v6 does not fit per-compartment `kappa_d`.

In strand-specific data:

- contained gDNA estimates feed regional density evidence;
- boundary gDNA estimates feed boundary imputation;
- contained and boundary RNA lower bounds exclude regions from background training;
- expressed captured and expressed off-target regions can still inform gDNA density because strand deconvolution separates gDNA from RNA.

In unstranded data:

- strand fields are unavailable;
- unspliced boundary counts are gDNA-compatible evidence, not definitive gDNA;
- uncertainty should stay broad when expression and capture are confounded.

### 3.2 Background off-target density

The background model estimates `rho_off`, the off-target gDNA density.

Initial background candidates:

```text
(intergenic OR intron-only)
AND no spliced mass
AND no strand-RNA evidence when stranded
AND not in the initial top-T contained-density exclusion set
```

The top-T exclusion is computed once from initial intergenic/intronic candidates and carried as an exclusion mask. It should not compound with later filters in a way that silently starves the background fit.

Fit target:

```text
D_bg(r) ~ Poisson(rho_off * L_contained(r))
rho_off ~ Gamma(alpha_floor, beta_floor)
```

Use strand-deconvolved contained gDNA counts when available. Use unspliced contained counts in unstranded mode.

No exact zero mode. `rho_off` can be very small but should remain finite and positive through a weak prior.

Outputs:

```python
@dataclass(frozen=True, slots=True)
class BackgroundModel:
    rho_off_alpha: float
    rho_off_beta: float
    rho_off_mean: float
    seed_mask: np.ndarray          # bool[R]
    top_t_exclusion_mask: np.ndarray
    n_seed_regions: int
    fit_status: str                # "ok" | "sparse" | "prior_only"
    flags: np.ndarray              # uint16[R]
```

`seed_mask` is important. The four-state solver uses it as a soft early-pass anchor for `background`, not as a permanent hard label.

### 3.3 Boundary-to-contained imputation

For each region side, boundary evidence estimates a local gDNA density and predicts contained gDNA:

```text
lambda_side | B_side ~ Gamma(alpha_prior + B_side, beta_prior + L_side)
D_contained_from_side ~ NegativeBinomial(
    r = alpha_prior + B_side,
    p = (beta_prior + L_side) / (beta_prior + L_side + L_contained)
)
```

Combine left and right sides in excess-Gamma space:

```text
alpha_excess = B_left + B_right
beta_excess  = L_left + L_right
```

Then predict contained gDNA with the same Negative Binomial marginal.

Outputs:

```python
@dataclass(frozen=True, slots=True)
class BoundaryLocalPosterior:
    alpha_excess: np.ndarray       # float64[R]
    beta_excess: np.ndarray        # float64[R]
    mu_local: np.ndarray           # float32[R]
    upper_local: np.ndarray        # float32[R]
    flags: np.ndarray              # uint16[R]
```

### 3.4 Capture enrichment shrinkage

The outside review is right that a global captured enrichment is not directly identifiable without labels. v6 therefore avoids treating it as truth.

Use a conservative shrink target:

```text
capture_enrichment_target >= 1
```

Rules:

- Initialize `capture_enrichment_target = 1`.
- Move it above 1 only when posterior captured evidence is strong and broad enough.
- Dampen updates.
- If evidence is weak or non-capture-like, it stays near 1 and the captured/offtarget distinction naturally collapses.

This is not a user-facing knob. It is an internal stabilizer to keep noisy local boundary ratios from directly entering `A_r`.

---

## 4. Boundary Sweeps

Boundary sweeps move boundary-derived gDNA evidence into adjacent ambiguous regions. They are especially important for internal exon cases like the T1/T2 giant-exon example.

### 4.1 Message representation

Use excess Gamma evidence only:

```text
alpha_excess >= 0
beta_excess  >= 0
```

The prior is not propagated. This prevents double-counting when transfer weights are high.

### 4.2 Boundary transfer weight

The transfer weight is evidence reliability, not a gDNA fraction.

```text
w_boundary in [0, 1]
```

Stranded data can use strand-deconvolved boundary gDNA evidence and boundary RNA lower bounds. Unstranded data should use boundary opportunity and compatible boundary mass, while treating spliced mass as a reason to reduce confidence, not as proof of no gDNA.

First implementation lives inside `boundary_sweep.py`:

```python
def compute_boundary_transfer_weight(
    boundaries: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate | None,
    background: BackgroundModel,
) -> np.ndarray:
    """Return float64[B] transfer weights. Terminal boundaries return 0."""
```

No separate `transfer_weight.py` module in v6. No hard nonzero floor. If a boundary has no useful information, zero transfer is acceptable.

### 4.3 Sequential scan

Drop the cumulative-product formula from v5. It is numerically fragile. Use a plain per-reference scan over sorted arrays:

```text
f_alpha[0] = local_alpha[0]
f_beta[0]  = local_beta[0]

f_alpha[i] = local_alpha[i] + w_boundary[i] * f_alpha[i - 1]
f_beta[i]  = local_beta[i]  + w_boundary[i] * f_beta[i - 1]
```

Here `w_boundary[i]` is the internal boundary between local region `i - 1` and local region `i`.

Reverse scan uses the same rule on reversed regions. The implementation may use a Python loop per reference slice at first. This is simpler and safer than clever vectorization. If profiling later identifies this as hot, optimize the scan then.

### 4.4 Avoiding local double-counting

The final combined evidence for region `i` should include:

```text
local evidence for i
forward evidence transmitted from left neighbors
reverse evidence transmitted from right neighbors
```

It should not count region `i`'s local evidence three times. Implement scans so `forward_from_left[i]` excludes `local[i]`, and `reverse_from_right[i]` excludes `local[i]`.

Final:

```text
alpha_total = local_alpha + forward_from_left_alpha + reverse_from_right_alpha
beta_total  = local_beta  + forward_from_left_beta  + reverse_from_right_beta
```

Then convert to contained gDNA predictive mean and upper tail using the same Negative Binomial predictive distribution.

---

## 5. Four-State Tensor As Log Bayes Factors

### 5.1 State priors

The four-state model needs weak anchors. Use `region_arrays.signature` and early background seeds to build `state_log_prior[R, 4]`.

Principles:

- Intergenic regions: prior favors `background` and `gdna_only_capture` over expressed states.
- Intron-only regions: prior favors `background`, but less strongly than intergenic if strand-RNA evidence exists.
- Exonic regions: weakly informative or flat prior unless spliced/strand evidence is present.
- `background.seed_mask`: add a soft `background` prior boost for the first 1-2 passes.

Do not hard-clamp states except for impossible rows such as non-informative stranded orientation if the strand term cannot be evaluated.

### 5.2 Log Bayes factor terms

Use log Bayes factors, not posterior probabilities, as additive tensor terms. This avoids double-counting priors.

```text
logBF_expression      shape [R, 2]    columns: not_expressed, expressed
logBF_capture         shape [R, 2]    columns: not_captured, captured
logBF_gdna_density    shape [R, 2]    columns: off_target_density, captured_density
logBF_strand          shape [R, 4]    state-specific when strand evidence exists
state_log_prior       shape [R, 4]
```

Assembly:

```python
L = state_log_prior.copy()
L += logBF_expression[:, STATE_IS_EXPRESSED]
L += logBF_capture[:, STATE_IS_CAPTURED]
L += logBF_gdna_density[:, STATE_IS_CAPTURED]
L += logBF_strand

log_norm = logsumexp(L, axis=1, keepdims=True)
p_states = np.exp(L - log_norm).astype(np.float32)
```

Terms:

```text
expression: spliced mass, strand RNA lower bound, contained excess over gDNA upper
capture: boundary-imputed density relative to off-target density
gDNA density: NB predictive under rho_off vs enriched density
strand: existing exact/approx strand posterior machinery, shared kappa_d
```

Reuse existing exact strand helpers for small counts. Do not resurrect duplicate grids inside the tensor path.

### 5.3 Output

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
    n_passes: int
    converged: bool
    flags: np.ndarray                 # uint16[R]
```

Convenience properties expose `p_background`, `p_captured`, and `p_expressed`.

---

## 6. Calibration Iteration

The iteration is a small calibration loop, not the transcript EM.

### 6.1 E-step

Given current model parameters:

```text
fit local boundary posterior
run boundary sweep
build log Bayes factor tensor
normalize to p_states
compute mu_gdna, upper_gdna, rna_lower, A_r
```

### 6.2 M-step

Use posterior state weights:

```text
background weight    = p_background
not-expressed weight = p_background + p_gdna_only_capture
captured weight      = p_gdna_only_capture + p_expressed_capture
off-target weight    = p_background + p_expressed_offtarget
```

Refit:

```text
rho_off                    from background-weighted gDNA evidence
kappa_d                    from high-confidence not-expressed evidence, stranded only
capture_enrichment_target  from captured-weighted gDNA density, damped toward previous value
```

Boundary local posteriors and boundary sweeps are recalculated in the next E-step.

### 6.3 Damping and stop rule

Use internal damping for scalar updates:

```text
theta_next = (1 - eta) * theta_current + eta * theta_hat
```

Start with `eta = 0.5`. Keep this internal until benchmarks suggest otherwise.

Stop when both are true:

```text
max_abs(p_states_next - p_states_current) < p_tol
relative_change(rho_off) < rho_tol
```

Also report a diagnostic `sum_log_evidence = sum(logsumexp(L_r))`, but do not claim formal ELBO monotonicity.

Default `max_calibration_passes = 5`.

---

## 7. Downstream EM Contract

### 7.1 Numerator

Use `RegionCalibration.mu_gdna` and `RegionCalibration.upper_gdna` for regional gDNA prior mass.

Boundary evidence and smooth priors should not be converted into brittle per-locus on/off decisions. The v6 path should run with a gDNA component enabled by default unless the user explicitly disables gDNA at the config level.

During migration, existing compatibility fields may still carry an `enable_gdna` column, but it should be derived conservatively and should not drive the v6 conceptual model.

### 7.2 Denominator

Keep the agreed first implementation:

```text
gdna_eff_len(locus) = unweighted_gdna_eff_len(locus) * bp_weighted_mean(A_r)
```

Equivalent interpretation:

```text
gdna_eff_len_corrected = raw_locus_gdna_eff_len * sum_r(bp_r * A_r) / sum_r(bp_r)
```

The FL-convolved exposure integral remains deferred until benchmarks demand it.

### 7.3 Migration safety

Do not break current quantification while v6 is being built.

- Add `RegionCalibration` as a new field on `CalibrationResult`.
- Derive old compatibility fields from `RegionCalibration` only after the v6 path is ready.
- Gate downstream use behind an internal config flag at first, for example `use_v6_calibration`.
- Do not update existing golden outputs until the benchmark suite demonstrates improvement.

---

## 8. Implementation Phases

```text
Phase I  -> Phase II -> Phase III -> Phase IV -> Phase V
regions     boundary    four-state   downstream   benchmark
background  imputation  solver       wiring       refit
```

### Phase I - Regions, boundaries, strand compartments, background

Scope:

- remove old `400` precision-cap path and old density exposure source;
- add `src/rigel/calibration/boundaries.py` with `BoundaryTable`;
- add exact boundary-table validation tests;
- extend [strand_deconv.py](src/rigel/calibration/strand_deconv.py) with compartment-aware posteriors using one shared `kappa_d`;
- add `background_model.py` with `BackgroundModel` and `seed_mask`.

Verification:

- `tests/test_boundaries.py`: exact reconstruction of current ledger slots.
- `tests/test_compartment_strand_deconv.py`: compartment sums agree with aggregate path where expected.
- `tests/test_background_model.py`: `rho_off` recovers synthetic background and excludes strand-RNA-positive introns.
- Existing tests remain green because this phase is additive.

### Phase II - Boundary imputation and sweeps

Scope:

- add `boundary_model.py` for Gamma-Poisson / NB local boundary-to-contained prediction;
- add `boundary_sweep.py` with transfer weights and simple sequential scans;
- keep all transfer-weight constants internal;
- test internal exon propagation with the T1/T2 fixture.

Verification:

- Sparse non-capture boundaries shrink to background.
- Strong captured boundaries produce larger local/sweep gDNA predictions.
- A region with no local boundary evidence can receive attenuated evidence from neighbors.
- The sweep equals local posterior when all transfer weights are zero.

### Phase III - Four-state tensor solver

Scope:

- add `latent_states.py`;
- implement annotation-derived `state_log_prior`;
- implement log Bayes factor terms;
- implement damped calibration loop;
- output `RegionCalibration`.

Verification:

- Four archetype synthetic fixture identifies the correct dominant state when evidence is informative.
- Ambiguous unstranded regions remain uncertain and flagged.
- State rows sum to one.
- Iteration converges or reports non-convergence within five passes.

### Phase IV - Downstream wiring and diagnostics

Scope:

- rewire [_orchestrator.py](src/rigel/calibration/_orchestrator.py) to produce `RegionCalibration`;
- update [exposure.py](src/rigel/calibration/exposure.py) to build exposure from `RegionCalibration.A_r`;
- update [prior.py](src/rigel/calibration/prior.py) to consume `mu_gdna`, `upper_gdna`, and `A_r`;
- update [pipeline.py](src/rigel/pipeline.py) behind `use_v6_calibration`;
- keep compatibility fields until benchmarks pass.

Verification:

- Existing golden outputs unchanged when `use_v6_calibration=False`.
- v6-on outputs are deterministic and get their own temporary tests.
- Region debug output is opt-in only.

### Phase V - Benchmarks and defaults

Scope:

- emit post-capture truth from simulation: transcript abundances and FL distributions after capture;
- benchmark unstranded/non-capture, unstranded/capture, stranded/non-capture, stranded/capture;
- compare v6-off and v6-on;
- tune internal defaults only when benchmark evidence demands it;
- decide whether to make v6 the default.

Verification:

- No regression in non-capture data.
- Stranded capture high-gDNA improves gDNA-to-RNA leakage.
- Unstranded capture reports appropriate uncertainty.
- Existing goldens update only after benchmark acceptance.

---

## 9. Diagnostics

Compact user-facing outputs:

```text
quant.feather:       mean_A_r, n_regions_captured, n_regions_expressed
gene_quant.feather:  mean_A_r, n_regions_captured, n_regions_expressed
loci.feather:        mean_A_r, state_mass summaries
summary.json:        pass-level calibration diagnostics
```

Optional region debug output, behind an explicit flag such as `--debug-regions`:

```text
region_id
p_background
p_gdna_only_capture
p_expressed_capture
p_expressed_offtarget
p_expressed
p_captured
mu_gdna
upper_gdna
rna_lower
A_r
gamma_r
boundary_local_mean
boundary_sweep_mean
strand_contained_gdna_mean
strand_boundary_left_gdna_mean
strand_boundary_right_gdna_mean
flags
```

Pass-level diagnostics:

```text
rho_off per pass
capture_enrichment_target per pass
max state probability shift
relative rho_off shift
sum_log_evidence diagnostic
n_background_seed_regions
n_top_t_excluded_regions
converged
n_passes
```

---

## 10. Minimal User-Facing Knobs

Expose only:

```text
gdna_density_confidence     existing confidence concept; keep old name for now
background_trim_fraction    one robustness control
max_calibration_passes      safety cap
use_v6_calibration          temporary internal/advanced migration flag
```

Do not expose transfer-weight constants, damping, capture shrinkage hyperparameters, or tensor penalties until benchmarks prove they need tuning.

---

## 11. Residual Issues

### R1. Boundary transfer remains approximate

The transfer weight is an evidence-reliability approximation. It should be tested on synthetic boundary patterns before adding sophistication.

### R2. Capture enrichment is weakly identified without probes

v6 handles this by shrinking toward 1 and only moving with strong evidence. This is not a full capture-target discovery model yet.

### R3. Unstranded captured/expression mixtures remain ambiguous

v6 predicts gDNA with uncertainty rather than disabling gDNA. It cannot fully separate all-captured/all-expressed unstranded regions without additional evidence.

### R4. Exact strand helpers must stay contained

If exact small-count enumeration is needed, keep it inside existing strand helpers or a small wrapper. Do not let region-wise loops leak into the tensor assembly.

### R5. BoundaryTable is transitional storage

The first implementation derives boundaries from current region-side slots. Native scanner boundary storage is deferred until profiling or code clarity justifies it.

### R6. Mappability remains deferred

Do not add a mappability hook in v6 unless the touched code naturally exposes one. The model should still be written so gDNA opportunity can later be replaced by mappability-corrected opportunity.

---

## 12. Short Version

v6 is v5 with the risky parts tightened:

```text
Regions + boundaries
-> one-library kappa_d and compartment strand posteriors
-> background rho_off with seed_mask
-> boundary-to-contained Gamma/NB imputation
-> simple sequential boundary sweeps over excess evidence
-> [R,4] log Bayes factor tensor with annotation priors
-> damped calibration loop
-> RegionCalibration
-> compatibility-gated downstream wiring
```

This stays lightweight: no new C++, no probe BED, no native boundary storage, no separate transfer-weight framework, no hard zero mode, and no broad user knob surface. The plan is now implementation-ready enough to start Phase I while leaving the larger modeling choices to benchmarks.
