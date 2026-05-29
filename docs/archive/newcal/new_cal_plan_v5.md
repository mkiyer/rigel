# New Calibration Plan v5

Date: 2026-05-26
Status: implementation-ready rewrite after v4 review
Supersedes: `new_cal_plan_v4.md`, `new_cal_plan_v3.md`, `new_cal_plan_v2.md`, `new_cal_plan_v1.md`
Primary target: strand-specific hybrid capture. Secondary target: non-capture and unstranded data through the same architecture.

v4 got the statistical model into the right shape: expression x capture is a continuous joint latent field, and boundary flux is local gDNA evidence rather than a threshold. v5 keeps that model but makes the implementation plan sharper. The main changes are:

1. Use compact region and boundary tables with biologically intuitive names.
2. Implement the four-state posterior as a vectorized log-likelihood tensor.
3. Collapse the eight v4 phases into four build phases plus a benchmark/refit phase.
4. Keep mathematical primitives stateless, immutable, and testable.
5. Avoid hard zero modes and avoid exposing a forest of user-facing knobs.

---

## 0. Review Of The v4 Critiques

### 0.1 Boundary sweeps: accept, with a simpler region/boundary layout

The critique is correct. Regions are already sorted, non-overlapping genomic intervals. Sorted arrays are enough. The useful data model is a one-dimensional boundary table:

```text
region i       = non-overlapping contiguous genomic interval
boundary i     = boundary immediately before local region i
boundary i + 1 = boundary immediately after local region i
```

For a reference with `N` regions, the boundary array has `N + 1` slots. The first and last slots are terminal boundaries. They make indexing simple, but they are not real crossing opportunities and carry zero crossing evidence. Internal boundary slot `i` is the boundary between local regions `i - 1` and `i`.

The prior boundary-storage note in [edge_centric_model.md](docs/edgecentric/edge_centric_model.md) makes the key point: today's left/right boundary slots are not duplicate data. They are complementary masses from the same crossing events. v5 therefore does **not** collapse them into one number. It exposes one conceptual boundary carrying two side-specific masses: mass landing in the region to the left and mass landing in the region to the right.

First implementation: pure Python `BoundaryTable` derived from existing `RegionArrays` + `RegionCountLedger`. No scanner schema change yet. A C++ accumulator migration can happen later if the boundary API proves itself.

### 0.2 Log-likelihood tensor: accept

The critique is exactly right. The four state classes should not become four branches of bespoke code. We should build orthogonal likelihood terms and combine them into:

```python
log_likelihoods = np.zeros((n_regions, 4), dtype=np.float64)
```

Then normalize once:

```python
log_p_states = log_likelihoods - logsumexp(log_likelihoods, axis=1, keepdims=True)
p_states = np.exp(log_p_states).astype(np.float32)
```

This keeps the E-step readable and gives us a natural place to inspect every evidence term.

### 0.3 Phase grouping: accept, with one correction

The proposed grouping into four build phases is better than v4's eight phases. The only correction is that [pipeline.py](src/rigel/pipeline.py) already implements `quant_from_buffer`. Phase IV should replace the current `calibration.fused_region_gdna` and `region_exposure` conduit with `RegionCalibration`, not remove a `NotImplementedError`.

### 0.4 Code aesthetics: accept, with practical boundaries

Rules for implementation:

- Mathematical objects are frozen dataclasses or stateless functions.
- Hot likelihood calculations are vectorized across regions.
- Loops over the four states and over reference slices are acceptable.
- Python loops over regions are not allowed in the calibration hot path unless isolated behind an exact small-count helper with tests and benchmarks.
- Diagnostic flags are explicit bit fields, not hidden side effects.

---

## 1. Core Statistical Contract

For each region `r`, estimate four posterior probabilities:

```text
state 0: background             = P(not expressed, not captured | data)
state 1: gdna_only_capture      = P(not expressed, captured     | data)
state 2: expressed_capture      = P(expressed,     captured     | data)
state 3: expressed_offtarget    = P(expressed,     not captured | data)
```

The state order is fixed and should live in code as constants:

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

`A_r` remains a denominator scalar:

```text
A_r = E[lambda_gdna(r) | data] / E[rho_off | data]
```

but the numerator is now a gDNA-only local density posterior informed by strand deconvolution and boundary imputation. It is not the old raw density posterior over gDNA + RNA mixture counts.

---

## 2. Data Model: Regions And Boundaries

### 2.1 Regions keep contained evidence

A region is one existing sorted region row:

```text
region_id = row in RegionArrays / DensityObservation
fields    = ref_id, start, end, signature, ts_class, length, contained counts
```

Regions keep contained evidence:

```text
contained_unspliced_pos
contained_unspliced_neg
contained_spliced_pos
contained_spliced_neg
contained_leff_gdna
```

### 2.2 Boundaries keep crossing evidence

A boundary slot exists before and after each region. For a reference with `N` local regions, there are `N + 1` local boundary slots:

```text
boundary 0       terminal boundary before region 0
boundary i       internal boundary between region i - 1 and region i
boundary N       terminal boundary after region N - 1
```

This gives a simple region-to-boundary relationship:

```text
left_boundary(region i)  = boundary i
right_boundary(region i) = boundary i + 1
```

In global arrays, the same relationship is recovered from per-reference offsets:

```text
local_i      = region_id - ref_region_offsets[ref_id]
left_b       = ref_boundary_offsets[ref_id] + local_i
right_b      = ref_boundary_offsets[ref_id] + local_i + 1
```

The boundary table does not need to store `left_region` or `right_region` index arrays. Those indices are derivable from offsets. This keeps the structure small and makes the indexing invariant obvious.

A pure Python `BoundaryTable` should expose:

```python
@dataclass(frozen=True, slots=True)
class BoundaryTable:
    boundary_pos: np.ndarray                  # int64[B]
    ref_region_offsets: np.ndarray            # int64[n_ref + 1]
    ref_boundary_offsets: np.ndarray          # int64[n_ref + 1]

    # Mass landing in the region to the left of this boundary.
    # Terminal left-of-reference boundaries carry zeros here.
    left_region_unspliced_pos: np.ndarray     # float32[B]
    left_region_unspliced_neg: np.ndarray
    left_region_spliced_pos: np.ndarray
    left_region_spliced_neg: np.ndarray

    # Mass landing in the region to the right of this boundary.
    # Terminal right-of-reference boundaries carry zeros here.
    right_region_unspliced_pos: np.ndarray    # float32[B]
    right_region_unspliced_neg: np.ndarray
    right_region_spliced_pos: np.ndarray
    right_region_spliced_neg: np.ndarray

    left_region_boundary_leff: np.ndarray     # opportunity for left region side
    right_region_boundary_leff: np.ndarray    # opportunity for right region side
    is_terminal: np.ndarray                   # bool[B]
```

The table is derived mechanically from the current ledger:

```text
internal boundary b between local regions i - 1 and i
left_region_*[b]  = region i - 1 boundary_right_* channels
right_region_*[b] = region i     boundary_left_* channels
```

Terminal boundaries are zero-filled and flagged. This preserves all current information while making boundary sweeps natural:

```text
region -> boundary -> region -> boundary -> region
```

### 2.3 Why not change the scanner now?

The current scanner/accumulator already records the needed information. A Python `BoundaryTable` gives us the API we want without C++ churn. Once boundary methods stabilize, a C++ boundary-centric accumulator can replace the backend without changing the calibration modules.

---

## 3. Evidence Products

### 3.1 Compartment-aware strand deconvolution

The current aggregate strand deconvolution is not enough for v5. We need contained and boundary channels separately:

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

For strand-specific libraries:

- contained gDNA estimates feed the regional density likelihood;
- boundary gDNA estimates feed the boundary imputation model;
- boundary RNA lower bounds prevent nascent RNA boundary flux from training the gDNA boundary model;
- regions with contained or boundary RNA lower bounds can be removed from background training.

For unstranded libraries, these fields are absent and the model falls back to gDNA-compatible counts with wider uncertainty.

### 3.2 Background off-target density

The background model estimates `rho_off`, the off-target gDNA density.

Initial background candidates:

```text
(intergenic OR intron-only)
AND no spliced mass
AND no strand-RNA evidence when stranded
AND not in the initial top-T contained-density exclusion set
```

The top-T exclusion set is computed once from initial intergenic/intronic candidates. It is a union exclusion, not a compounding filter.

Fit target:

```text
D_bg(r) ~ Poisson(rho_off * L_contained(r))
rho_off ~ Gamma(alpha_floor, beta_floor)
```

Use strand-deconvolved contained gDNA counts when available. Use unspliced contained counts in unstranded mode.

No exact zero mode. The posterior mean can become very small, but code should not branch on `rho_off == 0` to disable gDNA globally.

### 3.3 Boundary imputation

For each region side, boundary evidence estimates a local gDNA density and predicts contained gDNA:

```text
lambda_side | B_side ~ Gamma(alpha_prior + B_side, beta_prior + L_side)
D_contained_from_side ~ NegativeBinomial(
    r = alpha_prior + B_side,
    p = (beta_prior + L_side) / (beta_prior + L_side + L_contained)
)
```

The model returns:

```text
boundary_mean_gdna(r)
boundary_upper_gdna(r)
boundary_lambda_mean(r)
boundary_lambda_var(r)
p_boundary_captured(r)
```

In stranded data, `B_side` is the strand-deconvolved boundary gDNA count. In unstranded data, `B_side` is unspliced boundary mass, treated as gDNA-compatible rather than definitive gDNA.

### 3.4 Captured enrichment shrinkage

Per-region capture enrichment is noisy. Estimate a global captured component and shrink local values toward it:

```text
gamma_global = captured density / rho_off
log_gamma_r  = shrink(log(raw_gamma_r), log(gamma_global), precision_r)
A_r          = E[lambda_gdna(r)] / E[rho_off]
```

The exact shrinkage distribution can be decided in Phase II. The first candidate should be Gamma-Poisson empirical Bayes because it keeps conjugacy with the boundary count model. A log-normal shrinker is acceptable if the NB predictive tests are cleaner.

---

## 4. Boundary Sweeps As Vectorized 1D Message Passing

The sweep exists to move boundary-derived gDNA information into internal ambiguous regions such as the T1/T2 giant-exon case.

### 4.1 Message representation

Represent messages as excess Gamma evidence over the state prior:

```text
alpha_excess >= 0
beta_excess  >= 0
```

A local boundary posterior has:

```text
alpha_local = alpha_prior + B_boundary
beta_local  = beta_prior + L_boundary

alpha_excess_local = alpha_local - alpha_prior
beta_excess_local  = beta_local  - beta_prior
```

Combining evidence is addition in excess-parameter space:

```text
alpha_combined = alpha_prior + alpha_excess_a + alpha_excess_b
beta_combined  = beta_prior  + beta_excess_a  + beta_excess_b
```

### 4.2 Boundary transfer weights

Each internal boundary has a gDNA compatibility weight `w_boundary in [0, 1]`.

Stranded:

```text
w_boundary = deconvolved_boundary_gdna / max(total_boundary_compatible, eps)
```

Unstranded:

```text
w_boundary = unspliced_boundary / (unspliced_boundary + spliced_boundary + prior_mass)
```

This weight is not a hard truth. It controls how much density evidence crosses the boundary. Sparse or spliced-heavy boundaries naturally transmit less information.

### 4.3 Forward recurrence

For regions in one reference slice, the forward recurrence is:

```text
f[0] = local[0]
f[i] = local[i] + w_boundary[i] * f[i - 1]       for i >= 1
```

Here `w_boundary[i]` is the internal boundary between local region `i - 1` and local region `i`. Terminal boundary slots `0` and `N` are not used for transfer.

This recurrence can be vectorized per reference slice. Let:

```text
cp[i] = product_{k = 1..i} w_boundary[k]
f[i]  = cp[i] * cumulative_sum(local[i] / cp[i])
```

In implementation, use float64, reset at boundaries where `w_boundary` is effectively zero, and apply the formula separately for `alpha_excess` and `beta_excess`. The reverse pass is the same operation on reversed region arrays, using the boundary slots between each adjacent pair.

This is the important simplification: each reference can be handled as a plain one-dimensional array, without a per-region Python loop for the common case.

### 4.4 Combining local, forward, and reverse evidence

For each region:

```text
alpha_excess_total = alpha_local + alpha_forward + alpha_reverse
beta_excess_total  = beta_local  + beta_forward  + beta_reverse
```

then:

```text
lambda_sweep_mean = (alpha_prior + alpha_excess_total) / (beta_prior + beta_excess_total)
```

Convert this density posterior into contained gDNA predictive mean and upper tail with the same NB predictive distribution used by the local boundary model.

### 4.5 Scope of Phase II sweep

The first sweep should be deliberately simple:

- one-dimensional per reference and strand-compatible slice;
- only region and boundary tables;
- no transcript-aware sweep;
- no C++ implementation;
- no user-facing knobs beyond the existing confidence level.

The T1/T2 giant-exon fixture is the acceptance test.

---

## 5. Four-State Log-Likelihood Tensor

### 5.1 State templates

Use boolean state templates rather than class-specific branches:

```python
STATE_IS_EXPRESSED = np.array([False, False, True, True])
STATE_IS_CAPTURED = np.array([False, True, True, False])
```

Build evidence in two-column or four-column tensors:

```text
logL_expression      shape [R, 2]    columns: not_expressed, expressed
logL_capture         shape [R, 2]    columns: not_captured, captured
logL_gdna_density    shape [R, 2]    columns: offtarget_density, captured_density
logL_strand          shape [R, 4]    state-specific when strand evidence exists
log_prior_states     shape [R, 4]
```

Then assemble:

```python
log_likelihoods = log_prior_states.copy()
log_likelihoods += logL_expression[:, STATE_IS_EXPRESSED]
log_likelihoods += logL_capture[:, STATE_IS_CAPTURED]
log_likelihoods += logL_gdna_density[:, STATE_IS_CAPTURED]
log_likelihoods += logL_strand

log_p_states = log_likelihoods - logsumexp(log_likelihoods, axis=1, keepdims=True)
p_states = np.exp(log_p_states).astype(np.float32)
```

### 5.2 Orthogonal likelihood terms

Expression evidence:

```text
spliced mass                         penalizes not-expressed states
strand RNA lower bound               penalizes not-expressed states when stranded
contained excess over gDNA upper      supports expressed states
```

Capture evidence:

```text
boundary mixture posterior            supports captured vs not captured
shrunk local gamma                     supports captured density state
```

gDNA density evidence:

```text
NB predictive likelihood under rho_off
NB predictive likelihood under captured density
strand-deconvolved contained gDNA likelihood when stranded
```

Strand contrast evidence:

```text
BetaBinomial(k_sense | D_gdna, kappa_d) for gDNA component
Binomial or near-binomial RNA strand term for RNA component
```

The first implementation can use the existing exact/approx strand posterior helpers where needed, but the tensor assembly should remain the public shape of the solver.

### 5.3 State posterior outputs

The solver returns:

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
    flags: np.ndarray                 # uint16[R]
```

Convenience properties can expose:

```text
p_background
p_gdna_only_capture
p_expressed_capture
p_expressed_offtarget
p_expressed
p_captured
```

---

## 6. Calibration Iteration

The iteration is a small calibration EM, not the transcript EM.

### 6.1 E-step

Given current `rho_off`, `kappa_d`, boundary model, and gamma shrinker:

```text
build log-likelihood tensor
normalize to p_states
compute mu_gdna, upper_gdna, rna_lower, A_r
```

### 6.2 M-step

Use posterior state weights to refit model components:

```text
background weight       = p_background
not_expressed weight    = p_background + p_gdna_only_capture
captured weight         = p_gdna_only_capture + p_expressed_capture
offtarget weight        = p_background + p_expressed_offtarget
```

Refit:

```text
rho_off       from background-weighted gDNA evidence
kappa_d       from high-confidence not-expressed regions when stranded
boundary mix  from boundary evidence weighted by captured/offtarget probability
gamma_global  from captured-weighted gDNA density
```

### 6.3 Stop rule

Keep this simple in v5:

```text
max_calibration_passes = 5
stop when max absolute state-probability shift and relative rho_off shift are small
```

The exact tolerances should be internal defaults at first. Report pass-level diagnostics. Do not expose many knobs until benchmarks show a real need.

---

## 7. Downstream EM Contract

### 7.1 Numerator

Use `RegionCalibration.mu_gdna` and `RegionCalibration.upper_gdna` in prior assembly.

Boundary evidence is real gDNA-compatible evidence. Do not disable the gDNA component in a locus that has boundary evidence or positive `upper_gdna` simply because a hard state label is uncertain.

Longer term, the cleanest behavior is to always include a gDNA state with a tiny smooth prior. If the short-term code still needs an enable flag, use evidence, not a hard mode switch:

```text
enable_gdna = any(boundary evidence OR strand-deconvolved gDNA evidence OR upper_gdna > 0)
```

### 7.2 Denominator

Keep the agreed v3/v4 denominator for the first implementation:

```text
gdna_eff_len(locus) = unweighted_gdna_eff_len(locus) * bp_weighted_mean(A_r)
```

The key v5 improvement is that `A_r` is now shrunk and model-derived before it reaches this denominator. Short noisy regions should not inject raw boundary ratios into locus effective length.

The FL-convolved exposure integral remains a deferred improvement if benchmark evidence requires it.

### 7.3 Current pipeline wiring

[pipeline.py](src/rigel/pipeline.py) currently requires `calibration.fused_region_gdna`. Phase IV replaces that conduit with `calibration.region_calibration` while keeping a compatibility accessor only if needed for incremental migration.

---

## 8. Build Phases

The plan has four build phases plus a benchmark/refit phase.

```text
Phase I  -> Phase II -> Phase III -> Phase IV -> Phase V
strand      boundary    four-state   downstream   benchmark
background  imputation  solver       wiring       refit
```

### Phase I - Compartment strand and background

Scope:

- cleanup the old `400` precision cap and old density exposure path;
- add `BoundaryTable` as a derived Python table over existing region/ledger arrays;
- extend strand deconvolution to contained, boundary-left, and boundary-right compartments;
- create `background_model.py`;
- fit `rho_off` with top-T initial exclusion and strand-RNA exclusion when available.

Primary files:

- [strand_deconv.py](src/rigel/calibration/strand_deconv.py)
- `src/rigel/calibration/boundaries.py`
- `src/rigel/calibration/background_model.py`
- [density_observation.py](src/rigel/calibration/density_observation.py)

Verification gate:

- `tests/test_boundaries.py`: boundary reconstruction equals current left/right boundary ledger slots.
- `tests/test_compartment_strand_deconv.py`: compartment sums agree with aggregate behavior where expected.
- `tests/test_background_model.py`: capture and non-capture fixtures recover `rho_off` and exclude strand-RNA-positive introns.

### Phase II - Boundary imputation and contiguous sweeps

Scope:

- create `boundary_model.py` for Gamma-Poisson / NB boundary-to-contained prediction;
- fit off-target and captured boundary components;
- shrink `gamma_r` toward a global captured enrichment;
- create `boundary_sweep.py` implementing vectorized forward/backward scans over `BoundaryTable` reference slices;
- return local, forward, reverse, and combined boundary-imputed gDNA estimates.

Primary files:

- `src/rigel/calibration/boundary_model.py`
- `src/rigel/calibration/boundary_sweep.py`
- `tests/test_boundary_model.py`
- `tests/test_boundary_sweep.py`

Verification gate:

- sparse non-capture boundaries shrink to `rho_off`;
- high captured boundaries produce high `A_r` with shrinkage;
- no-gDNA libraries produce tiny but nonzero smooth priors;
- T1/T2 giant-exon fixture propagates outer boundary evidence inward with uncertainty decay.

### Phase III - Four-state tensor solver

Scope:

- implement `latent_states.py` around the `[R, 4]` log-likelihood tensor;
- implement vectorized expression, capture, gDNA density, and strand likelihood terms;
- implement calibration E/M iteration with posterior state weights;
- output `RegionCalibration`;
- refit `kappa_d` after each pass from high-confidence not-expressed regions when stranded.

Primary files:

- `src/rigel/calibration/latent_states.py`
- `src/rigel/calibration/_result.py`
- `tests/test_latent_states.py`
- `tests/test_calibration_iteration.py`

Verification gate:

- synthetic four-category fixture converges within five passes;
- expressed captured and expressed off-target regions are usable in stranded data;
- unstranded expressed classes remain uncertain rather than falsely zeroed;
- state probabilities sum to one and diagnostics identify ambiguous regions.

### Phase IV - Downstream wiring and diagnostics

Scope:

- rewire `_orchestrator.py`:

```text
FL models
-> strand summary/deconvolution when available
-> BoundaryTable
-> background model
-> boundary model
-> boundary sweeps
-> four-state calibration iteration
-> RegionCalibration
-> exposure/prior assembly
```

- rewrite `exposure.py` to read `RegionCalibration.A_r`;
- rewrite `prior.py` to read `RegionCalibration.mu_gdna` and `upper_gdna`;
- update [pipeline.py](src/rigel/pipeline.py) to require `RegionCalibration` instead of `fused_region_gdna`;
- add compact transcript/gene/locus diagnostics;
- add optional region-level calibration debug output.

Primary files:

- [src/rigel/calibration/_orchestrator.py](src/rigel/calibration/_orchestrator.py)
- [exposure.py](src/rigel/calibration/exposure.py)
- [prior.py](src/rigel/calibration/prior.py)
- [pipeline.py](src/rigel/pipeline.py)

Verification gate:

- capture-off golden outputs match within tolerance;
- capture-on no-gDNA does not invent substantial gDNA;
- capture-on high-gDNA reduces gDNA-to-RNA leakage;
- summary diagnostics include state masses by pass.

### Phase V - Benchmark and calibration refit

Scope:

- emit post-capture truth from simulation: transcript abundances and FL distributions after capture;
- benchmark four regimes:
  - unstranded, no capture;
  - unstranded, capture;
  - strand-specific, no capture;
  - strand-specific, capture;
- tune internal defaults only when benchmark evidence demands it;
- decide whether FL-convolved exposure integral is needed.

Verification gate:

- benchmarks compare against post-capture truth only;
- the T1/T2 fixture demonstrates boundary-sweep benefit;
- stranded capture is the first acceptance target;
- unstranded capture is reported with honest uncertainty.

---

## 9. Maintainability Rules

### 9.1 Statelessness

Model builders return immutable dataclasses. They do not mutate global state or hide convergence arrays inside objects.

Good:

```python
background = fit_background_model(observation, strand_channels, config)
```

Avoid:

```python
model.fit(...)
model.update_pass(...)
```

### 9.2 Explicit vectorization

All hot terms are arrays:

```text
[R]
[R, 2]
[R, 4]
[B]
```

Use `scipy.special`, `scipy.stats` vectorized functions, and NumPy scans. Do not write region-wise Python loops in likelihood assembly.

### 9.3 Transparent diagnostics

Every important fallback or approximation gets a flag:

```text
FLAG_LOW_BOUNDARY_OPPORTUNITY
FLAG_STRAND_UNINFORMATIVE
FLAG_BOUNDARY_SPARSE
FLAG_SWEEP_DAMPED
FLAG_STATE_AMBIGUOUS
FLAG_EXACT_POSTERIOR
FLAG_APPROX_POSTERIOR
```

Flags are bitwise arrays in `RegionCalibration.flags` and in intermediate debug tables.

### 9.4 Few user-facing knobs

Expose only:

```text
gdna_confidence
background_trim_fraction
max_calibration_passes
```

Everything else starts as an internal default or a named profile. Do not expose gamma shrinkage hyperparameters, sweep transfer priors, or tensor penalties unless benchmarks force it.

---

## 10. Diagnostics

Compact user-facing outputs:

```text
quant.feather:       mean_A_r, n_regions_captured, n_regions_expressed
gene_quant.feather:  mean_A_r, n_regions_captured, n_regions_expressed
loci.feather:        mean_A_r, n_regions_captured, state_mass summaries
summary.json:        pass-level calibration diagnostics
```

Optional region debug output:

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
boundary_forward_mean
boundary_reverse_mean
boundary_combined_mean
strand_contained_gdna_mean
strand_boundary_left_gdna_mean
strand_boundary_right_gdna_mean
flags
```

---

## 11. Residual Issues

### R1. Boundary sweeps are mathematically elegant but need empirical damping

The vectorized recurrence can propagate errors if a boundary transfer weight is too permissive. Start with conservative transfer weights and report sweep contribution separately from local boundary contribution.

### R2. Gamma shrinkage distribution remains a choice

The plan allows Gamma-Poisson EB or log-normal shrinkage. Phase II should pick the simpler one that passes NB predictive tests. Do not implement both unless benchmarks require it.

### R3. Unstranded capture remains weakly identifiable

v5 predicts gDNA rather than zeroing it, but all-captured/all-expressed unstranded regions can still be broad. Report uncertainty rather than pretending the state is resolved.

### R4. Vectorized exact posterior may still need a small-count helper

The tensor E-step should be vectorized. If exact enumeration over possible gDNA counts is still needed for small counts, quarantine it in a helper, cap it, flag it, and benchmark it.

### R5. BoundaryTable is a transitional API

The first implementation derives boundaries from current region-side slots. If the design proves stable, a later scanner migration can make boundaries native. The calibration code should not care which backend is used.

### R6. Mappability is still deferred

The model should be written so mappability-corrected effective length can replace raw gDNA opportunity later, but v5 does not block on it.

---

## 12. Short Version

v5 keeps v4's four biological state classes, but turns the plan into an implementable architecture:

```text
Regions + boundaries
-> compartment-aware strand deconvolution
-> background rho_off
-> boundary-to-contained gDNA imputation
-> vectorized forward/backward boundary sweeps
-> [R, 4] log-likelihood tensor
-> RegionCalibration
-> prior/exposure wiring into transcript EM
```

The plan is ambitious, but now the ambition is shaped into testable modules. The first implementation can stay pure Python/NumPy/SciPy, with no new C++ scanner work, while still giving us intuitive region/boundary semantics and the joint expression x capture posterior we need.
