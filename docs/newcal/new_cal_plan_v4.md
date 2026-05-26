# New Calibration Plan v4

Date: 2026-05-26
Status: rewrite after v3 review
Supersedes: `new_cal_plan_v3.md`, `new_cal_plan_v2.md`, `new_cal_plan_v1.md`
Primary target: strand-specific hybrid capture, while preserving a path for non-capture and unstranded data.

This is a rewrite, not a patch. The v3 document moved in the right direction by separating expression from capture, but it introduced a misleading three-pool abstraction and treated boundary flux as a one-shot enrichment test. v4 puts the four expression x capture state classes at the center and promotes boundary-to-contained gDNA imputation into a first-class model.

---

## 0. Critique Of v3

v3 has several major problems:

1. It introduced three pools named `A`, `B`, and `C`. That creates an extra layer of abstraction and hides the biology. The natural model has four named state classes: not expressed/not captured, not expressed/captured, expressed/captured, and expressed/not captured.
2. It said expressed captured regions are "never used to fit any parameter". That is false for strand-specific data. If the strand model can deconvolve gDNA and RNA in an expressed exon, that region can inform the captured gDNA surface.
3. It said expressed off-target regions are not a pool. That is also false. In strand-specific data, expressed off-target regions are useful because their deconvolved gDNA density should stay near `rho_off`.
4. It made the captured-not-expressed state too dependent on strand-specific exon screening. In unstranded data, captured-not-expressed regions can still be seeded by unspliced boundary evidence plus no spliced evidence, if we have a boundary imputation model.
5. It treated boundary flux mostly as a hard capture seed. Boundary flux must instead be modeled as local evidence for the gDNA density inside neighboring regions.
6. It proposed setting unidentifiable unstranded loci to no gDNA. That is wrong. If a region has boundary evidence, we must predict its gDNA from the appropriate density model, even when uncertainty is large.
7. It introduced too many knobs and too many hard mode switches. v4 should glide smoothly toward low-gDNA behavior through weak posterior mass, not through a magic zero trigger.
8. It did not fully use strand-specific data. Strand deconvolution should apply separately to contained, left-boundary, and right-boundary compartments, so density models can be fit on the gDNA portion of each mixture.

The right architecture is therefore:

```text
observed 12-channel ledger
  -> optional strand deconvolution by compartment
  -> background off-target density model
  -> boundary imputation model
  -> four-state posterior over expression x capture
  -> regional gDNA numerator and exposure denominator
  -> locus EM priors
```

---

## 1. Core Contract

For every region `r`, calibration estimates four state probabilities:

```text
p_background(r)          = P(not expressed, not captured | data)
p_gdna_only_capture(r)   = P(not expressed, captured     | data)
p_expressed_capture(r)   = P(expressed,     captured     | data)
p_expressed_offtarget(r) = P(expressed,     not captured | data)
```

These replace the v3 three-pool model. They are not hard labels. They are posterior probabilities that can be used as training weights during iteration and reported as diagnostics.

The two marginal probabilities are derived from those four values:

```text
p_expressed(r) = p_expressed_capture(r) + p_expressed_offtarget(r)
p_captured(r) = p_gdna_only_capture(r) + p_expressed_capture(r)
```

The downstream EM still only needs two regional products:

```text
mu_gdna(r)       expected gDNA fragments in region r
upper_gdna(r)    conservative upper credible bound on gDNA fragments in region r
A_r              gDNA exposure multiplier for region r
```

`A_r` is an exposure denominator term. It is computed from the gDNA density posterior, not from the old `rho_post / rho_ref` path that was dominated by the arbitrary `400` precision cap.

A useful interpretation is:

```text
A_r = E[lambda_gdna(r) | data] / E[rho_off | data]
```

This looks like a density ratio, but it is not the old density-ratio implementation. The numerator is a gDNA-only local density posterior informed by boundary imputation and strand deconvolution. It is not a posterior over the raw contained mixture of gDNA + RNA.

---

## 2. Evidence Products

### 2.1 The 12-channel ledger is the raw source

The existing regional accumulator already gives the right raw material:

```text
3 compartments: contained, boundary_left, boundary_right
2 splice states: unspliced, spliced
2 observed strands: pos, neg
```

The current helpers in [region_count_ledger.py](src/rigel/calibration/region_count_ledger.py) expose the needed totals. v4 should continue to keep these channels separate as long as possible.

### 2.2 Strand-specific data should be used by compartment

The current `RegionGdnaEstimate` in [strand_deconv.py](src/rigel/calibration/strand_deconv.py) summarizes a per-region gDNA posterior after aggregating counts. v4 needs a compartment-aware version:

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

This is essential. In strand-specific data, we are not forced to assume all boundary unspliced flux is gDNA. We can deconvolve contained and boundary channels separately. An intron with nascent RNA in either contained or boundary channels can be marked expressed and excluded from the background seed set.

### 2.3 Unstranded data has weaker evidence, not no evidence

In unstranded data:

- spliced counts are definite RNA evidence;
- unspliced boundary counts are gDNA-compatible evidence;
- contained unspliced counts are a mixture of gDNA, nascent RNA, and unspliced RNA;
- expressed captured regions cannot be cleanly seeded from strand balance.

The unstranded path therefore leans heavily on the boundary imputation model. It should not set gDNA to zero just because strand evidence is unavailable.

---

## 3. Four Named State Classes

These are state classes, not hard pools. The names should appear in diagnostics and tests.

| State class | Meaning | Primary role | Strand-specific seeding | Unstranded seeding |
|---|---|---|---|---|
| `background` | not expressed, not captured | seeds `rho_off` | intergenic/intronic, no spliced mass, no strand-RNA evidence, not in initial top-T contained-density exclusion set | intergenic/intronic, no spliced mass, not in initial top-T contained-density exclusion set |
| `gdna_only_capture` | not expressed, captured | seeds/validates captured gDNA density and `gamma` | no spliced mass; strand test consistent with pure gDNA; deconvolved gDNA density above `rho_off` | no spliced boundary mass; nonzero unspliced boundary evidence; boundary-imputed contained gDNA density above `rho_off` |
| `expressed_capture` | expressed, captured | informs captured gDNA surface when strand deconvolution is available | RNA evidence present; deconvolved gDNA density above `rho_off` | not seedable; posterior may still assign mass from boundary model |
| `expressed_offtarget` | expressed, not captured | validates off-target gDNA surface in RNA-rich regions | RNA evidence present; deconvolved gDNA density near `rho_off` | not seedable; posterior may still assign mass from spliced RNA plus weak boundary evidence |

Important consequences:

- There are four state classes because expression and capture are crossed factors.
- Strand-specific data expands the useful training set. It lets us use RNA-rich regions because gDNA and RNA can be deconvolved.
- Unstranded data cannot seed the expressed state classes with gDNA density, but it can still infer them probabilistically.
- The background state is not "all introns/intergenic forever". It is an iterative posterior state.

---

## 4. Background Off-Target Density Model

The background model estimates the depleted off-target gDNA density `rho_off`.

### 4.1 Initial background candidates

Initial candidates are regions satisfying:

```text
(intergenic OR intron-only)
AND no spliced mass
AND no strand-RNA evidence if strand-specific data is available
AND not in the initial top-T contained-density exclusion set
```

The top-T exclusion set is computed once from the initial intergenic/intronic candidates. It is a union-style exclusion, not a compounding filter. A region is excluded if it is in the initial top-T density tail OR if strand evidence later says it has RNA.

The intent of the top-T exclusion is to protect against:

- unannotated transcripts;
- high nascent RNA introns;
- negative-control probes deliberately placed in intergenic/intronic space;
- anomalous regions with strong local enrichment.

### 4.2 Fit target

Fit a Gamma-Poisson density model on gDNA-compatible counts and gDNA effective length:

```text
D_bg(r) ~ Poisson(rho_off * L_bg(r))
rho_off ~ Gamma(alpha_off, beta_off)
```

For strand-specific data, `D_bg(r)` should be the strand-deconvolved gDNA count, not the raw unspliced total. For unstranded data, `D_bg(r)` is the unspliced contained count after excluding any region with spliced RNA evidence.

### 4.3 No hard zero mode

`rho_off` should never become an exact branch condition. Use a weak, interpretable prior so the posterior mean can become very small but nonzero:

```text
E[rho_off | data] = (alpha_floor + sum D_bg) / (beta_floor + sum L_bg)
```

A Jeffreys-style half-fragment prior is reasonable for the first implementation. The important rule is that no code path should say "if rho_off == 0, disable gDNA everywhere". Low-gDNA libraries should glide toward tiny priors and broad uncertainty.

---

## 5. Boundary Imputation Model

This is the center of v4.

Boundary flux is the only broadly available local evidence that can predict gDNA inside exons in unstranded data. It remains valuable in strand-specific data because it gives an exposure-local signal orthogonal to contained counts.

### 5.1 Boundary observations

For each region `r`, keep left and right sides separate:

```text
B_left_gdna(r)     gDNA-compatible or deconvolved gDNA count crossing left boundary
B_right_gdna(r)    gDNA-compatible or deconvolved gDNA count crossing right boundary
L_left(r)          left boundary effective length
L_right(r)         right boundary effective length
L_contained(r)     contained effective length
```

In strand-specific data:

```text
B_side_gdna = strand-deconvolved gDNA boundary count for that side
B_side_rna_lower = strand-deconvolved RNA lower bound for that side
```

In unstranded data:

```text
B_side_gdna = boundary_side_unspliced_count
B_side_rna_lower = boundary_side_spliced_count
```

The unstranded assignment is deliberately conservative. It treats unspliced boundary flux as gDNA-compatible evidence, not as definitive gDNA truth.

### 5.2 Boundary-to-contained prediction

Given a side boundary count `B` and opportunity `L_b`, infer a local gDNA density `lambda_side` and predict contained gDNA inside the receiving region:

```text
lambda_side | B ~ Gamma(alpha_prior + B, beta_prior + L_b)
D_contained_from_side ~ NegativeBinomial(
    r = alpha_prior + B,
    p = (beta_prior + L_b) / (beta_prior + L_b + L_contained)
)
```

This gives both:

```text
mean_contained_gdna_from_boundary
upper_contained_gdna_from_boundary(confidence)
```

The upper tail is critical. It lets us ask: "How much contained mass could plausibly be gDNA, given the boundary flux?"

### 5.3 Off-target and captured boundary components

Hybrid capture makes boundary densities a mixture:

```text
off-target boundary component: lambda ~ Gamma(alpha_off, beta_off)
captured boundary component:   lambda ~ Gamma(alpha_cap, beta_cap)
```

The captured component can be parameterized as enrichment over background:

```text
lambda_cap = gamma * rho_off
```

`gamma` should be shrunk. Per-region boundary counts are noisy, especially for short exons and sparse boundaries. The first implementation should estimate:

```text
global_gamma     from high-confidence captured boundary observations
local_gamma_r    shrunk toward global_gamma with strength determined by boundary opportunity
```

This avoids propagating a single noisy boundary ratio directly into `A_r`.

### 5.4 Capture probability from boundary evidence

For each side, compute:

```text
P(side captured | B_side, L_side)
```

using the off-target and captured Gamma-Poisson predictive likelihoods. Then combine left and right sides to get a regional boundary-derived capture probability.

This replaces v3's hard Poisson-tail capture seed. A diagnostic p-value can still be reported, but the model should use posterior probabilities.

### 5.5 Combining boundary sides

Combine left and right boundary predictions as posterior evidence, not by taking a max:

```text
posterior_precision = precision_left + precision_right + precision_prior
mean = weighted posterior mean
upper = posterior predictive quantile
```

If one side is underpowered because the region is shorter than the fragment-length opportunity, the other side can still contribute. If both sides are sparse, the posterior shrinks toward the state-specific density prior.

---

## 6. Internal Exon gDNA Imputation By Sweeps

The user's T1/T2 example makes sense, but it should be implemented as probabilistic message passing rather than deterministic copying.

Example:

```text
T1+: (1000, 2000), (5000, 6000), (10000, 15000)
T2+: (1000, 11000)  # one giant exon
```

A giant exon can make internal regions ambiguous. Boundary evidence at the outer edges can still inform internal regions if we propagate local gDNA density through adjacent regions.

### 6.1 Region adjacency graph

Build an ordered adjacency graph per reference and strand-compatible component:

```text
node = region
edge = physical boundary between adjacent regions
edge observations = left/right boundary unspliced/spliced counts, strand split when available
```

This graph is independent of transcript IDs. It is geometry plus regional signatures.

### 6.2 Left-to-right and right-to-left passes

For each component:

```text
left-to-right pass:
  initialize message from the left external boundary
  for each region:
    combine incoming message with local boundary posterior
    store forward gDNA density posterior
    propagate through the right edge with uncertainty inflation

right-to-left pass:
  same procedure from the right external boundary
```

The message is not a single count. It is a distribution over local gDNA density or contained gDNA count:

```text
message = (mean_density, variance_density, effective_observations, flags)
```

### 6.3 Edge transfer rule

At an edge between regions `i` and `j`, estimate how much of the crossing flux is gDNA-compatible.

With strand-specific data:

```text
edge_gdna_fraction = deconvolved_gdna_boundary / total_boundary_compatible
```

With unstranded data:

```text
edge_gdna_fraction = unspliced_boundary / (unspliced_boundary + spliced_boundary + prior_mass)
```

This unstranded fraction is not proof of gDNA. It is a noisy local compatibility score. The uncertainty on the message should widen when counts are sparse or spliced mass is high.

### 6.4 Combining local and sweep evidence

For each region, combine three evidence sources:

```text
local boundary posterior
forward sweep posterior
reverse sweep posterior
```

using precision-weighted posterior combination or a product-of-experts approximation with variance floors.

This solves the internal-exon problem in the right spirit:

- outer exon boundaries anchor gDNA density;
- internal ambiguous regions inherit evidence with growing uncertainty;
- exon-exon boundaries are not discarded;
- spliced boundary mass reduces confidence in unstranded gDNA projection;
- strand-specific data can separate gDNA and RNA throughout the sweep.

This should live behind a separate module, likely `boundary_sweep.py`, so we can test it independently.

---

## 7. Four-State Inference Loop

The algorithm is iterative because the training sets improve as state probabilities improve.

### 7.1 Initialization

```text
1. Build RegionCountLedger and DensityObservation.
2. Train observed FL models as today.
3. If the library is informatively stranded:
     fit strand_summary
     build compartment-aware strand counts
     estimate kappa_d from initial gDNA-like regions
     deconvolve contained, boundary_left, and boundary_right channels
   Else:
     mark strand-derived evidence unavailable.
4. Fit initial background density from initial background candidates.
5. Fit initial boundary imputation model.
6. Run boundary sweeps to impute local gDNA density in ambiguous internal regions.
7. Initialize four state probabilities from the rules in Section 3.
```

### 7.2 Regional posterior update

For each region, evaluate the four state likelihoods:

```text
background:
  gDNA density near rho_off
  no RNA likelihood except weak nuisance mass

gdna_only_capture:
  captured gDNA density
  no RNA likelihood

expressed_capture:
  captured gDNA density
  RNA nuisance mass allowed

expressed_offtarget:
  off-target gDNA density
  RNA nuisance mass allowed
```

Evidence terms:

```text
spliced mass                  strong evidence for expressed states
strand RNA lower bound         strong evidence for expressed states when available
strand gDNA estimate           direct gDNA evidence when available
boundary imputed gDNA          local gDNA evidence in all modes
contained unspliced excess     RNA evidence after gDNA upper bound
```

The posterior returns:

```text
state probabilities
mu_gdna(r)
upper_gdna(r)
rna_lower(r)
A_r
local_gamma_r (shrunk)
flags
```

### 7.3 Training-set refresh

After each posterior update, refresh weighted training sets:

```text
background training weight       = p_background(r)
not-expressed training weight    = p_background(r) + p_gdna_only_capture(r)
captured training weight         = p_gdna_only_capture(r) + p_expressed_capture(r)
offtarget validation weight      = p_background(r) + p_expressed_offtarget(r)
```

For strand-specific data, `kappa_d` should be refit after each pass using high-confidence not-expressed regions. That means regions weighted by `p_background + p_gdna_only_capture`, with ambiguous strand classes excluded.

### 7.4 Iteration and damping

A practical first implementation:

```text
for pass in 1..max_passes:
  refit rho_off with background weights
  refit kappa_d with not-expressed weights if stranded
  refit boundary mixture and gamma shrinkage
  rerun boundary sweeps
  update four-state posteriors
  damp updates if state probabilities oscillate
  stop when rho_off, median A_r, and state probabilities stabilize
```

This does not need a complicated convergence theory. We should start with a small cap, e.g. 5 passes, report all pass-level diagnostics, and let benchmarks decide whether damping or a sharper stopping rule is needed.

---

## 8. Downstream EM Contract

### 8.1 Regional calibration table

Introduce one object that replaces the current density/fusion/exposure split:

```python
@dataclass(frozen=True, slots=True)
class RegionCalibration:
    p_background: np.ndarray
    p_gdna_only_capture: np.ndarray
    p_expressed_capture: np.ndarray
    p_expressed_offtarget: np.ndarray
    mu_gdna: np.ndarray
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    A_r: np.ndarray
    gamma_r: np.ndarray
    rho_off: float
    kappa_d: float | None
    flags: np.ndarray
```

`CalibrationResult` should carry this object plus the intermediate models for diagnostics.

### 8.2 Numerator

Locus-level gDNA prior count uses `mu_gdna` and `upper_gdna`, allocated by region geometry as today. Boundary evidence is objective evidence of gDNA-compatible fragments; it must not be discarded by an overzealous `enable_gdna` gate.

If the EM implementation still requires a boolean gDNA component switch in the short term, the safe rule is:

```text
enable gDNA when any region in the locus has boundary evidence, strand-deconvolved gDNA evidence,
or positive upper_gdna.
```

Longer term, the cleaner design is to always allow a gDNA component with a very small smooth prior when evidence is absent.

### 8.3 Denominator

`A_r` modulates the gDNA effective length. For v4, keep the bp-weighted regional exposure mean agreed in v3 review:

```text
gdna_eff_len(locus) = unweighted_gdna_eff_len(locus) * bp_weighted_mean(A_r over locus blocks)
```

The important v4 change is that `A_r` itself is shrunk and model-derived. Short, noisy regions should shrink toward the appropriate global captured or off-target density before they ever reach this denominator.

The FL-convolved exposure integral remains a deferred improvement if benchmarks show the bp-weighted approximation is insufficient.

---

## 9. Implementation Phases

### Phase 0 - Cleanup and naming

- Delete the `400` precision cap path: `DENSITY_PRIOR_MIN_CV`, `DENSITY_PRIOR_MAX_PRECISION`, `compute_beta_cap`, and related summary fields.
- Remove the old `rho_post / rho_ref` exposure path.
- Remove `density_max_exposure`.
- Collapse duplicate `boundary_left_leff` and `boundary_right_leff` fields into one side opportunity if no caller needs separate arrays.
- Avoid `Pool A/B/C` names. Use the four state class names everywhere.

### Phase 1 - Compartment-aware strand deconvolution

- Extend [strand_deconv.py](src/rigel/calibration/strand_deconv.py) to deconvolve contained, left-boundary, and right-boundary channels separately.
- Keep the existing aggregate deconvolution as a compatibility wrapper while tests migrate.
- Refit `kappa_d` from high-confidence not-expressed regions after each calibration pass.
- Add tests where intronic boundary flux is nascent RNA, not gDNA, and verify it is excluded from background training.

### Phase 2 - Background density model

- Create `background_model.py`.
- Fit `rho_off` from weighted background candidates.
- Implement initial top-T exclusion as a one-time mask.
- Use strand-deconvolved gDNA counts when available.
- Keep the zero-gDNA behavior smooth through a weak Gamma prior, not a hard branch.

### Phase 3 - Boundary imputation model

- Create `boundary_model.py`.
- Fit off-target and captured boundary density components.
- Predict contained gDNA from boundary evidence using the Gamma-Poisson / Negative Binomial predictive model.
- Estimate per-region `gamma_r` with shrinkage toward a global captured enrichment.
- Return mean and upper-tail contained gDNA predictions.
- Add unit tests for sparse boundaries, high boundaries, no-gDNA libraries, and short regions.

### Phase 4 - Boundary sweep imputation

- Create `boundary_sweep.py`.
- Build ordered region adjacency components.
- Implement left-to-right and right-to-left probabilistic message passes.
- Combine local, forward, and reverse boundary evidence.
- Test the T1/T2 giant-exon example as a dedicated fixture.

### Phase 5 - Four-state posterior

- Create or rewrite `latent_states.py` around the four state probabilities.
- Use boundary-imputed gDNA, strand-deconvolved gDNA, spliced evidence, and contained excess together.
- Output `RegionCalibration`.
- Add flags for state confidence and ambiguity.

### Phase 6 - Orchestrator and EM wiring

- Reorder `_orchestrator.py`: FL -> strand if available -> background -> boundary model -> boundary sweep -> four-state posterior -> exposure -> priors.
- Rewrite `exposure.py` to read `RegionCalibration.A_r`.
- Rewrite `prior.py` to read `RegionCalibration.mu_gdna` and `upper_gdna`.
- Preserve compatibility accessors for one minor version if needed.

### Phase 7 - Diagnostics, outputs, and docs

- Add region-level calibration diagnostics, at least in debug output.
- Add transcript/gene/locus diagnostics: `mean_A_r`, `n_regions_captured`, `n_regions_background`, `n_regions_expressed`, and state probability summaries.
- Add pass-level summary diagnostics for `rho_off`, `kappa_d`, `global_gamma`, state counts, and convergence.
- Update manual and parameter docs.

### Phase 8 - Benchmarking

- Emit post-capture truth from the simulator: transcript abundances and FL distributions after capture.
- Benchmark the four RNA-seq regimes:
  - unstranded, no capture;
  - unstranded, hybrid capture;
  - strand-specific, no capture;
  - strand-specific, hybrid capture.
- Prioritize strand-specific hybrid capture first.

---

## 10. Diagnostics To Expose

Region-level diagnostics should include:

```text
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
boundary_imputed_mean
boundary_imputed_upper
forward_sweep_mean
reverse_sweep_mean
strand_gdna_mean_contained
strand_gdna_mean_boundary_left
strand_gdna_mean_boundary_right
flags
```

Summary diagnostics should include:

```text
rho_off_mean
rho_off_interval
global_gamma_mean
global_gamma_interval
kappa_d
n_passes
state_probability_mass_by_pass
top_T_excluded_regions
strand_RNA_excluded_background_regions
boundary_model_component_weights
post_capture_truth_used_for_benchmarking
```

The user-facing outputs (`quant.feather`, `gene_quant.feather`, `loci.feather`) should stay compact. Put detailed region diagnostics behind a debug or calibration-report output.

---

## 11. Benchmarking Principles

Capture changes the library. Therefore benchmark capture data against post-capture truth, not against pre-capture simulator parameters.

The simulator must emit:

```text
truth_post_capture_transcript_abundance.tsv
truth_post_capture_gdna_fl.tsv
truth_post_capture_rna_fl.tsv
truth_region_capture_state.tsv
```

Acceptance should focus on:

- gDNA-to-RNA leakage in captured expressed regions;
- recovery of gDNA in captured not-expressed regions;
- preservation of low gDNA in off-target expressed regions;
- no regression in non-capture data;
- uncertainty behavior in unstranded capture;
- whether boundary sweeps improve internal exon gDNA prediction in the T1/T2 fixture.

---

## 12. Residual Issues After v4

### R1. Unstranded boundary flux still conflates gDNA and nascent RNA

In unstranded data, unspliced boundary flux is only gDNA-compatible. It can include nascent RNA. v4 handles this by reporting uncertainty and using spliced mass to weaken gDNA projection, but it cannot fully solve the ambiguity without strand evidence.

### R2. Boundary sweeps can propagate local errors

A bad boundary estimate can influence downstream internal regions. The sweep must carry variance and effective-observation counts so messages decay when evidence is weak. The T1/T2 fixture should include adversarial sparse-boundary cases.

### R3. Gamma shrinkage needs a concrete distribution

The plan says shrink `gamma_r` toward a global captured enrichment. The first implementation must choose a simple distribution, probably log-normal shrinkage or Gamma-Poisson empirical Bayes. This should be decided in Phase 3 with unit tests.

### R4. Top-T exclusion is still a heuristic

The initial top-T background exclusion is useful, but it is still a robustness heuristic. It should be computed once, reported, and kept out of the main probabilistic model as much as possible.

### R5. Compartment-aware strand deconvolution is new work

Current code aggregates strand counts across compartments. v4 needs per-compartment deconvolution, which is conceptually straightforward but test-heavy.

### R6. Smooth low-gDNA behavior needs care

Avoiding a hard zero mode is correct, but tiny nonzero priors can still create numerical weirdness. The implementation should use log-space likelihoods and explicit posterior intervals rather than if/else mode switches.

### R7. All-captured, all-expressed, unstranded loci remain weakly identifiable

v3 proposed disabling gDNA. v4 rejects that. The right behavior is to predict gDNA from the captured boundary/density model and expose high uncertainty. Accuracy will depend on how well the global captured component has been learned elsewhere.

### R8. Short-region boundary opportunity is weak

Very short regions have little boundary opportunity after FL clipping. Boundary evidence should then shrink strongly to the component prior. This is another reason `gamma_r` must be shrunk before it reaches `A_r`.

### R9. Mappability is still missing

A future mappability-corrected effective length would improve gDNA density estimation, especially for intergenic and intronic background. v4 should not block on this, but the density model should be written so mappability can replace raw opportunity later.

### R10. Knob control needs restraint

Expose as few user-facing parameters as possible. Prefer:

```text
gdna_confidence             existing confidence concept for upper bounds
background_trim_fraction    one robustness control
max_calibration_passes      mostly diagnostic/safety
```

Everything else should be learned, derived, or hidden behind documented profiles until benchmarks prove it needs exposure.

---

## 13. Non-Goals

- v4 does not try to make capture and non-capture RNA TPM comparable.
- v4 does not require a probe BED.
- v4 does not add an RNA exposure correction.
- v4 does not implement per-fragment exposure scoring.
- v4 does not solve unstranded capture perfectly; it gives the best available boundary-imputation route with explicit uncertainty.
- v4 does not require new C++ scanner features.

---

## 14. Short Version

v4 replaces v3's three pools with four posterior state classes:

```text
background
 gdna_only_capture
 expressed_capture
 expressed_offtarget
```

Strand-specific data can seed all four because it can deconvolve gDNA even in expressed regions. Unstranded data can seed background and gDNA-only capture, but expressed captured/offtarget regions remain posterior inferences.

The boundary model is no longer a threshold. It is a local gDNA imputation model that predicts contained gDNA from boundary crossing fragments, including internal exon regions via left-to-right and right-to-left sweeps. This is the mechanism that makes unstranded capture possible at all and improves strand-specific capture by adding local exposure evidence.
