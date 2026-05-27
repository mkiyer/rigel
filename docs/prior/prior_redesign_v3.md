# Prior Reformulation v3 - Mass-Conserving Unspliced Group RNA/gDNA Priors

Date: 2026-05-26
Status: implementation-ready plan after v2 review
Supersedes: `docs/prior/prior_redesign_v2.md`
Review input: `docs/prior/prior_redesign_v2_review.md`

## Executive Decision

Adopt the grouped RNA/gDNA prior, but revise v2 around one sharper contract:

```text
Only unspliced fragment mass feeds the calibration-derived EM prior.
For every region, prior_gdna_unspliced + prior_rna_unspliced = prior_unspliced_total.
Spliced fragments are excluded from the prior mass.
```

Spliced fragments already enter the EM as fragment evidence. They can inform calibration state probabilities and diagnostics, but they must not become prior pseudocounts. The prior is a mass-conserving deconvolution of unspliced regional evidence into gDNA and aggregate RNA.

The native EM prior remains grouped:

```text
one gDNA aggregate vs one aggregate containing all transcript-like RNA components
```

Within the RNA aggregate, EM remains fully responsible for partitioning mass across annotated mRNA and synthetic nRNA transcript components. There is no per-transcript RNA pseudocount and no split between mature and nascent RNA groups.

## Review Disposition

Accepted blocking critiques:

- B1: raw `m_g + m_r` as default budget is too hot. v3 uses a bounded prior budget.
- B2: v2 `mu_rna` formula can double-count and wrongly includes spliced mass. v3 replaces it with explicit unspliced mass deconvolution.
- B3: an external `enable_gdna` prior gate is unnecessary in the new paradigm. v3 removes it as a Python/calibration policy input and keeps only native structural candidate existence.
- B4: grouped SQUAREM needs explicit instrumentation and fallback criteria.

Accepted non-blocking critiques:

- N1: logit bias is odds-scale and should not be sold as a simple linear share shift.
- N3: if both posterior RNA mass and carried RNA mass are zero, `alpha_rna` is inactive for that update.
- N4: `assemble_priors()` needs an explicit prior config or `EMConfig` input and the pipeline call site must change.
- N5: golden/output migration must name the affected locus prior columns.
- N6: projection diagnostics should check allocated prior mass against observed unspliced mass.
- N7: VBEM warm start is count-space and must not be normalized before SQUAREM.
- N8: do not double-normalize by exposure; add diagnostics instead.

Rejected or reframed critique:

- N2 suggested splitting annotated RNA and synthetic nRNA if siphon worsens. v3 rejects that direction. All RNA remains one aggregate. In unstranded data, nRNA may intentionally enter the gDNA prior because strand information cannot separate it from gDNA. In strand-specific data, regional strand deconvolution can preserve some nRNA as RNA prior mass. If siphon is unacceptable, improve regional deconvolution, not within-EM RNA grouping.

## Core Model

For each `MultiLocus`, native EM has:

```text
components 0..T-1 = transcript-like RNA components
component T       = one gDNA component
```

The grouped prior acts only on:

```text
G = total gDNA component mass
R = sum of all transcript-like RNA component mass
```

It does not decide the distribution inside `R`.

After each E-step, define count-space component totals:

```text
n_t   = unambiguous_count_t + ambiguous_posterior_count_t
n_g   = unambiguous_count_g + ambiguous_posterior_count_g
N_rna = sum_t n_t
```

Calibration provides additive pseudocounts:

```text
a_g = alpha_gdna_add[locus]
a_r = alpha_rna_add[locus]
```

The grouped update target is:

```text
G = n_g + a_g_effective
R = N_rna + a_r_effective
S = G + R
```

When `N_rna > eps`, the within-RNA distribution is:

```text
q_t = n_t / N_rna
```

Then:

```text
theta_g = G / S
theta_t = (R / S) * q_t
```

Equivalently, the RNA pseudocount is dynamically distributed as:

```text
dynamic_rna_prior_t = a_r_effective * n_t / N_rna
```

Therefore, a transcript with no RNA evidence receives no RNA prior mass.

## Regional Prior-Mass Deconvolution

### New Calibration Object

Add a small mass-conserving prior deconvolution payload to the calibration path. This can be a standalone dataclass or fields on `RegionCalibration`.

Recommended fields:

```python
@dataclass(frozen=True, slots=True)
class PriorMassDeconvolution:
    unspliced_total: np.ndarray          # float32[R]
    gdna_unspliced_mean: np.ndarray      # float32[R]
    rna_unspliced_mean: np.ndarray       # float32[R]
    method: np.ndarray                   # uint8[R]
    precision: np.ndarray                # float32[R]
    flags: np.ndarray                    # uint16[R]
```

Required invariant:

```text
gdna_unspliced_mean[r] + rna_unspliced_mean[r] == unspliced_total[r]
```

within floating tolerance for every region.

`RegionCalibration` should carry this object or equivalent aligned arrays:

```python
prior_mass: PriorMassDeconvolution
```

Keep existing `mu_gdna`, `upper_gdna`, and `A_r` for calibration diagnostics and gDNA effective-length weighting. Do not use those fields directly as prior pseudocounts after v3 cutover.

### Prior-Owned Unspliced Total

The prior-owned total is the exact region-owned unspliced fragment mass from the fractional accumulator:

```text
U_r = contained_unspliced_total[r]
    + boundary_left_unspliced_total[r]
    + boundary_right_unspliced_total[r]
```

Equivalently, this is the unspliced compatible mass currently exposed as `DensityObservation.observed_compatible_count` when built from the ledger.

Spliced channels are excluded:

```text
contained_spliced_total
boundary_spliced_total
spliced_total
```

are not added to `U_r` and do not feed `alpha_rna_add`.

### Strand-Informative Deconvolution

When strand deconvolution is identifiable, use compartment strand posterior means.

For each region:

```text
D_r = contained_mean[r]
    + boundary_left_mean[r]
    + boundary_right_mean[r]

U_r = contained_total[r]
    + boundary_left_total[r]
    + boundary_right_total[r]

D_r = clip(D_r, 0, U_r)
R_r = U_r - D_r
```

where `contained_mean`, `boundary_left_mean`, and `boundary_right_mean` are gDNA means from `RegionGdnaChannelEstimate` and the totals are the matching unspliced compartment totals from `CompartmentStrandCounts` or `RegionCountLedger`.

This uses the precise deconvolution concept:

```text
400 unspliced fragments -> 100 gDNA + 300 RNA
```

The invariant is exact by construction:

```text
N_gdna + N_rna = N_unspliced_total
```

### Strand-Uninformative Deconvolution

When strand is unavailable or near-unstranded, there is no reliable way to separate nascent RNA from gDNA by strand balance. In this case the fallback is intentionally gDNA-permissive.

Use density/capture calibration only to choose the gDNA share, then force mass conservation:

```text
D_density = clip(region_calibration.mu_gdna, 0, U_r)
D_r       = D_density
R_r       = U_r - D_r
```

Because the unstranded calibration path has no strand RNA protection, intronic/nascent-like unspliced mass may enter `D_r`. This is intentional. For unstranded data, the gDNA component is allowed to absorb mass that is biologically nascent RNA but statistically indistinguishable from gDNA.

Do not add spliced evidence to `R_r` to counteract this. Spliced fragments already participate in EM likelihoods and should not become prior mass.

### Diagnostics For Deconvolution

Record per-region summaries:

```text
sum_unspliced_total
sum_gdna_unspliced_mean
sum_rna_unspliced_mean
max_abs_mass_conservation_error
method histogram
precision summary
```

Also record stratified summaries by strand-identifiable vs strand-uninformative regions.

## Projection To MultiLocus Priors

`assemble_priors()` should allocate the two deconvolved regional masses to `MultiLocus` rows using the existing geometry allocation machinery:

```text
gdna_expected_count[locus] = allocated sum of prior_mass.gdna_unspliced_mean
rna_expected_count[locus]  = allocated sum of prior_mass.rna_unspliced_mean
unspliced_total[locus]     = allocated sum of prior_mass.unspliced_total
```

Projection invariant:

```text
gdna_expected_count[locus] + rna_expected_count[locus]
    ~= unspliced_total[locus]
```

up to partial-coverage and unallocated-region diagnostics.

Add locus diagnostics:

```text
prior_unspliced_total
prior_mass_conservation_error
prior_allocated_fraction
unallocated_unspliced_count
```

The existing gDNA effective-length path remains:

```text
gdna_eff_len_unweighted = FL-PMF locus gDNA overlap effective length
gdna_em_exposure_weight = bp-weighted mean A_r over locus blocks
gdna_eff_len            = max(gdna_eff_len_unweighted * gdna_em_exposure_weight, 1.0)
```

This is not a second prior count. It is the gDNA component denominator in the likelihood. Do not divide or normalize `alpha_gdna_add` by `A_r` in Python.

Add a diagnostic only:

```text
gdna_prior_density = alpha_gdna_add / gdna_eff_len
```

For off-target loci, this should be roughly compatible with `rho_off` after accounting for prior strength and budget capping.

## Prior Budget And Operating Point

### Inputs

For each locus, after projection:

```text
m_g = gdna_expected_count[locus]
m_r = rna_expected_count[locus]
m   = m_g + m_r
```

The raw share is:

```text
w_raw = m_g / max(m, eps)
```

### Bounded Budget

Do not use `m_g + m_r` directly as the pseudocount budget. It is too hot for highly expressed loci.

Use a bounded budget that has two parts:

```text
balanced_budget = 2 * min(m_g, m_r)
edge_budget     = min(m, aggregate_prior_edge_count)
budget_raw      = max(balanced_budget, edge_budget)
budget_raw      = min(budget_raw, aggregate_prior_max_count)
budget          = min(aggregate_prior_strength * budget_raw,
                      aggregate_prior_max_count)
```

Default internal values after Phase 5 sentinel tuning:

```text
aggregate_prior_strength = 3.0
aggregate_prior_edge_count = 1000.0
aggregate_prior_max_count = 3000.0
```

Rationale:

- `2 * min(m_g, m_r)` follows the review's weaker-side constraint for mixed loci.
- `edge_budget` preserves useful one-sided priors when deconvolution says nearly all unspliced mass is RNA or nearly all is gDNA.
- `aggregate_prior_max_count` prevents calibration from injecting raw expression-scale pseudo-fragments into highly expressed loci.
- `aggregate_prior_strength` remains a smooth trust knob around this bounded scale.

Benchmark alternatives in Phase 5:

```text
budget_raw = 2 * min(m_g, m_r)
budget_raw = sqrt(m_g * m_r)
budget_raw = min(m, cap)
```

but do not implement raw `m_g + m_r` as the default.

### gDNA Bias

Bias operates on log odds, not directly on share:

```text
w = sigmoid(logit(clamp(w_raw, eps, 1 - eps)) + gdna_prior_logit_bias)
```

Default:

```text
gdna_prior_logit_bias = -6.0
```

Practical benchmark range:

```text
[-2.0, 2.0]
```

Interpretation:

- Positive values increase gDNA odds.
- Negative values decrease gDNA odds.
- Extreme `w_raw` values remain hard to move, which is desirable when deconvolution is confident.

Implementation note:

- The originally proposed starting point (`strength=1`, `edge_count=5`,
    `max_count=10`, `bias=0`) was too weak in Phase 5 zero-gDNA sentinels:
    structurally available gDNA siphoned RNA mass in pure mRNA, single-exon,
    wide-intron, and nRNA-double-counting baselines.
- The tuned defaults above intentionally provide strong one-sided RNA
    counter-prior mass and a conservative gDNA log-odds bias while preserving
    true-gDNA recovery in the measured sentinel subset.

The final prior masses are:

```text
alpha_gdna_add = budget * w
alpha_rna_add  = budget * (1 - w)
```

If `m <= eps`, set both to zero.

Do not expose `gdna_prior_logit_bias` as a polished public user knob until benchmark sweeps establish sensible defaults and ranges.

## Native Candidate Availability

v3 removes `enable_gdna` as a Python/calibration policy input.

Native code should derive a structural flag while extracting the locus sub-problem:

```text
has_gdna_candidate = any unit is unspliced and has finite gdna_log_lik
```

This is not a calibration gate. It is simply whether the gDNA component has any likelihood support in the sub-problem.

Native behavior:

```text
if !has_gdna_candidate:
    alpha_gdna_effective = 0
    alpha_rna_effective  = 0
    theta_gdna = 0
```

The aggregate prior only acts when there is an actual aggregate split to constrain. If no gDNA candidate exists, RNA is the only available group and no grouped RNA counter-prior is needed.

Migration note:

- During ABI transition, tests may keep an explicit `enable_gdna` array.
- Final production path should not require Python to pass `enable_gdna`.
- If the array remains for compatibility, native should validate it against the derived `has_gdna_candidate` in debug/profile builds and prefer the derived structural value.

## Native Grouped Update

### Effective Prior Inputs

Per locus:

```text
a_g_input = locus_alpha_gdna_add[locus]
a_r_input = locus_alpha_rna_add[locus]
```

Sanitize:

```text
a_g_input finite and >= 0
a_r_input finite and >= 0
```

Effective values for one update:

```text
if !has_gdna_candidate:
    a_g = 0
    a_r = 0
else:
    a_g = a_g_input
    a_r = a_r_input
```

Additional cold-start guard:

```text
if N_rna <= eps and carried_rna_mass <= eps:
    a_r = 0 for this update
```

Calibration cannot create expression on isoforms when neither the current E-step nor the carried state has RNA mass.

### Core Helper

Use one helper for MAP and VBEM count-space updates:

```cpp
apply_grouped_prior_update(raw_counts, carried_state, aggregate_prior, out)
```

Semantics:

```text
n_g = raw_counts[gdna] if has_gdna_candidate else 0
n_r = sum RNA raw_counts

if no gDNA candidate:
    a_g = 0
    a_r = 0
else:
    a_g = alpha_gdna_add
    a_r = alpha_rna_add

if n_r <= eps and carried_rna_mass <= eps:
    a_r = 0

G = n_g + a_g
R = n_r + a_r
```

RNA proportions:

```text
if n_r > eps:
    q_t = raw_counts[t] / n_r
elif carried_rna_mass > eps:
    q_t = carried_state[t] / carried_rna_mass
else:
    q_t is not used because a_r_effective = 0
```

Output in count space:

```text
out[gdna] = G
out[t]    = R * q_t if q_t exists else 0
```

For MAP, normalize `out` to theta. For VBEM, keep `out` unnormalized as alpha-like count mass.

### MAP Step

`grouped_map_step()`:

1. Compute E-step log weights from normalized theta and component `log_eff_len`.
2. Accumulate `em_totals`.
3. Build `raw_counts = unambig_totals + em_totals`.
4. Apply grouped prior update.
5. Normalize to `theta_new`.

### VBEM Step

`grouped_vbem_step()`:

1. Compute E-step log weights from `digamma(alpha) - digamma(sum alpha) - log_eff_len`.
2. Accumulate `em_totals`.
3. Build `raw_counts = unambig_totals + em_totals`.
4. Apply grouped prior update.
5. Return unnormalized `alpha_new`.

Remove the modeled `use_vbem ? 0.5 : 0.0` baseline. Keep only numerical floors for log and digamma stability.

### Warm Start

The warm start is count-space:

```text
warm_counts = unambiguous RNA totals + coverage-weighted ambiguous shares
```

Apply grouped prior update once to `warm_counts`.

MAP:

```text
theta_init = normalized grouped warm counts
```

VBEM:

```text
alpha_init = grouped warm counts
```

Do not normalize the VBEM warm start before SQUAREM. It is already in the count/alpha space expected by the VBEM fixed-point map.

## SQUAREM Plan

Grouped prior changes the fixed-point map. SQUAREM can still accelerate it, but it must be instrumented.

Add per-locus diagnostics when `emit_locus_stats` is enabled:

```text
squarem_step_scale_mean
squarem_step_scale_max
squarem_extrapolation_clamp_count
squarem_nonfinite_count
squarem_grouped_fallback_used
squarem_grouped_stabilization_fail_count
```

If a true acceptance/rejection line search is added, also record:

```text
squarem_extrapolation_attempts
squarem_extrapolation_accepts
squarem_acceptance_rate
```

Fallback rule:

```text
if grouped SQUAREM produces nonfinite state,
   or repeated stabilization failures,
   or convergence stalls past the configured iteration budget:
       rerun that locus with plain grouped EM
       set squarem_grouped_fallback_used = true
```

Phase 5 benchmarks must compare grouped-prior loci against unpriored loci for:

```text
iteration count
fallback rate
step scale distribution
clamp/nonfinite counts
```

If grouped SQUAREM acceptance collapses or fallback is common, use plain grouped EM for aggregate-prior loci until the acceleration path is repaired.

## Configuration

Initial internal config fields:

```python
aggregate_prior_strength: float = 3.0
aggregate_prior_edge_count: float = 1000.0
aggregate_prior_max_count: float = 3000.0
gdna_prior_logit_bias: float = -6.0
```

Validation:

```text
aggregate_prior_strength >= 0
aggregate_prior_edge_count >= 0
aggregate_prior_max_count >= 0
gdna_prior_logit_bias finite
aggregate_prior_edge_count <= aggregate_prior_max_count is recommended, not required
```

Thread these through `assemble_priors()` via either:

```python
assemble_priors(..., em_config: EMConfig)
```

or a smaller extracted prior config object. Update the call site in `pipeline.quant_from_buffer()`.

Do not expose CLI flags until Phase 5 establishes default behavior.

## PriorTable Contract

Extend `PriorTable` with:

```python
prior_unspliced_total: np.ndarray
rna_expected_count: np.ndarray
alpha_gdna_add: np.ndarray
alpha_rna_add: np.ndarray
prior_budget_raw: np.ndarray
prior_budget: np.ndarray
prior_gdna_share_raw: np.ndarray
prior_gdna_share_biased: np.ndarray
prior_mass_conservation_error: np.ndarray
gdna_prior_density: np.ndarray
```

Keep for one migration checkpoint:

```python
gdna_prior_count_em
```

as an alias to `alpha_gdna_add` in output diagnostics. Mark it deprecated in comments and tests.

`gdna_upper_count` remains diagnostic only unless a later design explicitly uses uncertainty.

## Output And Golden Migration

Affected output surfaces:

```text
loci.feather / loci.tsv locus prior columns
summary.json calibration.prior_table summaries
tests/golden/*_loci_df.tsv
tests/golden/*summary*.json if prior summaries are included
tests/test_golden_output.py expected schemas
```

Migration rule:

- One checkpoint may retain `gdna_prior_count` and `gdna_prior_count_em` aliases.
- Final schema should prefer `alpha_gdna_add`, `alpha_rna_add`, and `prior_unspliced_total`.
- Golden updates must be intentional and happen only after sentinel tests are understood.

## Phased Implementation Plan

### Phase 0 - Lock Math And Review Decisions

Scope: tests and documentation only.

Tasks:

1. Add pure Python tests for grouped update math with additive pseudocounts.
2. Test bounded budget formula:
   - large expressed RNA locus does not receive raw-count-scale prior mass,
   - one-sided RNA/gDNA cases still receive bounded prior mass,
   - `aggregate_prior_strength = 0` gives zero prior budget.
3. Test cold-start guard:
   - `N_rna == 0` and carried RNA mass is zero disables `a_r` for that update.
4. Test no-spliced-prior contract in a small synthetic region fixture.
5. Log Phase 0 start in `docs/newcal/new_cal_implementation_log.md`.

Exit criteria:

- The corrected v3 math is locked before native changes.

### Phase 1 - Add Mass-Conserving Unspliced Prior Deconvolution

Files:

```text
src/rigel/calibration/calibration_iteration.py
src/rigel/calibration/strand_deconv.py if additional mean fields/helpers are needed
src/rigel/calibration/_result.py
tests/test_calibration_iteration.py
tests/test_compartment_strand_deconv.py
```

Tasks:

1. Add `PriorMassDeconvolution` or equivalent fields to `RegionCalibration`.
2. Build `unspliced_total` from ledger unspliced compartments only.
3. For strand-informative regions, compute `gdna_unspliced_mean` from compartment gDNA posterior means and `rna_unspliced_mean = U - D`.
4. For strand-uninformative regions, compute `D = clip(mu_gdna, 0, U)` and `R = U - D`.
5. Enforce and test `D + R = U` per region.
6. Remove the v2 `mu_rna = spliced + residual` concept. Spliced mass is not prior mass.
7. Add summaries and method histograms.

Exit criteria:

- Calibration produces exact unspliced RNA/gDNA prior masses aligned to regions.

### Phase 2 - Project Paired Priors To Loci

Files:

```text
src/rigel/calibration/prior.py
src/rigel/config.py
src/rigel/pipeline.py
src/rigel/cli.py only if internal flags are exposed later
tests/test_calibration_prior.py
```

Tasks:

1. Add prior config fields to `EMConfig` or a nested prior config.
2. Update `assemble_priors()` signature and pipeline call site.
3. Allocate `prior_mass.gdna_unspliced_mean`, `prior_mass.rna_unspliced_mean`, and `prior_mass.unspliced_total` to `MultiLocus` rows.
4. Compute bounded budget and log-odds-biased share.
5. Fill `alpha_gdna_add` and `alpha_rna_add`.
6. Add mass-conservation and allocation diagnostics.
7. Keep `gdna_prior_count_em` as a deprecated alias for one checkpoint.

Exit criteria:

- Python produces paired bounded priors from unspliced mass only.

### Phase 3 - Native Grouped Solver And Candidate-Derived gDNA Availability

Files:

```text
src/rigel/native/em_solver.cpp
src/rigel/estimator.py
src/rigel/pipeline.py
tests/test_batch_em_impl.py
tests/test_native_gdna_eligibility.py or renamed grouped-prior native tests
```

Tasks:

1. Update native ABI to accept `locus_alpha_gdna_add` and `locus_alpha_rna_add`.
2. Remove Python-calibration `enable_gdna` from the final production ABI, or demote it to temporary compatibility only.
3. Derive `has_gdna_candidate` inside native locus extraction.
4. If `!has_gdna_candidate`, zero both effective priors and gDNA output.
5. Add cold-start guard for `a_r`.
6. Implement grouped MAP and VBEM steps.
7. Remove modeled per-component VBEM baseline.
8. Make VBEM warm start count-space and unnormalized.
9. Add SQUAREM diagnostics and fallback.
10. Rebuild native extension before testing.

Exit criteria:

- Native EM applies grouped priors inside every update step.
- No fixed RNA pseudocount reaches individual transcripts.
- gDNA candidate availability is structural, not calibration-gated.

### Phase 4 - Output Migration

Files:

```text
src/rigel/pipeline.py
src/rigel/estimator.py
tests/golden/*_loci_df.tsv
tests/golden/*summary*.json
tests/test_golden_output.py
```

Tasks:

1. Add output columns for paired prior diagnostics.
2. Keep old prior column aliases for one checkpoint if needed.
3. Update schema tests and goldens intentionally.
4. Document output column meaning in user docs after behavior stabilizes.

Exit criteria:

- Users can audit prior mass, budget, share, and density for every locus.

### Phase 5 - Validation And Benchmark Tuning

Required sentinel tests:

```text
pure mRNA zero-gDNA
single exon zero-gDNA
wide intron zero-gDNA
gDNA light/heavy
strand-specific intronic nRNA preservation
unstranded intronic nRNA siphon into gDNA
multi-isoform sparsity preservation
```

Required benchmark sweeps:

```text
aggregate_prior_strength in {0.0, 0.3, 1.0, 3.0}
aggregate_prior_edge_count in {1.0, 2.0, 5.0, 10.0}
aggregate_prior_max_count in {2.0, 5.0, 10.0, 20.0}
gdna_prior_logit_bias in {-2.0, -1.0, 0.0, 1.0, 2.0}
```

Report:

```text
false-positive gDNA in zero-gDNA data
false-negative gDNA in true-gDNA data
strand-specific nRNA retained as RNA
unstranded nRNA absorbed by gDNA
unsupported isoform mass
SQUAREM diagnostics and fallback rate
prior mass conservation errors
prior density vs rho_off diagnostics
```

The nRNA siphon benchmark is required, but it is not a trigger to split RNA into multiple EM prior groups. If siphon behavior is wrong, revise regional unspliced deconvolution or operating-point defaults.

Exit criteria:

- Defaults are chosen from measured false-positive/false-negative tradeoffs.
- Grouped prior improves zero-gDNA sentinels without destroying true-gDNA recovery.
- Spliced fragments never contribute to prior pseudocount mass.

## Acceptance Criteria

The redesign is acceptable only when all are true:

1. Region prior masses obey `gDNA + RNA = unspliced_total` per region.
2. Spliced fragments do not feed `alpha_gdna_add` or `alpha_rna_add`.
3. All transcript-like RNA remains one aggregate prior group.
4. Native gDNA availability is derived from actual gDNA likelihood candidates, not calibration posterior thresholds.
5. If no gDNA candidate exists, both grouped prior sides are inactive.
6. If no RNA evidence or carried RNA state exists, `alpha_rna` is inactive for that update.
7. Bounded budget prevents raw-count-scale prior hammers in highly expressed loci.
8. `aggregate_prior_strength = 0` behaves like unpriored EM apart from numerical floors.
9. SQUAREM grouped-prior behavior is instrumented, and fallback is rare or deliberately used.
10. Output diagnostics expose all prior mass and budget decisions.

## Implementation Checklist

Phase 0:

- [ ] Add grouped update math tests.
- [ ] Add bounded budget tests.
- [ ] Add unspliced-only prior tests.
- [ ] Add cold-start guard tests.

Phase 1:

- [ ] Add `PriorMassDeconvolution` or equivalent fields.
- [ ] Build unspliced-only `U_r` from ledger compartments.
- [ ] Implement strand-informative mass-conserving deconvolution.
- [ ] Implement strand-uninformative clipped-density fallback.
- [ ] Add region diagnostics.

Phase 2:

- [ ] Add prior config fields.
- [ ] Update `assemble_priors()` signature and pipeline call.
- [ ] Allocate gDNA/RNA/unspliced prior masses to loci.
- [ ] Compute bounded budget and alpha fields.
- [ ] Add PriorTable diagnostics.

Phase 3:

- [ ] Update native ABI for paired alpha fields.
- [ ] Derive native `has_gdna_candidate`.
- [ ] Implement grouped MAP/VBEM update.
- [ ] Remove modeled VBEM baseline.
- [ ] Add cold-start and no-gDNA-candidate guards.
- [ ] Add SQUAREM diagnostics/fallback.
- [ ] Rebuild native extension.

Phase 4:

- [ ] Update locus output diagnostics.
- [ ] Update schema/golden tests intentionally.

Phase 5:

- [ ] Run sentinel tests.
- [ ] Run budget/bias sweeps.
- [ ] Choose defaults from benchmark tradeoffs.
