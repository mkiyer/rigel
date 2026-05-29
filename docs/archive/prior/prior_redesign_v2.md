# Prior Reformulation v2 - Calibration-Derived Group RNA/gDNA Priors

Date: 2026-05-26
Status: implementation plan
Inputs reviewed:

- `docs/prior/prior_redesign_v1.md`
- `docs/prior/prior_redesign_gemini.md`
- `docs/newcal/gdna_prior_current_construction.md`

## Executive Decision

Adopt a grouped RNA/gDNA prior.

The prior should operate on the aggregate split:

```text
total gDNA mass vs total transcript-like RNA mass
```

It must not place positive calibration-derived pseudocounts on individual RNA transcripts. The EM should still decide how RNA mass is partitioned across annotated mRNA and synthetic nRNA transcript components.

The production target is:

```text
RegionCalibration
-> regional mu_gdna and mu_rna
-> per-MultiLocus alpha_gdna and alpha_rna
-> native grouped RNA/gDNA M-step
-> no per-transcript RNA prior floor
```

This plan uses additive pseudocount notation throughout:

```text
alpha_gdna >= 0 means extra gDNA pseudo-fragment mass
alpha_rna  >= 0 means extra aggregate RNA pseudo-fragment mass
```

This intentionally differs from standard Dirichlet concentration notation where MAP updates contain `alpha - 1`. Calibration emits expected counts, so the implementation should consume them as non-negative additive evidence. If someone wants to map this back to a standard Beta/Dirichlet density, these additive counts correspond to concentration parameters `alpha_standard = alpha_additive + 1` under MAP.

## Non-Negotiable Constraints

### Transcript-First Biology

Rigel internal inference uses transcripts and transcript-derived states as model units. Gene-level summaries are output convenience only.

Therefore:

- Phase 1 groups all transcript-like components in a locus into one aggregate RNA group.
- Do not introduce gene-level priors into internal EM.
- Do not group by `gene_id` as a production inference primitive.
- A future multi-group extension must use transcript-derived or region-derived groups, not gene abstractions, unless the project-level modeling rule changes.

This is the main correction to `prior_redesign_v1.md` Option C. The group-prior idea is right; gene grouping is not acceptable as an internal default for this codebase.

### No Per-Transcript RNA Pseudocounts

Calibration may say the locus is RNA-rich, but it does not identify which isoform should carry that RNA mass. Any fixed positive count on each transcript creates a false expression floor in large isoform sets.

The native update must preserve this invariant:

```text
If a transcript has zero RNA posterior count and zero carried RNA proportion,
the grouped RNA prior does not create expression for it.
```

### gDNA Eligibility Remains Technical

`enable_gdna` remains only:

```text
the locus has at least one unspliced unit with finite gDNA likelihood
```

It must not depend on `mu_gdna`, `upper_gdna`, `alpha_gdna`, `alpha_rna`, `A_r`, or posterior thresholds. The new prior controls mass smoothly; it does not reintroduce the old `gdna_upper > eps` gate.

### Separate Strength From Operating Point

Two concepts should stay distinct internally:

```text
aggregate_prior_strength: how much calibration constrains the EM
gdna_prior_logit_bias:    operating point, i.e. gDNA sensitivity/aggressiveness
```

The first is trust in calibration. The second moves the RNA/gDNA decision boundary. They can later be wrapped in one user-facing sensitivity control, but the implementation should keep them separate for benchmarking.

## What v1 And Gemini Got Right

Both parallel plans converge on the essential architecture:

- Treat RNA as a group.
- Put calibration pressure on the RNA-vs-gDNA boundary.
- Preserve within-RNA sparsity.
- Remove the current asymmetry where gDNA gets explicit calibration prior mass but RNA does not.
- Remove or neutralize the native VBEM `0.5` per-component baseline as a hidden modeling prior.

The Gemini derivation is useful theoretically because it writes the density:

```text
p(theta) proportional to theta_gdna^(a_g) * (sum RNA theta_t)^(a_r)
```

where this document uses `a_g` and `a_r` as additive pseudocount exponents, not standard concentrations. The v1 plan is useful operationally because it recognizes the required native M-step change and the need for a smooth user-tunable prior strength.

## Corrections To The Attached Plans

### Correction 1: Use Additive Pseudocounts, Not Raw Standard Dirichlet Alpha

The Gemini document uses the standard MAP expression with `alpha - 1`. That is mathematically correct if `alpha` means a standard Beta concentration. It is not the safest implementation parameterization for Rigel because calibration emits expected counts, and expected counts can be zero.

Rigel should store and pass:

```text
alpha_gdna_add >= 0
alpha_rna_add  >= 0
```

The grouped update then uses:

```text
G = n_gdna + alpha_gdna_add
R = n_rna  + alpha_rna_add
```

No `-1` appears in code. This avoids negative numerator cases, avoids special handling for standard concentrations below one, and makes `0.0` mean exactly no prior mass.

### Correction 2: Do Not Claim SQUAREM Is Automatically Unaffected

The grouped projection is nonlinear in the full component vector. It can still be used inside the fixed-point operator that SQUAREM accelerates, but we should not assume it commutes with SQUAREM extrapolation.

Implementation rule:

- Build grouped MAP/VBEM step functions as the actual fixed-point map.
- Let SQUAREM accelerate that map.
- Clamp invalid extrapolations exactly as the current solver already does.
- Add convergence and fallback tests.
- Benchmark mega-locus iteration counts before claiming performance neutrality.

### Correction 3: Multi-Group RNA Is A Later Extension

Phase 1 is one RNA group per `MultiLocus`. This solves the immediate zero-gDNA leakage and gDNA/RNA asymmetry without introducing another biological abstraction.

Later multi-group priors are allowed only if groups are transcript-derived or region-derived. Examples that may be acceptable later:

- annotated transcript-like RNA vs synthetic nRNA transcript-like RNA,
- disjoint transcript-span groups within a mega-locus,
- groups derived from connected-component structure after excluding cross-region multimapper bridges.

Gene-level priors are intentionally not part of this plan.

## Mathematical Contract

For a single `MultiLocus`, native EM has:

```text
components 0..T-1 = transcript-like RNA components
component T       = one gDNA component
```

After an E-step, define count-space totals:

```text
n_t    = unambiguous_count_t + ambiguous_posterior_count_t
n_g    = unambiguous_count_g + ambiguous_posterior_count_g
N_rna  = sum_t n_t
```

In practice `unambiguous_count_g` is normally zero because gDNA is a locus component appended for ambiguous/unspliced units, but the formula should not depend on that.

Calibration provides non-negative aggregate prior masses:

```text
a_g = alpha_gdna_add[locus]
a_r = alpha_rna_add[locus]
```

The grouped update target is:

```text
G = n_g + a_g
R = N_rna + a_r
S = G + R
```

The aggregate gDNA share is:

```text
p_g = G / S
p_r = R / S
```

The within-RNA distribution is data-driven:

```text
q_t = n_t / N_rna    if N_rna > eps
```

The updated component masses are:

```text
theta_g = p_g
theta_t = p_r * q_t
```

In unnormalized VBEM/count space, the same operator is:

```text
alpha_g_out = G
alpha_t_out = R * q_t
```

This is equivalent to dynamically distributing the aggregate RNA pseudocount in proportion to the current RNA evidence:

```text
dynamic_rna_prior_t = a_r * n_t / N_rna
```

So transcripts with no RNA evidence receive no RNA prior mass.

### Zero RNA Count Edge Case

When `N_rna <= eps`, the grouped RNA prior identifies aggregate RNA mass but not the isoform distribution. Do not uniformly distribute across all transcripts by default; that creates exactly the false isoform floor this redesign is meant to avoid.

Use a carried RNA distribution instead:

```text
q_t = previous RNA proportions if sum(previous RNA mass) > eps
q_t = warm-start RNA proportions if available
q_t = effective-length-weighted fallback only as a last resort, with a diagnostic flag
```

This makes the edge case explicit and testable. In normal loci, `N_rna > 0` because the locus exists through transcript-like candidates.

### Disabled gDNA Edge Case

When `enable_gdna == 0`, the gDNA component is technically absent. The grouped prior must not manufacture a gDNA component.

Native behavior should be:

```text
if !enable_gdna:
    a_g_effective = 0
    theta_g = 0
    all mass remains in RNA components
```

The Python prior table should still report the calibration-derived `alpha_gdna_add` for diagnostics, but native EM must ignore it when the component is technically unavailable.

## Calibration-Side Contract

### Add `mu_rna` To RegionCalibration

`RegionCalibration` currently carries `mu_gdna`, `upper_gdna`, `rna_lower`, `A_r`, and state probabilities. The new prior needs an expected RNA count aligned to the same region table:

```python
@dataclass(frozen=True, slots=True)
class RegionCalibration:
    p_states: np.ndarray
    mu_gdna: np.ndarray
    mu_rna: np.ndarray          # new: expected region-level transcript-like RNA count
    upper_gdna: np.ndarray
    rna_lower: np.ndarray
    A_r: np.ndarray
    ...
```

Initial definition for Phase 1:

```text
contained_total = observation.contained_count
spliced_total   = observation.spliced_count
p_expr          = region_calibration.p_expressed

rna_residual = max(contained_total - mu_gdna, 0)
mu_rna       = spliced_total + max(rna_lower, p_expr * rna_residual)
mu_rna       = min(mu_rna, spliced_total + contained_total)
```

Rationale:

- `mu_gdna` is the final calibration estimate after boundary sweep and state mixing.
- `contained_total - mu_gdna` is the natural contained RNA residual.
- `p_expr` prevents weak background residuals from becoming strong RNA priors in non-expressed regions.
- `rna_lower` preserves strand-derived RNA protection when available.
- `spliced_total` is unambiguous expression evidence at the region level and should support the aggregate RNA group.
- The upper bound prevents the RNA expectation from exceeding observed region-compatible count mass.

This formula is intentionally conservative and should be validated in Phase 1 tests. If benchmarks show it is too conservative or too aggressive, this is the correct layer to revise, not the native grouped prior math.

### Region Summaries And Diagnostics

Update calibration summaries to include:

```text
mu_rna summary stats
mu_total = mu_rna + mu_gdna summary stats
rna/gdna share summary stats
```

Also include a consistency diagnostic:

```text
sum(mu_rna + mu_gdna) vs sum(observation.contained_count + observation.spliced_count)
```

## Prior Projection Contract

`assemble_priors()` should project both regional means onto each `MultiLocus` using the same geometry allocation contract already used for `mu_gdna`.

For each locus:

```text
gdna_expected_count[locus] = allocated sum of region_calibration.mu_gdna
rna_expected_count[locus]  = allocated sum of region_calibration.mu_rna
```

Keep the current gDNA denominator path:

```text
gdna_eff_len_unweighted = FL-PMF locus gDNA overlap effective length
gdna_em_exposure_weight = bp-weighted mean A_r over locus blocks
gdna_eff_len            = max(gdna_eff_len_unweighted * gdna_em_exposure_weight, 1.0)
```

Do not create per-transcript RNA prior denominators. RNA transcript effective lengths already enter native EM through `log_eff_len[t]` in the E-step. The aggregate RNA prior is count-mass on `sum_t theta_t`, not a per-transcript exposure correction.

### Convert Expected Counts To Aggregate Prior Mass

Let:

```text
m_g = gdna_expected_count[locus]
m_r = rna_expected_count[locus]
m   = m_g + m_r
```

Define a count-scale prior budget:

```text
budget_raw = m
```

Phase 1 should implement this direct count-scale budget because it matches the current gDNA prior scale: the existing system already passes allocated `mu_gdna` as a physical count prior. The new system simply gives RNA the missing paired aggregate counter-prior.

Then apply global strength and optional caps:

```text
budget = aggregate_prior_strength * budget_raw
budget = min(budget, aggregate_prior_max_count) if cap is configured
```

Default proposal:

```text
aggregate_prior_strength = 1.0
aggregate_prior_max_count = None
```

This default is deliberately conservative with respect to the current implementation: gDNA keeps the same raw count scale, while RNA receives the missing paired count scale.

### gDNA Aggressiveness / Sensitivity

Apply operating-point bias on the share, not on the total budget:

```text
w_raw = m_g / max(m_g + m_r, eps)
w     = sigmoid(logit(clamp(w_raw, eps, 1 - eps)) + gdna_prior_logit_bias)

alpha_gdna_add = budget * w
alpha_rna_add  = budget * (1 - w)
```

Defaults:

```text
gdna_prior_logit_bias = 0.0
```

Interpretation:

- Negative bias protects RNA and reduces false-positive gDNA calls.
- Positive bias increases gDNA sensitivity.
- Bias does not change total prior strength, only the RNA/gDNA operating point.

This is preferable to multiplying `alpha_gdna` alone because it avoids unintentionally changing calibration trust at the same time as sensitivity.

### PriorTable Fields

Extend `PriorTable` with:

```python
alpha_gdna_add: np.ndarray
alpha_rna_add: np.ndarray
rna_expected_count: np.ndarray
prior_budget: np.ndarray
prior_gdna_share_raw: np.ndarray
prior_gdna_share_biased: np.ndarray
```

Keep current fields during migration:

```python
gdna_prior_count_em
gdna_expected_count
gdna_upper_count
gdna_eff_len
enable_gdna
```

During Phase 1/2 migration, `gdna_prior_count_em` can be retained as an alias/diagnostic equal to `alpha_gdna_add` or `gdna_expected_count`, but final native code should consume `alpha_gdna_add` and `alpha_rna_add` explicitly.

## Native EM Contract

### New Native Inputs

Replace or extend native `batch_locus_em_partitioned()` inputs:

```cpp
f64_1d locus_alpha_gdna_add
f64_1d locus_alpha_rna_add
u8_1d  locus_enable_gdna
f64_1d locus_gdna_eff_lens
```

The old `locus_gdna_prior_count` should be removed after compatibility tests are migrated.

### Group Prior Struct

Add a small native struct:

```cpp
struct AggregatePrior {
    double alpha_gdna = 0.0;
    double alpha_rna = 0.0;
    int gdna_idx = -1;
    bool enable_gdna = false;
};
```

Sanitize inputs once per locus:

```text
alpha_gdna = finite and >= 0 else error
alpha_rna  = finite and >= 0 else error
if !enable_gdna: alpha_gdna_effective = 0
```

### Replace Per-Component Prior Vector

Remove the current mode-aware per-component baseline:

```cpp
baseline = use_vbem ? 0.5 : 0.0
```

The grouped prior itself supplies all modeled prior mass. Numerical floors remain numerical only.

RNA components receive no fixed per-component pseudocount.

### Grouped MAP Step

Create a grouped MAP step that performs the current E-step, then applies the aggregate update:

```cpp
static void grouped_map_step(
    const double* theta,
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    const AggregatePrior& prior,
    double* em_totals,
    double* theta_new,
    int n_components,
    ...)
{
    // E-step identical except log weights use theta and log_eff_len.
    // raw_counts[c] = unambig_totals[c] + em_totals[c]
    apply_grouped_prior_update(
        raw_counts,
        theta,          // carried RNA proportions for zero-RNA edge case
        prior,
        theta_new,
        normalized=true);
}
```

`apply_grouped_prior_update()` is the core helper and should be unit-tested through a small native-exposed test hook if practical.

### Grouped VBEM Step

Create a grouped VBEM step with the same grouped update but unnormalized output:

```cpp
static void grouped_vbem_step(
    const double* alpha,
    ...,
    double* alpha_new)
{
    // E-step uses digamma(alpha) - digamma(sum alpha) - log_eff_len.
    // raw_counts[c] = unambig_totals[c] + em_totals[c]
    apply_grouped_prior_update(
        raw_counts,
        alpha,          // carried RNA proportions for zero-RNA edge case
        prior,
        alpha_new,
        normalized=false);
}
```

This removes the hidden Jeffreys asymmetry while preserving the existing MAP/VBEM distinction in the E-step.

### Core Helper Semantics

Pseudocode:

```cpp
static void apply_grouped_prior_update(
    const double* raw_counts,
    const double* carried_state,
    const AggregatePrior& prior,
    double* out,
    int n_components,
    bool normalized)
{
    int g = prior.gdna_idx;
    double n_g = (prior.enable_gdna && g >= 0) ? raw_counts[g] : 0.0;
    double n_r = 0.0;
    for each RNA component c:
        n_r += raw_counts[c];

    double a_g = prior.enable_gdna ? prior.alpha_gdna : 0.0;
    double a_r = prior.alpha_rna;
    double G = n_g + a_g;
    double R = n_r + a_r;

    // q_t from evidence if possible, otherwise carried RNA proportions.
    if (n_r > eps):
        q_t = raw_counts[t] / n_r;
    else:
        q_t = carried_state[t] / sum_carried_rna;

    out[g] = G;
    out[t] = R * q_t;

    if normalized:
        divide out by sum(out)
}
```

Important guards:

- Disabled gDNA output is exactly zero.
- Ineligible components stay zero.
- If `G + R == 0`, fall back to normalized eligible RNA carried state, or uniform over eligible RNA only as a last-resort diagnostic path.
- Never distribute `alpha_rna` uniformly across transcripts when there is RNA evidence.

### Warm Start

Warm start should be consistent with the grouped prior but not become a hidden gate.

Phase 1 native behavior:

1. Build the same coverage-weighted warm start from eligible components as today.
2. Apply `apply_grouped_prior_update()` once to the warm-start count vector using the warm-start vector itself as carried state.
3. Use that as `theta_init` for MAP or `alpha_init` seed for VBEM.

This gives calibration a smooth initial influence while keeping gDNA technically recoverable through likelihood if enabled.

Do not multiply gDNA rows by zero when `alpha_gdna == 0`; that is a gate by another name. If a later benchmark demands warm-start tempering, add it as a separate documented parameter with a nonzero floor.

### SQUAREM Integration

The current `run_squarem()` should call grouped step functions in place of `map_em_step()` and `vbem_step()` when aggregate priors are enabled.

For MAP:

- SQUAREM state remains normalized `theta`.
- After extrapolation, clamp negative entries and renormalize as today.
- Stabilization step calls `grouped_map_step()`.

For VBEM:

- SQUAREM state remains unnormalized count/alpha-like mass.
- After extrapolation, clamp to numerical epsilon only for eligible components.
- Stabilization step calls `grouped_vbem_step()`.
- Final `theta` is normalized `alpha_out`.

Add a safety fallback:

```text
if grouped SQUAREM produces nonfinite state or stalls for max iterations:
    rerun the locus with plain grouped EM for that locus
    increment a diagnostic counter
```

This fallback should be rare, but mega-loci deserve a guardrail.

## Phased Implementation Plan

### Phase 0 - Lock The Contract And Tests Before Native Changes

Scope: no behavior change.

Tasks:

1. Add pure Python tests for grouped update math using a small helper in tests.
2. Test additive pseudocount semantics:
   - `alpha_gdna = alpha_rna = 0` gives unpriored aggregate update.
   - `alpha_rna > 0, alpha_gdna = 0` shrinks gDNA share smoothly.
   - `alpha_gdna > 0, alpha_rna = 0` increases gDNA share smoothly.
   - zero-count transcripts do not receive RNA prior mass when `N_rna > 0`.
3. Add a checkpoint to `docs/newcal/new_cal_implementation_log.md` when Phase 0 starts.

Exit criteria:

- The expected grouped update is specified by tests before touching C++.

### Phase 1 - Surface Regional `mu_rna`

Files:

- `src/rigel/calibration/calibration_iteration.py`
- `src/rigel/calibration/_result.py`
- tests for calibration iteration/result summaries

Tasks:

1. Add `mu_rna` to `CalibrationStepResult` and `RegionCalibration`.
2. Compute `mu_rna` in `calibration_e_step()` using the Phase 1 formula above.
3. Validate shape, finiteness, and non-negativity in `RegionCalibration.__post_init__()`.
4. Add summary stats for `mu_rna` and `mu_rna + mu_gdna`.
5. Add tests where:
   - spliced evidence increases `mu_rna`,
   - pure gDNA evidence gives low `mu_rna`,
   - strand RNA lower bound protects `mu_rna`,
   - `mu_rna` does not exceed compatible observed count.

Exit criteria:

- Calibration result carries aligned regional RNA and gDNA expected counts.

### Phase 2 - Build Paired PriorTable Fields

Files:

- `src/rigel/calibration/prior.py`
- `src/rigel/config.py`
- `src/rigel/cli.py` if CLI exposure is desired immediately
- prior assembly tests

Tasks:

1. Extend `EMConfig` with internal defaults:

```python
aggregate_prior_strength: float = 1.0
aggregate_prior_max_count: float | None = None
gdna_prior_logit_bias: float = 0.0
```

2. Thread `em_config` into `assemble_priors()` or pass a small prior config object.
3. Allocate `mu_rna` to `MultiLocus` with the same geometry allocator used for `mu_gdna`.
4. Compute `alpha_gdna_add` and `alpha_rna_add` from expected counts, strength, optional cap, and logit bias.
5. Preserve current `enable_gdna` technical eligibility unchanged.
6. Add PriorTable summaries:

```text
sum_alpha_gdna_add
sum_alpha_rna_add
prior_budget stats
prior_gdna_share_raw stats
prior_gdna_share_biased stats
rna_expected_count stats
```

7. Keep old `gdna_prior_count_em` diagnostics for one migration checkpoint.

Exit criteria:

- Python can produce paired locus-level prior masses without changing native EM behavior yet.

### Phase 3 - Add Native Two-Group Prior Solver

Files:

- `src/rigel/native/em_solver.cpp`
- `src/rigel/estimator.py`
- `src/rigel/pipeline.py`
- native/batch EM tests

Tasks:

1. Update native ABI to accept `locus_alpha_gdna_add` and `locus_alpha_rna_add`.
2. Replace `compute_gdna_prior_and_warm_start()` with a grouped-prior initializer.
3. Add `AggregatePrior` and `apply_grouped_prior_update()`.
4. Add `grouped_map_step()` and `grouped_vbem_step()`.
5. Remove the modeled VBEM `0.5` baseline from aggregate-prior mode.
6. Keep numerical floors only for log/digamma stability.
7. Integrate grouped steps into SQUAREM with fallback diagnostics.
8. Rebuild native extension before running tests.

Exit criteria:

- Native EM consumes paired aggregate priors.
- RNA components receive no fixed per-transcript calibration pseudocount.
- Disabled gDNA remains exactly disabled.

### Phase 4 - Pipeline Output And Golden Migration

Files:

- `src/rigel/pipeline.py`
- estimator locus result metadata
- golden tests

Tasks:

1. Add locus diagnostics:

```text
alpha_gdna_add
alpha_rna_add
rna_expected_count
prior_budget
prior_gdna_share_raw
prior_gdna_share_biased
```

2. Deprecate or rename old columns:

```text
gdna_prior_count
gdna_prior_count_em
```

Recommended migration:

- Keep old names as aliases in one checkpoint if tests and downstream scripts need it.
- Final schema should prefer `alpha_gdna_add` / `alpha_rna_add`.

3. Update golden outputs intentionally after behavior is validated.

Exit criteria:

- Output diagnostics explain why each locus was pulled toward RNA or gDNA.

### Phase 5 - Focused Validation And Operating-Point Tuning

Required tests:

1. Prior projection tests:
   - geometry allocation preserves sums for both RNA and gDNA,
   - logit bias changes share but not total budget,
   - strength changes total budget but not raw share,
   - technical `enable_gdna` does not depend on prior share.

2. Native grouped EM tests:
   - one RNA component plus one gDNA component matches closed-form aggregate update,
   - many RNA components with one supported isoform preserves zero/near-zero unsupported isoforms,
   - zero-gDNA calibration with strong RNA prior suppresses gDNA without disabling it,
   - true-gDNA calibration with strong gDNA prior recovers gDNA,
   - `enable_gdna == 0` produces zero gDNA regardless of `alpha_gdna_add`.

3. Scenario sentinels:
   - pure mRNA baseline,
   - single exon clean baseline,
   - wide intron baseline,
   - gDNA light/heavy scenarios,
   - nRNA-heavy scenarios.

4. Benchmark sweeps:

```text
aggregate_prior_strength in {0.0, 0.1, 0.3, 1.0, 3.0}
gdna_prior_logit_bias in {-2.0, -1.0, 0.0, 1.0, 2.0}
```

Report:

- false-positive gDNA in zero-gDNA libraries,
- false-negative gDNA in true-gDNA libraries,
- mRNA/nRNA siphon behavior,
- unsupported isoform mass in multi-isoform loci,
- SQUAREM iterations and fallback counts.

Exit criteria:

- Pick defaults from measured tradeoffs, not from green tests alone.

### Phase 6 - Optional Multi-Group RNA Extension

Only after the two-group RNA/gDNA prior works.

Allowed group sources must be transcript-derived or region-derived, for example:

- one aggregate for annotated transcript-like RNA and one aggregate for synthetic nRNA,
- groups by transcript-derived connected subcomponents inside a mega-locus,
- groups by non-overlapping genomic span classes.

Not allowed as default internal inference:

- one group per gene,
- gene-name or gene-id priors,
- gene-level calibration summaries feeding EM.

Native generalization:

```text
(theta_gdna, sum group_0, sum group_1, ..., sum group_K)
```

with no within-group per-component pseudocounts.

This is a later design. Do not block Phase 1 on it.

## Configuration Proposal

Add to `EMConfig` first, CLI later after defaults are benchmarked:

```python
aggregate_prior_strength: float = 1.0
aggregate_prior_max_count: float | None = None
gdna_prior_logit_bias: float = 0.0
```

Validation:

```text
aggregate_prior_strength >= 0
aggregate_prior_max_count is None or >= 0
gdna_prior_logit_bias finite
```

Possible future CLI:

```text
--aggregate-prior-strength
--gdna-prior-bias
```

Possible user-friendly sensitivity wrapper later:

```text
gDNA sensitivity low    -> negative gdna_prior_logit_bias
gDNA sensitivity normal -> 0
gDNA sensitivity high   -> positive gdna_prior_logit_bias
```

Do not expose a polished public knob until Phase 5 benchmark sweeps choose defaults and ranges.

## Compatibility And Migration Notes

This is a native ABI change. After editing `src/rigel/native/em_solver.cpp`, rebuild with:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

Recommended migration sequence:

1. Add `mu_rna` and paired `PriorTable` diagnostics while old native input still exists.
2. Add native grouped inputs and switch the live pipeline in one clean checkpoint.
3. Refresh goldens intentionally.
4. Remove old `gdna_prior_count_em` native semantics after tests and downstream scripts are updated.

No production feature flag is planned. During development, tests may instantiate old/new helpers directly, but the live pipeline should cut over cleanly once the grouped solver is ready.

## Acceptance Criteria

The redesign is successful only if all of the following are true:

1. Zero-gDNA sentinel scenarios no longer leak large gDNA mass without reintroducing a hard eligibility gate.
2. True-gDNA scenarios still recover gDNA when likelihood and calibration support it.
3. Multi-isoform loci do not gain a broad false-positive transcript floor.
4. Native `enable_gdna` remains technical eligibility only.
5. `aggregate_prior_strength = 0` behaves like pure likelihood EM apart from numerical floors.
6. `gdna_prior_logit_bias` changes RNA/gDNA operating point without changing total prior budget.
7. SQUAREM remains stable on mega-loci, or fallback diagnostics prove rare and bounded fallback behavior.
8. Output diagnostics make every prior contribution auditable at region, locus, and native EM levels.

## Implementation Checklist

Phase 0:

- [ ] Add grouped-prior math tests before C++ changes.
- [ ] Lock additive pseudocount semantics in tests.

Phase 1:

- [ ] Add `mu_rna` to `CalibrationStepResult`.
- [ ] Add `mu_rna` to `RegionCalibration`.
- [ ] Add calibration summary diagnostics.

Phase 2:

- [ ] Add `aggregate_prior_strength`, `aggregate_prior_max_count`, and `gdna_prior_logit_bias` config fields.
- [ ] Allocate regional `mu_rna` to loci.
- [ ] Add `alpha_gdna_add` and `alpha_rna_add` to `PriorTable`.

Phase 3:

- [ ] Update Python/native ABI.
- [ ] Add grouped native update helper.
- [ ] Replace per-component prior baseline.
- [ ] Integrate grouped MAP/VBEM steps with SQUAREM.
- [ ] Rebuild native extension.

Phase 4:

- [ ] Update locus diagnostics and output schema.
- [ ] Update tests and goldens intentionally.

Phase 5:

- [ ] Run focused sentinel tests.
- [ ] Run benchmark sweeps over strength and gDNA bias.
- [ ] Choose defaults from measured false-positive/false-negative tradeoffs.

Phase 6:

- [ ] Decide whether transcript-derived multi-group RNA priors are needed.
- [ ] Do not add gene-level internal priors.
