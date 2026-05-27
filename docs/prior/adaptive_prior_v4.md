# Adaptive Prior v4 - Evidence-Gated Quantile/ESS Prior

Date: 2026-05-26
Status: implementation-ready design
Supersedes: `docs/prior/adaptive_prior_v3.md`

## 1. Purpose

Adaptive prior v3 has the right high-level idea: one user-facing operating
point, a posterior over the gDNA share of unspliced mass, and prior strength
learned from uncertainty rather than from hand-tuned constants.

This v4 plan keeps that contract but fixes the places where v3 can still be
brittle:

- boundary cases like certain no-gDNA evidence must create strong RNA support,
  not collapse to zero effective sample size;
- global empirical-Bayes shrinkage must not inject prior mass into loci with no
  prior-owned unspliced mass;
- global shrinkage should borrow strength from other loci, not reinforce a
  locus with its own evidence;
- density fallback uncertainty must use the full four-state posterior, not only
  a binary captured/off-target shortcut;
- strand and density confidence levels used internally by calibration must be
  separated from the single user-facing prior decision knob.

The result is still simple:

```text
per locus:  U_l, E[D_l], Var[D_l]
         -> local Beta-like posterior over phi_l = D_l / U_l
         -> leave-one-locus-out sample shrinkage, capped by local data content
         -> share_l = posterior quantile at rna_confidence
         -> grouped EM alphas for one gDNA group and one aggregate RNA group
```

No EM solver ABI change is required beyond the grouped prior interface already
implemented for prior redesign v3.

## 2. User Contract

### 2.1 One Public Knob

Expose exactly one prior-policy option:

```text
--rna-confidence q
```

Range:

```text
0 < q < 1
```

Default:

```text
q = 0.5
```

Meaning:

```text
Use the q-th posterior quantile of the locus gDNA share as the grouped-prior
gDNA share.
```

Higher values are more conservative for RNA reporting because they summarize
the posterior by a higher plausible gDNA share. Lower values are more sensitive
for RNA.

### 2.2 Internal Calibration Intervals Are Not User Prior Knobs

The existing `rna_lower_confidence` and `gdna_density_confidence` values are
currently used inside calibration for RNA lower-bound screens and gDNA upper
diagnostics. In v4 they are no longer ordinary user-facing quantification
knobs.

Implementation rule:

- add canonical `CalibrationConfig.rna_confidence = 0.5`;
- keep private calibration interval constants at their current conservative
  values for strand screens and upper-bound diagnostics;
- do not let those internal interval constants determine adaptive-prior share
  or ESS;
- for one release, accept old CLI/YAML names as deprecated aliases to
  `rna_confidence` only, with a warning that they no longer tune internal
  calibration screens.

This avoids pretending that a strand-screen interval and a prior decision
quantile are the same object.

### 2.3 Removed Public Knobs

Reject these keys in YAML and remove them from CLI help:

```text
aggregate_prior_strength
aggregate_prior_edge_count
aggregate_prior_max_count
gdna_prior_logit_bias
```

They were useful to stabilize v3, but they are not scientific inputs a user can
set from a sample without ground truth.

`aggregate_prior_max_count` survives only as a private safety cap named
`MAX_ESS` in the helper module. It is a numerical guard and review point, not a
model parameter.

## 3. Core Objects

For each locus `l`:

```text
U_l      = prior-owned unspliced mass allocated to the locus
D_l      = latent gDNA-owned part of U_l
phi_l    = D_l / U_l, defined only when U_l > 0
mu_D_l   = E[D_l | calibration]
var_D_l  = Var[D_l | calibration]
```

The grouped EM prior remains aggregate-only:

```text
alpha_gdna_add_l
alpha_rna_add_l
```

All annotated mRNA and synthetic nRNA components remain one RNA group for prior
purposes. The prior does not split mature RNA from nascent RNA.

The native structural gDNA flag remains a hard technical gate:

```text
has_gdna_candidate_l = any unspliced EM unit with finite gDNA likelihood
```

If this is false, both grouped alphas are exactly zero.

## 4. Local Posterior From Moments

### 4.1 Moment Inputs

Adaptive v4 needs `PriorMassDeconvolution` to carry:

```python
gdna_unspliced_var: np.ndarray
```

The invariant is count-scale:

```text
0 <= var_D_l <= mu_D_l * (U_l - mu_D_l)
```

This is the maximum variance for a bounded count `D_l in [0, U_l]` with mean
`mu_D_l`. Clamp violations at construction time and set a diagnostic flag.

### 4.2 Boundary-Aware Beta-Like Construction

v3's main mathematical bug is the line `tau = 0 if v_l == 0`. That erases the
strongest possible evidence: a point-mass posterior at `phi = 0` or `phi = 1`.

Use this rule instead:

```text
p_l = clamp(mu_D_l / U_l, 0, 1)              if U_l > EPS_MASS else undefined
v_l = clamp(var_D_l / U_l^2, 0, p_l*(1-p_l)) if variance is finite
n_l = min(U_l, MAX_ESS)

if U_l <= EPS_MASS:
    tau_l = 0
elif variance is missing or non-finite:
    tau_l = 0
elif v_l == 0:
    tau_l = n_l
else:
    tau_l = clamp(p_l * (1 - p_l) / v_l - 1, 0, n_l)

a_local_l = tau_l * p_l
b_local_l = tau_l * (1 - p_l)
```

Boundary consequences:

- `p_l = 0, v_l = 0, U_l > 0` gives `Beta(0, n_l)`, a degenerate zero-gDNA
  posterior with real RNA ESS.
- `p_l = 1, v_l = 0, U_l > 0` gives `Beta(n_l, 0)`, a degenerate all-gDNA
  posterior.
- interior `v_l = 0` gives a very concentrated posterior capped by observed
  prior-owned mass and `MAX_ESS`.
- missing variance gives no local ESS. The locus may still borrow from the
  sample profile, but only if it has `U_l > 0`.

No alpha floors are returned. Degenerate boundary Betas are handled by
short-circuiting the quantile function.

## 5. Global Shrinkage

v3's inverse-variance global Beta is elegant in the interior but fails at the
important sentinels: if every usable locus says `phi = 0` with zero variance,
then `psi * (1 - psi) = 0`, so the global ESS can collapse to zero exactly when
we want strong no-gDNA evidence.

v4 replaces that with pooled local evidence plus an agreement shrinkage factor.

### 5.1 Pool Local Evidence

Use only loci with `U_l > EPS_MASS` and `tau_l > 0`:

```text
A = sum_l a_local_l
B = sum_l b_local_l
T = A + B
```

If `T == 0`, the global profile is empty.

Otherwise:

```text
psi = A / T
```

### 5.2 Shrink by Cross-Locus Agreement

Let `p_l = a_local_l / tau_l` for usable loci. Use `tau_l` as the weight:

```text
v_obs = sum_l tau_l * (p_l - psi)^2 / T
v_mea = sum_l tau_l * v_l / T

if v_obs <= v_mea:
    agreement = 1
else:
    agreement = v_mea / v_obs

kappa_pool = min(T * agreement, MAX_ESS)
```

Boundary behavior is now correct:

- all usable loci certain `phi = 0`: `v_obs = v_mea = 0`, agreement is 1,
  `psi = 0`, and `kappa_pool = min(T, MAX_ESS)`;
- all usable loci certain `phi = 1`: same, but `psi = 1`;
- heterogeneous capture samples: `v_obs` exceeds measurement variance, so the
  global ESS shrinks automatically;
- no usable loci: global ESS is zero and no global prior is invented.

### 5.3 Leave-One-Locus-Out Contribution

Do not add a locus's own evidence back to itself through the global profile.
For each locus:

```text
T_loo = T - tau_l
A_loo = A - a_local_l
B_loo = B - b_local_l

if T_loo > 0:
    psi_loo = A_loo / T_loo
    kappa_loo = kappa_pool * (T_loo / T)
else:
    psi_loo = 0
    kappa_loo = 0
```

Then cap the global contribution by the remaining local data-content budget:

```text
remaining_l = max(0, min(U_l, MAX_ESS) - tau_l)
kappa_global_l = min(kappa_loo, remaining_l)

a_global_l = kappa_global_l * psi_loo
b_global_l = kappa_global_l * (1 - psi_loo)
```

This is the main v4 robustness improvement. Global shrinkage can guide an
ambiguous locus, but it cannot:

- affect a locus with `U_l <= EPS_MASS`;
- push the final ESS above the locus's own prior-owned unspliced mass;
- reinforce a single-locus sample by recycling that locus's evidence;
- override a locally certain locus whose `tau_l` already consumes its data
  budget.

## 6. Final Prior Summary

For loci with `U_l > EPS_MASS` and `has_gdna_candidate_l`:

```text
a_post_l = a_local_l + a_global_l
b_post_l = b_local_l + b_global_l
ess_l    = a_post_l + b_post_l
```

Decision summary:

```text
share_l = Q_q(Beta(a_post_l, b_post_l))
```

where `q = rna_confidence`.

Degenerate cases:

```text
ess_l == 0       -> share_l = 0, alphas = 0
a_post_l == 0    -> share_l = 0
b_post_l == 0    -> share_l = 1
otherwise        -> scipy.special.betaincinv(a_post_l, b_post_l, q)
```

Return:

```text
alpha_gdna_add_l = ess_l * share_l
alpha_rna_add_l  = ess_l * (1 - share_l)
```

If `has_gdna_candidate_l` is false, return both alphas as zero. An RNA prior is
not needed when there is no competing gDNA component in the native EM state.

## 7. Variance Propagation

### 7.1 Strand Deconvolution Variance

Correct file location: `RegionGdnaEstimate` and `RegionGdnaChannelEstimate`
live in `src/rigel/calibration/strand_deconv.py`, not `_result.py`.

Add to `RegionGdnaEstimate`:

```python
var_count: np.ndarray
```

Both exact and normal paths already compute `sd_r` for RNA count `R`. Since
`D = N - R` at fixed observed total `N`:

```text
Var(D) = Var(R) = sd_r^2
```

Add to `RegionGdnaChannelEstimate`:

```python
contained_var: np.ndarray
boundary_left_var: np.ndarray
boundary_right_var: np.ndarray
```

For strand-informative calibration:

```text
mu_D_region  = contained_mean + boundary_left_mean + boundary_right_mean
var_D_region = contained_var  + boundary_left_var  + boundary_right_var
```

This treats the three compartments as approximately independent. Record that
assumption in summary diagnostics.

Near-unstranded strand paths should not emit fake zero variance. If strand is
not identifiable, the production prior path should use density fallback
variance instead.

### 7.2 Density Fallback Variance From Four States

The density fallback must use the full posterior over the four calibration
states:

```text
0 background
1 gdna_only_capture
2 expressed_capture
3 expressed_offtarget
```

For each region `r` and state `s`, define a density-predicted gDNA count mean
`m_rs` and variance `v_rs`.

Off-target states (`background`, `expressed_offtarget`) use the background
Gamma-Poisson predictive:

```text
alpha = background.rho_off_alpha
beta  = background.rho_off_beta
L     = observation.contained_leff_r

m_off = alpha * L / beta
v_off = m_off + m_off^2 / alpha
```

Captured states (`gdna_only_capture`, `expressed_capture`) use the background
plus swept excess Gamma posterior:

```text
alpha = background.rho_off_alpha + sweep.alpha_excess_r
beta  = background.rho_off_beta  + sweep.beta_excess_r
L     = observation.contained_leff_r

m_cap = alpha * L / beta
v_cap = m_cap + m_cap^2 / alpha
```

Then apply the law of total variance over the four states:

```text
mu_D_r  = sum_s p_rs * m_rs
var_D_r = sum_s p_rs * (v_rs + (m_rs - mu_D_r)^2)
```

This collapses to the current captured/off-target mixture when the two
captured states share a predictive distribution and the two off-target states
share one. The explicit four-state form is still preferred because it preserves
state uncertainty, is easier to audit, and will remain correct if expression
state later changes the gDNA predictive.

### 7.3 Projection From Regions To Loci

Use the same geometry shares that allocate means. If region `r` contributes
share `s_lr` to locus `l`:

```text
U_l      = sum_r s_lr * U_r
mu_D_l   = sum_r s_lr * mu_D_r
var_D_l  = sum_r s_lr^2 * var_D_r
```

The squared-share variance rule is exact under independent region errors. It
also gives the right limiting behavior when one region is split across several
loci: the total allocated variance does not exceed the source region variance.

If a region lacks finite variance, mark the affected loci with
`PRIOR_NO_LOCAL_VARIANCE` and omit that region's variance from local ESS. Do
not invent interval-derived variance from an upper bound as a fallback unless a
future design explicitly models that interval-generating process.

## 8. Four-State Dry Run

### 8.1 Background / Off-Target Not Expressed

Expected signals:

- no spliced expression evidence;
- no capture excess;
- off-target background density only.

Behavior:

- `mu_D` is the background Gamma-Poisson expectation;
- `var_D` is broad when the background model is prior-only or sparse;
- local ESS is small unless actual observed unspliced mass and variance support
  a sharp estimate;
- if `U_l == 0`, no prior is emitted;
- in true zero-gDNA samples, many loci with point-mass or near-zero evidence
  pool into a strong global RNA-protective profile.

### 8.2 gDNA-Only Capture

Expected signals:

- capture/boundary excess;
- little or no expression evidence;
- no strand RNA lower-bound support.

Behavior:

- captured-state density raises `mu_D`;
- if variance is small relative to `U`, local ESS rises;
- `share_l` is high, especially at high `rna_confidence`;
- the gDNA prior does not collapse just because other off-target loci are
  clean, because strong local ESS consumes the local data budget before global
  shrinkage can dominate.

### 8.3 Expressed Capture

Expected signals:

- capture/boundary excess;
- expression evidence from spliced counts or strand RNA support;
- possible high unspliced RNA/nRNA mass.

Behavior:

- in stranded data, strand deconvolution supplies a direct `mu_D` and
  `var_D`, so expressed RNA can protect itself even in captured regions;
- in unstranded data, density may see capture but cannot fully distinguish
  captured gDNA from expressed unspliced RNA, so state uncertainty appears in
  `var_D` and reduces ESS;
- the prior should be cautious rather than falsely decisive. Ambiguous
  expressed-capture loci remain likelihood-driven.

### 8.4 Expressed Off-Target

Expected signals:

- expression evidence;
- no capture excess;
- background density only.

Behavior:

- density fallback gives low `mu_D`;
- strand evidence, if available, further concentrates the posterior toward
  RNA;
- local or global RNA prior protects against false gDNA siphoning in pure-RNA
  single-exon and wide-intron sentinels.

## 9. Library-Regime Dry Run

### 9.1 Stranded Whole-Transcriptome RNA-Seq

This is the easiest regime. Strand deconvolution supplies count-scale means and
variances for unspliced compartments.

Expected v4 behavior:

- zero-gDNA loci: `p_l` near 0, often low variance, strong RNA prior;
- true gDNA loci: `p_l` high with finite ESS, strong gDNA prior;
- mixed loci: posterior uncertainty lowers ESS instead of requiring a tuned
  logit bias.

### 9.2 Unstranded Whole-Transcriptome RNA-Seq

This is intrinsically harder because unspliced RNA and gDNA can be
observationally similar.

Expected v4 behavior:

- density fallback carries broad variance when evidence is weak;
- sample-wide no-gDNA evidence can still produce RNA protection through the
  boundary-aware global pool;
- true nRNA in unstranded data should not become a high-confidence gDNA prior
  unless density evidence is also precise;
- unresolved identifiability remains visible as low ESS, not hidden behind a
  tuned bias.

### 9.3 Stranded Capture RNA-Seq

Capture raises the stakes because density enrichment alone is not enough to
separate captured gDNA from captured expressed RNA.

Expected v4 behavior:

- strand evidence protects expressed captured RNA;
- gDNA-only captured regions receive high local gDNA share;
- cross-locus heterogeneity shrinks the global ESS, preventing capture samples
  from imposing one global gDNA rate everywhere.

### 9.4 Unstranded Capture RNA-Seq

This is the hardest supported regime.

Expected v4 behavior:

- four-state density variance is essential;
- expressed-capture ambiguity should reduce ESS rather than forcing either RNA
  or gDNA;
- strong sample-wide profiles can help only up to the locus's own `U_l` data
  budget;
- remaining failures should be diagnosed as identifiability limits, not tuned
  away with hidden constants.

## 10. Helper Module

Add:

```text
src/rigel/calibration/adaptive_prior.py
```

Public pure functions:

```python
def local_beta_from_moments(
    mu_D: np.ndarray,
    var_D: np.ndarray,
    U: np.ndarray,
    *,
    max_ess: float = MAX_ESS,
) -> LocalPriorMoments:
    """Return p, v, a_local, b_local, tau, flags."""

def global_profile_from_local(
    local: LocalPriorMoments,
    *,
    max_ess: float = MAX_ESS,
) -> GlobalPriorProfile:
    """Return pooled psi, kappa_pool, agreement, and leave-one-out arrays."""

def compute_grouped_prior_counts(
    *,
    mu_D: np.ndarray,
    var_D: np.ndarray,
    U: np.ndarray,
    has_gdna_candidate: np.ndarray,
    rna_confidence: float,
    max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    """Return grouped alphas and diagnostics."""
```

Dataclasses should be frozen and slot-based. The module has no I/O and no
dependency on pipeline objects.

Private constants:

```python
MAX_ESS = 3000.0
EPS_MASS = 1.0e-12
EPS_NUM = 1.0e-12
```

Only `MAX_ESS` is a policy safety cap. The epsilons are numerical guards.

## 11. Output Schema

### 11.1 `loci.feather`

Add or keep these columns:

```text
prior_gdna_share_mean       # local posterior mean p_l
prior_gdna_share_var        # local share variance v_l
prior_gdna_share_quantile   # final decision share_l
prior_ess_local             # tau_l
prior_ess_global            # kappa_global_l after LOO and local cap
prior_ess_final             # ess_l
prior_global_share          # psi_loo used for this locus
prior_global_agreement      # sample agreement scalar
alpha_gdna_add              # final gDNA grouped alpha
alpha_rna_add               # final aggregate-RNA grouped alpha
prior_flags                 # bitfield
prior_conflict_score        # post-EM diagnostic, NaN before EM result
```

Keep `gdna_prior_count_em` as an output compatibility alias for
`alpha_gdna_add` for one release. Remove `prior_gdna_share_biased` and
`prior_budget*` from new tests and documentation; if retained temporarily in
the dataframe, mark them deprecated and populate from the new columns only for
compatibility.

### 11.2 `summary.json`

Add:

```json
"prior_policy": {
  "name": "adaptive_quantile_v4",
  "rna_confidence": 0.5,
  "max_ess": 3000.0,
  "global_gdna_share": 0.04,
  "global_ess_pool": 1200.0,
  "global_agreement": 0.71,
  "n_loci_total": 12345,
  "n_loci_with_prior_mass": 11890,
  "n_loci_with_local_ess": 9120,
  "n_loci_structural_gated": 410,
  "flag_histogram": {"PRIOR_NO_LOCAL_VARIANCE": 453}
}
```

## 12. Flags

Use a single `prior_flags` bitfield:

```text
0x001 PRIOR_NO_LOCAL_VARIANCE       missing/non-finite var_D_l
0x002 PRIOR_VARIANCE_CLAMPED        var_D_l exceeded bounded-count ceiling
0x004 PRIOR_BOUNDARY_DEGENERATE     v_l == 0 and U_l > 0
0x008 PRIOR_STRUCTURAL_GATED        no native gDNA candidate
0x010 PRIOR_NO_UNSPLICED_MASS       U_l <= EPS_MASS
0x020 GLOBAL_PROFILE_EMPTY          no usable local evidence anywhere
0x040 GLOBAL_LOO_EMPTY              no other locus available for shrinkage
0x080 PRIOR_FINAL_ESS_CAPPED        local/global ESS hit local or MAX_ESS cap
```

Flags explain why a prior is weak or absent. They never change the math after
being set.

## 13. Implementation Phases

### Phase A - Helper Math Only

Add `adaptive_prior.py` and unit tests.

Required tests:

1. zero total gives zero local and final alphas;
2. missing variance gives zero local ESS and sets `PRIOR_NO_LOCAL_VARIANCE`;
3. `p=0, var=0, U>0` gives boundary `Beta(0, n)` and nonzero RNA ESS;
4. `p=1, var=0, U>0` gives boundary `Beta(n, 0)` and nonzero gDNA ESS;
5. increasing variance lowers local ESS;
6. local ESS never exceeds `min(U, MAX_ESS)`;
7. all-zero global profile gives `psi=0`, positive `kappa_pool`, agreement 1;
8. heterogeneous loci shrink global agreement below 1;
9. one-locus samples have empty leave-one-out global contribution;
10. global contribution never makes final ESS exceed `min(U, MAX_ESS)`;
11. `share(q)` is monotone in q;
12. structural gating zeroes both alphas exactly.

### Phase B - Calibration Variance

- Add `var_count` to `RegionGdnaEstimate` in `strand_deconv.py`.
- Add compartment variances to `RegionGdnaChannelEstimate`.
- Add `gdna_unspliced_var` to `PriorMassDeconvolution`.
- Implement four-state density fallback variance in `calibration_iteration.py`.
- Extend calibration summaries with variance stats and clamp histograms.

Required tests:

- exact strand posterior variance equals the PMF variance;
- normal strand path round-trips `sd_r^2`;
- near-unstranded calibration does not emit fake precise variance;
- Gamma-Poisson variance matches `mean + mean^2 / alpha`;
- four-state variance includes between-state uncertainty;
- bounded-count variance clamp sets a flag.

### Phase C - Projection And Prior Assembly

- Project variance through `_allocate_counts_by_geometry()` with squared shares.
- Replace `_compute_grouped_prior_counts()` body with the adaptive helper.
- Feed `has_gdna_candidate` into the helper rather than gating only after the
  helper returns.
- Extend `PriorTable` and locus output diagnostics.

Required tests:

- disjoint regions preserve variance sums;
- split regions use squared shares;
- `U_l == 0` receives no global prior;
- all-spliced loci remain structurally gated;
- zero-gDNA sentinel priors produce RNA alpha without logit bias;
- high-gDNA sentinel priors produce nonzero gDNA alpha where local evidence is
  informative.

### Phase D - Config And CLI

- Add `--rna-confidence`.
- Add `CalibrationConfig.rna_confidence` with validation `0 < q < 1`.
- Remove prior-strength/logit fields from `EMConfig` and config registry.
- Keep old confidence names for one release as deprecated aliases to
  `rna_confidence`, warning that internal calibration intervals are no longer
  user-tuned.
- Reject removed prior constants in YAML with a migration error.

### Phase E - Output And Acceptance Sweep

- Update `loci.feather`, summary JSON, docs, and goldens.
- Add a q-sweep over `q in {0.1, 0.5, 0.9}`.

For every locus with fixed calibration evidence and `prior_ess_final > 0`:

```text
share(q=0.1) <= share(q=0.5) <= share(q=0.9)
alpha_gdna_add non-decreasing in q
alpha_rna_add non-increasing in q
```

Run the sweep over the four library regimes:

```text
stranded whole-transcriptome
unstranded whole-transcriptome
stranded capture
unstranded capture
```

## 14. Acceptance Criteria

Implementation is accepted when:

1. full tests pass, apart from explicitly documented unrelated pre-existing
   failures;
2. `rigel quant --help` exposes `--rna-confidence` as the only prior operating
   point;
3. removed prior constants are rejected from YAML;
4. no code path references `gdna_prior_logit_bias`;
5. zero-gDNA sentinels recover RNA without user-tuned bias or edge counts;
6. true high-gDNA sentinels keep nonzero gDNA alpha where calibration evidence
   is informative;
7. unstranded nRNA-heavy sentinels show low ESS when evidence is intrinsically
   ambiguous rather than high-confidence false gDNA;
8. capture sentinels show lower global agreement than homogeneous
   non-capture samples when cross-locus gDNA shares are heterogeneous;
9. output diagnostics make every zero or weak prior explainable by flags;
10. the native EM grouped-prior behavior remains unchanged except for the new
    alpha values it receives.

## 15. Non-Goals

- Splitting annotated mRNA and nRNA into separate prior groups.
- Changing transcript-first internal inference.
- Auto-selecting `rna_confidence` from the data.
- Replacing the EM solver or SQUAREM implementation.
- Solving fundamentally unidentifiable unstranded capture cases by hidden
  tuning. v4 should expose ambiguity as low ESS and diagnostics.

## 16. Summary

Adaptive v4 keeps the one-knob quantile design but makes the uncertainty model
honest at the boundaries and in capture data.

The key simplifications are:

- use local moment-matched evidence, including degenerate boundary evidence;
- pool local evidence for the sample profile instead of fitting a fragile
  inverse-variance Beta at `psi=0` or `psi=1`;
- shrink the global profile by observed cross-locus disagreement;
- apply global shrinkage leave-one-locus-out and cap it by the locus's own
  prior-owned unspliced mass;
- propagate density variance through the full four-state calibration tensor.

This removes the v3 magic constants from the user surface without replacing
them with a new layer of trust multipliers. Variance and data content do the
work.