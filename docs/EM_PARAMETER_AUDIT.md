# Rigel EM Parameter Audit

This note describes the EM initialization and prior system as it exists in the
current code, not as it existed in older design notes.

It is intended to answer four questions:

1. What parameters affect the locus EM?
2. Where do they enter the pipeline?
3. Which values act as true priors versus warm starts versus eligibility gates?
4. Which parts are conceptually redundant or easy to misread?

## 1. Current locus model

For a locus with `T` transcripts and `N` unique nRNA spans, the current solver
uses:

- `T` mRNA components
- `N` nRNA components
- `1` gDNA component

Total components:

`T + N + 1`

This is implemented in the native batch solver and in per-locus extraction.

## 2. High-level flow

The current EM-related path is:

1. Scan BAM and accumulate transcript- and nRNA-level evidence.
2. Build ambiguous-fragment CSR with mRNA and nRNA candidates; gDNA is added per
   locus later.
3. Build loci as connected components over transcript compatibility.
4. Compute `nrna_init` from intronic strand excess.
5. Compute per-locus `gdna_init` from the hierarchical EB gDNA density model.
6. Compute per-nRNA Beta prior parameters `nrna_frac_alpha` and
   `nrna_frac_beta` from the hierarchical EB nascent-fraction model.
7. For each locus, extract a local sub-problem, gate ineligible components, build
   a coverage-weighted warm start and OVR prior, then run SQUAREM-accelerated
   MAP-EM or VBEM.
8. Optionally prune weak components and run one final redistribution step.

The important point is that `nrna_init` and `gdna_init` are no longer injected
as additive pseudo-counts into the EM count totals. They mainly control whether
some components are eligible at all.

## 3. Parameter inventory by role

### 3.1 Global EM behavior

- `prior_alpha`
  - Flat pseudo-count per eligible component.
  - Used in OVR prior construction, not as the only prior mass.
  - Default: `0.01`.

- `prior_gamma`
  - Scale factor for the OVR prior.
  - If `0`, the prior becomes flat `alpha` over eligible components.
  - Default: `1.0`.

- `iterations`
  - Maximum EM iteration budget.
  - SQUAREM divides the practical number of outer iterations by `3` because one
    SQUAREM iteration uses multiple EM evaluations.

- `convergence_delta`
  - Stopping threshold on the change in the normalized state.

- `mode`
  - `map` uses normalized mixture weights `theta`.
  - `vbem` uses Dirichlet parameters `alpha` and digamma expectations.

- `prune_threshold`
  - Post-EM pruning threshold.
  - A component with zero deterministic support is pruned if
    `(unambig + em) / alpha_out < prune_threshold`.
  - gDNA is never pruned.

- `confidence_threshold`
  - Only affects reporting of high-confidence RNA assignments.
  - Does not affect the EM fit itself.

### 3.2 nRNA initialization and gating

- `nrna_init`
  - Computed from intronic strand excess:
    `(sense_intronic - antisense_intronic) / (2 * SS - 1)`.
  - Clamped at `0`.
  - Forced to `0` when strand specificity is too weak or when the nRNA span has
    no intronic territory.
  - Current role: eligibility gate only. If `nrna_init == 0`, the local nRNA
    prior is zeroed and the component is effectively disabled.

This is a major conceptual change from an older design where `nrna_init` was a
direct additive seed in the M-step.

### 3.3 nRNA Beta prior hierarchy

The nRNA prior system computes a mean nascent fraction per unique nRNA span,
then converts it to a Beta prior passed into the native M-step.

- `nrna_frac_kappa_global`
  - Shrinkage strength for locus-strand estimates toward the global mean.
  - `None` means auto-estimate with Method of Moments.

- `nrna_frac_kappa_locus`
  - Shrinkage strength for per-nRNA estimates toward their locus-strand parent.
  - `None` means auto-estimate.

- `nrna_frac_kappa_nrna`
  - Final fixed effective sample size of the Beta prior delivered to the native
    solver.
  - This controls how strongly the hierarchical prior mean constrains the EM.
  - Default: `5.0`.

- `nrna_frac_mom_min_evidence_global`
  - Minimum evidence for a locus-strand group to participate in automatic
    `kappa_global` estimation.

- `nrna_frac_mom_min_evidence_locus`
  - Minimum evidence for a per-nRNA estimate to participate in automatic
    `kappa_locus` estimation.

- `nrna_frac_kappa_min`
- `nrna_frac_kappa_max`
- `nrna_frac_kappa_fallback`
- `nrna_frac_kappa_min_obs`
  - Shared Method-of-Moments controls used by both nRNA and gDNA EB layers.

The current mean estimator is hybrid:

- density component from exonic and intronic densities after subtracting global
  gDNA density
- strand component from sense minus antisense divided by `2 * SS - 1`
- blend weight `W = (2 * SS - 1)^2`

Final Beta parameters are:

- `alpha_n = eta_n * nrna_frac_kappa_nrna`
- `beta_n = (1 - eta_n) * nrna_frac_kappa_nrna`

### 3.4 gDNA EB hierarchy

The gDNA prior system computes a per-locus expected gDNA count from a density
model. It is hierarchical:

- global density
- reference-level density shrunk toward global
- locus density shrunk toward parent reference

Relevant knobs:

- `gdna_kappa_ref`
  - Shrinkage strength for reference toward global.
  - `None` means auto-estimate.

- `gdna_kappa_locus`
  - Shrinkage strength for locus toward reference.
  - `None` means auto-estimate.

- `gdna_mom_min_evidence_ref`
- `gdna_mom_min_evidence_locus`
  - Minimum evidence thresholds for automatic Method-of-Moments estimation.

The raw gDNA density estimate is also hybrid:

- strand-based density from unspliced sense and antisense imbalance
- intergenic background density
- blend weight `W = (2 * SS - 1)^2`

Final per-locus value is:

- `gdna_init = shrunk_density * locus_exonic_bp`

Current role of `gdna_init`:

- it does not enter `unambig_totals`
- it does not directly seed `theta`
- it gates whether the gDNA component is eligible at all

If `gdna_init == 0`, the locus gDNA prior is zero and the component is disabled.

### 3.5 gDNA strand-symmetry penalty

These parameters act inside the native M-step, not in the Python prior builder.

- `strand_symmetry_kappa`
  - Controls how strongly asymmetric gDNA posterior mass is penalized.
  - Values `<= 2` disable the penalty.

- `strand_symmetry_pseudo`
  - Pseudo-count used when estimating the gDNA sense fraction before computing
    the penalty.

Important nuance: the solver does not penalize all gDNA mass equally. It
protects the symmetric baseline `2 * min(sense, anti)` and discounts only the
asymmetric excess. This makes the penalty much more targeted than a naive
Beta prior on total gDNA mass.

## 4. What the warm start actually does

After a locus is extracted, the native solver computes two related vectors from
coverage weights over ambiguous units.

### 4.1 `theta_init`

Initialized as:

- deterministic mRNA row sums from `unambig_counts`
- plus one coverage-normalized share from each ambiguous unit across eligible
  components

This is a warm start, not a prior.

### 4.2 `prior`

For each eligible component:

- `prior[c] = prior_alpha + gamma * OVR_share[c]`

where `OVR_share[c]` is the component's normalized share of the one-virtual-read
coverage distribution.

So the operational regularizer is not just `alpha`; it is `alpha + OVR`.

## 5. What the MAP M-step actually does

The native solver does not treat all components independently when `n_nrna > 0`.

For each unique nRNA group:

1. Sum mature counts over all transcripts sharing that nRNA.
2. Combine with the nRNA count to get total group mass.
3. Update the nascent fraction with a MAP Beta posterior:
   `eta = (c_nrna + a - 1) / (c_group + a + b - 2)`.
4. Assign `eta * c_group` to the nRNA component.
5. Redistribute `(1 - eta) * c_group` back to the mature transcripts in
   proportion to their transcript-specific mass.

This means the Beta prior is not just a soft additive count on the nRNA slot.
It directly constrains the mRNA versus nRNA partition of each shared group.

## 6. Eligibility rules before EM

The following components can be disabled before EM starts.

### 6.1 gDNA disabled when

- `gdna_init == 0`

### 6.2 nRNA disabled when

- all transcripts sharing that nRNA are single-exon
- `nrna_init == 0`

Ineligible components receive:

- `prior = 0`
- `eligible = 0`
- `unambig_totals = 0`

and therefore do not receive OVR mass or warm-start mass.

## 7. Practical interpretation of the current knobs

If we strip the system down to first-order effects, the dominant EM controls are:

1. scoring log-likelihoods for each fragment/component pair
2. OVR warm start and OVR prior via `prior_alpha` and `prior_gamma`
3. nRNA eligibility via `nrna_init`
4. gDNA eligibility via `gdna_init`
5. nRNA mRNA split constraint via `nrna_frac_alpha` and `nrna_frac_beta`
6. gDNA strand-asymmetry discount via `strand_symmetry_kappa` and
   `strand_symmetry_pseudo`
7. pruning via `prune_threshold`

Everything else is upstream machinery used to estimate those quantities.

## 8. Important code-level clarifications

These points are easy to misunderstand when reading only older docs.

### 8.1 `gdna_init` is not a seed count anymore

It is better thought of as a component-on/component-off decision derived from a
hierarchical EB density estimate.

### 8.2 `nrna_init` is not the final nRNA prior strength

It gates component eligibility. The actual nRNA versus mRNA prior used in the
solver is the Beta prior built from `nrna_frac_alpha` and `nrna_frac_beta`.

### 8.3 `prior_alpha` is not the whole prior

The actual per-component prior mass in MAP mode is:

- flat `alpha`
- plus OVR mass scaled by `gamma`

for eligible components only.

### 8.4 gDNA and nRNA are regularized in very different ways

- gDNA is regularized by eligibility, OVR, and the strand-symmetry penalty.
- nRNA is regularized by eligibility, OVR, and a structured Beta prior on the
  group nascent fraction.

### 8.5 Pruning uses post-fit evidence relative to `alpha_out`

So pruning is not purely data-driven; it depends on the combined effective prior
and fitted counts.

## 9. Likely simplification targets

If the solver is to be redesigned from first principles, the parts most worth
questioning are:

1. Whether `gdna_init` and `nrna_init` should be pure eligibility gates, or true
   probabilistic initializers.
2. Whether OVR should be both a warm start and a prior source.
3. Whether the nRNA Beta prior should operate on the current group-partition
   parameterization or on a simpler latent abundance model.
4. Whether post-EM pruning is compensating for avoidable prior or initialization
   issues earlier in the pipeline.
5. Whether gDNA needs both hierarchical EB gating and a strand-asymmetry penalty,
   or whether one of those mechanisms can absorb the other.

## 10. Source locations

The current implementation is primarily in:

- `src/rigel/config.py`
- `src/rigel/pipeline.py`
- `src/rigel/locus.py`
- `src/rigel/priors.py`
- `src/rigel/estimator.py`
- `src/rigel/native/em_solver.cpp`
