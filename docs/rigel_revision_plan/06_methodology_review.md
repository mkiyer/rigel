# Methodology Review

This document records the assessment of the March 2026 methodology critique and
the revisions it motivates.

## 1. Overall Assessment

The critique is directionally correct and worth acting on.

It does not invalidate the theory-first redesign. It does expose three places
where the implementation plan was still too optimistic:

1. the internal gDNA parameterization needed to be revisited
2. the calibration bootstrap path was underspecified
3. the sparsity story needed to be stated more explicitly in the objective

## 2. Point-by-Point Assessment

### 2.1 Collapsed gDNA versus two internal gDNA components

Assessment:

- mathematically, the critique is right that a collapsed gDNA component with a
  free strand parameter is isomorphic to a two-component `g_plus/g_minus`
  representation
- architecturally, the critique is also right that `T + N + 2` may be easier to
  implement because it keeps the solver closer to a standard mixture update

Revision:

- keep one collapsed gDNA abundance in the public model and outputs
- adopt `T + N + 2` as the preferred implementation target

### 2.2 Chicken-and-egg calibration risk

Assessment:

- the critique is right that immediate whole-genome soft weighting is too
  ambitious as a first calibration pass
- the first implementation needs a calibration bootstrap that can be audited

Revision:

- start from a highly conservative gDNA-dominant seed set
- fit initial nuisance parameters from that seed set
- expand to softer regional weights only after the seed-based calibration is
  stable

### 2.3 Replacing heuristic gates with real sparsity

Assessment:

- the critique is right in principle: if the model uses non-sparse priors, hard
  pruning pressure has not really been removed
- the current repository already uses `prior_alpha = 0.01`, so Rigel is already
  numerically in a sparse-prior regime in code
- the theory docs were simply not explicit enough about this

Revision:

- make the redesign target explicit: optional components should use
  `0 < alpha_c < 1`
- preserve the distinction between true model sparsity and any later reporting
  threshold or numerical cleanup

### 2.4 Asymptotic behavior of the strand-symmetry prior

Assessment:

- the critique is correct that a fixed global `kappa_sym` is not a hard law of
  symmetry
- high-coverage loci can overwhelm the prior and express real strand imbalance

Interpretation:

- this is acceptable if `kappa_sym` is understood as an overdispersion prior,
  not a rigid physical constraint
- if the scientific stance changes and gDNA symmetry must remain nearly fixed
  even at very high coverage, then the prior family would need to change

Current recommendation:

- keep the current interpretation for now
- validate empirically whether high-coverage loci show plausible or pathological
  strand escape

## 3. Revised Implementation Stance

The revised near-term plan is:

1. build region partition and evidence extraction first
2. define a conservative gDNA-dominant seed set
3. fit initial gDNA nuisance parameters from that seed set
4. only then add soft weighting for calibration expansion
5. implement the main EM migration around the `T + N + 2` gDNA pair model while
  preserving collapsed public gDNA output

## 4. What Did Not Change

The critique does not change the core redesign goals.

We still want:

- separation of generative model from inference heuristics
- one collapsed public gDNA abundance output per locus
- empirical-Bayes calibration of gDNA nuisance parameters
- upstream nuisance learning rather than embedded ad hoc EM behavior