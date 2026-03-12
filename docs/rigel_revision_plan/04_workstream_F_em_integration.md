# Workstream F: Main EM Integration

This workstream should start only after the calibration bundle is stable.

## 1. Implementation Goal

Replace the current embedded gDNA heuristics inside the locus EM path with a
cleaner calibration-informed parameterization.

The target state is:

- two internal strand-specific gDNA components per locus
- one collapsed public gDNA abundance obtained by summing the gDNA pair
- one global calibration-derived `kappa_sym`
- no hidden reuse of initialization quantities as if they were direct abundance
  pseudo-counts

## 2. Existing Integration Points

The current EM plumbing is already localized enough to revise in stages.

Relevant touchpoints:

- `src/rigel/locus.py`
  - locus-level component construction
  - prior gating
  - current gDNA init semantics
- `src/rigel/estimator.py`
  - batch locus EM interface and result handling
- `src/rigel/native/em_solver.cpp`
  - actual M-step treatment of gDNA symmetry and priors
- `src/rigel/pipeline.py`
  - orchestration and prior preparation

## 3. Recommended Migration Strategy

### 3.1 First integrate the calibration bundle without changing solver math

Before introducing free `phi_l`, pass the new calibration outputs through the
pipeline and make them available to the locus builder and solver layer.

That gives a staging point where:

- old behavior can still be reproduced if needed
- diagnostics can compare old and new calibration values side by side

### 3.2 Then separate gDNA prior gating from abundance prior mass

The current audit already showed that `gdna_init` behaves mostly as a gating
signal. The redesign should make that explicit.

Recommended direction:

- replace `gdna_init` semantics with either:
  - a calibration-derived eligibility decision, or
  - a minimal explicit abundance prior unrelated to the old EB count proxy
- do not let a calibration-derived nuisance quantity masquerade as direct EM
  abundance evidence

### 3.3 Then migrate the solver to the gDNA pair model

Once calibration is stable, update the solver-side gDNA strand treatment so that
the EM explicitly carries `g_plus` and `g_minus` as ordinary mixture
components. Their total mass is reported as locus-level gDNA abundance, and
their implied strand fraction is regularized by the symmetric Beta prior
controlled by `kappa_sym`.

## 4. Step-by-Step Coding Plan

1. define the calibration bundle inputs expected by the locus EM path
2. thread the bundle through pipeline and estimator orchestration
3. expand the locus component layout to include `g_plus` and `g_minus`
4. simplify `gdna_init` and gDNA prior gating semantics in `locus.py`
5. update native solver logic to use the gDNA pair and symmetry prior
  structure
6. add regression tests comparing old and new behavior on pristine cases
7. add targeted tests for locus-level strand imbalance and gDNA-heavy cases

## 5. Main Risk

The main risk is trying to refactor solver math before the upstream calibration
outputs are trustworthy.

The safest rule is:

- no native EM changes until the region evidence and calibration diagnostics are
  already inspectable and stable