# Workstream E: Empirical-Bayes Calibration

This workstream should be implemented before the main EM refactor.

## 1. Implementation Goal

Convert region evidence plus gDNA-dominance weights into a calibrated nuisance
parameter bundle for the run.

The first calibration bundle should contain at least:

- `kappa_sym`
- a weighted gDNA fragment-length model
- diagnostics describing the amount and distribution of effective calibration
  mass

## 2. Why This Fits the Current Code Well

The current code already has localized calibration-like logic:

- `src/rigel/priors.py` contains concentration estimation utilities
- `src/rigel/locus.py` contains current gDNA density and EB helper logic
- `src/rigel/frag_length_model.py` already supports weighted observation
- `src/rigel/pipeline.py` already has a natural place where upstream summaries
  are computed before the locus EM

That means the redesign can initially add a parallel calibration product rather
than immediately replacing everything.

## 3. Recommended First Implementation

### 3.1 Weighted `kappa_sym` estimation

Implement a dedicated estimator for the symmetric Beta-Binomial concentration
parameter using:

- per-region plus counts
- per-region total counts
- per-region gDNA-dominance weights

The estimator should be separate from the current nRNA `estimate_kappa()` helper
because the statistical object is different.

Recommended output:

- point estimate `kappa_sym`
- effective weighted number of regions
- effective weighted total fragments
- optional fit diagnostics

### 3.2 Weighted gDNA fragment-length model

Use the same region weights to train a gDNA-specific fragment-length model from
calibration-eligible unspliced fragments.

This is easy to support in the current code because `FragmentLengthModel.observe`
already accepts a floating-point weight.

Recommended implementation:

- keep a dedicated weighted gDNA histogram
- finalize it with the existing fragment-length model machinery
- report effective total weight and summary moments

### 3.3 Calibration bundle object

Add a simple container for calibrated nuisance outputs.

Suggested fields:

- `kappa_sym`
- `gdna_frag_length_model`
- `effective_region_weight`
- `effective_fragment_weight`
- optional diagnostics tables

## 4. Step-by-Step Coding Plan

1. define a calibration bundle dataclass
2. implement weighted `kappa_sym` fitting from the region evidence table
3. implement weighted gDNA fragment-length fitting from eligible fragments
4. add calibration diagnostics and summaries in the pipeline
5. validate on synthetic data where true gDNA symmetry dispersion and fragment
   length are known

## 5. Concrete Code Touchpoints

Most likely touchpoints:

- `src/rigel/priors.py`
- `src/rigel/frag_length_model.py`
- `src/rigel/pipeline.py`
- possibly a new module such as `src/rigel/gdna_calibration.py`

## 6. Suggested Separation of Responsibilities

Keep the roles clean.

- region evidence extraction should not fit nuisance parameters
- purity weighting should not update EM abundances
- calibration should return a bundle, not mutate the estimator directly
- the later EM integration layer should consume the bundle explicitly

That separation will make the redesign much easier to test and revise.