# Purity Model Options

This document addresses the open question you raised: what should the first
practical purity approximation look like?

## 1. Recommendation from the Current Code

Based on the existing architecture and the methodology critique, the best first
practical approximation is not a full latent continuous model and not a global
unsupervised soft-weighting pass over all regions.

The strongest first implementation is a seeded approach:

1. define a highly conservative set of gDNA-dominant seed regions
2. fit initial nuisance parameters from that seed set
3. fit a two-state regional purity model whose posterior state probabilities are
   the calibration weights
4. only then consider richer soft expansions beyond the seed set

Concretely:

1. build the region evidence table first
2. define a strict seed subset with minimal RNA evidence
3. compute a small set of interpretable RNA-likeness features per region
4. fit a two-state mixture over gDNA-dominant versus RNA-contaminated regions
5. use posterior state probabilities as gDNA-dominance weights
6. treat that weight explicitly as an approximation to the broader latent
   purity target, not as the final theoretical model

## 2. Why This Is the Best First Step

It fits the current codebase because:

1. the evidence features are already close to observable
2. a conservative seed set gives a stable starting point for calibration
3. a simple score can then be diagnosed easily on synthetic and real runs
4. the same interface can later be replaced by a richer continuous latent model
   without changing downstream calibration consumers
5. it avoids prematurely hard-coding a mixture model before the empirical
   regional feature distributions are known

## 3. Suggested First Approximation

### 3.1 Region-level feature vector

For each region, compute a compact feature vector such as:

- `splice_rate`
- `unspliced_density`
- `strand_imbalance`
- `context_class`
- optional local exon-versus-intron contrast features

### 3.2 Conservative seed set

The first calibration pass should start from regions that satisfy very strict
gDNA-dominance criteria.

Examples of first-pass seed criteria:

- no local splice evidence
- compatible context class such as intergenic or unambiguous intronic
- low RNA-like strand asymmetry
- sufficient unspliced coverage for stable nuisance estimation

This does not need to be perfect. It needs to be conservative enough that the
initial calibration is not contaminated by evidently RNA-rich regions.

### 3.3 Two-state posterior weighting rule

Then define a conservative two-state model in which:

- one component represents gDNA-dominant regions
- one component represents RNA-contaminated regions
- splice, density, and strand features shift posterior mass toward one state or
   the other

The operational weight is then:

$$
w_r = \Pr(u_r = \mathrm{gDNA\text{-}dominant} \mid y_r)
$$

### 3.4 Preserve the right interface

Even in the first version, the interface should already be:

- input: region evidence table
- output: one `w_r in [0,1]` per region

That keeps the downstream calibration layer agnostic to how purity was
approximated.

## 4. Why Seed-Then-Expand Is Better Than Immediate Global Soft Weighting

The main reason is the chicken-and-egg problem.

If the very first gDNA calibration depends on weights that themselves depend on
untrusted nuisance parameters, then calibration quality becomes hard to audit.

A seed-then-expand strategy is cleaner:

1. fit an initial calibration on regions that are very likely to be gDNA-like
2. inspect the resulting nuisance parameters
3. optionally expand to softer weights after that initial fit is stable

## 5. Better Than a Hard Threshold

This approach is much better than a hard yes or no region filter because it:

1. uses information from weakly informative regions without overtrusting them
2. avoids brittle cutoffs that will vary across library types
3. provides immediate diagnostics when weights look wrong

## 6. Better Than a Full Mixture Model as the Very First Implementation

A full mixture model is still the cleaner theoretical target. It is just not
the best very first coding step.

Reasons:

1. it requires choosing likelihood families before seeing the empirical feature
   distributions
2. it adds identifiability and fitting complexity before the evidence tables are
   debugged
3. it becomes harder to tell whether failures come from evidence extraction or
   from the mixture model itself

## 7. Concrete Suggested Version 1

If I were implementing the first practical approximation in this repository, I
would do this:

1. compute raw per-region features
2. define a very conservative seed set for initial calibration
3. fit initial nuisance parameters on that seed set
4. fit a two-state model for gDNA-dominant versus RNA-contaminated regions
5. use posterior state probabilities as `w_r`
6. inspect the induced weights on synthetic datasets before fitting any more
   ambitious latent model

That is the most pragmatic option suggested by the current code.

## 8. Upgrade Path

Once the region evidence table and diagnostic plots exist, the two-state
weights can later be replaced by either:

1. a richer semiparametric weighting model
2. a continuous latent purity model with posterior mean weights

Because the downstream interface stays fixed, that upgrade path is clean.