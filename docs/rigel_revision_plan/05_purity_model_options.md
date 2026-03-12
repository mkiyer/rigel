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
3. only then introduce a transparent monotone scoring model that produces soft
   region weights for expansion beyond the seed set

Concretely:

1. build the region evidence table first
2. define a strict seed subset with minimal RNA evidence
3. compute a small set of interpretable RNA-likeness features per region
4. transform those features into a soft gDNA-dominance weight in `[0, 1]`
5. treat that weight explicitly as an approximation to posterior purity, not as
   the final theoretical model

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

### 3.3 Monotone soft weighting rule

Then define a conservative hand-specified score:

- more splice evidence lowers the weight
- stronger strand asymmetry lowers the weight
- strongly exon-enriched density lowers the weight
- intergenic and unambiguous-intronic contexts start from a higher baseline

This can be implemented as a logistic or clipped linear score with a small
number of interpretable coefficients.

### 3.4 Preserve the right interface

Even if the first version is heuristic, the interface should already be:

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
4. standardize features within broad context classes if necessary
5. define a soft score with explicit monotone signs
6. map the score through a logistic transform to get `w_r`
7. inspect the induced weights on synthetic datasets before fitting any more
   ambitious latent model

That is the most pragmatic option suggested by the current code.

## 8. Upgrade Path

Once the region evidence table and diagnostic plots exist, the score-based
weights can later be replaced by either:

1. a two-state mixture over RNA-rich versus gDNA-dominant regions
2. a continuous latent purity model with posterior mean weights

Because the downstream interface stays fixed, that upgrade path is clean.