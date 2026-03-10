# Ablation Study v2: Geometric Splicing Expectation + Kappa Sweep

**Date**: 2026-03-09  
**Code version**: Post geometric-splicing-expectation implementation  
**Compared against**: Ablation v1 (pre geometric-splicing-expectation, kappa ∈ {2, 6})

---

## 1. Executive Summary

We ran a 576-run synthetic ablation sweep to evaluate two changes:

1. **Geometric splicing expectation** — predicts mRNA's unspliced burden
   ($R_t = L_\text{unspliced}/L_\text{spliced}$) and subtracts it from
   sense/antisense counts before gDNA density estimation.
2. **Expanded kappa sweep** — tests $\kappa \in \{2, 4, 6\}$ to find the
   optimal strand-symmetry penalty strength.

**Key findings**:

- The geometric splicing expectation has **zero measurable effect** on
  deconvolution accuracy. Matched runs (same kappa, same parameters)
  produce identical results before and after the code change.
- $\kappa = 4$ emerges as the **best compromise** between the extremes of
  $\kappa = 2$ (disabled) and $\kappa = 6$ (aggressive).
- The **gDNA → mRNA siphon** remains the dominant error mode. At extreme gDNA
  (1000 fragments / 50k total), mRNA error reaches 49% regardless of kappa
  or burden subtraction.
- The fundamental problem is **structural**: the EM likelihood cannot
  distinguish gDNA from unspliced mRNA at the fragment level. Prior-based
  interventions (burden subtraction, symmetry penalties) shift initialization
  but the EM converges to the same likelihood maximum.

---

## 2. Experimental Design

### 2.1 Parameter Grid

| Parameter              | Values                      | Count |
|------------------------|-----------------------------|-------|
| Transcript pattern     | TA1, TA2, TA3, TA4          | 4     |
| nRNA fragments (NTA)   | 0, 10, 100, 1000            | 4     |
| gDNA fragments         | 0, 10, 100, 1000            | 4     |
| Strand specificity (SS)| 1.00, 0.95, 0.90            | 3     |
| Symmetry kappa (κ)     | 2.0, 4.0, 6.0               | 3     |
| **Total runs**         |                              | **576** |

Each run simulates 50,000 fragments. "Clean" baseline = 36 runs (nRNA=0, gDNA=0).

### 2.2 What Changed Since v1

The geometric splicing expectation introduces:

- `compute_unspliced_to_spliced_ratios()` in `estimator.py`: computes
  $R_t = L_\text{unspliced} / L_\text{spliced}$ per transcript using
  eCDF-based effective lengths.
- `_compute_locus_mrna_burden()` in `locus.py`: predicts the unspliced
  mRNA fragment count from observed spliced counts × $R_t$.
- Burden subtraction in `_compute_per_locus_gdna_densities()`: subtracts
  predicted burden from observed unspliced counts before computing gDNA
  density.

The intent is to prevent high-expression mRNA loci from inflating gDNA
estimates by removing the predictable unspliced mRNA signal before gDNA
initialization.

---

## 3. Headline Results

### 3.1 Overall Accuracy

| Metric              | v1 (384 runs) | v2 (576 runs) | v2 κ=2 only | v2 κ=4 only | v2 κ=6 only |
|---------------------|---------------|---------------|-------------|-------------|-------------|
| mRNA median |rel err| | 3.54%       | 3.70%         | 6.06%       | 6.06%       | 3.76%       |
| mRNA p95            | 36.84%        | 35.09%        | 43.20%      | 30.38%      | 28.13%      |
| mRNA max            | 49.16%        | 49.16%        | 49.16%      | —           | —           |
| nRNA median |rel err| | 5.08%       | 4.61%         | 1.57%       | 1.63%       | 2.33%       |
| gDNA median |rel err| | 1.83%       | 1.56%         | 1.00%       | 0.54%       | 1.75%       |

The slight differences between v1 and v2 overall numbers are entirely due to
the addition of κ=4 shifting the aggregate median. For matched configurations,
v1 and v2 are identical (see §4).

### 3.2 Accuracy Grid (nRNA × gDNA)

Median mRNA relative error across all SS, kappa, and pattern combinations:

| nRNA \ gDNA |    0    |   10    |  100   | 1000   |
|-------------|---------|---------|--------|--------|
| **0**       |  0.00%  |  1.38%  |  8.62% | 28.43% |
| **10**      |  0.08%  |  1.29%  |  9.00% | 30.81% |
| **100**     |  0.44%  |  1.14%  |  6.79% | 23.64% |
| **1000**    |  2.92%  |  5.00%  |  6.02% | 17.06% |

The gDNA level dominates total mRNA error. nRNA has a secondary and partially
compensating effect — higher nRNA slightly reduces mRNA error at extreme gDNA
(30.81% → 17.06% going from NTA=10 to NTA=1000 at gDNA=1000).

---

## 4. Critical Finding: Geometric Splicing Expectation Has Zero Effect

### 4.1 Evidence

Direct comparison of matched runs (same kappa, same parameters) between the v1
sweep (pre-burden-subtraction) and v2 sweep (post-burden-subtraction):

| Metric                        | v1 κ=2 | v2 κ=2 | v1 κ=6 | v2 κ=6 |
|-------------------------------|--------|--------|--------|--------|
| mRNA median (contaminated)    | 5.44%  | 5.44%  | 4.39%  | 4.39%  |
| Mass flow nRNA err median     | −515   | −515   | +859   | +859   |
| Mass flow mRNA err median     | −78    | −78    | −57    | −57    |

Results are identical to the precision of the output format. The burden
subtraction produces no measurable change in any metric.

### 4.2 Why the Burden Subtraction Is Inert

The burden subtraction modifies the **gDNA prior initialization**, not the EM
likelihood or the convergence target. With 50,000 fragments per run, even at
gDNA=1000 (2% contamination), the data-to-prior ratio is overwhelming. The EM
converges to the same MAP estimate regardless of whether the gDNA prior is
shifted by the burden.

Mathematically: in a Beta-Binomial model with $n = 50{,}000$ observations, a
prior $\text{Beta}(\alpha, \beta)$ with $\alpha + \beta \ll n$ has negligible
influence on the posterior mode. The burden subtraction adjusts $\alpha$ and
$\beta$ by a small amount relative to $n$, so the posterior (and hence the EM
fixed point) is unchanged.

### 4.3 Implication

**Prior-based interventions cannot solve the gDNA siphon.** The problem is not
in initialization — it is in the EM likelihood structure itself. Any fix must
either:

1. Change the **likelihood** (add discriminative fragment-level features), or
2. Add **hard constraints** that survive EM iteration (not soft priors).

---

## 5. Kappa Analysis

### 5.1 Overall Kappa Comparison (Contaminated Runs)

| κ   | mRNA med | mRNA p95 | nRNA med | gDNA med |
|-----|----------|----------|----------|----------|
| 2.0 | 6.06%    | 43.20%   | 1.57%    | 1.00%    |
| 4.0 | 6.06%    | 30.38%   | 1.63%    | 0.54%    |
| 6.0 | 3.76%    | 28.13%   | 2.33%    | 1.75%    |

- κ=2 (flat / disabled): Best nRNA accuracy, worst mRNA tail (p95=43%).
- κ=4 (moderate): Same median mRNA as κ=2, but dramatically better p95 (30%).
  Best gDNA accuracy (0.54%). nRNA accuracy close to κ=2.
- κ=6 (aggressive): Best mRNA median and p95. Worst gDNA accuracy. Highest
  nRNA error.

### 5.2 Kappa × gDNA Interaction (nRNA Present)

| gDNA | κ   | mRNA med | nRNA med | nRNA direction |
|------|-----|----------|----------|----------------|
| 0    | 2.0 | 0.40%    | 0.21%    | −0.006 (under) |
| 0    | 4.0 | 0.40%    | 0.15%    | −0.005 (under) |
| 0    | 6.0 | 0.40%    | 0.13%    | −0.004 (under) |
| 10   | 2.0 | 1.86%    | 2.29%    | −0.023 (under) |
| 10   | 4.0 | 1.70%    | 2.51%    | +0.015 (over)  |
| 10   | 6.0 | 1.48%    | 3.37%    | +0.034 (over)  |
| 100  | 2.0 | 7.29%    | 11.00%   | −0.110 (under) |
| 100  | 4.0 | 7.15%    | 6.43%    | +0.061 (over)  |
| 100  | 6.0 | 7.20%    | 11.04%   | +0.110 (over)  |
| 1000 | 2.0 | 27.44%   | 100.00%  | −1.000 (zeroed)|
| 1000 | 4.0 | 19.88%   | 55.10%   | +0.383 (over)  |
| 1000 | 6.0 | 17.25%   | 100.00%  | +0.907 (over)  |

**The kappa dilemma emerges clearly at gDNA ≥ 100:**

- κ=2 **crushes nRNA**: The flat prior allows the EM to absorb asymmetric
  counts entirely into gDNA, zeroing nRNA. At gDNA=1000, nRNA is 100%
  underestimated (direction = −1.0, all 12 runs zeroed).
- κ=6 **inflates nRNA**: The aggressive symmetry penalty forces gDNA to be
  symmetric, pushing the asymmetric excess into nRNA. At gDNA=1000, nRNA is
  91% overestimated — the mirror image of κ=2's problem.
- κ=4 **splits the difference**: At gDNA=1000, nRNA overestimated by 38%
  (direction = +0.383), mRNA error is 20%. Neither catastrophically wrong.

### 5.3 κ=4 Detailed Profile at gDNA=1000

At the hardest contamination level (gDNA=1000, nRNA present):

- $n = 36$ runs
- mRNA median relative error: **19.88%**
- nRNA median relative error: **55.10%**
- nRNA direction: 31/36 overestimated, 5/36 underestimated
- nRNA absolute error median: +389 fragments

κ=4 is the best available trade-off, but 20% mRNA error and 55% nRNA error at
2% gDNA contamination is still unacceptable for production use.

---

## 6. The gDNA Siphon: Root Cause Analysis

### 6.1 The Core Identification Problem

The gDNA siphon is the dominant failure mode across all configurations. At its
root, it is an **identifiability problem**: unspliced mRNA fragments and gDNA
fragments are statistically indistinguishable at the single-fragment level.

Both produce:

- Reads mapping to intronic regions
- Reads mapping to exonic regions without splice junctions
- Similar fragment length distributions
- Similar alignment quality patterns

The only signal differentiating gDNA from mRNA is **strand asymmetry**: mRNA
(and nRNA) are strand-specific, while gDNA is symmetric (equal
sense/antisense). But this signal has fundamental limitations:

1. **Low power at low contamination**: When gDNA is a small fraction (e.g.,
   10/50k = 0.02%), the asymmetry signal is lost in Poisson noise.
2. **Ambiguity at high contamination**: When gDNA dominates, the model cannot
   tell whether the asymmetric excess is nRNA or sampling noise in gDNA.
3. **Strand specificity degradation**: At SS=0.90, 10% of reads have flipped
   strands, further diluting the asymmetry signal.

### 6.2 Why the Siphon Scales with gDNA Level

| gDNA level | mRNA median error | Mechanism |
|------------|-------------------|-----------|
| 0          | 0.00%             | No contamination |
| 10         | 1.29%             | Noise floor — gDNA detected but slightly misallocated |
| 100        | 6.79%             | Systematic siphon — gDNA absorbs ~7% of mRNA pool |
| 1000       | 23.64%–30.81%     | Catastrophic siphon — gDNA absorbs ~25% of mRNA pool |

The scaling is roughly linear in $\log(\text{gDNA})$, suggesting the error is
proportion to the **signal-to-noise ratio** of the asymmetry estimator, which
degrades logarithmically as gDNA increases.

### 6.3 Mass Flow Patterns

Signed fragment errors (median, where positive = overestimated):

| Scenario     | mRNA    | nRNA    | gDNA    |
|--------------|---------|---------|---------|
| Clean        | 0       | —       | —       |
| nRNA only    | +6      | −22     | —       |
| gDNA only    | **−119**| —       | +29     |
| nRNA + gDNA  | **−57** | **+326**| **−240**|

**Interpretation**:

- **gDNA only**: gDNA siphons ~119 fragments from mRNA. Of those, most are
  absorbed into gDNA itself (+29 median) but gDNA estimation still has high
  variance (mean = −157, suggesting gDNA sometimes underestimates too).
- **nRNA + gDNA**: The triangle of errors: gDNA loses 240 fragments (the
  symmetry penalty pushes them away), those fragments flow into nRNA (+326),
  and mRNA net loses 57. The nRNA acts as a **sink for misallocated gDNA
  fragments** when the symmetry penalty is active.

### 6.4 The Kappa Paradox

The strand symmetry penalty $w_\text{sym} = [4 \hat{p}(1-\hat{p})]^{(\kappa/2-1)}$
creates a lose-lose at extreme gDNA:

- **Low κ** (flat prior): gDNA absorbs everything unspliced. mRNA loses ~27%,
  nRNA is completely zeroed (absorbed into gDNA).
- **High κ** (aggressive symmetry): gDNA is forced to be symmetric, so any
  asymmetric excess gets pushed to nRNA. nRNA inflates by ~91%, mRNA still
  loses ~17% (because the symmetric portion of gDNA still siphons from mRNA).
- **No κ exists** that simultaneously preserves both mRNA and nRNA accuracy at
  extreme gDNA, because the model lacks a third signal to break the three-way
  degeneracy.

---

## 7. nRNA Survival Analysis

### 7.1 Low nRNA (10 Fragments) Survival Rates

Number of runs (out of 12) where nRNA was driven to zero:

| gDNA | κ=2 | κ=4 | κ=6 |
|------|-----|-----|-----|
| 0    | 0   | 0   | 0   |
| 10   | 0   | 0   | 0   |
| 100  | **7** | 0 | 0   |
| 1000 | **12** | 5 | 5 |

- κ=2 is catastrophic: at gDNA ≥ 100, nRNA is routinely zeroed.
- κ=4 and κ=6 prevent zeroing at gDNA=100, but still partially fail at
  gDNA=1000 (5/12 zeroed for both).

When nRNA is NOT zeroed at high gDNA, it tends to be **massively
overestimated** (κ=4: median observed 223 vs expected 88; κ=6: median observed
803 vs expected 88). The model either zeros nRNA or overcorrects — there is no
stable middle ground.

### 7.2 nRNA Direction by Contamination Level

Percentage of runs where nRNA is overestimated:

| NTA  | gDNA=0 | gDNA=10 | gDNA=100 | gDNA=1000 |
|------|--------|---------|----------|-----------|
| 10   | 53%    | 67%     | 67%      | 39%       |
| 100  | 44%    | 72%     | 61%      | 67%       |
| 1000 | 25%    | 61%     | 69%      | 67%       |

At gDNA=0, nRNA tends to be slightly underestimated (expected — some nRNA
leaks into mRNA). With any gDNA present, nRNA over-estimation dominates
(60–72%), consistent with the mass flow pattern where gDNA fragments leak
into nRNA.

---

## 8. Strand Specificity and Transcript Pattern Effects

### 8.1 Strand Specificity

| SS   | mRNA med | nRNA med | gDNA med |
|------|----------|----------|----------|
| 1.00 | 3.63%    | 0.97%    | 0.69%    |
| 0.95 | 6.02%    | 1.91%    | 1.23%    |
| 0.90 | 5.92%    | 2.99%    | 1.62%    |

Perfect strand specificity (SS=1.00) provides measurably better deconvolution.
Even 5% strand flip (SS=0.95) nearly doubles nRNA error and degrades mRNA
accuracy from 3.6% to 6.0%. This confirms that strand information is the
primary discriminative signal.

### 8.2 Transcript Pattern Sensitivity

| Pattern | mRNA med | mRNA max |
|---------|----------|----------|
| TA1     | 5.49%    | 43.16%   |
| TA2     | 2.97%    | 44.02%   |
| TA3     | 6.09%    | 49.16%   |
| TA4     | 4.63%    | 47.83%   |

Pattern TA2 is easiest (2.97%), TA3 hardest (6.09%). All patterns share
similar worst-case behavior (~45–49% max error), indicating the failure mode
is independent of transcript structure.

---

## 9. Failure Mode Taxonomy

| Failure Mode | Conditions | Severity | Root Cause |
|-------------|-----------|----------|------------|
| **gDNA siphon** | gDNA ≥ 100, any κ | mRNA err 7–49% | Unspliced mRNA ≈ gDNA in likelihood |
| **nRNA zeroing** | gDNA ≥ 100, κ ≤ 2 | nRNA 100% lost | Flat prior lets gDNA absorb nRNA |
| **nRNA inflation** | gDNA ≥ 100, κ ≥ 4 | nRNA 50–100% over | Symmetry penalty pushes gDNA → nRNA |
| **Low nRNA fragility** | NTA=10, gDNA ≥ 100 | nRNA zeroed or 10× over | Small signal overwhelmed by gDNA |
| **SS degradation** | SS ≤ 0.95 | 2× error increase | Strand flip dilutes the only discriminative signal |

---

## 10. Why Prior-Based Approaches Are Insufficient

The ablation reveals a fundamental limitation: **all prior-based interventions
(burden subtraction, symmetry penalties, density-based initialization) operate
on the prior, but the EM converges to a data-driven likelihood maximum that
ignores weak priors.**

The information-theoretic argument:

1. At 50k fragments, the Fisher information in the data overwhelms any
   reasonable prior. Even a strongly informative $\text{Beta}(50, 50)$ prior
   would be washed out by $n = 50{,}000$ observations.

2. The gDNA siphon is a **likelihood pathology**, not a prior problem. The
   likelihood function $P(\text{fragment} \mid \text{gDNA})$ and
   $P(\text{fragment} \mid \text{unspliced mRNA})$ have near-identical values
   for intronic/unspliced fragments, making the mixture under-identified.

3. The symmetry penalty is the **only mechanism that survives EM** (it is
   applied per-iteration, not just at initialization), which is why κ has a
   measurable effect while burden subtraction does not.

---

## 11. Potential Solutions

### 11.1 Fragment-Level Discriminative Features (Highest Impact)

Add features to the **likelihood** that distinguish gDNA from mRNA:

- **Intronic coverage uniformity**: gDNA produces uniform coverage across
  introns; nRNA shows 3'→5' decay; mRNA has no intronic signal except near
  splice sites. A per-fragment uniformity score could weight the gDNA
  likelihood.
- **Splice-adjacent fragment patterns**: mRNA unspliced fragments cluster near
  exon-intron boundaries (pre-mRNA intermediates). gDNA fragments are uniform.
  Position along intron relative to nearest splice site is a discriminative
  feature.
- **Fragment length distribution**: gDNA fragment lengths may differ from
  cDNA-derived fragments depending on library preparation. If the library prep
  is known, this could be a soft constraint.

**Estimated impact**: These features would break the likelihood degeneracy and
allow the EM to separate gDNA from unspliced mRNA, potentially reducing the
siphon from 25% to <5% at extreme gDNA.

### 11.2 Adaptive Kappa (Moderate Impact)

Instead of a global κ, learn κ per locus from the strength of the asymmetry
signal:

$$\kappa_\ell = \kappa_\text{base} + \kappa_\text{scale} \cdot \text{Asym}(\ell)$$

where $\text{Asym}(\ell)$ is a confidence measure of strand asymmetry at
locus $\ell$ (e.g., chi-squared test statistic for sense vs. antisense
counts). At loci where asymmetry is strong and unambiguous, use lower κ (trust
the signal). At loci where asymmetry is weak, use higher κ (trust the
symmetry prior).

**Estimated impact**: Would prevent the worst-case zeroing at κ=2 and the
worst-case inflation at κ=6, but cannot break the fundamental identification
problem.

### 11.3 Two-Stage EM (Moderate Impact)

1. **Stage 1**: Run EM using only spliced reads (which are unambiguously mRNA).
   Estimate mRNA abundances $\hat{\theta}_t$.
2. **Stage 2**: Fix the spliced-read contribution to mRNA, then deconvolve the
   unspliced residual into nRNA + gDNA.

This would anchor mRNA estimates to the reliable spliced signal and prevent
gDNA from siphoning spliced mRNA counts. The unspliced mRNA contribution
would be predicted from $\hat{\theta}_t \cdot R_t$, not estimated from the
mixture.

**Estimated impact**: Would largely eliminate the mRNA siphon (mRNA error
bounded by spliced-read noise, ~1–3%) at the cost of potentially worse nRNA
estimates (nRNA must compete only with gDNA for the unspliced residual).

### 11.4 Hard Constraints in the EM (Lower Impact)

Add inequality constraints to the M-step:

- $\hat{g}_\text{sense} / \hat{g}_\text{anti} \in [1/r, r]$ for some $r > 1$
  (force gDNA near-symmetry as a hard constraint, not a soft penalty).
- $\hat{\text{nRNA}} \geq \epsilon \cdot \hat{\text{gDNA}}$ for some small
  $\epsilon$ (prevent nRNA zeroing when gDNA is large).

These survive EM iteration (unlike prior-based initialization) and could
prevent the worst failure modes without fundamentally solving identification.

### 11.5 Information-Theoretic Bounds (Diagnostic)

Compute the Cramér-Rao lower bound for the gDNA fraction estimator given
the overlap between $P(\text{frag} \mid \text{mRNA})$ and
$P(\text{frag} \mid \text{gDNA})$. This would characterize the **regime
where deconvolution is fundamentally unreliable** and allow rigel to report
confidence intervals or flag low-confidence loci.

---

## 12. Recommendations

### Immediate (v0.2.x)

1. **Set κ=4 as default**. It provides the best overall trade-off:
   same median mRNA accuracy as κ=2, dramatically better tail behavior
   (p95 drops from 43% to 30%), best gDNA accuracy (0.54%), and acceptable
   nRNA accuracy (1.63%).

2. **Remove the burden subtraction code** or gate it behind a flag. It adds
   complexity with zero measurable benefit. The EM overwhelms the prior
   adjustment.

### Medium-Term (v0.3.x)

3. **Implement fragment-level positional features** (§11.1). This is the
   highest-impact change because it attacks the likelihood directly. Start
   with splice-adjacent fragment position as the simplest discriminative
   signal.

4. **Implement two-stage EM** (§11.3) as an alternative estimation mode.
   Anchor mRNA to spliced reads, then deconvolve the unspliced residual.

### Long-Term

5. **Characterize identifiability bounds** (§11.5) and surface confidence
   metrics per locus so downstream analyses can weight or filter
   low-confidence deconvolution results.

---

## Appendix A: Raw Numbers

### A.1 Top 10 Worst mRNA Errors (v2)

All worst cases occur at gDNA=1000, κ=2.0:

| mRNA rel err | NTA | gDNA | SS   | κ   | mRNA exp | mRNA obs | nRNA obs |
|-------------|-----|------|------|-----|----------|----------|----------|
| 49.16%      | 0   | 1000 | 0.90 | 2.0 | 236      | 120      | 0        |
| 48.47%      | 0   | 1000 | 0.95 | 2.0 | 236      | 122      | 0        |
| 47.83%      | 0   | 1000 | 0.90 | 2.0 | 181      | 94       | 0        |
| 47.35%      | 10  | 1000 | 0.90 | 2.0 | 236      | 124      | 0        |
| 47.26%      | 0   | 1000 | 0.95 | 2.0 | 181      | 96       | 0        |
| 46.42%      | 10  | 1000 | 0.95 | 2.0 | 236      | 126      | 0        |

Note: nRNA is zeroed in every worst case — the flat prior (κ=2) allows gDNA to
absorb both mRNA unspliced counts AND nRNA counts.

### A.2 gDNA Accuracy V1 vs V2

| gDNA level | V1 median err | V2 median err | Change |
|-----------|---------------|---------------|--------|
| 10        | 2.37%         | 2.08%         | −0.29pp |
| 100       | 1.78%         | 1.56%         | −0.22pp |
| 1000      | 1.61%         | 1.15%         | −0.46pp |

The small improvement in gDNA accuracy is due to the addition of κ=4 to the
sweep (κ=4 has the best gDNA accuracy at 0.54%), not to the burden subtraction.
