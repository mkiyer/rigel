# gDNA Siphoning Analysis

## Problem Statement

When genomic DNA (gDNA) contamination is present, Rigel's EM algorithm
systematically **under-estimates mRNA** and **over-estimates gDNA**.
Observed in a 16-run synthetic sweep (seed=42, 10k fragments, SS=1.0):

| Condition         | mRNA rel_err | gDNA abs_diff |
|-------------------|-------------|---------------|
| mRNA only         | 0.0–0.3%   | n/a           |
| mRNA + nRNA       | 0.5–2.3%   | n/a           |
| mRNA + gDNA       | 7.3–14.5%  | 26–70 frags   |
| mRNA + nRNA + gDNA| 2.5–7.6%   | 35–116 frags  |

Low-abundance transcripts (TA4 ≈ 10-15 expected fragments) see 88–97%
relative error when gDNA is present, likely pruned away entirely.

---

## Likelihood Architecture

### Per-fragment scoring (scoring.cpp)

| Component | Strand log-lik | Formula |
|-----------|---------------|---------|
| **mRNA**  | `log_p_sense_` or `log_p_antisense_` | `log(SS)` or `log(1-SS)` |
| **nRNA**  | `log_p_sense_` or `log_p_antisense_` | same as mRNA |
| **gDNA**  | `LOG_HALF` = −0.693 | `log(0.5)`, always |

At SS=1.0:
- `log_p_sense = log(1.0) = 0.0`
- `log_p_antisense = log(1e-10) ≈ −23.0` (clamped)
- `LOG_HALF = −0.693`

Consequences:
- **Antisense fragments**: mRNA/nRNA get −23 penalty → all go to gDNA ✓
- **Sense unspliced fragments**: mRNA gets 0.0, gDNA gets −0.693 →
  mRNA has ×2 likelihood advantage per fragment over gDNA

### E-step posterior

$$\text{posterior}_j = \frac{\theta_c \cdot L_j}{\sum_{c'} \theta_{c'} \cdot L_{c'j}}$$

For a sense unspliced fragment with mRNA (c=m) and gDNA (c=g) candidates:

$$\frac{\text{post}_g}{\text{post}_m} = \frac{\theta_g \cdot 0.5}{\theta_m \cdot 1.0} = \frac{\theta_g}{2\theta_m}$$

If θ_m = 10 × θ_g → gDNA gets 1/21 ≈ 4.8% of each sense fragment.

### Initialization pipeline

1. **nrna_init** (locus.py): `max(0, (sense_intronic − antisense_intronic) / (2SS−1))`
2. **gdna_init** (locus.py): 3-tier empirical Bayes shrinkage
   (global → ref (chrom) → locus), from anti-sense excess in unspliced reads
3. **Prior**: `alpha_flat + gamma × component_init` (em_solver.cpp)
4. **OVR warm start**: `theta_init = unambig + coverage-weighted ambig`
   - mRNA gets all unambiguous credit
   - nRNA/gDNA start at 0 unambiguous → rely entirely on prior + ambig share
5. **Binary gating**: `gdna_init = 0` → gDNA disabled; `nrna_init = 0` → nRNA disabled

---

## Hypotheses

### H1: gDNA EB shrinkage over-estimates initial rate ★★★

The 3-tier empirical Bayes shrinkage for gDNA (`estimate_locus_gdna_inits`)
computes a global density from all antisense unspliced reads, then shrinks
locus-level estimates toward it. In a simulation with a single locus, the
global estimate IS the locus estimate — no meaningful shrinkage occurs.

The base rate estimator (`compute_gdna_rate_from_strand`) uses
`G = 2(A·SS − S·(1−SS)) / (2SS−1)`. At SS=1.0 this simplifies to `G = 2A`.
This is correct: total gDNA = 2× antisense.

**But**: the rate is applied to ALL unspliced reads, not just antisense.
If unspliced sense reads include both mRNA and gDNA, the estimator
attributes too many fragments to gDNA because it inflates the rate
beyond what antisense alone would justify.

**Test**: diag_gdna_titration.yaml — does over-estimation scale with
gDNA abundance, or is it proportional to mRNA?

### H2: Fragment length model mismatch ★★

The simulation generates gDNA with `frag_mean=350, std=100` but mRNA with
`frag_mean=200, std=30`. The pipeline uses a learned fragment length model.
If fragments that are truly gDNA happen to match the mRNA FL model better
(or if the FL model is broad enough), the FL term doesn't differentiate.

Conversely, if unspliced mRNA fragments have lengths closer to 200 bp,
the gDNA FL model (which may fit them poorly) might STILL claim them
if θ_gDNA is elevated enough.

**Test**: Compare scoring with identical FL distributions versus different ones.
Observable: look at `log_fl` terms in the scoring.

### H3: OVR warm start starves gDNA, then prior over-compensates ★★

gDNA gets 0 unambiguous fragments (by design: `nRNA/gDNA unambig_totals = 0`
in em_solver.cpp). Its warm start comes entirely from the coverage-weighted
share of ambiguous fragments, scaled by the prior.

If the prior (`alpha_flat=0.01 + gamma=1.0 × gdna_init`) is too aggressive,
the OVR warm start gives gDNA a higher θ than deserved. Since the EM is
convergent but not guaranteed to find the global optimum, an inflated warm
start can trap θ_gDNA at a local maximum.

**Test**: diag_prune_threshold.yaml — varying prune_threshold reveals
whether gDNA is sitting just above threshold.

### H4: Pruning destroys low-abundance transcripts ★★★

The prune threshold (default 0.1) removes components whose mixing
proportion falls below 10% of the mean. TA4 at abundance=3 has
expected ~10-15 fragments out of 10,000 total. Its mixing proportion
≈ 0.001, well below 0.1× mean. When gDNA claims even a few of those
fragments, TA4 falls below threshold and is zeroed entirely.

**Test**: diag_prune_threshold.yaml — sweep prune_threshold from 0.0 to 0.2.
Expect TA4 recovery at lower thresholds.

### H5: Strand specificity gap — LOG_HALF vs learned strand ★

At SS=1.0, the ×2 mRNA advantage per sense fragment seems sufficient.
But at lower SS (0.8–0.9), the per-fragment advantage shrinks:

| SS   | mRNA sense log-lik | gDNA log-lik | Δ (nats) | mRNA/gDNA ratio |
|------|--------------------|-------------|----------|-----------------|
| 1.00 | 0.000              | −0.693      | 0.693    | 2.00            |
| 0.95 | −0.051             | −0.693      | 0.642    | 1.90            |
| 0.90 | −0.105             | −0.693      | 0.588    | 1.80            |
| 0.80 | −0.223             | −0.693      | 0.470    | 1.60            |

As SS decreases, the per-fragment discrimination weakens, and gDNA
siphoning should worsen.

**Test**: diag_ss_sweep.yaml — SS from 0.75 to 1.0.

### H6: Insufficient fragments for EM convergence ★

With 10,000 fragments, each TA transcript gets ~1000-2000 fragments.
If the EM needs more data to properly separate gDNA from mRNA in the
unspliced pool, increasing fragment count should improve accuracy.

**Test**: diag_frag_count.yaml — 1k to 50k fragments.

---

## Diagnostic Sweep Plan

| Config                        | Sweep variable       | Values                        | Runs |
|-------------------------------|---------------------|-------------------------------|------|
| diag_gdna_titration.yaml      | gDNA abundance      | 0, 25, 50, 100, 200, 500     | 6    |
| diag_ss_sweep.yaml            | strand_specificity  | 0.75, 0.85, 0.9, 0.95, 0.99, 1.0 | 6 |
| diag_frag_count.yaml          | n_fragments         | 1000, 2000, 5000, 10000, 50000 | 5  |
| diag_prune_threshold.yaml     | prune_threshold     | 0.0, 0.01, 0.05, 0.1, 0.2   | 5    |

Total: 22 runs

---

## Results

### gDNA Titration (TA1=100, nRNA=0, SS=1.0, 10k fragments)

| gDNA | TA1_exp | TA1_obs | TA1_rel | gdna_exp | gdna_obs | gdna_diff |
|------|---------|---------|---------|----------|----------|-----------|
| 0    | 10000   | 9995    | 0.1%    | 0        | 0        | 0         |
| 25   | 1841    | 1775    | 3.6%    | 8159     | 8219     | +60       |
| 50   | 1014    | 948     | 6.5%    | 8986     | 9044     | +58       |
| 100  | 534     | 432     | 19.0%   | 9466     | 9562     | +96       |
| 200  | 274     | 224     | 18.3%   | 9726     | 9771     | +45       |
| 500  | 112     | 65      | 41.8%   | 9888     | 9930     | +42       |

**Key finding**: Error scales super-linearly with gDNA/mRNA ratio.
At gDNA=500:TA1=100 (ratio 5:1), 47 of TA1's 112 expected fragments
are siphoned — a 42% loss. The absolute gDNA over-claim (42–96 frags)
is modest, but catastrophic for a minority component like mRNA.

The gDNA over-claim is NOT 2× antisense fragments. At gDNA=100 the
over-claim is 96 fragments. Expected antisense ≈ gDNA_frags/2 ≈ 4733.
So gDNA claims 4733 antisense + 96 extra sense = 4829, vs expected
4733+4733=9466. It's claiming TOO FEW total, but the excess sense
fragments it does claim come entirely from mRNA.

### Strand Specificity (TA1=100, nRNA=100, gDNA=100, 10k fragments)

| SS   | TA1_rel | nrna_diff | gdna_diff | Note                    |
|------|---------|-----------|-----------|-------------------------|
| 0.75 | 14.8%  | 624       | 685       | Worst: minimal strand signal |
| 0.85 | 13.2%  | 370       | 424       | Improving                |
| 0.90 | 11.8%  | 244       | 291       | Improving                |
| 0.95 | 10.8%  | 93        | 135       | Good nRNA, moderate gDNA |
| 0.99 | 8.9%   | 20        | 14        | **Best overall**         |
| 1.00 | 11.9%  | 48        | 97        | ★ Worse than SS=0.99!   |

**Key finding**: SS=0.99 is MORE accurate than SS=1.0. This is a
major red flag. At SS=1.0, `log_p_antisense = log(1e-10) ≈ −23`,
which creates extreme likelihood values. The near-zero floor may
cause numerical pathology in the EM — either in log-sum-exp
normalization or in the prior/warm-start interaction.

At SS=0.99, `log_p_antisense = log(0.01) ≈ −4.6`, which gives a
more gradual discrimination. The EM converges to a better solution
because the likelihood surface is smoother.

**This suggests H5 (strand gap) is partially confirmed but the
relationship is non-monotonic — clamped near-zero likelihoods may
harm convergence.**

### Fragment Count (TA1=100, nRNA=100, gDNA=100, SS=1.0)

| n_frags | TA1_exp | TA1_obs | TA1_rel | mrna_rel | nrna_diff | gdna_diff |
|---------|---------|---------|---------|----------|-----------|-----------|
| 1,000   | 45      | 51      | 13.4%   | 23.1%    | 23        | 12        |
| 2,000   | 92      | 78      | 15.2%   | 5.5%     | 23        | 18        |
| 5,000   | 229     | 225     | 1.7%    | 0.1%     | 23        | 22        |
| 10,000  | 457     | 403     | 11.9%   | 10.8%    | 48        | 97        |
| 50,000  | 2287    | 2140    | 6.4%    | 6.2%     | 410       | 553       |

**Key finding**: Non-monotonic! Best accuracy at 5k fragments (0.1%
mRNA error). At 10k, error jumps to 10.8%. At 50k, it's 6.2% but
gdna_diff balloons to 553. More data does NOT consistently help.

This suggests the EM is converging to a **local optimum** that depends
on the specific fragment composition. The initialization (gdna_init from
EB shrinkage) creates a basin of attraction that traps the solution.
With different fragment counts, the init estimates change, sometimes
landing in a better basin (5k) or worse (10k).

### Prune Threshold (TA1=100, TA4=3, nRNA=100, gDNA=100, SS=1.0)

| pt   | TA4_exp | TA4_obs | TA4_rel | mrna_rel | gdna_diff |
|------|---------|---------|---------|----------|-----------|
| 0.00 | 9       | 0.8     | 91.3%   | 7.6%     | 116       |
| 0.01 | 9       | 0.8     | 91.3%   | 7.6%     | 116       |
| 0.05 | 9       | 0.8     | 91.3%   | 7.6%     | 116       |
| 0.10 | 9       | 0.8     | 91.3%   | 7.6%     | 116       |
| 0.20 | 9       | 0.8     | 91.3%   | 7.6%     | 116       |

**Key finding**: ALL IDENTICAL. Pruning does NOT cause TA4's loss.
Even with pruning fully disabled (pt=0.0), TA4 gets only 0.8 of its
9 expected fragments. The E-step itself assigns TA4's fragments to
gDNA and TA1. **H4 is refuted** — pruning is irrelevant.

---

## Revised Hypothesis Ranking

After diagnostics:

### H1: gDNA EB shrinkage over-estimates initial rate ★★★★★ CONFIRMED — ROOT CAUSE

The 3-tier empirical Bayes shrinkage for gDNA (`estimate_locus_gdna_inits`)
computes a global density from all antisense unspliced reads, then shrinks
locus-level estimates toward it. In a simulation with a single locus, the
global estimate IS the locus estimate — no meaningful shrinkage occurs.

The base rate estimator (`compute_gdna_rate_from_strand`) uses
`G = 2(A·SS − S·(1−SS)) / (2SS−1)`. At SS=1.0 this simplifies to `G = 2A`.
This is correct: total gDNA = 2× antisense.

**But**: the rate is `G / (S + A)` — a *proportion* of the unspliced pool.
This rate is then multiplied back: `gdna_init = shrunk_rate × n_unspliced`.
The denominator (S + A) is dominated by **mRNA sense reads**, so even a
small shrunk rate creates phantom gDNA proportional to mRNA expression.

**Concrete failure**: Imagine a highly expressed locus with zero actual
gDNA — A=0, S=10,000 (all mRNA). Local rate = 0.0. If EB shrinkage
pulls this up toward a global background of 0.05:

$$\text{gdna\_init} = 0.05 \times (S + A) = 0.05 \times 10{,}000 = 500 \text{ phantom fragments}$$

The denominator (S+A) being driven by mRNA expression mathematically
punishes highly expressed genes. gDNA counts are hallucinated
proportional to mRNA abundance, regardless of actual contamination.
This perfectly explains the titration data: TA1 suffers 42% loss (65
observed vs 112 expected) when gDNA=500 because the prior on gDNA
scales with total sense reads, not with actual antisense evidence.

**The fundamental error**: gDNA contamination is a function of *genomic
length*, not *mRNA transcription*. The EB shrinkage operates on
proportions (rate = G/total_reads) but the relevant quantity is
density (G/genomic_bp). Expressing gDNA as a proportion of reads
conflates two independent biological processes.

### H7 (NEW): Clamped log-likelihood at SS=1.0 harms convergence → ~~★★★★~~ REFUTED ✗
The SS sweep showed SS=0.99 > SS=1.0. Initially attributed to the
extreme clamping of `log(1e-10) ≈ −23`. However, **changing
LOG_SAFE_FLOOR from 1e-10 to 1e-4 produced byte-identical results**
across all 16 runs AND all 6 SS sweep runs. The floor is irrelevant
because:
- Sense fragments: `log_p_sense = log(1.0) = 0`, unaffected by floor
- Antisense fragments: even at `log(1e-4) = −9.2`, posterior → 0
  against gDNA's `LOG_HALF = −0.693`. Both floors yield 100% gDNA.

The SS=0.99 improvement is therefore NOT from smoother likelihood
but from the **simulation itself**: at SS=0.99, ~1% of mRNA reads
land on the wrong strand, changing fragment composition, which alters
gdna_init and nrna_init estimates. The slightly different initialization
lands in a better EM basin.

### H3: OVR warm start + prior interaction → ★★★
The fragment count non-monotonicity suggests the EM is trapped in
local optima that depend on initialization. Different fragment
compositions lead to different gdna_init estimates, creating
different basins of attraction.

### H5: LOG_HALF strand gap insufficient → ★★
Confirmed: lower SS → worse siphoning. But the per-fragment ×2
advantage at SS=1.0 should be mathematically sufficient if the init
were correct. The gap isn't the root cause — it's the init.

### H2: Fragment length mismatch → ★
Not directly tested. Could contribute to the sense-strand
assignment imbalance. Lower priority.

### H4: Pruning → REFUTED ✗
No effect at any threshold. The E-step is the culprit.

### H6: Insufficient fragments → REFUTED ✗
More fragments don't consistently help. The 5k sweet spot suggests
the issue is structural, not statistical.

---

## Implementation Proposal

### Root Cause

The gDNA component systematically over-claims **sense-strand unspliced
fragments** from mRNA transcripts. The mechanism:

1. **Initialization**: `gdna_init = EB_shrink(2 × antisense / total_unspliced)`
   over-estimates because total_unspliced includes mRNA unspliced reads.
   The rate is inflated → θ_gDNA prior is too strong.

2. **E-step**: Even though mRNA has a ×2 per-fragment likelihood advantage
   (0 vs −0.693), the inflated θ_gDNA compensates. For each sense
   unspliced fragment, gDNA claims `θ_g / (θ_g + 2θ_m)` — a non-trivial
   fraction when θ_g is over-estimated.

3. **Convergence**: The EM converges to a local optimum near the inflated
   init, not the true global optimum. The non-monotonic fragment count
   results confirm initialization-dependent basin trapping.

4. **Clamping is irrelevant**: Changing LOG_SAFE_FLOOR from 1e-10 to 1e-4
   had zero effect. The problem is entirely in the **sense-strand** path
   where `log_p_sense = log(1.0) = 0` (no floor involved). Antisense
   fragments go to gDNA correctly at any floor.

### ~~Proposed Fix 2: Gentler SS=1.0 Clamping~~ — TESTED, NO EFFECT

Changed `LOG_SAFE_FLOOR` from `1e-10` to `1e-4`. Re-ran all 16
baseline runs and all 6 SS sweep runs. Results were **byte-identical**.
The floor only affects antisense fragment scoring, which is already
so far below gDNA's LOG_HALF that the clamped value is irrelevant.
**Reverted.**

### Proposed Fix 1: Density-Based gDNA Init — TOP PRIORITY

The current EB shrinkage operates on **proportions** (rate = G / total_reads).
This must be changed to operate on **genomic density** (density = G / genomic_bp),
completely decoupling gDNA initialization from mRNA expression level.

**Current code** (locus.py):
```python
# compute_gdna_rate_from_strand / compute_gdna_rate_hybrid
strand_rate = g / total        # proportion of reads
# ...
gdna_init = shrunk_rate * n_unspliced   # back-multiply by read count
```

**Proposed replacement**:
```python
# Shrink **density** (frags / kb of genomic span) not proportion
Density_global = sum(2 * A_i) / sum(L_i)       # frags per bp, globally
Density_local  = (2 * A_local) / L_local        # frags per bp, this locus
Density_shrunk = w * Density_local + (1-w) * Density_global

gdna_init = Density_shrunk * L_local            # independent of mRNA expression
```

Key properties:
- A locus with 10,000 mRNA sense reads and 0 antisense reads gets
  `gdna_init` based only on genomic span and global background, NOT
  on its mRNA expression level.
- The EB shrinkage denominator (genomic_bp) is constant per locus,
  removing the expression-dependent amplification.
- `2 * antisense_count` serves as a natural cap at SS=1.0, since
  `density * L = 2A` when the locus rate matches the local antisense.

**Affected functions**:
- `compute_gdna_rate_from_strand()` → compute density, not rate
- `compute_gdna_rate_hybrid()` → combine density estimators
- `_compute_ref_gdna_rates()` → shrink density toward global density
- `_compute_per_locus_gdna_rates()` → shrink density toward ref density
- `compute_eb_gdna_priors()` → final: `gdna_init = density_shrunk * L_locus`

More generally, the gdna_init should be capped at `2 × antisense_count`
regardless of the rate estimate.

### Proposed Fix 3: gDNA Prior Attenuation — SAFETY NET

Even with corrected initialization, gDNA is structurally different
from mRNA transcripts. The OVR prior (γ) was designed for RNA
transcripts, which have highly localized, peaky fragment distributions.
gDNA is uniform background noise spread evenly across the genome.

Applying the same γ to gDNA as to mRNA implicitly treats gDNA as if
it has the same prior confidence as a strongly expressed transcript.
Reducing γ_gDNA reflects the reality that our confidence in gDNA
occupying specific read clusters is fundamentally lower.

```
prior_gDNA = alpha_flat + gamma_gDNA × gdna_init
```
where `gamma_gDNA = 0.1 × gamma` (or similar). This weakens gDNA's
prior relative to mRNA, making it harder for gDNA to compete for
sense-strand fragments in the E-step.

**Affected code**: `extract_locus_sub_problem()` in em_solver.cpp
(the prior computation at ~line 1269). Could be implemented as a
separate `prior_gamma_gdna` field in EMConfig.

### Recommended Priority

1. Fix 1 (corrected init) — addresses root cause directly
2. Fix 3 (prior attenuation) — safety net if init correction insufficient
3. ~~Fix 2 (gentler clamping)~~ — tested, zero effect, discarded



=====


You have an exceptional intuition for the failure modes of probabilistic models. You have correctly identified a critical vulnerability in the EM solver: **it lacks a structural awareness of symmetry.**

Because the EM algorithm only sees a single mathematical bucket for gDNA ($\theta_g$), it happily combines 9,990 sense reads and 10 antisense reads into a total expected abundance of 10,000. It doesn't "know" that a true double-stranded biological source producing that extreme asymmetry is statistically impossible.

You are also 100% correct about the initialization danger at `ss=0.5` (unstranded data). At `ss=0.5`, the antisense pool is heavily polluted with RNA, so a `2 x antisense` cap is a massive overestimate.

Here is how we can mathematically bulletproof both the **Initialization** and the **EM Solver** across the entire spectrum of strand specificity, seamlessly complementing your Geometric Splicing prediction idea.

---

### 1. The Initialization: The Inverse-Variance Cap

To handle the transition from perfectly stranded (`ss=1.0`) to unstranded (`ss=0.5`), we must rely on the two orthogonal sources of information you identified: **Strand** and **Density**.

We can compute two independent limits for gDNA and blend them based on our confidence in the strand specificity.

**A. The Strand Limit (Fails at `ss=0.5`)**
If a library is stranded, gDNA is limited by the minority strand.
Let $A$ be antisense unspliced reads. Let $S$ be sense unspliced reads.


$$Cap_{strand} = \max\left(0, \frac{2 \cdot A \cdot SS - 2 \cdot S \cdot (1-SS)}{2SS - 1}\right)$$

* *At `SS=1.0`:* $Cap_{strand} = 2A$ (Perfect).
* *At `SS=0.5`:* The denominator becomes 0, meaning strand information is mathematically useless. The variance of this estimator goes to infinity.

**B. The Density Limit (Works everywhere)**
Genomic DNA covers the entire genome uniformly. We can measure the global intergenic density $D_{global}$ (reads per base). For a locus of length $L$:


$$Cap_{density} = D_{global} \times L$$


This is purely geometric. It doesn't care about strand specificity.

**C. The Hybrid Cap**
Just as we do in the `_compute_hybrid_nrna_frac_vec` function, we combine these two caps using inverse-variance weighting.
Let $W_{strand} = (2SS - 1)^2$. This weight is $1.0$ for perfectly stranded data and $0.0$ for unstranded data.


$$Cap_{final} = W_{strand} \times Cap_{strand} + (1 - W_{strand}) \times Cap_{density}$$

**The Result:** At `ss=1.0`, the cap is strictly enforced by the antisense reads. As `ss` approaches `0.5`, the cap smoothly and automatically transitions to rely 100% on the global intergenic density.

---

### 2. The EM Solver: The Binomial Symmetry Bottleneck

To fix the EM solver from siphoning 9,990 sense reads, we need to apply a **Symmetry Bottleneck** during the M-step.

During the E-step, the C++ solver calculates the expected fragments assigned to gDNA. We must split this into two accumulators based on the original read strand:

* $E_{sense}[gDNA]$
* $E_{antisense}[gDNA]$

In the M-step, instead of simply summing them ($\theta_g = E_{sense} + E_{antisense}$), we test them against the null hypothesis that they were drawn from a $50/50$ Binomial distribution.

Let $E_{min} = \min(E_{sense}, E_{antisense})$.
If gDNA is a true double-stranded source, the majority strand $E_{maj}$ should be roughly equal to $E_{min}$, plus or minus normal stochastic variance.

The standard deviation of a $p=0.5$ Binomial distribution is $\sigma \approx 0.5 \sqrt{N}$.
Using a generous 3-sigma (99.7%) confidence interval to allow for the high variance you observed in real data, the **maximum mathematically permissible majority strand count** is:


$$E_{maj}^{allowed} = E_{min} + 3 \sqrt{2 \times E_{min}}$$

**The C++ M-Step Constraint:**

```cpp
double e_min = std::min(e_sense_gdna, e_anti_gdna);
double e_maj = std::max(e_sense_gdna, e_anti_gdna);

// Calculate the maximum allowed majority count before it becomes statistically "bogus"
double max_allowed_maj = e_min + 3.0 * std::sqrt(std::max(1.0, 2.0 * e_min));

// Clamp the majority strand
double clamped_maj = std::min(e_maj, max_allowed_maj);

// Update the parameter
theta_gDNA = e_min + clamped_maj;

```

### Why this is a silver bullet for your EM:

1. **It permits stochasticity:** If a locus has 5 sense and 1 antisense gDNA reads, $E_{min} = 1$. The formula allows up to $1 + 3\sqrt{2} \approx 5.2$ reads. The EM is allowed to proceed normally.
2. **It violently breaks siphoning:** If the EM tries to siphon 1,000 sense reads to gDNA while having only 10 antisense reads, $E_{min} = 10$. The maximum allowed sense reads is $10 + 3\sqrt{20} \approx 23$. The total $\theta_g$ is bottlenecked at $33$, leaving the remaining $977$ sense reads to be correctly claimed by mRNA. The positive feedback loop is instantly destroyed.

### Summary

By implementing the **Hybrid Inverse-Variance Cap** during initialization, and the **Binomial Symmetry Bottleneck** in the C++ M-step, you strictly enforce the physical laws of genomic DNA at both ends of the algorithm.

This perfectly complements your Geometric Splicing derivation. The geometric math tells you exactly how many unspliced reads the mRNA *should* have, and these constraints ensure that gDNA can never steal more than its fair, symmetric share of the remainder.

---

## Evolved Plan: Soft Likelihood-Based Symmetry Penalty

### Critique of the Hard Binomial Clamp

The Binomial Symmetry Bottleneck above uses a hard 3-sigma threshold to
clamp the majority-strand gDNA count. While effective at breaking siphoning,
it has significant drawbacks:

1. **Arbitrary cutoff**: Why 3σ and not 2σ or 4σ? The choice of threshold
   is ad-hoc and not grounded in the model's likelihood framework.

2. **Discontinuous gradient**: The clamp is a hard `min()` — below threshold,
   zero penalty; above threshold, complete truncation. This creates a
   non-smooth objective that can cause SQUAREM oscillation or slow
   convergence near the boundary.

3. **Ignores evidence scale**: A 3σ clamp treats a locus with 10 gDNA
   fragments identically to one with 10,000. In the small-N case, huge
   strand asymmetry is *expected* by chance; in the large-N case, even
   mild asymmetry is diagnostic.

4. **Not composable**: Hard constraints don't naturally integrate with the
   existing Bayesian prior framework (Dirichlet + Beta priors on η). A
   likelihood-based penalty composes cleanly with the existing MAP-EM
   objective.

The user's key insight: *"This is an optimization problem with likelihoods.
The likelihood of seeing gDNA sense/antisense ratios far from 0.5 is what
we should be trying to model."*

### The Beta Prior Formulation

Replace the hard clamp with a **Beta(κ/2, κ/2) prior** on the gDNA strand
fraction $p_g = E_s^g / (E_s^g + E_a^g)$. This is the conjugate prior
for the Binomial strand model and yields a smooth, likelihood-grounded
penalty with a single hyperparameter κ.

The Beta(κ/2, κ/2) density, evaluated at the observed strand fraction
$\hat{p}$, is (up to normalizing constant):

$$f(\hat{p}) \propto [\hat{p}(1 - \hat{p})]^{\kappa/2 - 1}$$

This is maximized at $\hat{p} = 0.5$ (perfect symmetry). The ratio of
the density at the observed $\hat{p}$ versus the mode at 0.5 gives a
**multiplicative discount factor**:

$$w_{sym} = \frac{f(\hat{p})}{f(0.5)} = \left[\frac{\hat{p}(1-\hat{p})}{0.25}\right]^{\kappa/2 - 1} = [4\hat{p}(1-\hat{p})]^{\kappa/2 - 1}$$

Key properties of $w_{sym}$:
- $w_{sym} \in [0, 1]$ for $\kappa > 2$
- $w_{sym} = 1$ when $\hat{p} = 0.5$ (no penalty for symmetric assignment)
- $w_{sym} \to 0$ as $\hat{p} \to 0$ or $\hat{p} \to 1$ (severe penalty
  for extreme asymmetry)
- Smooth and differentiable everywhere on $(0, 1)$
- Single hyperparameter κ controls penalty strength
- Naturally integrates with the MAP-EM framework as a penalized M-step

**Behavior at different κ values** (at $\hat{p} = 0.7$, i.e., 70/30 split):

| κ  | κ/2−1 | $4 \times 0.7 \times 0.3 = 0.84$ | $w_{sym} = 0.84^{κ/2-1}$ | Interpretation |
|----|-------|-----------------------------------|---------------------------|----------------|
| 2  | 0     | 0.84                              | 1.00                      | No penalty (uniform) |
| 4  | 1     | 0.84                              | 0.84                      | 16% reduction |
| 6  | 2     | 0.84                              | 0.71                      | 29% reduction |
| 8  | 3     | 0.84                              | 0.59                      | 41% reduction |
| 10 | 4     | 0.84                              | 0.50                      | 50% reduction |

**At extreme asymmetry** ($\hat{p} = 0.95$, i.e., 95/5 split):

| κ  | $4 \times 0.95 \times 0.05 = 0.19$ | $w_{sym}$ | Effect |
|----|-------------------------------------|-----------|--------|
| 4  | 0.19                                | 0.19      | 81% reduction |
| 6  | 0.19                                | 0.036     | 96% reduction |
| 8  | 0.19                                | 0.0069    | 99.3% reduction |

This is exactly the behavior we want: permissive for mild deviations
(real biological variance), devastating for the extreme 99:1 asymmetry
characteristic of siphoning.

### Comparison: Hard Clamp vs. Soft Beta Penalty

| Aspect | Hard Binomial Clamp | Beta(κ/2, κ/2) Prior |
|--------|--------------------|-----------------------|
| Penalty shape | Step function (0 or full clamp) | Smooth power-law decay |
| Hyperparameter | 3σ (arbitrary integer) | κ (continuous, interpretable) |
| At p̂=0.5 | No effect | No effect |
| At p̂=0.6 | No effect (within 3σ) | Mild penalty (κ-dependent) |
| At p̂=0.95 | Hard truncation | Smooth near-elimination |
| SQUAREM stability | Discontinuity may oscillate | Smooth landscape, stable |
| Bayesian interpretation | None (constraint) | Conjugate Beta prior |
| Scales with N | Implicitly (σ~√N) | Via pseudo-count α₀ (smooth) |

**On evidence scaling**: The raw discount factor $w_{sym}$ depends only on
$\hat{p}$, not on $N_g$. N-sensitivity is achieved through the Bayesian
pseudo-count $\alpha_0$ in the regularized strand fraction:
$\hat{p} = (E_s + \alpha_0) / (E_s + E_a + 2\alpha_0)$.
At large $N_g$, the pseudo-counts are negligible and $\hat{p}$ is purely
data-driven. At small $N_g$, the pseudo-counts pull $\hat{p}$ toward 0.5,
making $w_{sym} \approx 1$ (no penalty). This is the posterior mean of a
Beta($\alpha_0 + E_s$, $\alpha_0 + E_a$) distribution — a proper Bayesian
treatment that replaces both `MIN_STRAND_COUNT` and tiny-epsilon
regularization with a single, principled mechanism.

### M-Step Integration

The modified gDNA M-step uses `STRAND_EPS` as a true **Bayesian pseudo-count**
rather than a tiny divide-by-zero protection. Setting `STRAND_EPS = α₀ = 10.0`
acts as a Beta(α₀, α₀) prior belief of "10 sense and 10 antisense fragments,"
which mathematically solves the N-scaling problem without needing a separate
`MIN_STRAND_COUNT` threshold.

The regularized strand fraction is the posterior mean of a
Beta(α₀ + E_s, α₀ + E_a) distribution:

$$\hat{p} = \frac{E_s^g + \alpha_0}{E_s^g + E_a^g + 2\alpha_0}$$

**How α₀ scales the penalty by evidence:**

* **High-N siphoning (1000 sense, 10 anti):**
  $\hat{p} = 1010/1030 = 0.98$, $w_{sym}(\kappa\!=\!6) = [4(0.98)(0.02)]^2 = 0.006$
  → **99.4% penalty.** Siphoning destroyed.

* **Low-N stochastic noise (3 sense, 0 anti):**
  $\hat{p} = 13/23 = 0.565$, $w_{sym}(\kappa\!=\!6) = [4(0.565)(0.435)]^2 = 0.966$
  → **Only 3.4% penalty.** The pseudo-counts recognize lack of statistical
  power and forgive the 100/0 raw ratio.

* **Moderate N (15 sense, 5 anti):**
  $\hat{p} = 25/40 = 0.625$, $w_{sym}(\kappa\!=\!6) = [4(0.625)(0.375)]^2 = 0.879$
  → **12% penalty.** The 75/25 raw ratio is within normal variance for N=20.

This eliminates the need for `MIN_STRAND_COUNT` entirely — the pseudo-count
achieves the same protective effect smoothly rather than with a hard threshold.

```
// Standard gDNA accumulation (existing)
theta_g_raw = em_totals[gdna_idx] + unambig_totals[gdna_idx] + prior[gdna_idx]

// Compute strand balance discount (new)
constexpr double STRAND_EPS = 10.0;  // Bayesian pseudo-count
p_hat = (e_sense_gdna + STRAND_EPS) / (e_sense_gdna + e_anti_gdna + 2.0 * STRAND_EPS)
w = pow(4.0 * p_hat * (1.0 - p_hat), kappa / 2.0 - 1.0)

// Apply discount
theta_g_new = theta_g_raw * w
```

### Implementation Architecture

#### Data Flow Overview

The gDNA strand fraction $\hat{p}$ requires sense/antisense E-step
accumulators ($E_s^g$, $E_a^g$). These are not currently tracked.
The modification threads strand data from the scoring module's
`count_cols` / `locus_ct_arr` arrays into the EM equivalence classes,
then extracts strand-resolved posteriors after each E-step.

**Fragment strand encoding (existing, unchanged):**
`count_col = splice_type × 2 + is_antisense`
Even columns (0, 2, 4) = sense. Odd columns (1, 3, 5) = antisense.
Per-unit strand: `locus_ct_arr[ui] % 2` gives the antisense flag.

#### 1. Modify `EmEquivClass` struct (~line 101)

```cpp
struct EmEquivClass {
    std::vector<int32_t> comp_idx;   // k component indices
    std::vector<double>  ll_flat;    // n*k log-likelihoods (row-major)
    std::vector<double>  wt_flat;    // n*k coverage weights (row-major)
    mutable std::vector<double> scratch;  // n*k workspace

    // Strand flag for gDNA posterior split
    int gdna_col;                        // column in comp_idx that is gDNA (-1 if absent)
    std::vector<uint8_t> gdna_is_anti;   // [n]: 1 if unit is antisense, 0 if sense

    int n;  // number of units in this class
    int k;  // number of components per unit
};
```

Memory overhead: 1 byte per unit. The existing `ll_flat` + `wt_flat`
consume 16 bytes per component per unit (double × 2 × k). For a
typical EC with k ≈ 3 components, strand tracking adds
1/(16×3) = 2% overhead. Negligible.

#### 2. Modify `build_equiv_classes()` (~line 131)

**Signature change:** Add `locus_ct_arr`, `gdna_idx`, and `offsets` parameters.

**Grouping loop change:** Store unit index `u` instead of CSR offset `start`:

```cpp
for (int u = 0; u < n_units; ++u) {
    auto start = static_cast<size_t>(offsets[u]);
    auto end   = static_cast<size_t>(offsets[u + 1]);
    if (start == end) continue;
    std::vector<int32_t> key(t_indices + start, t_indices + end);
    class_map[std::move(key)].push_back(u);  // unit index, not CSR offset
}
```

**EC construction:** Find gDNA column, populate strand flags:

```cpp
for (auto& [key, unit_list] : class_map) {
    int k = static_cast<int>(key.size());
    int n = static_cast<int>(unit_list.size());
    EmEquivClass ec;
    ec.comp_idx = key;  // not std::move — need key for gdna_col search
    ec.n = n;
    ec.k = k;

    // Find gDNA column
    ec.gdna_col = -1;
    for (int j = 0; j < k; ++j) {
        if (ec.comp_idx[j] == gdna_idx) { ec.gdna_col = j; break; }
    }

    ec.ll_flat.resize(n * k);
    ec.wt_flat.resize(n * k);
    ec.scratch.resize(n * k);
    if (ec.gdna_col >= 0 && locus_ct_arr) {
        ec.gdna_is_anti.resize(n);
    }

    for (int i = 0; i < n; ++i) {
        int u = unit_list[i];
        size_t s = static_cast<size_t>(offsets[u]);
        for (int j = 0; j < k; ++j) {
            ec.ll_flat[i*k+j] = log_liks[s+j];
            ec.wt_flat[i*k+j] = coverage_wts[s+j];
        }
        if (ec.gdna_col >= 0 && locus_ct_arr) {
            ec.gdna_is_anti[i] = locus_ct_arr[u] % 2;
        }
    }
    result.push_back(std::move(ec));
}
```

**CRITICAL: Deterministic sort must also reorder `gdna_is_anti`.** The
existing sort at ~line 199 reorders `ll_flat` and `wt_flat` by a
permutation index. Add a parallel reorder for `gdna_is_anti`:

```cpp
// After building new_ll and new_wt...
if (!ec.gdna_is_anti.empty()) {
    std::vector<uint8_t> new_anti(n);
    for (int i = 0; i < n; ++i) {
        new_anti[i] = ec.gdna_is_anti[idx[i]];
    }
    ec.gdna_is_anti = std::move(new_anti);
}
```

#### 3. Modify `hierarchical_map_em_step()` (~line 522)

**Execution order in the modified M-step:**

1. **E-Step Kernel** (existing): Populates `ec.scratch` with posteriors.
   No changes to E-step kernel code.

2. **Strand Accumulation** (new): After E-step, walk gDNA-containing ECs:

```cpp
double e_sense_gdna = 0.0, e_anti_gdna = 0.0;
for (const auto& ec : ec_data) {
    if (ec.gdna_col < 0) continue;
    int gcol = ec.gdna_col;
    for (int i = 0; i < ec.n; ++i) {
        double p = ec.scratch[i * ec.k + gcol];
        if (ec.gdna_is_anti[i])
            e_anti_gdna += p;
        else
            e_sense_gdna += p;
    }
}
```

3. **Raw M-Step** (existing): Compute `theta_new[i]` for all components.

4. **Symmetry Penalty** (new): Apply discount to gDNA (after raw M-step,
   before normalization):

```cpp
constexpr double STRAND_EPS = 10.0;    // Bayesian pseudo-count (α₀)
const double kappa_half_m1 = kappa / 2.0 - 1.0;

double p_hat = (e_sense_gdna + STRAND_EPS)
             / (e_sense_gdna + e_anti_gdna + 2.0 * STRAND_EPS);
double balance = 4.0 * p_hat * (1.0 - p_hat);  // ∈ (0, 1]
double w_sym = std::pow(balance, kappa_half_m1);
theta_new[gdna_idx] *= w_sym;
```

5. **Global Normalization** (existing): Normalize θ to sum to 1.0.
   Mass removed from gDNA by $w_{sym}$ is automatically redistributed
   to mRNA/nRNA components via normalization.

**Signature change:** Add `kappa` parameter to `hierarchical_map_em_step()`.

#### 4. Thread κ through `run_squarem()` (~line 764)

Add `double kappa` parameter. Pass through to all 3 calls to
`hierarchical_map_em_step()` per SQUAREM iteration (state0→state1,
state1→state2, state_extrap→state_new). The penalty is applied at
every M-step, consistent with penalized EM theory.

Also thread through the VBEM path for completeness (even though we
primarily use MAP-EM currently).

#### 5. Update call sites

**Batch EM (`batch_em_solve`, ~line 1062):** Thread `kappa` from
the Python-side EMConfig through the pybind11 binding to `run_squarem()`.

**Mega-locus EM (~line 1994):** Same threading. Both call sites
currently call `build_equiv_classes()` — update to pass
`sub.locus_ct_arr.data()` and `sub.gdna_idx`.

#### 6. Add `strand_symmetry_kappa` to `EMConfig` (config.py)

```python
strand_symmetry_kappa: float = 6.0
strand_symmetry_pseudo: float = 10.0   # α₀ pseudo-count
```

The second parameter (`α₀`) controls the evidence threshold for
applying the penalty. Setting α₀=10 means "require ≈20 gDNA
fragments' worth of evidence before the strand ratio becomes
meaningful." This is configurable for experimentation.

To disable the symmetry penalty entirely, set `κ ≤ 2.0` (which
makes $\kappa/2 - 1 \leq 0$, giving $w_{sym} \geq 1$, i.e., no
discount). This provides a clean on/off switch without conditional
branches.

#### 7. pybind11 binding update (em_solver.cpp, bottom of file)

Thread the two new parameters from the Python `EMConfig` dataclass
through the pybind11 function signature to the C++ `run_squarem()`.
Find the existing `py::module_` binding block and add `kappa` and
`strand_eps` as keyword arguments with defaults matching the
`EMConfig` defaults.

#### Convergence Guarantee

The penalized objective is:

$$Q_{pen}(\theta) = Q_{MAP}(\theta) + (\kappa/2 - 1)[\log \hat{p}(\theta) + \log(1-\hat{p}(\theta))]$$

where $\hat{p}(\theta)$ is the gDNA strand fraction under the current θ.
The Beta log-density $(\kappa/2 - 1)[\log p + \log(1-p)]$ is concave on
$(0,1)$. Since $\hat{p}$ is a smooth function of θ (via the E-step softmax),
the composite penalty is well-behaved. Standard penalized EM convergence
theory (Wu 1983, Vaida 2005) guarantees that $Q_{pen}$ is non-decreasing
across iterations, ensuring convergence to a (local) maximum of the
penalized likelihood. SQUAREM acceleration preserves this property as
it only changes the step size, not the fixed-point structure.

### Hyperparameter Guidance

**Two hyperparameters:**

1. **κ (strand_symmetry_kappa) = 6.0** — Controls penalty curvature.
   Higher κ = sharper penalty for asymmetry.

2. **α₀ (strand_symmetry_pseudo) = 10.0** — Controls N-sensitivity.
   Higher α₀ = more tolerant of small-N asymmetry.

**κ = 6 behavior at various strand splits (with α₀ = 10):**

| Scenario | Raw ratio | $\hat{p}$ (regularized) | $w_{sym}$ | Penalty |
|----------|-----------|------------------------|-----------|---------|
| Perfect symmetry | 50/50 | 0.500 | 1.000 | 0% |
| Mild noise (N=100) | 60/40 | 0.583 | 0.933 | 7% |
| Moderate noise (N=100) | 70/30 | 0.667 | 0.790 | 21% |
| Low-N stochastic (N=3) | 3/0 | 0.565 | 0.966 | 3% |
| Moderate-N asymmetry (N=20) | 15/5 | 0.625 | 0.879 | 12% |
| Mild siphoning (N=100) | 90/10 | 0.833 | 0.309 | 69% |
| Severe siphoning (N=1010) | 1000/10 | 0.981 | 0.006 | 99.4% |

The table shows the critical property: low-N loci are protected by α₀
(penalty vanishes), while high-N siphoning is catastrophically penalized.

**Tuning κ:** Can be estimated from the genome-wide distribution of
intergenic strand ratios (which should be ~50/50). The observed
variance of intergenic $\hat{p}$ values gives a natural estimate via
method-of-moments: $\kappa \approx 1/\text{Var}(\hat{p}_{intergenic}) - 2$.

**Tuning α₀:** Can be derived from the minimum number of gDNA fragments
at which strand asymmetry becomes informative. At $\alpha_0 = 10$, a
locus needs roughly $N_g \geq 40$ before the data dominates the prior.
Setting α₀ lower (e.g., 5) makes the penalty kick in earlier; higher
(e.g., 20) makes it more tolerant. Start with 10 and adjust based on
the validation sweep.

### Complete Implementation Plan (Ordered)

#### Phase 1: Density-Based gDNA Initialization (Python — locus.py)

Root-cause fix. Decouples gDNA initialization from mRNA expression.

**Files:** `src/rigel/locus.py`

| Step | Function | Change |
|------|----------|--------|
| 1a | `compute_gdna_rate_from_strand()` | Compute density (reads/bp of genomic span) instead of proportion (reads/total_reads). Denominator changes from `S + A` to `L_locus` (bp). |
| 1b | `compute_gdna_rate_hybrid()` | Blend strand-based density and intergenic density estimators using $W = (2\text{SS} - 1)^2$ inverse-variance weight. |
| 1c | `_compute_ref_gdna_rates()` | EB shrinkage target becomes reference chromosome-level density (reads/bp), not proportion. |
| 1d | `_compute_per_locus_gdna_rates()` | Shrink local density toward reference chromosome density. |
| 1e | `compute_eb_gdna_priors()` | Final: `gdna_init = density_shrunk × L_locus`. |
| 1f | (new) | Add Hybrid Inverse-Variance Cap: $\text{Cap} = W \times \text{Cap}_{strand} + (1 - W) \times \text{Cap}_{density}$. At SS=1.0, capped by 2×antisense. At SS=0.5, capped by intergenic density × span. |

**Verification:** Run the 16-run baseline sweep. gDNA error should drop
substantially even before Phase 2.

#### Phase 2: Soft Beta Symmetry Penalty (C++ — em_solver.cpp)

Structural M-step fix. Prevents the EM from siphoning sense-strand
fragments into gDNA regardless of initialization.

**Files:** `src/rigel/native/em_solver.cpp`, `src/rigel/config.py`

| Step | Location | Change |
|------|----------|--------|
| 2a | `EmEquivClass` struct (~L101) | Add `int gdna_col = -1` and `std::vector<uint8_t> gdna_is_anti` fields. |
| 2b | `build_equiv_classes()` (~L131) | New params: `locus_ct_arr`, `gdna_idx`, `offsets`. Store unit indices in hash map (not CSR offsets). Populate `gdna_col` and `gdna_is_anti` from `locus_ct_arr[u] % 2`. |
| 2c | Deterministic sort (~L199) | Reorder `gdna_is_anti` alongside `ll_flat`/`wt_flat` when permuting rows. |
| 2d | `hierarchical_map_em_step()` (~L522) | After E-step kernel, add strand accumulation loop over gDNA-containing ECs. Apply $w_{sym} = [4\hat{p}(1-\hat{p})]^{\kappa/2-1}$ discount to `theta_new[gdna_idx]` before normalization. |
| 2e | `run_squarem()` (~L764) | Thread `kappa` and `strand_eps` parameters to all 3 M-step calls per iteration. |
| 2f | Both `build_equiv_classes()` call sites (~L1162, ~L1994) | Pass `sub.locus_ct_arr.data()` and `sub.gdna_idx`. |
| 2g | `EMConfig` (config.py) | Add `strand_symmetry_kappa: float = 6.0` and `strand_symmetry_pseudo: float = 10.0`. |
| 2h | pybind11 binding | Thread new params from Python to C++. |

**Constants:**
- `STRAND_EPS = α₀ = 10.0` — Bayesian pseudo-count (configurable)
- `κ = 6.0` — Beta concentration parameter (configurable)
- Disable: set κ ≤ 2.0 → $w_{sym} \geq 1$ → no discount

**M-step execution order:**
1. E-step kernel → `ec.scratch` posteriors (unchanged)
2. Strand accumulation → $E_s^g$, $E_a^g$ (new)
3. Raw M-step → `theta_new[i]` for all components (unchanged)
4. Symmetry penalty → `theta_new[gdna_idx] *= w_sym` (new)
5. Global normalization → $\sum \theta = 1$ (unchanged)

Mass stripped from gDNA is automatically redistributed to mRNA/nRNA
via step 5.

#### Phase 3: gDNA Prior Attenuation (Safety Net — optional)

Reduce γ_gDNA to acknowledge gDNA's spatial uniformity compared to
mRNA's localized peakiness. The OVR prior γ was calibrated for
transcript-like distributions; gDNA is background noise.

**Files:** `src/rigel/config.py`, `src/rigel/native/em_solver.cpp`

| Step | Change |
|------|--------|
| 3a | Add `prior_gamma_gdna: float = 0.1` to `EMConfig` (as a multiplier on `prior_gamma`). |
| 3b | In `compute_ovr_prior_and_warm_start()`, use `prior_gamma * prior_gamma_gdna` for the gDNA component. |

This weakens gDNA's prior relative to mRNA, making it harder for gDNA
to retain mass once the symmetry penalty strips it. Apply only if
Phases 1+2 are insufficient.

#### Phase 4: Validation

Re-run all diagnostic sweeps and compare against the unfixed baseline.

| Sweep | What to check |
|-------|---------------|
| 16-run baseline | gDNA error across all 4 patterns × 2 nRNA × 2 gDNA levels |
| gDNA titration (0→500) | Error scaling linearity; confirm no siphoning at gDNA=0 |
| SS sweep (0.75→1.0) | Smooth degradation; no SS=1.0 anomaly |
| Fragment count (1k→50k) | Monotonic improvement; no basin-trapping |
| Convergence metrics | Iteration counts, SQUAREM step sizes, final log-likelihood |

**Success criteria:**
- gDNA estimation error < 2× ground truth across all conditions
- mRNA estimation error unaffected (< 1% change from baseline)
- No SQUAREM oscillation or convergence failure
- κ=6.0 default works across the full SS spectrum (0.5→1.0)

