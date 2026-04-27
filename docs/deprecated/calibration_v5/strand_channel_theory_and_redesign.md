# Strand calibration: theoretical foundation, current implementation, failure modes, and redesign

**Date:** 2026-04-23
**Status:** Design / research. No code changes.
**Priority:** Highest — blocks the NB density work in
[density_nb_model_plan.md](density_nb_model_plan.md) and
[channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md).
**Scope:** The per-region strand likelihood in
[`src/rigel/calibration/_em.py`](../../src/rigel/calibration/_em.py)
(functions `_strand_llr`, `_strand_llr_betabinom`, `_estimate_kappa_marginal`,
`_aggregate_strand_z`) and how it is fused with the density and FL channels.

---

## 1. Clean-slate theoretical model

### 1.1 Setup

For each genomic region *i*, let

- $n_i$ = number of **unspliced** fragments overlapping region *i*.
- $k_i$ = number of those fragments whose read orientation is **sense**
  relative to the annotated transcriptional strand $s_i \in \{+1, -1\}$.
  (Spliced fragments are excluded from the strand channel: they are a
  hard RNA anchor and carry no information about whether a region is
  "not expressed".)
- $\mathrm{SS} \in (0.5, 1]$ is the library strand specificity
  estimated from uniquely-mapped spliced fragments (trained upstream,
  fixed for the calibration EM).

We want to infer, for each region with a known $s_i$, the posterior

$$\gamma_i \;=\; \Pr(z_i = G \mid k_i, n_i, \ldots)$$

where $z_i = G$ ("not expressed") means all fragments in region *i*
are gDNA-origin, and $z_i = R$ ("expressed") means there is a
non-zero RNA contribution.

### 1.2 The per-fragment sense probability

Let $\rho_i \in [0, 1]$ be the **RNA fraction** of unspliced
fragments in region *i*: $\rho_i = 0$ under $G$ by definition,
$\rho_i \in (0, 1]$ under $R$. The per-fragment probability of
observing the sense orientation is a mixture over origin:

$$p_i \;=\; \rho_i \cdot \mathrm{SS} \;+\; (1-\rho_i) \cdot \tfrac{1}{2}$$

This is the fundamental identity of the strand channel. It says:

- Pure gDNA: $\rho = 0 \Rightarrow p = 0.5$. Reads partition 50/50.
- Pure RNA: $\rho = 1 \Rightarrow p = \mathrm{SS}$. Rare-strand
  fraction equals $1 - \mathrm{SS}$.
- Mixtures interpolate **linearly** between these two endpoints.

Crucially, $p_i$ alone is a bijection to $\rho_i$ given SS:
$\rho_i = (p_i - 0.5)/(\mathrm{SS} - 0.5)$. The strand channel is
therefore a direct, closed-form estimator of $\rho_i$ — that is its
superpower.

### 1.3 Overdispersion

If fragments within a region were strand-independent Bernoulli, we
would have $k_i \sim \mathrm{Binom}(n_i, p_i)$. In practice the
Binomial is too tight:

- PCR duplicates tie reads together.
- Capture bias, mappability, and repeat-family overlap cluster
  fragments that share orientation.
- Random misannotation of $s_i$ at the region-boundary level flips
  entire runs of fragments.

The standard Binomial-with-overdispersion model is Beta-Binomial:

$$k_i \mid n_i, \pi_i \sim \mathrm{Binom}(n_i, \pi_i), \quad
\pi_i \mid p_i, \kappa_i \sim \mathrm{Beta}(\kappa_i p_i, \kappa_i (1-p_i))$$

with concentration $\kappa_i > 0$. The mean is still $p_i$; the
variance inflates by a factor $(n_i + \kappa_i)/(\kappa_i + 1)$
relative to Binomial. As $\kappa_i \to \infty$ we recover the
Binomial.

Two different overdispersions operate on the two classes:

- $\kappa_G$: dispersion around $p = 0.5$ driven by coverage
  clustering and capture bias on gDNA fragments.
- $\kappa_R$: dispersion around $p \approx \mathrm{SS}$ driven by
  antisense transcription, misannotated gene strand, and the spread
  in $\rho$ across expressed regions.

Treating them as one shared parameter (current code) is a modelling
choice, not a theoretical requirement — see §3.

### 1.4 Per-region likelihoods

Under the Beta-Binomial model:

$$P(k_i \mid n_i, z_i = G) \;=\; \mathrm{BB}(k_i; n_i, \kappa_G \cdot \tfrac{1}{2}, \kappa_G \cdot \tfrac{1}{2})$$

For the R class, because $\rho_i$ is itself unknown, we must place a
prior $f(\rho)$ on $\rho \in (0, 1]$ and marginalise:

$$P(k_i \mid n_i, z_i = R) \;=\; \int_0^1 \mathrm{BB}\!\left(k_i; n_i, \kappa_R \cdot p(\rho), \kappa_R \cdot (1-p(\rho))\right) f(\rho)\,d\rho$$

where $p(\rho) = \rho \mathrm{SS} + (1-\rho)\tfrac{1}{2}$.

The strand LLR is then

$$\mathrm{LLR}^{\text{strand}}_i \;=\; \log P(k_i \mid n_i, G) \;-\; \log P(k_i \mid n_i, R)$$

and fuses additively with density and FL LLRs to form $\mathrm{logit}(\gamma_i)$.

### 1.5 Two theoretical consequences that must not be overlooked

**(a) Strand is asymmetric in what it can conclude.**

- When observed $\hat p_i = k_i/n_i$ is close to 0.5 with
  statistically convincing $n_i$, the only $\rho$ consistent with
  the data is $\rho \approx 0$. This is a **high-confidence G
  classification** — regardless of coverage, density, or library
  type. The user's 49/51 and 4912/5123 examples fall here: both are
  overwhelming evidence for $z = G$, and sample size only tightens
  the conclusion.
- When observed $\hat p_i$ is close to SS, it is consistent with
  $\rho$ near 1. This is a high-confidence R classification, but
  only for **RNA-dominant** regions.
- When $\hat p_i$ is intermediate (say, 0.6 with $n_i = 20$), the
  region has "some" RNA but we cannot pin $\rho$. The LLR is
  modest, which is correct.

The symmetry of the Binomial / BetaBinomial LLR around these two
regimes is already in the theory — nothing magical about the boundary
conditions. What matters is that the **R-class likelihood must be
integrated over $\rho$**, otherwise it behaves as if the R class only
admits pure RNA, which is wrong.

**(b) Strand can never distinguish "pure gDNA" from "gDNA with a
trace of RNA."**

If $\rho_i$ is very small, $p_i \approx 0.5 + \rho_i(\mathrm{SS} - 0.5)$.
To separate $\rho = 0$ from $\rho = \rho_0 > 0$ with confidence at the
strand channel's per-region noise level, we need

$$n_i \gtrsim \frac{1}{\big(\rho_0 (\mathrm{SS} - 0.5)\big)^2} \cdot \text{(variance factor including }\kappa\text{)}$$

For $\mathrm{SS} = 0.9$ and $\rho_0 = 0.01$ ("1% RNA contamination"):
$n_i \gtrsim 2.5 \times 10^5$. Strand **cannot** deliver this
distinction at realistic coverage. This is a fundamental
identifiability bound — the density channel, not strand, is what
distinguishes "pure gDNA" from "mostly gDNA + trace RNA". Any strand
model that pretends otherwise is fitting noise.

The practical upshot is that strand is the **strong-G authority** and
the **strong-R authority**, and density is the **weak-R authority**
(it has to detect trace expression above the gDNA floor). The three
channels are complementary, not redundant.

---

## 2. What the current code does

Source of truth: `src/rigel/calibration/_em.py`.

### 2.1 Library-level z-test gate

[`_aggregate_strand_z`](../../src/rigel/calibration/_em.py) pools sense
and antisense counts across **all eligible regions with a known
`tx_strand`** and computes a one-sided Binomial z-score testing
$H_0: p = 0.5$ vs $H_1: p > 0.5$ at the library-pool level. If
$z \ge 3$ (default `strand_z_threshold=3.0`), the strand channel is
**enabled**; otherwise it is silently **disabled** and contributes
zero LLR everywhere.

### 2.2 Binomial mode (`strand_llr_mode="binomial"`, the default path)

[`_strand_llr`](../../src/rigel/calibration/_em.py):

$$\mathrm{LLR}^{\text{bin}}_i \;=\; k_i \log\!\tfrac{0.5}{\mathrm{SS}} + (n_i - k_i) \log\!\tfrac{0.5}{1-\mathrm{SS}}$$

applied only to regions with $n_i \ge 2$ and $s_i \ne 0$. The
combinatorial $\binom{n}{k}$ cancels. A noise floor `max(ε_bio, ε_CI)`
is applied to $1 - \mathrm{SS}$ to avoid blow-up near perfect
strandedness.

This is exactly $\log \mathrm{Binom}(k; n, 0.5) - \log \mathrm{Binom}(k; n, \mathrm{SS})$
— a two-point likelihood ratio with **no overdispersion** and **no
$\rho$ prior**.

### 2.3 Beta-Binomial mode (`strand_llr_mode="betabinom"`)

[`_strand_llr_betabinom`](../../src/rigel/calibration/_em.py):

$$\mathrm{LLR}^{\text{BB}}_i \;=\; \log \mathrm{BB}(k_i; n_i, \tfrac{\kappa}{2}, \tfrac{\kappa}{2}) - \log \mathrm{BB}(k_i; n_i, \kappa \mathrm{SS}, \kappa(1-\mathrm{SS}))$$

with a **single shared** $\kappa$ across G and R, estimated each
EM iteration by
[`_estimate_kappa_marginal`](../../src/rigel/calibration/_em.py). The
estimator maximises the 2-class mixture marginal log-likelihood via
golden-section search on $\kappa \in [0.01, 500]$.

Key property: when one component has less than one effective region
of weight ($\sum \gamma_i < 1$ or $\sum (1-\gamma_i) < 1$), the
estimator resets $\gamma \to 0.5$ and fits a **single symmetric
Beta-Binomial** — label-free.

### 2.4 Fusion

[`_e_step`](../../src/rigel/calibration/_em.py) fuses channels with
unit weights:

```python
log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft
log_odds = np.clip(log_odds, -500, 500)
gamma[soft] = 1 / (1 + exp(-log_odds))
```

No per-channel weighting, no per-region reliability weighting, no
Fisher-information scaling. The 500-clip is a numerical safety rail,
not a modelling choice.

### 2.5 What the current model effectively assumes

- R class is **pure RNA** ($p = \mathrm{SS}$). No integration over
  $\rho$.
- G and R have **identical overdispersion** ($\kappa_G = \kappa_R$).
- Overdispersion is either off (Binomial mode) or a single scalar
  shared across the genome.
- Strand is a **binary gate**: fully on (library-level
  $z \ge 3$) or fully off.
- Regions with $n_i < 2$ or ambiguous strand contribute zero LLR —
  fine in principle, but silently folds them into density's opinion.

---

## 3. Failure modes and why strand is weaker than it should be

### 3.1 Per-read information is bounded by the KL divergence

Under the Binomial LLR, the expected signal per read is

$$\mathbb{E}_{k \mid G}\!\left[\tfrac{1}{n}\mathrm{LLR}\right] = \mathrm{KL}(0.5 \,\|\, \mathrm{SS}) = \tfrac{1}{2}\log\!\tfrac{0.5}{\mathrm{SS}} + \tfrac{1}{2}\log\!\tfrac{0.5}{1-\mathrm{SS}}$$

For $\mathrm{SS} = 0.81$: 0.243 nats/read. For $\mathrm{SS} = 0.95$:
0.664 nats/read. For $\mathrm{SS} = 0.99$: 1.54 nats/read. This is
**intrinsic** to the two-hypothesis strand model.

The density channel, by contrast, has unbounded LLR per read —
$|\log(\lambda_G E / \mu)|$ grows without bound as the hypothesised
rates diverge. Even a **correctly specified** density channel
dominates unit-weight summation at moderate SS. This is not a bug;
it reflects that density genuinely carries more information per
read *if and only if the model is correct*.

### 3.2 Binomial mode ignores overdispersion → inflated strand LLR on G-class outliers

When $\kappa_G < \infty$, a G-class region with a moderate 60/40 split
(say $n = 100$, $k = 60$) is entirely plausible — the
Beta-Binomial CI is wide. But Binomial mode reports

$$\mathrm{LLR}^{\text{bin}} = 60 \log(0.5/\mathrm{SS}) + 40 \log(0.5/(1-\mathrm{SS}))$$

which at $\mathrm{SS}=0.81$ gives $-29.5 + 38.7 = +9.2$ — "class G
by a mile." Consistent with truth, but with **fragile** confidence:
add a bit more sense bias (say 65/35) and the LLR flips to $+5.0$;
at 70/30 it's $+0.8$; at 75/25 it's $-3.4$ (wrongly class R).

A correctly overdispersed model (BB with $\kappa_G$ fitted from
actual G-class regions) gives a much smoother transition. The
Binomial mode is effectively infinitely-confident-in-p=0.5 under G,
which is wrong at any realistic $\kappa_G$.

### 3.3 Beta-Binomial mode shares $\kappa$ across G and R → wrong shape under either

RNA-dominant regions tend to have higher $\kappa$ (tighter clustering
around SS — fragments within a well-expressed gene agree on
orientation). gDNA regions tend to have lower $\kappa$ (more
heterogeneity from capture and mappability artefacts).

Forcing a single $\kappa$ pins a compromise value. The practical
effect: G-class LLR is too tight (because $\kappa$ got pulled high
by RNA regions), making the channel **over-confident** on G
classifications when a region happens to observe $\hat p$ slightly
away from 0.5 just from $\kappa_G$-scale gDNA noise.

This is the mirror image of (3.2): Binomial mode is over-confident on
G because $\kappa$ is infinite; shared-$\kappa$ BB mode can still be
over-confident on G because $\kappa$ is too large for the G class.

### 3.4 R class modelled as pure RNA → LLR strongly negative on mixed regions

A region with $\rho_i = 0.3$, $\mathrm{SS} = 0.81$, $n_i = 100$,
$k_i = 0.5 \cdot 0.7 \cdot 100 + 0.81 \cdot 0.3 \cdot 100 = 59.3$
fragments sense. Under current models:

- $\log P(k | G) = \log \mathrm{BB}(59; 100, \kappa/2, \kappa/2)$
  — near the peak, modest.
- $\log P(k | R) = \log \mathrm{BB}(59; 100, \kappa \cdot 0.81, \kappa \cdot 0.19)$
  — far from the peak (expected ~81), small.

LLR is strongly positive (class G), but **the region does contain
RNA**. The density and FL channels are supposed to rescue it, but
when density is mis-specified (hybrid capture) or when the FL
channel is slow to start, the misclassification sticks.

With a $\rho$-prior-integrated R likelihood, this region's $P(k|R)$
would be dominated by the $\rho \approx 0.3$ slice of the prior,
aligning with the data and giving a moderate LLR near 0 — correctly
expressing that strand alone can't adjudicate.

### 3.5 Binary gate at library-level z = 3 is all-or-nothing

Libraries with borderline SS (e.g. stranded with high gDNA fraction
$\Rightarrow$ small pooled $z$) can fail the aggregate test despite
carrying very strong per-region strand evidence. Once disabled,
strand contributes **zero** to every region — so the channel with
the cleanest theoretical grounding silently stops contributing, even
for regions where local $n_i$ is massive.

Conversely, once enabled, all regions get strand LLR regardless of
how reliable each region's local $\hat p$ is. There is no
per-region credibility.

### 3.6 No explicit per-region G-certainty output

The user's natural question — "for a given region, how confident
are we that this is gDNA-only?" — cannot be answered by the current
code without reconstructing the posterior from the LLR sum. And
that posterior is contaminated by whatever density is saying. We
lack a **strand-only, per-region, model-principled** answer.

This is not just an interface issue. It's what makes strand "look
weak": the user sees a γ that's been blended with a bad density
LLR, and concludes "strand is too soft." In reality, strand is not
soft — it's being drowned.

### 3.7 $\kappa$ is fit on the whole mixture, biased by whichever class is mis-specified

`_estimate_kappa_marginal` maximises the mixture marginal over both
G and R components, γ-weighted. When R is mis-specified (pure-RNA
assumption), the estimator pushes $\kappa$ toward values that
accommodate the mis-specification — typically by inflating variance
(smaller $\kappa$) so that mixed-$\rho$ regions look "plausible"
under the R branch. That smaller $\kappa$ then also widens the G
branch, reducing G-class LLR magnitude.

Net effect: **mis-specification on R leaks into reduced strand
decisiveness on G** via the shared $\kappa$.

---

## 4. Proposed redesign

The goal is a strand likelihood that (a) correctly integrates over
$\rho$ under R, (b) separately parameterises $\kappa_G$ and
$\kappa_R$, (c) surfaces a per-region strand-only posterior that the
fusion layer can use as a strong prior, (d) degrades gracefully on
low-$n$ regions instead of via a library-level gate.

### 4.1 Core likelihood (class-separated Beta-Binomial + ρ prior on R)

**G class (pure gDNA, one parameter $\kappa_G$):**

$$P(k_i \mid n_i, G, \kappa_G) \;=\; \mathrm{BB}(k_i; n_i, \tfrac{\kappa_G}{2}, \tfrac{\kappa_G}{2})$$

**R class (gDNA+RNA mixture, parameters $\kappa_R$ and a $\rho$ prior):**

$$P(k_i \mid n_i, R, \kappa_R, \theta_\rho) \;=\; \int_0^1 \mathrm{BB}\!\left(k_i; n_i, \kappa_R p(\rho), \kappa_R (1-p(\rho))\right) f(\rho \mid \theta_\rho)\,d\rho$$

with $p(\rho) = \rho \mathrm{SS} + (1-\rho) \tfrac{1}{2}$.

**Choice of $f(\rho)$.** Two defensible options:

- $\rho \sim \mathrm{Beta}(a_\rho, b_\rho)$. Two parameters; fit from
  the data.
- $\rho \sim \mathrm{Uniform}(0, 1)$. Zero parameters; a deliberately
  flat prior that says "under R, RNA fraction is unknown."

I recommend **starting with $\mathrm{Uniform}(0,1)$** as the default:
it is the principled non-informative choice, it has no free
parameters to fit (avoiding another identifiability surface inside
the EM), and it correctly encodes "the R class makes no claim about
how RNA-dominant a region is." If empirical validation on the
VCaP titration shows the $\rho$ distribution is sharply peaked near
$\rho = 1$ (expressed genes are usually RNA-dominant), a 2-param
Beta is a straightforward future upgrade.

**Numerical evaluation.** The $\rho$-integral has no closed form.
A 17-node Gauss-Legendre quadrature on $[0, 1]$ is more than enough
(the integrand is smooth for $n_i \le 10^4$; for bigger $n_i$, switch
to a Laplace approximation around the mode of
$\log \mathrm{BB} + \log f$). Cost: $O(N \cdot 17)$ per EM iter,
comparable to the existing 21-node Gauss-Hermite for the density
channel.

### 4.2 Separate overdispersions, γ-weighted soft MLE (no hard subsets)

We want $\kappa_G \ne \kappa_R$ without introducing brittle
cutoffs like "regions with $\gamma > 0.95$". The standard EM
treatment already gives us this: each M-step maximises the
γ-weighted expected complete-data log-likelihood,

$$Q(\kappa_G, \kappa_R) \;=\; \sum_i \gamma_i \log P(k_i \mid n_i, G, \kappa_G) \;+\; \sum_i (1 - \gamma_i) \log P(k_i \mid n_i, R, \kappa_R)$$

The two terms decouple in their parameters — $\kappa_G$ is fit on
$\{k_i, n_i\}$ with **soft weights** $\gamma_i$, and $\kappa_R$ is fit
with weights $1 - \gamma_i$, independently. Each fit is a
1-D golden-section search (as today for the shared $\kappa$) on
the weighted log-likelihood of its own class. No cutoff, no
percentile, no seed count.

This is fully parameter-free in the same sense as the current
shared-$\kappa$ estimator: it uses whatever soft classification the
previous E-step produced. When EM is initialized neutrally
($\gamma_i = 0.5$ everywhere), both fits see identical data and
recover identical $\kappa$'s. As γ separates in subsequent
iterations, the two $\kappa$'s differentiate naturally —
self-consistently, via EM.

**Why the EM escapes a neutral init.** The two classes have
structurally different expected rates ($p = 0.5$ under G, $p \ne 0.5$
under R after integrating ρ over the upper half of $[0, 1/2, \mathrm{SS}]$).
Even at $\gamma_i = 0.5$ and identical $\kappa$'s, the LLR is
non-zero wherever the observed $k_i/n_i$ differs from 0.5, which
drives γ separation on the very first E-step.

**Hard anchors are not cutoffs — they are physical impossibilities.**
The one "hard subset" we do keep is the existing spliced-evidence
initializer: regions with $k^s > 0$ get $\gamma_i := 0$ at init.
This is not a threshold-based classification — spliced fragments
cannot arise from gDNA by construction (they cross exon-exon
junctions). It is a truth-based anchor, identical in spirit to
"regions with zero mappable bp are ineligible." We do **not** use
it as the κ_R fit subset; κ_R uses the full γ-weighted data.

**Initialization.** Start $\kappa_G = \kappa_R = \kappa_0$ where
$\kappa_0$ is one scalar (e.g. 20, or even fit once on the whole
γ=0.5 pool before EM starts). The EM carries them self-consistently
from there. There is no "subset size guard" because the soft MLE
always has the full N regions contributing — low-γ regions just
count little toward $\kappa_G$, and vice versa.

**Tiny-test-case robustness.** Because both $\kappa$'s see the full
γ-weighted data, small test cases with only a handful of regions
fit without collapsing. The old "need N > 50 anchor regions" failure
mode disappears.

### 4.3 Strand-only per-region posterior

Expose — alongside the LLR — a **strand-only posterior**:

$$\tilde\gamma^{\text{strand}}_i \;=\; \mathrm{sigmoid}\!\left(\log\tfrac{\pi_0}{1-\pi_0} + \mathrm{LLR}^{\text{strand}}_i\right)$$

where $\pi_0$ is the library-level G prior (initially 0.5). This is
a by-product of the LLR computation — essentially free — and lets
downstream code answer "strand alone says this region is 99%
gDNA" without re-running anything.

The fusion layer can then optionally use $\tilde\gamma^{\text{strand}}_i$
as a **strong prior** for density/FL, rather than summing LLRs
(the sequential-Bayes approach from
[channel_fusion_hybrid_capture.md §4.3 C.1](channel_fusion_hybrid_capture.md)).

### 4.4 Remove the binary library-level gate

The aggregate $z \ge 3$ test collapses a per-region signal onto a
library-level flag. Replace with a per-region reliability weight or
just let the LLR speak for itself:

- A region with $n_i < 5$ contributes an LLR magnitude $\le 1.2$
  nats even under perfect SS. Small automatically.
- A region with $n_i = 5000$ contributes, under the clean BB model
  with $\kappa_G = 20$, an LLR bounded by the BB KL divergence times
  $n_i$ with variance factor $\tfrac{n_i + \kappa_G}{\kappa_G + 1}$,
  which for typical $\kappa_G$ still scales linearly in $n$ at the
  observed rate. Big automatically.

The library-level z-test is trying to protect against "the library
isn't stranded" — but if SS is trained upstream to $\approx 0.5$,
the BB-under-R likelihood collapses onto the BB-under-G likelihood
and $\mathrm{LLR} \to 0$ anyway. The gate is redundant.

**Concrete replacement:** drop the `strand_used` flag. Always compute
the strand LLR. Trust the SS estimate. If SS was trained on
uniquely-mapped spliced fragments and came back at 0.52, the strand
LLR will be nearly zero everywhere, as it should be.

### 4.5 Fusion: unit-weight LLR sum, nothing else

$$\mathrm{logit}(\gamma_i) = \log\!\tfrac{\pi_0}{1-\pi_0} + \mathrm{LLR}^{\text{strand}}_i + \mathrm{LLR}^{\text{density}}_i + \mathrm{LLR}^{\text{FL}}_i$$

Under correctly-specified channels this is the optimal fusion.
The current "strand is too weak" symptom is density being wrong,
not strand being wrong.

The temptation to add an interim veto rail ("clamp $\gamma_i$ when
strand says definite-G") is real, but we resist it. Rails mask
mis-specification instead of exposing it, and they would survive
past the density fix as a permanent heuristic debt. The explicit
plan is:

1. Rebuild strand (this doc) → strand likelihoods become
   theoretically correct. Density may still be mis-specified; some
   hybrid-capture regions will still be misclassified. That is
   diagnostic information, not a bug to paper over.
2. Rebuild density (NB plan, deferred) → density likelihoods become
   theoretically correct. Fusion becomes optimal by construction.
3. No rails, no clips, no gates, no thresholds at any point.

The only "safety" code that remains is the existing numerical
$\log$-odds clip at $\pm 500$ to prevent `exp` overflow — that is
an IEEE-754 concern, not a modelling choice.

### 4.6 Summary of changes vs current code

| Element | Current | Proposed |
|---|---|---|
| G likelihood | $\mathrm{Binom}(n, 0.5)$ or $\mathrm{BB}(\kappa/2, \kappa/2)$ | $\mathrm{BB}(\kappa_G/2, \kappa_G/2)$ |
| R likelihood | $\mathrm{Binom}(n, \mathrm{SS})$ or $\mathrm{BB}(\kappa \mathrm{SS}, \kappa(1-\mathrm{SS}))$ | $\int \mathrm{BB}(\kappa_R p(\rho), \kappa_R(1-p(\rho))) f(\rho)\,d\rho$ |
| $\rho$ prior | Implicit $\rho = 1$ | $\mathrm{Uniform}(0,1)$ default; Beta upgrade later |
| Dispersion sharing | Shared $\kappa$ or none | Separate $\kappa_G$, $\kappa_R$ |
| $\kappa$ estimation | Mixture marginal MLE, γ-weighted | γ-weighted soft MLE per class (no hard subsets) |
| Channel gate | Library-level $z \ge 3$ binary | None; per-region LLR magnitudes self-attenuate |
| Per-region strand posterior | Not exposed | $\tilde\gamma^{\text{strand}}_i$ surfaced |
| Fusion | Unit-weight LLR sum | Unit-weight LLR sum, no rails/gates/clips |

---

## 5. Expected behaviour on the test libraries

`/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna{00..80}m/rigel/annotated.bam`

A note on magnitudes: the per-region strand LLR is **genuinely
bounded** by the KL divergence between the two class distributions.
Under a correctly-overdispersed BB with moderate $\kappa_G$, a
near-50/50 region's LLR grows only slowly with $n$ (roughly
$\tfrac{1}{2}\log n$ rather than the linear scaling one gets under a
mis-specified Binomial). This is a feature, not a bug: strand is
honestly reporting the information content of $\hat p$ near 0.5.

**User examples (back-of-envelope, $\mathrm{SS} = 0.81$, $\rho \sim U(0,1)$, $\kappa_G = \kappa_R = 20$):**

- Region with 49 antisense / 51 sense ($n=100$, $\hat p = 0.51$):
  $\mathrm{LLR}^{\text{strand}} \approx +1$. Weak G preference —
  correct, because $n$ is too small to rule out modest $\rho$. The
  current Binomial code gives $\approx +0.3$.

- Region with 4912 antisense / 5123 sense ($n = 10035$, $\hat p = 0.5105$):
  $\mathrm{LLR}^{\text{strand}} \approx +3$ to $+5$. Moderate G
  preference — much stronger than the small-$n$ case, bounded by the
  class overlap under Uniform($\rho$) prior. The current Binomial
  code gives $\approx +2.2$.

The user's intuition that these cases are "overwhelming gDNA"
presumes near-Binomial behaviour ($\kappa_G \to \infty$). At the
realistic overdispersion levels we'll actually fit from the data
($\kappa_G$ in the tens for hybrid capture is a defensible guess),
strand is genuinely less decisive than intuition suggests on
near-50/50 regions. **This is why we need density and FL to operate
correctly** — strand carries much of the signal but not all of it.

**The bigger win** is on mid-$\hat p$ regions that current models
mis-classify. A region with $\hat p = 0.30$ at $n = 200$ (strongly
anti-sense-biased gDNA fluctuation):
- Current Binomial: $\mathrm{LLR} \approx -30$ (pushed hard to R
  because $k/n < 0.5$ and the Binomial $p=0.5$ null has no variance
  room for $\hat p = 0.3$).
- Proposed BB($\kappa_G = 20$): $\mathrm{LLR} \approx +2$. Correctly
  recognized as a plausible gDNA fluctuation.

And on **hybrid-capture intronic hotspots** ($\hat p$ near 0.5 at
very high $n$): post-redesign the strand LLR remains modest-positive
(a few nats) rather than flipping sign under Binomial mis-fit. This
doesn't by itself overpower a mis-specified density LLR of $-50$ nats
— but it puts the strand channel on solid ground so that when
density is fixed (separately), the fusion becomes correct.

---

## 6. Plan (no code in this turn)

1. **User review** of §1 (theory) and §4.1–4.2 (the two core modelling
   changes). Specifically confirm:
   - (a) $\rho \sim \mathrm{Uniform}(0,1)$ as the default R-class
     prior, Beta upgrade deferred.
   - (b) Separate $\kappa_G$ / $\kappa_R$ with class-anchored
     estimation.
   - (c) Remove the library-level z-gate entirely.
   - (d) Add the per-region $\tilde\gamma^{\text{strand}}$ output
     for diagnostics only. **No veto rail in the fusion layer**
     (consistent with §4.5: fusion is a pure unit-weight LLR sum).

2. **If approved**, write an implementation plan (separate doc)
   covering the new kernels (`_strand_llr_bb_rho`, `_estimate_kappa_G`,
   `_estimate_kappa_R`), integration with `_e_step` / `_m_step`,
   and surfacing of $\tilde\gamma^{\text{strand}}$ in `EMFit` as a
   diagnostic field. Plus unit tests (Binomial-limit regression,
   user's 49/51 and 4912/5123 examples, pure-R class recovery,
   pure-G class recovery).

3. **Validate** on the VCaP titration:
   - Per-region strand-only γ vs current density-fused γ; disagreement
     analysis on regions where strand says "definite G" but density
     says "R".
   - Sensitivity of fit $\kappa_G$, $\kappa_R$ across gDNA spike levels.
   - Overall `quant.gdna_total` change (should move toward
     strand-inferred totals, which we already know is the right
     answer on VCaP).

4. **Only after strand redesign lands**, return to the NB density
   work per [density_nb_model_plan.md](density_nb_model_plan.md).

---

## 7. References in this repository

- [`src/rigel/calibration/_em.py`](../../src/rigel/calibration/_em.py) — current strand code (sections 2.1–2.4 above)
- [`docs/calibration/channel_fusion_hybrid_capture.md`](channel_fusion_hybrid_capture.md) — parent analysis, §4.3 C.1 sequential-Bayes fusion
- [`docs/calibration/density_nb_model_plan.md`](density_nb_model_plan.md) — density channel plan, deferred until strand is fixed
- [`docs/calibration/SESSION_HANDOFF_2026-04-23.md`](SESSION_HANDOFF_2026-04-23.md) — prior session context
