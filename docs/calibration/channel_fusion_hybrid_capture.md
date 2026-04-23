# Channel fusion in calibration: problem statement, current failure, and theoretically-grounded alternatives

**Date:** 2026-04-22
**Status:** Research / design proposal. No code changes. Supersedes the "per-kind stratification" idea in `strand_first_hybrid_capture_plan.md`.
**Scope:** Rethink how the three calibration channels (density / strand / gDNA-FL) are fused. Do not assume uniform λ_G. Do not use arbitrary thresholds.

---

## 1. Problem statement

Rigel's calibration EM (`_em.py::run_em`) infers, for each genomic region *i*, a posterior probability γ_i ∈ [0,1] that the region's unspliced fragments are gDNA (class G) rather than expressed RNA (class R = gDNA+RNA mixture). The posterior is fused from three conditionally-independent channels, each computed as a per-region log-likelihood ratio (LLR) under the two classes:

1. **Density (count) channel.** `ℓ_count(k_i, E_i | λ_G, μ_R, σ_R²)`. Assumes a single global per-bp gDNA rate λ_G, and that expressed regions have a LogNormal per-bp excess rate μ ~ LogNormal(μ_R, σ_R²).
2. **Strand channel.** `ℓ_strand(k_sense_i, n_i | ss, κ)` for oriented regions. Under class G, sense ≈ antisense (Binomial / BetaBinomial with mean 0.5). Under class R, sense ≈ SS (≈0.81 in the VCaP library, ≈1.0 in dUTP libraries with no antisense contamination).
3. **gDNA FL channel.** `ℓ_fl(FL observations in region i | f_gDNA(x), f_RNA(x))`. Compares the per-region FL histogram against the γ-weighted gDNA FL distribution vs a (1−γ)-weighted RNA FL distribution, both re-estimated each iteration.

The E-step currently combines them additively with **unit weights**:

```python
# src/rigel/calibration/_em.py, _e_step
log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft
gamma = 1 / (1 + exp(-log_odds))
```

And the M-step updates λ_G via a γ-weighted Poisson MLE, μ_R/σ_R from "hard-spliced anchor regions" only (`k^s > 0 & eligible`), and π/π_soft from γ-mean over eligible soft regions.

### 1.1 The hybrid-capture pathology

In hybrid-capture RNA-seq (VCaP titration, and many clinical panels):

- **λ_G is not spatially uniform.** The per-bp gDNA rate varies 100–140× between on-target and off-target regions (e.g., `λ_G_exonic / λ_G_intergenic ≈ 139` at dna80m).
- **The density channel's two-component Poisson-LogNormal mixture is mis-specified.** It assumes one λ_G shared across all regions. Forced to produce a single λ_G, the EM settles on a compromise value (roughly the geometric mean of the on/off-target rates). That value is too high for intergenic regions (where `λ_G_cal · E_i >> observed k_i` → clip-dominated) and too low for intronic regions (where `λ_G_cal · E_i << observed k_i` → almost all evidence goes to class R).
- **The strand channel tells a clean, simple story.** With SS=0.81, the rare-strand fraction has signature means 0.19 (RNA) vs 0.5 (gDNA). The separation is 0.31 — roughly 6σ per 100-read region. Across the VCaP titration, strand-inferred gDNA tracks quant.gdna_total to within 1.1–1.5×, with the strand channel correctly identifying intronic regions as 94% gDNA at dna10m.
- **The FL channel is invariant to strand specificity but has slow start-up.** It is disabled in iter 0 (no gDNA FL model yet), and it contributes information only when gDNA and RNA FL distributions differ noticeably (typical: gDNA FL has wider tails and a different mode).

The density channel is **confidently wrong** inside expressed loci under hybrid capture: it pushes γ_i → 0 (class R) because `k_i / E_i` looks "too high to be gDNA" at the uniform λ_G. But the strand channel correctly reports γ_i → 1 because the rare fraction is ≈0.5. They fight, and addition with unit weights gives the density channel a numerical advantage (it can produce LLRs of ±100 per region with modest read counts, while the strand channel's LLR is capped by the finite Kullback–Leibler divergence between `Binom(n, 0.5)` and `Binom(n, ss)`).

### 1.2 Why unit-weight LLR summation is not the right answer

Under the **ideal generative model**, summing LLRs is optimal — LLR additivity is what "conditional independence" means, and no reweighting is needed. The failure is not a weighting issue; it is a **model-specification issue**: the density channel is asserting a likelihood derived from a model that the data systematically violates. Summing a high-confidence wrong likelihood with a low-confidence right likelihood gives high-confidence wrong posteriors.

The right conceptual fix is one of:

(a) **Fix the density model so its likelihoods are well-calibrated** (model λ_G non-uniformity explicitly).
(b) **Provide targets to the density model** (targets-BED), restricting density inference to on-target regions where uniformity is roughly valid.
(c) **Down-weight the density channel when it is unreliable** — but based on a data-reliability signal, not an arbitrary threshold. This is equivalent to making the variance of the density likelihood proportional to its mis-specification, i.e., widening its effective σ.
(d) **Replace the density channel with a richer model that accommodates local rate variation** (hierarchical λ_G; neighborhood smoothing; spatially-adaptive background).

---

## 2. Current implementation (what's actually in the code)

### 2.1 The three channels

**Density LLR** [`_count_llr_poisson_ln`](src/rigel/calibration/_em.py#L134):

$$\mathrm{LLR}_{\text{count},i} = \log \text{Poisson}(k_i;\;\lambda_G E_i) - \log \int_0^\infty \text{Poisson}(k_i;\;(\lambda_G+\mu)E_i)\,\text{LogNormal}(\mu;\;\mu_R,\sigma_R^2)\,d\mu$$

with the integral approximated by 21-node Gauss-HermiteE quadrature. This is applied to **every soft eligible region** unconditionally.

**Strand LLR** ([`_strand_llr`](src/rigel/calibration/_em.py#L210) / [`_strand_llr_betabinom`](src/rigel/calibration/_em.py#L295)):

$$\mathrm{LLR}_{\text{strand},i} = \log \text{Binom}(k_{\text{sense},i};\;n_i, 0.5) - \log \text{Binom}(k_{\text{sense},i};\;n_i, \text{ss})$$

Gated by an aggregate Wald z-test on the pooled sense fraction: **disabled entirely** unless the whole-library z ≥ 3.0 (default), and for each region is zero unless `n_unspliced ≥ 2` and `tx_strand ≠ 0`. BetaBinomial variant adds a shared κ estimated iteratively.

**FL LLR** ([`_fl_llr`](src/rigel/calibration/_em.py#L421)):

$$\mathrm{LLR}_{\text{fl},i} = \sum_{x \in \text{FLs at region }i} \log f_{\text{gdna}}(x) - \log f_{\text{rna}}(x)$$

Disabled in iter 0. From iter 1+, `f_gdna` is a γ-weighted histogram of per-region observed FLs; `f_rna` is a (1−γ)-weighted histogram.

**Fusion:**

$$\mathrm{logit}(\gamma_i) = \log \frac{\pi_{\text{soft}}}{1-\pi_{\text{soft}}} + \mathrm{LLR}_{\text{count},i} + \mathrm{LLR}_{\text{strand},i} + \mathrm{LLR}_{\text{fl},i}$$

**No weights.** Hard clipping of `log_odds` to [−500, +500]. No per-channel reliability, no variance modeling, no mis-specification tempering. This is vanilla "sum the LLRs" under the assumption that the models are jointly correctly specified.

### 2.2 The M-step

$$\hat\lambda_G = \frac{\sum_i \gamma_i k^u_i}{\sum_i \gamma_i E_i}\quad\text{(over eligible soft regions)}$$

This is the Poisson MLE for a single-λ model given weighted assignments. **The M-step is not what's broken** — given correct γ, it produces a correct single-λ estimate. The problem is that under hybrid capture, no single λ is correct, so producing the "right" single-λ estimate still yields wrong per-locus priors γ_ℓ.

### 2.3 Prior channel weighting (historical)

Old docs ([docs/METHODS.md](docs/METHODS.md) §5.1, [docs/mappability/stress_test_findings.md](docs/mappability/stress_test_findings.md)) describe a prior system using:

$$W = (2\cdot\text{SS} - 1)^2 \in [0,1]$$

$$\hat\gamma_\ell = W \cdot \hat\gamma_{\text{strand},\ell} + (1-W) \cdot \hat\gamma_{\text{density},\ell}$$

This was grounded as "inverse-variance weighting" but the derivation is heuristic — it's the squared strand-specificity, not the inverse variance of anything. The **current v4 code does not use this** (verified: no `w_strand`, `w_density`, or `(2*ss-1)**2` anywhere in src/rigel/). It was removed when calibration was rewritten into the two-component EM in `_em.py`.

### 2.4 The unified-proposal direction (unimplemented)

[`docs/mappability/unified_proposal.md`](docs/mappability/unified_proposal.md) proposed a joint log-likelihood approach: all three channels as terms in one objective, "self-weighting by Fisher information". The density-channel replacement (two-component Poisson mixture MLE) **was implemented** and is what's in `_em.py` today. The unmapped-read channel and the "Fisher-information self-weighting" fusion were **never implemented**. The current code just sums LLRs.

### 2.5 Summary of the gap

The current system is in a half-finished state:
- Density channel: new two-component Poisson-LogNormal EM. **Misspecified under hybrid capture.**
- Strand channel: upgraded to BetaBinomial with learned κ. **Working well.**
- FL channel: γ-weighted histogram comparison. **Working, but late-starting.**
- Fusion: **unit-weight sum of LLRs.** No Fisher-info weighting; no awareness that one of the three channels may be systematically wrong.

---

## 3. Why the current implementation fails on hybrid-capture data

### 3.1 The density LLR is systematically wrong

Under hybrid capture, for an on-target intronic region with `E_i = 2000 bp`, `k_i = 40`:
- "True" on-target gDNA rate `λ_G_on = 2e-3`. So `λ_G_on · E_i = 4.0`.
- Observed `k_i = 40`, which is consistent with `Poisson(rate = 20e-3 · E_i)` — class R with modest μ.
- Fit λ_G_cal (a compromise value, ≈ 8e-4): `λ_G_cal · E_i = 1.6`.
- `ll_G = 40·log(1.6) − 1.6 − lgamma(41)` ≈ very negative.
- `ll_R` under LogNormal(μ_R ≈ -4.3, σ_R ≈ 1.5) integrates much higher.
- `LLR_count_i ≈ -20`. **High confidence, wrong direction.**

The strand channel for the same region, if `n_sense=20`, `n_anti=20`:
- `ll_G = 20·log(0.5) + 20·log(0.5) = -27.7`
- `ll_R = 20·log(0.19) + 20·log(0.81) = -37.4`
- `LLR_strand_i = +9.7`. **Moderate confidence, correct direction.**

Sum: `-20 + 9.7 = -10.3` → γ_i ≈ 3e-5. **Wrong.**

### 3.2 Why the LLR magnitudes differ so much

Per-region Fisher information ratios:
- Density LLR per region: `∝ k_i · (log(λ_G)/log(μ_R))²` — unbounded in k_i. Scales with **read count**.
- Strand LLR per region: `≤ n_i · KL(0.5 || ss) ≈ n_i · 0.23` for ss=0.81 — bounded by the KL divergence between the two Binomial parameters. Scales with **oriented read count** and **strand separation**.

At ss=0.81 the strand LLR per oriented read is ~0.23 nats. The density LLR per read — under a mis-specified model — can be 1–3 nats. So **every read that enters both channels gets counted 5–10× more strongly in density than in strand**, despite strand being the correct channel in hybrid capture. Unit-weight summation is effectively weighted *against* strand.

### 3.3 Where Fisher information would actually help

Under a **correctly specified** model, LLR per region equals that region's contribution to the Fisher-information-weighted sum. The asymmetry in 3.2 is not a weighting pathology under a correct model — it reflects that density contains more information per read when the model is right. The asymmetry becomes pathological only when density is wrong.

Any fusion scheme that simply rescales unit weights will **not fix** this: the density channel is still wrong, just at smaller amplitude.

---

## 4. Candidate solutions

I group these by theoretical cleanliness.

### 4.1 CLASS A — Fix the density model (most principled)

#### A.1. Targets-BED mode (user's suggested "purest solution")

Accept `--targets targets.bed`. Inside the calibration:

- Set `E_i = 0` for off-target regions in the density channel. Off-target regions contribute zero density likelihood (but can still contribute strand & FL LLRs if oriented).
- Fit a single λ_G over on-target regions only, where uniformity is approximately valid.
- Report separately the off-target intergenic rate for QC.

**Why it works:** Eliminates the non-uniformity. Density is now a valid Poisson-LogNormal mixture over a subset of the genome where its assumptions hold.

**Limitations:**
- Requires user to provide targets BED. Many datasets don't have one.
- Targets are often "aspirational" — real probe capture has off-target bleed. Hard boundaries may leak signal.
- Doesn't help unstranded hybrid-capture libraries directly (but strand channel is useless there anyway, so density-only with targets-BED is the only viable path).

#### A.2. Hierarchical λ_G with region-adjacency smoothing

Model log λ_G as a spatial Gaussian process or MRF over the genome. Each region *i* has its own log `λ_{G,i}` drawn from a smooth prior. On-target regions will get a high local rate, off-target a low rate, learned from the data.

**Why it works:** Lets the data dictate the spatial structure of λ_G. No user-provided BED required; no hard categorization.

**Limitations:**
- Much more expensive. Turns a 1-D Poisson MLE into a Laplace-approximated GP over ~700k regions.
- Additional hyperparameters (smoothing length, prior variance).
- Identifiability risk: the GP can learn to "follow" the expressed regions, which is exactly wrong — we want it to follow the capture pattern, not the expression pattern.

Conclusion: attractive in principle; too heavy for a fix-this-iteration solution.

#### A.3. Per-region self-calibration via orthogonal signals

Use the **strand + FL** channels to identify a clean subset of gDNA-confident regions (γ > 0.95), then fit λ_G from that subset alone, per local neighborhood. This is essentially a "strand-anchored local λ_G".

**Why it works:** Sidesteps density mis-specification by letting strand define "gDNA" first, then letting density's M-step only run on that subset. Per-locus priors then come from "local λ_G × mbp".

**Limitations:**
- Chicken-and-egg: in iter 0, strand is gated off until z ≥ 3, and γ is uniform π_soft. No gDNA-confident subset exists until the EM has converged somewhat.
- Could work as an **EM re-initialization strategy**: run one EM pass, use γ_final to define a strand-confident G subset, re-fit λ_G locally from that subset, then iterate again.

### 4.2 CLASS B — Model density unreliability directly (variance-based fusion)

Instead of down-weighting by a heuristic, **make the density LLR variance-aware** so it self-attenuates when it's wrong.

#### B.1. Dispersion-inflated density channel (Quasi-Poisson / Negative Binomial)

Replace `Poisson(k_i | λ_G E_i)` with `NegBin(k_i | λ_G E_i, φ)` where φ is a learned dispersion. In hybrid capture, the on/off-target mixture produces huge overdispersion at the genome level — φ ≫ 1 is empirically required to fit the data. A well-fit φ naturally widens the density likelihoods for both classes and **reduces the density LLR magnitude**.

**Why it works:** This is how Fisher-info self-weighting actually looks in practice. An over-dispersed likelihood has less information per observation. When the model is mis-specified by unobserved spatial structure, fitting a dispersion parameter absorbs that mis-specification into a larger effective variance, producing smaller LLR magnitudes.

**Analytical intuition:** For a NegBin, the LLR between two rate hypotheses scales with `k / (k/μ + φ)`. As φ → ∞, every region's LLR → 0. The fit will find a φ such that the residuals between predicted and observed regional counts are consistent with the NB — which in hybrid-capture data will inflate φ dramatically.

**Pros:** One extra parameter. Works alongside existing strand + FL channels unchanged. Fully automatic. Theoretically principled — it is Bayesian inference under a **well-specified** NB mixture.

**Cons:** NB overdispersion absorbs any unaccounted structure — not only hybrid-capture targeting, but also e.g. mappability variation, coverage bias, low-complexity regions. This could under-weight density even when it would have been informative.

**Rough validation plan:** Fit NB on the VCaP titration. Expect φ ≈ 1 (standard Poisson) for whole-genome WGS input, φ ≈ 20–100 for MCTP hybrid-capture. Check that the density LLR per region drops commensurately.

#### B.2. Dispersion-inflation learned per-region neighborhood

Same as B.1 but with a locally-learned φ_i over sliding genomic windows. Catches spatially-heterogeneous mis-specification (capture hotspots vs off-target cold zones) without requiring full hierarchical modeling.

**Pros:** Strictly more expressive than B.1.
**Cons:** More parameters; harder to test.

#### B.3. Channel-level temperature learned from agreement

Treat the three channels as noisy experts and fit per-channel temperatures `T_density, T_strand, T_fl`:

$$\text{logit}(\gamma_i) = \text{prior} + \tfrac{\mathrm{LLR}_c}{T_c} + \tfrac{\mathrm{LLR}_s}{T_s} + \tfrac{\mathrm{LLR}_f}{T_f}$$

Fit {T_c, T_s, T_f} to maximize the marginal likelihood (via a held-out calibration set, or via an outer EM). If density is systematically inconsistent with strand & FL, T_density → ∞ automatically.

**Why it works:** Operationalizes the Fisher-info self-weighting from the unified proposal, without requiring us to analytically derive Fisher information for a mis-specified model.

**Pros:** Theoretically defensible (Dawid-Sebastiani calibration; temperature scaling is well-studied in classification). Channel-level granularity is usually sufficient.
**Cons:** Needs a held-out set or an optimization loop. Can be identifiability-brittle if two channels are highly correlated.

### 4.3 CLASS C — Strand-anchored posteriors (user's "strand-first" direction, done right)

The core observation: **the strand channel in strand-specific data is almost-sufficient alone**. The real question is how to combine it with the others *without thresholds*.

#### C.1. Strand-posterior-as-prior (sequential Bayes)

Structure the inference so the strand channel sets the *prior*, and density/FL *update* it.

1. **Step 1 (per region, data-driven prior):** Compute `π_strand_i = P(G | strand observations at i)` from the strand LLR plus a broad strand-specificity-derived prior.
   - For low-coverage regions (few oriented reads) this defaults back toward π_soft, contributing little.
   - For high-coverage regions, this is a confident per-region prior.
2. **Step 2 (density/FL update):** Use `π_strand_i` in place of the global π_soft inside the E-step. Density and FL then move γ_i around this per-region base rate.

**Why it works:** When strand is very confident (high coverage, high SS), it dominates via the prior. When strand is weak (low coverage, low SS), it falls back to the global prior and density takes over — smoothly, by virtue of Binomial's finite KL.

**Why this avoids arbitrary thresholds:** The transition between "strand dominates" and "density dominates" is governed by the relative sharpness of the strand posterior vs density likelihood — no `ss ≥ 0.75` condition.

**Why this is not equivalent to just summing LLRs:** In standard sequential Bayes it *is* identical (LLR additivity is unchanged by which channel you treat as "prior"). However, we can additionally **reduce the density LLR's influence** on γ_i by evaluating it against `π_strand_i`-implied (λ_G, μ_R, σ_R) that have been *fit only over regions strand-identified as class G*.

This is essentially **profile likelihood on the density channel**: the density LLR stops trying to adjudicate class for every region and instead only reports whether a region's (k, E) is consistent with the strand-determined class. For off-target intergenic regions, this will make the density LLR nearly zero (because λ_G is fit to match the strand-identified gDNA regions, which exist across both on- and off-target zones).

**Implementation sketch:**
1. E-step: compute γ only from strand + FL LLRs first.
2. M-step: fit λ_G, μ_R, σ_R using strand+FL γ as weights.
3. Update γ with density included; but the density LLR contribution is now well-calibrated against a λ_G that matches the strand-identified gDNA subset.
4. Iterate.

This is **a different EM ordering**, not a weighting scheme.

#### C.2. Per-region variance-weighted fusion (the theoretically-clean version of "W·strand + (1−W)·density")

For each region i, compute per-channel **log-odds + variance** estimates:

- Strand: `l_s_i ± σ_s_i`. σ_s_i from the Binomial/BetaBinomial SE at (n_i, γ_i).
- Density: `l_d_i ± σ_d_i`. σ_d_i from the Poisson SE (or NB if B.1).
- FL: `l_f_i ± σ_f_i`. σ_f_i from the bootstrap variance of the histogram-comparison statistic.

Fuse by inverse-variance weighting:

$$l_i = \frac{l_s/\sigma_s^2 + l_d/\sigma_d^2 + l_f/\sigma_f^2}{1/\sigma_s^2 + 1/\sigma_d^2 + 1/\sigma_f^2}$$

This is **the optimal linear estimator** (BLUE) under the assumption that each channel's log-odds estimate is unbiased and Gaussian.

**Pros:** Theoretically clean. No thresholds. Naturally reduces density's influence when density is highly uncertain (low read count) without needing SS-based rules.
**Cons:** Assumes unbiasedness, which hybrid-capture violates for density. If density is biased with small variance, inverse-variance weighting gives it *more* weight, not less. So this does not fix the hybrid-capture issue unless combined with B.1.

### 4.4 CLASS D — Detect-then-adapt

Run a **channel-agreement diagnostic** at the start of calibration. If strand and density disagree systematically across regions, conclude that density is mis-specified and switch to a strand-first path automatically.

**Concrete diagnostic:** Compute `γ_strand_only` (strand+FL, no density) and `γ_density_only`. If the Pearson correlation across regions is < 0.5 (or equivalently, if their disagreement rate exceeds what would be expected from independent noisy estimates of the same γ), flag density as untrustworthy and use `γ_strand_only` exclusively.

**Pros:** Fails gracefully; doesn't require the user to provide targets. Purely diagnostic, no model changes.
**Cons:** Needs a threshold after all (the "< 0.5" cutoff). But it's a *diagnostic* threshold, not a fusion threshold — it only determines which fusion path to take, not how to weight channels.

---

## 5. Recommended direction

The user's stated "purest solution" (targets-BED, §4.1 A.1) is correct: it fixes the density channel's mis-specification directly. Every other option in Class A is heavier. But targets-BED requires user input.

For the automatic (no-BED) path, the theoretically-cleanest fix is **Class B.1 (Negative Binomial density channel with learned dispersion φ)**:

- It's a one-parameter extension to the existing Poisson-LogNormal mixture.
- It's Bayesian inference under a **correctly-specified NB mixture**, so no heuristic weights.
- It **automatically** down-weights density when the data is over-dispersed (hybrid capture) and **automatically** keeps density at full strength when data is well-approximated by Poisson (WGS-like).
- It interacts cleanly with the existing strand and FL channels via simple LLR addition.
- It requires no thresholds, no user input, no region classification.

**This is the right answer for the user's stated concern.** The hybrid-capture failure is exactly "unmodeled spatial variance in λ_G manifesting as over-dispersion in the density channel" — and the NB dispersion parameter is the minimal correct model for that.

**Recommended implementation layering:**

1. **NB density channel (B.1) — default automatic path.** Absorbs hybrid-capture mis-specification via learned φ. Falls back to standard Poisson behavior when φ → 1.
2. **Targets-BED mode — power-user path.** Keep Poisson density but restrict to on-target regions. Cleaner when the BED is known; not required.
3. **Strand-anchored re-initialization (C.1) — diagnostic path.** As a convergence improvement: after a first pass, re-fit (λ_G, μ_R, σ_R, φ) using only strand-identified class-G regions. This tightens the density model in exchange for one more EM pass.

The heuristic `W = (2·SS−1)²` blending, explicit channel temperatures, and per-kind λ_G stratification are all **not needed** once the density channel's variance correctly reflects its mis-specification.

### 5.1 Validation plan

For the VCaP titration (8 libraries, 0–80M gDNA spike):
1. Measure learned φ across the titration. Expect φ → 1 for dna00m (pure RNA, Poisson-like) and φ >> 1 for dna10m+ (capture-enrichment dominated).
2. Compare per-region density LLR magnitudes pre- and post-NB. Expect reduction of 5–10× inside transcribed loci.
3. Compare final γ_i to strand-inferred γ_strand_i. Expect correlation to improve from ~0.5 (current) to >0.9 (with NB).
4. Compare quant.gdna_total to strand-decomposed total. Expect ratio to drop from current 1.1–1.5× to ~1.0–1.1×.
5. No regression on synthetic benchmarks (uniform λ_G, no hybrid capture). Expect φ ≈ 1 and identical results to current code.

### 5.2 Out of scope

- The "nRNA siphon" I flagged in the prior document is not a real problem in hybrid capture — the intronic reads really are ~94% gDNA, and the locus EM's final gDNA allocation is within 1.5× of strand-inferred truth. Fixing calibration's per-locus γ_ℓ should resolve any residual 1.5× overshoot automatically. No scorer or locus-EM changes required.

---

## 6. Open questions for the user

1. **Accept B.1 (NB density channel) as the primary direction?** It's the theoretically-cleanest fix that doesn't require user input. I'd also want to implement §4.4 (D) as a sanity diagnostic so if NB still disagrees with strand we know early.
2. **Do we want targets-BED mode as a parallel feature?** It's a separate (maybe simpler) code path. Worth having for users who do have kit BEDs.
3. **Is there any appetite for §4.1 A.3 (strand-anchored λ_G reinitialization)** as a belt-and-suspenders addition? It's cheap and robust but adds iterations.
4. **Should we retire the `(2·SS−1)²` weighting from [METHODS.md](docs/METHODS.md)?** It's documented but not implemented. Either re-implement it in a theoretically-clean form (inverse-variance, per-region) or delete from docs.

---

## 7. References (in this repo)

- Calibration EM: [src/rigel/calibration/_em.py](src/rigel/calibration/_em.py), especially `_e_step` (fusion) and `_m_step` (λ_G update).
- Orchestration: [src/rigel/calibration/_calibrate.py](src/rigel/calibration/_calibrate.py).
- Region partition: [src/rigel/index.py](src/rigel/index.py) `build_region_table`.
- Historical W = (2·SS−1)² weighting: [docs/METHODS.md](docs/METHODS.md) §5.1.
- Unified proposal (density replaced; fusion never implemented): [docs/mappability/unified_proposal.md](docs/mappability/unified_proposal.md).
- Stress-test findings (SS-gated blending, P₁₀ zero collapse): [docs/mappability/stress_test_findings.md](docs/mappability/stress_test_findings.md).
- This session's diagnostic data: `/scratch/.../runs/human_optionb_diag/dna*m/region_counts.feather` (8 libs).
- Superseded by this doc: [docs/calibration/strand_first_hybrid_capture_plan.md](docs/calibration/strand_first_hybrid_capture_plan.md), [docs/calibration/root_cause_locus_em_nrna_siphon.md](docs/calibration/root_cause_locus_em_nrna_siphon.md).
