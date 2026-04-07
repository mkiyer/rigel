# gDNA vs nRNA Identifiability at Low Strand Specificity

## Problem Statement

At SS=0.50 (unstranded data) or very low strand specificity (SS≤0.65), the EM cannot distinguish gDNA contamination from nascent RNA (nRNA). This document analyzes why, and proposes solutions.

## Current Scoring Architecture

| Signal | mRNA | nRNA | gDNA |
|--------|------|------|------|
| **Strand** | `log(SS)` sense / `log(1-SS)` anti | same as mRNA | `LOG_HALF` always |
| **FL model** | RNA FL (Dirichlet prior) | RNA FL (same) | gDNA FL (Dirichlet prior) |
| **Splicing** | spliced → full LL | unspliced only | unspliced only, `gdna_log_sp=-∞` if spliced |
| **Bias correction** | `-log(exonic_len - frag + 1)` | `-log(genomic_span - frag + 1)` | `-log(locus_span - frag + 1)` |

nRNA components are synthetic single-exon transcripts (genomic span = exonic length). They use the **same** scoring function as mRNA — same RNA FL model, same strand model. Only during pruning are they separated into independent pools.

## Degeneracy at SS=0.50

When SS=0.50, `log(SS) = log(0.5) = LOG_HALF`:

- **nRNA strand LL = `LOG_HALF`** (sense and antisense equally likely)
- **gDNA strand LL = `LOG_HALF`** (always)
- **Strand is identical** — zero discrimination

FL models: with the Dirichlet prior and SS=0.50, calibration uses the pure density pathway (blending weight `w = (2·SS−1)² = 0`). The gDNA FL model is trained from low-density regions. If gDNA contamination is genuinely low, low-density regions are mostly very-low-expression RNA — so the gDNA FL model approximates the RNA FL model. **FL discrimination is weak to zero.**

That leaves **bias correction** as the sole per-fragment differentiator:

- nRNA: `-log(transcript_genomic_span - frag + 1)`
- gDNA: `-log(locus_span - frag + 1)`

**For a single-transcript locus**: `locus_span ≈ transcript_genomic_span` → bias correction is identical → **complete degeneracy**. The EM cannot separate gDNA from nRNA. The fragment posterior is determined entirely by the prior ratio `α_gDNA / α_nRNA`.

**For a multi-transcript locus** where transcripts span different regions: `locus_span > any single transcript's genomic_span`, so gDNA has a larger effective length and gets a small per-fragment penalty. But this is often marginal (e.g., 50kb vs 45kb → ~0.1 log-unit difference).

## Why the Prior is Too Weak

The prior sets `α_gDNA = γ × C` and `α_RNA = (1−γ) × C` with `C = 1.0`. Even if `γ = 0.01` (1% contamination), the prior is:

- `α_gDNA = 0.01`
- `α_RNA = 0.99`

That's 1 pseudo-observation split across all components. With N=2000 real fragments, the prior has 0.05% influence. The EM simply ignores it.

## The Core Insight: gDNA is Global but Modeled Locally

The calibration correctly estimates `λ_G` (gDNA density per bp) from global data — cross-locus information. But then this global estimate is converted to a per-locus prior with `C=1.0`, losing all its statistical weight.

Inside a highly-expressed locus with abundant nRNA:

- Local unspliced density = `λ_nRNA + λ_G` (high)
- Globally calibrated `λ_G` is low
- But the weak prior can't enforce the global constraint → the local gDNA component is free to absorb nRNA fragments

The gDNA component acts as a "sponge" for any fragments that are equally plausible under both gDNA and nRNA likelihoods.

## When IS gDNA Separable from nRNA?

1. **Stranded data (SS > 0.7)**: Anti-sense fragments are preferentially gDNA. Strand signal breaks the degeneracy.

2. **Different FL distributions**: If gDNA genuinely has a different FL profile (e.g., FFPE degradation → shorter fragments), the FL model discriminates. Requires sufficient gDNA data for calibration.

3. **Spliced fragments**: Definitively mRNA (not gDNA or nRNA). Once mRNA is separated, the spliced-to-unspliced ratio constrains the mRNA-nRNA split, indirectly constraining the gDNA residual.

4. **Multi-transcript loci with diverse spans**: gDNA's `locus_span` exceeds any single nRNA transcript's span → effective length difference → per-fragment likelihood difference.

5. **Cross-locus density consistency** (NOT currently used): gDNA should produce approximately uniform density everywhere. A locus with 100× the global unspliced density is almost certainly nRNA, not gDNA.

---

## Proposed Solutions

### Idea 1: Scale C Proportional to Locus Fragment Count

Instead of `C = 1.0`, set `C = N_locus` (the number of fragments in the locus):

$$\alpha_{\text{gDNA}} = \gamma \times N_{\text{locus}}, \quad \alpha_{\text{RNA}} = (1 - \gamma) \times N_{\text{locus}}$$

This makes the prior competitive with the data. With 2000 fragments and `γ = 0.01`, the prior says "expect 20 gDNA fragments" with strength proportional to the data size. The EM can still deviate if the evidence is strong, but for degenerate components (gDNA ≈ nRNA), the prior holds.

**Pro**: Simple, uses existing architecture. The calibrated fraction γ is well-estimated from thousands of regions.

**Con**: If γ is wrong (poorly calibrated), the error propagates proportionally to N. Also makes the prior dominate even when the data has genuine information (e.g., a locus with a true local spike in gDNA contamination — rare but possible near CNVs).

### Idea 2: Fix gDNA at Calibrated Level (Not a Free EM Parameter)

Instead of letting the EM find gDNA abundance, set `θ_gDNA = λ_G × locus_span` as a fixed constant, exclude gDNA from EM iteration, and subtract its expected contribution before allocating the rest to mRNA/nRNA.

**Pro**: Eliminates the degeneracy entirely. Uses the global calibration as a hard constraint.

**Con**: Can't capture local gDNA variation. If a region has genuinely higher gDNA (e.g., amplified region), it will be attributed to RNA.

### Idea 3: Informative Prior with ESS from Calibration Confidence

The calibration estimates γ from N_regions with N_total fragments. The ESS of this estimate is roughly `N_total / variance(γ)`. Use this as the prior strength:

$$C = \text{ESS}_{\text{calibration}}, \quad \alpha_{\text{gDNA}} = \gamma \times C$$

This is a principled Empirical Bayes approach: the prior's strength reflects how confident we are in the calibrated fraction. With millions of fragments across thousands of regions, C could be in the hundreds or thousands — strong enough to constrain gDNA but not hard-fixed.

**Pro**: Principled, data-driven, no magic numbers. Works on a spectrum between "very weak prior" (poor calibration) and "near-fixed" (strong calibration).

**Con**: Need to correctly estimate the calibration ESS. The calibration is a regression, not a simple proportion — extracting a meaningful ESS requires care.

### Idea 4: Explicitly Differentiate gDNA Effective Length

Force the gDNA component to use a larger effective length than nRNA by extending `locus_span` to include flanking intergenic sequence (e.g., ±5kb). Since gDNA genuinely can produce fragments from anywhere in the genome (including intergenic regions adjacent to the locus), a wider span is physically motivated.

This creates a per-fragment likelihood asymmetry:

```
bias_correction(gDNA) = -log(span_ext - frag + 1) < bias_correction(nRNA) = -log(span - frag + 1)
```

Each nRNA fragment gets a higher per-fragment probability under nRNA than gDNA because nRNA's effective length is smaller. Over many fragments, this accumulates.

**Pro**: Physically motivated, creates discrimination even in the degenerate case.

**Con**: The amount of flanking extension is somewhat arbitrary. Small extensions give minimal discrimination. Large extensions may over-penalize gDNA. Also, if we extend gDNA's span, fragments outside the transcript body but within the extended region have zero nRNA likelihood and would be forced to gDNA — which is correct and actually helps.

### Idea 5: Two-Stage Decomposition

Stage 1: Separate spliced (mRNA) from unspliced (nRNA + gDNA). This is well-posed — splicing is a definitive signal.

Stage 2: Decompose the unspliced pool into nRNA + gDNA using the calibrated density as a *fixed constraint*:

```
N_gDNA_locus = λ_G × locus_span
N_nRNA_locus = N_unspliced - N_gDNA_locus
```

This avoids the identifiability problem by not trying to solve it with EM. Instead, it uses the global calibration deterministically for the degenerate component pair.

**Pro**: Clean conceptual separation of what's identifiable (splicing) from what isn't (gDNA vs nRNA for unspliced).

**Con**: Loses the joint estimation benefit — mRNA/nRNA assignments are correlated through shared fragments, and a two-stage approach doesn't account for this.

### Idea 6: Hierarchical gDNA Model

Model gDNA density as a single global parameter `λ_G` shared across all loci, estimated jointly:

$$\text{E-step}: \quad P(\text{gDNA} \mid \text{fragment}) \propto \lambda_G \times f(\text{fragment} \mid \text{gDNA})$$

$$\text{M-step}: \quad \lambda_G = \frac{\sum_{\text{loci}} N_{\text{gDNA, locus}}}{\sum_{\text{loci}} \text{span}_{\text{locus}}}$$

The key insight: gDNA is the ONLY component with a global parameter (density is genome-wide constant). mRNA and nRNA are per-transcript. This asymmetry provides identification — the model can distinguish gDNA from nRNA because gDNA's per-locus abundance is constrained to be proportional to locus span (global density × local span), while nRNA abundance varies freely by expression level.

**Pro**: Statistically correct. Uses all cross-locus information. Naturally handles the identifiability problem.

**Con**: Requires fundamentally restructuring the EM to have a global M-step for `λ_G` across all loci. Currently each locus runs independently. Would need an outer EM loop: (1) run per-locus EM with fixed `λ_G`, (2) update `λ_G` from across-locus gDNA totals, (3) repeat.

---

## Assessment

**Ideas 1–2** are quick fixes but have real limitations (Idea 1: γ error propagation; Idea 2: no local variation).

**Idea 3** (ESS from calibration confidence) is the most natural extension of the current Empirical Bayes framework. The calibration already estimates γ — encoding its confidence as prior strength is the principled next step.

**Idea 4** (effective length differentiation) is worth exploring because it creates likelihood discrimination independent of the prior. It's complementary to prior-based approaches.

**Idea 6** (hierarchical model) is the correct theoretical answer but the biggest architectural change. An outer-inner EM loop is feasible — the inner per-locus EM already converges fast, and the outer `λ_G` update is just a ratio of sums. But it's a significant implementation.

A practical middle ground may be **Idea 3 + Idea 4**: use an informative prior (ESS from calibration) to constrain gDNA globally, plus effective length differentiation to give the EM some likelihood signal even in the degenerate case. Together, these should make gDNA vs nRNA separation robust even at SS=0.50 for loci where nRNA >> gDNA.
